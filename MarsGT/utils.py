import anndata as ad
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from collections import Counter
from tqdm import tqdm
import math
import scanpy as sc
from sklearn.metrics import accuracy_score
from sklearn.metrics.cluster import normalized_mutual_info_score

def subgraph(graph, seed, n_neighbors, node_sele_prob):
    total_matrix_size = 1 + np.cumprod(n_neighbors).sum()  # Number of nodes in the subgraph
    picked_nodes = {seed}  # One node in the batch
    last_layer_nodes = {seed}

    # Number of nodes selected in each layer. Initially, only the seed node is selected.
    to_pick = 1
    for n_neighbors_current in n_neighbors:  # Current layer neighbors
        to_pick = to_pick * n_neighbors_current
        neighbors = graph[list(last_layer_nodes), :].nonzero()[1]  # Find neighbors of last_layer_nodes

        neighbors_prob = node_sele_prob[list(neighbors)]
        neighbors = list(set(neighbors))  # Make all nodes from the last layer part of the neighbors set
        n_neigbors_real = min(
            to_pick,
            len(neighbors))  # Handle the case where the required number of neighbors is less than the actual number of neighbors
        if len(neighbors_prob) == 0:
            continue
        last_layer_nodes = set(
            np.random.choice(neighbors, n_neigbors_real, replace=False,
                             p=softmax(neighbors_prob)))  # Select non-repeated nodes from neighbors
        picked_nodes |= last_layer_nodes  # Update picked_nodes as last_layer_nodes ∪ picked_nodes
    indices = list(sorted(picked_nodes - {seed}))
    return indices


# ---------------------------------------------------------------------------------------------------------------
# Device-side feature rows.
# Every batch needs dense rows of RNA_matrix / ATAC_matrix (features x ALL cells). The original code densified them
# on the CPU and copied ~0.8 GB per batch to the GPU at full dataset width. Instead, keep one sparse (CSR) copy of
# each matrix on the device and densify the selected rows there. Values are identical (float32 counts).
#   GPU_FEATURES = 'auto' -> use the device cache when the device is CUDA (default)
#                  True   -> always (also on CPU; mainly for testing)
#                  False  -> original CPU path
# ---------------------------------------------------------------------------------------------------------------
GPU_FEATURES = 'auto'
_ROW_CACHE = {}


class DeviceRows:
    """Sparse CSR copy of a (features x cells) scipy matrix on `device`; .rows(idx) returns dense float32 rows."""

    def __init__(self, M, device):
        import scipy.sparse as _sp
        # The original M[idx,].todense() is column-major (Fortran order) when M is CSC. Matmul rounding depends on
        # memory layout, so return the same layout to reproduce the original arithmetic exactly.
        self.fortran = getattr(M, 'format', None) == 'csc'
        M = _sp.csr_matrix(M, dtype=np.float32, copy=True)
        M.sum_duplicates(); M.sort_indices()
        self.shape, self.device = M.shape, device
        self.indptr = torch.as_tensor(M.indptr.astype(np.int64), device=device)
        self.indices = torch.as_tensor(M.indices.astype(np.int32), device=device)
        self.data = torch.as_tensor(M.data, device=device)

    def rows(self, idx):
        idx = torch.as_tensor(np.asarray(list(idx), dtype=np.int64), device=self.device)
        if self.fortran:   # (n, C) view with strides (1, n): same layout as the original, no extra memory
            out = torch.zeros((self.shape[1], idx.numel()), dtype=torch.float32, device=self.device).t()
        else:
            out = torch.zeros((idx.numel(), self.shape[1]), dtype=torch.float32, device=self.device)
        start = self.indptr[idx]
        cnt = self.indptr[idx + 1] - start
        total = int(cnt.sum())
        if total == 0:
            return out
        row = torch.repeat_interleave(torch.arange(idx.numel(), device=self.device), cnt)
        first = torch.repeat_interleave(torch.cumsum(cnt, 0) - cnt, cnt)
        pos = torch.repeat_interleave(start, cnt) + (torch.arange(total, device=self.device) - first)
        out[row, self.indices[pos].long()] = self.data[pos]
        return out

    def nbytes(self):
        return sum(t.numel() * t.element_size() for t in (self.indptr, self.indices, self.data))


def _use_device_rows(device):
    if GPU_FEATURES == 'auto':
        return torch.device(device).type == 'cuda'
    return bool(GPU_FEATURES)


def dense_rows(M, idx, device):
    """Dense float32 tensor of rows `idx` of scipy matrix M, on `device`. Same values as the original
    torch.tensor(np.array(M[list(idx),].todense()), dtype=torch.float32).to(device)."""
    if not _use_device_rows(device):
        return torch.tensor(np.array(M[list(idx), ].todense()), dtype=torch.float32).to(device)
    key = (id(M), str(device))
    hit = _ROW_CACHE.get(key)
    if hit is None or hit[0] is not M or hit[1] != (M.shape, M.nnz):
        hit = _ROW_CACHE[key] = (M, (M.shape, M.nnz), DeviceRows(M, device))
    return hit[2].rows(idx)


def clear_feature_cache():
    """Free the device copies of RNA/ATAC matrices (e.g. between datasets)."""
    _ROW_CACHE.clear()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


SAMPLING_MODES = ('upstream', 'aligned', 'paper')


def feature_probs(counts, totals, mode, is_rna, quantile=0.25):
    """Per-cell sampling probabilities over that cell's nonzero features.

    counts : 1-D array, the cell's counts for its nonzero features (genes or peaks)
    totals : 1-D array, each of those features' total count across all cells
    mode   : 'aligned' -> upstream weights (RNA: log(x+1) after zeroing x<5; ATAC: log(x+1)) passed
                          through softmax, i.e. proportional to x+1 (or 1), but kept aligned to features
             'paper'   -> Wang et al. 2024 (Nat Commun 15:338, doi:10.1038/s41467-023-44570-8) eqs. 2-5:
                          Prop = x_{i,k0} / sum_k x_{i,k}   (share of the feature's total found in this cell)
                          RNA: keep genes with x > a, a = `quantile` of the cell's nonzero counts (falls back
                               to x >= a if that leaves nothing); Prob = Prop / sum(Prop over kept genes)
                          ATAC: Prob = Prop / sum(Prop)
    Returns probabilities (same length as counts; 0 for filtered features).
    """
    counts = np.asarray(counts, dtype=np.float64)
    if counts.size == 0:
        return counts
    if mode == 'aligned':
        w = np.log(np.where(counts < 5, 0, counts) + 1) if is_rna else np.log(counts + 1)
        w = np.exp(w - w.max())
        return w / w.sum()
    if mode == 'paper':
        prop = counts / np.asarray(totals, dtype=np.float64)
        if is_rna:
            a = np.quantile(counts, quantile)
            keep = counts > a
            if not keep.any():
                keep = counts >= a
            prop = np.where(keep, prop, 0.0)
        return prop / prop.sum()
    raise ValueError(f'unknown sampling mode {mode!r}; choose from {SAMPLING_MODES}')


def _sample_cell(matrix_csc, totals, node, n, mode, is_rna, quantile):
    col = matrix_csc[:, node]
    idx, counts = col.indices, col.data
    if idx.size == 0:
        return []
    p = feature_probs(counts, totals[idx], mode, is_rna, quantile)
    k = min(n, int((p > 0).sum()))
    return sorted(np.random.choice(idx, k, replace=False, p=p).tolist())


def _batch_select_weighted(RNA_matrix, ATAC_matrix, neighbor, cell_size, sampling, quantile):
    if len(neighbor) != 1:
        raise ValueError('aligned/paper sampling supports a single hop only, e.g. neighbor=[20]')
    n = int(neighbor[0])
    R, A = RNA_matrix.tocsc(), ATAC_matrix.tocsc()
    R.sort_indices(); A.sort_indices()
    r_tot = np.asarray(R.sum(axis=1)).ravel()
    a_tot = np.asarray(A.sum(axis=1)).ravel()
    node_ids = np.random.choice(R.shape[1], size=R.shape[1], replace=False)
    n_batch = math.ceil(node_ids.shape[0] / cell_size)
    indices_ss, dic = [], {}
    for i in tqdm(range(n_batch)):
        cells = node_ids[i * cell_size:(i + 1) * cell_size]
        g_all, p_all = set(), set()
        for node in cells:
            g = _sample_cell(R, r_tot, node, n, sampling, True, quantile)
            p = _sample_cell(A, a_tot, node, n, sampling, False, quantile)
            dic[node] = {'g': g, 'p': p}
            g_all.update(g); p_all.update(p)
        indices_ss.append({'gene_index': sorted(g_all), 'peak_index': sorted(p_all), 'cell_index': cells})
    return indices_ss, node_ids, dic


def batch_select_whole(RNA_matrix, ATAC_matrix, neighbor=[20], cell_size=30, sampling='upstream', quantile=0.25):
    """Partition cells into batches and sample each cell's gene/peak neighbors.

    sampling : 'upstream' (default) -> original MarsGT 0.2.1 code, unchanged. Note its probabilities are
                                       misaligned with features (list(set(...)) reorders the features after
                                       the probabilities were computed), so weights are effectively shuffled.
               'aligned'            -> the upstream weighting with the alignment bug fixed.
               'paper'              -> the specificity-based scheme described in the MarsGT paper (eqs. 2-5).
    quantile : 'paper' only; gene filter threshold as a quantile of the cell's nonzero counts (paper: first quartile).
    """
    if sampling not in SAMPLING_MODES:
        raise ValueError(f'unknown sampling mode {sampling!r}; choose from {SAMPLING_MODES}')
    if sampling != 'upstream':
        print(f'Partitioning cells into batches ({sampling} sampling).')
        return _batch_select_weighted(RNA_matrix, ATAC_matrix, neighbor, cell_size, sampling, quantile)
    print('We are currently in the process of partitioning the data into batches. Kindly wait for a moment, please.')
    node_ids = np.random.choice(RNA_matrix.shape[1], size=RNA_matrix.shape[1], replace=False)
    n_batch = math.ceil(node_ids.shape[0] / cell_size)
    indices_ss = []

    RNA_matrix1 = RNA_matrix
    dic = {}
    for i in tqdm(range(n_batch)):
        gene_indices_all = []
        peak_indices_all = []
        if i < n_batch:
            for index, node in enumerate(node_ids[i * cell_size:(i + 1) * cell_size]):
                rna_ = RNA_matrix1[:, node].todense()
                rna_[rna_ < 5] = 0
                gene_indices = subgraph(RNA_matrix.transpose(), node, neighbor, np.squeeze(np.array(np.log(rna_ + 1))))
                peak_indices = subgraph(ATAC_matrix.transpose(), node, neighbor,
                                        np.squeeze(np.array(np.log(ATAC_matrix[:, node].todense() + 1))))
                dic[node] = {'g': gene_indices, 'p': peak_indices}
                gene_indices_all = gene_indices_all + gene_indices
                peak_indices_all = peak_indices_all + peak_indices
            node_indices_all = node_ids[i * cell_size:(i + 1) * cell_size]
        else:
            for index, node in enumerate(node_ids[i * cell_size:]):
                rna_ = RNA_matrix1[:, node].todense()
                rna_[rna_ < 5] = 0
                gene_indices = subgraph(RNA_matrix.transpose(), node, neighbor,
                                        np.squeeze(np.array(np.log(rna_[:, node].todense() + 1))))
                peak_indices = subgraph(ATAC_matrix.transpose(), node, neighbor,
                                        np.squeeze(np.array(np.log(ATAC_matrix[:, node].todense() + 1))))
                dic[node] = {'g': gene_indices, 'p': peak_indices}
                gene_indices_all = gene_indices_all + gene_indices
                peak_indices_all = peak_indices_all + peak_indices
            node_indices_all = node_ids[i * cell_size:]

        gene_indices_all = list(set(gene_indices_all))
        peak_indices_all = list(set(peak_indices_all))
        h = dict()
        h['gene_index'] = gene_indices_all
        h['peak_index'] = peak_indices_all
        h['cell_index'] = node_indices_all
        indices_ss.append(h)
    return indices_ss, node_ids, dic


def softmax(x):
    return (np.exp(x) / np.exp(x).sum())


class LabelSmoothing(nn.Module):
    """NLL loss with label smoothing.
    """

    def __init__(self, smoothing=0.0):
        """Constructor for LabelSmoothing module.
        :param smoothing: Label smoothing factor
        """
        super(LabelSmoothing, self).__init__()
        self.confidence = 1.0 - smoothing
        self.smoothing = smoothing

    def forward(self, x, target):
        logprobs = torch.nn.functional.log_softmax(x, dim=-1)
        nll_loss = -logprobs.gather(dim=-1, index=target.unsqueeze(1))
        nll_loss = nll_loss.squeeze(1)
        smooth_loss = -logprobs.mean(dim=-1)
        loss = self.confidence * nll_loss + self.smoothing * smooth_loss
        return loss.mean()


def initial_clustering(RNA_matrix, custom_n_neighbors=None, n_pcs=40, custom_resolution=None, use_rep=None):
    print(
        '\tWhen the number of cells is less than or equal to 500, it is recommended to set the resolution value to 0.2.')
    print('\tWhen the number of cells is within the range of 500 to 5000, the resolution value should be set to 0.5.')
    print('\tWhen the number of cells is greater than 5000, the resolution value should be set to 0.8.')

    def segment_function(x):
        if x <= 500:
            return 0.2, 5
        elif x <= 5000:
            return 0.5, 10
        else:
            return 0.8, 15

    adata = ad.AnnData(RNA_matrix.transpose().astype('float32'))

    # If the user did not provide a custom resolution or n_neighbors value, use the values calculated by segment_function
    if custom_resolution is None or custom_n_neighbors is None:
        resolution, n_neighbors = segment_function(adata.shape[0])
    else:
        resolution = custom_resolution
        n_neighbors = custom_n_neighbors

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Use the user-provided embedding if available, otherwise use n_pcs
    if use_rep is not None:
        adata.obsm['use_rep']=use_rep
        sc.pp.neighbors(adata, use_rep='use_rep', n_neighbors=n_neighbors)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.leiden(adata, resolution)
    return adata.obs['leiden']


def purity_score(y_true, y_pred):
    """Purity score

    Args:
        y_true (np.ndarray): n*1 matrix, true labels
        y_pred (np.ndarray): n*1 matrix, predicted clusters

    Returns:
        float: Purity score
    """
    # Create a matrix to store the majority-voted labels
    y_voted_labels = np.zeros(y_true.shape)

    # Sort the labels
    # Some labels might be missing, e.g., a set {0,2} where 1 is missing
    # First, find the unique labels and then map them to an ordered set
    # E.g., {0,2} should be mapped to {0,1}
    labels = np.unique(y_true)
    ordered_labels = np.arange(labels.shape[0])
    for k in range(labels.shape[0]):
        y_true[y_true == labels[k]] = ordered_labels[k]
    y_true = np.array(y_true, dtype='int64')

    # Update the unique labels
    labels = np.unique(y_true)

    # Set the number of bins to n_classes + 2 so that we can compute the actual
    # class occurrences between two consecutive bins
    # The larger bin is excluded: [bin_i, bin_i+1[
    bins = np.concatenate((labels, [np.max(labels) + 1]), axis=0)

    for cluster in np.unique(y_pred):
        hist, _ = np.histogram(y_true[y_pred == cluster], bins=bins)
        # Find the most frequent label in the cluster
        winner = np.argmax(hist)
        y_voted_labels[y_pred == cluster] = winner

    y_true = np.array(y_true, dtype='int8')
    y_voted_labels = np.array(y_voted_labels, dtype='int8')
    return accuracy_score(y_true, y_voted_labels), y_true


def Entropy(pred_label, true_label):
    e = 0
    for k in set(pred_label):
        en = 0
        pred_k = Counter(pred_label)[k]
        index_pred_k = pred_label == k
        for j in set(true_label):
            true_j = Counter(true_label)[j]
            intersection_kj = (true_label[index_pred_k] == j).sum()
            p = np.array(intersection_kj) / np.array(pred_k)
            if p != 0:
                en += np.log(p) * p
        e = e + en * pred_k / true_label.shape[0]
    return abs(e)
