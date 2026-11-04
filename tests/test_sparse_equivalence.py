"""Vectorized peak_Sparse / gene_Sparse must match the original loop implementations exactly.

Run:  python -m pytest tests/  (or: python tests/test_sparse_equivalence.py)
"""
import itertools
import torch

from MarsGT.egrn import peak_Sparse, gene_Sparse


# --- reference: upstream MarsGT 0.2.1 loop implementations (legacy constructor swapped for sparse_coo_tensor) ---
def peak_Sparse_ref(peak_feature, gene_feature, device):
    peak_row = torch.where(peak_feature != 0)[0].to(device)
    peak_col = torch.where(peak_feature != 0)[1].to(device)
    peak_data = peak_feature[torch.where(peak_feature != 0)]
    peak_combine_row = torch.zeros((peak_row.shape[0] * gene_feature.shape[0]), dtype=torch.int32)
    peak_combine_col = torch.zeros((peak_row.shape[0] * gene_feature.shape[0]), dtype=torch.int32)
    peak_combine_data = torch.zeros((peak_row.shape[0] * gene_feature.shape[0]))
    non_count = 0
    for peak_num in range(peak_feature.shape[0]):
        peak_sum = len(peak_row[peak_row == peak_num])
        peak_init = torch.zeros((gene_feature.shape[0], peak_sum), dtype=torch.int32).to(device)
        for id in range(gene_feature.shape[0]):
            peak_init[id,] = id + peak_num * gene_feature.shape[0]
        n = peak_sum * gene_feature.shape[0]
        peak_combine_row[non_count:non_count + n] = (peak_row[peak_row == peak_num] + peak_init - peak_num).reshape(1, -1)
        peak_combine_col[non_count:non_count + n] = (peak_col[peak_row == peak_num]).repeat(gene_feature.shape[0])
        peak_combine_data[non_count:non_count + n] = (peak_data[peak_row == peak_num]).repeat(gene_feature.shape[0])
        non_count += n
    return torch.sparse_coo_tensor(torch.vstack((peak_combine_row, peak_combine_col)).long(), peak_combine_data,
                                   torch.Size((gene_feature.shape[0] * peak_feature.shape[0], peak_feature.shape[1]))).to(device)


def gene_Sparse_ref(peak_feature, gene_feature, device):
    gene_row = torch.where(gene_feature != 0)[0].to(device)
    gene_col = torch.where(gene_feature != 0)[1].to(device)
    gene_data = gene_feature[torch.where(gene_feature != 0)]
    gene_init = torch.zeros((peak_feature.shape[0], len(gene_row)), dtype=torch.int32).to(device)
    for peak_id in range(peak_feature.shape[0]):
        gene_init[peak_id] = gene_row + peak_id * gene_feature.shape[0]
    gene_combine_row = gene_init.reshape(1, -1)[0]
    gene_combine_col = gene_col.repeat(1, peak_feature.shape[0])[0]
    gene_combine_data = gene_data.repeat(peak_feature.shape[0])
    return torch.sparse_coo_tensor(torch.vstack((gene_combine_row, gene_combine_col)), gene_combine_data,
                                   torch.Size((gene_feature.shape[0] * peak_feature.shape[0], gene_feature.shape[1]))).to(device)


def _features(G, P, C, density, seed):
    g = torch.Generator().manual_seed(seed)
    gf = torch.poisson(torch.full((G, C), 2.0), generator=g) * (torch.rand(G, C, generator=g) < density)
    pf = torch.poisson(torch.full((P, C), 1.0), generator=g) * (torch.rand(P, C, generator=g) < density)
    if P > 1: pf[0] = 0          # empty peak row
    if G > 1: gf[-1] = 0         # empty gene row
    return pf, gf


SIZES = [(1, 1, 1), (5, 4, 3), (40, 37, 30), (120, 150, 30), (465, 590, 30)]   # last ~ a real MarsGT batch
DENSITIES = [0.0, 0.05, 0.3, 1.0]


def test_vectorized_matches_reference():
    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    for (G, P, C), dens in itertools.product(SIZES, DENSITIES):
        pf, gf = (t.to(dev) for t in _features(G, P, C, dens, seed=G * 1000 + P))
        rp, rg = peak_Sparse_ref(pf, gf, dev), gene_Sparse_ref(pf, gf, dev)
        vp, vg = peak_Sparse(pf, gf, dev), gene_Sparse(pf, gf, dev)
        assert rp.shape == vp.shape and rg.shape == vg.shape, (G, P, C, dens)
        assert torch.equal(rp.to_dense(), vp.to_dense()), ('peak_Sparse', G, P, C, dens)
        assert torch.equal(rg.to_dense(), vg.to_dense()), ('gene_Sparse', G, P, C, dens)
        # the quantity MarsGT actually uses downstream
        assert torch.equal(rg.mul(rp).to_dense(), vg.mul(vp).to_dense()), ('product', G, P, C, dens)


if __name__ == '__main__':
    test_vectorized_matches_reference()
    print('vectorized peak_Sparse/gene_Sparse identical to reference')
