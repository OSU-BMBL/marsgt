"""Sub-sampling modes of batch_select_whole.

Run:  python -m pytest tests/   (or: python tests/test_sampling.py)
"""
import numpy as np
import scipy.sparse as sp

from MarsGT.utils import batch_select_whole, feature_probs, subgraph


def _toy(seed=0, n_genes=300, n_peaks=800, n_cells=90):
    rng = np.random.default_rng(seed)
    def m(r, c, d, lam):
        return sp.random(r, c, density=d, format='csc', random_state=rng,
                         data_rvs=lambda n: rng.poisson(lam, n) + 1).astype(np.float32)
    return m(n_genes, n_cells, 0.2, 4), m(n_peaks, n_cells, 0.05, 1)


def test_upstream_mode_is_unchanged():
    """sampling='upstream' must reproduce the original 0.2.1 loop (reimplemented here verbatim) draw for draw."""
    R, A = _toy()

    def original(RNA_matrix, ATAC_matrix, neighbor, cell_size):
        node_ids = np.random.choice(RNA_matrix.shape[1], size=RNA_matrix.shape[1], replace=False)
        n_batch = int(np.ceil(node_ids.shape[0] / cell_size)); out = []
        for i in range(n_batch):
            g_all, p_all = [], []
            for node in node_ids[i * cell_size:(i + 1) * cell_size]:
                rna_ = RNA_matrix[:, node].todense(); rna_[rna_ < 5] = 0
                g_all += subgraph(RNA_matrix.transpose(), node, neighbor, np.squeeze(np.array(np.log(rna_ + 1))))
                p_all += subgraph(ATAC_matrix.transpose(), node, neighbor,
                                  np.squeeze(np.array(np.log(ATAC_matrix[:, node].todense() + 1))))
            out.append((sorted(set(g_all)), sorted(set(p_all)), list(node_ids[i * cell_size:(i + 1) * cell_size])))
        return out

    np.random.seed(7); ref = original(R, A, [20], 30)
    np.random.seed(7); new, _, _ = batch_select_whole(R, A, neighbor=[20], cell_size=30)
    got = [(sorted(b['gene_index']), sorted(b['peak_index']), list(b['cell_index'])) for b in new]
    assert got == ref


def test_upstream_probabilities_are_misaligned_and_aligned_fixes_it():
    """One gene with a huge count among many count-1 genes. Its intended weight dominates, so with n=1
    'aligned' should pick it almost always; 'upstream' attaches that weight to a different gene."""
    n_genes, n_cells = 3000, 2
    rows = np.arange(0, 3000, 5)                         # 600 expressed genes, spread over the index range
    data = np.ones(rows.size, np.float32); big = rows[-1]; data[-1] = 1000
    R = sp.csc_matrix((np.r_[data, data], (np.r_[rows, rows], np.r_[np.zeros(rows.size), np.ones(rows.size)])),
                      shape=(n_genes, n_cells))
    A = sp.csc_matrix(np.ones((10, n_cells), np.float32))
    hits = {'upstream': 0, 'aligned': 0}
    for mode in hits:
        for s in range(200):
            np.random.seed(s)
            idx, _, dic = batch_select_whole(R, A, neighbor=[1], cell_size=2, sampling=mode)
            hits[mode] += sum(dic[c]['g'] == [big] for c in dic)
    # intended probability: 1001 / (1001 + 599) ~= 0.63 per draw (400 draws)
    assert hits['aligned'] > 0.5 * 400, hits
    assert hits['upstream'] < 0.1 * 400, hits


def test_paper_probabilities_match_equations():
    # cell k0 counts for its nonzero genes, and those genes' totals across all cells
    x = np.array([1., 2., 6., 10.])
    tot = np.array([1., 20., 6., 100.])                 # gene 0 and 2 are exclusive to this cell
    a = np.quantile(x, 0.25)                            # = 1.75 -> gene 0 (x=1) filtered out
    prop = x / tot
    expected = np.where(x > a, prop, 0) / np.where(x > a, prop, 0).sum()
    assert np.allclose(feature_probs(x, tot, 'paper', is_rna=True), expected)
    assert feature_probs(x, tot, 'paper', is_rna=True)[0] == 0
    # ATAC: no filter
    assert np.allclose(feature_probs(x, tot, 'paper', is_rna=False), prop / prop.sum())
    # specificity: an exclusive gene outranks a ubiquitous gene with higher count
    assert feature_probs(x, tot, 'paper', is_rna=True)[2] > feature_probs(x, tot, 'paper', is_rna=True)[3]


def test_paper_quartile_fallback_when_all_counts_equal():
    x = np.ones(5); p = feature_probs(x, np.full(5, 3.), 'paper', is_rna=True)
    assert np.allclose(p, 0.2)


def test_weighted_modes_cover_every_cell_once_and_respect_neighbor():
    R, A = _toy(seed=3)
    for mode in ('aligned', 'paper'):
        np.random.seed(0)
        idx, nodes_id, dic = batch_select_whole(R, A, neighbor=[15], cell_size=30, sampling=mode)
        cells = np.concatenate([b['cell_index'] for b in idx])
        assert sorted(cells.tolist()) == list(range(R.shape[1]))
        assert all(len(v['g']) <= 15 and len(v['p']) <= 15 for v in dic.values())
        for c, v in dic.items():                       # sampled features must be expressed/accessible in that cell
            assert set(v['g']) <= set(R[:, c].indices) and set(v['p']) <= set(A[:, c].indices)


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print('ok', name)
