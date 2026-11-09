"""Device-side feature rows (MarsGT.utils.dense_rows / DeviceRows) must give exactly the original values,
and switching them on must not change training.

Run:  python -m pytest tests/   (or: python tests/test_device_rows.py)
"""
import numpy as np
import scipy.sparse as sp
import torch

import MarsGT.utils as U


def _counts(r, c, density, lam, seed):
    rng = np.random.default_rng(seed)
    return sp.random(r, c, density=density, format='csc', random_state=rng,
                     data_rvs=lambda n: rng.poisson(lam, n) + 1).astype(np.float32)


def _original(M, idx, device):
    return torch.tensor(np.array(M[list(idx), ].todense()), dtype=torch.float32).to(device)


def test_rows_identical_to_original():
    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    rng = np.random.default_rng(0)
    for (r, c, dens) in [(50, 40, 0.0), (300, 1000, 0.05), (3000, 2000, 0.15), (20000, 500, 0.03)]:
        M = _counts(r, c, dens, 3, seed=r)
        M[min(3, r - 1)] = 0                          # an empty row
        M.eliminate_zeros()
        rows = U.DeviceRows(M, dev)
        for idx in [[0], list(rng.choice(r, min(r, 470), replace=False)), [1, 1, 2], list(range(r))[:5], []]:
            a, b = _original(M, idx, dev) if idx else torch.zeros((0, c)), rows.rows(idx)
            assert a.shape == b.shape and torch.equal(a.cpu(), b.cpu()), (r, c, dens, len(idx))
            if idx:   # same memory order (row- vs column-major); strides of size-1 dims are irrelevant
                assert (a.is_contiguous(), a.t().is_contiguous()) == (b.is_contiguous(), b.t().is_contiguous()), \
                    ('layout', a.stride(), b.stride())


def test_dense_rows_switch_and_cache():
    M = _counts(200, 300, 0.1, 2, seed=1); idx = [5, 7, 100]
    U.GPU_FEATURES = False; a = U.dense_rows(M, idx, 'cpu')
    U.GPU_FEATURES = True;  b = U.dense_rows(M, idx, 'cpu'); c = U.dense_rows(M, idx, 'cpu')
    assert torch.equal(a, b) and torch.equal(b, c) and len(U._ROW_CACHE) == 1
    M2 = M.copy(); M2[5, 0] = 99.0                    # a different (modified) matrix must not hit the stale cache
    assert U.dense_rows(M2, idx, 'cpu')[0, 0] == 99.0
    U.clear_feature_cache(); U.GPU_FEATURES = 'auto'


def test_training_unchanged():
    """Stages 1 + 2 + prediction with device rows on vs off: bitwise-identical weights and embeddings
    (CPU, 1 thread; device rows reproduce the original values AND memory layout, so the arithmetic is the same)."""
    import pandas as pd
    import MarsGT.marsgt_model as mm
    torch.set_num_threads(1)
    rng = np.random.default_rng(1); N, G, P = 400, 800, 4000
    R, A = _counts(G, N, 0.15, 4, 2), _counts(P, N, 0.03, 1, 3)
    GP = sp.random(G, P, density=0.002, format='csr', random_state=4).astype(np.float32); ini = rng.integers(0, 6, N)
    np.random.seed(0); idx, nodes_id, _ = U.batch_select_whole(R, A, neighbor=[15], cell_size=30)
    names = (pd.DataFrame([f'g{i}' for i in range(G)]), pd.DataFrame([f'p{i}' for i in range(P)]))

    def run(flag):
        U.GPU_FEATURES = flag; U.clear_feature_cache()
        torch.manual_seed(0); np.random.seed(0)
        nd = mm.NodeDimensionReduction(R, A, idx, ini, n_hid=48, n_heads=4, n_layers=2, labsm=0.1, lr=5e-4, wd=0.1,
                                       device=torch.device('cpu'), epochs=1)
        gnn, _, _, _, h = nd.train_model(n_batch=3)
        m = mm.MarsGT(gnn=gnn, h=h, labsm=0.1, n_hid=48, n_batch=3, device=torch.device('cpu'), lr=5e-4, wd=0.1)
        g = m.train_model(indices=idx, RNA_matrix=R, ATAC_matrix=A, Gene_Peak=GP, ini_p1=ini)
        res = mm.MarsGT_pred(R, A, GP, egrn=False, MarsGT_gnn=g, indices=idx[:3], nodes_id=nodes_id[:90], cell_size=30,
                             device=torch.device('cpu'), gene_names=names[0], peak_names=names[1])
        return [p.detach().clone() for p in g.parameters()], res['cell_embedding']

    w_off, e_off = run(False)
    w_on, e_on = run(True)
    U.GPU_FEATURES = 'auto'; U.clear_feature_cache()
    assert all(torch.equal(a, b) for a, b in zip(w_off, w_on))
    assert np.array_equal(e_off, e_on)


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print('ok', name)
