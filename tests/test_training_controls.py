"""Training controls used for early stopping must not change the math.

- one call with N epochs  ==  N calls with 1 epoch on the same objects (same optimizers)     [bitwise]
- an eval-mode prediction between epochs leaves training unchanged                           [bitwise]
- eval-mode predictions are repeatable; the model's train/eval mode is restored afterwards
- per-epoch loss history is recorded for both stages

CPU, single thread (multithreaded CPU and GPU kernels are not bitwise deterministic).
Run:  python -m pytest tests/   (or: python tests/test_training_controls.py)
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

import MarsGT.utils as U
import MarsGT.marsgt_model as mm

DEV = torch.device('cpu')


def _data(seed=1, N=300, G=600, P=3000):
    rng = np.random.default_rng(seed)
    def counts(r, c, d, lam, s):
        g = np.random.default_rng(s)
        return sp.random(r, c, density=d, format='csc', random_state=g, data_rvs=lambda n: g.poisson(lam, n) + 1).astype(np.float32)
    R, A = counts(G, N, 0.15, 4, seed + 1), counts(P, N, 0.03, 1, seed + 2)
    GP = sp.random(G, P, density=0.003, format='csr', random_state=seed + 3).astype(np.float32)
    ini = rng.integers(0, 5, N)
    np.random.seed(0)
    idx, nodes, _ = U.batch_select_whole(R, A, neighbor=[12], cell_size=30, sampling='paper')
    names = (pd.DataFrame([f'g{i}' for i in range(G)]), pd.DataFrame([f'p{i}' for i in range(P)]))
    return R, A, GP, ini, idx, nodes, names


def _train(D, e1, e2, stepwise=False, monitor=False):
    R, A, GP, ini, idx, nodes, names = D
    torch.manual_seed(0); np.random.seed(0)
    nd = mm.NodeDimensionReduction(R, A, idx, ini, n_hid=32, n_heads=4, n_layers=2, labsm=0.1, lr=5e-4, wd=0.1,
                                   device=DEV, epochs=1 if stepwise else e1)
    for _ in range(e1 if stepwise else 1):
        gnn, _, _, _, h = nd.train_model(n_batch=len(idx))
    m = mm.MarsGT(gnn=gnn, h=h, labsm=0.1, n_hid=32, n_batch=len(idx), device=DEV, lr=5e-4, wd=0.1,
                  num_epochs=1 if stepwise else e2)
    preds = []
    for _ in range(e2 if stepwise else 1):
        g = m.train_model(indices=idx, RNA_matrix=R, ATAC_matrix=A, Gene_Peak=GP, ini_p1=ini)
        if monitor:
            preds.append(mm.MarsGT_pred(R, A, GP, egrn=False, MarsGT_gnn=g, indices=idx, nodes_id=nodes, cell_size=30,
                                        device=DEV, gene_names=names[0], peak_names=names[1])['pred_label'])
    return [p.detach().clone() for p in g.parameters()], nd.history, m.history, preds, g


def test_controls():
    torch.set_num_threads(1)
    D = _data()
    w_ref, h1, h2, _, _ = _train(D, 2, 3)
    w_step, h1s, h2s, _, _ = _train(D, 2, 3, stepwise=True)
    w_mon, _, _, preds, g = _train(D, 2, 3, stepwise=True, monitor=True)
    assert all(torch.equal(a, b) for a, b in zip(w_ref, w_step)), 'stepwise epochs differ from one multi-epoch call'
    assert all(torch.equal(a, b) for a, b in zip(w_ref, w_mon)), 'eval-mode monitoring changed training'
    assert len(h1) == 2 and len(h2) == 3 and h1 == h1s and h2 == h2s
    assert all(np.isfinite(v) for e in h2 for v in e.values())
    assert g.training, 'train mode not restored after prediction'

    R, A, GP, ini, idx, nodes, names = D
    kw = dict(egrn=False, MarsGT_gnn=g, indices=idx, nodes_id=nodes, cell_size=30, device=DEV,
              gene_names=names[0], peak_names=names[1])
    a = mm.MarsGT_pred(R, A, GP, **kw)['pred_label']; b = mm.MarsGT_pred(R, A, GP, **kw)['pred_label']
    assert np.array_equal(a, b), 'eval-mode predictions not repeatable'
    torch.manual_seed(5); c = mm.MarsGT_pred(R, A, GP, eval_mode=False, **kw)['pred_label']
    torch.manual_seed(6); d = mm.MarsGT_pred(R, A, GP, eval_mode=False, **kw)['pred_label']
    assert not np.array_equal(c, d), 'eval_mode=False should reproduce the original stochastic (dropout) prediction'
    assert g.training


if __name__ == '__main__':
    test_controls(); print('ok test_controls')
