"""Opt-in options that bring the code closer to the paper (Wang et al. 2024, doi:10.1038/s41467-023-44570-8).

- num_relations=4 (paper: "four graph relations to pass the message"): peak<->cell edges now pass messages,
  so ATAC values affect stage 1 and the peak-source projections are trained.
- MarsGT_pred(pred_graph='training'): prediction uses exactly the graph built in training.
- MarsGT(cosine_loss=True): stage-2 loss includes the paper's within-cluster cosine term (eq. 12), on the current batch.

Run:  python -m pytest tests/   (or: python tests/test_paper_options.py)
"""
import numpy as np
import torch

import MarsGT.marsgt_model as mm
from test_training_controls import _data

DEV = torch.device('cpu')
KW = dict(n_hid=32, n_heads=4, n_layers=2, labsm=0.1, lr=5e-4, device=DEV)


def _stage1(D, num_relations, ATAC=None, wd=0.1, epochs=1):
    R, A, GP, ini, idx, nodes, names = D
    torch.manual_seed(0); np.random.seed(0)
    nd = mm.NodeDimensionReduction(R, A if ATAC is None else ATAC, idx, ini, wd=wd, epochs=epochs,
                                   num_relations=num_relations, **KW)
    gnn, _, _, _, h = nd.train_model(n_batch=len(idx))
    return nd, gnn, h


def test_four_relations_use_atac():
    torch.set_num_threads(1)
    D = _data(); A = D[1]
    A2 = A.copy(); A2.data = np.random.default_rng(9).uniform(1, 50, A2.data.size).astype(np.float32)
    shared = lambda g: {n: p.detach().clone() for n, p in g.named_parameters() if '.2.' not in n and not n.endswith('.2')}
    for nrel, expect_same in ((2, True), (4, False)):
        _, g1, _ = _stage1(D, nrel); _, g2, _ = _stage1(D, nrel, ATAC=A2)
        a, b = shared(g1), shared(g2)
        same = all(torch.equal(a[n], b[n]) for n in a)
        assert same == expect_same, f'num_relations={nrel}: ATAC values {"should not" if expect_same else "should"} affect stage 1'
    # peak-source key/value projections are trained only with 4 relations (weight decay off to isolate gradients)
    for nrel, expect_trained in ((2, False), (4, True)):
        torch.manual_seed(0); np.random.seed(0)
        R, A, GP, ini, idx, nodes, names = D
        nd = mm.NodeDimensionReduction(R, A, idx, ini, wd=0.0, epochs=1, num_relations=nrel, **KW)
        w0 = {n: p.detach().clone() for n, p in nd.gnn.named_parameters() if 'k_linears.2' in n or 'v_linears.2' in n}
        nd.train_model(n_batch=len(idx))
        w1 = {n: p.detach() for n, p in nd.gnn.named_parameters() if n in w0}
        trained = any(not torch.equal(w0[n], w1[n]) for n in w0)
        assert trained == expect_trained, (nrel, trained)


def test_prediction_graph_matches_training():
    torch.set_num_threads(1)
    D = _data(); R, A, GP, ini, idx, nodes, names = D
    nd, gnn, h = _stage1(D, 4)
    m = mm.MarsGT(gnn=gnn, h=h, labsm=0.1, n_hid=32, n_batch=len(idx), device=DEV, lr=5e-4, wd=0.1)
    seen = []
    orig = gnn.forward
    def spy(node_feature, node_type, edge_index, edge_type):
        seen.append((node_type.clone(), edge_index.clone(), edge_type.clone()))
        return orig(node_feature, node_type, edge_index, edge_type)
    gnn.forward = spy
    m.train_model(indices=idx, RNA_matrix=R, ATAC_matrix=A, Gene_Peak=GP, ini_p1=ini)
    train_b0 = seen[0]; seen.clear()
    kw = dict(egrn=False, MarsGT_gnn=gnn, indices=idx, nodes_id=nodes, cell_size=30, device=DEV,
              gene_names=names[0], peak_names=names[1])
    mm.MarsGT_pred(R, A, GP, pred_graph='training', **kw); pred_b0 = seen[0]; seen.clear()
    mm.MarsGT_pred(R, A, GP, **kw); released_b0 = seen[0]
    gnn.forward = orig
    assert all(torch.equal(a, b) for a, b in zip(train_b0, pred_b0)), "pred_graph='training' differs from the training graph"
    assert set(train_b0[2].unique().tolist()) == {0, 1, 2, 3}
    assert set(released_b0[2].unique().tolist()) == {0, 1}, 'released prediction graph should have only relations 0 and 1'
    assert released_b0[1].shape[1] < train_b0[1].shape[1], 'released prediction graph should have fewer edges'


def test_cosine_loss():
    torch.set_num_threads(1)
    D = _data(); R, A, GP, ini, idx, nodes, names = D
    out = {}
    for cos in (False, True):
        nd, gnn, h = _stage1(D, 2)
        m = mm.MarsGT(gnn=gnn, h=h, labsm=0.1, n_hid=32, n_batch=len(idx), device=DEV, lr=5e-4, wd=0.1, cosine_loss=cos)
        g = m.train_model(indices=idx, RNA_matrix=R, ATAC_matrix=A, Gene_Peak=GP, ini_p1=ini)
        out[cos] = ([p.detach().clone() for p in g.parameters()], m.history[-1])
    assert 'cos_sim' in out[True][1] and 'cos_sim' not in out[False][1]
    assert np.isfinite(out[True][1]['cos_sim']) and out[True][1]['cos_sim'] > 0
    assert not all(torch.equal(a, b) for a, b in zip(out[False][0], out[True][0])), 'cosine term had no effect'


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print('ok', name)
