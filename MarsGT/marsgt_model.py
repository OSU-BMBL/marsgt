import torch
import torch.nn as nn
import torch.nn.functional as F
from .conv import *
from .utils import *
from .egrn import *


class GNN_from_raw(nn.Module):
    def __init__(self, in_dim, n_hid, num_types, num_relations, n_heads, n_layers, dropout=0.2, conv_name='hgt',
                 prev_norm=True, last_norm=True):
        super(GNN_from_raw, self).__init__()
        self.gcs = nn.ModuleList()
        self.num_types = num_types
        self.in_dim = in_dim
        self.n_hid = n_hid
        self.adapt_ws = nn.ModuleList()
        self.drop = nn.Dropout(dropout)
        self.embedding1 = nn.ModuleList()

        # Initialize MLP weight matrices
        for ti in range(num_types):
            self.embedding1.append(nn.Linear(in_dim[ti], 256))

        for t in range(num_types):
            self.adapt_ws.append(nn.Linear(256, n_hid))

        # Initialize graph convolution layers
        for l in range(n_layers - 1):
            self.gcs.append(
                GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=prev_norm))
        self.gcs.append(
            GeneralConv(conv_name, n_hid, n_hid, num_types, num_relations, n_heads, dropout, use_norm=last_norm))

    def encode(self, x, t_id):
        h1 = F.relu(self.embedding1[t_id](x))
        return h1

    def forward(self, node_feature, node_type, edge_index, edge_type):
        node_embedding = []
        for t_id in range(self.num_types):
            node_embedding += list(self.encode(node_feature[t_id], t_id))

        node_embedding = torch.stack(node_embedding)
        # Initialize result matrix
        res = torch.zeros(node_embedding.size(0), self.n_hid).to(node_feature[0].device)

        # Process each node type
        for t_id in range(self.num_types):
            idx = (node_type == int(t_id))
            if idx.sum() == 0:
                continue
            # Update result matrix
            res[idx] = torch.tanh(self.adapt_ws[t_id](node_embedding[idx]))

        # Apply dropout to the result matrix
        meta_xs = self.drop(res)
        del res

        # Iterate through graph convolution layers and update result matrix
        for gc in self.gcs:
            meta_xs = gc(meta_xs, node_type, edge_index, edge_type)

        return meta_xs


class Net(nn.Module):
    def __init__(self, dim_in, dim_out):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(dim_in, dim_out)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        return x


class NodeDimensionReduction(nn.Module):
    def __init__(self, RNA_matrix, ATAC_matrix, indices, ini_p1, n_hid, n_heads,
                 n_layers, labsm, lr, wd, device, num_types=3, num_relations=2, epochs=1):
        super(NodeDimensionReduction, self).__init__()
        self.RNA_matrix = RNA_matrix
        self.ATAC_matrix = ATAC_matrix
        self.indices = indices
        self.ini_p1 = ini_p1
        self.in_dim = [RNA_matrix.shape[0], RNA_matrix.shape[1], ATAC_matrix.shape[1]]
        self.n_hid = n_hid
        self.num_types = num_types
        self.num_relations = num_relations
        self.n_heads = n_heads
        self.n_layers = n_layers
        self.labsm = labsm
        self.lr = lr
        self.wd = wd
        self.device = device
        self.epochs = epochs

        self.LabSm = LabelSmoothing(self.labsm)

        self.gnn = GNN_from_raw(in_dim=self.in_dim,
                                n_hid=self.n_hid,
                                num_types=self.num_types,
                                num_relations=self.num_relations,
                                n_heads=self.n_heads,
                                n_layers=self.n_layers,
                                dropout=0.3).to(self.device)

        self.optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=0.5, patience=5,
                                                                    verbose=True)

    def train_model(self, n_batch):
        print('The training process for the NodeDimensionReduction model has started. Please wait.')
        for epoch in tqdm(range(self.epochs)):
            for batch_id in np.arange(n_batch):
                gene_index = self.indices[batch_id]['gene_index']
                cell_index = self.indices[batch_id]['cell_index']
                peak_index = self.indices[batch_id]['peak_index']
                gene_feature = self.RNA_matrix[list(gene_index),]
                cell_feature = self.RNA_matrix[:, list(cell_index)].T
                peak_feature = self.ATAC_matrix[list(peak_index),]
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)
                peak_feature = torch.tensor(np.array(peak_feature.todense()), dtype=torch.float32).to(self.device)

                node_feature = [cell_feature, gene_feature, peak_feature]
                gene_cell_sub = self.RNA_matrix[list(gene_index),][:, list(cell_index)]
                peak_cell_sub = self.ATAC_matrix[list(peak_index),][:, list(cell_index)]
                # gene_cell_edge_index = torch.LongTensor([np.nonzero(gene_cell_sub)[0]+gene_cell_sub.shape[1],np.nonzero(gene_cell_sub)[1]]).to(device)
                # peak_cell_edge_index = torch.LongTensor([np.nonzero(peak_cell_sub)[0]+gene_cell_sub.shape[0]+gene_cell_sub.shape[1],np.nonzero(peak_cell_sub)[1]]).to(device)
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(
                    np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)
                peak_cell_edge_index1 = list(
                    np.nonzero(peak_cell_sub)[0] + gene_cell_sub.shape[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(peak_cell_sub)[1])
                peak_cell_edge_index2 = list(np.nonzero(peak_cell_sub)[1]) + list(
                    np.nonzero(peak_cell_sub)[0] + gene_cell_sub.shape[0] + gene_cell_sub.shape[1])
                peak_cell_edge_index = torch.LongTensor([peak_cell_edge_index1, peak_cell_edge_index2]).to(self.device)

                edge_index = torch.cat((gene_cell_edge_index, peak_cell_edge_index), dim=1)
                node_type = torch.LongTensor(np.array(
                    list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))) + list(
                        np.ones(len(peak_index)) * 2))).to(self.device)
                # edge_type = torch.LongTensor(np.array(list(np.zeros(gene_cell_edge_index.shape[1]))+list(np.ones(peak_cell_edge_index.shape[1]) ))).to(device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(
                    np.ones(np.nonzero(gene_cell_sub)[1].shape[0])) + list(
                    2 * np.ones(np.nonzero(peak_cell_sub)[0].shape[0])) + list(
                    3 * np.ones(np.nonzero(peak_cell_sub)[1].shape[0])))).to(self.device)
                l = torch.LongTensor(np.array(self.ini_p1)[[cell_index]]).to(self.device)
                node_rep = self.gnn.forward(node_feature, node_type,
                                            edge_index,
                                            edge_type).to(self.device)
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]
                peak_emb = node_rep[node_type == 2]

                decoder1 = torch.mm(gene_emb, cell_emb.t())
                decoder2 = torch.mm(peak_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)
                peak_cell_sub = torch.tensor(np.array(peak_cell_sub.todense()), dtype=torch.float32).to(self.device)

                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')

                logp_x2 = F.log_softmax(decoder2, dim=-1)
                p_y2 = F.softmax(peak_cell_sub, dim=-1)

                loss_kl2 = F.kl_div(logp_x2, p_y2, reduction='mean')
                loss_kl = loss_kl1 + loss_kl2
                loss_cluster = self.LabSm(cell_emb, l)
                lll = 0
                g = [int(i) for i in l]
                for i in set([int(k) for k in l]):
                    h = cell_emb[[True if i == j else False for j in g]]
                    ll = F.cosine_similarity(h[list(range(h.shape[0])) * h.shape[0],],
                                             h[[v for v in range(h.shape[0]) for i in range(h.shape[0])]]).mean()
                    lll = ll + lll
                loss = loss_cluster - lll
                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()
        print('The training for the NodeDimensionReduction model has been completed.')
        return self.gnn, cell_emb, gene_emb, peak_emb, h


class MarsGT(nn.Module):
    def __init__(self, gnn, h, labsm, n_hid, n_batch, device, lr, wd, num_epochs=1):
        super(MarsGT, self).__init__()
        self.lr = lr
        self.wd = wd
        self.gnn = gnn
        self.h = h
        self.n_hid = n_hid
        self.n_batch = n_batch
        self.device = device
        self.num_epochs = num_epochs
        self.net = Net(2 * self.n_hid, self.n_hid).to(self.device)
        self.gnn_optimizer = torch.optim.AdamW(self.gnn.parameters(), lr=self.lr, weight_decay=self.wd)
        self.net_optimizer = torch.optim.AdamW(self.net.parameters(), lr=1e-2)
        self.labsm = labsm
        self.LabSm = LabelSmoothing(self.labsm)

    def forward(self, indices, RNA_matrix, ATAC_matrix, Gene_Peak, ini_p1):
        cluster_l = list()
        cluster_kl_l = list()
        sim_l = list()
        nmi_l = list()
        ini_p1 = np.array(ini_p1)

        for epoch in range(self.num_epochs):
            for batch_id in tqdm(np.arange(self.n_batch)):
                gene_index = indices[batch_id]['gene_index']
                cell_index = indices[batch_id]['cell_index']
                peak_index = indices[batch_id]['peak_index']
                gene_feature = RNA_matrix[list(gene_index),]
                cell_feature = RNA_matrix[:, list(cell_index)].T
                peak_feature = ATAC_matrix[list(peak_index),]
                gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(self.device)
                cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(self.device)
                peak_feature = torch.tensor(np.array(peak_feature.todense()), dtype=torch.float32).to(self.device)
                node_feature = [cell_feature, gene_feature, peak_feature]
                gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
                peak_cell_sub = ATAC_matrix[list(peak_index),][:, list(cell_index)]
                gene_cell_edge_index1 = list(np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(gene_cell_sub)[1])
                gene_cell_edge_index2 = list(np.nonzero(gene_cell_sub)[1]) + list(
                    np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1])
                gene_cell_edge_index = torch.LongTensor([gene_cell_edge_index1, gene_cell_edge_index2]).to(self.device)
                peak_cell_edge_index1 = list(
                    np.nonzero(peak_cell_sub)[0] + gene_cell_sub.shape[0] + gene_cell_sub.shape[1]) + list(
                    np.nonzero(peak_cell_sub)[1])
                peak_cell_edge_index2 = list(np.nonzero(peak_cell_sub)[1]) + list(
                    np.nonzero(peak_cell_sub)[0] + gene_cell_sub.shape[0] + gene_cell_sub.shape[1])
                peak_cell_edge_index = torch.LongTensor([peak_cell_edge_index1, peak_cell_edge_index2]).to(self.device)
                # gene_cell_edge_index = torch.LongTensor([np.nonzero(gene_cell_sub)[0]+gene_cell_sub.shape[1],np.nonzero(gene_cell_sub)[1]]).to(device)
                # peak_cell_edge_index = torch.LongTensor([np.nonzero(peak_cell_sub)[0]+gene_cell_sub.shape[0]+gene_cell_sub.shape[1],np.nonzero(peak_cell_sub)[1]]).to(device)
                edge_index = torch.cat((gene_cell_edge_index, peak_cell_edge_index), dim=1)
                node_type = torch.LongTensor(np.array(
                    list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))) + list(
                        np.ones(len(peak_index)) * 2))).to(self.device)
                edge_type = torch.LongTensor(np.array(list(np.zeros(np.nonzero(gene_cell_sub)[0].shape[0])) + list(
                    np.ones(np.nonzero(gene_cell_sub)[1].shape[0])) + list(
                    2 * np.ones(np.nonzero(peak_cell_sub)[0].shape[0])) + list(
                    3 * np.ones(np.nonzero(peak_cell_sub)[1].shape[0])))).to(self.device)

                # edge_type = torch.LongTensor(np.array(list(np.zeros(gene_cell_edge_index.shape[1]))+list(np.ones(peak_cell_edge_index.shape[1]) ))).to(device)
                l = torch.LongTensor(np.array(ini_p1)[[cell_index]]).to(self.device)
                # l2 = torch.LongTensor(label[[cell_index]])

                # l = torch.LongTensor(ini_p1)[[cell_index]].to(device)

                node_rep = self.gnn.forward(node_feature, node_type,
                                            edge_index,
                                            edge_type).to(self.device)
                cell_emb = node_rep[node_type == 0]
                gene_emb = node_rep[node_type == 1]
                peak_emb = node_rep[node_type == 2]

                decoder1 = torch.mm(gene_emb, cell_emb.t())
                decoder2 = torch.mm(peak_emb, cell_emb.t())
                gene_cell_sub = torch.tensor(np.array(gene_cell_sub.todense()), dtype=torch.float32).to(self.device)
                peak_cell_sub = torch.tensor(np.array(peak_cell_sub.todense()), dtype=torch.float32).to(self.device)

                logp_x1 = F.log_softmax(decoder1, dim=-1)
                p_y1 = F.softmax(gene_cell_sub, dim=-1)

                loss_kl1 = F.kl_div(logp_x1, p_y1, reduction='mean')

                logp_x2 = F.log_softmax(decoder2, dim=-1)
                p_y2 = F.softmax(peak_cell_sub, dim=-1)

                loss_kl2 = F.kl_div(logp_x2, p_y2, reduction='mean')

                loss_kl = loss_kl1 + loss_kl2

                lll2 = 0
                g = [int(i) for i in l]
                for i in set([int(k) for k in l]):
                    ll2 = F.cosine_similarity(self.h[list(range(self.h.shape[0])) * self.h.shape[0],], self.h[
                        [v for v in range(self.h.shape[0]) for i in range(self.h.shape[0])]]).mean()
                    lll2 = ll2 + lll2

                loss_cluster = self.LabSm(cell_emb, l)

                m = range(peak_emb.shape[0])
                gene_emb_enh = gene_emb[list(range(gene_emb.shape[0])) * peak_emb.shape[0]]
                peak_emb_enh = peak_emb[[v for v in m for i in range(gene_emb.shape[0])]]

                net_input = torch.cat((gene_emb_enh, peak_emb_enh), 1)
                gene_peak_cluster_pre = self.net(net_input)
                m = range(peak_feature.shape[0])
                # gene_feature_ori = gene_feature[list(range (gene_feature.shape[0]))*peak_emb.shape[0]]
                # peak_feature_ori = peak_feature[[v for v in m for i in range(gene_feature.shape[0])]]
                peak_feature_ori = peak_Sparse(peak_feature[:,cell_index],gene_feature[:,cell_index],self.device)
                gene_feature_ori = gene_Sparse(peak_feature[:,cell_index],gene_feature[:,cell_index],self.device)
                # g1p1 g2p1 g3p1
                gene_peak_ori = gene_feature_ori.mul(peak_feature_ori)
                gene_peak_1 = Gene_Peak[list(gene_index),][:, list(peak_index)].reshape(
                    len(gene_index) * len(peak_index), 1)
                row_ind = torch.Tensor(list(gene_peak_1.tocoo().row) * len(cell_index))
                col_ind = torch.Tensor([v for v in range(len(cell_index)) for i in range(len(gene_peak_1.tocoo().row))])
                data = torch.Tensor(list(gene_peak_1.tocoo().data) * len(cell_index))
                a = torch.sparse.FloatTensor(torch.vstack((row_ind, col_ind)).long(), data,
                                             (gene_peak_1.shape[0], len(cell_index))).to(self.device)
                gene_peak_cell = gene_peak_ori.mul(a)
                # the original in cell level (peak * gene * peak_gene)
                gene_peak_cell = gene_peak_ori.mul(a)
                # the original in cell cluster level (peak * gene * peak_gene)
                gene_peak_cell_cluster = torch.mm(gene_peak_ori.mul(a), cell_emb)

                logp_x3 = F.log_softmax(gene_peak_cluster_pre, dim=-1)
                p_y3 = F.softmax(gene_peak_cell_cluster, dim=-1)

                loss_net = F.kl_div(logp_x3, p_y3, reduction='mean')

                loss = loss_net + loss_cluster + loss_kl  # - lll2
                # loss = loss_cluster  + loss_kl - lll + loss_net  #+ loss_S_R#+ loss_net
                self.gnn_optimizer.zero_grad()
                self.net_optimizer.zero_grad()
                loss.backward()
                self.gnn_optimizer.step()
                self.net_optimizer.step()
        return self.gnn

    def train_model(self, indices, RNA_matrix, ATAC_matrix, Gene_Peak, ini_p1):
        self.train()
        print('The training process for the MarsGT model has started. Please wait.')
        Mars_gnn = self.forward(indices, RNA_matrix, ATAC_matrix, Gene_Peak, ini_p1)
        print('The training for the MarsGT model has been completed.')
        return Mars_gnn


def MarsGT_pred(RNA_matrix, ATAC_matrix, RP_matrix, egrn, MarsGT_gnn, indices, nodes_id, cell_size, device, gene_names,
                peak_names):
    n_batch = math.ceil(nodes_id.shape[0] / cell_size)
    embedding = []
    l_pre = []
    MarsGT_result = {}
    with torch.no_grad():
        for batch_id in tqdm(range(n_batch)):
            gene_index = indices[batch_id]['gene_index']
            cell_index = indices[batch_id]['cell_index']
            peak_index = indices[batch_id]['peak_index']
            gene_feature = RNA_matrix[list(gene_index),]
            cell_feature = RNA_matrix[:, list(cell_index)].T
            peak_feature = ATAC_matrix[list(peak_index),]
            gene_feature = torch.tensor(np.array(gene_feature.todense()), dtype=torch.float32).to(device)
            cell_feature = torch.tensor(np.array(cell_feature.todense()), dtype=torch.float32).to(device)
            peak_feature = torch.tensor(np.array(peak_feature.todense()), dtype=torch.float32).to(device)
            node_feature = [cell_feature, gene_feature, peak_feature]
            gene_cell_sub = RNA_matrix[list(gene_index),][:, list(cell_index)]
            peak_cell_sub = ATAC_matrix[list(peak_index),][:, list(cell_index)]
            gene_cell_edge_index = torch.LongTensor(
                [np.nonzero(gene_cell_sub)[0] + gene_cell_sub.shape[1], np.nonzero(gene_cell_sub)[1]]).to(device)
            peak_cell_edge_index = torch.LongTensor(
                [np.nonzero(peak_cell_sub)[0] + gene_cell_sub.shape[0] + gene_cell_sub.shape[1],
                 np.nonzero(peak_cell_sub)[1]]).to(device)
            edge_index = torch.cat((gene_cell_edge_index, peak_cell_edge_index), dim=1)
            node_type = torch.LongTensor(np.array(
                list(np.zeros(len(cell_index))) + list(np.ones(len(gene_index))) + list(
                    np.ones(len(peak_index)) * 2))).to(device)
            edge_type = torch.LongTensor(np.array(
                list(np.zeros(gene_cell_edge_index.shape[1])) + list(np.ones(peak_cell_edge_index.shape[1])))).to(
                device)
            node_rep = MarsGT_gnn.forward(node_feature, node_type,
                                          edge_index,
                                          edge_type).to(device)
            cell_emb = node_rep[node_type == 0]
            gene_emb = node_rep[node_type == 1]
            peak_emb = node_rep[node_type == 2]

            # If the device is CUDA, copy the tensor to CPU memory
            if device.type == "cuda":
                cell_emb = cell_emb.cpu()
            # It is now safe to convert the tensor to a NumPy array
            embedding.append(cell_emb.detach().numpy())

            cell_pre = list(cell_emb.argmax(dim=1).detach().numpy())
            l_pre.extend(cell_pre)

    cell_embedding = np.vstack(embedding)
    cell_clu = np.array(l_pre)

    if egrn:
        final_egrn_df = egrn_calculate(cell_clu, nodes_id, RNA_matrix, ATAC_matrix, RP_matrix, gene_names, peak_names)
        MarsGT_result = {'pred_label': cell_clu, 'cell_embedding': cell_embedding, 'egrn': final_egrn_df}
        return MarsGT_result
    else:
        MarsGT_result = {'pred_label': cell_clu, 'cell_embedding': cell_embedding}
        return MarsGT_result