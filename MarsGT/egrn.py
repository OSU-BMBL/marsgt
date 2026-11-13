import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import scipy.sparse as sp
import anndata as ad
from collections import Counter
import pandas as pd


def peak_Sparse(peak_feature, gene_feature, device):
    """Vectorized. For every nonzero (p, c, v) of peak_feature and every gene g, emit (p*G + g, c, v).
    Shape (G*P, C); identical to the original loop version."""
    G = gene_feature.shape[0]
    peak_row, peak_col = torch.nonzero(peak_feature, as_tuple=True)
    peak_data = peak_feature[peak_row, peak_col]
    g = torch.arange(G, device=peak_feature.device)
    rows = (peak_row.unsqueeze(1) * G + g.unsqueeze(0)).reshape(-1)
    cols = peak_col.unsqueeze(1).expand(-1, G).reshape(-1)
    vals = peak_data.unsqueeze(1).expand(-1, G).reshape(-1)
    return torch.sparse_coo_tensor(torch.stack((rows, cols)), vals,
                                   (G * peak_feature.shape[0], peak_feature.shape[1])).to(device)


def gene_Sparse(peak_feature, gene_feature, device):
    """Vectorized. For every nonzero (g, c, v) of gene_feature and every peak p, emit (p*G + g, c, v).
    Shape (G*P, C); identical to the original loop version."""
    G, P = gene_feature.shape[0], peak_feature.shape[0]
    gene_row, gene_col = torch.nonzero(gene_feature, as_tuple=True)
    gene_data = gene_feature[gene_row, gene_col]
    p = torch.arange(P, device=gene_feature.device)
    rows = (p.unsqueeze(1) * G + gene_row.unsqueeze(0)).reshape(-1)
    cols = gene_col.repeat(P)
    vals = gene_data.repeat(P)
    return torch.sparse_coo_tensor(torch.stack((rows, cols)), vals,
                                   (G * P, gene_feature.shape[1])).to(device)


def gene_Peak_Score(pred_label,nodes_id,RNA_matrix,ATAC_matrix,RP_matrix,gene_names):
    y = ad.AnnData(RNA_matrix.transpose().astype('float32'))
    y.var_names=gene_names[0]
    y.obs['pred'] = np.array(pred_label,dtype='int64')
    y.obs['pred'] = y.obs['pred'].astype("category")
    # Store the prediction results as a sparse matrix, where rows represent cells and columns represent cluster labels. A value of 1 at the corresponding position indicates that the cell belongs to that cluster.
    cell_emb_binary = sp.coo_matrix((np.ones(len(y.obs['pred'])),(np.array(range(len(y.obs['pred']))),list(y.obs['pred']))))
    cell_emb_binary.todense()

    adata_atac_all = ATAC_matrix[:,nodes_id]
    adata_rna_all = RNA_matrix[:,nodes_id]

    ATAC_matrix_ct = adata_atac_all*(cell_emb_binary)
    RNA_matrix_ct = adata_rna_all*(cell_emb_binary)
    RNA_matrix_ct = RNA_matrix_ct/np.sum(cell_emb_binary,axis=0)
    ATAC_matrix_ct = ATAC_matrix_ct/np.sum(cell_emb_binary,axis=0)

    rna_matrix = RNA_matrix_ct
    m = range(rna_matrix.shape[0])
    gene_peak = RP_matrix
    gp = gene_peak.reshape(gene_peak.shape[1]*rna_matrix.shape[0],1).todense()
    gene_emb_enh = ATAC_matrix_ct[list(range (ATAC_matrix_ct.shape[0]))*rna_matrix.shape[0]]
    peak_emb_enh = RNA_matrix_ct[[v for v in m for i in range(ATAC_matrix_ct.shape[0])]]
    egrn = np.multiply(gene_emb_enh, peak_emb_enh)
    egrn = np.multiply(egrn,gp)
    return egrn


def small_gene_Peak_Score(pred_label,RNA_matrix,ATAC_matrix,RP_matrix,gene_names,peak_names,choice_cluser):  
    y = ad.AnnData(RNA_matrix.transpose().astype('float32'))
    y.var_names=gene_names[0]
    y.obs['pred'] = np.array(pred_label,dtype='int64')
    cell_emb_binary = sp.coo_matrix((np.ones(len(y.obs['pred'])),(np.array(range(len(y.obs['pred']))),list(y.obs['pred']))))
    cell_emb_binary.todense()

    adata_atac_all = ATAC_matrix
    adata_rna_all = RNA_matrix

    ATAC_matrix_ct = adata_atac_all*(cell_emb_binary)
    RNA_matrix_ct = adata_rna_all*(cell_emb_binary)
    
    RNA_matrix_ct = RNA_matrix_ct/np.sum(cell_emb_binary,axis=0)
    ATAC_matrix_ct = ATAC_matrix_ct/np.sum(cell_emb_binary,axis=0)

    rna_matrix = RNA_matrix_ct
    m = range(rna_matrix.shape[0])
    gene_peak = RP_matrix
    gp = gene_peak.reshape(gene_peak.shape[1]*rna_matrix.shape[0],1).todense()
    gene_emb_enh = ATAC_matrix_ct[list(range (ATAC_matrix_ct.shape[0]))*rna_matrix.shape[0]]
    peak_emb_enh = RNA_matrix_ct[[v for v in m for i in range(ATAC_matrix_ct.shape[0])]]
    egrn = np.multiply(gene_emb_enh, peak_emb_enh)
    egrn = np.multiply(egrn,gp)
    gn = np.array(gene_names[0])[[v for v in m for i in range(ATAC_matrix.shape[0])]]
    pn = np.array(peak_names[0])[list(range (ATAC_matrix.shape[0]))*RNA_matrix.shape[0]]
    data = np.array(np.squeeze(egrn[:,choice_cluser]))[0]
    gene_peak_score_Dict = {'gene': gn, 'peak': pn, 'score':data}
    gene_peak_score_df = pd.DataFrame(data = gene_peak_score_Dict)
    gene_peak_score = gene_peak_score_df[gene_peak_score_df['score'] > 0]
    gene_peak_score = gene_peak_score.drop_duplicates(['gene','peak'])
    gene_peak_score = gene_peak_score.sort_values(by="score",ascending=False)
    
    return gene_peak_score_df,gene_peak_score


def egrn_calculate(pred_label,nodes_id,RNA_matrix,ATAC_matrix,RP_matrix,gene_names,peak_names,threshold=0):
    print('We are currently performing calculations for EGRN. Please bear with us as this process will take approximately around 10 minutes.')
    
    try:
        egrn = gene_Peak_Score(pred_label,nodes_id,RNA_matrix,ATAC_matrix,RP_matrix,gene_names)
        label_num = len(Counter(pred_label))
        print('class_num:',label_num)
        gn = np.array(gene_names[0])[[v for v in range(RNA_matrix.shape[0]) for i in range(ATAC_matrix.shape[0])]]
        pn = np.array(peak_names[0])[list(range (ATAC_matrix.shape[0]))*RNA_matrix.shape[0]]
        egrn_df = pd.DataFrame(columns=['gene', 'peak', 'score', 'class'])
        print('We are currently conducting filtering and categorization operations.')
        for i in range(label_num):
            data = np.array(np.squeeze(egrn[:,i]))[0]
            gene_peak_score_Dict = {'gene': gn, 'peak': pn, 'score':data}
            gene_peak_score_df = pd.DataFrame(data = gene_peak_score_Dict)
            gene_peak_score = gene_peak_score_df[gene_peak_score_df['score'] >threshold]
            gene_peak_score = gene_peak_score.drop_duplicates(['gene','peak'])
            gene_peak_score = gene_peak_score.sort_values(by="score",ascending=False)
            gene_peak_score['class'] = i
            egrn_df = pd.concat([egrn_df, gene_peak_score])
    
    except:
        try:
            print('An error occurred with gene_Peak_Score. We will try using small_gene_Peak_Score instead.')
            egrn_df = pd.DataFrame(columns=['gene', 'peak', 'score', 'class'])
            Total_name = list(range(len(Counter(pred_label))))
            for i in range(len(Counter(pred_label))):
                Total_gene_peak_score_df,Total_gene_peak_score = small_gene_Peak_Score(pred_label,RNA_matrix,ATAC_matrix,RP_matrix,gene_names,peak_names,choice_cluser=i)
                Total_gene_peak_score_df = Total_gene_peak_score_df[Total_gene_peak_score_df['score'] > threshold]
                Total_gene_peak_score_df['class'] = Total_name[i]
                egrn_df = pd.concat([egrn_df, Total_gene_peak_score_df.sort_values(by="score",ascending=False)])
        except:
            print('You are encountering memory shortage issues, please follow the tutorial instructions to address it.')
            egrn_df = pd.DataFrame(columns=['gene', 'peak', 'score', 'class'])
            pass
            
    return egrn_df