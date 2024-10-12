import numpy as np
import pandas as pd
import sys
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import encode
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as pyg_nn
import torch_geometric.utils as pyg_utils
import torch
import torch.optim as optim

from torch_geometric.data import DataLoader
import torch_geometric.transforms as T
from torch_geometric.data import InMemoryDataset
from torch_geometric.data import Data

input_fasta = sys.argv[1]
output_file = sys.argv[2]
reverse = sys.argv[3]
thread = int(sys.argv[4])

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

HIDDEN_DIM = 3
PNODE_DIM = 3
PNODE_NUM = 4096
FNODE_NUM = 64
NODE_NUM = 4160
BATCH_SIZE = 32
LR = 1e-4
GCN_HIDDEN_DIM = 128
CNN_HIDDEN_DIM = 64
FC_HIDDEN_DIM = 100
GCN_LAYER_NUM = 2
DROP_RATE = 0.2


class GNN_Model(nn.Module):
    def __init__(self):
        super(GNN_Model, self).__init__()
        self.gcn_dim = GCN_HIDDEN_DIM
        self.cnn_dim = CNN_HIDDEN_DIM
        self.fc_dim = FC_HIDDEN_DIM
        self.num_layers = GCN_LAYER_NUM
        self.dropout = DROP_RATE

        self.pnode_d = nn.Linear(PNODE_NUM * PNODE_DIM, PNODE_NUM * HIDDEN_DIM)
        self.fnode_d = nn.Linear(FNODE_NUM, FNODE_NUM * HIDDEN_DIM)

        self.convs_1 = nn.ModuleList()
        for l in range(self.num_layers):
            if l == 0:
                self.convs_1.append(pyg_nn.SAGEConv((HIDDEN_DIM, HIDDEN_DIM), self.gcn_dim))
            else:
                self.convs_1.append(pyg_nn.SAGEConv((self.gcn_dim, self.gcn_dim), self.gcn_dim))

        self.convs_2 = nn.ModuleList()
        for l in range(self.num_layers):
            if l == 0:
                self.convs_2.append(pyg_nn.SAGEConv((self.gcn_dim, HIDDEN_DIM), self.gcn_dim))
            else:
                self.convs_2.append(pyg_nn.SAGEConv((self.gcn_dim, self.gcn_dim), self.gcn_dim))

        self.lns = nn.ModuleList()
        for l in range(self.num_layers-1):
            self.lns.append(nn.LayerNorm(self.gcn_dim))

        self.conv1 = nn.Conv1d(in_channels=self.gcn_dim, out_channels=64, kernel_size=8)
        self.conv2 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=8)
        self.conv3 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=8)
        self.d1 = nn.Linear(4075 * 64, 100)
        self.d2 = nn.Linear(100, 2)


    def forward(self, data):
        x_f = data.x_src
        x_p = data.x_dst
        edge_index_forward = data.edge_index[:,::2]
        edge_index_backward = data.edge_index[[1, 0], :][:,1::2]

        # reserve primary nodes dim
        x_p = torch.reshape(x_p, (-1, PNODE_NUM * PNODE_DIM))
        x_p = self.pnode_d(x_p)
        x_p = torch.reshape(x_p, (-1, HIDDEN_DIM))

        # change feature nodes dim
        x_f = torch.reshape(x_f, (-1, FNODE_NUM))
        x_f = self.fnode_d(x_f)
        x_f = torch.reshape(x_f, (-1, HIDDEN_DIM))

        for i in range(self.num_layers):
            x_p = self.convs_1[i]((x_f, x_p), edge_index_forward)
            x_p = F.relu(x_p)
            x_p = F.dropout(x_p, p=self.dropout, training=self.training)
            x_f = self.convs_2[i]((x_p, x_f), edge_index_backward)
            x_f = F.relu(x_f)
            x_f = F.dropout(x_f, p=self.dropout, training=self.training)
            if not i == self.num_layers - 1:
                x_p = self.lns[i](x_p)
                x_f = self.lns[i](x_f)
        x = torch.reshape(x_p, (-1, self.gcn_dim, PNODE_NUM))
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.conv3(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)

        x = x.flatten(start_dim = 1)
        x = self.d1(x)
        x = F.relu(x)

        logits = self.d2(x)
        out = F.softmax(logits, dim=1)

        return out


class BipartiteData(Data):
    def __inc__(self, key, value):
        if key == 'edge_index':
            return torch.tensor([[self.x_src.size(0)], [self.x_dst.size(0)]])
        else:
            return super(BipartiteData, self).__inc__(key, value)


def make_edge():
    edge = []
    for i in range(4096):
        a = i // 64
        b = i % 64
        edge.append([a, i])
        edge.append([b, i])
    edge = np.array(edge).T

    return edge


def generate_model_input(fasta_file, K, thread):
    seq_list = []
    data_id = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_list.append(str(seq_record.seq))
        data_id.append(seq_record.id)
    
    pool = Pool(thread)
    partial_encode_seq = partial(encode.matrix_encoding, K=K)
    feature = np.array(pool.map(partial_encode_seq, seq_list))
    pool.close()
    pool.join()

    return data_id, feature


if __name__ == "__main__":
    data_id, input_data = generate_model_input(input_fasta, 3, thread)
    edge = make_edge()
    pnode_feature = input_data.reshape(-1, 3, 4096)
    pnode_feature = np.moveaxis(pnode_feature, 1, 2)

    zero_layer = input_data.reshape(-1, 3, 64, 64)[:, 0, :, :]
    fnode_feature = np.sum(zero_layer, axis=2).reshape(-1, 64, 1)
    
    data_list = []
    for i in range(pnode_feature.shape[0]):
        edge_index = torch.tensor(edge, dtype=torch.long)  # edge_index should be long type
        x_p = torch.tensor(pnode_feature[i, :, :], dtype=torch.float)
        x_f = torch.tensor(fnode_feature[i, :, :], dtype=torch.float)
        
        data = BipartiteData(x_src=x_f, x_dst=x_p, edge_index=edge_index, num_nodes=None)
        data_list.append(data)

    loader = DataLoader(data_list, batch_size=100, shuffle=False, follow_batch=['x_src', 'x_dst'])
    
    model = torch.load("/home/ruohawang2/12.Phage_assem/phage_scoring/GCNFrame/GCN_model_retrained.pt")
    model.eval()
    if reverse:
        f_out = open(output_file, "w")
        count = 0
        for data in loader:
            with torch.no_grad():
                data = data.to(device)
                pred = model(data)
                # pred = np.array(pred.argmax(dim=1).cpu())
                pred = np.array(pred[:, 1].cpu())
                for each in pred:
                    if count != 0:
                        f_out.write("\n")
                    f_out.write(str(data_id[count]) + "\t")
                    f_out.write(str(each) + "\n")
                    f_out.write(str(data_id[count]) + "'" + "\t")
                    f_out.write(str(each))
                    count += 1
        f_out.close()
    
    else:
        f_out = open(output_file, "w")
        count = 0
        for data in loader:
            with torch.no_grad():
                data = data.to(device)
                pred = model(data)
                # pred = np.array(pred.argmax(dim=1).cpu())
                pred = np.array(pred[:, 1].cpu())
                for each in pred:
                    if count != 0:
                        f_out.write("\n")
                    f_out.write(str(data_id[count]) + "\t")
                    f_out.write(str(each))
                    count += 1
        f_out.close()


