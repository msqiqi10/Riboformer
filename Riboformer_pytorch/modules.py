# modules for the model - PyTorch version

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, Dropout, Conv1d, Conv2d, BatchNorm1d, BatchNorm2d
from typing import Optional
import math


class ConvTower(nn.Module):
    def __init__(self, 
                 func: str, 
                 filters: list[int],
                 kernel_size: int,
                 padding: str = "same", 
                 activation: Optional[str] = None):
        super(ConvTower, self).__init__()
        
        self.func = Conv1d if func == '1D' else Conv2d
        self.batch_norm = BatchNorm1d if func == '1D' else BatchNorm2d
        self.activation = activation
        
        # Calculate padding for 'same' padding
        if padding == "same":
            self.padding = kernel_size // 2
        else:
            self.padding = 0
        
        layers = []
        for i, f in enumerate(filters):
            if func == '1D':
                conv_layer = Conv1d(
                    in_channels=1 if i == 0 else filters[i-1],
                    out_channels=f,
                    kernel_size=kernel_size,
                    padding=self.padding
                )
                bn_layer = BatchNorm1d(f)
            else:  # 2D
                conv_layer = Conv2d(
                    in_channels=1 if i == 0 else filters[i-1],
                    out_channels=f,
                    kernel_size=kernel_size,
                    padding=self.padding
                )
                bn_layer = BatchNorm2d(f)
            
            layers.append(conv_layer)
            layers.append(bn_layer)
            
            if activation:
                if activation == 'relu':
                    layers.append(nn.ReLU())
                elif activation == 'gelu':
                    layers.append(nn.GELU())
                    
        self.nn = Sequential(*layers)
            
    def forward(self, inputs):
        return self.nn(inputs)


class TransformerBlock(nn.Module):
    def __init__(self, 
                 embed_dim, 
                 num_heads, 
                 ff_dim, 
                 dropout_rate=0.1,
                 **kwargs):
        super(TransformerBlock, self).__init__()
        
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.ff_dim = ff_dim
        
        self.mha = nn.MultiheadAttention(
            embed_dim=embed_dim,
            num_heads=num_heads,
            dropout=dropout_rate,
            batch_first=True
        )
        
        self.mlp = Sequential(
            Linear(embed_dim, ff_dim),
            nn.ReLU(),
            Linear(ff_dim, embed_dim)
        )
        
        self.layernorm1 = nn.LayerNorm(embed_dim, eps=1e-6)
        self.layernorm2 = nn.LayerNorm(embed_dim, eps=1e-6)
        self.dropout = nn.Dropout(dropout_rate)

    def forward(self, inputs, training=True):
        # multi-head attention layer
        residual = inputs
        inputs = self.layernorm1(inputs)
        mha_output, weights = self.mha(inputs, inputs, inputs)
        mha_output = self.dropout(mha_output) if training else mha_output
        mha_output = residual + mha_output
        
        # multi-layer perceptron layer
        residual = mha_output
        output = self.layernorm2(mha_output)
        output = self.mlp(output)
        output = self.dropout(output) if training else output
        output = residual + output
        
        return output, weights


class TokenAndPositionEmbedding(nn.Module):
    def __init__(self, 
                 maxlen, 
                 vocab_size, 
                 embed_dim,
                 **kwargs):
        super(TokenAndPositionEmbedding, self).__init__()
        
        self.maxlen = maxlen
        self.vocab_size = vocab_size
        self.embed_dim = embed_dim
        
        self.token_emb = nn.Embedding(vocab_size, embed_dim)
        self.pos_emb = nn.Embedding(maxlen, embed_dim)
        
        # Initialize position indices
        self.register_buffer('position_ids', torch.arange(maxlen).expand((1, -1)))
        self.position_ids: torch.Tensor

    def forward(self, inputs):
        seq_len = inputs.size(1)
        position_ids = self.position_ids[:, :seq_len]
        
        token_embeddings = self.token_emb(inputs)
        position_embeddings = self.pos_emb(position_ids)
        
        return token_embeddings + position_embeddings
