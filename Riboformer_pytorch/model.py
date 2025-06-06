# main model structure - PyTorch version

import torch
import torch.nn as nn
import torch.nn.functional as F

from modules import ConvTower, TransformerBlock, TokenAndPositionEmbedding


class Riboformer(nn.Module):

    def __init__(self, configs):
        super().__init__()
        
        self.wsize = configs.wsize
        self.embed_dim = configs.embed_dim
        
        self.embedding_layer = TokenAndPositionEmbedding(configs.wsize, configs.vocab_size, configs.embed_dim)
        
        self.conv_tower1 = ConvTower('2D', [32,32,32,32,32], 5, activation=configs.activation)
        self.conv_tower2 = ConvTower('1D', [32,32,32,32,8], 9, activation=configs.activation)
        
        self.transformer_block1 = TransformerBlock(configs.embed_dim, configs.num_heads, configs.mlp_dim)
        self.transformer_block2 = TransformerBlock(configs.embed_dim, configs.num_heads, configs.mlp_dim)
        
        # Linear layer to project conv2 output to embed_dim for transformer compatibility
        self.conv2_to_embed = nn.Linear(8, configs.embed_dim)  # 8 is the last filter size in conv_tower2
        
        # Calculate the flattened size after conv operations
        # This will need to be adjusted based on actual tensor sizes
        self.head1 = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(configs.dropout_rate),
            nn.Linear(self._get_conv1_output_size(), 32),
            self._get_activation(configs.activation)
        )
        
        self.head2 = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(configs.dropout_rate),
            nn.Linear(self._get_conv2_output_size(), 32),
            self._get_activation(configs.activation)
        )
        
        self.final_dense = nn.Linear(32, 1)
        if configs.activation and configs.activation != 'linear':
            self.final_activation = self._get_activation(configs.activation)
        else:
            self.final_activation = None
    
    def _get_activation(self, activation):
        if activation == 'relu':
            return nn.ReLU()
        elif activation == 'gelu':
            return nn.GELU()
        elif activation == 'sigmoid':
            return nn.Sigmoid()
        elif activation == 'tanh':
            return nn.Tanh()
        else:
            return nn.Identity()
    
    def _get_conv1_output_size(self):
        # Calculate output size dynamically
        # After conv operations, the spatial dimensions are preserved due to 'same' padding
        # After transformer: (batch, wsize, embed_dim), then flattened
        return self.wsize * self.embed_dim
    
    def _get_conv2_output_size(self):
        # Calculate output size dynamically  
        # After conv_tower2, then transformer: (batch, wsize, embed_dim), then flattened
        # Note: The transformer output maintains embed_dim, not the conv filter size
        return self.wsize * self.embed_dim
    
    def forward(self, inputs):
        seq, exp = inputs

        # mRNA sequence branch
        x = self.embedding_layer(seq)  # (batch, wsize, embed_dim)
        x = x.unsqueeze(-1)  # Add channel dimension: (batch, wsize, embed_dim, 1)
        x = x.permute(0, 3, 1, 2)  # PyTorch conv2d expects (batch, channels, height, width)
        
        x = self.conv_tower1(x)  # (batch, 32, wsize, embed_dim)
        x = torch.mean(x, dim=1)  # Reduce channel dimension: (batch, wsize, embed_dim)
        
        x, weights1 = self.transformer_block1(x, training=self.training)
        x = self.head1(x)
        
        # control experiment branch  
        y = exp.unsqueeze(1)  # Add channel dimension: (batch, 1, wsize)
        
        y = self.conv_tower2(y)  # (batch, 8, wsize)
        y = y.permute(0, 2, 1)  # Convert to (batch, wsize, 8) for projection
        
        # Project to embed_dim for transformer compatibility
        y = self.conv2_to_embed(y)  # (batch, wsize, embed_dim)
        
        y, weights2 = self.transformer_block2(y, training=self.training)
        y = self.head2(y)

        # Element-wise multiplication and final prediction
        outputs = self.final_dense(x * y)
        
        if self.final_activation:
            outputs = self.final_activation(outputs)
        
        return outputs
