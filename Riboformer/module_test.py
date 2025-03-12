import unittest
import tensorflow as tf
import numpy as np
from modules import ConvTower, TransformerBlock, TokenAndPositionEmbedding

class TestModules(unittest.TestCase):
    def setUp(self):
        self.batch_size = 32
        self.seq_length = 100
        self.channels = 16
        self.embed_dim = 32
        self.num_heads = 4
        self.ff_dim = 64
        self.vocab_size = 1000

    def test_conv_tower_1d(self):
        # Test 1D ConvTower
        conv_tower = ConvTower('1D', [32, 64], kernel_size=3)
        inputs = tf.random.normal((self.batch_size, self.seq_length, self.channels))
        outputs = conv_tower(inputs)
        self.assertEqual(outputs.shape, (self.batch_size, self.seq_length, 64))

    def test_conv_tower_2d(self):
        # Test 2D ConvTower
        conv_tower = ConvTower('2D', [32, 64], kernel_size=3)
        inputs = tf.random.normal((self.batch_size, self.seq_length, self.seq_length, self.channels))
        outputs = conv_tower(inputs)
        self.assertEqual(outputs.shape, (self.batch_size, self.seq_length, self.seq_length, 64))

    def test_transformer_block(self):
        transformer = TransformerBlock(
            embed_dim=self.embed_dim,
            num_heads=self.num_heads,
            ff_dim=self.ff_dim
        )
        inputs = tf.random.normal((self.batch_size, self.seq_length, self.embed_dim))
        outputs, attention_weights = transformer(inputs, training=True)
        
        # Test output shapes
        self.assertEqual(outputs.shape, (self.batch_size, self.seq_length, self.embed_dim))
        self.assertEqual(attention_weights.shape, (self.batch_size, self.num_heads, self.seq_length, self.seq_length))

    def test_token_position_embedding(self):
        embedding_layer = TokenAndPositionEmbedding(
            maxlen=self.seq_length,
            vocab_size=self.vocab_size,
            embed_dim=self.embed_dim
        )
        inputs = tf.random.uniform((self.batch_size, self.seq_length), 
                                 dtype=tf.int32, 
                                 maxval=self.vocab_size)
        outputs = embedding_layer(inputs)
        
        # Test output shape
        self.assertEqual(outputs.shape, (self.batch_size, self.seq_length, self.embed_dim))

    def test_layer_configs(self):
        # Test get_config() for TransformerBlock
        transformer = TransformerBlock(
            embed_dim=self.embed_dim,
            num_heads=self.num_heads,
            ff_dim=self.ff_dim
        )
        config = transformer.get_config()
        self.assertEqual(config['embed_dim'], self.embed_dim)
        self.assertEqual(config['num_heads'], self.num_heads)
        self.assertEqual(config['mlp_dim'], self.ff_dim)

        # Test get_config() for TokenAndPositionEmbedding
        embedding_layer = TokenAndPositionEmbedding(
            maxlen=self.seq_length,
            vocab_size=self.vocab_size,
            embed_dim=self.embed_dim
        )
        config = embedding_layer.get_config()
        self.assertEqual(config['vocab_size'], self.vocab_size)
        self.assertEqual(config['max_len'], self.seq_length)
        self.assertEqual(config['embed_dim'], self.embed_dim)

if __name__ == '__main__':
    unittest.main()