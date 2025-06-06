#!/usr/bin/env python3
"""
PyTorch Model Test Script
Tests the Riboformer PyTorch model with sample data
"""

import torch
import numpy as np
from config import Config
from model import Riboformer
from modules import ConvTower, TransformerBlock, TokenAndPositionEmbedding

def test_modules():
    """Test individual modules"""
    print("Testing PyTorch modules...")
    
    # Test TokenAndPositionEmbedding
    print("Testing TokenAndPositionEmbedding...")
    embed = TokenAndPositionEmbedding(maxlen=40, vocab_size=64, embed_dim=8)
    test_input = torch.randint(0, 64, (2, 40))
    embed_output = embed(test_input)
    print(f"Embedding output shape: {embed_output.shape}")
    assert embed_output.shape == (2, 40, 8), f"Expected (2, 40, 8), got {embed_output.shape}"
    
    # Test ConvTower 1D
    print("Testing ConvTower 1D...")
    conv1d = ConvTower('1D', [32, 32, 8], 9, activation='relu')
    test_input_1d = torch.randn(2, 1, 40)
    conv1d_output = conv1d(test_input_1d)
    print(f"Conv1D output shape: {conv1d_output.shape}")
    
    # Test ConvTower 2D
    print("Testing ConvTower 2D...")
    conv2d = ConvTower('2D', [32, 32, 32], 5, activation='relu')
    test_input_2d = torch.randn(2, 1, 40, 8)
    conv2d_output = conv2d(test_input_2d)
    print(f"Conv2D output shape: {conv2d_output.shape}")
    
    # Test TransformerBlock
    print("Testing TransformerBlock...")
    transformer = TransformerBlock(embed_dim=8, num_heads=2, ff_dim=64)
    test_input_tf = torch.randn(2, 40, 8)
    tf_output, weights = transformer(test_input_tf)
    print(f"Transformer output shape: {tf_output.shape}")
    assert tf_output.shape == (2, 40, 8), f"Expected (2, 40, 8), got {tf_output.shape}"
    
    print("All module tests passed!")

def test_model():
    """Test the complete Riboformer model"""
    print("\nTesting complete Riboformer model...")
    
    # Create config
    config = Config()
    
    # Create model
    model = Riboformer(config)
    model.eval()
    
    # Create sample data
    batch_size = 4
    seq_data = torch.randint(0, config.vocab_size, (batch_size, config.wsize))
    exp_data = torch.randn(batch_size, config.wsize)
    
    print(f"Input sequence shape: {seq_data.shape}")
    print(f"Input expression shape: {exp_data.shape}")
    
    # Forward pass
    with torch.no_grad():
        output = model((seq_data, exp_data))
    
    print(f"Model output shape: {output.shape}")
    print(f"Model output: {output}")
    
    assert output.shape == (batch_size, 1), f"Expected ({batch_size}, 1), got {output.shape}"
    
    print("Model test passed!")

def test_model_training_mode():
    """Test model in training mode"""
    print("\nTesting model in training mode...")
    
    config = Config()
    model = Riboformer(config)
    model.train()
    
    # Create sample data
    batch_size = 2
    seq_data = torch.randint(0, config.vocab_size, (batch_size, config.wsize))
    exp_data = torch.randn(batch_size, config.wsize)
    
    # Forward pass
    output = model((seq_data, exp_data))
    
    # Test backpropagation
    loss = torch.mean(output ** 2)  # Simple loss for testing
    loss.backward()
    
    print(f"Training mode output shape: {output.shape}")
    print(f"Loss: {loss.item()}")
    
    # Check if gradients are computed
    has_gradients = any(p.grad is not None for p in model.parameters())
    assert has_gradients, "No gradients computed"
    
    print("Training mode test passed!")

def test_model_compatibility():
    """Test model compatibility with different input sizes"""
    print("\nTesting model compatibility...")
    
    config = Config()
    model = Riboformer(config)
    model.eval()
    
    # Test with different batch sizes
    for batch_size in [1, 3, 16]:
        seq_data = torch.randint(0, config.vocab_size, (batch_size, config.wsize))
        exp_data = torch.randn(batch_size, config.wsize)
        
        with torch.no_grad():
            output = model((seq_data, exp_data))
        
        assert output.shape == (batch_size, 1), f"Failed for batch_size {batch_size}"
        print(f"Batch size {batch_size}: OK")
    
    print("Compatibility tests passed!")

if __name__ == "__main__":
    print("Starting PyTorch Riboformer tests...")
    
    try:
        test_modules()
        test_model()
        test_model_training_mode()
        test_model_compatibility()
        
        print("\n" + "="*50)
        print("All tests passed successfully!")
        print("PyTorch Riboformer model is working correctly.")
        print("="*50)
        
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
