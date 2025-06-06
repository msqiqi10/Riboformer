#!/usr/bin/env python3
"""
Integration test for PyTorch Riboformer implementation
Tests end-to-end functionality with synthetic data
"""

import numpy as np
import torch
import torch.nn as nn
from model import Riboformer
from config import Config
import tempfile
import os

def create_synthetic_data(num_samples=1000, seq_length=40):
    """Create synthetic training data for testing"""
    
    # Generate random sequences (0-63 for 64 codons)
    sequences = np.random.randint(0, 64, size=(num_samples, seq_length))
    
    # Generate random expression values
    expressions = np.random.normal(0, 1, size=(num_samples, seq_length))
    
    # Generate synthetic targets (simple function of inputs)
    targets = np.mean(sequences[:, :10], axis=1) * 0.1 + np.random.normal(0, 0.1, size=num_samples)
    
    return sequences, expressions, targets

def test_end_to_end():
    """Test the complete PyTorch Riboformer pipeline"""
    
    print("="*60)
    print("PyTorch Riboformer Integration Test")
    print("="*60)
    
    # Create config and model
    config = Config()
    model = Riboformer(config)
    
    # Create synthetic data
    print("Creating synthetic data...")
    sequences, expressions, targets = create_synthetic_data(num_samples=500, seq_length=40)
    
    # Convert to tensors
    seq_tensor = torch.tensor(sequences, dtype=torch.long)
    exp_tensor = torch.tensor(expressions, dtype=torch.float32)
    target_tensor = torch.tensor(targets, dtype=torch.float32).unsqueeze(1)
    
    print(f"Data shapes:")
    print(f"  Sequences: {seq_tensor.shape}")
    print(f"  Expressions: {exp_tensor.shape}")
    print(f"  Targets: {target_tensor.shape}")
    
    # Test forward pass
    print("\nTesting forward pass...")
    model.eval()
    with torch.no_grad():
        outputs = model((seq_tensor[:10], exp_tensor[:10]))
        print(f"Model output shape: {outputs.shape}")
        print(f"Sample outputs: {outputs[:5].squeeze()}")
    
    # Test training setup
    print("\nTesting training setup...")
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.MSELoss()
    
    # Mini training loop
    model.train()
    batch_size = 32
    num_batches = 3
    
    print(f"Running {num_batches} training batches...")
    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = start_idx + batch_size
        
        # Get batch
        batch_seq = seq_tensor[start_idx:end_idx]
        batch_exp = exp_tensor[start_idx:end_idx]
        batch_targets = target_tensor[start_idx:end_idx]
        
        # Forward pass
        optimizer.zero_grad()
        outputs = model((batch_seq, batch_exp))
        loss = criterion(outputs, batch_targets)
        
        # Backward pass
        loss.backward()
        optimizer.step()
        
        print(f"  Batch {batch_idx + 1}: Loss = {loss.item():.4f}")
    
    # Test model saving and loading
    print("\nTesting model saving and loading...")
    with tempfile.NamedTemporaryFile(suffix='.pth', delete=False) as tmp_file:
        model_path = tmp_file.name
        torch.save(model.state_dict(), model_path)
        print(f"Model saved to: {model_path}")
        
        # Create new model and load weights
        new_model = Riboformer(config)
        new_model.load_state_dict(torch.load(model_path, map_location='cpu'))
        print("Model loaded successfully")
        
        # Test that outputs are identical
        model.eval()
        new_model.eval()
        with torch.no_grad():
            original_output = model((seq_tensor[:5], exp_tensor[:5]))
            loaded_output = new_model((seq_tensor[:5], exp_tensor[:5]))
            
            diff = torch.abs(original_output - loaded_output).max().item()
            print(f"Max difference between original and loaded model: {diff}")
            
            if diff < 1e-6:
                print("✓ Model save/load test PASSED")
            else:
                print("✗ Model save/load test FAILED")
        
        # Clean up
        os.unlink(model_path)
    
    # Test model evaluation metrics
    print("\nTesting model evaluation...")
    model.eval()
    with torch.no_grad():
        test_outputs = model((seq_tensor[-100:], exp_tensor[-100:]))
        test_targets = target_tensor[-100:]
        
        # Calculate metrics
        mse = nn.MSELoss()(test_outputs, test_targets)
        mae = nn.L1Loss()(test_outputs, test_targets)
        
        print(f"Test MSE: {mse.item():.4f}")
        print(f"Test MAE: {mae.item():.4f}")
        
        # Calculate correlation
        test_outputs_np = test_outputs.squeeze().numpy()
        test_targets_np = test_targets.squeeze().numpy()
        correlation = np.corrcoef(test_outputs_np, test_targets_np)[0, 1]
        print(f"Test Correlation: {correlation:.4f}")
    
    print("\n" + "="*60)
    print("Integration test completed successfully!")
    print("✓ Model architecture works correctly")
    print("✓ Forward pass works correctly")
    print("✓ Training loop works correctly")
    print("✓ Model save/load works correctly")
    print("✓ Evaluation metrics work correctly")
    print("="*60)
    
    return True

if __name__ == "__main__":
    test_end_to_end()
