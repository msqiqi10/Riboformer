# Riboformer PyTorch

This is a PyTorch implementation of the Riboformer model, refactored from the original TensorFlow/Keras version.

## Files

- `config.py` - Model configuration parameters
- `modules.py` - PyTorch implementations of model components (ConvTower, TransformerBlock, TokenAndPositionEmbedding)
- `model.py` - Main Riboformer model architecture
- `training.py` - Training script for the model
- `transfer.py` - Script for making predictions with pre-trained models
- `impact_score.py` - Script for calculating sequence impact scores
- `model_test.py` - Test script to verify model functionality
- `data_processing.py` - Data processing utilities (compatible with original)
- `codon_table.json` - Codon mapping table

## Key Changes from TensorFlow Version

### Architecture
- Converted from Keras to PyTorch nn.Module
- Updated tensor operations and dimensions to PyTorch conventions
- Maintained the same model architecture and functionality

### Modules
- **ConvTower**: Converted Conv1D/Conv2D layers with proper PyTorch syntax
- **TransformerBlock**: Implemented using PyTorch's MultiheadAttention
- **TokenAndPositionEmbedding**: Used PyTorch Embedding layers with proper position encoding

### Model Structure
- Input handling adapted for PyTorch tensor format (batch_first=True)
- Proper tensor reshaping and permutations for convolution operations
- Maintained compatibility with original data format

### Training
- Implemented PyTorch-style training loop with DataLoader
- Added support for GPU/CPU device selection
- Included learning rate scheduling and model checkpointing

## Usage

### Testing the Model
```bash
python model_test.py
```

### Training
```bash
python training.py -e 50 -b 256 -l 0.001 --save -i GSE119104_Mg_buffer
```

### Making Predictions
```bash
python transfer.py -i GSE119104_Mg_buffer -m model_name.pth
```

### Calculating Impact Scores
```bash
python impact_score.py -i GSE139036_disome -m yeast_disome.pth
```

## Requirements

- PyTorch >= 1.8.0
- NumPy
- tqdm
- BioPython (for data processing)
- BCBio-GFF (for data processing)

## Data Compatibility

The PyTorch version maintains full compatibility with the original data processing pipeline and data formats. You can use the same datasets and preprocessing scripts from the original TensorFlow version.

## Model Checkpoints

Models are saved as `.pth` files containing:
- `model_state_dict`: Model parameters
- `model_config`: Configuration object
- `optimizer_state_dict`: Optimizer state (for resuming training)

## Device Support

The implementation automatically detects and uses GPU acceleration when available, falling back to CPU when necessary. You can also explicitly specify the device using the `--device` argument.
