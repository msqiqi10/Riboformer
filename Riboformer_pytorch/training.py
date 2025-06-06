import os
import argparse
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from torch.utils.tensorboard import SummaryWriter

from config import Config
from model import Riboformer

def main():
    
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

    test = parser.add_argument_group(title='Test parameters',
                                     description='Parameters for testing model')
    test.add_argument('-e', '--epoch', default=10, type=int, 
                      help='epoch number for model training')
    test.add_argument('-b', '--batch', default=256, type=int, 
                      help='batch size for model training')
    test.add_argument('-s', '--split', default=0.7, type=float,
                      help='proportion for the training dataset')
    test.add_argument('-l', '--learning', default=0.001, type=float, 
                      help='learning rate for model training')
    test.add_argument('-save', '--save', action='store_true', 
                      help='Save the model and training results')
    test.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda'],
                      help='Device to use for training')

    run = parser.add_argument_group(title='Prediction parameters',
                                    description='Parameters for Prediction using model')
    run.add_argument('-o', '--output_model', default='cm_mg_model2', help='Output PyTorch model')
    run.add_argument('-i', '--input_folder', default='GSE119104_Mg_buffer', help='Input data folder')

    args = parser.parse_args()
    print(args)
    
    # Set device
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")
    
    # initiate the model & the model configs
    model_config = Config()
    
    model = Riboformer(model_config)
    model = model.to(device)
    
    wsize = model_config.wsize

    print("--------------------------------------------------\nData loading")
    
    # You'll need to implement data loading here
    # This is a placeholder for the actual data loading logic
    # The data loading should be compatible with the existing data_processing.py
    
    # Example data loading (you'll need to adapt this to your actual data)
    # X_seq, X_exp, y = load_data(args.input_folder)  # Implement this function
    # 
    # # Convert to PyTorch tensors
    # X_seq = torch.tensor(X_seq, dtype=torch.long)
    # X_exp = torch.tensor(X_exp, dtype=torch.float32)
    # y = torch.tensor(y, dtype=torch.float32)
    
    # # Create data loaders
    # dataset = TensorDataset(X_seq, X_exp, y)
    # train_size = int(args.split * len(dataset))
    # val_size = len(dataset) - train_size
    # train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    # 
    # train_loader = DataLoader(train_dataset, batch_size=args.batch, shuffle=True)
    # val_loader = DataLoader(val_dataset, batch_size=args.batch, shuffle=False)
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=args.learning)
    
    # Optional: Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=3, factor=0.5)
    
    # Training loop
    print(f"Starting training for {args.epoch} epochs...")
    
    # You'll need to uncomment and adapt this training loop once you have data loading implemented
    # for epoch in range(args.epoch):
    #     model.train()
    #     train_loss = 0.0
    #     
    #     for batch_idx, (seq_batch, exp_batch, target_batch) in enumerate(train_loader):
    #         seq_batch = seq_batch.to(device)
    #         exp_batch = exp_batch.to(device)
    #         target_batch = target_batch.to(device)
    #         
    #         optimizer.zero_grad()
    #         outputs = model((seq_batch, exp_batch))
    #         loss = criterion(outputs.squeeze(), target_batch)
    #         loss.backward()
    #         optimizer.step()
    #         
    #         train_loss += loss.item()
    #     
    #     # Validation
    #     model.eval()
    #     val_loss = 0.0
    #     with torch.no_grad():
    #         for seq_batch, exp_batch, target_batch in val_loader:
    #             seq_batch = seq_batch.to(device)
    #             exp_batch = exp_batch.to(device)
    #             target_batch = target_batch.to(device)
    #             
    #             outputs = model((seq_batch, exp_batch))
    #             loss = criterion(outputs.squeeze(), target_batch)
    #             val_loss += loss.item()
    #     
    #     train_loss /= len(train_loader)
    #     val_loss /= len(val_loader)
    #     
    #     print(f'Epoch [{epoch+1}/{args.epoch}], Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}')
    #     
    #     # Step the scheduler
    #     scheduler.step(val_loss)
    
    # Save model if requested
    if args.save:
        torch.save({
            'model_state_dict': model.state_dict(),
            'model_config': model_config,
            'optimizer_state_dict': optimizer.state_dict(),
        }, f'{args.output_model}.pth')
        print(f"Model saved as {args.output_model}.pth")

def load_data(input_folder):
    """
    Load data from the specified folder.
    This function should be implemented to load your specific data format.
    Should return sequences, expressions, and targets as numpy arrays.
    """
    # Placeholder - implement actual data loading logic here
    # You can use the functions from data_processing.py
    pass

if __name__ == "__main__":
    main()
