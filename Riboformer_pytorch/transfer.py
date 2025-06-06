import os
import numpy as np
import argparse
import torch
from tqdm import tqdm

from modules import TransformerBlock, TokenAndPositionEmbedding
from config import Config
from model import Riboformer


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Script for making predictions using a pre-trained PyTorch model',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    run = parser.add_argument_group(title='Prediction Parameters',
                                    description='Parameters for prediction using the model')
    run.add_argument('-i', '--input_folder', help='Input data folder')
    run.add_argument('-m', '--model_file', help='Model file (.pth)')
    run.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda'],
                     help='Device to use for prediction')
    args = parser.parse_args()

    # Set device
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")

    # Initialize model parameters
    model_config = Config()

    print("--------------------------------------------------")
    parpath = os.path.dirname(os.getcwd())

    # Load the pre-trained model
    model = Riboformer(model_config)
    
    if args.model_file.endswith('.pth'):
        checkpoint = torch.load(os.path.join(parpath, 'models', args.model_file), 
                               map_location=device)
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        # Try to load as a simple state dict
        model.load_state_dict(torch.load(os.path.join(parpath, 'models', args.model_file), 
                                        map_location=device))
    
    model.to(device)
    model.eval()

    inputpath = os.path.join(parpath, 'datasets', args.input_folder)
    all_files = os.listdir(inputpath)
    xc_files = [f for f in all_files if f.endswith('xc.txt')]

    if len(xc_files) < 1:
        print("Input data not found. Please use the data_processing.py to generate the input dataset.")
        return

    if len(xc_files) > 1:
        print("Multiple input datasets exist.")
    
    # Data loading progress bar
    with tqdm(total=len(xc_files), desc="Data Loading", unit="file") as pbar:
        x_c = None
        for xc_file in xc_files:
            xc_path = os.path.join(inputpath, xc_file)
            xc_data = np.loadtxt(xc_path, delimiter="\t")
            if x_c is None:
                x_c = xc_data
            else:
                x_c = np.concatenate((x_c, xc_data))
            pbar.update(1)

    if x_c is None:
        print("No data loaded. Exiting.")
        return

    # Data preprocessing (same as original)
    x_c[:, :40] = x_c[:, 0:40] / 100 - 5
    x_c[:, 40] = x_c[:, 40] / 100

    print("--------------------------------------------------")

    # Model prediction progress bar
    with torch.no_grad():
        with tqdm(total=len(x_c), desc="Model Prediction", unit="data") as pbar:
            y_pred = []
            batch_size = 200000
            for i in range(0, len(x_c), batch_size):
                # Prepare batch data
                seq_batch = torch.tensor(x_c[i:i + batch_size, -40:], 
                                       dtype=torch.long, device=device)
                exp_batch = torch.tensor(x_c[i:i + batch_size, :40], 
                                       dtype=torch.float32, device=device)
                
                # Make predictions
                batch_pred = model((seq_batch, exp_batch))
                y_pred.extend(batch_pred.cpu().numpy())
                pbar.update(len(seq_batch))

    # Save predictions
    np.savetxt(os.path.join(parpath, 'datasets', args.input_folder, 'model_prediction.txt'),
               y_pred, delimiter="\t")
    print(f"Predictions saved to {os.path.join(parpath, 'datasets', args.input_folder, 'model_prediction.txt')}")


if __name__ == "__main__":
    main()
