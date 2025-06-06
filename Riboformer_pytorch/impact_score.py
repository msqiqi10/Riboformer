# for the calculation of impact scores - PyTorch version
import os
import numpy as np
import argparse
import torch

from modules import TransformerBlock, TokenAndPositionEmbedding
from config import Config
from model import Riboformer

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    run = parser.add_argument_group(title='Prediction parameters',
                                    description='Parameters for Prediction '
                                                'using model')
    run.add_argument('-i', '--input_folder', default='GSE139036_disome', help='Input data folder')
    run.add_argument('-m', '--model_file', default='yeast_disome.pth', help='Model file')
    run.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda'],
                     help='Device to use for prediction')
    args = parser.parse_args()
    
    # Set device
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")
        
    # model parameters
    model_config = Config()
    
    print("--------------------------------------------------")
    print("Data loading")
    parpath = os.path.dirname(os.getcwd())

    # Load model
    model = Riboformer(model_config)
    
    if args.model_file.endswith('.pth'):
        checkpoint = torch.load(os.path.join(parpath, 'models', args.model_file), 
                               map_location=device)
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
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

    x_c = np.loadtxt(os.path.join(inputpath, xc_files[0]), delimiter="\t")
    print(f"Input data size is {len(x_c)}")

    x_c[:, :40] = x_c[:,0:40]/100 - 5
    x_c[:, 40] = x_c[:,40]/100

    # load ribosome pausing sites
    indices = np.loadtxt(os.path.join(parpath, 'datasets', args.input_folder, 
                                     'pause_indices.txt'), delimiter="\t").astype('int')[:, 0]
    print(len(indices), "ribosome pausing sites.")

    # Make initial predictions
    with torch.no_grad():
        seq_tensor = torch.tensor(x_c[:, -40:], dtype=torch.long, device=device)
        exp_tensor = torch.tensor(x_c[:, :40], dtype=torch.float32, device=device)
        y_pred = model((seq_tensor, exp_tensor)).cpu().numpy()
    
    print("--------------------------------------------------")
    print("Calculating sequence impact score")

    N_plot = len(indices)
    pos_start, pos_end = 0, 30
    window_size = 10
    N_rand = 50
    y_rand_mean = np.ones([N_plot, pos_end - pos_start])

    codon_list = np.append(np.arange(58), [60, 61, 63])
    
    # Prepare batch data for efficient processing
    batch_size = 1000  # Process in batches to manage memory
    all_results = []

    # calculate the sequence impact scores
    for i, idx in enumerate(indices):
        for j in range(pos_start, pos_end):
            
            # generate the random codons
            codon_rand = np.random.randint(low=0, high=61, size=[N_rand, window_size])
            switch_index = lambda t: codon_list[int(t)]
            vfunc = np.vectorize(switch_index)
            codon_rand = vfunc(codon_rand)
            
            # make predictions and record results
            x_codon_rand = np.tile(x_c[idx, -40:], (N_rand, 1))
            x_codon_rand[:, j:j + window_size] = codon_rand
            x_level = np.tile(x_c[idx, :40], (N_rand, 1))

            # Store for batch processing
            batch_data = {
                'seq': x_codon_rand,
                'exp': x_level,
                'original_pred': y_pred[idx],
                'idx': (i, j)
            }
            all_results.append(batch_data)

        if i % 1000 == 0:
            print(f"Finished preparing {i} pause sites.")
    
    print(f"Finished preparing {len(indices)} pause sites.")
    
    print("--------------------------------------------------")
    print("Model prediction.")
    
    # Process in batches
    y_rand_changes = np.zeros((N_plot, pos_end - pos_start))
    
    with torch.no_grad():
        for batch_idx in range(0, len(all_results), batch_size):
            batch_end = min(batch_idx + batch_size, len(all_results))
            batch_items = all_results[batch_idx:batch_end]
            
            # Combine all sequences and expressions in this batch
            batch_seqs = []
            batch_exps = []
            batch_info = []
            
            for item in batch_items:
                batch_seqs.append(item['seq'])
                batch_exps.append(item['exp'])
                batch_info.extend([(item['original_pred'], item['idx'])] * len(item['seq']))
            
            # Convert to tensors
            batch_seqs = np.vstack(batch_seqs)
            batch_exps = np.vstack(batch_exps)
            
            seq_tensor = torch.tensor(batch_seqs, dtype=torch.long, device=device)
            exp_tensor = torch.tensor(batch_exps, dtype=torch.float32, device=device)
            
            # Make predictions
            batch_preds = model((seq_tensor, exp_tensor)).cpu().numpy().squeeze()
            
            # Calculate changes and aggregate
            start_idx = 0
            for item in batch_items:
                end_idx = start_idx + N_rand
                item_preds = batch_preds[start_idx:end_idx]
                original_pred = item['original_pred']
                i, j = item['idx']
                
                # Calculate mean change for this position
                mean_change = np.mean(item_preds - original_pred)
                y_rand_changes[i, j - pos_start] = mean_change
                
                start_idx = end_idx
            
            if batch_idx % 10000 == 0:
                print(f"Processed {batch_idx} batches...")
    
    print("--------------------------------------------------")
    print("Finishing calculation and saving results.")
    np.savetxt(os.path.join(parpath, 'datasets', args.input_folder, 'SIS.txt'),
               y_rand_changes, delimiter="\t")

if __name__ == "__main__":
    main()
