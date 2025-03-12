
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from riboformer_utils import *
# Read the DNA sequence and the GFF file
sequence_iter = fasta_iter("/home/zzz0054/chen_data/Riboformer/datasets/AA/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta")
seq_0 = sequence_iter.__next__()[1]
# for header, seq in sequence_iter:
#     print("Genome seq length:", len(seq))
    
annotations = np.genfromtxt("/home/zzz0054/chen_data/Riboformer/datasets/AA/gff.csv", delimiter="\t")

pause_score_benchmark = {}
L_score_benchmark0 = {}
L_score_benchmark1 = {}


def process_experiment(densities, zc, ypred, annotations, seq, pause_scores_aa):
    """Process an experiment by calculating pause scores and returning the results."""
    # densities = read_gene_densities(file_paths['data_folder'], file_paths[exp_name], suffixes)
    # zc = np.loadtxt(file_paths['data_folder'] + file_paths[exp_name] + file_io['z_index'], delimiter="\t")
    # ypred = np.loadtxt(file_paths['data_folder'] + 'Sym_' + file_paths[exp_name] + file_io['y_pred'], delimiter="\t")

    # Input 6 -> E site, 3 -> P site, 0 -> A site
    pause_scores_1 = get_pause_score(0, 0, densities, annotations, seq, ypred, zc, np.arange(len(annotations)), pred=1)
    pause_scores_0 = get_pause_score(0, 0, densities, annotations, seq, ypred, zc, np.arange(len(annotations)), pred=0)

    pause_scores_aa = {k: pause_scores_aa[k] for k in sorted(pause_scores_aa)}
    pause_scores_total = []

    for k in pause_scores_aa.keys():
        codon_0_scores = [np.append([], pause_scores_0[s]) for s in pause_scores_aa[k]]
        codon_1_scores = [np.append([], pause_scores_1[s]) for s in pause_scores_aa[k]]
        pause_scores_total.append([np.mean(np.concatenate(codon_0_scores)), np.mean(np.concatenate(codon_1_scores))])

    H_codon_list = ['CTA', 'CTG', 'CTC', 'CTT', 'TTA', 'TTG']
    l_scores_0 = [np.mean(pause_scores_0[s]) for s in H_codon_list]
    l_scores_1 = [np.mean(pause_scores_1[s]) for s in H_codon_list]

    return np.array(pause_scores_total).transpose(), l_scores_0, l_scores_1


densities = read_gene_densities_new("/home/zzz0054/chen_data/Riboformer/datasets/AA/AuxinA0FootprintR4_DNA_f/AuxinA0FootprintR4_DNA_f_1_var.wig", "/home/zzz0054/chen_data/Riboformer/datasets/AA/AuxinA0FootprintR4_DNA_r/AuxinA0FootprintR4_DNA_r_1_var.wig")

zc = np.loadtxt("/home/zzz0054/chen_data/Riboformer/datasets/AA/zc.txt", delimiter="\t")
ypred = np.loadtxt("/home/zzz0054/chen_data/Riboformer/datasets/AA/model_prediction.txt", delimiter="\t")
idx = np.loadtxt('/home/zzz0054/chen_data/Riboformer/datasets/AA/test_idx.txt', delimiter="\t").astype(int)
zc = zc[idx]
zc.shape


densities.shape, zc.shape, ypred.shape, annotations.shape


ypred.shape


# Pearson correlation coefficient
np.corrcoef(ypred[:, 0], ypred[:, 1])[0,1]


pause_scores, l_scores_0, l_scores_1 = process_experiment(densities, zc, ypred[:,1], annotations, seq_0, pause_scores_aa)

# the index int(z_c[p, 1]) * 3 + 1 - 30 + (read_a_site - a_site) of r_pro_2 goes out of bounds, cannot be debugged,
# r_pro_2 shape is defined by dwig[start - 1 + a_site:end + a_site, 1].shape, which is 49909, 51210

# here p is m descriped m beblow, from data_processing.py
# i: Index of the gene or data segment being processed (from an outer loop).
# m: Position within the data sequence RD1 for the current gene or segment.

print('done')



