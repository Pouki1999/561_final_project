from __future__ import print_function, division
import os
import pandas as pd
from skimage import io, transform
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader

class DNA_dataset(Dataset):
    def __init__(self, csv_file, dir):
        """
        Args:
          dir (string): Directory with all the DNA sequence pairs.
          each file should be a .txt or fasta file with a DNA sequence inside
          csv_file (string): Path to the csv file with all pair-wise labels.
            1st line of the file has to be title because data at line 0 is skipped
            each line should be, 'read#, read#, boolean'
        """
        self.path_DNA_pairs = dir
        self.labels = pd.read_file(csv_file)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        """
        get the 0-based idx-th element in the dataset
        returns a numpy array containing both sequences and the label for the pair of DNA sequences at that index

        DNA_dataset.DNA_pairs[0] = [ sequence1 = ACCTGTTACCC, sequence2 = GGTCAACACTTTA, label = False ]
        DNA_dataset.DNA_pairs[1] = [ sequence1 = AGTGACTGGAC, sequence2 = ACTGGACTCAACC, label = True ]
         .
         .
         .
        DNA_dataset.DNA_pairs[-1] = [ sequence1 = ACTGGATCAACC, sequence2 = AGTGACTGGACA, label = True ]
        """
        assert idx >= 0
        to_return = list()

        line_of_interest = str(self.labels.iloc(idx)).split()
        # line of interest will look something like this :
        # ['title_line', 'True', 'Name:', '(3,', '1),', 'dtype:', 'object'] where we
        # want to extract True, 3, and 1 as done using the re module below
        import re
        seq1_num = re.sub("[^0-9]", "", line_of_interest[3])
        seq2_num = re.sub("[^0-9]", "", line_of_interest[4])
        with open(os.path.join(self.path_DNA_pairs, '/pair_', seq1_num, '.txt'), 'r') as f:
            seq1 = f.readlines()
        with open(os.path.join(self.path_DNA_pairs, '/pair_', seq2_num, '.txt'), 'r') as f:
            seq2 = f.readlines()
        to_return.append(seq1)
        to_return.append(seq2)
        to_return.append(line_of_interest[1])
        return np.array(to_return)