# 561 Final Project Fall 2021
## Overlap discovery with Pacific Biosciences sequencing reads

***

#### Séraphin Bassas
#### Georgi Maslyankov
#### Emma Rousseau


***

##### McGill University
##### Dec. 10th 2021

Groud_truth.py: File used to generate individual fasta file and to blast reads on the first chromosome. Also, the file contained the code
                used to choose the reads that will be composing the dataset (as explained in the paper), but due to a awkward use of git 
                and an unfortunate merge, that part of the file was lost.

minHash.py : File containing helper functions for the kmers finding inside the reads and their hashing
hashing_all_reads.py: File used to produce test_dataset.csv and train_dataset.csv by producing the sketch for each read and dividing them in 20% test set and 80% train set 
dataset.csv: Main dataset. Contains all the pairs of reads that were selected and a boolean value indicating if they overlap or not on chromosome 1
test_dataset.csv: sketches for all reads in the test set
train_dataset.csv: sketches for all reads in the train set

Model_training.py: Code for the convolutional neural network. Disclaimer: the code was run on the Google colab. We are not sure if it is formatted correctly 
                    for it to be run elsewhere.



    
