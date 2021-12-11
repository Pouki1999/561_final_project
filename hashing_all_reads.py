import pandas as pd
import os
import minHash
import kmer_approach_joint as kmer_file
import numpy as np

if __name__ == '__main__':

    df = pd.read_csv('dataset.csv')
    print(df)

    train_percentage = 0.8

    idx_list = []
    for i in range(df.shape[0]):
        idx_list.append(i)

    seed = 1
    np.random.seed(seed=seed)
    train_idx = []
    test_idx = []
    for idx in idx_list:
        if np.random.uniform(0.0, 1.0) < train_percentage:
            train_idx.append(idx)
        else:
            test_idx.append(idx)
    print(len(train_idx))
    print(len(test_idx))
    #print(train_idx)
    #print(test_idx)

    new_train = []
    for idx in train_idx:
        if np.random.uniform(0.0, 1.0) <= 0.1:
            new_train.append(idx)

    new_test = []
    for idx in test_idx:
        if np.random.uniform(0.0, 1.0) <= 0.1:
            new_test.append(idx)


    print(len(new_train))
    print(len(new_test))

    print(str(df['read_1'][idx]) + ',' + str(df['read_2'][idx]) + ',' + str(df['overlap'][idx]))

    with open('test_dataset.csv', 'w') as test_file:
        test_file.write('read_1,read_2,overlap' + '\n')
        for idx in new_test:
            test_file.write(str(df['read_1'][idx]) + ',' + str(df['read_2'][idx]) + ',' + str(df['overlap'][idx]) + '\n')

    with open('train_dataset.csv', 'w') as train_file:
        train_file.write('read_1,read_2,overlap' + '\n')
        for idx in new_train:
            train_file.write(str(df['read_1'][idx]) + ',' + str(df['read_2'][idx]) + ',' + str(df['overlap'][idx]) + '\n')


    """
    
    train_idx = list()
    test_idx = list()
    for idx in idx_list:
        if np.random.uniform(0.0, 1.0) < train_percentage:
            train_idx.append(idx)
        else:
            test_idx.append(idx)
    
    label_df = pd.read_csv('dataset.csv', dtype=str)
    #get all reads in dataset
    available_reads = []
    for i in range(int(label_df.shape[0]/10)):
        print('try')
        if label_df['read_1'][i] not in available_reads:
            print(i)
            available_reads.append(label_df['read_1'][i])
        if label_df['read_2'][i] not in available_reads:
            print(i)
            available_reads.append(label_df['read_2'][i])

    max_kmers = 0
    kmer_dict_list = []

    k = 16
    num_hash_functions = 200
    #determine max number of kmers in a read (to find c value for hash functions)
    for read in available_reads:
        print(read)
        seq = ''
        with open(os.path.join('reads', 'read_{}.fa'.format(read)), 'r') as file:
            for line in file.readlines():
                if line[0] != '>':
                    seq += line.replace('\n', '')

        kmer_dict = minHash.get_kmer_list(seq, k)
        kmer_dict_list.append(kmer_dict)

        if len(kmer_dict) > max_kmers:
            max_kmers = len(kmer_dict)

    print('max is:', max_kmers)
    sketch_list = []

    # determine c (next biggest prime and parameters a and b for each hash function)
    c = minHash.nextprime(max_kmers)
    if c < num_hash_functions:
        c = minHash.nextprime(num_hash_functions)
    hash_func_params = minHash.generate_hash_functions(num_hash_functions, c)

    os.makedirs('sketches', exist_ok=True)

    for i, read in enumerate(available_reads):

        kmer_dict = kmer_dict_list[i]

        hash_set = []
        cur_values = list(kmer_dict.values())

        hash_set.append(cur_values)

        for params in hash_func_params:
            cur_values = minHash.get_hash_set(params[0], params[1], c, cur_values)
            hash_set.append(cur_values)

        #print(hash_set_1)

        sketch = minHash.get_hash_sketch(hash_set)

        print(str(list(sketch)))

        outfile = open(os.path.join('sketches', 'sketch_{}.txt'.format(read)), 'w')
        outfile.write(str(list(sketch)))
        outfile.close()

    """


    """
        kmers_count1, kmers_pos1 = {}, {}
        kmers_count1, kmers_pos1 = kmer_file.count_n_locate_kmers(read=seq1, read_id=read1,
                                                             kmer_count=kmers_count1,
                                                             kmer_locations=kmers_pos1)
        kmers_count2, kmers_pos2 = kmer_file.count_n_locate_kmers(read=seq2, read_id=read2,
                                                             kmer_count=kmers_count1,
                                                             kmer_locations=kmers_pos1)

        solid_kmers, solid_kmers_pos = kmer_file.find_solid_kmers(kmers_count2, kmers_pos2)
        print(solid_kmers)
        print(solid_kmers_pos)

        
        kmers = list(solid_kmers.keys()).copy()
        for kmer in kmers:
            reads_present = []
            for pos in solid_kmers_pos[kmer]:
                if pos[0] not in reads_present:
                    reads_present.append(pos[0])
            if len(reads_present) < 2:
                del solid_kmers_pos[kmer]
                del solid_kmers[kmer]
        print(solid_kmers)
        
        kmers_read_1 = []
        kmers_read_2 = []
        kmers = list(solid_kmers.keys()).copy()
        for kmer in kmers:
            reads_present = []
            for pos in solid_kmers_pos[kmer]:
                if pos[0] not in reads_present:
                    reads_present.append(pos[0])
            if len(reads_present) == 2:
                kmers_read_1.append(kmer)
                kmers_read_2.append(kmer)
            else:
                if read1 in reads_present:
                    kmers_read_1.append(kmer)
                elif read2 in reads_present:
                    kmers_read_2.append(kmer)
        print(kmers_read_1)
        print(kmers_read_2)

        solid_kmers_per_pair.append((kmers_read_1, kmers_read_2))

        """












