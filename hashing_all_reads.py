import pandas as pd
import os
import minHash
import kmer_approach_joint as kmer_file

if __name__ == '__main__':

    label_df = pd.read_csv('dataset.csv', dtype=str)
    available_reads = []
    for i in range(label_df.shape[0]):
        print(i)
        if label_df['read_1'][i] not in available_reads:
            available_reads.append(label_df['read_1'][i])
        if label_df['read_2'][i] not in available_reads:
            available_reads.append(label_df['read_2'][i])

    print(available_reads)

    dataset = {}
    dataset['read_1'] = []
    dataset['read_2'] = []
    dataset['overlap'] = []
    dataset['vector_1'] = []
    dataset['vector_2'] = []

    solid_kmers_per_pair = []

    for read in available_reads:
        read1 = label_df['read_1'][i]
        read2 = label_df['read_2'][i]
        seq1 = ''
        with open(os.path.join('reads', 'read_{}.fa'.format(read1)), 'r') as file:
            for line in file.readlines():
                if line[0] != '>':
                    seq1 += line.replace('\n', '')
        seq2 = ''
        with open(os.path.join('reads', 'read_{}.fa'.format(read2)), 'r') as file:
            for line in file.readlines():
                if line[0] != '>':
                    seq2 += line.replace('\n', '')
        print(len(seq1))
        print(len(seq2))
        k = 16
        num_hash_functions = 200

        kmer_dict_1 = minHash.get_kmer_list(seq1, k)
        kmer_dict_2 = minHash.get_kmer_list(seq2, k)

        max_len = max(len(kmer_dict_1.keys()), len(kmer_dict_2.keys()))
        c = minHash.nextprime(max_len)
        if c < num_hash_functions:
            c = minHash.nextprime(num_hash_functions)
        hash_func_params = minHash.generate_hash_functions(num_hash_functions, c)
        hash_set_1, hash_set_2 = [], []
        cur_values_1, cur_values_2 = list(kmer_dict_1.values()), list(kmer_dict_2.values())

        hash_set_1.append(cur_values_1)
        hash_set_2.append(cur_values_2)

        for params in hash_func_params:
            cur_values_1 = minHash.get_hash_set(params[0], params[1], c, cur_values_1)
            hash_set_1.append(cur_values_1)
            cur_values_2 = minHash.get_hash_set(params[0], params[1], c, cur_values_2)
            hash_set_2.append(cur_values_2)

        #print(hash_set_1)
        #print(hash_set_2)

        sketch_1 = minHash.get_hash_sketch(hash_set_1)
        sketch_2 = minHash.get_hash_sketch(hash_set_2)

        print(list(sketch_1))
        print(list(sketch_2))




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












