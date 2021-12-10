import pandas as pd
import os
import minHash
import kmer_approach_joint as kmer_file

if __name__ == '__main__':

    label_df = pd.read_csv('dataset.csv', dtype=str)
    available_reads = []

    dataset = {}
    dataset['read_1'] = []
    dataset['read_2'] = []
    dataset['overlap'] = []
    dataset['vector_1'] = []
    dataset['vector_2'] = []

    solid_kmers_per_pair = []

    for i in range(0, 20000, 300):
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

        """
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
        """
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














