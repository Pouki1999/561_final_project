# kmer_count = {'AAA':1, 'AAT':3, "AAC":2, "AAG":2, 'ATA':2, 'ATT':3, 'ATC':2, 'ATG':5, 'ACA':3, 'ACT':2, 'ACC':4, 'ACG':3,
#          'CAA':6, 'CAT':5, 'CAC':1, 'CAG':9}

read_0 = 'AATACGATCGCGGGATAT'
read_1 = 'CGAATCGACTTTAACCCA'

#####
#this counting and locating step is performed once for each read in the whole read set in the very begining
#we then obtain one very large hashtable that contains the positions of the kmers in all reads
#####
def count_n_locate_kmers(read, read_id, k, kmer_count={}, kmer_locations={}, spacing=1):
    # TODO: adapt method to have customizable length spacing between each kmer
    """
    :param read: string of nucleotides making up a DNA read
    :param read_id: string or int identifying which of the two reads is being processed
    :param k: int kmer length
    :param kmer_count: dict key is a kmer that has been found and value is that kmer's frequency
    :param kmer_locations: dict key is a kmer that has been found and value is a dict with that kmer's locations in all
    reads that contain it
    :param spacing: number of nucleotides shifted between each kmer
    :return: two dicts giving the frequency and locations in the reads of each kmer
    """
    i = 0
    # iterate once per amount of kmers in the input sequence, given the spacing
    while i < (len(read) - (k * spacing) + 1):
        # getting a kmer
        curr = read[i:i + k]
        if curr in kmer_count:
            kmer_count[curr] += 1
            if read_id not in kmer_locations[curr]:
                kmer_locations[curr] = {read_id:[i]}
            else:
                kmer_locations[curr][read_id].append(i)
        else:
            kmer_count[curr] = 1
            kmer_locations[curr] = {}
            kmer_locations[curr] = {read_id:[i]}
        i += spacing
    return kmer_count, kmer_locations


def find_solid_kmers(k_mers, kmer_locations, theta=0.9):
    """
    :param k_mers: dict key is a kmer that has been found and value is that kmer's frequency
    :param kmer_locations: dict key is a kmer that has been found and value is a list of that kmer's locations in both
    reads
    :param theta: cutoff value that determines how big f_max gets
    :return: two dicts giving the frequency and locations in the reads of each kmer that is deemed solid f_min < frequency(kmer) < f_max
    """
    # removing all the k_mers that have frequency 1
    no_singletons = {k: k_mers[k] for k in k_mers if k_mers[k] > 1}
    sum_h = len(no_singletons)
    sum_f, f_max = 0, 2
    # determining the value of f_max
    while sum_f <= int(sum_h * theta):
        sum_f += len({k for k in no_singletons if no_singletons[k] == f_max})
        f_max += 1
    # removing all kmers that don't fall in the interval [f_min; f_max] and their locations
    solid_kmers = {k: no_singletons[k] for k in no_singletons if no_singletons[k] <= f_max}
    solid_kmer_locations = {k: kmer_locations[k] for k in solid_kmers}
    # return 2, f_max, solid_kmers
    return solid_kmers, solid_kmer_locations


def build_common_kmer_set(read_A, read_B, read_B_id, solid_kmers, k, spacing=1):
    common_set = {}
    i = 0
    while i < (len(read_A) - (k * spacing) + 1):
        curr = read_A[i:i + k]
        if read_B_id in solid_kmers[curr].keys() and solid_kmers[curr][read_B_id] <= 2 :
            common_set.append(curr)
        i+=spacing
    return common_set


kmer_size = 3
read_set = [ read_0, read_1]
kmer_count, kmer_location = {}, {}
for i in range(len(read_set)):
    kmer_count, kmer_locations = count_n_locate_kmers(read_0, str(i), kmer_size, kmer_count, kmer_locations)
# kmer_count_B, kmer_locations_B = count_n_locate_kmers(read_1, '1', kmer_size)
solid_kmers, solid_kmer_locations = find_solid_kmers(kmer_count, kmer_locations)
# solid_kmers_B, solid_kmer_locations_B = find_solid_kmers(kmer_count_B, kmer_locations_B)


# testing whether the kmer locating process succeeded
# kmer = 'CGA'
# print(kmer_count[kmer] == len(kmer_locations[kmer]))
# print("Truth: CGA, predicted:")
# for i in range(len(kmer_locations[kmer])):
#     print(read_0[kmer_locations[kmer][i]:kmer_locations[kmer][i] + kmer_size])


