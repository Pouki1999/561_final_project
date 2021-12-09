import numpy as np
from random import randint
import math

test_array = [
    [22, 52, 12, 6],
    [33, 25, 64, 82],
    [91, 14, 42, 1],
]

xorshift_seed = 30  # The initial seed should go here
def xorshift():
    global xorshift_seed
    xorshift_seed ^= xorshift_seed << 13
    xorshift_seed ^= xorshift_seed >> 17
    xorshift_seed ^= xorshift_seed << 5
    xorshift_seed %= int("ffffffff", 16)  # The modulus limits it to a 32-bit number
    return xorshift_seed
# for i in range(5):
#   # print(xorshift())
#   print(hash("CATTCC"))
# print(len(str(hash("CATTCC"))))

def get_kmer_list( read, k, spacing=1):
    i = 0
    kmer_dict = {}
    # iterate once per amount of kmers in the input sequence, given the spacing
    while i < int(math.ceil( ( len(read) - k + 1) / spacing ) ):
        curr = read[i:i + k]
        kmer_dict[curr] = hash(str(curr))
        i += spacing
    return kmer_dict

##taken from https://www.codegrepper.com/code-examples/python/next+prime+number+in+python
#finds the next prime number greater than the input number
def nextprime(n):
    prime = 0
    n += 1
    for i in range(2, int(n ** 0.5) + 2):
        if n % i == 0:
            prime = 0
            break
        else:
            prime = 1
    if prime == 1:
        return n
    else:
        return nextprime(n)
#inspired from https://stackoverflow.com/questions/24676237/generating-random-hash-functions-for-lsh-minhash-algorithm
def generate_hash_functions(how_many, limit):
    a = []
    b = []
    for i in range(how_many - 1):
        num = randint(0, limit)
        while num in a:
            num = randint(0, limit)
        a.append(num)
        num = randint(0, limit)
        while num in b:
            num = randint(0, limit)
        b.append(num)
    return list( [ list(x) for x in zip(a, b) ] )

#returns the hash set of the input set of kmer hash values
def get_hash_set(a, b, c, kmers):
    hash_set = []
    for kmer in kmers:
        hash_set.append( ( (a*kmer)+b) %c )
    return hash_set

#returns the hash sketch of the input set of hash sets for a sequence
def get_hash_sketch(hash_sets):
    return np.amin(np.array(hash_sets), axis=1)

def main():
    #params
    read_0 = 'AATACGATCGCGGGATAT'
    read_1 = 'CGAATCGACTTTAACCCA'
    kmer_size = 16
    num_hash_functions = 200

    #reading all the kmers in each read and storing their hashed values in a python dictionary
    kmer_dict_0 = get_kmer_list(read_0, kmer_size)
    kmer_dict_1 = get_kmer_list(read_1, kmer_size)

    #finding the appropriate c for hash function of the form h(x) = ((a*x)+b)%c
    #where c is the next largest prime greater than the number of k-mers in the longest of the two reads
    #where a and b are random integers in the interval [0,c]
    max_len = max( len(kmer_dict_0.keys()), len(kmer_dict_1.keys()) )
    c = nextprime( max_len )
    if c < num_hash_functions:
        c = nextprime(num_hash_functions)
    #finding the a's and b's for the above mentioned hash functions
    hash_func_params = generate_hash_functions(num_hash_functions, c)
    # print("hash_func_params=", hash_func_params )
    hash_set_0, hash_set_1 = [], []
    cur_values_0, cur_values_1 = list(kmer_dict_0.values()), list( kmer_dict_1.values() )
    hash_set_0.append(cur_values_0)
    hash_set_1.append(cur_values_1)
    for params in hash_func_params:
        cur_values_0 = get_hash_set(params[0], params[1], c, cur_values_0)
        hash_set_0.append( cur_values_0 )
        cur_values_1 = get_hash_set(params[0], params[1], c, cur_values_1)
        hash_set_1.append( cur_values_1 )
    sketch_0 = get_hash_sketch(hash_set_0)
    sketch_1 = get_hash_sketch(hash_set_1)

    # print( kmer_dict_0.keys())
    # print(list([x for x in zip(kmer_dict_0.keys(),kmer_dict_0.values())]))

if __name__ == '__main__':
    main()