from Bio import SeqIO

if __name__ == '__main__':
    c = 0
    for record in SeqIO.parse('readsMappingToChr1.fa', 'fasta'):
        print(record.id)
        c += len(record)

    print(c)

    for record in SeqIO.parse('GCA_000227135.2_ASM22713v2_genomic.fna', 'fasta'):
        print(record.id)




