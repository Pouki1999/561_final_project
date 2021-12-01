import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
if __name__ == '__main__':

    reads = []
    os.makedirs('reads', exist_ok=True)
    os.makedirs('chromosomes', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    n = 0
    for record in SeqIO.parse('readsMappingToChr1.fa', 'fasta'):
        reads.append(record)
        if not os.path.isfile(os.path.join('reads', 'read_' + str(n) + '.fa')):
            SeqIO.write(record, open(os.path.join('reads', 'read_' + str(n) + '.fa'), 'w'), 'fasta')
            print('create file for read ' + str(n))
        n += 1
    n = 0
    os.makedirs('chromosomes', exist_ok=True)
    for record in SeqIO.parse('GCA_000227135.2_ASM22713v2_genomic.fna', 'fasta'):
        if record.id == 'FR799588.2_Leishmania_donovani,_chromosome_1':
            chr1 = record
        if not os.path.isfile(os.path.join('chromosomes', 'chromosome_' + str(n) + '.fa')):
            SeqIO.write(record, open(os.path.join('chromosomes', 'chromosome_' + str(n) + '.fa'), 'w'), 'fasta')
            print('create file for chromosome ' + str(n))
        n += 1

    blastn_path = r"C:\Program Files\NCBI\blast-2.12.0+\bin\blastn.exe"
    subject_path = r"C:\Users\gogom\PycharmProjects\561_final_project\chromosomes\chromosome_0.fa"
    for file in os.listdir('reads'):
        read_name = file.split('.')[0]
        query_path = r"C:\Users\gogom\PycharmProjects\561_final_project\reads\{}.fa".format(read_name)
        out_path = r"C:\Users\gogom\PycharmProjects\561_final_project\results\{}.xml".format(read_name)
        if not os.path.isfile(r"C:\Users\gogom\PycharmProjects\561_final_project\results\{}.xml".format(read_name)):
            cl = NcbiblastnCommandline(cmd=blastn_path, query=query_path, subject=subject_path, out=out_path, outfmt=5)
            cl()

    hsp_per_read = {}
    for result in os.listdir('results'):
        handle = open(os.path.join('results', result), 'r')
        blast_record = NCBIXML.read(handle)
        best = float('inf')
        for hsp in blast_record.alignments[0].hsps:
            if hsp.expect < best:
                hsp_per_read[result.split('.')[0]] = hsp
                best = hsp.expect
    print(hsp_per_read)

    for name, hsp in hsp_per_read.items():
        print(name)
        print((hsp.query_start, hsp.query_end))
        print(hsp.query_end - hsp.query_start)
        print((hsp.sbjct_start, hsp.sbjct_end))
        print(hsp.sbjct_end - hsp.sbjct_start)















