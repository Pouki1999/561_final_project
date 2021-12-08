import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

results_dir = 'new_results1'

blastn_path = r"/Users/seraphinbassas/Documents/Mcgill 2021-2022/COMP-561/ncbi-blast-2.12.0+/bin/blastn"
subject_path = r"/Users/seraphinbassas/PycharmProjects/561_final_project/chromosomes/chromosome_0.fa"
for file in os.listdir('reads'):
    read_name = file.split('.')[0]
    query_path = r"/Users/seraphinbassas/PycharmProjects/561_final_project/reads/{}.fa".format(read_name)
    out_path = r"/Users/seraphinbassas/PycharmProjects/561_final_project/{}/{}.xml".format(results_dir ,read_name)
    if not os.path.isfile(r"C:/Users/seraphinbassas/PycharmProjects/561_final_project/{}/{}.xml".format(results_dir, read_name)):
        cl = NcbiblastnCommandline(cmd=blastn_path, query=query_path, subject=subject_path, out=out_path, outfmt=5, gapopen=-2, gapextend=-1)
        cl()

hsp_per_read = {}
# directories = os.listdir('results')
# handle = open(os.path.join('results',directories[0]), 'r')
# blast_record = NCBIXML.read(handle)
# best = float('inf')
for result in os.listdir(results_dir):
    handle = open(os.path.join(results_dir, result), 'r')
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
