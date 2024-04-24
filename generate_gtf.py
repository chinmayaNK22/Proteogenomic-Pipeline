from itertools import islice
import pd_peptidegroup_parser
import read_fasta_file
import argparse
import os

parser = argparse.ArgumentParser(description='''Generates GTF file for peptides identified from Sixframe translated protein sequence database search in Proteome Discoverer''')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='PeptideGroups output from Proteome Discoverer')
parser.add_argument('proteome_fasta', metavar='-f', type=str, nargs='+', help='Proteome database used in the first step (FASTA format)')
parser.add_argument('sixframe_fasta', metavar='-sf', type=str, nargs='+', help='Sixframe translated proteome database used (FASTA format)')

args = parser.parse_args()

##infile = 'M_avium_Mtb_H37RvRa_search_unassigned_6frame_search_071721_PeptideGroups.txt'
##sixframe_fasta = 'sixframe/Mavium_hominissuis_GCF_000007865.1_ASM786v1_genomic_sixframe.fasta'
##fasta = 'Mavium_hominissuis_GCF_000007865.1_ASM786v1_protein.fasta'

def fetch_pos(peptide, fasta):
    chromosome = {}
    for rows in read_fasta_file.read_fasta(os.path.join(fasta)):
        try:
            if rows[0].split(' ')[1] not in chromosome:
                chromosome[rows[0].split(' ')[1]] = [rows[0] + '@' + rows[1]]
            else:
                chromosome[rows[0].split(' ')[1]].append(rows[0] + '@' + rows[1])
        except:
            if rows[0] not in chromosome:
                chromosome[rows[0]] = [rows[0] + '@' + rows[1]]
            else:
                chromosome[rows[0]].append(rows[0] + '@' + rows[1])
            
    for k, v in chromosome.items():
        start_pos = ""
        end_pos = ""
        strand = ""
        for j in v:
            if peptide in j.split('@')[-1]:
                if 'f' in j.split('@')[0].split(' ')[0].split('#')[1]:
                    start_pos = int(j.split('@')[0].split(' ')[0].split('#')[-1].split(':')[0]) + (j.split('@')[1].index(peptide) * 3)
                    end_pos = int(len(peptide)*3) + int(start_pos) - 1
                    strand = '+'
                elif 'f' in j.split('@')[0].split(' ')[0].split('#')[1].split('_')[-1][0]:
                    start_pos = int(j.split('@')[0].split(' ')[0].split('#')[-1].split('-')[0].split('_')[1]) + (j.split('@')[1].index(peptide) * 3)
                    end_pos = int(len(peptide)*3) + int(start_pos) - 1
                    strand = '+'
                if 'r' in j.split('@')[0].split(' ')[0].split('#')[1]:
                    end_pos = int(j.split('@')[0].split(' ')[0].split('#')[-1].split(':')[-1]) - (j.split('@')[1].index(peptide) * 3) + 3
                    start_pos = int(end_pos) - int(len(peptide)*3) + 1
                    strand = '-'
                elif 'r' in j.split('@')[0].split(' ')[0].split('#')[1].split('_')[-1][0]:
                    end_pos = int(j.split('@')[0].split(' ')[0].split('#')[-1].split('-')[-1].split('_')[0]) - (j.split('@')[1].index(peptide) * 3) + 3
                    start_pos = int(end_pos) - int(len(peptide)*3) + 1
                    strand = '-'

        #print (peptide, str(start_pos), str(end_pos), strand, k)
        return peptide, str(start_pos), str(end_pos), strand, k

output = []
def generate_gtf(pep_file, fasta, sixframe_fasta):
    a = pd_peptidegroup_parser.get_seq_idx(os.path.join(pep_file))
    c = 0
    dicts = ''.join(rows[1] for rows in read_fasta_file.read_fasta(os.path.join(fasta)))
    with open(os.path.join(pep_file)) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            peptide = split_i[a].strip('"').split('.')[1]
            if peptide not in dicts:
                pep, start_pos, end_pos, strand, chromosome = fetch_pos(peptide, sixframe_fasta)
                if len(start_pos) > 0 and len(end_pos) > 0:
                    c += 1
                    #print ('NC_002944.2','Custom_script','CDS',  start_pos, end_pos, '.', strand, '0', 'gene ' + '"Mavium_' + str(c) + '";transcript "' + pep + '"')
                    output.append([chromosome,'Custom_script','CDS',  start_pos, end_pos, '.', strand, '0', 'gene ' + '"Mavium_' + str(c) + '";transcript "' + pep + '"'])
                

generate_gtf(args.infile[0], args.proteome_fasta[0], args.sixframe_fasta[0]) 

outfile = "{0}.gtf".format(args.infile[0].rstrip('.txt'))
with open(outfile, 'w') as outf:
    outf.writelines('\t'.join(i) + '\n' for i in output)
