#!/usr/bin/python
import regex as re
import os
import sys

GENOME = sys.argv[1]

############################################################
# Find all PAM sites
# we use a trick to minimize memory- by searching only
# the current string, and directly outputing the results
# INPUTS:
# 1. Genome fasta (to search)
############################################################


def live_fasta(fafile, OUT, log):
    '''reads a fasta file into dict of strings'''

    # every 100,000 nt update the % completeness
    file_len = len(open(fafile).readlines())
    live_write(['lines to read:', str(file_len)], log)

    cur_name = ''
    cur_seq = ''
    N = 0

    PAMheader = ['cas9.site', 'start', 'end', 'fragment', 'PAM', 'strand']
    live_write([PAMheader], OUT+'_PAMs.tsv')

    with open(fafile, 'r') as file:
        for line in file:
            if line[0] == '>':
                cur_name = line[1:].strip()
                cur_seq = ''
                N = 0
            else:
                cur_seq = cur_seq[-22:] + line.strip()
                # no hit check, for speed
                pams = live_sites(cur_seq, cur_name, N=N)
                live_write(pams, OUT+'_PAMs.tsv')
                N += len(line.strip())

            if (N % 100000) <= 80:
                live_write([100*N/(file_len*len(line.strip())),
                            '% complete'], log)
    return()


# nt to convert from forward to reverse
old_chars = "ACGT"
replace_chars = "TGCA"


def live_sites(genome, name, N=0):
    '''takes a genome, and finds all the potential sites
    for cas9 binding- saving index'''
    L = len(genome)
    PAMs = []

    # Find all sequences with a valid or alternate PAM site
    for PAM in re.finditer('([ATGC]{20})[ATCG][AG]G', genome,
                           overlapped=True, flags=re.IGNORECASE):
        region = (N+PAM.span()[0], N+PAM.span()[1])
        new_pams = (PAM.group(), *region, name, 'N'+PAM.group()[-2:].upper(), '+')
        PAMs.append(new_pams)

    # translate sequence
    tab = str.maketrans(old_chars, replace_chars)
    rev_comp = genome.translate(tab)[::-1]

    for PAM in re.finditer('([ATGC]{20})[ATCG][AG]G', rev_comp,
                           overlapped=True, flags=re.IGNORECASE):
        region = (N+L-PAM.span()[1], N+L-PAM.span()[0])
        new_pams = (PAM.group(), *region, name, 'N'+PAM.group()[-2:].upper(), '-')
        PAMs.append(new_pams)

    return(PAMs)


def live_write(tuple_list, filename):
    '''add rows to end of tsv table'''
    with open(filename, 'a') as file:
        for row in tuple_list:
            file.write('\t'.join([str(r) for r in row])+'\n')


############################################################
# RUN
############################################################

OUT = GENOME.split('/')[-1].split('.')[0]
if not os.path.exists(OUT):
    os.mkdir(OUT)

log = open(OUT+'/log.PAM_search.txt', 'a')

live_write(['making file:', OUT+'/'+'PAMs_'+GENOME.split('/')[-1]], log)

live_fasta(GENOME, OUT+'/'+OUT, log)
