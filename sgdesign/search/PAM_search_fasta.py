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


def live_fasta(fafile, out, log):
    '''reads a fasta file into dict of strings'''
    log = open(log, 'a')
    out = open(out, 'a')

    # every 100,000 nt update the % completeness
    file_len = len(open(fafile).readlines())
    live_write(log, ['lines to read:', str(file_len)])

    cur_name = ''
    cur_seq = ''
    N = 0

    PAMheader = ['cas9.site', 'start', 'end', 'fragment', 'PAM', 'strand']
    live_write(out, PAMheader)

    file = open(fafile, 'r')
    for line in file:
        if line[0] == '>':
            cur_name = line[1:].strip()
            cur_seq = ''
            N = 0
        else:
            cur_seq = cur_seq[-22:] + line.strip()
            # no hit check, for speed
            pams = live_sites(cur_seq, cur_name, N=N)

            for pam in pams:
                live_write(out, pam)

            N += len(line.strip())

        if (N % 100000) <= 80:
            live_write([100*N/(file_len*len(line.strip())), '% complete'], log)
            log.flush()

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


def live_write(out, row):
    '''writes row to the end of output table'''
    out.write('\t'.join([str(r) for r in row]) + '\n')


############################################################
# RUN
############################################################

out = GENOME.split('/')[-1].split('.')[0]
if not os.path.exists(out):
    os.mkdir(out)

log = out + '/' + '/log.PAM_search.txt'
out = out + '/' + 'PAMS_' + out + '.tsv'

live_fasta(GENOME, out, log)
