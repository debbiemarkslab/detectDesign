#!/usr/bin/python
import pandas as pd
import numpy as np
import glob
import sys
import os

GENOMEFASTA = sys.argv[1]
FEATURETSV = sys.argv[2]
OUTFILE = sys.argv[3]
WINDOW = sys.argv[4:]

############################################################
# Given a list of genome coordinates (tsv: chromosome, position)
# and genome (fasta), returns a table with the seqs receiving mappings.
#
# INPUTS:
# 1. genome (fasta)
# 2. feature table (tsv)
# 3. desired output filename (string)
# 4. extra region to keep before features (int)
# 5. extra region to keep after features (int)
# OUTPUT:
# 1. seq regions mapped to nearby genes (tsv)
############################################################

# retrieving sequences according to features ###############


def get_feature_seqs(fasta, features, outfile, chromosomes, log):
    '''extracts sequence regions from genome, maps to feature info
    it is necessary that the features appear in order (by start, and end)
    at any moment, we need to sequence and coordinate range of the prev feature
    '''
    Lgenome = file_len(fasta)
    Lfeats = len(features)
    live_write(log, [Lgenome, 'lines in genome'])
    live_write(log, [Lfeats, 'features to extract'])
    live_write(log, ['reporting progress every 1,000 features'])

    out = open(outfile, 'a')
    gen = open(fasta, 'r')

    c = ''
    i = 0  # current feature
    j = 1  # current bp index
    coll = ''  # current string collected

    for _, row in features.iterrows():
        s = row['region.start']
        e = row['region.end']

        # if next feature is a new chromosome, skip up to it
        if chromosomes and (row['region.chromosome'] != c):
            c = row['region.chromosome']
            skip_to_chromosome(gen, c)
            coll = ''
            j = 1

        # line-by-line skip up to start, collecting string from there
        # stops once end is passed, keeps full line, gets last index
        # retain previous sequence if start appears within it

        if e > j:
            coll = coll * (s <= j)
            add, j = read_from_to(s, e, gen, j)
            coll = coll + add

        # sorted by start, it's safe to drop the left-hand tail
        # retain right-hand tail in case of overlapping features

        drop, seq, rem = within(s, e, coll, j)
        coll = seq + rem

        live_write(out, list(row.values) + [seq])

        i += 1
        if (i % 1000) == 0:
            live_write(log, [100*i/Lfeats, '% of features in chr. complete'])
            live_write(log, [j, 'bp processed in chr.'])


def read_from_to(a, b, openfile, i):
    '''skips up to line including coord a- collects lines up til b is passed
    keeps whole lines so results in excess before and after, to be trimmed '''
    subseq = ''
    while i < b:
        newline = openfile.readline()
        newline = newline.replace('\n', '').strip()
        i += len(newline)
        if i >= a:
            subseq += newline
    return(subseq, i)


def within(a, b, text, c):
    '''gets region of string *ending at c*, in range [a, b)'''
    l = len(text)
    ar = np.max([0, a - (c - l)])
    br = np.min([0, b - c])

    # seq region starts within string
    if (ar < l):
        return(text[:ar], text[ar:br], text[br:])
    else:
        return('', '', '')


def skip_to_chromosome(openfile, c):
    '''moves up to desired chromosome'''
    newline = openfile.readline()
    while (newline[0] != '>') or (c != get_chromosome(newline).iloc[0]):
        newline = openfile.readline()


def live_write(out, row):
    '''writes row to the end of output table'''
    out.write('\t'.join([str(r) for r in row]) + '\n')


def file_len(filename):
    '''gets length of text file'''
    with open(filename, 'r') as file:
        for l, line in enumerate(file):
            continue
    return(l)


# processing the feature table #############################


def load_features(featuretable, window, log):
    '''reads and sorts the feature table (tsv)
    for extra space around the specified features, specify:
        window = [ - {int before}, + {int after} ]'''

    features = pd.read_csv(featuretable, sep='\t')
    sorton = ['region.start', 'region.end']

    if 'end' in features.keys():
        live_write(log, ['"end" detected, windows will be "start" -> "end"'])
        B = 'end'
    else:
        live_write(log, ['feature windows centering around "start"'])
        B = 'start'

    if 'chromosome' in features.keys():
        features.loc[:, 'region.chromosome'] = get_chromosome(features)
        sorton = ['region.chromosome'] + sorton

    features.loc[:, 'region.start'] = features['start'].apply(
        lambda x: int(x) + window[0])

    features.loc[:, 'region.end'] = features[B].apply(
        lambda x: int(x) + window[1])

    features = features.sort_values(sorton, ascending=True)
    return(features)


def get_chromosome(table):
    '''extract chromosome number as integer, from various common formats'''
    if isinstance(table, str):
        table = pd.DataFrame({'chromosome': [table.replace('\n', '')]})

    v = table['chromosome']
    y = v.iloc[0].lower()

    if isinstance(y, int):
        return(v)

    v = v.apply(lambda x: x.lower().replace(' ', ''))
    v = v.apply(lambda x: x.split(',')[0])

    for c in ['chromosome', 'chr']:
        if c in y:
            v = v.apply(lambda x: x.split(c)[-1])

    v = v.apply(int)
    return(v)

# for genomes split into 1 fasta / chromosome #############################


def find_chromosomes(regex):
    '''given a search string for multiple chromosomes,
    find the corresponding chromosomes from header'''
    cs = {}
    for filename in glob.iglob(regex):

        with open(filename, 'r') as file:
            newline = file.readline()

            while newline[0] != '>':
                newline = file.readline()

            c = get_chromosome(newline).iloc[0]

        cs[c] = filename

    return(cs)


# RUN #############################
log = open('log.' + OUTFILE.split('/')[-1], 'a')

WINDOW = [int(WINDOW[0]), int(WINDOW[1])]
features = load_features(FEATURETSV, WINDOW, log)
header = list(features.keys()) + ['sequence']

if os.path.isfile(OUTFILE):
    live_write(log, [OUTFILE, 'already exists- appending rows to file'])
else:
    with open(OUTFILE, 'a') as out:
        live_write(out, header)

chromosomes = ('chromosome' in features.keys())

if chromosomes and ('*' in GENOMEFASTA):
    clist = find_chromosomes(GENOMEFASTA)
    jobs = list(clist.keys())
    jobs.sort()

    for c in jobs:
        live_write(log, [c, clist[c]])

    for c in jobs:
        live_write(log, ['running chromosome', c, 'file:', clist[c]])
        cfeat = features[features['region.chromosome'].eq(c)]
        get_feature_seqs(clist[c], cfeat, OUTFILE, chromosomes, log)

else:
    get_feature_seqs(GENOMEFASTA, features, OUTFILE, chromosomes, log)
