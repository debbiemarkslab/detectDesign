# Last updated 05/27/2021
# Last modified
#   
# Last added
#   import plotly.graph_objects as go
#   import plotly.express as px
#   import datetime
#   import unittest
#   import os
#
#   def convert_theta(data, data_max)
#   def strand_mod(strand)
#   def viz_xform(score, max_score=8)
#   def detectDesign_pipeline(PAIR_SHARE_STRAND, INCLUDE_TARGET, HAMMING_MAX, 
#       SEED_SIZE, PROXIMAL_FULL, PROXIMAL_SEED, PROXIMAL_MIN, PROXIMAL_MAX, 
#       SMALLSEED_MAX, SMALLSEED_SIZE, DATE, TARGET_CSV, OFF_TARGET_CSV, 
#       TARGET_INPUT_FILE, previous_genome_folder, total_candidate_table, 
#       target_col, full_col, seed_col):

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import seaborn as sns
import pandas as pd
import regex as re
import numpy as np
import datetime
import unittest
import sys
import os

from collections import Counter
from ast import literal_eval
from copy import deepcopy
from glob import iglob


# 1. import genomes ######################################################


def read_seq_gb(gb_file: str) -> str:
    """Returns the sequences in a Genbank file"""
    with open(gb_file, "r") as file:
        gb_seqs = []
        in_origin = False
        for line in file:
            # removes surrounding whitespace
            line = line.replace("\n", "").strip()
            if line == "ORIGIN":
                gb_seq = ""
                in_origin = True

            if line == "//":
                gb_seqs.append(gb_seq)
                in_origin = False

            if in_origin:
                gb_seq += "".join(line.split(" ")[1:])

    return gb_seqs[0]


def read_seq_fa(fa_file: str) -> str:
    """Returns the information from a FASTA file"""
    # Read file
    with open(fa_file, "r") as open_file:
        data = open_file.read()

    # Separate sequence portion of the file
    fa_list = data.split("\n")[1:]
    fa_seq = "".join(fa_list)

    return fa_seq


def read_seq_file(seq_file: str) -> str:
    """Returns the information from a FASTA or genbank file"""
    if seq_file.endswith((".fa", ".FA", ".fasta", ".FASTA")):
        return read_seq_fa(seq_file)

    elif seq_file.endswith((".gb", ".gbk", ".genbank")):
        return read_seq_gb(seq_file)

    else:
        print("Invalid file type, please use '.fa' or '.gb'")
        return ""


def one_hot_seq(sequence: str) -> list:
    """ Converts a sequence from ACTG into one hot format"""
    one_hot_dict = {
        "a": [1, 0, 0, 0],
        "c": [0, 1, 0, 0],
        "g": [0, 0, 1, 0],
        "t": [0, 0, 0, 1],
    }

    one_hot = []
    for bp in sequence:
        one_hot += one_hot_dict[bp.lower()]

    return np.array(one_hot)


def translate_seq_to_regex(sequence_to_translate: str) -> str:
    """Translates a sequence to the regex equivalent"""
    IUPAC_dict = {
        "W": "[AT]",
        "S": "[CG]",
        "M": "[AC]",
        "K": "[GT]",
        "R": "[AG]",
        "Y": "[CT]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACTG]",
    }

    translated_sequence = ""

    # Translate every character if it is within the dictionary
    for base in sequence_to_translate:
        if (base in IUPAC_dict) or (base in ["A", "C", "G", "T"]):
            translated_sequence += IUPAC_dict.get(base.upper(), base.upper())
        else:
            raise ValueError(
                "Sequence contains a character not in the base dictionary!"
            )

    return translated_sequence

def live_write(out, row):
    '''writes row to the end of output table'''
    out.write(','.join([str(r) for r in row]) + '\n')


def convert_to_strand(x):
    ''' Convert 1 and -1 to + and - strand denoters'''
    if x == 1:
        return '+'
    return '-'


# 2. find PAM sites ######################################################


def find_guides(
    genome: str, genome_name: str = "N_A", PAM: str = "NGG"
) -> pd.DataFrame:
    '''Find all guides, the sites within a genome with the specified PAM sequence'''
    PAM_length = len(PAM)
    guide_sites = []
    L = len(genome)

    # Translate PAM to regex
    PAM_regex = translate_seq_to_regex(PAM)

    # Translate sequence to get reverse complement sequence
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    rev_comp = genome.translate(tab)[::-1]

    # Record sites on each strand
    for sequence, strand in [(genome, "+"), (rev_comp, "-")]:
        sites = re.finditer(
            "([ATGC]{20})" + PAM_regex, sequence, overlapped=True, flags=re.IGNORECASE
        )

        for site in sites:
            site_start, site_end = site.span()
            site_seq = site.group()

            if strand == "+":
                i_start = site_start
                i_end = site_end - PAM_length
            if strand == "-":
                i_start = L - site_end + PAM_length - 1
                i_end = L - site_start - 1

            guide_sites.append(
                {
                    "Guide": site_seq[0:20],
                    "PAM": site_seq[-PAM_length:].upper(),
                    "Genome": genome_name,
                    "Strand": strand,
                    "Start": i_start,
                    "End": i_end,
                }
            )

    return pd.DataFrame(guide_sites)


def find_guides_multiple_pams(genome: str, genome_name: str, PAMs: list) -> list:
    """For each of the given PAMs, find the guides in the given genome"""
    
    all_guides = []
    for PAM in PAMs:
        guides = find_guides(genome, genome_name, PAM)
        all_guides.append(guides)
    all_guides = pd.concat(all_guides).drop_duplicates()
    all_guides = all_guides.reset_index(drop=True)
    return all_guides
    	

def PAMs_near_peaks(peaksfile, pamsfile, outfile, window=200, start_col='start'):
    ''' Finds PAMs near peaks given peaks and PAMs files'''

    peaks = pd.read_csv(peaksfile, sep='\t').sort_values(start_col, ascending=True)
    peaks_index = list(peaks[start_col])

    pams = open(pamsfile, 'r')
    out = open(outfile, 'w')
    current_peak = 0
    for i, pam in enumerate(pams):
        if i == 0:
            print(pam.strip('\n').split(','))
            live_write(out, pam.strip('\n').split(',') + list(peaks.columns))
            continue
        N = int(pam.split(',')[5])

        while N > peaks_index[current_peak]+window:
            current_peak += 1
            if current_peak >= len(peaks_index):
                break
        else:
            if (N >= peaks_index[current_peak] - window) and (N <= peaks_index[current_peak] + window):
                #print(pam.split(','))
                #print(list(peaks.iloc[current_peak]))
                live_write(out, pam.strip('\n').split(',') + list(peaks.iloc[current_peak]))
            continue
        break
    pams.close()
    out.close()

    



# 3. find PAM site pairs ######################################################


def all_v_all_distances(sites, start_col="Start", stranded=False, fill=-1):
    """Find the distances between each of the sites in the sites vector"""
    start_mat = np.tile(np.array(sites[start_col]), (len(sites), 1))
    dist_mat = np.abs(start_mat - start_mat.T)

    if stranded:
        strand_mat = np.tile(np.array(sites["Strand"] == "+"), (len(sites), 1))
        strand_mat = strand_mat == strand_mat.T
        dist_mat[~strand_mat] = fill

    return dist_mat


def extract_pairs(pair_coords, sites):
    """Combines pair coordinates and the sites into rows of (pair_1, pair_2)"""
    if len(pair_coords) < 1:
        return []

    i, j = zip(*pair_coords)
    cols = sites.keys()

    sites_i = sites.loc[list(i)].rename(
        columns=dict(zip(cols, [c + "_1" for c in cols]))
    )

    sites_j = sites.loc[list(j)].rename(
        columns=dict(zip(cols, [c + "_2" for c in cols]))
    )

    pairs = pd.concat(
        [sites_i.reset_index(drop=True), sites_j.reset_index(drop=True)], axis=1
    )

    return pairs


def pair_sites(sites, max_dist=200, min_dist=23, pos="Start", gen="Genome"):
    """Checks for sites neighboring one another (between max and min apart),
    and within the same genome or same strand (set gen = 'Strand')
    returns the extracted pairs of the pairs between min and max distances"""
    coords = sites.sort_values(pos)
    pairs = []
    nearby = []
    curr = (0, np.nan, np.nan)

    for i in coords.index:
        nearby.append(curr)
        curr = (i, coords.at[i, pos], coords.at[i, gen])

        for k, vs in enumerate(nearby[::-1]):
            pair_dist = curr[1] - vs[1]
            same_gen = curr[2] == vs[2]

            if pair_dist > max_dist:
                break
            elif (pair_dist >= min_dist) and same_gen:
                pairs.append((curr[0], vs[0]))

        nearby = nearby[-(k + 1) :]

    return extract_pairs(pairs, sites)


def find_target_pairs(target_sites, TARGET_NAME, MAX_DIST, MIN_DIST, TARGET_REGION=False, PAIR_SHARE_STRAND=True, TARGET_CSV=False, OUTPUT_FOLDER='.'):
    '''Find the target pairs with NGG pam, implicit with 'NGG' in target_sites'''
    target_pairs = pair_sites(target_sites, MAX_DIST, MIN_DIST)
    if len(target_pairs) <= 0:
        return []

    target_pairs.loc[:, 'Pair_Dist'] = np.abs(target_pairs['Start_2'] - target_pairs['Start_1'])
    target_pairs.loc[:, 'Shared_Strand'] = (target_pairs['Strand_2'] == target_pairs['Strand_1'])

    dist_cond = target_pairs.Pair_Dist.between(MIN_DIST, MAX_DIST)
    
    candidate_table = target_pairs[(dist_cond) & (target_pairs.Shared_Strand == PAIR_SHARE_STRAND)]
    
    if TARGET_CSV == True:
        if TARGET_REGION == False:
            candidate_table.to_csv(OUTPUT_FOLDER + TARGET_NAME + '.csv')
        else:
            candidate_table.to_csv(OUTPUT_FOLDER + TARGET_NAME + '_' + str(TARGET_REGION[0]) + '_' + str(TARGET_REGION[1]) + '.csv')

    return candidate_table


def split_pairs(target_pairs):
    
    # Stack the guide table
    cols = list(target_pairs.keys())
    cols_1 = [k for k in cols if "_1" in k]
    cols_2 = [k for k in cols if "_2" in k]

    target_list = []

    for i, row in target_pairs.iterrows():
        target_1 = target_pairs.loc[[i], cols_1].rename(
            columns={k: k[:-2] for k in cols_1}
        )
        target_2 = target_pairs.loc[[i], cols_2].rename(
            columns={k: k[:-2] for k in cols_2}
        )

        target_list.append(target_1)
        target_list.append(target_2)

    target_guides_pd = pd.concat(target_list)
    target_guides_pd = target_guides_pd.reset_index()
    return target_guides_pd


def off_target_analysis(
    target_pairs,
    genome,
    hamming_max=8,
    seed_max=3,
    seed_size=8,
    off_target_pams=["NGG", "NAG"],
):
    """ Compute the potential binding sites for each given target"""

    # Load in off-target sites
    sites = find_guides_multiple_pams(genome[0], genome[1], off_target_pams)

    cols = list(target_pairs.keys())
    cols_1 = [k for k in cols if "_1" in k]
    cols_2 = [k for k in cols if "_2" in k]

    offtarget_sites = {}
    offtarget_pairs = {}

    for i, row in target_pairs.iterrows():
        print(i, "computing sites in ", genome[1])
        target_1 = target_pairs.loc[[i], cols_1].rename(
            columns={k: k[:-2] for k in cols_1}
        )
        target_2 = target_pairs.loc[[i], cols_2].rename(
            columns={k: k[:-2] for k in cols_2}
        )
        target_pair = pd.concat([target_1, target_2]).reset_index(drop=True)

        offtarget_sites[i] = find_hamming(
            sites, target_pair, hamming_max, seed_max, seed_size
        )
        offtarget_sites[i] = offtarget_sites[i][
            offtarget_sites[i].Seed_Mism.between(0.0, seed_max)
        ]
        offtarget_sites[i] = offtarget_sites[i][
            offtarget_sites[i].Full_Mism.between(0.0, hamming_max)
        ]
        print(len(offtarget_sites[i]), "site")

        print(i, "computing pairs in ", genome[1])
        offtarget_pairs[i] = pair_sites(offtarget_sites[i])
        print(len(offtarget_pairs[i]), "pairs")

    return offtarget_sites, offtarget_pairs


def find_proximal_perpair(offtargets_pd, MAX_DIST, MIN_DIST):
    ''' Find the proximal off-targets per Pair_idx'''
    total_proximal_sites = []

    for curr_pair_idx in set(offtargets_pd['Pair_idx']):
        curr_offtargets_pd = offtargets_pd[offtargets_pd['Pair_idx'] == curr_pair_idx]
        curr_proximal_sites = pair_sites(curr_offtargets_pd, MAX_DIST, MIN_DIST, gen='Genome')
        if len(curr_proximal_sites) > 0:
            total_proximal_sites.append(curr_proximal_sites)
            
    if len(total_proximal_sites) > 0:
        total_proximal_pd = pd.concat(total_proximal_sites)
    else:
        total_proximal_pd = pd.DataFrame(columns=offtargets_pd.columns)
        
    return total_proximal_pd

def zip_values(spread_pd, col_cut, indices=[]):
    ''' Zip the values, used for prep into side_by_side bar '''
    zipped_list = []
    for i in spread_pd.index:
        if len(indices) == 0:
            zipped_list.append(list(zip(spread_pd.columns[col_cut:], spread_pd.iloc[i][col_cut:])))
        else:
            zipped_list.append(list(zip(indices, spread_pd.loc[i])))
        
    return zipped_list

def expand_list_to_cols(df, col, nan=None, fix_int=False):
    ''' Takes nested list and creates columns'''
    exp_df = deepcopy(df)
    new_cols = pd.DataFrame(list(df[col].apply(dict)))
    new_cols = new_cols.rename(columns = {k: col+str(int(k)) for k in new_cols.keys()})
    if nan is not None:
        new_cols[new_cols.isna()] = nan
    if fix_int:
        new_cols = new_cols.astype(int)
    exp_df = pd.concat([exp_df, new_cols], axis=1)
    return exp_df

def side_by_side_bar(xys, titles=None ,fig=None, ax=None, space=0.15, ticks=None):
    ''' xys should look like [(x1, counts1), (x2, counts2), ..]), [...]..] 
        zipped_data = zip_values(curr_pair_counts_pd[filter_col], 0, split_num)
        side_by_side_bar(zipped_data, curr_pair_counts_pd['Pair_idx'], 
                             ticks=range(int(split_num[-1])+1), space=0.15)'''
    if fig is None:
        fig, ax = plt.subplots()
    N = len(xys)
    offset = (1/N)# * (1+space)
    for i, xy in enumerate(xys):
        x, counts = zip(*xy)
        ax.bar(np.array(x).astype(np.int)+offset*i, counts, width=(1/N)*(1-space))
    if titles is not None:
        ax.legend(titles, loc='center left', bbox_to_anchor=(1, 0.5))
    if ticks is not None:
        ax.set_xticks(ticks)
        



# 4. score sites versus target ######################################################


def find_hamming(
    search_sites: pd.DataFrame,
    guides: pd.DataFrame,
    hamming_max: int,
    seed_max: int,
    seed_size: int = 8,
) -> pd.DataFrame:
    """Find hits in the search sites that match the seed criteria"""
    one_hot_seed_size = 4 * seed_size

    # Add one hot conversion of guide sequence to data frame for testing
    if not "One_Hot" in guides:
        guides.loc[:, "One_Hot"] = guides["Guide"].apply(lambda x: one_hot_seq(x[-20:]))

    if not "One_Hot" in search_sites:
        search_sites.loc[:, "One_Hot"] = search_sites["Guide"].apply(
            lambda x: one_hot_seq(x[-20:])
        )

    if not "Seed_One_Hot" in search_sites:
        search_sites.loc[:, "Seed_One_Hot"] = search_sites["Guide"].apply(
            lambda x: one_hot_seq(x[-seed_size:])
        )

    # Loop through guides and score vs all sites
    scores = []
    for i in guides.index:
        score = deepcopy(search_sites)
        target = guides.at[i, "One_Hot"]
        score.loc[:, "Full_Mism"] = np.count_nonzero(np.stack(score["One_Hot"].to_numpy()) - target, axis=1) / 2
        score.loc[:, "Seed_Mism"] = np.count_nonzero(np.stack(score["Seed_One_Hot"].to_numpy()) - target[-one_hot_seed_size:], axis=1) / 2
        score = score.drop(columns=["One_Hot"])
        scores.append(score.assign(Target_Guide=guides.at[i, "Guide"]))

    scores = pd.concat(scores).reset_index(drop=True)
    scores = scores[scores['Full_Mism'].le(hamming_max)]
    scores = scores[scores['Seed_Mism'].le(seed_max)]

    return scores

def get_mism_counts(off_targets, target_col, full_col, seed_col):
    ''' Get the individual counts for the off-targets provided:
        counts_pd = get_mism_counts(offtargets_total_pd, target_col, full_col, seed_col)
        counts_pd = expand_list_to_cols(counts_pd, 'Seed_Mism', nan=0, fix_int=True)
        counts_pd = expand_list_to_cols(counts_pd, 'Full_Mism', nan=0, fix_int=True)'''
    
    counts_table = []
    for curr_pair_idx in set(off_targets['Pair_idx']):
        curr_pair_table = off_targets[off_targets['Pair_idx'] == curr_pair_idx]
        for curr_guide in set(curr_pair_table[target_col]):
            curr_guide_table = curr_pair_table[curr_pair_table[target_col] == curr_guide]
            for curr_genome in set(curr_guide_table['Genome']):
                table = curr_guide_table[curr_guide_table['Genome'] == curr_genome]
                counts_table.append([curr_pair_idx, curr_genome, curr_guide, Counter(table[full_col]), Counter(table[seed_col])])

    counts_pd = pd.DataFrame(counts_table, columns=['Pair_idx', 'Genome', 'Guide', 'Full_Mism', 'Seed_Mism'])
    return counts_pd

# 5. tabulated outputs and plots ######################################################

def guide_counts(off_targets, target_col='Target_Guide', full_col='Full_Mism', seed_col='Seed_Mism'):
    """Calculate the total and seed compare value counts"""
    counts_table = []
    for guide in set(off_targets[target_col]):
        table = off_targets[off_targets[target_col] == guide]
        counts_table.append([guide, sorted(Counter(table[full_col]).items()), sorted(Counter(table[seed_col]).items())])
    
    return counts_table


def gene_pairs(gene_sites, exact_dist=62, min_dist=23, max_dist=200):
    """Numbers of exact pair targets in a gene and potential 'off-targets'
    with less exact pairing."""
    gene_table = []
    target_table = []

    for k, sites in gene_sites.items():
        target_sites = sites[sites["PAM"].apply(lambda x: x[1:] == "GG")]
        target_pairs = pair_sites(target_sites)
        target_pairs.loc[:, "Distance"] = abs(
            target_pairs["End_1"] - target_pairs["End_2"]
        )
        target_pairs.loc[:, "Shared_Strand"] = (
            target_pairs["Strand_2"] == target_pairs["Strand_1"]
        )

        # Grab target pairs and summarize
        if target_pairs.empty == False:
            target_table.append(target_pairs)
            num_exact = np.sum(target_pairs["Distance"] == exact_dist)
            num_nearby = np.sum(
                (target_pairs["Distance"] >= min_dist)
                & (target_pairs["Distance"] <= max_dist)
            )
        else:
            num_exact = 0
            num_nearby = 0

        gene_table.append(
            {
                "target_gene": k,
                "num_nearby_pairs": num_nearby,
                "num_exact_pairs": num_exact,
            }
        )

    gene_table = pd.DataFrame(gene_table)
    target_table = pd.concat(target_table)

    return [gene_table, target_table]

def convert_theta(data, data_max):
    ''' Convert position in genome to theta values for polar plotting.'''
    return (data/data_max) * 360

def strand_mod(strand):
    ''' Convert strand to positive or negative number 
    for polar plotting reverse strand as negative'''
    if strand == '-':
        return -1
    return 1

def viz_xform(score, max_score=8):
    ''' Modify the score (typically seed match) for visualization in plotly, 
        emphasize closest to max'''
    return 10**(2*score/max_score - 1)


# Default pipeline ##########################

def detectDesign_pipeline(PAIR_SHARE_STRAND, INCLUDE_TARGET, HAMMING_MAX, SEED_SIZE, PROXIMAL_FULL, PROXIMAL_SEED, PROXIMAL_MIN, PROXIMAL_MAX, SMALLSEED_MAX, SMALLSEED_SIZE, DATE, TARGET_CSV, OFF_TARGET_CSV, TARGET_INPUT_FILE, previous_genome_folder, total_candidate_table, target_col, full_col, seed_col):
    ''' Run the whole detectDesign pipeline, creating output folders and files for results
    '''
    target_input = pd.read_csv(TARGET_INPUT_FILE)
    for curr_target_file in set(target_input['Target_Genome']):
        print(curr_target_file)

        # Load in target genome
        target_pd = target_input[target_input['Target_Genome']
                                 == curr_target_file].reset_index()
        TARGET_NAME = target_pd['Genome_Name'][0]
        if curr_target_file.endswith('.csv'):
            target_sites = pd.read_csv(curr_target_file, index_col=0)
        else:
            target_seq = read_seq_file(curr_target_file)
            target_sites = find_guides_multiple_pams(
                target_seq, TARGET_NAME, ['NGG'])

        # Go through target input rows
        for i, row in target_pd.iterrows():
            curr_target_genome, curr_target_name, curr_name, curr_range, curr_offtarget_folder, curr_pair_range = row[
                1], row[2], row[3], row[4], row[5], row[6]
            print(
                curr_target_genome,
                curr_target_name,
                curr_name,
                curr_range,
                curr_offtarget_folder)
            
            # Load in off-target genome(s)

            # Collect all the off-target genome sites
            offtarget_sites_list = [] 
            for curr_file in iglob(curr_offtarget_folder + '*'):
                curr_offtarget_name = curr_file.split('/')[-1].split('.')[-2]
                if curr_file.endswith('.csv') == True:
                    offtarget_sites = pd.read_csv(curr_file, index_col=0)
                else:
                    offtarget_seq = read_seq_file(curr_file)
                    offtarget_sites = find_guides_multiple_pams(offtarget_seq, curr_offtarget_name, ['NGG'])
                offtarget_sites_list.append(offtarget_sites)

            if INCLUDE_TARGET == True:
                offtarget_sites_list.append(target_sites)
            # Combine offtarget genome sites into one file
            offtarget_genome_sites_pd = pd.concat(offtarget_sites_list, sort=True)

            # Save name of loaded genome folder to save computation
            previous_genome_folder = curr_offtarget_folder

            # Set target region names and output folders
            if isinstance(curr_range, str):
                TARGET_REGION = literal_eval(curr_range)
            else:
                TARGET_REGION = [0, len(target_seq)]


            # Set output folder and file names
            OUTPUT_FOLDER = './results/' + curr_target_name + \
                '/' + DATE + '_' + curr_name + '_targets/'
            PREFIX = OUTPUT_FOLDER + TARGET_NAME
            AUX_FOLDER = OUTPUT_FOLDER + '/aux/'
            AUX_PREFIX = AUX_FOLDER + TARGET_NAME
            COUNTS_OUTFILE = AUX_PREFIX + '_offtargets_indv_summary_' + DATE + '.csv'
            COUNTS_PAIR_OUTFILE = AUX_PREFIX + '_offtargets_pair_summary_' + DATE + '.csv'
            COUNTS_PAIR_SUMMARY_OUTFILE = PREFIX + '_offtargets_summary_' + DATE + '.csv'
            PROXIMAL_OUTFILE = AUX_PREFIX + '_offtargets_proximal_fullmism' + \
                str(PROXIMAL_FULL) + '_seedlen' + str(SEED_SIZE) + '_seedmism_' + str(PROXIMAL_SEED) + '_' + DATE + '.csv'
            PROXIMAL_SMALLSEED_OUTFILE = AUX_PREFIX + '_offtargets_smallseed_proximal_fullmism_seedlen' + \
                str(SMALLSEED_SIZE) + '_seedmism_' + str(SMALLSEED_MAX) + '_' + DATE + '.csv'
            FIGURE_OUTFILE = AUX_PREFIX + '_pair_summary_' + DATE + '.png'
            FIGURE_SEED_OUTFILE = AUX_PREFIX + '_seed_pair_summary_' + DATE + '.png'

            # Create output folder if it doesn't exist
            if not os.path.exists(AUX_FOLDER):
                os.makedirs(AUX_FOLDER)

            if isinstance(curr_pair_range, int):
                MIN_DIST, MAX_DIST = [curr_pair_range] * 2 # Turn int into [int, int] for formatting
            else:
                # Select range of distances if exact dist not given
                MIN_DIST, MAX_DIST = literal_eval(curr_pair_range) # Turn '[int, int]' into [int, int]


            # Find all on-target sgRNA pairs
            target_sites_region = target_sites[target_sites['Start'].between(
                TARGET_REGION[0], TARGET_REGION[1])]
            candidate_table = find_target_pairs(
                target_sites_region,
                TARGET_NAME,
                MAX_DIST,
                MIN_DIST,
                TARGET_REGION,
                PAIR_SHARE_STRAND,
                TARGET_CSV,
                AUX_FOLDER)
            # If there are no candidate pairs continue
            if len(candidate_table) == 0:
                print('No candidate pairs in ' + str(curr_name))
                continue


            # Find off-targets for individual sgRNA for all genomes
            split_candidates = split_pairs(candidate_table)
            offtargets_total_pd = find_hamming(
                offtarget_genome_sites_pd,
                split_candidates.drop_duplicates('Guide'),
                HAMMING_MAX,
                SEED_SIZE,
                SEED_SIZE)
            offtargets_total_pd = offtargets_total_pd.merge(
                split_candidates[['Guide', 'index']], left_on='Target_Guide', right_on='Guide')
            offtargets_total_pd = offtargets_total_pd.rename(
                columns={'Guide_x': 'Guide', 'index': 'Pair_idx'})


            # Create summary counts
            # Get individual counts
            counts_pd = get_mism_counts(
                offtargets_total_pd, target_col, full_col, seed_col)

            # Get and sort the counts for on-target paired output and plots
            counts_pd = expand_list_to_cols(
                counts_pd, 'Seed_Mism', nan=0, fix_int=True)
            counts_pd = expand_list_to_cols(
                counts_pd, 'Full_Mism', nan=0, fix_int=True)

            # Find the columns for our given dataframe and sort it by name and
            # number
            n = pd.DataFrame(counts_pd.columns[5:], columns=['Name'])
            n['Pre'] = n['Name'].apply(lambda x: str(x).split('Mism')[0])
            n['Num'] = n['Name'].apply(lambda x: int(str(x).split('Mism')[1]))
            sorted_columns = list(counts_pd.columns[0:3]) + list(
                n.sort_values(['Pre', 'Num'], ascending=[False, True])['Name'])
            counts_pd = counts_pd[sorted_columns]
            counts_pd.to_csv(COUNTS_OUTFILE)

            x = counts_pd.groupby(['Pair_idx', 'Genome']).sum()
            y = counts_pd.groupby(['Pair_idx', 'Genome'])['Guide'].apply(
                lambda x: list(x.drop_duplicates()))
            pair_table_pd = pd.concat([y.apply(lambda x: x[0]).rename('Guide1'),
                                       y.apply(lambda x: x[1]).rename('Guide2'),
                                       x], axis=1).reset_index()
            pair_table_pd = pair_table_pd[sorted_columns[0:2] +
                                          ['Guide1', 'Guide2'] + sorted_columns[3:]]
            pair_table_pd.to_csv(COUNTS_PAIR_OUTFILE)


            # Full_Mism and Seed_Mism paired summary plots
            pair_counts_pd = pair_table_pd.reset_index()  # Reset index to undo grouping

            # Plot Full_Mism and Seed_Mism count plots for on-target pairs in each
            # genome
            for curr_genome in set(pair_counts_pd['Genome']):
                curr_pair_counts_pd = pair_counts_pd[pair_counts_pd['Genome']
                                                     == curr_genome]
                for plot_col in ['Full_Mism', 'Seed_Mism']:
                    filter_col = [
                        col for col in sorted_columns if str(col).startswith(plot_col)]
                    split_num = [x.split('Mism')[-1] for x in filter_col]
                    zipped_data = zip_values(
                        curr_pair_counts_pd[filter_col], 0, split_num)
                    curr_FIG_OUTFILE = AUX_FOLDER + plot_col + \
                        '_summary_' + str(curr_genome) + '_' + DATE + '.png'
                    side_by_side_bar(zipped_data,
                                     curr_pair_counts_pd['Pair_idx'],
                                     ticks=range(int(split_num[-1]) + 1),
                                     space=0.15)
                    plt.title(plot_col + ' in ' + curr_genome)
                    plt.ylabel('Counts')
                    plt.xlabel(plot_col)
                    plt.savefig(curr_FIG_OUTFILE, bbox_inches='tight')
                    plt.close()

            # Find off-targets filtered for pair prep
            total_proximal_pd = find_proximal_perpair(
                offtargets_total_pd[offtargets_total_pd['Seed_Mism'].le(0)], PROXIMAL_MAX, PROXIMAL_MIN)
            cols_for_drop = total_proximal_pd.columns[[('One_Hot' in x) or (
                'Unnamed' in x) or ('_y_' in x) for x in total_proximal_pd.columns]]
            total_proximal_pd = total_proximal_pd.drop(columns=cols_for_drop)
            total_proximal_pd = total_proximal_pd.drop_duplicates()
            total_proximal_pd.to_csv(PROXIMAL_OUTFILE)

            # Create proximal summary counts
            # Get the individual sgRNA proximal counts
            # Get and sort the proximal counts for on-target paired sgRNA output and plots
            # Small seed calculations
            offtargets_smallseed_total_pd = find_hamming(
                offtarget_genome_sites_pd.drop(
                    columns=['Seed_One_Hot']),
                split_candidates.drop_duplicates('Guide'),
                HAMMING_MAX,
                SMALLSEED_MAX,
                SMALLSEED_SIZE)
            offtargets_smallseed_total_pd = offtargets_smallseed_total_pd.merge(
                split_candidates[['Guide', 'index']], left_on='Target_Guide', right_on='Guide')
            offtargets_smallseed_total_pd = offtargets_smallseed_total_pd.rename(
                columns={'Guide_x': 'Guide', 'index': 'Pair_idx'})
            total_smallseed_proximal_pd = find_proximal_perpair(
                offtargets_smallseed_total_pd[offtargets_smallseed_total_pd['Seed_Mism'].le(0)], PROXIMAL_MAX, PROXIMAL_MIN)
            cols_for_drop = total_smallseed_proximal_pd.columns[[('One_Hot' in x) or (
                'Unnamed' in x) or ('_y_' in x) for x in total_smallseed_proximal_pd.columns]]
            total_smallseed_proximal_pd = total_smallseed_proximal_pd.drop(
                columns=cols_for_drop)
            total_smallseed_proximal_pd = total_smallseed_proximal_pd.drop_duplicates()
            total_smallseed_proximal_pd.to_csv(PROXIMAL_SMALLSEED_OUTFILE)

            if len(total_smallseed_proximal_pd) > 0:
                for curr_genome in set(total_smallseed_proximal_pd['Genome_1']):
                    print(curr_name, curr_genome)
                    curr_total_smallseed_proximal_pd = total_smallseed_proximal_pd[
                        total_smallseed_proximal_pd['Genome_1'] == curr_genome]

                    # If any, plot the small seed proximal matches
                    if len(curr_total_smallseed_proximal_pd) > 0:
                        # Fill any missing data with zeros, carry through filled
                        # data
                        s = curr_total_smallseed_proximal_pd.groupby(
                            ['Pair_idx_1'])
                        objects = list(s.groups.keys())
                        data_s_pd = s.count()[['Start_1']].reset_index()
                        fill_zeros_s = set(
                            objects) - set(curr_total_smallseed_proximal_pd['Pair_idx_1'])
                        zeros_s_pd = pd.DataFrame(list(
                            zip(fill_zeros_s, [0] * len(fill_zeros_s))), columns=['Pair_idx_1', 'Start_1'])
                        data_s_pd = data_s_pd.append(
                            zeros_s_pd).sort_values('Pair_idx_1')
                        # Save ordered and filled data
                        data_s = list(data_s_pd['Start_1'])
                        s_pos = np.arange(len(objects))
                        # width:s_pos/2, height:3
                        plt.figure(figsize=(len(s_pos) / 2, 3))
                        plt.bar(s_pos, data_s, align='center', alpha=0.75)
                        name_s = '\n Proximal perfect seed counts per pair in ' + str(curr_genome) + ' for ' + str(curr_name) + '\n Seed ' + str(
                            SMALLSEED_SIZE) + 'bp match: ' + str(SMALLSEED_MAX) + ', Dist: [' + str(PROXIMAL_MIN) + ',' + str(PROXIMAL_MAX) + ']'
                    
                        # If any, plot the large seed proximal matches
                        if len(total_proximal_pd) > 0:
                            curr_total_proximal_pd = total_proximal_pd[total_proximal_pd['Genome_1'] == curr_genome]
                            if len(curr_total_proximal_pd) > 0:
                                # Fill any missing data with zeros, carry through
                                # filled data
                                y = curr_total_proximal_pd.groupby(['Pair_idx_1'])
                                print(y)
                                data_y_pd = y.count()[['Start_1']].reset_index()
                                fill_zeros = set(
                                    objects) - set(curr_total_proximal_pd['Pair_idx_1'])
                                zeros_pd = pd.DataFrame(
                                    list(zip(fill_zeros, [0] * len(fill_zeros))), columns=['Pair_idx_1', 'Start_1'])
                                data_y_pd = data_y_pd.append(
                                    zeros_pd).sort_values('Pair_idx_1')
                                # Save ordered and filled data
                                data_y = list(data_y_pd['Start_1'])
                                # width:s_pos/2, height:3
                                print(data_y)
                                plt.bar(s_pos, data_y, align='center', alpha=0.75)
                            else:
                                print(
                                    'Empty proximal counts ' +
                                    str(SEED_SIZE) +
                                    ' seed in ' +
                                    str(curr_genome))
                        else:
                            print(
                                'Empty proximal counts ' +
                                str(SEED_SIZE) +
                                ' seed in ' +
                                str(curr_genome))
                    else:
                        print(
                            'Empty proximal counts, small seed ' +
                            str(SMALLSEED_SIZE) +
                            ' in ' +
                            str(curr_genome))
                        continue

                    

                    plt.xticks(s_pos, objects)
                    plt.ylabel('Counts')
                    plt.xlabel('Pair idx')
                    name_y = '\n Seed ' + str(SEED_SIZE) + 'bp match: ' + str(
                        PROXIMAL_SEED) + ', Dist: [' + str(PROXIMAL_MIN) + ',' + str(PROXIMAL_MAX) + ']'
                    plt.title(name_s + name_y)
                    plt.legend([str(SMALLSEED_SIZE) + 'bp seed',
                                str(SEED_SIZE) + 'bp seed'])

                    curr_FIG_OUTFILE = OUTPUT_FOLDER + 'proximal_summary_' + \
                        str(curr_genome) + '_' + DATE + '.png'
                    plt.savefig(curr_FIG_OUTFILE, bbox_inches='tight')
                    plt.close()
            else:
                print(
                    'No off-target proximal pairs with ' +
                    str(SMALLSEED_SIZE) +
                    ' seed in ' +
                    str(curr_name))

            if len(total_smallseed_proximal_pd) > 0:
                count_fullseed = 'Prox_' + str(SEED_SIZE) + 'Seed'
                count_smallseed = 'Prox_' + str(SMALLSEED_SIZE) + 'Seed'
                counts_ss = total_smallseed_proximal_pd.groupby(
                    ['Pair_idx_1', 'Genome_1']).count()['Start_1'].reset_index()
                counts_ss.columns = ['Pair_idx', 'Genome', count_smallseed]
                pair_table_pd = pd.merge(
                    pair_table_pd, counts_ss, on=[
                        'Pair_idx', 'Genome'], how='left')
                counts_total_cols = [
                    'Pair_idx',
                    'Genome',
                    'Guide1',
                    'Guide2',
                    count_smallseed,
                    'Seed_Mism0',
                    'Seed_Mism1',
                    'Seed_Mism2']

                if len(total_proximal_pd) > 0:
                    counts = total_proximal_pd.groupby(['Pair_idx_1', 'Genome_1']).count()[
                        'Start_1'].reset_index()
                    counts.columns = ['Pair_idx', 'Genome', count_fullseed]

                    pair_table_pd = pd.merge(
                        pair_table_pd, counts, on=[
                            'Pair_idx', 'Genome'], how='left')

                    counts_total_cols = [
                        'Pair_idx',
                        'Genome',
                        'Guide1',
                        'Guide2',
                        count_smallseed,
                        count_fullseed,
                        'Seed_Mism0',
                        'Seed_Mism1',
                        'Seed_Mism2']

                counts_total_cols = counts_total_cols + \
                    list(n[n['Pre'] == 'Full_'].sort_values(['Pre', 'Num'], ascending=[False, True])['Name'])[0:3]
                pair_table_pd[counts_total_cols[4:]] = pair_table_pd[counts_total_cols[4:]].fillna(
                    0).astype(int)
                pair_table_pd[counts_total_cols].to_csv(
                    COUNTS_PAIR_SUMMARY_OUTFILE)


if __name__ == "__main__":
    raise NotImplementedError("View example notebook for usecase")
