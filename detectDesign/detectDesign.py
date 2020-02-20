# Last updated 02/20/2020
# Last modified
#   Update .get_value() to .at[] for version change
#   From find_target_pairs removed PAIR_EXACT and PAIR_DIST
# Last added
#   import literal_eval
#   convert_to_strand(x)
#   get_mism_counts(off_targets, target_col, full_col, seed_col)
#   side_by_side_bar(xys, titles=None ,fig=None, ax=None, space=0.15, ticks=None)
#   find_proximal_perpair(offtargets_pd, MAX_DIST, MIN_DIST)
#   zip_values(spread_pd, col_cut, indices=[])
#   expand_list_to_cols(df, col, nan=None, fix_int=False)

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import regex as re
import numpy as np
import sys

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


# Default pipeline ##########################

if __name__ == "__main__":
    raise NotImplementedError("View example notebook for usecase")
