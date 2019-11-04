import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import regex as re
import numpy as np
import sys
from glob import iglob
from copy import deepcopy


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


# 2. find PAM sites ######################################################


def find_guides(
    genome: str, genome_name: str = "N_A", PAM: str = "NGG"
) -> pd.DataFrame:
    """Find all guides, the sites within a genome with the specified PAM sequence"""
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
    if len(pair_coords) < 2:
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
        curr = (i, coords[pos].get_value(i), coords[gen].get_value(i))

        for k, vs in enumerate(nearby[::-1]):
            pair_dist = curr[1] - vs[1]
            same_gen = curr[2] == vs[2]

            if pair_dist > max_dist:
                break
            elif (pair_dist >= min_dist) and same_gen:
                pairs.append((curr[0], vs[0]))

        nearby = nearby[-(k + 1) :]

    return extract_pairs(pairs, sites)


def off_target_analysis(
    target_pairs,
    genome,
    hamming_max=8,
    seed_max=3,
    seed_size=8,
    off_target_pams=["NGG", "NAG"],
):
    """ Compute the potential binding sites for each given target"""

    # Load in Subtilis sites
    sites = find_guides_multiple_pams(genome[0], genome[1], off_target_pams)

    cols = list(target_pairs.keys())
    cols_1 = [k for k in cols if "_1" in k]
    cols_2 = [k for k in cols if "_2" in k]

    offtarget_sites = {}
    offtarget_pairs = {}

    for i, row in target_pairs.iterrows():
        print(i, "computing sites")
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

        print(i, "computing pairs")
        offtarget_pairs[i] = pair_sites(offtarget_sites[i])
        print(len(offtarget_pairs[i]), "pairs")

    return offtarget_pairs


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
        target = guides["One_Hot"].get_value(i)
        score.loc[:, "Full_Mism"] = (score["One_Hot"] - target).apply(
            np.count_nonzero
        ) / 2
        score.loc[:, "Full_Mism"] = score["Full_Mism"].apply(np.sum)
        score.loc[:, "Seed_Mism"] = (
            score["Seed_One_Hot"] - target[-one_hot_seed_size:]
        ).apply(np.count_nonzero) / 2
        score.loc[:, "Seed_Mism"] = score["Seed_Mism"].apply(np.sum)
        score = score.drop(columns=["One_Hot"])
        scores.append(score.assign(Target_Guide=guides["Guide"].get_value(i)))
    scores = pd.concat(scores).reset_index(drop=True)
    return scores


# 5. tabulated outputs and plots ######################################################


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
