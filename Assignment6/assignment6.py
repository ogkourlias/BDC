#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -f [INPUT VCF FILE] -r [VCF COMPARISON FILE]
        -ih [INPUT HEADERS TEXT FILE] -ch [COMPARISON HEADER TEXT FILE]
        -chr [CHROMOSOME TEXT FILE] -n [CHUNKSIZE (Variant amount per run)]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import tabix
import gzip
import numpy as np
import time
from scipy.stats import spearmanr, pearsonr


# FUNCTIONS
def arg_parse():
    """
    Function to parse arguments.

    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="VcfCompare")
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("-r", "--reference", type=str)
    parser.add_argument("-ih", "--input_header", type=str)
    parser.add_argument("-ch", "--comp_header", type=str)
    parser.add_argument("-chr", "--chr", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-n", "--chunk", type=int)
    return parser.parse_args()


def open_tbi(input_file):
    """
    Tabix file opener function.
    :input_file: Input VCF file.
    :return: Opened tabix file.
    """
    return tabix.open(input_file)


def compare(chr, start, stop, recs_i, r_handle, headers, chunksize):
    """
    Takes chromosome entry from nextflow process and compare record entries between two files.
    i in var name = VCF File containing all records for which you want to find hits.
    c in var name = vCF File to be compard with.
    :chr: chromosome number in (chr1) format.
    :i_handle: Input VCF file handle.
    :r_handle: Comparison VCF file handle.
    :input_headers_r: Input VCF header file handle.
    :comp_headers_r: Comparison VCF header file handle.
    :re
    turn: Pandas dataframes containing records with intersectioned positions.
    """
    # Open and query tabix files.
    comp_r = open_tbi(r_handle)
    i_records = recs_i
    results = []

    # Query the chunk in the comparison file for a generator.
    r_records = comp_r.query(chr, start, stop)
    # For comparison record in the generator.
    for r_record in r_records:
        # For each retrieved record check whether there's a match in the input record list.
        for i_record in i_records:
            if (
                i_record[1] == r_record[1]  # If CHR is the same.
                and i_record[3] == r_record[3]  # If ref is the same.
                and i_record[4] == r_record[4]  # If alt is the same.
            ):
                results.append(
                    record_extract(i_record, r_record, list(headers.values()))
                )

    # Return resultss as list of strings.
    return results


def header_subset(input_headers_r, comp_headers_r):
    """Looks for header intersections between files and retrieves
    the relevant indexes.

    Args:
        input_headers_r (path): File path for input headers.
        comp_headers_r (path): File path for input headers.

    Returns:
        dictionary: Dictionary of indexes corresponding to intersecting headers between vcf files.
    """

    # Initialise an empty dictionary.
    head_idx = {}

    # For each input header, check whether there is a match in comparison headers.
    for header_i in np.genfromtxt(input_headers_r, delimiter="\t", dtype=str):
        for i, header_r in enumerate(
            np.genfromtxt(comp_headers_r, delimiter="\t", dtype=str)
        ):  # If match is found, save the index in the dictionary.
            if header_i == header_r:
                head_idx[header_i] = i
                break

    return head_idx


def header_subset_gtex(input_headers_r, comp_headers_r):
    """Looks for header intersections between files and retrieves
    the relevant indexes.

    Args:
        input_headers_r (path): File path for input headers.
        comp_headers_r (path): File path for input headers.

    Returns:
        dictionary: Dictionary of indexes corresponding to intersecting headers between vcf files.
    """
    # Initialise an empty dictionary.
    head_idx = {}
    headers = np.genfromtxt(input_headers_r, delimiter="\t", dtype=str)
    # For each input header, check whether there is a match in comparison headers.
    with open(
        "/444",
        "r",
    ) as tbl:
        for line in tbl:
            line_split = line.split(",")
            srr = line_split[0]
            gtex = line_split[28]
            for i, header in enumerate(headers):
                # If match is found, save the index in the dictionary.
                header = header.replace("-splitreads", "")
                if header == srr:
                    headers[i] = gtex

    # For each input header, check whether there is a match in comparison headers.
    for header_i in headers:
        for i, header_r in enumerate(
            np.genfromtxt(comp_headers_r, delimiter="\t", dtype=str)
        ):  # If match is found, save the index in the dictionary.
            if header_i in header_r:
                head_idx[header_i] = i
                break

    return head_idx


def process_file(chr, i_handle, r_handle, headers, output, chunksize=1000):
    """Retrieves the start and stop positions of relevant vcf
    selection.

    Args:
        chr (string): Chromosome string.
        i_handle (path): Input vcf file handle.

    Returns:
        string, string: number strings of start and stop range.
    """
    # Get start and stop positions of chr selection.
    with open(output, "w") as fo:
        # Write headers to output file.
        fo.write(
            "var_id\t"
            + "is_indel\t"
            + "sample_size\t"
            + "nrhoma_i\t"
            + "nrhets_i\t"
            + "nrhomb_i\t"
            + "nrhoma_r\t"
            + "nrhets_r\t"
            + "nrhomb_r\t"
            + "matching\t"
            + "match_ratio\t"
            + "maf_i\t"
            + "maf_r\t"
            + "hwe_i\t"
            + "hwe_r\t"
            + "p_corr\t"
            + "s_corr\n"
        )
        # Open outputfile first for output writing.
        with gzip.open(i_handle, "rt") as f:
            # Create a list with N entries where N = chunksize.
            recs = []
            # Initialise an iteration and index tracker variable.
            for line in f:
                # If line contains relevant variant information.
                if line.startswith(chr):
                    # Split the variant line on tabs.
                    recs.append(line.split("\t"))

                if len(recs) == chunksize:
                    results = compare(
                        chr,
                        int(recs[0][1]) - 1000,
                        int(recs[-1][1]) + 1000,
                        recs,
                        r_handle,
                        headers,
                        chunksize,
                    )
                    recs = []
                    # Write the results to output file.
                    fo.write("".join(results))

            if len(recs) != 0:
                results = compare(
                    chr,
                    int(recs[0][1]) - 1000,
                    int(recs[-1][1]) + 1000,
                    recs,
                    r_handle,
                    headers,
                    chunksize,
                )
                recs = []
                # Write the results to output file.
                fo.write("".join(results))

    return 0


def record_extract(i_record, r_record, head_idx):
    """Extracts the relevant intersection information from two tabix records.

    Args:
        i_record (list): list containing record entries for input file.
        r_record (list): list containing record entries for comparison file.
        head_idx (dict): header indexes

    Returns:
        string: string to be written to output file.
    """
    # Get a column value for each row if header/column index is present in both files.
    r_record = [entry for i, entry in enumerate(r_record) if i in head_idx or i > 8]

    # sample_size and indel variables.
    sample_size = 0
    is_indel = 0

    # Alele count trackers for input.
    nrhoma_i = 0
    nrhets_i = 0
    nrhomb_i = 0

    # Alele count trackers for reference file.
    nrhoma_r = 0
    nrhets_r = 0
    nrhomb_r = 0

    # Genotype lists for use in pearson corr.
    gt_list_r = []
    gt_list_i = []

    # Count of mathcing gts
    matching = 0

    # Construct identifier string for the record/row.
    var_id = f"{i_record[0]}:{i_record[1]}:{i_record[3]}_{i_record[4]}"

    # Selection is same length because of column index selection.
    # Zip to iterate over both record values and track index.
    for i, (entry_i, entry_r) in enumerate(zip(i_record, r_record)):
        match i:  # Match the index value.
            case 0 | 1 | 2 | 5 | 6 | 7:  # Values on these indexes are not relevant.
                continue

            case 3 | 4:  # Check whether ref or alt contains an indel.
                if len(entry_i) > 1:
                    is_indel = 1

            case 8:  # Initialse the format column indicators/headers.
                val_dict_i = {}
                val_dict_r = {}

                # Store format indication/header into a dictionary with empty values.
                for ind_i, ind_r in zip(entry_i.split(":"), entry_r.split(":")):
                    val_dict_i[ind_i] = "./."
                    val_dict_r[ind_r] = "./."

            case other:  # Sample rows.
                # Assign values in format order.
                for key_i, key_r, val_i, val_r in zip(
                    val_dict_i, val_dict_r, entry_i.split(":"), entry_r.split(":")
                ):
                    val_dict_i[key_i] = val_i
                    val_dict_r[key_r] = val_r

                # Save genotypes in seperate var.
                gt_i = val_dict_i["GT"]
                gt_r = val_dict_r["GT"]

                # Split and save genotype numbers for sum.
                # Example: "1/0" becomes list: [1,0]
                vals_i = gt_i.replace("|", "/").split("/")
                vals_r = gt_r.replace("|", "/").split("/")

                # Check whether genotypes are the same between files.
                # Also check whether they are not missing genotypes.
                if vals_i[0] != "." and vals_r[0] != ".":
                    sample_size += 1
                    gt_i = int(vals_i[0]) + int(vals_i[1])
                    gt_r = int(vals_r[0]) + int(vals_r[1])
                    gt_list_i.append(gt_i)
                    gt_list_r.append(gt_r)

                    match (gt_i):
                        case 0:
                            nrhoma_i += 1
                        case 1:
                            nrhets_i += 1
                        case 2:
                            nrhomb_i += 1

                    match (gt_r):
                        case 0:
                            nrhoma_r += 1
                        case 1:
                            nrhets_r += 1
                        case 2:
                            nrhomb_r += 1

                    if gt_r == gt_i:
                        matching += 1

    maf_i = calc_maf(nrhoma_i, nrhets_i, nrhomb_i)
    maf_r = calc_maf(nrhoma_r, nrhets_r, nrhomb_r)

    # Get HWE
    hwe_i = calc_hwe(nrhoma_i, nrhets_i, nrhomb_i)
    hwe_r = calc_hwe(nrhoma_r, nrhets_r, nrhomb_r)

    gt_var_i = np.var(gt_list_i)
    gt_var_r = np.var(gt_list_r)

    if gt_var_i != 0 and gt_var_r != 0 and len(gt_list_i) > 1 and len(gt_list_r) > 1:
        p_corr, p_pval = pearsonr(gt_list_i, gt_list_r)
        s_corr, pval = spearmanr(gt_list_i, gt_list_r)
    else:
        p_corr, s_corr = "nan", "nan"

    if matching != 0 and sample_size != 0:
        match_ratio = matching / sample_size
    else:
        match_ratio = "nan"

    # Return a constructed row string for immediate output writing.
    return (
        f"{var_id}\t{is_indel}\t{sample_size}\t"
        f"{nrhoma_i}\t{nrhets_i}\t{nrhomb_i}\t"
        f"{nrhoma_r}\t{nrhets_r}\t{nrhomb_r}\t"
        f"{matching}\t{match_ratio}\t"
        f"{maf_i}\t{maf_r}\t"
        f"{hwe_i}\t{hwe_r}\t"
        f"{p_corr}\t{s_corr}\n"
    )


def calc_maf(nrhoma, nrhets, nrhomb):
    sample_size = nrhoma + nrhomb + nrhets
    if sample_size == 0:
        return 0
    nrB = (nrhomb * 2) + nrhets
    nrAlleles = sample_size * 2
    MAF = nrB / nrAlleles
    if MAF > 0.5:
        return 1 - MAF  # reference allele is the minor allele
    return MAF  # alternate allele is the minor allele


def calc_hwe(obs_hom1, obs_hets, obs_hom2):
    obs_homc = obs_hom1
    if obs_hom1 < obs_hom2:
        obs_homc = obs_hom2
    obs_homr = obs_hom2
    if obs_hom1 < obs_hom2:
        obs_homr = obs_hom1

    rare_copies = 2 * obs_homr + obs_hets
    l_genotypes = obs_hets + obs_homc + obs_homr

    if l_genotypes == 0:
        return -1

    het_probs = [0] * (rare_copies + 1)
    mid = int(rare_copies) * int(2 * l_genotypes - rare_copies) / int(2 * l_genotypes)
    mid = int(mid)
    if mid % 2 != rare_copies % 2:
        mid += 1
    mid = int(mid)

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    het_probs[int(mid)] = 1.0
    sum = het_probs[int(mid)]

    curr_hets = int(mid)
    while curr_hets > 1:
        het_probs[curr_hets - 2] = (
            het_probs[curr_hets]
            * curr_hets
            * (curr_hets - 1.0)
            / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        )
        sum += het_probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = int(mid)
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    while curr_hets <= (rare_copies - 2):
        het_probs[curr_hets + 2] = (
            het_probs[curr_hets]
            * 4.0
            * curr_homr
            * curr_homc
            / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        sum += het_probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    i = 0
    while i <= rare_copies:
        het_probs[i] /= sum
        i += 1

    p_hwe = 0.0
    i = 0
    while i <= rare_copies:
        if het_probs[i] <= het_probs[obs_hets]:
            p_hwe += het_probs[i]
        i += 1

    if p_hwe > 1:
        p_hwe = 1
    return p_hwe


def calc_hwe_old(nrhoma, nrhet, nrhomb):
    """Calculates HWE deviation using chi-squared test.

    Args:
        nrhoma (int): homozygote for ref count.
        nrhets (_type_): heterozygote count.
        nrhomb (_type_): homozygote for alt count.

    Returns:
        int: chi-square result of HWE deviation.
    """
    # Init vars.
    # AA, Aa, aa.
    dom, het, rec = nrhomb, nrhet, nrhoma
    pop = rec + dom + het

    # Alele freqs.
    p = (2 * dom + het) / (2 * (dom + het + rec))
    q = 1 - p

    # Expected HW.
    exp_dom = p**2 * pop
    exp_het = 2 * p * q * pop
    exp_rec = q**2 * pop

    # Chi-squared for deviation.
    hwe_p = (
        (((dom - exp_dom) ** 2) / exp_dom)
        + (((het - exp_het) ** 2) / exp_het)
        + (((rec - exp_rec) ** 2) / exp_rec)
    )
    return hwe_p


# MAIN
def main(args):
    """Main function"""
    # Get args.
    start = time.time()
    args = arg_parse()
    # Get headers
    headers = header_subset(args.input_header, args.comp_header)
    # Perform comparison and write output.
    process_file(args.chr, args.file, args.reference, headers, args.output, args.chunk)
    end = time.time()
    with open("timings.csv", "a") as f:
        f.write(f"{end - start},{start},{end}\n")
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
