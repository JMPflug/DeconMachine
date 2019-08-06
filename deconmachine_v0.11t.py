#!/usr/bin/env python3

import argparse
import os
import re
import sys
import shutil
import shlex
import subprocess
from collections import Counter
from ete3 import NCBITaxa
from Bio import SeqIO
from Bio import SearchIO
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import plotly
import time
import gzip
import binascii


ncbi = NCBITaxa()


def arguments():
    parser = argparse.ArgumentParser(description=" ")

    # Files

    parser.add_argument("-assembly", "-i", required=True,
                        help="",
                        type=str)

    parser.add_argument("-cov_file", "-c", required=True,
                        help="",
                        type=str)

    parser.add_argument("-ref", "-r", required=True,
                        help="",
                        type=str)

    parser.add_argument("-dmnd_db", "-d", required=False,
                        help="", default=None,
                        type=str)
    # Parameters
    parser.add_argument("-sim_score", "-s", required=True,
                        help="", default=80,
                        type=float)

    parser.add_argument("-min_cov", "-m", required=True,
                        help="", default=10,
                        type=float)

    parser.add_argument("-cov_factor", "-f", required=True,
                        help="", default=0.5,
                        type=float)

    parser.add_argument("-exclude", "-e", default="2",
                        help="NCBI taxon IDs to exclude.",
                        type=str)

    parser.add_argument("-include", "-n", default=None,
                        help="NCBI taxon IDs to include.",
                        type=str)

    parser.add_argument("-threads", "-t", help="Number of threads to use",
                        default=2, type=int)
    # Output

    parser.add_argument("-outfile", "-o", required=True,
                        help="", default="decon_sequences.fa",
                        type=str)

    parser.add_argument("-rejected", "-j", required=False,
                        help="", default="rejected.out.fa",
                        type=str)

    parser.add_argument("-outdir", required=False,
                        help="", default=None,
                        type=str)

    parser.add_argument("-contam_out", default="contam.out.fasta",
                        help="Output FASTA containing sequences identified as originating comtaminate taxa.",
                        type=str)

    args = parser.parse_args()

    if args.cov_factor > 1 and args.cov_factor < 0:
         raise ValueError(
             "The value for coverage factor must be between 0 and 1")

    assembly = args.assembly
    cov_file = args.cov_file
    ref = args.ref
    dmnd_db = args.dmnd_db

    sim_score = args.sim_score
    min_cov = args.min_cov
    cov_factor = args.cov_factor
    include_ids = args.include
    exclude_ids = args.exclude
    threads = args.threads

    outfile = args.outfile
    rejected = args.rejected
    out_dir = args.out_dir
    contam_out = args.contam_out
    if out_dir is None:
        out_dir = "decon_out_{}".format(assembly.split(".")[0])



def args_test_quick():
    global assembly, cov_file, ref, dmnd_db, sim_score, min_cov, cov_factor, taxon, outfile, rejected, out_dir, threads, include_ids, exclude_ids, use_uniref

    assembly = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/dev/assemblies/Siagona-LIB0446-DNA6666-scaffolds_head_A.fa"]
    cov_file = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/dev/coverage_files/Siagona-DNA6666-LIB0446_A.covs"]
    ref = "/Users/MaddisonLab/Documents/JMP/DeconMachine/Asaphidion_Probes_mini.fasta"
    # dmnd_db = "/Users/MaddisonLab/Documents/JMP/DeconMachine/uniref50_plus_taxonomy.dmnd"
    dmnd_db = "/Users/MaddisonLab/Documents/JMP/DeconMachine/swissprot_plus_taxonomy.dmnd"

    use_uniref = False
    sim_score = 80
    min_cov = 10
    cov_factor = 0.2
    include_ids = None
    exclude_ids = 2
    threads = 4

    outfile = "decon.out.fa"
    rejected = "rejected.out.fa"

    contam_out = "contam.out.fasta"


def args_test_quick2():
    global assembly, cov_file, ref, dmnd_db, sim_score, min_cov, cov_factor, taxon, outfile, rejected, out_dir, threads, include_ids, exclude_ids, use_uniref

    assembly = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/SiagonaCLC_head.fas"]
    cov_file = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/SiagonaCLC_head.covstats"]
    ref = "/Users/MaddisonLab/Documents/JMP/DeconMachine/Asaphidion_Probes_mini.fasta"
    # dmnd_db = "/Users/MaddisonLab/Documents/JMP/DeconMachine/uniref50_plus_taxonomy.dmnd"
    dmnd_db = "/Users/MaddisonLab/Documents/JMP/DeconMachine/swissprot_plus_taxonomy.dmnd"

    use_uniref = False
    sim_score = 80
    min_cov = 10
    cov_factor = 0.2
    include_ids = None
    exclude_ids = 2
    threads = 4

    outfile = "decon.out.fa"
    rejected = "rejected.out.fa"

    contam_out = "contam.out.fasta"

def args_test():
    global assembly, cov_file, ref, dmnd_db, sim_score, min_cov, cov_factor, taxon, outfile, rejected, out_dir, threads, include_ids, exclude_ids, use_uniref

    # assembly = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/dev/assemblies/Siagona-LIB0446-scaffolds_head_A.fa"]
    # cov_file = ["/Users/MaddisonLab/Documents/JMP/DeconMachine/dev/coverage_files/Siagona-LIB0446_A.covs"]
    ref = "/Users/maddisonlab/Documents/James/DeconMachine/Asaphidion_Probes_mini.fasta"
    dmnd_db = "/Users/maddisonlab/Documents/James/DeconMachine/uniref50_plus_taxonomy.dmnd"
    # dmnd_db = "/Users/maddisonlab/Documents/James/DeconMachine/swissprot_plus_taxonomy.dmnd"

    sim_score = 80
    min_cov = 10
    cov_factor = 0.2
    include_ids = None
    exclude_ids = 2
    threads = 6
    use_uniref = True

    outfile = "decon.out.fa"
    rejected = "rejected.out.fa"

    contam_out = "contam.out.fasta"

    # assembly, cov_file = read_inputs("/Users/maddisonlab/Documents/James/DeconMachine/test_set_in.txt",
    #             "/Users/maddisonlab/Documents/James/DeconMachine/test_set_covs.txt")
    #
    assembly, cov_file = read_inputs("/Users/maddisonlab/Documents/James/DeconMachine/test_set_in.txt",
                "/Users/maddisonlab/Documents/James/DeconMachine/test_set_in_cov.txt")


def read_inputs(assemblies_in, covs_in):
    assemblies_list = []
    covs_list = []
    with open(assemblies_in) as assemblies:
        for line in assemblies:
            assemblies_list.append(line.strip())
    with open(covs_in) as covs:
        for line in covs:
            covs_list.append(line.strip())

    return assemblies_list, covs_list


def process_cov(cov_file):
    cov_file = cov_file
    df = pd.read_csv(cov_file, sep='\t')
    df.rename(columns={"#ID": "Contig"}, inplace=True)
    df.set_index('Contig', inplace=True)
    df['Status'] = 0
    df.drop(
        ['Ref_GC', 'Covered_bases', 'Covered_percent', 'Plus_reads',
         'Minus_reads', 'Median_fold', 'Read_GC', 'Std_Dev'],
          axis=1, inplace=True)

    return df


def read_fasta(fasta, sanitize=None):
    if sanitize is None:
        sanitize = False
    if fasta is not None:
        print("Loading {}".format(fasta))
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    if sanitize is True:
        for seq in record_dict:
            seq_id = record_dict[seq].id
            if [s for s in ["=", "|", " "] if s in seq_id]:
                raise ValueError(
                    'FASTA file cannot contain spaces or nonstandard characters:\n{}\n{}'.format(fasta ,seq_id))

    return record_dict


def seqs_to_df(df_in, seqs):
    df = df_in.copy()
    for seq in seqs:
        sequence = seqs[seq].seq
        seq_id = seqs[seq].id
        df.loc[[seq_id]]


def min_cov_filter(df_in, min_cov, status):
    df = df_in.copy()
    cov_mean = df['Avg_fold'].mean()
    df.loc[df.Avg_fold <= min_cov, 'Status'] = status

    return df


def prepare_exonerate(out_dir, ref_dict):
    ''' '''
    out_dir = os.path.join(os.getcwd(), out_dir)
    os.makedirs(out_dir, exist_ok=True)

    to_return = []
    for ref in ref_dict:
        ref_seq = ref_dict[ref]
        ref_seq.id = ref_seq.id.split("-")[1]
        ref_dir = os.path.join(out_dir, ref_seq.id)
        os.makedirs(ref_dir, exist_ok=True)
        ref_seq_path = os.path.join(ref_dir, ref_seq.id) + "_reference.fa"
        SeqIO.write(ref_seq, ref_seq_path, "fasta")
        exon_out_name = "{}_exonerate.{}.txt".format(ref_seq.id, "iter1")
        exon_out = os.path.join(ref_dir, exon_out_name)
        out_pair = (ref_seq_path, exon_out)
        to_return.append(out_pair)

    return to_return


def generate_exonerate_cmd(assembly, param_tuple, hits, score_min=None, model=None):
    if score_min == None:
        score_min = "300"
    if score_min == None:
        score_min = "est2genome"
    commands = []
    for param in param_tuple:
        ref_seq = param[0]
        exon_out = param[1]
        ref_seq_path = ref_seq
        command = _generate_exonerate_cmd(
            ref_seq_path, assembly, hits, exon_out, score_min)
        commands.append(command)

    return commands


def _generate_exonerate_cmd(query, target, hits, out_file, score_min, model="est2genome"):
    command = "while true; do exonerate {} {} -m {} -n {} --score {} --extensionthreshold 1000 --ryo \">%qi %ti_TrimLen=%tal_ExScr=%s_Cord=%tab_%tae\\n%ts\" --verbose 0 --showalignment no --showvulgar no > {} && break ; sleep 10; done".format(
        query, target, model, hits, score_min, out_file)
    return command


def run_exonerate(commands, threads, final=False):
    """"""
    if final:
        exon_cmd_sh = "exonerate_final.sh"
    else:
        exon_cmd_sh = "exonerate_round1.sh"
    exonerate_path = os.path.join(out_dir, exon_cmd_sh)
    with open(exonerate_path, 'w') as f:
        for cmd in commands:
            f.write(cmd + "\n")
    parallel_cmd = "parallel --jobs {} < {}".format(threads, exonerate_path)

    #parallel_cmd = shlex.split(parallel_cmd)
    subprocess.run(parallel_cmd, shell=True)

    #subprocess.run([parallel_cmd], shell=True)


def tag_coverage(df_in, seq_dict):
    df = df_in.copy()
    for seq in [seq_dict[seq] for seq in seq_dict]:
        seq_info = df.loc[seq.id]
        seq.id = "{}|Cov={}_Len={}".format(seq.id, seq_info.Avg_fold, seq_info.Length)
        seq.description = ""
    return seq_dict


def write_filtered_fasta(df_in, seq_dict, assembly, out_type, out_dir):
    df = df_in.copy()
    out_name = os.path.splitext(os.path.basename(assembly))[0]
    filtered_out_path = os.path.join(
        out_dir, out_name + "_" + out_type + ".fa")
    df_filtered = df[df.Status == 0]
    seqs_passing = df_filtered.index.tolist()
    seqs_to_write = [
        seq for seq in seq_dict if seq_dict[seq].name in seqs_passing]
    write_dict = {x: seq_dict[x] for x in seqs_to_write}
    write_list = [v for v in write_dict.values()]
    SeqIO.write(write_list, filtered_out_path, "fasta")

    return filtered_out_path


def read_taxon_file(exclude_ids, include_ids, dmnd_out):
    """ """
    ex_tax_list = []
    in_tax_list = []
    # Enhancement: Change this to a tuple with the seq_id and
    # taxon_id_temp so the taxon hit can be specified in final df
    contam_seqs = []
    to_exclude = populate_descendants(exclude_ids, "exc")

    if include_ids is None:
        to_include = None
    else:
        to_include = populate_descendants(include_ids, "inc")

    with open(dmnd_out, 'r') as taxonfile:
        for line in taxonfile:
            line_split = line.split("\t")
            gene_id = line_split[0]
            if line_split[2] == "\n" or line_split[2] == "":
                line_split[2] = "0"
            if ";" not in line_split[2]:
                taxon_id = int(line_split[2])
            else:
                taxon_id_temp = line_split[2].split(";")
                taxon_id = min(list(map(int, taxon_id_temp)))
            if to_include is None:
                if taxon_id in to_exclude:
                    ex_tax_list.append(get_taxon_string(taxon_id))
                    contam_seqs.append(gene_id)
                elif taxon_id not in to_exclude:
                    in_tax_list.append(get_taxon_string(taxon_id))
            else:
                if taxon_id in to_include and taxon_id not in to_exclude:
                    in_tax_list.append(get_taxon_string(taxon_id))
                else:
                    ex_tax_list.append(get_taxon_string(taxon_id))
                    contam_seqs.append(gene_id)
    # print counts of genes removed due to taxonomy filtering
    ex_tax_counter = Counter(ex_tax_list)
    in_tax_counter = Counter(in_tax_list)
    print("Taxonomy screening:")
    print("Excluded {} sequences in {} taxa (top 10 shown):".format(
        len(ex_tax_list), len(ex_tax_counter)))
    for taxon, count in ex_tax_counter.most_common(10):
        print("\t{}\t{}".format(taxon, count))
    print("Included {} sequences in {} taxa (top 10 shown):".format(
        len(in_tax_list), len(in_tax_counter)))
    for taxon, count in in_tax_counter.most_common(10):
        print("\t{}\t{}".format(taxon, count))

    return contam_seqs


def populate_descendants(taxids, mode):
    """Calls ete3 get_descendant_taxa and returns a list containing the taxon
    IDs of all descendant taxa. Used for both including and excluding"""
    desc_tax = []
    if mode == "inc":
        message = "Including: "
    else:
        message = "Excluding: "
    if isinstance(taxids, list):
        temp_taxa = []
        for taxon in taxids:
            temp_taxa = ncbi.get_descendant_taxa(
                taxon, intermediate_nodes=True)
            desc_tax.append(int(taxon))
            desc_tax.extend(temp_taxa)
            print("{}{} {}".format(message, taxon, get_taxon_string(taxon)))
    else:
        desc_tax = ncbi.get_descendant_taxa(taxids, intermediate_nodes=True)
        desc_tax.append(int(taxids))
        print("{}{} {}".format(message, taxids, get_taxon_string(taxids)))

    return set(desc_tax)


def get_taxon_string(taxon_id):
    """Converts ete3 taxon dict to string"""

    taxid2name = ncbi.get_taxid_translator([taxon_id])
    taxkeys = list(taxid2name.keys())
    # Taxon ID 0 causes errors in other functions
    if taxon_id != 0:
        return taxid2name[taxkeys[0]]
    else:
        return "Unknown"


def run_diamond(threads, dmnd_db, assembly, out_dir, outfile="dmnd.out"):
    outfile = os.path.join(out_dir, outfile)
    diamond_cmd = "diamond blastx -p {} -d {} -q {} -o {} -b 0.6 --bin 4\
        --max-target-seqs 1 --outfmt 6 qseqid evalue staxids sseqid".format(
        threads, dmnd_db, assembly, outfile)
    print(diamond_cmd)
    diamond_cmd = shlex.split(diamond_cmd)
    subprocess.run(diamond_cmd, shell=False)
    return outfile


def update_status(df_in, seq_ids, status):
    """Changes status of sequences in a list of
    sequence ID in a dataframe to status"""
    df = df_in.copy()
    for seq in seq_ids:
        if df.loc[seq, 'Status'] == 0:
            df.loc[df.index == seq, 'Status'] = status
        else:
            print("SKIPPING {}".format(seq))

    return df


def add_contam_status(df_in, seq_ids, status):
    """Changes status of sequences in a list of
    sequence ID in a dataframe to status"""
    df = df_in.copy()
    for seq in seq_ids:
        df.loc[df.index == seq, 'Status'] = status

    return df


def get_exonerate_files(path, round_id):
    exonerate_files = []
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            for subentry in os.scandir(entry):
                if round_id in subentry.name:
                    exonerate_files.append(subentry.path)

    return exonerate_files


def parse_exonerate(exonerate_files, assembly):
    print("\tLoading sequences")
    locus_list = []
    locusdf_list = []
    print("Done loading {}".format(assembly))

    for out_file in exonerate_files:
        hits_dict = {}
        if os.stat(out_file).st_size == 0:
            continue
        for seq_record in SeqIO.parse(out_file, "fasta"):
            description = seq_record.description.split(" ")[1]
            cov_match = float(
                re.search("\|Cov=([0-9.]*)_", description).group(1))
            len_match = int(
                re.search("_Len=([0-9]*)_", description).group(1))
            trimlen_match = int(
                re.search("_TrimLen=([0-9]*)_", description).group(1))
            exscore_match = int(
                re.search("_ExScr=([0-9]*)_", description).group(1))

            new_ids = seq_record.description.split(" ")
            contig_id = new_ids[1].split("|")[0]
            new_id = "|".join(new_ids)
            hits_dict = {'SeqID': new_id, 'Contig':contig_id, 'Ref': seq_record.id,
                    'Cov': cov_match, 'Len': len_match, 'TrimLen': trimlen_match, 'ExonerateScore' : exscore_match,
                    'Sequence': seq_record, 'TopCov': False, 'TopLen': False}
            locus_list.append(hits_dict)

    all_hits_df = pd.DataFrame(locus_list)

    df = pd.DataFrame(locus_list)
    df.set_index('SeqID', inplace=True)

    return df, all_hits_df


def df_to_list(df_in, ids_only=None, seq_col="Sequence"):
    df = df_in.copy()
    if ids_only is None:
        ids_only = False
    seq_list = df['Sequence'].tolist()
    if ids_only is True:
        seq_list = [s.description.split()[1].split("|")[0] for s in seq_list]
    return seq_list


def rename_func(seqid, sequence):
    sequence.id = seqid
    return sequence


def seqdf_to_fasta_uniref(df_in, out_dir, assembly, out_tag="hits", seq_col="Sequence"):
    df = df_in.copy()
    out_name = os.path.splitext(os.path.basename(assembly))[0]
    seq_list = [rename_func(y, x) for x, y in zip(df[seq_col], df.index)]

    # >Asaphidion_yukonense-AU0731_Ortho6776 36650_950_28138|cov=63.0979_sim=63.00_len=947
    hitseq_filename = "{}.{}.uniref.fa".format(out_name, out_tag)
    hitseq_path = os.path.join(out_dir, hitseq_filename)
    SeqIO.write(seq_list, hitseq_path, "fasta")
    return hitseq_path


def seqdf_to_fasta(df_in, out_dir, assembly, out_tag="hits", seq_col="Sequence"):
    df = df_in.copy()
    out_name = os.path.splitext(os.path.basename(assembly))[0]
    seq_list = df_to_list(df)
    hitseq_filename = "{}.{}.fa".format(out_name, out_tag)
    hitseq_path = os.path.join(out_dir, hitseq_filename)
    SeqIO.write(seq_list, hitseq_path, "fasta")
    return hitseq_path


def open_gz_fasta(fasta, return_id_dict=None):
    if return_id_dict is None:
        return_id_dict = False
    with open(fasta, 'rb') as test_f:
        is_gzip = binascii.hexlify(test_f.read(2)) == b'1f8b'
    if not return_id_dict:
        if is_gzip:
            with gzip.open(fasta, "rt") as f:
                record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        else:
            record_dict = read_fasta(fasta, sanitize=True)
    if return_id_dict:
        record_dict = {}
        if is_gzip:
            with gzip.open(fasta, "rt") as f:
                for record in SeqIO.parse(f, "fasta"):
                    record_dict[record.id] = record.description
        else:
            for record in SeqIO.parse(fasta, "fasta"):
                record_dict = read_fasta(fasta, sanitize=True)
    return record_dict


def parse_uniprot(dmnd_results, record_dict, out_dir):
    new_dmnd = []
    with open(dmnd_results, "r") as f:
        for line in f:
            dmnd_split = line.strip().split("\t")
            contig_id = dmnd_split[0]
            evalue = dmnd_split[1]
            init_taxid = dmnd_split[2]
            header = dmnd_split[3]
            taxid = taxid_from_uniprot(record_dict[header])
            if taxid != init_taxid:
                print("TaxID Mismatch\nOriginal: {} New: {}".format(init_taxid, taxid))
            new_dmnd_line = "{}\t{}\t{}\t{}".format(contig_id, evalue, taxid, header)
            new_dmnd.append(new_dmnd_line)
    dmnd_out = os.path.join(out_dir, "uniref_dmnd.out")
    with open(dmnd_out, "w") as g:
        for entry in new_dmnd:
            g.write(entry + "\n")
    return dmnd_out


def taxid_from_uniprot(header):
    taxid = re.search(" TaxID=([0-9.]*) ", header).group(1)
    return taxid

def write_final_fasta(df_in, seq_dict, assembly, out_type, out_dir):
    df = df_in.copy()
    out_name = os.path.splitext(os.path.basename(assembly))[0]
    filtered_out_path = os.path.join(
        out_dir, out_name + "_" + out_type + ".fa")
    df_filtered = df[df.Status == 0]
    seqs_passing = df_filtered.index.tolist()
    seqs_to_write = [
        seq for seq in seq_dict if seq_dict[seq].name in seqs_passing]
    write_dict = {x: seq_dict[x] for x in seqs_to_write}
    write_list = [v for v in write_dict.values()]

    SeqIO.write(write_list, filtered_out_path, "fasta")
    return filtered_out_path


def main(assembly, cov_file, uniref_id_dict, min_cov):

    print("Beginning main")
    df = process_cov(cov_file)

    # logger.info(df.head())

    print("Load assembly")
    assembly_dict = read_fasta(assembly, sanitize=True)

    print("Load reference")
    ref_dict = read_fasta(ref)

    print("Prep directories")
    for_exonerate = prepare_exonerate(out_dir, ref_dict)

    print("Min coverage filter")
    cov_df = min_cov_filter(df, min_cov, "VeryLowCov") # Filtering to increse speed of later steps

    print("Tagging coverage")
    assembly_dict = tag_coverage(cov_df, assembly_dict)

    print("Initial filtering")
    filtered_assembly = write_filtered_fasta(cov_df, assembly_dict, assembly, "filtered", out_dir) # Write new fasta without contigs failing minimum filter

    ########## EXONERATE BLOCK ##########
    print("Generating Exonerate commands")
    commands = generate_exonerate_cmd(filtered_assembly, for_exonerate, 10)

    print("Reading filtered assembly")
    assembly_filtered_dict = read_fasta(filtered_assembly)

    print("Running Exonerate")
    run_exonerate(commands, threads)

    print("Getting Exonerate results")
    exonerate_files = get_exonerate_files(out_dir, ".iter1.")

    print("Processing Exonerate results")
    exonerate_results_df, exonerate_TopCov_df = parse_exonerate(exonerate_files, filtered_assembly) # Include all hits for calculating preliminary coverage

    print("Updating database for Exonerate results")
    exonerate_for_update = df_to_list(exonerate_results_df, ids_only = True)

    print("Updating database with Exonerate hits")
    match_df = update_status(cov_df, exonerate_for_update, "Hit")

    ########## TAXON ANNOTATION BLOCK ##########

    print("Taxon Annotation")
    if use_uniref is False:
        print("Writing hit sequences")
        hit_seqs = seqdf_to_fasta(exonerate_results_df, out_dir, assembly)

        dmnd_results = run_diamond(threads, dmnd_db, hit_seqs, out_dir)

    if use_uniref is True:
        hit_seqs_uni = seqdf_to_fasta_uniref(exonerate_results_df, out_dir, assembly)

        init_dmnd_results = run_diamond(threads, dmnd_db, hit_seqs_uni, out_dir)

        # uniref_id_dict = open_gz_fasta("/Users/MaddisonLab/Documents/JMP/DeconMachine/uniref50.fasta.gz", return_id_dict=True)

        dmnd_results = parse_uniprot(init_dmnd_results, uniref_id_dict, out_dir)

    print("Checking for taxonomy matches")
    contam_ids = read_taxon_file([2, 4751], [2759], dmnd_results)

    print("Updating database with clean sequences")
    clean_df = add_contam_status(match_df, contam_ids, "Contam")

    hit_mean_cov = exonerate_results_df['Cov'].mean()
    hit_cov_sdev = exonerate_results_df['Cov'].std()

    hit_mean_cov = clean_df[clean_df.Status == "Hit"]['Avg_fold'].std()
    hit_mean_cov = clean_df[(clean_df.Status == "Hit") | (clean_df.Status != "Contam")]['Avg_fold'].mean()


    # cov_df = min_cov_filter(clean_df[clean_df.Avg_fold > min_cov], hit_mean_cov * hit_cov_sdev * 2, "LowCov")

    print("Updating database with coverage filter")
    cov_df = min_cov_filter(clean_df, hit_mean_cov * cov_factor, "LowCov")
    # cov_df = min_cov_filter(clean_df, min_cov, "LowCov")

    passing_hits = cov_df[cov_df.Status == "Hit"].index.tolist()

    # Get passing_hits from exonerate_results_df, write to file
    passing_seqs = exonerate_results_df[exonerate_results_df['Contig'].isin(passing_hits)].Sequence.tolist()

    pass_name = os.path.splitext(os.path.basename(assembly))[0]

    for seq in passing_seqs:
        seq.id = seq.id.replace("Asaphidion_yukonense-", pass_name + "|")
        # seq.description = seq.description.split()[1]

    # Write single multiFASTA
    pass_path = os.path.join(out_dir, pass_name + "_final.fa")
    SeqIO.write(passing_seqs, pass_path, "fasta")

    # Write individual FASTAs
    print("Writing final sequences")
    single_fastas = os.path.join(out_dir, "final_sequences")

    specimen_name = os.path.splitext(os.path.basename(assembly))[0]
    # if re.search("DNA[0-9]{4}$", assembly) is not None:
    #     dna_number = re.search("(DNA[0-9]{4})", assembly).group(1)
    #     dna_rep = "{}_LIB".format(dna_number)
    #     specimen_name = re.sub(r'LIB', dna_rep, assembly)

    os.makedirs(single_fastas, exist_ok=True)
    for seq in passing_seqs:
        finn = "{}_{}"

        seq_path = os.path.join(single_fastas, "{}.txt".format(seq.id, seq.description))
        # seq.id = "{}_{}|{}".format(specimen_name, seq.id, seq.description.split("|")[0])
        # seq.id = "{}|#|{}|#|{}".format(specimen_name, seq.id, seq.description)

        # seq.description = ''
        SeqIO.write(seq, seq_path, "fasta")


    ########## GRAPHING BLOCK ##########

    print("Drawing graph")

    ex_to_plot = exonerate_results_df[exonerate_results_df['Cov'] > min_cov]

    hit_mean_cov = ex_to_plot['Cov'].mean()
    hit_cov_sdev = ex_to_plot['Cov'].std()

    sdev_lineL = hit_mean_cov - (2 * hit_cov_sdev)
    sdev_lineR = hit_mean_cov - (2 * hit_cov_sdev)
    fig = px.histogram(ex_to_plot, title=os.path.basename(assembly).replace(".fa",""),
            x='Cov').update_layout(shapes=[go.layout.Shape(type="line", ysizemode="scaled",
            x0=exonerate_results_df.Cov.mean(), y0=0, x1=exonerate_results_df.Cov.mean(), y1=1,
            xref="x", yref="paper", line=dict(color="Red", width=1)
            )],

            annotations=[dict(x=hit_mean_cov,
            y=1,
            xref='x',
            yref='paper',
            text="Mean = {:,.0f}".format(hit_mean_cov),
            showarrow=True,
            arrowhead=7,
            ax=1,
            ay=1,
            )])

    fig2 = px.histogram(ex_to_plot, title=os.path.basename(assembly).replace(".fa",""),
            x='Len').update_layout(shapes=[go.layout.Shape(type="line", ysizemode="scaled",
            x0=exonerate_results_df.Len.mean(), y0=0, x1=exonerate_results_df.Cov.mean(), y1=1,
            xref="x", yref="paper", line=dict(color="Red", width=1)
            )])

    fig3 = px.scatter(cov_df, title=os.path.basename(assembly).replace(".fa",""),
               x='Avg_fold', log_x=True, log_y=True, y="Length", marginal_x="box", color="Status",
               color_continuous_scale=px.colors.sequential.Viridis, render_mode="webgl")


    fig_name = "{}.pdf".format(os.path.basename(assembly).replace(".fa",""))
    fig.write_image(os.path.join(out_dir, fig_name))

    fig_name2 = "{}_Len.pdf".format(os.path.basename(assembly).replace(".fa",""))
    fig2.write_image(os.path.join(out_dir, fig_name2))

    fig_name3 = "{}_scatter.pdf".format(os.path.basename(assembly).replace(".fa",""))
    fig3.write_image(os.path.join(out_dir, fig_name3))

    no_micro_assembly = write_final_fasta(clean_df, assembly_dict, assembly, "unscored_hits", out_dir)

    return cov_df, match_df, exonerate_results_df


def write_single_fastas(passing_seqs, assembly):
    single_fastas = os.path.join(out_dir, "scored_sequences")

    specimen_name = os.path.splitext(os.path.basename(assembly))[0]
    # if re.search("DNA[0-9]{4}$", assembly) is not None:
    #     dna_number = re.search("(DNA[0-9]{4})", assembly).group(1)
    #     dna_rep = "{}_LIB".format(dna_number)
    #     specimen_name = re.sub(r'LIB', dna_rep, assembly)

    os.makedirs(single_fastas, exist_ok=True)
    for seq in passing_seqs:
        print("")
        print(seq)
        print("id")
        print(seq.id)
        print(seq.description)
        print("des")
        print("")
        finn = "{}|{}|{}".format(specimen_name, seq.name, seq.id)
        seq.id = finn
        seq_path = os.path.join(single_fastas, "{}.fa".format(finn))
        print(seq_path)
        seq.description = ''
        SeqIO.write(seq, seq_path, "fasta")


def bbb(row):
    seq = row['Sequence']
    score = format(row['Score'], '.4f')
    seq_id = row['Sequence'].id
    seq_des = row['Sequence'].description.replace(" ","_").split("|")[1]
    print("BBBSeqDes = {}".format(seq_des))
    if row['Tag'] is not False:
        tag = row['Tag']
        seq.id = "{}_Score={}_{}".format(seq_des, score, tag, seq_des)
    else:
        # tag = "none"
        # seq.id = "{}_Score={}".format(seq_id, score)
        seq.id = "{}_Score={}".format(seq_des, score)

    return seq


def score(df_in):
    exdf = df_in.copy()
    exdf['Tag'] = False
    exdf['Score'] = 0

    exdf.loc[exdf.groupby(['Ref'])['Cov'].idxmax(), 'TopCov'] = True
    exdf.loc[exdf.groupby(['Ref'])['Len'].idxmax(), 'TopLen'] = True
    exdf.loc[exdf.groupby(['Ref'])['Len'].idxmax(), 'TopLen'] = True

    exdf.loc[(exdf['TopCov'] == True) & (exdf['TopLen'] == True), 'Tag'] = "BOTH"
    exdf.loc[(exdf['TopCov'] == True) & (exdf['TopLen'] == False), 'Tag'] = "HIGHCOV"
    exdf.loc[(exdf['TopCov'] == False) & (exdf['TopLen'] == True), 'Tag'] = "LONGEST"


    exdf['Score'] = (exdf['Cov'] / exdf.groupby('Ref')['Cov'].transform('max'))\
        * (exdf['Len'] / exdf.groupby('Ref')['Len'].transform('max'))

    exdf['Sequence'] = exdf.apply(bbb, axis=1)
    return exdf


# args_test()
# args_test_quick()
args_test_quick2()
# Siagona_test()
# Trachysarus_test()


times = []
if use_uniref:
    print("Loading UniRef50")
    start = time.time()
    uniref_id_dict = open_gz_fasta("uniref50.fasta.gz", return_id_dict=True)
    end = time.time()
    run_time = end - start
    times.append("UniRef loaded in seconds {}".format(run_time))
else:
    uniref_id_dict = None


for i in range(0, len(assembly)):
    out_dir = "decon_out_{}".format(
        os.path.splitext(os.path.basename(assembly[i]))[0])
    print(out_dir)
    start = time.time()
    cov_df, match_df, exonerate_results_df = main(assembly[i], cov_file[i], uniref_id_dict, min_cov)
    exdf = score(exonerate_results_df)
    seqdf_to_fasta(exdf, out_dir, assembly[i], out_tag="scored_hits", seq_col="Sequence")
    scored_seqs = df_to_list(exdf)
    write_single_fastas(scored_seqs, assembly[i])
    end = time.time()
    run_time = end - start
    times.append("{} seconds {}".format(run_time, assembly[i]))

for t in times:
    print(t)









# *$*$*$*$*$*$*$*$*$
# exonerate_TopCov_df # Get only TopCov hit_seqs
# second coverage filtering using hit_mean_cov
# write to file


# TODO:
# Make reference optional
# Use arbitrary coverage cutoff for now
# Plot histogram of coverages for all loci matching references


# Change log
# 0.5
# Moved Diamond annotation to after exonerate to increase speed

# 0.4
# Taxon cleaning now accepts multiple NCBI taxon IDs












#
