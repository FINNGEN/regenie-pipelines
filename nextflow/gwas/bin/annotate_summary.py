#!/usr/bin/env python3
"""Annotate a pheweb-format munged summary stat file with FinnGen variant
annotations (gene, most-severe consequence, rsid, enrichment) via tabix lookup.

Emits two files next to the CWD:
    <pheno>_summary.txt   — variants with pval < summary_pval_thresh
    <pheno>_coding.txt    — variants with pval < coding_pval_thresh and a coding
                             consequence (transcript_ablation, splice_*, stop_*,
                             frameshift_variant, start_lost, inframe_*,
                             missense_variant, protein_altering_variant)

Usage:
    annotate_summary.py <input_file> <finngen_annotation_gz> <summary_pval_thresh> <coding_pval_thresh>
"""
import gzip
import os
import sys
from collections import namedtuple

import pysam


CODING_GROUPS = frozenset([
    "transcript_ablation",
    "splice_donor_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
])

FG_REQ_COLS = [
    "#variant",
    "gene_most_severe",
    "most_severe",
    "rsid",
    "EXOME_enrichment_nfe",
    "GENOME_enrichment_nfe",
]


FGAnnotation = namedtuple(
    "FGAnnotation",
    ["gene", "consequence", "rsid", "EXOME_enrichment_nfe", "GENOME_enrichment_nfe"],
)


def get_header(reader, path):
    with reader(path, "rt") as f:
        return f.readline().strip("\n").split("\t")


def get_fg_annotation(
    iterator,
    variant,
    gene_idx,
    consequence_idx,
    rsid_idx,
    exome_enr_idx,
    genome_enr_idx,
    variant_idx,
):
    for v in iterator:
        cols = v.strip("\n").split("\t")
        if cols[variant_idx] == variant:
            return FGAnnotation(
                cols[gene_idx],
                cols[consequence_idx],
                cols[rsid_idx],
                cols[exome_enr_idx],
                cols[genome_enr_idx],
            )
    return FGAnnotation("", "", "", "", "")


def main(fname, finngen_annotation_file, summary_threshold, coding_threshold):
    pheno = os.path.splitext(os.path.basename(fname))[0]
    output_summary_name = pheno + "_summary.txt"
    output_coding_name = pheno + "_coding.txt"

    fg_tabix = pysam.TabixFile(finngen_annotation_file, parser=None)
    fg_header = get_header(gzip.open, finngen_annotation_file)
    fg_idx = {v: i for (i, v) in enumerate(fg_header)}

    if not all(a in fg_header for a in FG_REQ_COLS):
        raise Exception("Not all columns present in FinnGen annotation! Aborting...")

    var_idx = fg_idx["#variant"]
    gene_idx = fg_idx["gene_most_severe"]
    cons_idx = fg_idx["most_severe"]
    rsid_idx = fg_idx["rsid"]
    exome_enr_idx = fg_idx["EXOME_enrichment_nfe"]
    genome_enr_idx = fg_idx["GENOME_enrichment_nfe"]

    with gzip.open(fname, "rt") as file:
        with open(output_summary_name, "w") as summary_outfile, open(output_coding_name, "w") as coding_outfile:
            header = file.readline().strip("\n").split("\t")
            header_idx = {v: i for (i, v) in enumerate(header)}
            pval_idx = header_idx["pval"]
            cid, pid, rid, aid = (
                header_idx["#chrom"],
                header_idx["pos"],
                header_idx["ref"],
                header_idx["alt"],
            )

            header.extend([
                "rsid",
                "gene_most_severe",
                "most_severe",
                "EXOME_enrichment_nfe",
                "GENOME_enrichment_nfe",
                "phenotype",
            ])
            summary_outfile.write("\t".join(header) + "\n")
            coding_outfile.write("\t".join(header) + "\n")

            for line in file:
                line_columns = line.strip("\n").split("\t")
                pvalue = float(line_columns[pval_idx])
                if pvalue < coding_threshold or pvalue < summary_threshold:
                    cpra = (
                        line_columns[cid],
                        int(float(line_columns[pid])),
                        line_columns[rid],
                        line_columns[aid],
                    )
                    variant = "{}:{}:{}:{}".format(cpra[0], cpra[1], cpra[2], cpra[3])
                    fg_c = (
                        cpra[0]
                        .replace("chr", "")
                        .replace("X", "23")
                        .replace("Y", "24")
                        .replace("XY", "25")
                        .replace("MT", "26")
                        .replace("M", "26")
                    )
                    fg_iter = fg_tabix.fetch(fg_c, cpra[1] - 1, cpra[1])
                    fg_a = get_fg_annotation(
                        fg_iter,
                        variant,
                        gene_idx,
                        cons_idx,
                        rsid_idx,
                        exome_enr_idx,
                        genome_enr_idx,
                        var_idx,
                    )

                    line_columns.extend([
                        fg_a.rsid,
                        fg_a.gene,
                        fg_a.consequence,
                        fg_a.EXOME_enrichment_nfe,
                        fg_a.GENOME_enrichment_nfe,
                        pheno,
                    ])

                    if pvalue < summary_threshold:
                        summary_outfile.write("\t".join(line_columns) + "\n")

                    if pvalue < coding_threshold and fg_a.consequence in CODING_GROUPS:
                        coding_outfile.write("\t".join(line_columns) + "\n")

    print("summary created")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage: annotate_summary.py <input_file> <finngen_annotation_gz> <summary_pval_thresh> <coding_pval_thresh>",
            file=sys.stderr,
        )
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
