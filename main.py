# TBC
from glob import glob
import argparse
import math
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import os
from pathlib import Path
import seaborn as sns
import statistics

from argparse import ArgumentParser

argp = ArgumentParser()
argp.add_argument("--input-dir", type=str, required=True, help="Path to the input directory containing results")
argp.add_argument("--use-alignment-score", action="store_true", default=False)
argp.add_argument("--plot", action="store_true", default=False)
args = argp.parse_args()


def main():
    taxtr = TaxTranslator()

    for abundance_path_str in glob(
        f"{args.input_dir}/results/*rel-abundance.tsv"
    ):
        abundance_path = Path(abundance_path_str)
        abundance_fname = abundance_path.name
        sample_name = (
            abundance_fname.replace("_downsampled.fastq_rel-abundance.tsv", "")
        )
        readassmt_path = abundance_path_str.replace(
            "rel-abundance", "read-assignment-distributions"
        )
        alignments_pathobj = Path(
            abundance_path_str.replace("_rel-abundance.tsv", "_emu_alignments.sam")
        )
        results_dir = abundance_path.parent

        #print(
        #    "--------------------------------------------------------------------------------"
        #)
        #print(f"Processing {sample_name} ...")
        #print(
        #    "--------------------------------------------------------------------------------"
        #)

        df_abundance_unsorted = pd.read_csv(abundance_path_str, sep="\t")
        df_abundance = df_abundance_unsorted.sort_values("abundance", ascending=False)
        #import ipdb; ipdb.set_trace()

        df_reads = load_read_file(readassmt_path, df_abundance)
        colnames_taxids = df_reads.columns

        # ================================================
        # Start of alignment score calculation
        # ================================================
        if args.use_alignment_score:
            align_file = pysam.AlignmentFile(str(alignments_pathobj))

            alns_all = {}
            for aln in align_file:
                readid = aln.query_name  # Ex: b9bb144e-eb53-4509-8931-5f4477444a48
                refname = aln.reference_name  # Ex: 562:emu_db:23853
                taxid = str(aln.reference_name).split(":")[0]  # Ex: 562
                refid = aln.reference_id  # Ex: 23853

                if aln.is_secondary or aln.is_supplementary:
                    # We don't count these
                    continue

                if taxid not in alns_all:
                    alns_all[taxid] = []
                alns_all[taxid].append(aln)

            aln_infos = []
            i = 1
            for taxid in colnames_taxids:
                if taxid in alns_all:
                    alns = alns_all[taxid]
                    alns_cnt = len(alns)
                    identities, coverages = collect_distribution(alns)
                    median_id = statistics.median(identities)
                    median_cov = statistics.median(coverages)
                    abundance = float(df_abundance[df_abundance["tax_id"] == taxid]["abundance"].values[0])
                    taxon = taxtr.taxid_to_label(taxid)
                    aln_infos.append(
                        {
                            "sample": sample_name,
                            "abundance": abundance,
                            "taxid": taxid,
                            "taxon": taxon,
                            "aligned_reads": alns_cnt,
                            "median_identity": median_id,
                            "median_coverage": median_cov,
                            "identities": identities,
                            "coverages": coverages,
                        }
                    )
                    if args.plot:
                        _, axes = plt.subplots(
                            nrows=2,
                            figsize=(12, 7),
                            squeeze=False,
                        )
                        ax = sns.histplot(
                            data=identities,
                            bins=100,
                            binrange=(0, 1),
                            fill=True,
                            edgecolor="white",
                            linewidth=1,
                            legend=True,
                            ax=axes[0,0]
                        )
                        title = f"Med %Id | Sp. {i}: {sample_name} - {taxon} ({abundance:.3f})"
                        ax.set_title(title)
                        ax = sns.histplot(
                            data=coverages,
                            bins=100,
                            binrange=(0, 1),
                            fill=True,
                            edgecolor="white",
                            linewidth=1,
                            legend=True,
                            ax=axes[1,0]
                        )
                        title = f"Med. cov | Sp. {i}: {sample_name} - {taxon} ({abundance:.3f})"
                        ax.set_title(title)

                        outdir=f"plots/{sample_name}"
                        os.makedirs(outdir, exist_ok=True)

                        png_path = f"{outdir}/{sample_name}-{i:03d}-{taxon.replace(' ', '_')}.png"
                        print(f"Saving raw alignment plot to {png_path} ...")
                        plt.savefig(png_path)

                        i += 1
                        plt.close()

            for ai in aln_infos:
                print(f"{ai['sample']}\t{ai['taxid']}\t{ai['abundance']:.5f}\t{ai['aligned_reads']}\t{ai['median_identity']:.3f}\t{ai['median_coverage']:.3f}\t{ai['taxon']}")

            continue
        # ================================================
        # End of alignment score calculation
        # ================================================

        handles, labels = None, None
        selected_cols = colnames_taxids[0:10]
        n_cols = len(selected_cols)
        subplot_height = 1.37
        fig, axes = plt.subplots(
            nrows=n_cols,
            figsize=(12, subplot_height * n_cols),
            squeeze=False,
        )
        for i, colname in enumerate(selected_cols):
            ax = axes[i, 0]

            # Keep only reads which has some abundance at all for this taxa
            df_reads_for_taxa = df_reads[df_reads[colname].notna()]

            for col in df_reads_for_taxa[selected_cols].columns:
                sns.histplot(
                    data=df_reads_for_taxa[selected_cols],
                    x=col,
                    bins=10,
                    binrange=(0, 1),
                    fill=True,
                    edgecolor="white",
                    linewidth=1,
                    ax=ax,
                    label=col,
                    alpha=0.67,
                    legend=True,
                )
            # sns.kdeplot(
            #        df_reads_for_taxa.iloc[:, 0:10], fill=True
            # )
            ax.set_xlim(0, 1.1)
            title = colname
            ax.set_title(f"{title}")
            ax.set_xlabel("Assignment probability", fontsize=8)
            ax.set_ylabel("Number of reads per bin", fontsize=8)

            # Capture legend handles and labels from first plot only
            if i == 0:
                handles, labels = ax.get_legend_handles_labels()

            # Remove individual legends from each subplot
            legend = ax.get_legend()
            if legend is not None:
                legend.remove()

        # Add single legend outside the plots (to the right)
        fig.legend(
            handles,
            labels,
            loc="upper right",
            bbox_to_anchor=(1.0, 0.5),
            fontsize=9,
            frameon=True,
        )

        plt.tight_layout(rect=[0, 0, 0.7, 0.99])

        png_path = f"{readassmt_path}_hist.png"
        print(f"Saving figure to {png_path} ...")
        plt.savefig(png_path)

        plt.close()


def collect_distribution(alns):
    identities = []
    coverages = []
    for aln in alns:
        identity, coverage = get_align_stats(aln)
        identities.append(identity)
        coverages.append(coverage)
    return identities, coverages


def load_read_file(readassmt_path, df_abundance):
    df_reads = pd.read_csv(readassmt_path, sep="\t", header=0)
    colnames_sorted = [cn for cn in df_abundance["tax_id"] if cn in df_reads.columns]
    df_reads = df_reads[colnames_sorted]
    return df_reads


def get_align_stats(alignment):
    """
    Return list of inquired cigar stats (I,D,S,X) for alignment
    """
    cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]

    insertions = cigar_stats[1]
    deletions = cigar_stats[2]
    soft_clips = cigar_stats[4]

    query_len = alignment.query_length
    aln_len = alignment.query_alignment_length

    nm = alignment.get_tag("NM") if alignment.has_tag("NM") else None

    matches = 0
    mismatches = 0
    if nm is not None:
        mismatches = nm - insertions - deletions
        matches = aln_len - insertions - mismatches

    # We can not normalize over query or reference length, as these
    # won't contain either insertions or deletions
    divisor = matches + mismatches + insertions + deletions

    identity = 0
    if divisor > 0:
        identity = matches / divisor

    coverage = 0
    if query_len > 0:
        coverage = aln_len / query_len

    return identity, coverage


class TaxTranslator(object):
    def __init__(self, taxonomy_path="taxonomy.tsv"):
        self.taxdf = pd.read_csv(taxonomy_path, sep="\t", dtype=str).set_index("tax_id")
        self.taxid_to_label_mapping = {
            tax_id: self.get_best_tax_label(row) for tax_id, row in self.taxdf.iterrows()
        }

    def taxid_to_label(self, taxid):
        if taxid in self.taxid_to_label_mapping:
            return self.taxid_to_label_mapping[taxid]
        return taxid

    def translate_taxids_in_df_columns(self, df):
        df_cols_orig = df.columns.tolist()
        new_headers = [
            self.taxid_to_label(col.strip()) if col.strip().isdigit() else col
            for col in df_cols_orig
        ]
        df.columns = new_headers
        return df

    def get_best_tax_label(self, row):
        """Return the best available taxonomic label from left to right."""
        for level in [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "clade",
            "superkingdom",
            "subspecies",
            "species subgroup",
            "species group",
        ]:
            val = row.get(level, "")
            if pd.notna(val) and str(val).strip() != "":
                return val.strip()
        return "Unknown"


if __name__ == "__main__":
    main()
