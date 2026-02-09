# TBC
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    for abundance_path in glob(
        "data/16s_ont_251203_lmg-20251204-235001/results/*rel-abundance.tsv"
    ):
        abundance_fname = abundance_path.split("/")[-1]
        print(f"Processing {abundance_fname} ...")

        df_abundance = pd.read_csv(abundance_path, sep="\t")

        df_abundance = df_abundance.sort_values("abundance", ascending=False)
        df_abundance["tax_id"]

        readassmt_path = abundance_path.replace(
            "rel-abundance", "read-assignment-distributions"
        )

        df_reads = pd.read_csv(readassmt_path, sep="\t", header=0)
        colnames_sorted = [
            cn for cn in df_abundance["tax_id"] if cn in df_reads.columns
        ]
        df_reads = df_reads[colnames_sorted]
        df_reads = translate_taxids(df_reads)

        colnames_translated = df_reads.columns

        handles, labels = None, None
        # import ipdb; ipdb.set_trace()
        selected_cols = colnames_translated[0:10]
        n_cols = len(selected_cols)
        subplot_height = 2
        fig, axes = plt.subplots(
            nrows=n_cols,
            figsize=(12, subplot_height * n_cols),
            squeeze=False,
        )
        base_title = (
            abundance_fname.replace("_downsampled.fastq_rel-abundance.tsv", "") + ": "
        )
        for i, colname in enumerate(selected_cols):
            ax = axes[i, 0]

            # Keep only reads which has some abundance at all for this taxa
            df_reads_for_taxa = df_reads[df_reads[str(colname)].notna()]

            sns.histplot(
                df_reads_for_taxa.iloc[:, 0:15],
                bins=10,
                fill=True,
                edgecolor="white",
                linewidth=1,
                ax=ax,
            )
            # sns.kdeplot(
            #        df_reads_for_taxa.iloc[:, 0:10], fill=True
            # )
            ax.set_xlim(0, 1.1)
            title = base_title + colname
            ax.set_title(f"{title}")
            ax.set_xlabel("Assignment probability", fontsize=8)
            ax.set_ylabel("Number of reads per bin", fontsize=8)

            # Capture legend handles and labels from first plot only
            if i == 0:
                handles, labels = ax.get_legend_handles_labels()

            # Remove individual legends from each subplot
            if ax.get_legend() is not None:
                ax.get_legend().remove()

        # Add single legend outside the plots (to the right)
        fig.legend(
            handles,
            labels,
            loc="center right",
            bbox_to_anchor=(1.0, 0.5),
            fontsize=9,
            frameon=True,
        )

        plt.tight_layout(rect=[0, 0, 0.85, 0.99])

        plot_path = f"{readassmt_path}_hist.png"
        plt.savefig(plot_path)
        plt.close()


def translate_taxids(df, taxonomy_path="taxonomy.tsv"):
    original_headers = df.columns.tolist()

    taxdf = pd.read_csv(taxonomy_path, sep="\t", dtype=str).set_index("tax_id")
    taxid_to_label = {
        tax_id: get_best_tax_label(row) for tax_id, row in taxdf.iterrows()
    }

    new_headers = [
        taxid_to_label.get(col.strip(), col) if col.strip().isdigit() else col
        for col in original_headers
    ]
    df.columns = new_headers
    return df


def get_best_tax_label(row):
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
