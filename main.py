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

        #import ipdb; ipdb.set_trace()
        for colname in colnames_translated:
            # Keep only reads which has some abundance at all for this taxa
            df_reads_for_taxa = df_reads[df_reads[str(colname)].notna()]

            fig = plt.figure(figsize=(8, 5))
            #fig = sns.histplot(
            #        df_reads_for_taxa.iloc[:, 0:15], bins=20, fill=True, edgecolor="white", linewidth=1,
            #)
            fig = sns.kdeplot(
                    df_reads_for_taxa.iloc[:, 0:10], fill=True
            )
            fig.set_xlim(0, 1.1)

            title = abundance_fname.replace("_downsampled.fastq_rel-abundance.tsv", "") + ": " + colname
            plt.title(f"{title}")

            colname_path = colname.replace(" ", "_")
            plot_path = f"{readassmt_path}_{colname_path}.png"
            plt.savefig(plot_path)

            plt.clf()


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
