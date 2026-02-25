"""Plot transcripts.parquet outcomes."""

from pathlib import Path
from typing import Tuple

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd


def count_expression(df, gene_list=["KRT5", "MRC1", "16S"]):
    """Count number of transcripts per gene and rank them."""
    results = []

    for gene in gene_list:
        subset = df[df["feature_name"] == gene]
        total_counts = len(subset)
        high_qv_counts = len(subset[subset["qv"] > 20])
        percent = (high_qv_counts / total_counts) * 100
        results.append(
            {
                "gene": gene,
                "total_counts": total_counts,
                "high_qv_counts": high_qv_counts,
                "percent_high_qv": percent,
                "ROI": subset["ROI"].iloc[0] if not subset.empty else "N/A",
                "Batch": subset["Batch"].iloc[0] if not subset.empty else "N/A",
            }
        )

    summary_df = pd.DataFrame(results)
    summary_df.to_excel("gene_transcript_summary.xlsx", index=False)
    return summary_df


def plot_transcript_ranking(
    df: pd.DataFrame,
    gene_of_interest: str = "16S",
    top_n: int = 10,
    qv_threshold: float = 20,
    out_dir: [Path] = None,
    donor_name: [str] = None,
    figsize: Tuple[float, float] = (8, 5),
) -> pd.DataFrame:
    """Filter transcripts by quality and generate a ranking scatter plot.

    Args:
    df : pd.DataFrame
        Input dataframe containing at least 'feature_name' and 'qv'.
    gene_of_interest : str
        Gene to highlight on the plot.
    top_n : int
        Number of top genes to highlight and label.
    qv_threshold : float
        Minimum qv value to retain transcripts.
    out_dir : Optional[Path]
        Directory to save the figure. If None, figure is not saved.
    donor_name : Optional[str]
        Donor/sample name used for filename if saving.
    figsize : tuple
        Figure size.

    Returns:
    pd.DataFrame
        Ranked gene count table.
    """
    required_cols = {"feature_name", "qv"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Filter high-quality transcripts
    df_filtered = df.loc[df["qv"] > qv_threshold].copy()
    print(f"Transcripts after qv > {qv_threshold} filter: {len(df_filtered):,}")

    # Count transcripts per gene & rank
    gene_counts = (
        df_filtered.groupby("feature_name", observed=True)
        .size()
        .reset_index(name="total_transcripts")
        .sort_values("total_transcripts", ascending=False)
        .reset_index(drop=True)
    )
    gene_counts["rank"] = gene_counts.index + 1

    total_genes = len(gene_counts)
    print(f"Total unique genes: {total_genes:,}")

    if total_genes == 0:
        raise ValueError("No genes remaining after filtering.")

    top_genes = gene_counts.head(top_n)
    gene_row = gene_counts.loc[gene_counts["feature_name"] == gene_of_interest]

    print(
        f"\nTop {top_n} genes:\n"
        f"{top_genes[['rank', 'feature_name', 'total_transcripts']].to_string(index=False)}"
    )

    if not gene_row.empty:
        r = gene_row.iloc[0]
        print(
            f"\n{gene_of_interest}: rank {int(r['rank'])} / {total_genes} | "
            f"{int(r['total_transcripts']):,} transcripts"
        )

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # All genes — grey
    ax.scatter(
        gene_counts["rank"],
        gene_counts["total_transcripts"],
        color="#AAAAAA",
        s=15,
        alpha=0.6,
        zorder=2,
        label="All genes",
    )

    # Top N — blue
    ax.scatter(
        top_genes["rank"],
        top_genes["total_transcripts"],
        color="#0279EE",
        s=60,
        zorder=4,
        label=f"Top {top_n} genes",
    )

    # Gene of interest — orange
    if not gene_row.empty:
        ax.scatter(
            gene_row["rank"],
            gene_row["total_transcripts"],
            color="#FF9400",
            s=80,
            zorder=5,
            label=gene_of_interest,
        )

    # Label top N
    x_offset = total_genes * 0.01
    y_max = gene_counts["total_transcripts"].max()

    for _, row in top_genes.iterrows():
        ax.annotate(
            row["feature_name"],
            xy=(row["rank"], row["total_transcripts"]),
            xytext=(row["rank"] + x_offset, row["total_transcripts"]),
            fontsize=7.5,
            color="#0279EE",
            va="center",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    # Label gene of interest if not in top N
    if not gene_row.empty and gene_row.iloc[0]["rank"] > top_n:
        r = gene_row.iloc[0]
        ax.annotate(
            gene_of_interest,
            xy=(r["rank"], r["total_transcripts"]),
            xytext=(
                r["rank"] + x_offset,
                r["total_transcripts"] + y_max * 0.03,
            ),
            fontsize=8.5,
            color="#FF9400",
            fontweight="bold",
            arrowprops=dict(
                arrowstyle="->",
                color="#FF9400",
                lw=1.2,
            ),
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ax.set_xlabel(
        f"Gene rank (1 = most transcripts, qv > {qv_threshold})",
        fontsize=11,
    )
    ax.set_ylabel("Total transcript count", fontsize=11)
    ax.set_title(
        "Transcript Count Ranking — Xenium 10X",
        fontsize=13,
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()

    if out_dir is not None and donor_name is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / f"{donor_name}_transcript_ranking.png"
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"\nSaved figure to: {output_path}")

    plt.show()

    return gene_counts


# HPC
BASE_DIR = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1"
)
# LOCAL
# BASE_DIR = Path("/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/")

XENIUM_DIR = BASE_DIR / "live/Sara_Patti/009_ST_Xenium/data/xenium_raw"
OUTPATH_DIR = BASE_DIR / "live/Sara_Patti/009_ST_Xenium/data/out/transcripts/"

# Load data
DATA = "/live/Sara_Patti/009_ST_Xenium/data/out/transcripts/xenium_transcripts_combined.csv"
all_transcripts = pd.read_csv(BASE_DIR / DATA)
