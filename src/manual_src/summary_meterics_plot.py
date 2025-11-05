"""Script to read and visualize the combined metrics summary CSV file."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Setup
base_dir = Path(
    "/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/data/out_data/summary_metrics/"
)

out_dir = base_dir / "figs"
out_dir.mkdir(parents=True, exist_ok=True)

# Load combined dataframe
df = pd.read_excel(base_dir / "csv/combined_run_dfs.xlsx")
print(df.columns)

# Color palettes
condition_color_dict = {"PM08": "#BC3C29", "IPF": "#E18727", "COPD": "#0072B5"}
treatment_color_dict = {"Sham": "#FFDC91", "Treamtment": "#20854E"}

# Global style
sns.set_theme(style="white")  # removes grid lines
plt.rcParams.update(
    {
        "figure.figsize": (6, 6),
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.titlepad": 12,
        "axes.labelpad": 8,
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "figure.autolayout": True,
    }
)

# Data prep
df["Condition"] = pd.Categorical(
    df["Condition"], categories=["PM08", "COPD", "IPF"], ordered=True
)
df["norm_transcript"] = (
    df["total_high_quality_decoded_transcripts"] / df["total_cell_area"]
)
df["norm_cells"] = df["num_cells_detected"] / df["total_cell_area"]


# Helper function
def plot_violin(data, x, y, palette, title, ylabel, filename, order=None, ylim=None):
    """Plot violin plots.

    Args:
        data (_type_): _description_
        x (_type_): _description_
        y (_type_): _description_
        palette (_type_): _description_
        title (_type_): _description_
        ylabel (_type_): _description_
        filename (_type_): _description_
        order (_type_, optional): _description_. Defaults to None.
        ylim (_type_, optional): _description_. Defaults to None.
    """
    plt.figure()
    sns.violinplot(
        x=x,
        y=y,
        data=data,
        inner="point",
        palette=palette,
        order=order,
    )
    plt.title(title)
    plt.xlabel(x)
    plt.ylabel(ylabel)
    if ylim:
        plt.ylim(ylim)
    plt.savefig(filename)
    plt.close()


# Plots
plots = [
    (
        "total_cell_area",
        "Total Cell Area by Condition",
        "Total Cell Area (um²)",
        None,
        "total_cell_area_by_condition.png",
    ),
    (
        "norm_transcript",
        "High Quality Decoded Transcripts Normalized by Cell Area (um²)",
        "Normalized Transcripts (Mean Q >20)",
        None,
        "norm_high_Q_transcript_condition.png",
    ),
    (
        "fraction_transcripts_decoded_q20",
        "Fraction of Transcripts Decoded (Mean Q >20) by Condition",
        "Fraction of Transcripts Decoded (Mean Q >20)",
        (0, 1),
        "fraction_transcripts_decoded_q20_condition.png",
    ),
    (
        "norm_cells",
        "Number of Cells Detected Normalized by Cell Area (um²) by Condition",
        "Normalized Number of Cells Detected (Cells/um²)",
        None,
        "norm_cells_detected_condition.png",
    ),
    (
        "fraction_empty_cells",
        "Fraction of Empty Cells by Condition",
        "Fraction of Empty Cells",
        (-0.001, 0.2),
        "fraction_empty_cells_condition.png",
    ),
    (
        "median_genes_per_cell",
        "Median Genes per Cell by Condition",
        "Median Genes per Cell",
        None,
        "median_genes_per_cell_condition.png",
    ),
    (
        "median_transcripts_per_cell",
        "Median Transcripts per Cell by Condition",
        "Median Transcripts per Cell",
        None,
        "median_transcripts_per_cell_condition.png",
    ),
]

for y, title, ylabel, ylim, filename in plots:
    plot_violin(
        data=df,
        x="Condition",
        y=y,
        palette=condition_color_dict,
        title=title,
        ylabel=ylabel,
        filename=out_dir / filename,
        order=["PM08", "COPD", "IPF"],
        ylim=ylim,
    )

print("Plots saved successfully.")
