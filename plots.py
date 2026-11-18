import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

#============Probability of fixation as a function of age/tMRCA=========#


def plot_fixation(df, filename, n_bins=100, color="blue"):
    data = df[["tMRCA", "Mutation_Age", "fixed"]].dropna().copy()
    fig, axes = plt.subplots(2, 1, figsize=(8, 10))
    fig.suptitle(filename)
    for ax in axes:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:_.0f}"))
    for ax, x_column, x_label in zip(axes, ["tMRCA", "Mutation_Age"], ["Carrier tMRCA (generations)", "Mutation Age (generations)"]):
        data["bin"] = pd.qcut(data[x_column], q=n_bins, duplicates="drop")
        summary = data.groupby("bin", observed=True).agg(x_mean=(x_column, "mean"), fixation_fraction=("fixed", "mean"))
        ax.plot(summary["x_mean"], summary["fixation_fraction"], marker="o", color=color)
        ax.set_xlabel(x_label)
        ax.set_ylim(0, 1)
        ax.grid(True)
    axes[0].set_ylabel("Fraction fixed in the human sample")
    for ax in fig.axes:
        line = ax.lines[0]
        x, y = line.get_data()

        for xi, yi in zip(x[-2:], y[-2:]):
            ax.annotate(
                f"{yi:.3f}",
                xy=(xi, yi),
                xytext=(-8, -20),
                textcoords="offset points",
                ha="right",
                fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8),
            )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")




#=========Distribution of age/TMRCA of all fixed mutations=========#


def plot_distribution(df_fixed, name, color):
    fig, axes = plt.subplots(2, 1, figsize=(8, 10))
    fig.suptitle(name)
    for ax in axes:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:_.0f}"))
        
    axes[0].hist(df_fixed["Mutation_Age"].dropna(), bins=30, edgecolor="black", color=color)
    axes[0].set_xlabel("Mutation Age (generations)")
    axes[0].set_ylabel("Number of fixed mutations")
    axes[0].set_title("Age distribution of fixed mutations")

    axes[1].hist(df_fixed["tMRCA"].dropna(), bins=30, edgecolor="black", color=color)
    axes[1].set_xlabel("Carrier tMRCA (generations)")
    axes[1].set_ylabel("Number of fixed mutations")
    axes[1].set_title("tMRCA distribution of fixed mutations")

    
    fig.tight_layout()
    fig.savefig(f"{name}_distributions.png",dpi=300,bbox_inches="tight")

df_mut = pd.read_csv("finite_mutations_1000.csv")
df_mut_inf = pd.read_csv("variants_infinite_1000.csv")
plot_distribution(df_mut[df_mut["fixed"] == True], "finite-sites model", "blue")
plot_distribution(df_mut_inf[df_mut_inf["fixed"] == True], "infinite-sites model", "green")
# #use the same human-carrier mutations for both panels
df_human = df_mut[df_mut["tMRCA"].notna()].copy()
df_human_inf = df_mut_inf[df_mut_inf["tMRCA"].notna()].copy()
plot_fixation(df_human, "finite_fixation_1000.png", color="blue")
plot_fixation(df_human_inf, "infinite_fixation.png", color="green")
