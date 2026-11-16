
import msprime
import pandas as pd
import csv
import numpy as np
import sys
print(sys.executable)
import matplotlib.pyplot as plt

#===========================#
# Demographic model - based on the simulations in Deep coalescent history of the hominin lineage
# - ancestral: 20,000
# 6 Mya (split: 200,000 generations) - ancestral population splits into proto1 and proto2 
# - proto1: 20,000
# - proto2: 20,000
# 3 Mya (admix: 100,000 generations) - human lineage forms as an 80:20% mixture of proto1 and proto2 
#                                    - chimpanzee lineage forms as a 50:50% mixture
# - human: 20,000
# - chimpanzee: 20,000
#
# more parameters:
# - generation time: 30 years
# - mutation rate: 1.25e-08
# - seq length: 2 Gb
#===========================#

def find_tmrca(tree, carriers):
    curr_mrca = carriers[0]
    for carrier in carriers[1:]:
        curr_mrca = tree.mrca(curr_mrca, carrier)
    
    return tree.time(curr_mrca)

# def fixation_check(genotypes):
#     ref = genotypes[0]
#     if ref == 0:
#         return False
#     for g in genotypes:
#         if g != ref:
#             return False
#     return True



N = 20_000
SAMPLES = 100
SPLIT_TIME = 200_000
ADMIX_TIME = 100_000

demography = msprime.Demography()
demography.add_population(name="Chimpanzee", initial_size=N)
demography.add_population(name="Human", initial_size=N)
demography.add_population(name="proto1", initial_size=N)
demography.add_population(name="proto2", initial_size=N)
demography.add_population(name="ANC", initial_size=N)

demography.add_population_split(time=SPLIT_TIME, derived=["proto1", "proto2"], ancestral="ANC")
demography.add_admixture(time=ADMIX_TIME, derived="Human", ancestral=["proto1", "proto2"], proportions=[0.8, 0.2])
demography.add_admixture(time=ADMIX_TIME, derived="Chimpanzee", ancestral=["proto1", "proto2"], proportions=[0.5, 0.5])

demography.sort_events()


#-----Simulate ancestry------#

ts = msprime.sim_ancestry(
    samples={"Chimpanzee": SAMPLES, "Human": SAMPLES},
    demography=demography,
    ploidy=2,
    sequence_length=2_000_000,
    recombination_rate=1.1600523082788902e-08,
    random_seed=90,
    record_migrations=True,
)


# ============infinite sites assumption============= #
# Using discrete_genome=False means that the mutation model will conform 
# to the classic infinite sites assumption (each mutation occurs at a distinct position).
# Adapted from the msprime documentation: https://tskit.dev/msprime/docs/stable/mutations.html

mts = msprime.sim_mutations(ts, rate=1.25e-08, random_seed=32, discrete_genome=False)

#samples returns an array of the sample *node IDs*
chimp_samples = np.array(ts.samples(population=0))
human_samples = np.array(ts.samples(population=1))

i = 0
df_var = pd.DataFrame(columns=["Site_ID", "Position", "Site_Alleles", "Observed_Alleles_in_Humans", "Human_Genotypes", "fixed", "tMRCA", "Mutation_Age"])
for var in mts.variants(samples=human_samples):
    human_genotypes = var.genotypes
    alleles_set_in_human = set(human_genotypes)
    human_carriers = human_samples[human_genotypes == 1]
    fixed = len(human_carriers) == len(human_samples)
    mut = var.site.mutations[0]
    mutation_age = mut.time
    if len(human_carriers) == 0:
        tmrca = np.nan
    else:
        tmrca = find_tmrca(mts.at(var.site.position), human_carriers)
    df_var.loc[i] = [var.site.id, var.site.position, var.alleles, np.array(var.alleles)[sorted(alleles_set_in_human)], human_genotypes, fixed, tmrca, mutation_age]
    i += 1
df_var.to_csv(f"variants_infinite_{SAMPLES}.csv", index=False)


def plot_fixation(df, filename, n_bins=100, color="blue"):
    data = df[["tMRCA", "Mutation_Age", "fixed"]].dropna().copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, x_column, x_label in zip(axes, ["tMRCA", "Mutation_Age"], ["Carrier tMRCA (generations)", "Mutation Age (generations)"]):
        data["bin"] = pd.qcut(data[x_column], q=n_bins, duplicates="drop")
        summary = data.groupby("bin", observed=True).agg(x_mean=(x_column, "mean"), fixation_fraction=("fixed", "mean"))
        ax.plot(summary["x_mean"], summary["fixation_fraction"], marker="o", color=color)
        ax.set_xlabel(x_label)
        ax.set_ylim(0, 1)
        ax.grid(True)
    axes[0].set_ylabel("Fraction fixed in the human sample")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

# Keep mutations with human carriers for both fixation plots.
# df_human_inf = df_var[df_var["tMRCA"].notna()].copy()
# plot_fixation(df_human_inf, "infinite_fixation.png", color="green")

