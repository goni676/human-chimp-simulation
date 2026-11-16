
import msprime
import pandas as pd
import csv
import numpy as np
import sys
import matplotlib.pyplot as plt

#===========================#
# Demographic model - based on the simulations in "Deep coalescent history of the hominin lineage"
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


N = 20_000
SPLIT_TIME = 200_000
ADMIX_TIME = 100_000
SAMPLES = 100

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



# generic mutation example
# ts = msprime.sim_ancestry(100, sequence_length=100, random_seed=1234)
# mts = msprime.sim_mutations(ts, rate=0.01, random_seed=5678)

mts = msprime.sim_mutations(ts, rate=1.25e-08, random_seed=32)

# ts.samples returns an array of the sample *node IDs*
chimp_samples = np.array(ts.samples(population=0))
human_samples = np.array(ts.samples(population=1))



# create mutations DF
i = 0
df_mut = pd.DataFrame(columns=["Site_ID", "Mutation_ID", "fixed", "Inherited_State", "Derived_State", "Mutation_Age", "tMRCA", "Human_Carriers"])
for site in mts.sites():
    tree = mts.at(site.position)
    for mut in site.mutations:
        carriers = set(tree.samples(mut.node))
        for mut_2 in site.mutations:
            if mut_2.parent == mut.id:
                carriers -= set(tree.samples(mut_2.node))
        human_carriers = [int(c) for c in carriers & set(ts.samples(population=1))]
        fixed = len(human_carriers) == len(human_samples)
        if len(human_carriers) == 0:
            tmrca = np.nan
        elif len(human_carriers) == 1:
            tmrca = 0
        else:
            tmrca = find_tmrca(tree, human_carriers)
        df_mut.loc[i] = [site.id, mut.id, fixed, mut.inherited_state, mut.derived_state, mut.time, tmrca, human_carriers]
        i += 1
df_mut.to_csv(f"mutations_by_site_{SAMPLES}.csv", index=False)

# create variants DF
i = 0
df_var = pd.DataFrame(columns=["Site_ID", "Position", "Site_Alleles", "Observed_Alleles_in_Humans", "Human_Genotypes"])
for var in mts.variants(samples=human_samples):
    human_genotypes = var.genotypes
    alleles_set_in_human = set(human_genotypes)
    df_var.loc[i] = [var.site.id, var.site.position, var.alleles, np.array(var.alleles)[sorted(alleles_set_in_human)], human_genotypes]
    i += 1

df_var.to_csv(f"variants_by_site_{SAMPLES}.csv", index=False)


# plot
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

df_mut = pd.read_csv("mutations_by_site_100.csv")


#============Probability of fixation as a function of age/tMRCA=========#
# Keep mutations with human carriers for both fixation plots.
# df_human = df_mut[df_mut["tMRCA"].notna()].copy()
# plot_fixation(df_human, "finite_fixation.png", color="blue")




#plots
#Fixation Probability vs. Mutation Age
# l = df_mut["Mutation_Age"].tolist()
# s = pd.qcut(l, q=20, duplicates="drop")
# d_fixed = {interval: [] for interval in s.categories}
# d_age = {interval: [] for interval in s.categories}

# for i in range(len(l)):
#     d_fixed[s[i]].append(df_mut["fixed"].tolist()[i])
#     d_age[s[i]].append(df_mut["Mutation_Age"].tolist()[i])


# x = []
# y = []

# for interval in d_fixed.keys():
#     x.append(np.mean(d_age[interval]))
#     y.append(np.mean(d_fixed[interval]))

# plt.plot(x, y, marker="o")
# plt.xlabel("Mutation Age (generations)")
# plt.ylabel("Probability of Fixation")
# plt.ylim(0, 1)
# plt.grid(True)
# plt.savefig("Fixation Probability vs. Mutation Age.png", dpi=300, bbox_inches="tight")
# plt.close()


# #Fixation Probability vs. tMRCA
# l = df_mut["tMRCA"].tolist()
# s = pd.qcut(l, q=20, duplicates="drop")
# d_fixed = {interval: [] for interval in s.categories}
# d_tMRCA = {interval: [] for interval in s.categories}

# for i in range(len(l)):
#     d_fixed[s[i]].append(df_mut["fixed"].tolist()[i])
#     d_tMRCA[s[i]].append(df_mut["tMRCA"].tolist()[i])


# x = []
# y = []

# for interval in d_fixed:
#     if len(d_fixed[interval]) > 0:
#         x.append(np.mean(d_tMRCA[interval]))
#         y.append(np.mean(d_fixed[interval]))

# plt.plot(x, y, marker="o")
# plt.xlabel("Carrier tMRCA (generations)")
# plt.ylabel("Probability of Fixation")
# plt.ylim(0, 1)
# plt.grid(True)
# plt.savefig("Fixation Probability vs. tMRCA.png", dpi=300, bbox_inches="tight")
# plt.close()

# print("DONE")