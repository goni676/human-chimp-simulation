import pandas as pd
df = pd.read_csv("finite_mutations_1000.csv")
flipping_ids = set()

for site_id, group in df.groupby("Site_ID"):
    mutations = group.to_dict("records")

    for old in mutations:
        for young in mutations:
            if (
                old["Mutation_Age"] > young["Mutation_Age"]
                and old["Inherited_State"] != old["Derived_State"]
                and old["Inherited_State"] == young["Derived_State"]
                and old["Derived_State"] == young["Inherited_State"]
            ):
                flipping_ids.update([
                    old["Mutation_ID"],
                    young["Mutation_ID"]
                ])

df_flipping = df[df["Mutation_ID"].isin(flipping_ids)].copy()
df_flipping = df_flipping.sort_values(["Site_ID", "Mutation_Age"],ascending=[True, False])

df_flipping.to_csv("flipping_mutations.csv", index=False)