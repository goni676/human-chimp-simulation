from pathlib import Path

import pandas as pd


def summarize_model(df):
    mutations_per_site = df.groupby("Site_ID").size()
    fixed = df["fixed"].astype(str).str.lower().eq("true")
    return [len(df), len(mutations_per_site), int((mutations_per_site > 1).sum()), int(fixed.sum()), fixed.mean() * 100]


def compare_models(finite_file, infinite_file):
    metrics = ["Mutation events", "Sites", "Sites with multiple mutations", "Mutations fixed in the human sample", "Fixed mutations (%)"]
    finite = pd.read_csv(finite_file)
    infinite = pd.read_csv(infinite_file)
    # The finite table has one row per mutation; the infinite-sites table
    # has one row per site, which contains exactly one mutation in this model.
    return pd.DataFrame({"Metric": metrics, "Finite sites": summarize_model(finite), "Infinite sites": summarize_model(infinite)})


if __name__ == "__main__":
    folder = Path(__file__).resolve().parent
    comparison = compare_models(folder / "mutations_by_site_100.csv", folder / "variants_infinite_100.csv")
    comparison.to_csv(folder / "model_comparison.csv", index=False)
    display = comparison.copy()
    for column in ["Finite sites", "Infinite sites"]:
        display[column] = [f"{value:,.0f}" if i < 4 else f"{value:.3f}%" for i, value in enumerate(comparison[column])]
    print(display.to_string(index=False))
