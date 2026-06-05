import re

import matplotlib.pyplot as plt
import pandas as pd

tables = [
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-3/ten-images/plots/avg_img_CV_250x3500x500_bin1x1_125_10_stitched_3deadef8/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-4/plots/avg_img_CV_250x3500x500_bin1x1_125_10_stitched_ad841c6d/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-5/plots/avg_img_CV_250x3500x500_bin1x1_125_10_stitched_213b37c2/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-6/plots/avg_img_CV_250x3500x500_bin1x1_125_10_stitched_dbec0ad2/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-7/plots/cds_avg_img_L2_250x3500x500_1x1_L2_124_12_stitched_0392c058/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-8/plots/cds_avg_img_L2_250x3500x500_1x1_L2_125_10_stitched_836d8173/extension_summary.csv",
    "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-9/plots/cds_avg_img_L2_250x3500x500_1x1_L2_125_10_stitched_30fff068/extension_summary.csv",
]

# Quantities to plot: CSV column name -> (axis label, plot title)
quantities = {
    "gain_adu_per_e": ("Gain (ADU/e-)", "Gain vs VR"),
    "noise_e": ("Noise (e-)", "Noise vs VR"),
    "nonlinearity_at_500_e": ("Nonlinearity at 500 e- (e-)", "Nonlinearity at 500 e- vs VR"),
}

# Read every table, keyed by the VR value parsed from its path.
records = []
for table in tables:
    vr = -int(re.search(r"/VR-(\d+)/", table).group(1))
    df = pd.read_csv(table).sort_values("ext")
    df["VR"] = vr
    records.append(df)

data = pd.concat(records, ignore_index=True).sort_values(["ext", "VR"])
extensions = sorted(data["ext"].unique())

print(data.sort_values(["VR", "ext"]).to_string(index=False))

output_dir = "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies"

# Pastel color per extension (more saturated so the lines aren't washed out).
pastel = ["#ff82be", "#96d4fd", "#ffc165", "#d29ce8", "#59e45e", "#76d7c4"]
ext_colors = {ext: pastel[i % len(pastel)] for i, ext in enumerate(extensions)}

for column, (ylabel, title) in quantities.items():
    plt.figure()
    for ext in extensions:
        sub = data[data["ext"] == ext]
        plt.plot(sub["VR"], sub[column], marker="o", color=ext_colors[ext], label=f"EXT{ext}")
    plt.gca().invert_xaxis()
    plt.xlabel("VR (V)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{column}_vs_VR.jpeg", dpi=350)

plt.show()

tables_VR_6 = ["/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-6/plots/avg_img_CV_250x3500x500_bin1x1_125_10_stitched_dbec0ad2/extension_summary.csv",
               "/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR_6_Diego_params/plots/cds_avg_img_L2_250x3500x500_1x1_L2_125_11_stitched_ab2690ae/extension_summary.csv"]

# Read both tables, keyed by VR and the delta(V) label parsed from the path.
delta_labels = {"VR-6": "Standard", "VR_6_Diego_params": "Increased"}
records = []
for table in tables_VR_6:
    vr = -int(re.search(r"/VR[-_](\d+)", table).group(1))
    folder = re.search(r"/VR_studies/(VR[-_][^/]+)/", table).group(1)
    df = pd.read_csv(table).sort_values("ext")
    df["VR"] = vr
    df["delta(V)"] = delta_labels[folder]
    records.append(df)

data_VR_6 = pd.concat(records, ignore_index=True)

print(data_VR_6.sort_values(["delta(V)", "ext", "VR"]).to_string(index=False))