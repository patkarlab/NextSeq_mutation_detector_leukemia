import pandas as pd
import matplotlib.pyplot as plt
import statistics

df = pd.read_csv("Reference_combpanel.cnn", sep="\t")

mean_depth = df["depth"].mean()
median_depth = df["depth"].median()
std_depth = df["depth"].std()
plus_2sd = mean_depth + 2 * std_depth
minus_2sd = mean_depth - 2 * std_depth

plt.figure(figsize=(12, 6))
plt.scatter(df.index, df["depth"], s=10, color='teal')

for i, (idx, row) in enumerate(df.iterrows()):
    if row["depth"] > plus_2sd or row["depth"] < minus_2sd:
        plt.annotate(
            row["gene"],
            (idx, row["depth"]),
            textcoords="offset points",
            xytext=(0, 5),
            ha='center',
            fontsize=8,
            rotation=30
        )

plt.axhline(mean_depth, color='black', linestyle='-', label=f'Mean: {mean_depth:.2f}')
plt.axhline(median_depth, color='grey', linestyle='--', label=f'Median: {median_depth:.2f}')
plt.axhline(plus_2sd, color='red', linestyle='-', label=f'+2 SD: {plus_2sd:.2f}')
plt.axhline(minus_2sd, color='red', linestyle='-', label=f'-2 SD: {minus_2sd:.2f}')

plt.title("Depth Plot with Statistics")
plt.xlabel("Index")
plt.ylabel("Depth")
plt.legend(loc="upper right", bbox_to_anchor=(1.15, 1))
plt.tight_layout()

plt.savefig("depth_plot_with_stats.pdf", format="pdf")


outliers = df[(df["depth"] > plus_2sd) | (df["depth"] < minus_2sd)]
print(outliers)

outliers.to_csv("depth_outliers_above_below_2SD.tsv", sep="\t", index=False)

