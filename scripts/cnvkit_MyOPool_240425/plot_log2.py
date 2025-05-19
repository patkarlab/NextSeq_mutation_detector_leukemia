import pandas as pd
import matplotlib.pyplot as plt
import statistics

df = pd.read_csv("Reference_combpanel.cnn", sep="\t")

mean_log2 = df["log2"].mean()
median_log2 = df["log2"].median()
std_log2 = df["log2"].std()
plus_2sd = mean_log2 + 2 * std_log2
minus_2sd = mean_log2 - 2 * std_log2

plt.figure(figsize=(12, 6))
plt.scatter(df.index, df["log2"], s=10, color='teal')

for i, (idx, row) in enumerate(df.iterrows()):
    if row["log2"] > plus_2sd or row["log2"] < minus_2sd:
        plt.annotate(
            row["gene"],
            (idx, row["log2"]),
            textcoords="offset points",
            xytext=(0, 5),
            ha='center',
            fontsize=8,
            rotation=30
        )

plt.axhline(mean_log2, color='black', linestyle='-', label=f'Mean: {mean_log2:.2f}')
plt.axhline(median_log2, color='grey', linestyle='--', label=f'Median: {median_log2:.2f}')
plt.axhline(plus_2sd, color='red', linestyle='-', label=f'+2 SD: {plus_2sd:.2f}')
plt.axhline(minus_2sd, color='red', linestyle='-', label=f'-2 SD: {minus_2sd:.2f}')

plt.title("Log2 Plot with Statistics")
plt.xlabel("Index")
plt.ylabel("Log2")
plt.legend(loc="upper right", bbox_to_anchor=(1.15, 1))
plt.tight_layout()

plt.savefig("log2_plot_with_stats.pdf", format="pdf")

