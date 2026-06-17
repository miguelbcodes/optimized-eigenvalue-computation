import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

ALGO_STYLES = {
    "hessenberg_qr": {"color": "steelblue", "label": "Hessenberg + QR (general)"},
    "tridiagonal_qr": {"color": "darkorange", "label": "Tridiagonal + QR (symmetric)"},
}

def read_csv(filename):
    rows = []
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "algorithm": row["algorithm"],
                "N": int(row["N"]),
                "time_ms": float(row["time_ms"]),
                "memory_kb": int(row["memory_kb"]),
            })
    return rows

def aggregate(rows, algo):
    filtered = [r for r in rows if r["algorithm"] == algo]
    if not filtered:
        return None, None, None, None
    sizes = np.array([r["N"] for r in filtered])
    times = np.array([r["time_ms"] for r in filtered])
    mems = np.array([r["memory_kb"] for r in filtered])
    unique = np.unique(sizes)
    mean_t = np.array([times[sizes == n].mean() for n in unique])
    std_t = np.array([times[sizes == n].std() for n in unique])
    mean_m = np.array([mems[sizes == n].mean() for n in unique])
    std_m = np.array([mems[sizes == n].std() for n in unique])
    counts = np.array([np.sum(sizes == n) for n in unique])
    return unique, (mean_t, std_t), (mean_m, std_m), counts

def fit_power_law(x, y, min_x=0):
    mask = (x > min_x) & (y > 0)
    if mask.sum() < 2:
        return None, None
    coeffs = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return coeffs[0], np.exp(coeffs[1])

def plot_metric(algos_data, ylabel, title, fit_min_x, outfile):
    fig, ax = plt.subplots(figsize=(8, 5))

    for algo, (sizes, means, stds, counts) in algos_data.items():
        style = ALGO_STYLES.get(algo, {"color": "gray", "label": algo})
        has_multiple = np.any(counts > 1)

        if has_multiple:
            ax.errorbar(sizes, means, yerr=stds, fmt="o", color=style["color"],
                        capsize=4, markersize=6, zorder=5)
        else:
            ax.scatter(sizes, means, color=style["color"], s=60, zorder=5)

        exponent, coeff = fit_power_law(sizes, means, min_x=fit_min_x)
        if exponent is not None:
            x_fit = np.linspace(sizes.min(), sizes.max(), 200)
            y_fit = coeff * x_fit ** exponent
            ax.plot(x_fit, y_fit, "--", color=style["color"],
                    label=f"{style['label']} — O(N^{exponent:.2f})")

    ax.legend(fontsize=10)
    ax.set_xlabel("Matrix size (N)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")

def report_gains(rows, general="hessenberg_qr", optimized="tridiagonal_qr"):
    sizes_g, time_g, mem_g, _ = aggregate(rows, general)
    sizes_o, time_o, mem_o, _ = aggregate(rows, optimized)
    if sizes_g is None or sizes_o is None:
        return

    common = np.intersect1d(sizes_g, sizes_o)
    if common.size == 0:
        return

    print(f"\nPerformance gain ({optimized} vs {general}):")
    print(f"{'N':>6} {'time gain':>12} {'memory gain':>12}")
    for n in common:
        t_g = time_g[0][sizes_g == n][0]
        t_o = time_o[0][sizes_o == n][0]
        m_g = mem_g[0][sizes_g == n][0]
        m_o = mem_o[0][sizes_o == n][0]
        time_gain = (t_g - t_o) / t_g * 100 if t_g > 0 else float("nan")
        mem_gain = (m_g - m_o) / m_g * 100 if m_g > 0 else float("nan")
        print(f"{int(n):>6} {time_gain:>11.1f}% {mem_gain:>11.1f}%")

def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "benchmark_results.csv"
    rows = read_csv(filename)

    algorithms = sorted(set(r["algorithm"] for r in rows))

    time_data = {}
    mem_data = {}
    for algo in algorithms:
        sizes, time_stats, mem_stats, counts = aggregate(rows, algo)
        if sizes is not None:
            time_data[algo] = (sizes, time_stats[0], time_stats[1], counts)
            mem_data[algo] = (sizes, mem_stats[0], mem_stats[1], counts)

    plot_metric(time_data, "Time (ms)", "Scalability: N x Execution Time",
                0, "scalability_time.png")
    plot_metric(mem_data, "Peak memory (KB)", "Scalability: N x Peak Memory",
                50, "scalability_memory.png")

    for algo in algorithms:
        sizes, time_stats, mem_stats, counts = aggregate(rows, algo)
        if sizes is not None and np.any(counts > 1):
            print(f"\n[{algo}] Samples per N: {dict(zip(sizes.astype(int), counts.astype(int)))}")

    report_gains(rows)

if __name__ == "__main__":
    main()
