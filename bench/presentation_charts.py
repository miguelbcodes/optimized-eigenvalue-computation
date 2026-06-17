"""Generate didactic charts for the presentation (apresentacao.html).

Produces FOUR PNGs in a consistent, Portuguese-labeled style:

  Phase 2 (general algorithm only — does NOT spoil the optimization):
    - scalability_time_general.png
    - scalability_memory_general.png

  Phase 3 (general vs. optimized comparison):
    - scalability_time_comparison.png
    - scalability_memory_comparison.png

Reuses read_csv/aggregate/fit_power_law from scalability_analysis.py so the
methodology (mean +/- std, power-law fit) matches the report's charts.

Usage:  python3 bench/presentation_charts.py [benchmark_results.csv]
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scalability_analysis import read_csv, aggregate, fit_power_law

# Palette tied to the deck: slate = general, green = optimized ("the better one").
STYLES = {
    "hessenberg_qr": {"color": "#4a5158", "label": "Hessenberg + QR (geral)"},
    "tridiagonal_qr": {"color": "#2f7d4f", "label": "Tridiagonal + QL (simétrico)"},
}


def plot(series, ylabel, title, fit_min_x, outfile):
    """series: list of (algo, sizes, means, stds, counts)."""
    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    for algo, sizes, means, stds, counts in series:
        style = STYLES.get(algo, {"color": "gray", "label": algo})
        if np.any(counts > 1):
            ax.errorbar(sizes, means, yerr=stds, fmt="o", color=style["color"],
                        capsize=4, markersize=6, zorder=5)
        else:
            ax.scatter(sizes, means, color=style["color"], s=55, zorder=5)

        exponent, coeff = fit_power_law(sizes, means, min_x=fit_min_x)
        if exponent is not None:
            x_fit = np.linspace(sizes.min(), sizes.max(), 200)
            ax.plot(x_fit, coeff * x_fit ** exponent, "--", color=style["color"],
                    linewidth=1.8, label=f"{style['label']} — O(N^{exponent:.2f})")

    ax.legend(fontsize=10, frameon=True, framealpha=0.95)
    ax.set_xlabel("Tamanho da matriz (N)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold", color="#23272b")
    ax.grid(True, alpha=0.3)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def collect(rows, algo):
    sizes, time_stats, mem_stats, counts = aggregate(rows, algo)
    if sizes is None:
        return None, None
    time_series = (algo, sizes, time_stats[0], time_stats[1], counts)
    mem_series = (algo, sizes, mem_stats[0], mem_stats[1], counts)
    return time_series, mem_series


def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "benchmark_results.csv"
    rows = read_csv(filename)

    gen_t, gen_m = collect(rows, "hessenberg_qr")
    opt_t, opt_m = collect(rows, "tridiagonal_qr")
    if gen_t is None:
        sys.exit("No 'hessenberg_qr' rows found in CSV.")

    # --- Phase 2: general algorithm only ---
    plot([gen_t], "Tempo (ms)",
         "Escalabilidade — Algoritmo Geral (Hessenberg + QR)",
         0, "scalability_time_general.png")
    plot([gen_m], "Pico de memória (KB)",
         "Escalabilidade — Algoritmo Geral (Hessenberg + QR)",
         50, "scalability_memory_general.png")

    # --- Phase 3: comparison (only if the optimized series is present) ---
    if opt_t is not None:
        plot([gen_t, opt_t], "Tempo (ms)",
             "Comparação — Geral vs. Otimizado (simétrico)",
             0, "scalability_time_comparison.png")
        plot([gen_m, opt_m], "Pico de memória (KB)",
             "Comparação — Geral vs. Otimizado (simétrico)",
             50, "scalability_memory_comparison.png")
    else:
        print("Optimized series absent — skipped comparison charts.")


if __name__ == "__main__":
    main()
