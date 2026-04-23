import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).parent
BED = ROOT / "alignments" / "windows_1kb.bed"
BAMS = [ROOT / "alignments" / f"assembly{i}.bam" for i in range(1, 14)]
OUT = ROOT / "alignments" / "coverage_heatmap.png"

# samtools bedcov returns summed base depth per window per BAM
cmd = ["samtools", "bedcov", str(BED), *map(str, BAMS)]
res = subprocess.run(cmd, capture_output=True, text=True, check=True)

rows = [line.split("\t") for line in res.stdout.strip().split("\n")]
starts = np.array([int(r[1]) for r in rows])
# columns 3..15 are per-BAM summed base depths; divide by window size for mean depth
window_sizes = np.array([int(r[2]) - int(r[1]) for r in rows])
depths = np.array([[int(v) for v in r[3:]] for r in rows], dtype=float)
mean_depth = depths / window_sizes[:, None]

# matrix: samples (rows) x windows (cols)
mat = mean_depth.T

fig, ax = plt.subplots(figsize=(16, 5))
im = ax.imshow(
    mat,
    aspect="auto",
    cmap="viridis",
    norm=LogNorm(vmin=0.1, vmax=max(mat.max(), 1)),
    interpolation="nearest",
    extent=[starts[0] / 1e6, (starts[-1] + 1000) / 1e6, 13.5, 0.5],
)
ax.set_yticks(range(1, 14))
ax.set_yticklabels([f"assembly{i}" for i in range(1, 14)])
ax.set_xlabel("NC_007795.1 position (Mb)")
ax.set_ylabel("MRSA isolate")
ax.set_title("Per-1kb mean coverage vs. NCTC 8325 (log scale)")
cbar = fig.colorbar(im, ax=ax, pad=0.01)
cbar.set_label("Mean depth")
fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"wrote {OUT}  shape={mat.shape}")
