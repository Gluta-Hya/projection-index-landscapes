"""Gaya gambar bersama untuk naskah (palet tervalidasi, tinta teks, garis tipis)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BIRU, ORANYE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, REDUP = "#0b0b0b", "#52514e", "#898781"
GRID, SUMBU, PERMUKAAN = "#e1e0d9", "#c3c2b7", "#ffffff"
RAMP_BIRU = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP_BIRU = LinearSegmentedColormap.from_list("biru", RAMP_BIRU)
LEBAR = 6.3   # inci, lebar teks A4 dengan margin 2,5 cm

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7, "legend.frameon": False,
    "axes.edgecolor": SUMBU, "axes.linewidth": 0.8, "axes.labelcolor": TINTA2, "axes.titlecolor": TINTA,
    "xtick.color": TINTA2, "ytick.color": TINTA2, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "grid.linestyle": "-", "axes.axisbelow": True, "lines.linewidth": 1.4,
    "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round", "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    "axes.spines.top": False, "axes.spines.right": False, "mathtext.fontset": "dejavusans",
})

def label_panel(ax, huruf, x=-0.13, y=1.04):
    ax.text(x, y, huruf, transform=ax.transAxes, fontsize=9, fontweight="bold", color=TINTA, va="bottom")

from matplotlib.ticker import ScalarFormatter
class _Koma(ScalarFormatter):
    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")

def koma_semua(fig):
    """Desimal koma (gaya Indonesia) di semua sumbu linear 2D, termasuk colorbar."""
    for ax in fig.axes:
        if hasattr(ax, "zaxis"):
            continue
        for sumbu, skala in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if skala == "linear":
                sumbu.set_major_formatter(_Koma())
