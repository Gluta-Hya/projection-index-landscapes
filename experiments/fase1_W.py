"""
FASE 1: Uji ketahanan temuan W*.

Hipotesis dari run pertama:
  Fitur bimodalitas lanskap I(u) tidak bergantung pada m dan tau secara
  terpisah, tetapi pada panjang jendela  W = (m-1) * tau.

Yang diuji di sini:
  1. Banyak realisasi (kondisi awal berbeda)  -> pita ketidakpastian
  2. Dua panjang deret                         -> efek ukuran sampel
  3. m = 2..7                                  -> rentang dimensi lebih lebar
  4. Ringkasan KONTINU (rata-rata I3 atas arah) selain fraksi berambang 0.555
  5. Surrogate acak-fase per realisasi         -> null
  6. Skala waktu alami Lorenz (autokorelasi, mutual information, periode
     spektral) sebagai pembanding W*

Keluaran:
  fase1_hasil.csv     semua angka mentah
  fase1_ringkasan.csv W* per (m, panjang, realisasi)
  fase1_W.png         kurva vs W, dengan pita +-1 simpangan baku
  Ringkasan W* dicetak di terminal.

Jalankan:  python fase1_W.py
Kebutuhan: numpy scipy matplotlib  (sudah terpasang dari run sebelumnya)
Waktu: sekitar 20 sampai 40 menit di laptop biasa (progres dicetak per realisasi).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import csv
import time

DT = 0.01
RNG = np.random.default_rng(2026)

# ---------------- pengaturan ----------------
N_REALISASI = 10
PANJANG_T = (75.0, 150.0)          # detik data setelah transien dibuang
M_LIST = (2, 3, 4, 5, 6, 7)
W_GRID = np.round(np.arange(0.08, 1.29, 0.08), 2)   # detik
K_ARAH = 600                         # banyak arah acak per m
SUB = 3                              # ambil tiap SUB sampel saat menghitung indeks
AMBANG_BIMODAL = 0.555


# ---------------- data ----------------
def lorenz_x(T, buang=20.0, ic=None):
    f = lambda t, s: [10 * (s[1] - s[0]), s[0] * (28 - s[2]) - s[1], s[0] * s[1] - 8 / 3 * s[2]]
    tt = np.arange(0, T + buang, DT)
    ic = ic if ic is not None else RNG.normal([1, 1, 20], 5)
    sol = solve_ivp(f, (0, T + buang), ic, t_eval=tt, rtol=1e-8, atol=1e-8)
    return sol.y[0, int(buang / DT):]


def surrogate_acak_fase(x):
    X = np.fft.rfft(x - x.mean())
    fase = RNG.uniform(0, 2 * np.pi, len(X))
    fase[0] = 0
    if len(x) % 2 == 0:
        fase[-1] = 0
    return np.fft.irfft(np.abs(X) * np.exp(1j * fase), n=len(x)) + x.mean()


def embed(x, m, tau):
    N = len(x) - (m - 1) * tau
    V = np.stack([x[(m - 1 - k) * tau:(m - 1 - k) * tau + N] for k in range(m)], axis=1)
    return (V - V.mean(0)) / V.std(0)


def arah(m, K):
    U = RNG.standard_normal((K, m))
    return U / np.linalg.norm(U, axis=1, keepdims=True)


# ---------------- indeks (tervektorisasi atas semua arah) ----------------
G = lambda y: np.log(np.cosh(y))
EGN = G(np.random.default_rng(0).standard_normal(400_000)).mean()


def indeks_semua_arah(V, U):
    Y = V[::SUB] @ U.T                       # (n, K)
    Y = (Y - Y.mean(0)) / Y.std(0)
    n = Y.shape[0]
    I1 = (G(Y).mean(0) - EGN) ** 2           # negentropi log cosh
    g = (Y ** 3).mean(0)                     # skewness
    k = (Y ** 4).mean(0) - 3                 # excess kurtosis
    I3 = (g ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    return I1, I3


def ringkas(V, U):
    I1, I3 = indeks_semua_arah(V, U)
    return dict(I3_rata=I3.mean(), I3_fraksi=(I3 > AMBANG_BIMODAL).mean(),
                I1_rata=I1.mean(), I1_maks=I1.max())


# ---------------- skala waktu alami ----------------
def skala_waktu(x):
    xc = x - x.mean()
    n = len(xc)
    F = np.fft.rfft(xc, 2 * n)
    ac = np.fft.irfft(np.abs(F) ** 2)[:n]
    ac /= ac[0]
    nol = np.argmax(ac < 0) * DT

    def mi(lag, bins=32):
        h, _, _ = np.histogram2d(x[:-lag], x[lag:], bins)
        p = h / h.sum()
        px, py = p.sum(1), p.sum(0)
        nz = p > 0
        return (p[nz] * np.log(p[nz] / np.outer(px, py)[nz])).sum()
    lags = np.arange(1, 150)
    mis = np.array([mi(l) for l in lags])
    i = np.argmax((mis[1:-1] < mis[:-2]) & (mis[1:-1] < mis[2:])) + 1
    mi_min = lags[i] * DT

    # Periode putaran di dalam satu sayap: median selang antar puncak lokal |x|.
    # (Puncak spektrum tidak dipakai: ia didominasi perpindahan sayap yang lambat.)
    pk, _ = find_peaks(np.abs(x))
    periode = np.median(np.diff(pk)) * DT
    return dict(autokorelasi_nol=nol, mi_minimum=mi_min, periode_sayap=periode)


# ---------------- eksperimen ----------------
def jalankan():
    U_per_m = {m: arah(m, K_ARAH) for m in M_LIST}
    baris, skala = [], []
    t0 = time.time()
    for r in range(N_REALISASI):
        x_penuh = lorenz_x(max(PANJANG_T))
        skala.append(skala_waktu(x_penuh))
        for T in PANJANG_T:
            x = x_penuh[:int(T / DT)]
            s = surrogate_acak_fase(x)
            for m in M_LIST:
                for W in W_GRID:
                    tau = int(round(W / ((m - 1) * DT)))
                    if tau < 1:
                        continue
                    W_nyata = (m - 1) * tau * DT
                    for nama, sig in (("lorenz", x), ("surrogate", s)):
                        d = ringkas(embed(sig, m, tau), U_per_m[m])
                        d.update(realisasi=r, T=T, m=m, tau=tau, W=round(W_nyata, 3), sistem=nama)
                        baris.append(d)
        print(f"realisasi {r + 1}/{N_REALISASI} selesai ({time.time() - t0:.0f} detik)")
    return baris, skala


def W_bintang(baris, kunci="I3_rata"):
    """W* = W yang memaksimalkan kurva (dihaluskan ringan) per (realisasi, T, m)."""
    out = []
    for r in sorted({b["realisasi"] for b in baris}):
        for T in PANJANG_T:
            for m in M_LIST:
                d = sorted([b for b in baris if b["realisasi"] == r and b["T"] == T
                            and b["m"] == m and b["sistem"] == "lorenz"], key=lambda b: b["W"])
                if len(d) < 3:
                    continue
                Ws = np.array([b["W"] for b in d])
                y = uniform_filter1d(np.array([b[kunci] for b in d]), 3)
                out.append(dict(realisasi=r, T=T, m=m, W_bintang=Ws[np.argmax(y)], puncak=y.max()))
    return out


def gambar(baris, skala):
    kolom = [("I3_rata", "rata-rata $I_3$ atas arah"),
             ("I3_fraksi", "fraksi arah bimodal ($I_3>0.555$)"),
             ("I1_rata", "rata-rata $I_1$ (negentropi)")]
    warna = plt.cm.viridis(np.linspace(0, 0.9, len(M_LIST)))
    fig, axs = plt.subplots(len(PANJANG_T), len(kolom), figsize=(17, 5 * len(PANJANG_T)), squeeze=False)
    sk = {k: np.mean([s[k] for s in skala]) for k in skala[0]}
    for i, T in enumerate(PANJANG_T):
        for j, (k, judul) in enumerate(kolom):
            ax = axs[i, j]
            for c, m in zip(warna, M_LIST):
                for sistem, gaya, alpha in (("lorenz", "-", 0.25), ("surrogate", ":", 0.0)):
                    d = [b for b in baris if b["T"] == T and b["m"] == m and b["sistem"] == sistem]
                    Ws = sorted({b["W"] for b in d})
                    mu = np.array([np.mean([b[k] for b in d if b["W"] == w]) for w in Ws])
                    sd = np.array([np.std([b[k] for b in d if b["W"] == w]) for w in Ws])
                    ax.plot(Ws, mu, gaya, color=c, lw=1.6 if sistem == "lorenz" else 1,
                            label=f"m={m}" if sistem == "lorenz" else None)
                    if alpha:
                        ax.fill_between(Ws, mu - sd, mu + sd, color=c, alpha=alpha)
            for (nama, nilai), ls in zip(sk.items(), ("--", "-.", (0, (1, 1)))):
                ax.axvline(nilai, color="k", ls=ls, lw=0.9, label=f"{nama} = {nilai:.2f}s" if j == 0 else None)
            ax.set_xlim(0, W_GRID.max() + 0.05)
            ax.set_title(f"{judul}, T = {T:.0f} s")
            ax.set_xlabel("W = (m-1)τ  (detik)")
            if j == 0:
                ax.legend(fontsize=8, ncol=2)
    fig.suptitle("Garis tebal: Lorenz (pita = ±1 sd antar realisasi).  Garis titik-titik: surrogate.", y=1.0)
    plt.tight_layout()
    plt.savefig("fase1_W.png", dpi=120, bbox_inches="tight")


def simpan(nama, baris):
    with open(nama, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(baris[0].keys()))
        w.writeheader()
        w.writerows(baris)


if __name__ == "__main__":
    baris, skala = jalankan()
    simpan("fase1_hasil.csv", baris)
    wb = W_bintang(baris)
    simpan("fase1_ringkasan.csv", wb)
    gambar(baris, skala)

    print("\nSkala waktu alami (rata-rata antar realisasi):")
    for k in skala[0]:
        v = [s[k] for s in skala]
        print(f"  {k:18s} {np.mean(v):.3f} ± {np.std(v):.3f} s")
    print("\nW* (puncak rata-rata I3), rata-rata ± sd antar realisasi:")
    for T in PANJANG_T:
        print(f"  T = {T:.0f} s")
        for m in M_LIST:
            v = [d["W_bintang"] for d in wb if d["T"] == T and d["m"] == m]
            p = [d["puncak"] for d in wb if d["T"] == T and d["m"] == m]
            if v:
                print(f"    m={m}:  W* = {np.mean(v):.3f} ± {np.std(v):.3f} s   puncak = {np.mean(p):.3f}")
    print("\nTersimpan: fase1_hasil.csv, fase1_ringkasan.csv, fase1_W.png")
