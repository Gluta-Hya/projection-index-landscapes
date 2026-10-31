"""
FASE 3c: Apakah bentuk gelombang nyata Lorenz menjelaskan kenaikan c_m = W*/P?

Hipotesis (Larry): gelombang Lorenz di dalam satu sayap tidak sinusoidal
(memiliki harmonik), dan harmonik itu menaikkan c_m, terutama untuk m besar.

Tangga model. Setiap anak tangga menambah SATU ciri Lorenz yang diukur langsung:
  A  cos           : x = s (1 + a cos),         a/c = median Lorenz, konstan
  B  harmonik 1-2  : bentuk Lorenz dipotong ke harmonik 1 dan 2 saja
  C  bentuk rata2  : bentuk rata-rata Lorenz (semua harmonik), a/c konstan
  D  + a/c acak    : bentuk rata-rata, a/c per putaran diundi dari distribusi Lorenz
                     (tanpa ekor atas yang membuat x memotong nol; lihat peringatan di bangun())
  E  putaran nyata : setiap putaran memakai bentuk DAN a/c dari satu putaran Lorenz
                     nyata (diundi), jadi variasi bentuk antar putaran ikut masuk
Semua anak tangga memakai perpindahan sayap yang sama (lama tinggal geometrik,
rata-rata D = 2.2 putaran), dan perpindahan/pergantian putaran terjadi di titik
|x| minimum, seperti model Fase 3.

Yang TIDAK ada di model mana pun: lintasan perpindahan sayap itu sendiri
(lompatan seketika di model, ekskursi panjang di Lorenz) dan variasi periode.
Selisih antara E dan Lorenz asli adalah batas atas pengaruh kedua hal itu.

Bentuk gelombang diukur dari Lorenz rho=28: segmen antara dua minimum |x(t)|
berturut-turut yang x-nya tidak berganti tanda (satu putaran di dalam sayap),
lalu |x| dinormalisasi:  g = (|x| - c) / a,  c = rata-rata, a = setengah jangkauan.

Keluaran: fase3c_hasil.csv, fase3c_ringkasan.csv, fase3c_cm.png, fase3c_bentuk.png
Jalankan:  python fase3c_gelombang.py   (perlu fase2_sistem.py; fase2_ringkasan.csv opsional)
Waktu: sekitar 10 menit.
"""
import csv
import os
import time
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import fase2_sistem as F   # embed, arah, ringkas, daftar_tau, puncak_parabola yang SAMA dengan Fase 2

RNG = np.random.default_rng(1618)
NPH = 200                       # titik fase per putaran untuk bentuk gelombang
DT = 0.01                       # 100 sampel per putaran pada model (P = 1)
N_SAMPEL = 60_000
N_PUTARAN = int(N_SAMPEL * DT) + 2
D_SAYAP = 2.2
M_LIST = (2, 3, 4, 5, 6)
N_REALISASI = 3
W_MAKS, W_LANGKAH = 1.6, 0.03
K_ARAH = 500


# ---------------- 1. ukur bentuk gelombang Lorenz ----------------
def ukur_lorenz():
    dt = 0.005
    f = lambda t, s: [10 * (s[1] - s[0]), s[0] * (28 - s[2]) - s[1], s[0] * s[1] - 8 / 3 * s[2]]
    simpan_n = F.N_SAMPEL
    F.N_SAMPEL = 200_000                                  # sekitar 1300 putaran
    S = F._flow(f, [1.0, 1.0, 20.0], dt, 30.0)
    F.N_SAMPEL = simpan_n
    x = S[0]
    # Satu putaran = dari satu minimum |x| ke minimum |x| berikutnya, sehingga sambungan
    # antar putaran jatuh di minimum |x|, sama seperti titik pergantian putaran di model.
    # (Memotong di maksimum z menimbulkan patahan buatan di tengah putaran.)
    pm, _ = find_peaks(-np.abs(x))
    bentuk, rasio = [], []
    fase = np.linspace(0, 1, NPH, endpoint=False)
    for a, b in zip(pm[:-1], pm[1:]):
        seg = x[a:b + 1]
        if not (np.all(seg > 0) or np.all(seg < 0)):
            continue                                      # putaran yang berpindah sayap dibuang
        g = np.interp(fase, np.linspace(0, 1, len(seg)), np.abs(seg))
        c, amp = g.mean(), (g.max() - g.min()) / 2
        g = (g - c) / amp
        bentuk.append(np.roll(g, NPH // 2))               # minimum (ujung segmen) di theta = pi
        rasio.append(amp / c)
    return np.array(bentuk), np.array(rasio)


def normalisasi(g):
    g = g - g.mean()
    g = g / ((g.max() - g.min()) / 2)
    return np.roll(g, NPH // 2 - int(np.argmin(g)))


def potong_harmonik(g, n):
    G = np.fft.rfft(g)
    G[n + 1:] = 0
    return normalisasi(np.fft.irfft(G, n=NPH))


# ---------------- 2. model ----------------
def tanda_sayap():
    s, tanda = [], 1.0
    while len(s) < N_PUTARAN:
        s += [tanda] * RNG.geometric(1.0 / D_SAYAP)
        tanda = -tanda
    return np.array(s[:N_PUTARAN])


def deret(s, bentuk_per_putaran, rasio_per_putaran):
    """x = s_r (1 + a_r g_r(theta)); putaran r berganti di theta = pi (minimum |x|)."""
    theta = 2 * np.pi * np.arange(N_SAMPEL) * DT
    r = np.floor((theta + np.pi) / (2 * np.pi)).astype(int)
    k = np.round((theta % (2 * np.pi)) / (2 * np.pi) * NPH).astype(int) % NPH
    return s[r] * (1.0 + rasio_per_putaran[r] * bentuk_per_putaran[r, k])


def bangun(jenis, s, G_rata, G_12, bentuk_lor, rasio_lor, ac_med):
    n = N_PUTARAN
    if jenis == "A cos":
        g = normalisasi(np.cos(2 * np.pi * np.arange(NPH) / NPH))
        return deret(s, np.tile(g, (n, 1)), np.full(n, ac_med))
    if jenis == "B harmonik 1-2":
        return deret(s, np.tile(G_12, (n, 1)), np.full(n, ac_med))
    if jenis == "C bentuk rata2":
        return deret(s, np.tile(G_rata, (n, 1)), np.full(n, ac_med))
    if jenis == "D + a/c acak":
        # Hanya a/c yang tidak membuat x memotong nol di dalam sayap (1 + a min g > 0).
        # PERINGATAN: sekitar 40% nilai terbuang, karena di Lorenz putaran ber-a/c besar
        # punya lembah yang lebih dangkal (bentuk dan amplitudo berkorelasi). Model D
        # mengabaikan korelasi itu, jadi hasilnya bias ke bawah; model E yang benar.
        sah = rasio_lor[1 + rasio_lor * G_rata.min() > 0.02]
        return deret(s, np.tile(G_rata, (n, 1)), RNG.choice(sah, n))
    if jenis == "E putaran nyata":
        idx = RNG.integers(0, len(bentuk_lor), n)
        return deret(s, bentuk_lor[idx], rasio_lor[idx])
    raise ValueError(jenis)


JENIS = ("A cos", "B harmonik 1-2", "C bentuk rata2", "D + a/c acak", "E putaran nyata")


# ---------------- 3. eksperimen ----------------
def kurva_dan_puncak(x, m, U):
    W_target = np.arange(W_LANGKAH, W_MAKS + 1e-9, W_LANGKAH)
    xs, ys = [], []
    for tau in F.daftar_tau(m, DT, W_target):
        xs.append((m - 1) * tau * DT)
        ys.append(F.ringkas(F.embed(x, m, tau), U)["I3_rata"])
    wb, pk, tepi = F.puncak_parabola(xs, ys)
    return xs, ys, wb, pk, tepi


def lorenz_acuan(path="fase2_ringkasan.csv", sistem="lorenz28"):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh) if r["sistem"] == sistem]
    return [np.mean([float(r["Wbintang_per_P"]) for r in rows if int(r["m"]) == m]) for m in M_LIST]


def simpan(nama, rows):
    with open(nama, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    t0 = time.time()
    bentuk_lor, rasio_lor = ukur_lorenz()
    G_rata = normalisasi(bentuk_lor.mean(0))
    G_12 = potong_harmonik(G_rata, 2)
    ac_med = float(np.median(rasio_lor))
    H = np.abs(np.fft.rfft(G_rata))[1:7]
    print(f"Lorenz: {len(bentuk_lor)} putaran di dalam sayap, a/c median = {ac_med:.3f}")
    print(f"Model D membuang {np.mean(1 + rasio_lor * G_rata.min() <= 0.02):.1%} nilai a/c (memotong nol)")
    print("Amplitudo harmonik 1..6 relatif ke harmonik 1: " + "  ".join(f"{h:.3f}" for h in H / H[0]))
    print(f"({time.time() - t0:.0f} detik)\n")

    U = {m: F.arah(m, K_ARAH) for m in M_LIST}
    baris, ringkas = [], []
    for r in range(N_REALISASI):
        s = tanda_sayap()                                  # semua anak tangga berbagi barisan sayap
        for jenis in JENIS:
            x = bangun(jenis, s, G_rata, G_12, bentuk_lor, rasio_lor, ac_med)
            for m in M_LIST:
                xs, ys, wb, pk, tepi = kurva_dan_puncak(x, m, U[m])
                baris += [dict(jenis=jenis, realisasi=r, m=m, W_per_P=round(w, 6), I3_rata=v)
                          for w, v in zip(xs, ys)]
                ringkas.append(dict(jenis=jenis, realisasi=r, m=m, c_m=wb, puncak=pk, di_tepi_grid=tepi))
            print(f"realisasi {r + 1}/{N_REALISASI}  {jenis:16s} ({time.time() - t0:.0f} detik)")
    simpan("fase3c_hasil.csv", baris)
    simpan("fase3c_ringkasan.csv", ringkas)
    lor = lorenz_acuan()

    def stat(jenis):
        return [(np.mean(v), np.std(v)) for v in
                ([d["c_m"] for d in ringkas if d["jenis"] == jenis and d["m"] == m] for m in M_LIST)]

    # ---------- gambar: bentuk gelombang ----------
    fig, ax = plt.subplots(figsize=(8, 4))
    fase = np.arange(NPH) / NPH
    for g in bentuk_lor[RNG.choice(len(bentuk_lor), 40, replace=False)]:
        ax.plot(fase, g, color="0.85", lw=0.6)
    ax.plot(fase, normalisasi(np.cos(2 * np.pi * fase)), "k--", label="cos (model A)")
    ax.plot(fase, G_12, color="tab:orange", label="harmonik 1-2 (model B)")
    ax.plot(fase, G_rata, color="tab:red", lw=2, label="bentuk rata-rata Lorenz (model C)")
    ax.set_xlabel("fase putaran (minimum |x| di 0.5)")
    ax.set_ylabel("g = (|x| - c) / a")
    ax.set_title("Bentuk gelombang di dalam sayap (abu-abu: 40 putaran nyata)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("fase3c_bentuk.png", dpi=120)

    # ---------- gambar: c_m ----------
    fig, ax = plt.subplots(figsize=(8, 5))
    warna = plt.cm.viridis(np.linspace(0, 0.9, len(JENIS)))
    for c, jenis in zip(warna, JENIS):
        mu, sd = zip(*stat(jenis))
        ax.errorbar(M_LIST, mu, yerr=sd, fmt="o-", color=c, capsize=3, label=jenis)
    if lor:
        ax.plot(M_LIST, lor, "ks--", lw=2, label="Lorenz ρ=28 asli")
    ax.set_xlabel("m")
    ax.set_ylabel("c_m = W*/P")
    ax.set_xticks(M_LIST)
    ax.set_title("Tangga model: seberapa jauh setiap ciri Lorenz menutup selisih?")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("fase3c_cm.png", dpi=120)

    print("\nc_m = W*/P (rata-rata antar realisasi);  kemiringan = c_6 - c_3")
    if lor:
        print(f"  {'Lorenz asli':16s} " + "  ".join(f"m{m}={v:.3f}" for m, v in zip(M_LIST, lor))
              + f"   kemiringan={lor[4] - lor[1]:+.3f}")
    for jenis in JENIS:
        st = stat(jenis)
        tepi = any(d["di_tepi_grid"] for d in ringkas if d["jenis"] == jenis)
        print(f"  {jenis:16s} " + "  ".join(f"m{m}={mu:.3f}" for m, (mu, _) in zip(M_LIST, st))
              + f"   kemiringan={st[4][0] - st[1][0]:+.3f}" + ("   (ada puncak di tepi grid!)" if tepi else ""))
    print("\nTersimpan: fase3c_hasil.csv, fase3c_ringkasan.csv, fase3c_cm.png, fase3c_bentuk.png")
