"""
FASE 3b: Apakah amplitudo spiral menjelaskan kenaikan c_m = W*/P terhadap m?

Hipotesis (Larry):
  Di dalam satu sayap, amplitudo osilasi tumbuh A_{j+1} = lambda * A_j,
  dan lintasan berpindah sayap saat amplitudo mencapai A_crit.

Catatan rancangan: memaksa lambda dan lama-tinggal independen membuat amplitudo
awal jatuh mendekati nol (lambda = 2, D = 8 -> separuh putaran beramplitudo < 0.05),
jadi model menjadi tidak fisis. Dalam model Larry, lambda, rentang amplitudo, dan
lama tinggal memang saling terkait. Karena itu yang diuji adalah PERBANDINGAN
DENGAN KONTROL ACAK, yang menjaga segalanya sama kecuali urutan pertumbuhan:

  Eksperimen 1 (amplitudo konstan): a/c = 0.30 ... 0.90, perpindahan ala Lorenz
      -> berapa besar pengaruh TINGKAT amplitudo saja terhadap c_m?
  Eksperimen 2 (spiral vs acak), untuk lambda = 1.25, 1.6, 2.0:
      spiral : amplitudo awal ~ U(0.3, 0.5), A_j = min(A_0 lambda^j, A_crit),
               pindah sayap setelah putaran yang mencapai A_crit = 0.9
      acak   : barisan sayap s(t) IDENTIK, multiset amplitudo tiap sayap IDENTIK,
               tetapi urutan amplitudo di dalam tiap sayap diacak
      -> kalau spiral dan acak memberi c_m yang sama, maka yang penting adalah
         keragaman amplitudo, bukan pertumbuhan teraturnya.

Model: x(t) = s(t) * ( 1 + A(t) cos(2 pi t / P) ), P = 1, pindah sayap di fase pi.
Keluaran:
  fase3b_hasil.csv, fase3b_ringkasan.csv, fase3b_cm.png
Jalankan:  python fase3b_spiral.py   (perlu fase2_sistem.py; fase2_ringkasan.csv opsional)
Waktu: sekitar 10 sampai 20 menit.
"""
import csv
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import fase2_sistem as F   # embed, arah, ringkas, daftar_tau, puncak_parabola yang SAMA dengan Fase 2

RNG = np.random.default_rng(2718)
DT = 0.01                       # 100 sampel per putaran
N_SAMPEL = 60_000               # 600 putaran
N_PUTARAN = int(N_SAMPEL * DT) + 2
M_LIST = (2, 3, 4, 5, 6)
N_REALISASI = 3
W_MAKS, W_LANGKAH = 1.6, 0.03
K_ARAH = 500

D_KONSTAN = 2.2                 # rata-rata lama di sayap untuk eksperimen 1
AC_LIST = (0.30, 0.45, 0.60, 0.75, 0.90)
LAMBDA_LIST = (1.25, 1.6, 2.0)
A0_RENTANG = (0.3, 0.5)
A_CRIT = 0.9


# ---------------- pembangkit barisan per-putaran ----------------
def putaran_konstan(a, D):
    """Kembalikan (s, A) per putaran: amplitudo tetap a, lama tinggal geometrik rata-rata D."""
    s, A, tanda = [], [], 1.0
    while len(s) < N_PUTARAN:
        K = RNG.geometric(1.0 / D)
        s += [tanda] * K
        A += [a] * K
        tanda = -tanda
    return np.array(s[:N_PUTARAN]), np.array(A[:N_PUTARAN])


def putaran_spiral(lam):
    """Model Larry. Kembalikan (s, A, daftar_sayap) agar kontrol acak bisa memakai sayap yang sama."""
    sayap, n, tanda = [], 0, 1.0
    while n < N_PUTARAN:
        amp = [RNG.uniform(*A0_RENTANG)]
        while amp[-1] < A_CRIT:
            amp.append(min(amp[-1] * lam, A_CRIT))
        sayap.append((tanda, np.array(amp)))
        n += len(amp)
        tanda = -tanda
    return sayap


def dari_sayap(sayap, acak=False):
    # Catatan: deret dipotong di N_PUTARAN, jadi sayap TERAKHIR bisa terpotong berbeda
    # antara spiral dan acak. Semua sayap lain identik isinya (hanya urutan yang beda).
    s, A = [], []
    for tanda, amp in sayap:
        amp = RNG.permutation(amp) if acak else amp
        s += [tanda] * len(amp)
        A += list(amp)
    return np.array(s[:N_PUTARAN]), np.array(A[:N_PUTARAN])


def deret(s, A):
    theta = 2 * np.pi * np.arange(N_SAMPEL) * DT
    r = np.floor((theta + np.pi) / (2 * np.pi)).astype(int)   # putaran berganti di theta = pi
    return s[r] * (1.0 + A[r] * np.cos(theta))


# ---------------- kurva dan puncak ----------------
def kurva_dan_puncak(x, m, U):
    W_target = np.arange(W_LANGKAH, W_MAKS + 1e-9, W_LANGKAH)
    xs, ys = [], []
    for tau in F.daftar_tau(m, DT, W_target):                    # tiap tau unik sekali
        xs.append((m - 1) * tau * DT)
        ys.append(F.ringkas(F.embed(x, m, tau), U)["I3_rata"])
    wb, pk, tepi = F.puncak_parabola(xs, ys)
    return xs, ys, wb, pk, tepi


def lama_tinggal_rata(s):
    sw = np.flatnonzero(np.diff(s) != 0)
    return np.diff(sw).mean() if len(sw) > 2 else np.nan


def lorenz_acuan(path="fase2_ringkasan.csv", sistem="lorenz28"):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["sistem"] == sistem]
    return [np.mean([float(r["Wbintang_per_P"]) for r in rows if int(r["m"]) == m]) for m in M_LIST]


def jalankan():
    U = {m: F.arah(m, K_ARAH) for m in M_LIST}
    konfig = [("konstan", f"a/c={a:.2f}", a) for a in AC_LIST]
    konfig += [(jenis, f"λ={lam}", lam) for lam in LAMBDA_LIST for jenis in ("spiral", "acak")]
    baris, ringkas = [], []
    t0 = time.time()
    for r in range(N_REALISASI):
        sayap_per_lam = {lam: putaran_spiral(lam) for lam in LAMBDA_LIST}   # spiral dan acak berbagi sayap
        for jenis, label, par in konfig:
            if jenis == "konstan":
                s, A = putaran_konstan(par, D_KONSTAN)
            else:
                s, A = dari_sayap(sayap_per_lam[par], acak=(jenis == "acak"))
            x = deret(s, A)
            D = lama_tinggal_rata(s)
            for m in M_LIST:
                xs, ys, wb, pk, tepi = kurva_dan_puncak(x, m, U[m])
                baris += [dict(jenis=jenis, label=label, realisasi=r, m=m, W_per_P=round(w, 6), I3_rata=v)
                          for w, v in zip(xs, ys)]
                ringkas.append(dict(jenis=jenis, label=label, realisasi=r, m=m, c_m=wb, puncak=pk,
                                    di_tepi_grid=tepi, D_terukur=D, A_rata=A.mean(), A_sd=A.std()))
            print(f"realisasi {r + 1}/{N_REALISASI}  {jenis:8s} {label:9s}  D={D:.2f}  "
                  f"A rata-rata={A.mean():.2f}  ({time.time() - t0:.0f} detik)")
    return baris, ringkas, konfig


def simpan(nama, rows):
    with open(nama, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def stat(ringkas, jenis, label, kunci="c_m"):
    out = []
    for m in M_LIST:
        v = [d[kunci] for d in ringkas if d["jenis"] == jenis and d["label"] == label and d["m"] == m]
        out.append((np.mean(v), np.std(v)))
    return out


if __name__ == "__main__":
    baris, ringkas, konfig = jalankan()
    simpan("fase3b_hasil.csv", baris)
    simpan("fase3b_ringkasan.csv", ringkas)
    lor = lorenz_acuan()

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    warna = plt.cm.viridis(np.linspace(0, 0.9, len(AC_LIST)))
    for c, a in zip(warna, AC_LIST):
        mu, sd = zip(*stat(ringkas, "konstan", f"a/c={a:.2f}"))
        axs[0].errorbar(M_LIST, mu, yerr=sd, fmt="o-", color=c, capsize=3, label=f"a/c = {a:.2f}")
    axs[0].set_title(f"Eksperimen 1: amplitudo konstan (D = {D_KONSTAN})")
    warna = plt.cm.plasma(np.linspace(0, 0.8, len(LAMBDA_LIST)))
    for c, lam in zip(warna, LAMBDA_LIST):
        for jenis, gaya in (("spiral", "o-"), ("acak", "x:")):
            mu, sd = zip(*stat(ringkas, jenis, f"λ={lam}"))
            axs[1].errorbar(M_LIST, mu, yerr=sd, fmt=gaya, color=c, capsize=3, label=f"{jenis} λ={lam}")
    axs[1].set_title("Eksperimen 2: spiral (garis) vs urutan acak (titik-titik)")
    for ax in axs:
        if lor:
            ax.plot(M_LIST, lor, "ks--", lw=2, label="Lorenz ρ=28 asli")
        ax.set_xlabel("m")
        ax.set_xticks(M_LIST)
        ax.legend(fontsize=7, ncol=2)
    axs[0].set_ylabel("c_m = W*/P")
    plt.tight_layout()
    plt.savefig("fase3b_cm.png", dpi=120)

    print("\nc_m = W*/P (rata-rata ± sd antar realisasi);  kemiringan = c_6 - c_3")
    if lor:
        print("  Lorenz ρ=28 asli:   " + "  ".join(f"m{m}={v:.3f}" for m, v in zip(M_LIST, lor))
              + f"   kemiringan={lor[4] - lor[1]:+.3f}")
    for jenis, label, _ in konfig:
        st = stat(ringkas, jenis, label)
        D = np.mean([d["D_terukur"] for d in ringkas if d["jenis"] == jenis and d["label"] == label])
        tepi = any(d["di_tepi_grid"] for d in ringkas if d["jenis"] == jenis and d["label"] == label)
        print(f"  {jenis:8s} {label:9s} D={D:4.2f}  " + "  ".join(f"m{m}={mu:.3f}" for m, (mu, _) in zip(M_LIST, st))
              + f"   kemiringan={st[4][0] - st[1][0]:+.3f}" + ("   (ada puncak di tepi grid!)" if tepi else ""))
    print("\nTersimpan: fase3b_hasil.csv, fase3b_ringkasan.csv, fase3b_cm.png")
