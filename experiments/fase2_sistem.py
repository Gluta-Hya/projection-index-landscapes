"""
FASE 2: Apakah W*/P ~ 0.7 berlaku di sistem lain?

Hipotesis dari Fase 1 (Lorenz rho=28):
  (H1) Untuk m >= 3, panjang jendela W = (m-1)*tau yang memaksimalkan
       rata-rata I3 atas semua arah berada di sekitar W* ~ 0.7 P,
       dengan P = periode satu putaran osilasi.
  (H2) Untuk m = 2, rata-rata I3 punya lembah tajam di W ~ P.

Sistem yang diuji:
  lorenz28, lorenz35, lorenz45 : sistem sama, periode berbeda
  rossler                      : bentuk attractor berbeda, satu osilasi dominan
  henon                        : peta diskret tanpa putaran (kontrol negatif)

Definisi P = periode satu putaran, diukur dari KEADAAN PENUH sistem
(bukan dari x saja), memakai cara baku untuk tiap sistem:
  Lorenz : median selang antar maksimum lokal z(t)  (dasar "peta Lorenz")
  Rössler: median selang antar maksimum lokal x(t)  (x berputar sekali per putaran)
  Hénon  : tidak punya putaran; P = 1 langkah, sehingga W/P = W dalam langkah.
Catatan: penaksir dari x saja ternyata tidak seragam. Puncak x melewatkan
sekitar sepertiga putaran Lorenz di sayap kiri, sedangkan puncak |x|
menghitung putaran Rössler dua kali. Karena itu P diambil dari keadaan penuh.
Analisis lanskap tetap HANYA memakai x(t).

Keluaran:
  fase2_hasil.csv      semua angka mentah
  fase2_ringkasan.csv  W*/P per (sistem, realisasi, m) dan lembah m=2
  fase2_kurva.png      rata-rata I3 dan I1 vs W/P, satu baris per sistem
  fase2_Wbintang.png   W*/P per sistem dan m, dengan galat antar realisasi
  Ringkasan dicetak di terminal.

Jalankan:  python fase2_sistem.py
Kebutuhan: numpy scipy matplotlib (sudah terpasang)
"""
import csv
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

RNG = np.random.default_rng(2027)

# ---------------- pengaturan ----------------
N_REALISASI = 5
N_SAMPEL = 20_000            # panjang deret per realisasi (setelah transien)
M_LIST = (2, 3, 4, 5, 6)
K_ARAH = 500                  # banyak arah acak per m
N_INDEKS = 4_000              # banyak titik yang dipakai menghitung indeks
W_PER_P_MAKS = 2.0            # grid W/P untuk sistem kontinu
W_PER_P_LANGKAH = 0.06
HENON_W_MAKS = 24             # grid W (langkah) untuk Hénon
AMBANG_BIMODAL = 0.555


# ---------------- sistem ----------------
def _flow(f, ic, dt, buang):
    """Kembalikan seluruh keadaan (3, N_SAMPEL) setelah transien dibuang."""
    T = (N_SAMPEL - 1) * dt + buang
    t = np.arange(0, T + dt / 2, dt)
    sol = solve_ivp(f, (0, t[-1]), ic, t_eval=t, method="RK45", rtol=1e-9, atol=1e-9)
    if not sol.success:
        raise RuntimeError(sol.message)
    S = sol.y[:, int(round(buang / dt)):][:, :N_SAMPEL]
    if S.shape[1] != N_SAMPEL:
        raise RuntimeError("panjang deret tidak sesuai")
    return S


def _median_selang_puncak(sinyal, dt):
    pk, _ = find_peaks(sinyal)
    return np.median(np.diff(pk)) * dt


def lorenz(rho):
    def gen():
        f = lambda t, s: [10 * (s[1] - s[0]), s[0] * (rho - s[2]) - s[1], s[0] * s[1] - 8 / 3 * s[2]]
        S = _flow(f, RNG.normal([1, 1, 20], 5), dt=0.005, buang=30.0)
        return S[0], _median_selang_puncak(S[2], 0.005)      # x, P dari maksimum z
    return gen, 0.005


def rossler():
    def gen():
        f = lambda t, s: [-s[1] - s[2], s[0] + 0.2 * s[1], 0.2 + s[2] * (s[0] - 5.7)]
        S = _flow(f, RNG.normal([1, 1, 0], 0.5), dt=0.05, buang=300.0)
        return S[0], _median_selang_puncak(S[0], 0.05)       # x, P dari maksimum x
    return gen, 0.05


def henon():
    def gen():
        x, y = RNG.normal(0.1, 0.05), RNG.normal(0.1, 0.05)
        out = np.empty(N_SAMPEL)
        for i in range(1000 + N_SAMPEL):
            x, y = 1 - 1.4 * x * x + y, 0.3 * x
            if i >= 1000:
                out[i - 1000] = x
        return out, 1.0                                       # P = 1 langkah
    return gen, 1.0


SISTEM = {
    "lorenz28": lorenz(28.0),
    "lorenz35": lorenz(35.0),
    "lorenz45": lorenz(45.0),
    "rossler": rossler(),
    "henon": henon(),
}
DISKRET = {"henon"}


# ---------------- alat ----------------
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


def daftar_tau(m, dt, W_target):
    """tau (dalam sampel) unik, >= 1, untuk setiap target W."""
    t = np.unique(np.round(np.asarray(W_target) / ((m - 1) * dt)).astype(int))
    return t[t >= 1]


# ---------------- indeks ----------------
_Z = np.random.default_rng(0).standard_normal(400_000)
EGN = np.mean(np.log(np.cosh(_Z)))           # E[log cosh(nu)], nu ~ N(0,1)


def ringkas(V, U):
    sub = max(1, len(V) // N_INDEKS)
    Y = V[::sub] @ U.T
    Y = (Y - Y.mean(0)) / Y.std(0)
    n = Y.shape[0]
    # log cosh aman di sini: Y distandarkan, |Y| <= sqrt(n) jauh di bawah batas overflow (~710)
    I1 = (np.log(np.cosh(Y)).mean(0) - EGN) ** 2
    Y2 = Y * Y                                  # perkalian jauh lebih cepat daripada Y**3, Y**4
    g = (Y2 * Y).mean(0)
    k = (Y2 * Y2).mean(0) - 3
    I3 = (g ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    return dict(I3_rata=I3.mean(), I3_fraksi=(I3 > AMBANG_BIMODAL).mean(), I1_rata=I1.mean())


def puncak_parabola(xs, ys):
    """Lokasi maksimum dengan interpolasi parabola lewat 3 titik (grid tak seragam boleh)."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    i = int(np.argmax(ys))
    if i == 0 or i == len(xs) - 1:
        return xs[i], ys[i], True             # di tepi grid: tidak bisa dihaluskan
    a, b, c = np.polyfit(xs[i - 1:i + 2], ys[i - 1:i + 2], 2)
    if a >= 0:
        return xs[i], ys[i], False
    xv = -b / (2 * a)
    xv = min(max(xv, xs[i - 1]), xs[i + 1])
    return xv, a * xv * xv + b * xv + c, False


# ---------------- eksperimen ----------------
def jalankan():
    U_per_m = {m: arah(m, K_ARAH) for m in M_LIST}
    baris, info = [], {}
    t0 = time.time()
    for nama, (gen, dt) in SISTEM.items():
        P_list = []
        grid = None
        for r in range(N_REALISASI):
            x, P = gen()
            if not np.all(np.isfinite(x)):
                raise RuntimeError(f"{nama}: deret tidak hingga (lintasan lepas dari attractor)")
            P_list.append(P)
            if grid is None:   # grid dibangun sekali per sistem dari realisasi pertama
                if nama in DISKRET:
                    grid = np.arange(1, HENON_W_MAKS + 1) * dt
                else:
                    grid = np.arange(W_PER_P_LANGKAH, W_PER_P_MAKS + 1e-9, W_PER_P_LANGKAH) * P
            s = surrogate_acak_fase(x)
            for m in M_LIST:
                for tau in daftar_tau(m, dt, grid):
                    W = (m - 1) * tau * dt
                    for sumber, sig in (("asli", x), ("surrogate", s)):
                        d = ringkas(embed(sig, m, tau), U_per_m[m])
                        d.update(sistem=nama, realisasi=r, sumber=sumber, m=m, tau=int(tau),
                                 W=round(W, 6), P=round(P, 6), W_per_P=round(W / P, 6))
                        baris.append(d)
            print(f"{nama:9s} realisasi {r + 1}/{N_REALISASI}  P = {P:.4f}  "
                  f"({time.time() - t0:.0f} detik)")
        info[nama] = dict(P_rata=np.mean(P_list), P_sd=np.std(P_list), dt=dt)
    return baris, info


def ringkasan(baris):
    out = []
    for nama in SISTEM:
        for r in range(N_REALISASI):
            for m in M_LIST:
                d = sorted([b for b in baris if b["sistem"] == nama and b["realisasi"] == r
                            and b["m"] == m and b["sumber"] == "asli"], key=lambda b: b["W"])
                if len(d) < 3:
                    continue
                xs = [b["W_per_P"] for b in d]
                ys = [b["I3_rata"] for b in d]
                wb, pk, tepi = puncak_parabola(xs, ys)
                row = dict(sistem=nama, realisasi=r, m=m, P=d[0]["P"],
                           Wbintang_per_P=wb, Wbintang=wb * d[0]["P"], puncak=pk, di_tepi_grid=tepi,
                           lembah_per_P=np.nan, lembah_nilai=np.nan)
                if m == 2 and nama not in DISKRET:
                    # lembah: minimum I3_rata di 0.6 <= W/P <= 1.6
                    j = [i for i, x in enumerate(xs) if 0.6 <= x <= 1.6]
                    if j:
                        i = min(j, key=lambda i: ys[i])
                        row["lembah_per_P"], row["lembah_nilai"] = xs[i], ys[i]
                out.append(row)
    return out


def gambar_kurva(baris, info):
    warna = plt.cm.viridis(np.linspace(0, 0.9, len(M_LIST)))
    kolom = [("I3_rata", "rata-rata $I_3$"), ("I1_rata", "rata-rata $I_1$")]
    fig, axs = plt.subplots(len(SISTEM), 2, figsize=(14, 3.6 * len(SISTEM)), squeeze=False)
    for i, nama in enumerate(SISTEM):
        for j, (k, judul) in enumerate(kolom):
            ax = axs[i, j]
            for c, m in zip(warna, M_LIST):
                for sumber, gaya in (("asli", "-"), ("surrogate", ":")):
                    d = [b for b in baris if b["sistem"] == nama and b["m"] == m and b["sumber"] == sumber]
                    taus = sorted({b["tau"] for b in d})
                    xs = np.array([np.mean([b["W_per_P"] for b in d if b["tau"] == t]) for t in taus])
                    mu = np.array([np.mean([b[k] for b in d if b["tau"] == t]) for t in taus])
                    sd = np.array([np.std([b[k] for b in d if b["tau"] == t]) for t in taus])
                    ax.plot(xs, mu, gaya, color=c, lw=1.6 if sumber == "asli" else 1,
                            label=f"m={m}" if sumber == "asli" else None)
                    if sumber == "asli":
                        ax.fill_between(xs, mu - sd, mu + sd, color=c, alpha=0.2)
            if nama not in DISKRET:
                ax.axvline(0.7, color="r", ls="--", lw=0.9, label="W/P = 0.7" if j == 0 else None)
                ax.axvline(1.0, color="k", ls=":", lw=0.9, label="W/P = 1" if j == 0 else None)
            P = info[nama]["P_rata"]
            ax.set_title(f"{nama}: {judul}   (P = {P:.3f})")
            ax.set_xlabel("W (langkah)" if nama in DISKRET else "W / P")
            if j == 0:
                ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Garis tebal: data asli (pita = ±1 sd antar realisasi).  Titik-titik: surrogate.", y=1.0)
    plt.tight_layout()
    plt.savefig("fase2_kurva.png", dpi=110, bbox_inches="tight")


def gambar_wbintang(rk):
    fig, ax = plt.subplots(figsize=(9, 5))
    nama_list = [s for s in SISTEM if s not in DISKRET]   # Hénon tidak punya P: tidak digambar di sini
    lebar = 0.8 / len(M_LIST)
    warna = plt.cm.viridis(np.linspace(0, 0.9, len(M_LIST)))
    for j, (m, c) in enumerate(zip(M_LIST, warna)):
        mu = [np.mean([d["Wbintang_per_P"] for d in rk if d["sistem"] == s and d["m"] == m]) for s in nama_list]
        sd = [np.std([d["Wbintang_per_P"] for d in rk if d["sistem"] == s and d["m"] == m]) for s in nama_list]
        xpos = np.arange(len(nama_list)) - 0.4 + lebar * (j + 0.5)
        ax.errorbar(xpos, mu, yerr=sd, fmt="o", color=c, capsize=3, label=f"m={m}")
    ax.axhline(0.7, color="r", ls="--", lw=1, label="0.7")
    ax.set_xticks(np.arange(len(nama_list)))
    ax.set_xticklabels(nama_list)
    ax.set_ylabel("W* / P")
    ax.set_title("Lokasi puncak rata-rata $I_3$, dinormalisasi oleh periode P")
    ax.legend(fontsize=8, ncol=3)
    plt.tight_layout()
    plt.savefig("fase2_Wbintang.png", dpi=120)


def simpan(nama, rows):
    with open(nama, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    baris, info = jalankan()
    simpan("fase2_hasil.csv", baris)
    rk = ringkasan(baris)
    simpan("fase2_ringkasan.csv", rk)
    gambar_kurva(baris, info)
    gambar_wbintang(rk)

    print("\nRINGKASAN (rata-rata ± sd antar realisasi)")
    for nama in SISTEM:
        i = info[nama]
        print(f"\n{nama}:  P = {i['P_rata']:.4f} ± {i['P_sd']:.4f}   (dt = {i['dt']})")
        for m in M_LIST:
            d = [x for x in rk if x["sistem"] == nama and x["m"] == m]
            if not d:
                continue
            v = np.array([x["Wbintang_per_P"] for x in d])
            p = np.array([x["puncak"] for x in d])
            tepi = sum(x["di_tepi_grid"] for x in d)
            if tepi and nama in DISKRET:
                catatan = "   (puncak di tepi grid: tidak ada puncak di dalam, sesuai harapan kontrol)"
            elif tepi:
                catatan = f"   PERINGATAN: {tepi} realisasi puncaknya di tepi grid"
            else:
                catatan = ""
            satuan = "W* (langkah)" if nama in DISKRET else "W*/P"
            print(f"  m={m}:  {satuan} = {v.mean():.3f} ± {v.std():.3f}   puncak I3 = {p.mean():.3f}{catatan}")
        lem = np.array([x["lembah_per_P"] for x in rk if x["sistem"] == nama and x["m"] == 2])
        lem = lem[np.isfinite(lem)]
        if len(lem):
            print(f"  m=2 lembah di W/P = {lem.mean():.3f} ± {lem.std():.3f}")
        sur = [b["I3_rata"] for b in baris if b["sistem"] == nama and b["sumber"] == "surrogate"]
        print(f"  surrogate: rata-rata I3 = {np.mean(sur):.3f}  (rentang {np.min(sur):.3f} sampai {np.max(sur):.3f})")
    print("\nCatatan: untuk proses Gaussian, I3 = 1/3 = 0.333. Surrogate sistem yang hampir periodik")
    print("(misalnya Rössler) wajar berada di atas 0.333, karena sinyal hampir-sinus sendiri bimodal.")
    print("\nTersimpan: fase2_hasil.csv, fase2_ringkasan.csv, fase2_kurva.png, fase2_Wbintang.png")
