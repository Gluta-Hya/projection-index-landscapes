"""
FASE 3: Model mainan "osilator dua sayap" untuk menjelaskan mekanisme.

Model (P = 1, jadi W langsung bernilai W/P):
    x(t) = s(t) * ( c + a cos(2 pi t) ),     s(t) in {-1, +1}
s(t) berganti tanda hanya di fase theta = pi (titik x paling dekat nol),
dengan lama tinggal di satu sayap berdistribusi geometrik (rata-rata D putaran).

Tiga varian:
  sinus           : cos(2 pi t) saja                     -> harus datar di 2/3
  dua_sayap_jarang: D = 40 putaran (perpindahan langka)  -> harus cocok rumus analitik
  dua_sayap_lorenz: D = 2.2 putaran (seperti Lorenz rho=28)

Rumus analitik untuk m = 2 (mengabaikan pasangan yang melintasi perpindahan sayap):
    arah u = (cos al, sin al),  psi = 2 pi W/P
    rho(al) = (a/c)^2 * (1 + sin(2al) cos psi) / (1 + sin(2al))
    b(rho)  = 2 (2 + rho)^2 / (3 rho^2 + 24 rho + 8)
    kurva analitik = rata-rata b(rho(al)) atas al

Keluaran:
  fase3_m2.png       kurva m=2: analitik, model mainan, dan Lorenz asli (jika fase2_hasil.csv ada)
  fase3_Wbintang.png W*/P vs m: model mainan vs Lorenz
  fase3_hasil.csv    angka mentah
Jalankan:  python fase3_model_mainan.py      (perlu fase2_sistem.py di folder yang sama)
Waktu: sekitar 5 sampai 10 menit.
"""
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import fase2_sistem as F   # memakai embed, arah, ringkas yang SAMA dengan Fase 2

RNG = np.random.default_rng(31)
DT = 0.01                  # 100 sampel per putaran
A_PER_C = 0.6
M_LIST = (2, 3, 4, 5, 6)
W_GRID = np.round(np.arange(0.02, 2.001, 0.02), 3)
K_ARAH = 500


def dua_sayap(n_sampel, D, c=1.0, a=A_PER_C):
    n_putaran = int(n_sampel * DT) + 2
    s = np.empty(n_putaran)
    tanda = 1.0
    for i in range(n_putaran):
        s[i] = tanda
        if RNG.random() < 1.0 / D:
            tanda = -tanda
    theta = 2 * np.pi * np.arange(n_sampel) * DT
    putaran_ke = np.floor((theta + np.pi) / (2 * np.pi)).astype(int)   # ganti putaran di theta = pi
    return s[putaran_ke] * (c + a * np.cos(theta))


def analitik_m2(W, a_per_c=A_PER_C, n=20000):
    al = np.linspace(0, np.pi, n, endpoint=False)
    s2 = np.sin(2 * al)
    ok = (1 + s2) > 1e-9                     # buang satu arah degenerat (ukuran nol)
    rho = a_per_c ** 2 * (1 + s2[ok] * np.cos(2 * np.pi * W)) / (1 + s2[ok])
    return np.mean(2 * (2 + rho) ** 2 / (3 * rho ** 2 + 24 * rho + 8))


def kurva(x, m, U):
    hasil = []
    for W in W_GRID:
        tau = max(1, int(round(W / ((m - 1) * DT))))
        hasil.append((tau * (m - 1) * DT, F.ringkas(F.embed(x, m, tau), U)["I3_rata"]))
    W_nyata, y = map(np.array, zip(*hasil))
    return W_nyata, y


def lorenz_dari_fase2(path="fase2_hasil.csv", sistem="lorenz28"):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["sistem"] == sistem and r["sumber"] == "asli"]
    out = {}
    for m in M_LIST:
        d = [r for r in rows if int(r["m"]) == m]
        taus = sorted({int(r["tau"]) for r in d})
        xs = np.array([np.mean([float(r["W_per_P"]) for r in d if int(r["tau"]) == t]) for t in taus])
        ys = np.array([np.mean([float(r["I3_rata"]) for r in d if int(r["tau"]) == t]) for t in taus])
        out[m] = (xs, ys)
    return out


if __name__ == "__main__":
    U = {m: F.arah(m, K_ARAH) for m in M_LIST}
    model = {
        "sinus": np.cos(2 * np.pi * np.arange(60_000) * DT),
        "dua_sayap_jarang": dua_sayap(200_000, D=40),
        "dua_sayap_lorenz": dua_sayap(60_000, D=2.2),
    }
    baris, kurva_semua = [], {}
    for nama, x in model.items():
        for m in (M_LIST if nama == "dua_sayap_lorenz" else (2,)):
            Wn, y = kurva(x, m, U[m])
            kurva_semua[(nama, m)] = (Wn, y)
            baris += [dict(model=nama, m=m, W_per_P=w, I3_rata=v) for w, v in zip(Wn, y)]
            print(f"{nama:17s} m={m} selesai")
    an_W = np.linspace(0, 2, 401)
    an_y = np.array([analitik_m2(w) for w in an_W])
    lor = lorenz_dari_fase2()

    # ---------- uji kecocokan ----------
    Wn, y = kurva_semua[("sinus", 2)]
    print(f"\nSinus: I3 rata-rata = {y.mean():.4f} (teori 2/3 = 0.6667), rentang {y.min():.4f} sampai {y.max():.4f}")
    Wn, y = kurva_semua[("dua_sayap_jarang", 2)]
    awal = Wn <= 0.5
    beda = np.abs(y[awal] - np.interp(Wn[awal], an_W, an_y))
    print(f"Dua sayap jarang vs analitik, W/P <= 0.5: selisih maksimum {beda.max():.4f}")
    print(f"Analitik: maksimum {an_y.max():.4f} di W/P = {an_W[np.argmax(an_y[:201])]:.3f};  "
          f"nilai di W/P = 0 dan 1: {an_y[0]:.4f}, {an_y[200]:.4f}")

    # ---------- gambar 1: m = 2 ----------
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(an_W, an_y, "k-", lw=2, label="analitik (tanpa perpindahan sayap)")
    for nama, gaya in (("sinus", ":"), ("dua_sayap_jarang", "--"), ("dua_sayap_lorenz", "-")):
        Wn, y = kurva_semua[(nama, 2)]
        ax.plot(Wn, y, gaya, lw=1.6, label=f"model: {nama}")
    if lor:
        ax.plot(*lor[2], "-", color="purple", lw=2, alpha=0.7, label="Lorenz ρ=28 asli (Fase 2)")
    ax.axhline(1 / 3, color="gray", lw=0.8, ls=":", label="Gaussian (1/3)")
    for v in (0.5, 1.0, 1.5):
        ax.axvline(v, color="gray", lw=0.6, ls="--")
    ax.set_xlabel("W / P")
    ax.set_ylabel("rata-rata $I_3$ atas arah")
    ax.set_title("m = 2: rumus analitik, model mainan, dan Lorenz")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("fase3_m2.png", dpi=120)

    # ---------- gambar 2: W*/P vs m ----------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    wb_model = [kurva_semua[("dua_sayap_lorenz", m)][0][np.argmax(kurva_semua[("dua_sayap_lorenz", m)][1])]
                for m in M_LIST]
    ax.plot(M_LIST, wb_model, "o-", label="model dua sayap (D = 2.2)")
    if lor:
        wb_lor = [lor[m][0][np.argmax(lor[m][1])] for m in M_LIST]
        ax.plot(M_LIST, wb_lor, "s-", label="Lorenz ρ=28 asli (Fase 2, titik grid)")
    ax.set_xlabel("m")
    ax.set_ylabel("W* / P")
    ax.set_xticks(M_LIST)
    ax.set_title("Apakah model mainan mereproduksi kenaikan W*/P terhadap m?")
    ax.legend()
    plt.tight_layout()
    plt.savefig("fase3_Wbintang.png", dpi=120)

    print("\nW*/P per m:")
    for i, m in enumerate(M_LIST):
        extra = f"   Lorenz asli = {wb_lor[i]:.2f}" if lor else ""
        print(f"  m={m}: model = {wb_model[i]:.2f}{extra}")
    with open("fase3_hasil.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(baris[0].keys()))
        w.writeheader()
        w.writerows(baris)
    print("\nTersimpan: fase3_m2.png, fase3_Wbintang.png, fase3_hasil.csv")
