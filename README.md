# Projection-index landscapes of delay-embedded chaotic attractors

Code, data and manuscript sources for

> Inandra, A. F. *Projection-index landscapes of delay-embedded chaotic attractors: period scaling and its mechanism.* Manuscript submitted to *Dynamical Systems: An International Journal*.

A scalar time series is lifted to $\mathbb{R}^m$ by delay embedding, and the whole landscape of a projection index (mainly Sarle's bimodality coefficient) over the sphere of directions is studied as a function of the window length $W = (m-1)\tau$. For the Lorenz family the direction-averaged landscape peaks at $W^* \approx c_m P$, where $P$ is the rotation period, and a two-wing oscillator model explains the mechanism.

Author: Athaillah Farrel Inandra, Universitas Pelita Harapan (athaillahinandra@gmail.com, ORCID 0009-0005-5729-5457).

## Repository layout

```
experiments/   simulation scripts, one per phase of the study
data/          CSV output of every phase (the numbers behind all tables and figures)
analysis/      figure scripts, statistics and verification scripts
manuscript/    manuscript sources (Indonesian and English), build scripts, figures, .docx drafts
```

Code comments and variable names are in Indonesian. A short glossary: *fase* = phase, *hasil* = raw results, *ringkasan* = summary, *gambar* = figure, *naskah* = manuscript, *cek* = check, *bangun* = build, *lanskap* = landscape, *sayap* = wing, *putaran* = rotation, *realisasi* = realization, *arah* = direction, *benih* = random seed.

## Phases, outputs and where they appear in the paper

| Script (`experiments/`) | Output (`data/`) | Paper |
|---|---|---|
| `invarians_lanskap.py` | `hasil_invarians.csv` | exploratory invariance run that motivated the window-length summary (Section 2.2) |
| `fase1_W.py` | `fase1_hasil.csv`, `fase1_ringkasan.csv` | stability of the peak at $\rho = 28$ (Section 3.2) |
| `fase2_sistem.py` | `fase2_hasil.csv`, `fase2_ringkasan.csv` | period scaling, Rössler and Hénon (Sections 3.2 and 3.3, Table 1, Figure 2) |
| `fase3_model_mainan.py` | `fase3_hasil.csv` | two-wing model, no-switch formula and the valley at $W = P$ (Sections 4.3 and 4.4, Figure 3c) |
| `fase3b_spiral.py` | `fase3b_hasil.csv`, `fase3b_ringkasan.csv` | constant amplitude and growing-spiral models with shuffled control (Section 4.5, Table 2) |
| `fase3c_gelombang.py` | `fase3c_hasil.csv`, `fase3c_ringkasan.csv` | measured Lorenz waveform and the model ladder A to E (Section 4.5, Table 3, Figure 4) |

`fase2_sistem.py` also provides the shared routines (embedding, random directions, index computation, parabolic peak location) that the later phases import. Settings and random seeds for every phase are listed in Appendix B of the manuscript.

## Data provenance

The CSV files for the invariance run and for Phases 1, 2, 3 and 3c were produced on the author's machine (Windows). The Phase 3b files come from a complete rerun on Linux; on the rows both runs share, the two machines agree to within $3 \times 10^{-16}$. All toy-model results depend only on fixed seeds and reproduce to floating-point round-off. Results that integrate the Lorenz equations reproduce only statistically, because round-off differences are enough to separate chaotic trajectories (Section 5.2 of the manuscript).

## Reproducing the results

Install the requirements (Python 3.10 or newer):

```
pip install -r requirements.txt
```

Run the experiments from inside `experiments/`. Each script writes its CSV files to the current folder; copy them to `data/` if you want the figures to use your own run. Approximate run times on a laptop: Phase 1, 20 to 40 minutes; Phase 3, 5 to 10 minutes; Phase 3b, 10 to 20 minutes; Phase 3c, about 10 minutes. Phase 2 is the longest.

```
cd experiments
python fase1_W.py
python fase2_sistem.py
python fase3_model_mainan.py
python fase3b_spiral.py
python fase3c_gelombang.py
```

The figure, statistics and verification scripts (written for the Indonesian draft, so they still call the ratio $R^2/C^2$ `rho`; the English manuscript writes it as $\eta$) can be run from any folder; they read `data/` and write into `manuscript/`.

```
python analysis/gambar1.py      # Figure 1, Indonesian labels  -> manuscript/gambar/
python analysis/fig_en_1.py     # Figure 1, English labels     -> manuscript/gambar_en/
                                # likewise gambar2-4.py and fig_en_2-4.py
python analysis/statistik.py    # the numbers quoted in Sections 3 and 4, recomputed from data/
python analysis/cek_angka.py    # table values, rounding, gap percentages and residuals
python analysis/cek_simbolik.py # SymPy: phase moments, Proposition 2, m = 2 case of Proposition 1, Proposition 3 and Corollary 2
python analysis/cek_numerik.py  # numerical check of Corollary 1 and the threshold a/c = 0.4573
python analysis/cek_campuran.py # Proposition 3 against direct computation on the simulated series
python analysis/cek_m3_halus.py # no-switch peak at tau = P/2 for m = 3 and 4 (Appendix A.4)
python analysis/cek_sitasi.py   # every citation has a reference entry and vice versa
```

`analysis/cek_m3.py` is a coarser version of `cek_m3_halus.py`. `analysis/gaya.py` holds the shared plotting style.

## Building the manuscript

The manuscripts are written in Pandoc Markdown and converted to Word with native, editable equations. This needs [pandoc](https://pandoc.org) (version 3.1.3 was used) and `python-docx`.

```
cd manuscript
python bangun.py      # naskah.md    -> naskah.docx    (Indonesian draft)
python bangun_en.py   # naskah_en.md -> naskah_en.docx (English submission draft)
```

The `.docx` files in `manuscript/` are the drafts built from these sources.

## Software versions used

Python 3.13, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, SymPy 1.14.0, scikit-learn (for `invarians_lanskap.py` only), python-docx 1.2.0 and pandoc 3.1.3. The author's own runs used the versions installed on his machine, which may differ.

## License

To be added by the author before the repository is made public.
