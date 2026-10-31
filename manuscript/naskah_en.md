---
lang: en-US
---

::: {custom-style="PaperTitle"}
Projection-index landscapes of delay-embedded chaotic attractors: period scaling and its mechanism
:::

::: {custom-style="AuthorName"}
Athaillah Farrel Inandra
:::

::: {custom-style="AuthorAffiliation"}
Universitas Pelita Harapan, Tangerang, Indonesia

Corresponding author. Email: athaillahinandra@gmail.com · ORCID: 0009-0005-5729-5457
:::

::: {custom-style="AbstractHeading"}
Abstract
:::

::: {custom-style="AbstractBody"}
Principal component analysis and projection pursuit seek a single most informative direction in high-dimensional data. We reverse this viewpoint: a scalar time series is lifted to a higher dimension by delay embedding, and the whole landscape of a projection index over the sphere of directions becomes the object of study. Averaging Sarle's bimodality coefficient over directions as a function of the window length $W = (m-1)\tau$, we find that for the Lorenz system at $\rho$ = 28, 35 and 45 the bimodality-maximizing window obeys $W^{∗} \approx c_m P$, where $P$ is the rotation period and $c_m$ rises from 0.48 ($m = 2$) to 0.75 ($m = 6$), varying by at most 0.014 across $\rho$. The pattern is absent from phase-randomized surrogates and from the Hénon map, and fails for the Rössler system. A two-wing oscillator model explains the mechanism: same-wing delay vectors project onto a single shifted sinusoid whose bimodality depends on one ratio, with minimum 7/15, and even rare wing switches depress the landscape at $W = P$ through fourth-moment outliers. A ladder of models attributes 77 to 84 percent of the remaining gap to the second harmonic of the in-wing waveform and to rotation-to-rotation variability.
:::

::: {custom-style="KeywordLine"}
Keywords: delay embedding; projection pursuit; bimodality coefficient; Lorenz attractor; surrogate data; embedding window

2020 Mathematics Subject Classification: Primary 37M10; Secondary 37D45, 37M05, 62H25
:::

# 1. Introduction

Principal component analysis (PCA) projects high-dimensional data onto the few directions of largest variance (Pearson 1901). Projection pursuit generalizes this idea: every unit direction $u$ receives a score $I(u)$ that measures how "interesting" the projected data are, usually how far their distribution departs from a Gaussian, and the direction with the highest score is sought (Friedman and Tukey 1974; Huber 1985; Jones and Sibson 1987). The Gaussian serves as the uninteresting reference because most low-dimensional projections of high-dimensional data are close to Gaussian (Diaconis and Freedman 1984). In both frameworks the target is a single peak of $I(u)$.

This paper starts from the opposite question. Instead of reducing the dimension, a scalar time series $x(t)$ is lifted to $\mathbb{R}^m$ by delay embedding (Packard et al. 1980), and the whole landscape $I(u)$ over the sphere of directions $S^{m-1}$ is studied, not only its peak. A linear lift would add no information, because the data would remain in a subspace of the original dimension. Delay embedding is different because each coordinate carries the series at a different time. The theorem of Takens (1981) and its extensions (Sauer, Yorke, and Casdagli 1991) guarantee that for $m$ large enough the attractor is reconstructed from a single scalar observable up to a diffeomorphism. Since the reconstruction is defined only up to such a change of coordinates, a natural question arises for any statistic computed on it: which features of the landscape are properties of the attractor, and which are artifacts of the embedding coordinates $\tau$ and $m$?

The ingredients of this idea are individually familiar. Combining delay embedding with PCA gives singular spectrum analysis (Broomhead and King 1986), and PCA preceded by a nonlinear lift is kernel PCA (Schölkopf, Smola, and Müller 1998). Because every projection of a delay vector is the output of a finite impulse response (FIR) filter, maximizing non-Gaussianity over directions is equivalent to single-channel independent component analysis (Davies and James 2007). Index functions have been studied as objects in their own right for the purpose of optimization, through their smoothness and squintability (Laa and Cook 2020; Zhang et al. 2026). A function over all directions in $S^{m-1}$ is also at the heart of the persistent homology transform (Turner, Mukherjee, and Boyer 2014), although the quantity measured there is topological rather than distributional. On the reconstruction side, standard rules choose $\tau$ and $m$ separately (Fraser and Swinney 1986; Kennel, Brown, and Abarbanel 1992), whereas Kugiumtzis (1996) argued that the meaningful parameter is the window length $\tau_w = (m-1)\tau$ and recommended a window of at least one mean orbital period for estimating the correlation dimension.

To the author's knowledge, no study has treated the landscape of a distribution-shape index on a chaotic attractor as an object that depends on the embedding parameters. The paper makes three contributions.

1. An empirical scaling law $W^{∗} \approx c_m P$ for the Lorenz family: the peak of the direction-averaged landscape is set by the rotation period $P$, with a constant $c_m$ that depends on $m$ but hardly on $\rho$ (Section 3).
2. Closed-form landscapes for a two-wing oscillator model, including the effect of wing switching (Sections 4.1 to 4.4).
3. A ladder of models that separates the contributions of waveform shape, rotation-to-rotation variability and amplitude growth to the rise of $c_m$ with $m$ (Section 4.5).

Section 2 describes the framework and methods. Section 5 discusses the results in relation to previous work, reproducibility and limitations, and Section 6 closes with open questions. Proofs are collected in Appendix A and numerical details in Appendix B.

# 2. Framework and methods

## 2.1. Delay embedding as a bank of filters

For a scalar series $x(t)$, delay $\tau$ and dimension $m$, the delay vector and its window length are

::: {custom-style="Equation"}
$$v(t) = \left(x(t),\, x(t-\tau),\, \dots,\, x(t-(m-1)\tau)\right)^\top \in \mathbb{R}^m, \qquad W = (m-1)\tau.$$
:::

Each coordinate is normalized to zero mean and unit standard deviation. For a unit direction $u = (u_0, \dots, u_{m-1})^\top \in S^{m-1}$, the projection is

::: {custom-style="Equation"}
$$y_u(t) = u^\top v(t) = \sum_{k=0}^{m-1} u_k\, x(t-k\tau).$$
:::

Equation (2) shows that every direction is an FIR filter of unit norm and length $W$. Since $y_{-u} = -y_u$ and the indices used here are invariant under reflection, directions need only be considered modulo sign.

## 2.2. Bimodality index and the mean landscape

The main index is Sarle's bimodality coefficient with its finite-sample correction (Pfister et al. 2013),

::: {custom-style="Equation"}
$$I_3(u) = \frac{\gamma_u^2 + 1}{\kappa_u + \dfrac{3(N-1)^2}{(N-2)(N-3)}},$$
:::

where $\gamma_u$ is the skewness and $\kappa_u$ the excess kurtosis of $y_u$, and $N$ is the number of points. Its reference values are 1/3 for a Gaussian, 5/9 for the uniform distribution (the conventional bimodality threshold), 2/3 for the arcsine distribution (a pure sinusoid) and 1 for a two-point distribution. For a symmetric distribution and large $N$, $I_3 \to \mu_2^2/\mu_4$, where $\mu_k$ is the $k$th central moment. As a second check we use the negentropy approximation $I_1(u) = \left(\mathbb{E}[\log\cosh y_u] - \mathbb{E}[\log\cosh \nu]\right)^2$ with $\nu \sim N(0,1)$ (Hyvärinen 1998).

The landscape is summarized by its average over directions,

::: {custom-style="Equation"}
$$\overline{I}_3(W) = \mathbb{E}_{u \sim \operatorname{Unif}(S^{m-1})}\left[I_3(u)\right], \qquad W^{∗} = \arg\max_W \overline{I}_3(W), \qquad c_m = \frac{W^{∗}}{P}.$$
:::

The average is estimated with 500 to 600 random directions drawn uniformly on the sphere. The peak location $W^{∗}$ is taken from the grid maximum and refined by a parabola through the three surrounding points; peaks at the edge of the grid are flagged and not interpreted. This continuous summary was chosen after two alternatives, the fraction of directions with $I_3 > 5/9$ and the number of high-valued "islands" of directions, proved dependent on the threshold and less sensitive.

## 2.3. Systems, period and null model

Three systems are used: the Lorenz (1963) system with $\sigma = 10$, $\beta = 8/3$ and $\rho \in \{28, 35, 45\}$; the Rössler (1976) system with its standard parameter values 0.2, 0.2 and 5.7; and the Hénon (1976) map with parameters 1.4 and 0.3. The flows are integrated with an explicit Runge-Kutta (4,5) scheme with relative and absolute tolerances of $10^{-8}$ to $10^{-9}$ (Appendix B). Time in both flows is dimensionless and all durations below are in these natural time units. Only the $x$ component is used to build the landscape.

The period $P$ is defined as the median interval between successive local maxima of $z(t)$ for the Lorenz system and of $x(t)$ for the Rössler system; for the Hénon map, $P$ is one iteration. $P$ is deliberately measured from the full state, because estimators based on $x(t)$ alone are unreliable: maxima of $x$ miss about a quarter of the Lorenz rotations, while maxima of $|x|$ count every Rössler rotation twice.

The null model is phase-randomized surrogate data (Theiler et al. 1992). The power spectrum is preserved and the Fourier phases are randomized, so the linear autocorrelation is unchanged while nonlinear structure is destroyed. One surrogate is generated for each realization, and its landscape is computed in the same way.

## 2.4. Sequence of experiments

The experiments were run in stages. Phase 1 tests the stability of the peak location for the Lorenz system at $\rho = 28$ with respect to initial conditions and data length. Phase 2 tests whether the peak follows the period across three values of $\rho$, and for the Rössler and Hénon systems. Phases 3 to 3c use toy models to explain the mechanism. Full settings are listed in Appendix B.

# 3. Empirical results

## 3.1. The landscape of a single embedding

Figure 1 shows the Lorenz attractor ($\rho = 28$), the observed series $x(t)$, its reconstruction for $m = 3$ and $\tau = 0.12$, and the landscape $I_3(u)$ on $S^2$. The landscape is far from flat. It ranges from 0.19 to 0.75, and about 25 percent of the sphere has $I_3 > 5/9$. High-valued regions come in antipodal pairs because $I(u) = I(-u)$.

In the direction of the maximum the projection is bimodal, with two sharp spikes at its edges (Figure 1e). Seen from that direction the attractor appears edge-on as an S-shaped curve, and the spikes arise at the folds, where many points pile up at one projected value. In the direction of the minimum the projection is sharply peaked at the center with heavy tails (Figure 1f).

Different indices agree on where there is no structure, but not on where structure is strongest. The Spearman rank correlation between the $I_1$ and $I_3$ landscapes is 0.84, yet their maximizing directions are 41 degrees apart. The analysis that follows therefore uses the direction-averaged summary and always checks both indices.

![Figure 1. (a) Lorenz attractor with $\rho = 28$. (b) The series $x(t)$, the only data used. (c) Delay reconstruction with $m = 3$ and $\tau = 0.12$. (d) Landscape $I_3(u)$ on $S^2$ in spherical coordinates; the black line is the contour $I_3 = 5/9$, the star marks the maximum and the circle the minimum. (e, f) Distribution of the projection in the directions of the maximum and the minimum, compared with a Gaussian.](gambar_en/fig1_landscape.png){width=16cm}

## 3.2. Period scaling in the Lorenz family

At $\rho = 28$ the peak location of $\overline{I}_3(W)$ is stable with respect to data length and initial conditions. Over 10 initial conditions and two series lengths (75 and 150), $W^{∗}$ differs by at most 0.008 between lengths, on a grid with step 0.08. The peak of $\overline{I}_1$ falls at the same $W$, within one grid step, for every $m$ from 2 to 7. Surrogates give $\overline{I}_3$ between 0.31 and 0.37 (mean 0.345), close to the Gaussian value 1/3, with no peak.

To test whether the peak follows the period, the experiment was repeated at $\rho$ = 28, 35 and 45 with five realizations per system and with the $W$ axis normalized by the respective period. The periods are $P$ = 0.724 ± 0.016, 0.627 ± 0.005 and 0.540 ± 0.007. Figure 2 shows that the curves for the three values of $\rho$ nearly collapse after normalization, and Table 1 summarizes the peak locations.

For $m = 2$ and $m = 3$, the values of $c_m$ differ by at most 0.002 across $\rho$. For $m = 4$ to 6 the differences are 0.008 to 0.014 and decrease consistently as $\rho$ increases, slightly exceeding the standard deviation across realizations (0.002 to 0.006). The scaling law

$$W^{∗} \approx c_m P$$

therefore holds nearly, but not exactly. The constant $c_m$ itself increases with $m$, from 0.48 at $m = 2$ to 0.75 at $m = 6$.

: Table 1. Values of $c_m = W^{∗}/P$ (mean ± standard deviation over five realizations) in the Lorenz family. Differences quoted in the text are computed before rounding.

| $\rho$ | *P* | $m = 2$ | $m = 3$ | $m = 4$ | $m = 5$ | $m = 6$ |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 28 | 0.724 ± 0.016 | 0.478 ± 0.004 | 0.626 ± 0.005 | 0.660 ± 0.003 | 0.686 ± 0.004 | 0.748 ± 0.006 |
| 35 | 0.627 ± 0.005 | 0.479 ± 0.003 | 0.627 ± 0.004 | 0.655 ± 0.003 | 0.681 ± 0.003 | 0.744 ± 0.003 |
| 45 | 0.540 ± 0.007 | 0.479 ± 0.002 | 0.625 ± 0.003 | 0.647 ± 0.002 | 0.676 ± 0.003 | 0.741 ± 0.004 |

For $m = 2$ the curve also has a sharp valley at *W*/*P* = 1.03 ± 0.02, 1.07 ± 0.01 and 1.05 ± 0.02 for the three values of $\rho$. At $\rho = 28$ it drops to 0.344, practically the surrogate level. The valley is no longer visible at $m = 4$ (Figure 2b).

The peak location is not recovered by standard rules for choosing $\tau$. The first minimum of the mutual information (Fraser and Swinney 1986) lies at $\tau$ = 0.164 ± 0.005, which for $m = 3$ gives $W$ = 0.33, whereas the landscape peaks at $W^{∗} = c_3 P$ = 0.45. The first zero of the autocorrelation (2.4 ± 0.8) is not even stable across realizations.

![Figure 2. Direction-averaged $I_3$ as a function of $W/P$ for (a) $m = 2$ and (b) $m = 4$ at three values of $\rho$. Bands show ±1 standard deviation across realizations, the dashed gray line is the surrogate ($\rho = 28$), and the horizontal line is the Gaussian value 1/3. (c) $c_m = W^{∗}/P$ against $m$, with error bars of ±1 standard deviation over five realizations.](gambar_en/fig2_scaling.png){width=16cm}

## 3.3. Limits of scope: the Rössler and Hénon systems

For the Rössler system ($P$ = 5.96), $\overline{I}_3(W)$ oscillates with several peaks of nearly equal height. As a result the location of the global maximum jumps between values of $m$: $W^{∗}/P$ equals 0.80, 1.67, 1.19 and 1.63 for $m = 2$ to 5, and falls at the edge of the grid for $m = 6$. The Rössler spectrum is dominated by a single narrow peak, so its surrogate retains an almost sinusoidal oscillation and gives $\overline{I}_3$ = 0.45, well above 1/3. For systems of this kind the comparison must be made against their own surrogates.

For the Hénon map the maximum always lies at the smallest window, one iteration per coordinate, with no interior peak. This is expected for a system without rotation. The scaling law is therefore claimed only for attractors with a two-wing geometry like that of the Lorenz system.

# 4. Mechanism: a two-wing oscillator model

## 4.1. The model

Two features of the Lorenz dynamics are visible in $x(t)$: oscillation within a wing and irregular switching between wings. The following model retains both without chaos:

::: {custom-style="Equation"}
$$x(t) = s(t)\left[c + a\cos\theta(t)\right], \qquad \theta(t) = \frac{2\pi t}{P}, \qquad s(t) \in \{-1, +1\}.$$
:::

The sign $s(t)$ changes only at the phase $\theta = \pi$, where $|x|$ is minimal, and the dwell time in a wing is geometrically distributed with mean $D$ rotations. The condition $a/c < 1$ ensures that $x$ does not cross zero within a wing. At $\rho = 28$ the Lorenz trajectory switches wings roughly once every 2.2 rotations, so $D = 2.2$ is used as the Lorenz-like value. In this model the delay between coordinates becomes a phase shift $\psi = 2\pi\tau/P$. A pair or a vector whose coordinates all lie in the same wing is called same-wing.

## 4.2. Projections of same-wing vectors

::: {custom-style="TheoremHeading"}
Proposition 1 (form of the projection).
:::

::: {custom-style="TheoremBody"}
If all coordinates of $v(t)$ lie in the same wing, then for every direction $u \in S^{m-1}$
:::

::: {custom-style="Equation"}
$$y_u = s\left[C + R\cos(\theta - \delta)\right], \qquad C = c\sum_{k=0}^{m-1} u_k, \qquad R = a\,|Z|, \qquad \delta = -\arg Z, \qquad Z = \sum_{k=0}^{m-1} u_k\, e^{-ik\psi}.$$
:::

The proof (Appendix A.2) uses only $\cos(\theta - k\psi) = \operatorname{Re}\, e^{i(\theta - k\psi)}$: a sum of sinusoids of equal frequency is a single sinusoid. The proposition holds for any $m$, and all dependence on the direction $u$ and the delay $\tau$ is collected in two numbers, $C$ and $R$.

## 4.3. The bimodality coefficient as a function of one variable

::: {custom-style="TheoremHeading"}
Proposition 2 (bimodality function).
:::

::: {custom-style="TheoremBody"}
Let $\theta$ be uniform on [0, 2π), let $s = \pm 1$ with equal probability independently of $\theta$, and let $C \neq 0$. With $\eta = R^2/C^2$, the projection in Proposition 1 satisfies $\mathbb{E}[y] = 0$, $\mu_2 = C^2 + R^2/2$ and $\mu_4 = C^4 + 3C^2R^2 + 3R^4/8$, so that
:::

::: {custom-style="Equation"}
$$b(\eta) = \frac{\mu_2^2}{\mu_4} = \frac{2(2+\eta)^2}{3\eta^2 + 24\eta + 8}.$$
:::

::: {custom-style="TheoremBody"}
The function $b$ decreases from 1 (as $\eta \to 0$) to its minimum 7/15 at $\eta = 8/3$ and then increases towards 2/3 (as $\eta \to \infty$), since
:::

::: {custom-style="Equation"}
$$b'(\eta) = \frac{8(\eta+2)(3\eta-8)}{(3\eta^2 + 24\eta + 8)^2}.$$
:::

Both limits are meaningful. The limit $\eta \to 0$ is the two-point distribution $\pm C$, which is perfectly bimodal, whereas $\eta \to \infty$ is a pure sinusoid with an arcsine distribution. Figure 3a plots $b(\eta)$ together with the reference values.

::: {custom-style="TheoremHeading"}
Corollary 1 (landscape for *m* = 2 without switching).
:::

::: {custom-style="TheoremBody"}
For $m = 2$ and $u = (\cos\alpha, \sin\alpha)$, one has $C^2 = c^2(1 + \sin 2\alpha)$ and $R^2 = a^2(1 + \sin 2\alpha \cos\psi)$ with $\psi = 2\pi W/P$, so that
:::

::: {custom-style="Equation"}
$$\eta(\alpha, W) = \left(\frac{a}{c}\right)^2 \frac{1 + \sin 2\alpha\, \cos(2\pi W/P)}{1 + \sin 2\alpha}.$$
:::

::: {custom-style="TheoremBody"}
Consequently the mean landscape $L(W) = \mathbb{E}_\alpha\left[b(\eta(\alpha, W))\right]$ depends on $W$ only through $\cos(2\pi W/P)$. It is periodic with period $P$, symmetric about $W = P/2$ and stationary at $W = P/2$.
:::

At $W = P$, and also as $W \to 0$, one has $\eta = (a/c)^2$ for every direction, so that $L = b\left((a/c)^2\right)$. At $W = P/2$, $\eta$ sweeps all of [0, ∞) as $\alpha$ varies. Numerically, for $0 < a/c < 1$ the point $W = P/2$ is the global maximum of the landscape on $(0, P)$ if and only if $a/c > 0.4573$. Below this threshold the maximum moves to the ends $W \to 0$ and $W \to P$ (Figure 3b).

The geometric explanation is simple. At $W = P$, and as $W \to 0$, the two coordinates of a same-wing pair are equal, so the points lie on the diagonal and every projection is merely a rescaled copy of $x$ itself. At $W = P/2$ a same-wing pair has the form $(c + a\cos\theta,\; c - a\cos\theta)$, a segment perpendicular to the diagonal, one segment per wing. Directions near the diagonal collapse each segment almost to a point and see two tight clusters, directions near the antidiagonal see a single arcsine distribution, and directions in between pass through the minimum 7/15. If $a/c$ is small, $x$ itself is already close to a two-point distribution, and its rescaled copy at $W = P$ is hard to beat. If $a/c$ is large, each wing of $x$ is spread out and $b\left((a/c)^2\right)$ falls towards 18/35 ≈ 0.51 as *a*/*c* → 1, so rotating the segments at $W = P/2$ gives many directions a more bimodal view.

Simulations support all three results. The model with rare switching (*D* = 40, *a*/*c* = 0.6) matches the formula for $W \le P/2$ to within 0.008 (Figure 3c). A pure sinusoid gives 0.6657 for every $W$, which is the expected value of 2/3 after Sarle's finite-sample correction for about 4000 points. The threshold in $a/c$ is also consistent with simulations that include switching: with *D* = 2.2, the values *a*/*c* = 0.30 and 0.45 produce no interior peak for any $m$, whereas *a*/*c* ≥ 0.60 produces a peak near $P/2$ for $m = 2$, with *c*~2~ between 0.45 and 0.49.

## 4.4. Wing switching and the valley at *W* = *P*

The no-switch formula has no valley at $W = P$ (black curve in Figure 3c), whereas the simulations and the Lorenz system do. The valley comes from pairs that straddle a wing switch.

::: {custom-style="TheoremHeading"}
Proposition 3 (mixture at *W* = *P*).
:::

::: {custom-style="TheoremBody"}
For $m = 2$, an opposite-wing pair has the form of Proposition 1 with $u_1$ replaced by $-u_1$, that is, $C' = c(u_0 - u_1)$ and $R'^2 = a^2(1 - \sin 2\alpha\cos\psi)$. At $W = P$ every pair straddles exactly one rotation boundary, so the proportion of same-wing pairs is $p = 1 - 1/D$, independent of $\theta$. With $q = \sin 2\alpha$, the moments of the two kinds of pairs mix linearly and
:::

::: {custom-style="Equation"}
$$b_{\text{mix}}(\alpha) = b\left((a/c)^2\right)\, \frac{\left[1 + (2p-1)q\right]^2}{1 + 2(2p-1)q + q^2}.$$
:::

::: {custom-style="TheoremHeading"}
Corollary 2 (directions without same-wing signal).
:::

::: {custom-style="TheoremBody"}
In the two diagonal directions,
:::

::: {custom-style="Equation"}
$$b_{\text{mix}}(3\pi/4) = (1-p)\, b\left((a/c)^2\right), \qquad b_{\text{mix}}(\pi/4) = p\, b\left((a/c)^2\right).$$
:::

Corollary 2 explains why even rare switching has a large effect. In the direction $\alpha = 3\pi/4$ all same-wing pairs project exactly onto zero because $C = R = 0$. The only variance comes from the rare opposite-wing pairs, which become outliers and dominate the fourth moment. The rarer the switching, the smaller $b$ becomes in that direction, and the average over $\alpha$ recovers its no-switch value only very slowly. For *a*/*c* = 0.6 the no-switch value is 0.654, while the mixture average is 0.483, 0.512 and 0.578 for *p* = 0.9, 0.975 and 0.999. A single switch per thousand rotations still lowers the mean landscape by more than 10 percent.

With the switching proportions actually realized in the simulations (0.4268 for *D* = 2.2 and 0.0235 for *D* = 40), equation (10) gives 0.4629 and 0.5131. Direct computation on the simulated series with 4000 uniform directions and all points gives 0.4633 and 0.5187. The remaining difference for *D* = 40 is unsurprising, since that simulation contains only about 47 switches in 2000 rotations. The Lorenz valley near $W = P$ (Section 3.2) is therefore the signature of wing switching in the landscape.

![Figure 3. (a) The function $b(\eta)$ of Proposition 2 with the reference values for the two-point (1), sinusoidal (2/3), uniform (5/9) and Gaussian (1/3) distributions. (b) The $m = 2$ landscape without switching: mean $b$ at $W = P/2$ and at $W \to 0$ or $W = P$ as a function of $a/c$; $P/2$ is the global maximum when *a*/*c* > 0.457. (c) The no-switch formula (black) compared with simulations of the two-wing model (*a*/*c* = 0.6) for *D* = 40 and *D* = 2.2. Diamonds are the predictions of Proposition 3 at $W = P$ with the realized switching proportions.](gambar_en/fig3_analytic.png){width=16cm}

## 4.5. The rise of *c*~*m*~ with *m*

For $m \ge 3$, Propositions 1 and 2 still hold, and the symmetry argument of Corollary 1 carries over. Because $u$ is real, $|Z|$ is invariant under $\psi \to -\psi$, so the no-switch landscape is symmetric about $\tau = P/2$ and stationary there, that is, at $W = (m-1)P/2$ (Appendix A.4). Numerically this point is the global maximum for *a*/*c* between 0.6 and 1.0 at $m = 3$ and $m = 4$, so the no-switch formula predicts *c*~3~ = 1 and *c*~4~ = 1.5. The sinusoidal model with Lorenz-like switching instead peaks at 0.43*P* to 0.54*P* for $m = 3$ and *a*/*c* between 0.6 and 1.0 (Table 2 and model A in Table 3), while the Lorenz system gives *c*~3~ = 0.626. For $m \ge 3$, wing switching therefore determines the peak location rather than merely correcting it, and $c_m$ cannot yet be derived in closed form. The rise of $c_m$ is instead studied with a ladder of models, each of which adds one directly measured feature of the Lorenz dynamics.

The amplitude level affects the slope. In the sinusoidal model with constant amplitude, the slope $c_6 - c_3$ rises from 0.035 to 0.069 as *a*/*c* increases from 0.6 to 0.9, but $c_3$ stays below 0.53, far from the Lorenz value 0.626 (Table 2).

The growing-spiral hypothesis is not supported. In the Lorenz system the oscillation amplitude grows within a wing before the trajectory is ejected to the other wing. The spiral model mimics this with $A_{j+1} = \lambda A_j$, starting from a random amplitude between 0.3 and 0.5 and continuing until a critical amplitude of 0.9 triggers a switch. As a control, the order of the amplitudes within each wing is shuffled while the wing sequence and the set of amplitudes are preserved. Relative to the control, the spiral lowers $c_3$ by 0.009 to 0.015, away from the Lorenz value, and changes $c_6$ by no more than 0.004. The slightly larger slope of the spiral model comes from a lower $c_3$, not from a higher $c_6$.

: Table 2. Peak locations for the constant-amplitude model and for the spiral model with its shuffled control (means over three realizations). Column *D* is the measured mean dwell time in rotations. The column $c_6 - c_3$ and the differences quoted in the text are computed before rounding.

| Model | Parameter | *D* | $c_3$ | $c_6$ | $c_6 - c_3$ |
|:--|:--:|:--:|:--:|:--:|:--:|
| Constant amplitude | *a*/*c* = 0.60 | 2.21 | 0.426 | 0.461 | +0.035 |
| Constant amplitude | *a*/*c* = 0.75 | 2.19 | 0.488 | 0.539 | +0.051 |
| Constant amplitude | *a*/*c* = 0.90 | 2.29 | 0.528 | 0.598 | +0.069 |
| Spiral | *λ* = 1.25 | 5.17 | 0.475 | 0.537 | +0.063 |
| Shuffled control | *λ* = 1.25 | 5.17 | 0.489 | 0.536 | +0.047 |
| Spiral | *λ* = 1.6 | 3.26 | 0.479 | 0.538 | +0.059 |
| Shuffled control | *λ* = 1.6 | 3.26 | 0.490 | 0.534 | +0.044 |
| Spiral | *λ* = 2.0 | 2.76 | 0.484 | 0.541 | +0.056 |
| Shuffled control | *λ* = 2.0 | 2.76 | 0.494 | 0.539 | +0.045 |
| Lorenz | *ρ* = 28 | about 2.2 | 0.626 | 0.748 | +0.122 |

The in-wing waveform was measured directly from Lorenz trajectories. Each rotation that does not switch wings is cut from one minimum of $|x|$ to the next and normalized as $g = (|x| - c)/a$, where $c$ is the mean and $a$ the half-range of $|x|$ over that rotation (Figure 4a). Over five trajectories with different initial conditions, each contributing about 813 same-wing rotations, the second harmonic of the mean shape is 0.171 ± 0.003 of the fundamental, the third harmonic is 0.021 ± 0.001, and the median $a/c$ is 1.00 ± 0.01, with an interquartile range of about 0.75 to 1.21. The asymmetry is vertical rather than temporal: the peaks of $|x|$ are higher and sharper, and the troughs wider and shallower. A temporal asymmetry could not affect the mean landscape in any case, because time reversal only permutes the embedding coordinates, and the average over uniformly distributed directions is invariant under such a permutation.

The model ladder (Table 3 and Figure 4b) uses the same switching (*D* = 2.2) and the median Lorenz value of $a/c$, and each rung adds one feature. Model A uses a cosine; model B uses the Lorenz shape truncated to harmonics 1 and 2; model C uses the full mean shape; model D adds a random $a/c$ per rotation; and model E uses the shape and $a/c$ of randomly drawn real Lorenz rotations.

The second harmonic is the main carrier. Model B closes 61 percent of the gap in $c_3$ between model A and the Lorenz system, and 43 to 52 percent for $m = 4$ to 6. Model C, which contains all harmonics, does no better than B. Rotation-to-rotation variability adds to the effect at larger $m$: model E closes 77 to 84 percent of the gap for every $m$ from 3 to 6.

Model D cannot be used. In the Lorenz system, rotations with large $a/c$ have shallower normalized troughs (correlation 0.62 ± 0.02 between $a/c$ and the minimum of $g$). Pairing the mean shape with a random $a/c$ makes about 41 percent of the rotations cross zero, so those values must be discarded and the model becomes biased. Shape and amplitude are tied together, and only model E preserves that tie.

The residual gap between model E and the Lorenz system is 0.018 to 0.025. It estimates the combined effect of two features absent from every model: the actual switching trajectory, which the models replace by an instantaneous jump, and the variability of the period between rotations. Two independent runs with different Lorenz trajectories gave *c*~6~ = 0.728 and 0.714 for model E, so the residual is known only roughly.

: Table 3. The model ladder: $c_m = W^{∗}/P$ (means over three realizations; the largest standard deviation across realizations is 0.008). Model D is invalid because it ignores the coupling between shape and amplitude (see text). The column $c_6 - c_3$ and the percentages and residuals quoted in the text are computed before rounding.

| Model | $c_2$ | $c_3$ | $c_4$ | $c_5$ | $c_6$ | $c_6 - c_3$ |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| A, cosine | 0.495 | 0.543 | 0.549 | 0.585 | 0.619 | +0.076 |
| B, harmonics 1 and 2 | 0.491 | 0.594 | 0.606 | 0.629 | 0.676 | +0.083 |
| C, mean shape | 0.488 | 0.587 | 0.601 | 0.626 | 0.673 | +0.086 |
| D, plus random *a*/*c* | 0.468 | 0.566 | 0.567 | 0.603 | 0.638 | +0.071 |
| E, real rotations | 0.461 | 0.608 | 0.635 | 0.669 | 0.728 | +0.120 |
| Lorenz, *ρ* = 28 | 0.478 | 0.626 | 0.660 | 0.686 | 0.748 | +0.122 |

![Figure 4. (a) In-wing waveform: 40 real Lorenz rotations (gray), the mean shape (black, model C), its truncation to harmonics 1 and 2 (orange, model B) and a cosine (blue, model A). The small jump at phase 0.5 is the rotation boundary, where two successive minima differ in depth. (b) $c_m$ against $m$ for models A, B and E (mean ± standard deviation over three realizations) and for the Lorenz system at $\rho = 28$ (five realizations). Model C nearly coincides with B, and model D is invalid; both are listed in Table 3.](gambar_en/fig4_ladder.png){width=16cm}

# 5. Discussion

## 5.1. Relation to previous work

The results agree with the argument of Kugiumtzis (1996) that the window length, rather than $\tau$ and $m$ separately, is the meaningful parameter, and that the window is tied to the mean orbital period. The direction of the effect, however, differs. Kugiumtzis recommended $\tau_w \ge \tau_p$ as a lower bound for correlation-dimension estimation, whereas the bimodality landscape here peaks at $W^{∗} = c_m P$ with $c_m$ between 0.48 and 0.75, before one full period. The window that maximizes distributional structure in projections is shorter than the window recommended for dimension estimation, because the two criteria answer different questions.

From a dynamical-systems viewpoint, the landscape is a coordinate-dependent observable of the reconstruction, yet its features are tied to distinct ingredients of the flow. The peak location is set by an intrinsic time scale, the mean rotation period, together with the in-wing waveform, while the valley at $W = P$ is set by the switching statistics through Proposition 3. The landscape therefore reflects, in separate features, the oscillation within a wing and the symbolic sequence of wings that standard embedding diagnostics treat together.

The filter view of equation (2) places this work close to single-channel independent component analysis (Davies and James 2007), which seeks one filter $u$ of maximal non-Gaussianity. The object here is different: the average over all filters as a function of window length, rather than the optimal filter. Likewise, studies of index landscapes in projection pursuit (Laa and Cook 2020; Zhang et al. 2026) use smoothness and squintability to assess how easily an index can be optimized, not how the landscape depends on the embedding parameters of a dynamical system. The warning of Bickel, Kur, and Nadler (2018) that projection pursuit can find spurious structure in high-dimensional Gaussian data is also relevant; averaging over directions and comparing with surrogates are used precisely to avoid that trap.

## 5.2. Reproducibility

Results involving only toy models are reproducible almost exactly, because all random numbers come from fixed seeds. Runs of Phases 3 and 3b on two machines with different operating systems (Windows and Linux) gave index values that differ by at most $3 \times 10^{-16}$, the level of floating-point round-off. Results that integrate the Lorenz equations are reproducible only statistically. Round-off differences suffice to separate trajectories: perturbing the initial condition by $10^{-14}$ changes the number of same-wing rotations from 824 to 830, and a run on another machine gave 817. All Lorenz quantities are therefore reported as means and standard deviations over several trajectories, and the uncertainty of model E is estimated from two runs with different trajectories.

## 5.3. Limitations

Several limitations should be noted.

- The scaling law was tested on a single family of attractors, and the Rössler results show that it is not universal.
- Only one family of indices was used, namely the shape of one-dimensional projection distributions, with negentropy as a check.
- The resolution of the $W$ grid limits the precision of $W^{∗}$ before interpolation. The Phase 1 grid is coarse (0.08), so part of the zero spread across realizations there is a grid artifact.
- Proposition 3 is exact only at $W = P$. At other windows the probability that a pair straddles a switch depends on the phase $\theta$.
- The values of $c_m$ for $m \ge 3$ have not been derived analytically, and the small drift of $c_m$ with $\rho$ for $m \ge 4$ remains unexplained.
- The model treats a wing switch as an instantaneous jump, whereas in the Lorenz system a switch is a long excursion that passes near zero.

# 6. Conclusions

The bimodality landscape of the delay-embedded Lorenz attractor has structure governed by the dynamics rather than by the choice of coordinates: its direction-averaged peak follows $W^{∗} \approx c_m P$ at three values of $\rho$. A two-wing oscillator model explains its three main features. The peak at $P/2$ for $m = 2$ arises from in-wing geometry and is the global maximum only when *a*/*c* > 0.457; the valley at $W = P$ arises from wing switching, and its mixture formula matches simulation; and the rise of $c_m$ with $m$ comes mainly from the second harmonic of the waveform and from rotation-to-rotation variability.

Four questions remain open.

1. A closed-form landscape with wing switching for arbitrary $W$ and $m \ge 3$, which would yield $c_m$ analytically.
2. Qualitative changes of the landscape as $m$ crosses the Takens bound. The low-dimensional projection argument of Diaconis and Freedman (1984) suggests that such changes are gradual.
3. The map of disagreement between indices as an object in its own right.
4. Other families of attractors, and experimental data.

# Acknowledgements

The author used Claude, a large language model developed by Anthropic, accessed in October 2026, as a research and writing assistant. It was used to draft simulation code that the author ran and reproduced, to check the derivations numerically and symbolically, to extend the symmetry argument to arbitrary $m$ in Appendix A.4, to search the literature, to write the plotting scripts, and to draft and translate the manuscript text. All figures were produced by deterministic Python and Matplotlib scripts from the simulation output; no generative AI tool was used to create or alter images. The author formulated the research question and hypotheses, derived Propositions 1 to 3, ran all simulations, revised the text, and takes full responsibility for the content of the article.

# Disclosure statement

No potential conflict of interest was reported by the author.

# Funding

The author received no specific funding for this work.

# Data availability statement

The Python scripts for all experiments (invarians_lanskap.py, fase1_W.py, fase2_sistem.py, fase3_model_mainan.py, fase3b_spiral.py and fase3c_gelombang.py), the plotting scripts and the CSV outputs will be deposited in a public repository before publication [repository DOI to be added], and are available from the author on request in the meantime.

# References

::: {custom-style="ReferenceEntry"}
Bickel, P. J., G. Kur, and B. Nadler. 2018. "Projection Pursuit in High Dimensions." *Proceedings of the National Academy of Sciences* 115 (37): 9151-9156. https://doi.org/10.1073/pnas.1801177115.

Broomhead, D. S., and G. P. King. 1986. "Extracting Qualitative Dynamics from Experimental Data." *Physica D* 20 (2-3): 217-236.

Davies, M. E., and C. J. James. 2007. "Source Separation Using Single Channel ICA." *Signal Processing* 87 (8): 1819-1832. https://doi.org/10.1016/j.sigpro.2007.01.011.

Diaconis, P., and D. Freedman. 1984. "Asymptotics of Graphical Projection Pursuit." *The Annals of Statistics* 12 (3): 793-815.

Fraser, A. M., and H. L. Swinney. 1986. "Independent Coordinates for Strange Attractors from Mutual Information." *Physical Review A* 33 (2): 1134-1140.

Friedman, J. H., and J. W. Tukey. 1974. "A Projection Pursuit Algorithm for Exploratory Data Analysis." *IEEE Transactions on Computers* C-23 (9): 881-890.

Hénon, M. 1976. "A Two-Dimensional Mapping with a Strange Attractor." *Communications in Mathematical Physics* 50 (1): 69-77.

Huber, P. J. 1985. "Projection Pursuit." *The Annals of Statistics* 13 (2): 435-475.

Hyvärinen, A. 1998. "New Approximations of Differential Entropy for Independent Component Analysis and Projection Pursuit." In *Advances in Neural Information Processing Systems 10*, edited by M. I. Jordan, M. J. Kearns, and S. A. Solla, 273-279. Cambridge, MA: MIT Press.

Jones, M. C., and R. Sibson. 1987. "What Is Projection Pursuit?" *Journal of the Royal Statistical Society: Series A* 150 (1): 1-36.

Kennel, M. B., R. Brown, and H. D. I. Abarbanel. 1992. "Determining Embedding Dimension for Phase-Space Reconstruction Using a Geometrical Construction." *Physical Review A* 45 (6): 3403-3411.

Kugiumtzis, D. 1996. "State Space Reconstruction Parameters in the Analysis of Chaotic Time Series: The Role of the Time Window Length." *Physica D* 95 (1): 13-28. https://doi.org/10.1016/0167-2789(96)00054-1.

Laa, U., and D. Cook. 2020. "Using Tours to Visually Investigate Properties of New Projection Pursuit Indexes with Application to Problems in Physics." *Computational Statistics* 35: 1171-1205. https://doi.org/10.1007/s00180-020-00954-8.

Lorenz, E. N. 1963. "Deterministic Nonperiodic Flow." *Journal of the Atmospheric Sciences* 20 (2): 130-141.

Packard, N. H., J. P. Crutchfield, J. D. Farmer, and R. S. Shaw. 1980. "Geometry from a Time Series." *Physical Review Letters* 45 (9): 712-716.

Pearson, K. 1901. "On Lines and Planes of Closest Fit to Systems of Points in Space." *Philosophical Magazine* 2 (11): 559-572.

Pfister, R., K. A. Schwarz, M. Janczyk, R. Dale, and J. B. Freeman. 2013. "Good Things Peak in Pairs: A Note on the Bimodality Coefficient." *Frontiers in Psychology* 4: 700. https://doi.org/10.3389/fpsyg.2013.00700.

Rössler, O. E. 1976. "An Equation for Continuous Chaos." *Physics Letters A* 57 (5): 397-398.

Sauer, T., J. A. Yorke, and M. Casdagli. 1991. "Embedology." *Journal of Statistical Physics* 65 (3-4): 579-616.

Schölkopf, B., A. Smola, and K.-R. Müller. 1998. "Nonlinear Component Analysis as a Kernel Eigenvalue Problem." *Neural Computation* 10 (5): 1299-1319.

Takens, F. 1981. "Detecting Strange Attractors in Turbulence." In *Dynamical Systems and Turbulence, Warwick 1980*, edited by D. Rand and L.-S. Young, 366-381. Lecture Notes in Mathematics 898. Berlin: Springer.

Theiler, J., S. Eubank, A. Longtin, B. Galdrikian, and J. D. Farmer. 1992. "Testing for Nonlinearity in Time Series: The Method of Surrogate Data." *Physica D* 58 (1-4): 77-94.

Turner, K., S. Mukherjee, and D. M. Boyer. 2014. "Persistent Homology Transform for Modeling Shapes and Surfaces." *Information and Inference* 3 (4): 310-344. https://doi.org/10.1093/imaiai/iau011.

Zhang, H. S., D. Cook, N. Langrené, and J. W. Y. Leung. 2026. "Squintability and Other Metrics for Assessing Projection Pursuit Indexes, and Guiding Optimization Choices." *Journal of Computational and Graphical Statistics*. https://doi.org/10.1080/10618600.2026.2614076.
:::

# Appendix A. Proofs

This appendix collects the derivations of Propositions 1 to 3 and Corollaries 1 and 2. Coordinates are indexed from zero, so that $u = (u_0, \dots, u_{m-1})$ and coordinate $k$ is $x(t - k\tau)$.

## A.1. Moments of a uniform phase

For $\theta$ uniform on [0, 2π), the identities $\cos^2\theta = (1 + \cos 2\theta)/2$, $\cos^3\theta = (3\cos\theta + \cos 3\theta)/4$ and $\cos^4\theta = (3 + 4\cos 2\theta + \cos 4\theta)/8$ give

::: {custom-style="EquationA"}
$$\mathbb{E}[\cos\theta] = \mathbb{E}[\cos^3\theta] = 0, \qquad \mathbb{E}[\cos^2\theta] = \frac{1}{2}, \qquad \mathbb{E}[\cos^4\theta] = \frac{3}{8}.$$
:::

Since $\theta - \delta$ (modulo $2\pi$) is also uniform, the same moments hold for $X = \cos(\theta - \delta)$ for any $\delta$.

## A.2. Proof of Proposition 1

If all coordinates lie in the same wing, coordinate $k$ equals $s\,[c + a\cos(\theta - k\psi)]$. Hence

::: {custom-style="EquationA"}
$$y_u = \sum_{k=0}^{m-1} u_k\, s\left[c + a\cos(\theta - k\psi)\right] = s\left[c\sum_{k=0}^{m-1} u_k + a\sum_{k=0}^{m-1} u_k \cos(\theta - k\psi)\right].$$
:::

Using $\cos x = \operatorname{Re}\, e^{ix}$,

::: {custom-style="EquationA"}
$$\sum_{k=0}^{m-1} u_k \cos(\theta - k\psi) = \operatorname{Re}\left[e^{i\theta} \sum_{k=0}^{m-1} u_k e^{-ik\psi}\right] = \operatorname{Re}\left[e^{i\theta} Z\right].$$
:::

Write $Z = |Z|\, e^{i\phi}$. Then $\operatorname{Re}[e^{i\theta} Z] = |Z| \cos(\theta + \phi) = |Z| \cos(\theta - \delta)$ with $\delta = -\phi = -\arg Z$. With $C = c\sum_k u_k$ and $R = a|Z|$ we obtain $y_u = s[C + R\cos(\theta - \delta)]$. ∎

For $m = 2$ and $u = (\cos\alpha, \sin\alpha)$ one has $C^2 = c^2(\cos\alpha + \sin\alpha)^2 = c^2(1 + \sin 2\alpha)$ and

::: {custom-style="EquationA"}
$$|Z|^2 = (\cos\alpha + \sin\alpha\cos\psi)^2 + (\sin\alpha\sin\psi)^2 = 1 + \sin 2\alpha\cos\psi.$$
:::

## A.3. Proof of Proposition 2

Write $X = \cos(\theta - \delta)$. Since $s$ is independent of $\theta$, $\mathbb{E}[s] = 0$ and $s^2 = s^4 = 1$, we have $\mathbb{E}[y] = \mathbb{E}[s]\,\mathbb{E}[C + RX] = 0$. By A.1 the odd powers of $X$ vanish in expectation, so that $\mu_2 = \mathbb{E}[(C + RX)^2] = C^2 + R^2/2$ and

::: {custom-style="EquationA"}
$$\mu_4 = \mathbb{E}\left[(C + RX)^4\right] = C^4 + 6C^2R^2\cdot\frac{1}{2} + R^4\cdot\frac{3}{8} = C^4 + 3C^2R^2 + \frac{3R^4}{8}.$$
:::

The distribution of $y$ is symmetric, so its skewness vanishes, and for large $N$ the correction term in equation (3) tends to 3, so that $I_3 \to 1/(\kappa + 3) = \mu_2^2/\mu_4$. With $R^2 = \eta C^2$,

::: {custom-style="EquationA"}
$$b(\eta) = \frac{C^4(1 + \eta/2)^2}{C^4(1 + 3\eta + 3\eta^2/8)} = \frac{2(2+\eta)^2}{3\eta^2 + 24\eta + 8}.$$
:::

With numerator $n(\eta) = (1 + \eta/2)^2$ and denominator $d(\eta) = 1 + 3\eta + 3\eta^2/8$,

::: {custom-style="EquationA"}
$$n'd - nd' = \left(1 + \frac{\eta}{2}\right)\left[d - \left(1 + \frac{\eta}{2}\right)\left(3 + \frac{3\eta}{4}\right)\right] = \left(1 + \frac{\eta}{2}\right)\left(\frac{3\eta}{4} - 2\right),$$
:::

which is equivalent to equation (8). For $\eta \ge 0$ the factor $\eta + 2$ is positive, so $b'$ has the sign of $3\eta - 8$: $b$ decreases on [0, 8/3) and increases afterwards. Hence $\eta = 8/3$ is the global minimum on [0, ∞), with $b(8/3) = 2(14/3)^2/(64/3 + 72) = 7/15$. The limits are $b(0) = 1$ and $b(\eta) \to 2/3$ as $\eta \to \infty$. ∎

## A.4. Proof of Corollary 1

From A.2, $\eta = R^2/C^2 = (a/c)^2 (1 + \sin 2\alpha\cos\psi)/(1 + \sin 2\alpha)$ with $\psi = 2\pi W/P$ for $m = 2$, so $L(W)$ depends on $W$ only through $\cos\psi$. Its derivative with respect to $W$ contains the factor $\sin(2\pi W/P)$, which vanishes at $W = P/2$, so this point is stationary. At $W = P$, $\cos\psi = 1$ and $\eta = (a/c)^2$ for every $\alpha$. At $W = P/2$, $\cos\psi = -1$ and

::: {custom-style="EquationA"}
$$\eta = \left(\frac{a}{c}\right)^2 \frac{1 - q}{1 + q}, \qquad q = \sin 2\alpha.$$
:::

The value $q$ lies in (−1, 1], since $q = -1$ gives $C = 0$, which Proposition 2 excludes. For every $k \ge 0$ the value $q = (1 - k)/(1 + k)$ lies in this interval and gives $\eta = k\,(a/c)^2$, so $\eta$ sweeps all of [0, ∞). The threshold 0.4573 was obtained by solving $L(P/2) = b\left((a/c)^2\right)$ numerically (Brent's method, averaging over $2 \times 10^5$ values of $\alpha$). A sweep of *a*/*c* from 0.30 to 1.00 found no other local maximum in $(0, P/2)$. For *a*/*c* well above 1, outside the valid range of the model, the ordering changes again.

The symmetry argument holds for any $m$. By Proposition 1, $C$ does not depend on $\tau$, and since $u$ is real, $Z(-\psi) = \overline{Z(\psi)}$, so $|Z|$ is invariant under $\psi \to -\psi$. Together with $2\pi$-periodicity in $\psi$, the no-switch landscape satisfies $L(\tau) = L(P - \tau)$ and is therefore stationary at $\tau = P/2$, that is, at $W = (m-1)P/2$. A numerical sweep with 40,000 random directions shows that this point is the global maximum for *a*/*c* between 0.6 and 1.0 at $m = 3$ and $m = 4$.

## A.5. Proof of Proposition 3

For an opposite-wing pair, $v = s\left(c + a\cos\theta,\; -[c + a\cos(\theta - \psi)]\right)^\top$, so that

::: {custom-style="EquationA"}
$$y' = s\left[c(u_0 - u_1) + a\left(u_0\cos\theta - u_1\cos(\theta - \psi)\right)\right].$$
:::

This is the form of Proposition 1 with $u_1$ replaced by $-u_1$. Step A.2 then gives $C'^2 = c^2(1 - \sin 2\alpha)$ and $R'^2 = a^2(1 - \sin 2\alpha\cos\psi)$. At $W = P$, $\psi = 2\pi$ and $\cos\psi = 1$. With $q = \sin 2\alpha$, $K_2 = c^2 + a^2/2$ and $K_4 = c^4 + 3c^2a^2 + 3a^4/8$, Proposition 2 gives

::: {custom-style="EquationA"}
$$\mu_2^{\text{same}} = K_2(1+q), \qquad \mu_4^{\text{same}} = K_4(1+q)^2, \qquad \mu_2^{\text{opp}} = K_2(1-q), \qquad \mu_4^{\text{opp}} = K_4(1-q)^2,$$
:::

where the labels refer to same-wing and opposite-wing pairs. The pair $(x(t), x(t-P))$ always straddles exactly one rotation boundary, and a switch occurs at that boundary with probability $1/D$, independently of $\theta$. Both components have zero mean, so the central moments of the mixture are weighted averages. For the second moment, $\mu_2 = p\,\mu_2^{\text{same}} + (1-p)\,\mu_2^{\text{opp}} = K_2\left[1 + (2p-1)q\right]$, and for the fourth

::: {custom-style="EquationA"}
$$\mu_4 = K_4\left[p(1+q)^2 + (1-p)(1-q)^2\right] = K_4\left[1 + 2(2p-1)q + q^2\right].$$
:::

The ratio $\mu_2^2/\mu_4$ gives equation (10), since $K_2^2/K_4 = b\left((a/c)^2\right)$ by Proposition 2. For $q = -1$ the numerator is $[1 - (2p-1)]^2 = 4(1-p)^2$ and the denominator is $4(1-p)$, so that $b_{\text{mix}} = (1-p)\,b\left((a/c)^2\right)$. For $q = 1$ the same computation gives $b_{\text{mix}} = p\,b\left((a/c)^2\right)$. ∎

For example, with *D* = 40 one has *p* = 39/40, so that $b_{\text{mix}}(3\pi/4) = b\left((a/c)^2\right)/40$, which for *a*/*c* = 0.6 equals 0.654/40 ≈ 0.0164.

# Appendix B. Numerical details

: Table B1. Settings of each phase. The number of random directions is 600 in Phase 1 and 500 in the other phases.

| Phase | System or model | Data per realization | Window grid | Realizations | Seed |
|:--|:--|:--|:--|:--:|:--:|
| 1 | Lorenz, *ρ* = 28 | 75 and 150 time units, step 0.01 | *W* = 0.08 to 1.28, step 0.08 | 10 | 2026 |
| 2 | Lorenz (three *ρ*), Rössler, Hénon | 20,000 samples; step 0.005 (Lorenz), 0.05 (Rössler) | *W*/*P* = 0.06 to 2.0, step 0.06; Hénon 1 to 24 iterations | 5 | 2027 |
| 3 | Two-wing, *a*/*c* = 0.6 | 200,000 (*D* = 40) and 60,000 (*D* = 2.2) samples, 100 per rotation | *W*/*P* = 0.02 to 2.0, step 0.02 | 1 | 31 |
| 3b | Constant amplitude, spiral, shuffled control | 60,000 samples | *W*/*P* = 0.03 to 1.6, step 0.03 | 3 | 2718 |
| 3c | Model ladder A to E | 60,000 samples | *W*/*P* = 0.03 to 1.6, step 0.03 | 3 | 1618 |

From Phase 2 onwards the index is computed on about 4000 points per projection (a regular subsample); in Phase 1 every third point is used. Random directions are generated as normalized standard Gaussian vectors, so they are uniformly distributed on the sphere, and the same set of directions is used for every $\tau$ within an experiment. The discarded transients are 20 (Phase 1), 30 (Lorenz, Phase 2) and 300 (Rössler) time units, and 1000 iterations for the Hénon map. Lorenz initial conditions are drawn from a normal distribution with mean (1, 1, 20) and standard deviation 5, Rössler initial conditions from mean (1, 1, 0) and standard deviation 0.5, and Hénon initial conditions from mean 0.1 and standard deviation 0.05 in each coordinate. Integration uses `solve_ivp` (RK45) from SciPy with tolerances $10^{-8}$ in Phase 1 and $10^{-9}$ elsewhere. Skewness and excess kurtosis are computed from uncorrected sample moments; for about 4000 points the difference from bias-corrected versions in equation (3) is negligible.

In Phase 1 the peak is the maximum after a three-point moving average. From Phase 2 onwards it is located by three-point parabolic interpolation, and peaks at the edge of the grid are flagged. The waveform in Phase 3c is measured on a single Lorenz trajectory at $\rho = 28$ of length 1000 time units (step 0.005) starting from (1, 1, 20), with 200 phase points per rotation. The waveform statistics in Section 4.5 are computed from five trajectories of that length: this one and four with random initial conditions. All computations use Python with NumPy, SciPy and Matplotlib.
