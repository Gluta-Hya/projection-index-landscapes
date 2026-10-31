---
lang: id-ID
---

::: {custom-style="JudulNaskah"}
Lanskap Indeks Proyeksi pada Atraktor Hasil *Embedding* Tunda: Skala Periode dan Mekanismenya
:::

::: {custom-style="JudulInggris"}
Projection-Index Landscapes of Delay-Embedded Chaotic Attractors: Period Scaling and Its Mechanism
:::

::: {custom-style="Penulis"}
Athaillah Farrel Inandra
:::

::: {custom-style="Afiliasi"}
Universitas Pelita Harapan, Tangerang, Indonesia

athaillahinandra@gmail.com · ORCID 0009-0005-5729-5457
:::

::: {custom-style="StatusDraf"}
Draf untuk umpan balik, 31 Oktober 2026
:::

::: {custom-style="JudulAbstrak"}
Abstrak
:::

::: {custom-style="Abstrak"}
Analisis komponen utama dan *projection pursuit* mencari satu arah pandang terbaik pada data berdimensi tinggi. Naskah ini membalik arah itu: deret waktu skalar diangkat ke ruang berdimensi lebih tinggi melalui *embedding* tunda, lalu seluruh lanskap indeks pola di atas bola arah dipelajari sebagai objek. Indeks yang dipakai adalah koefisien bimodalitas Sarle, dan lanskap diringkas sebagai rata-ratanya atas semua arah, sebagai fungsi panjang jendela $W = (m-1)\tau$. Pada sistem Lorenz dengan $\rho$ = 28, 35 dan 45, jendela yang memaksimalkan bimodalitas memenuhi $W^{∗} \approx c_m P$, dengan $P$ periode satu putaran. Konstanta $c_m$ naik dari 0,48 ($m = 2$) ke 0,75 ($m = 6$) dan berbeda paling banyak 0,014 antarnilai $\rho$, sementara $P$ berubah dari 0,72 s ke 0,54 s. Pola ini tidak muncul pada data *surrogate* acak-fase, tidak berlaku sebagai satu konstanta pada sistem Rössler, dan tidak ada pada peta Hénon. Mekanismenya dijelaskan oleh model osilator dua sayap. Tiga proposisi analitik menunjukkan bahwa proyeksi vektor tunda sesayap selalu berupa satu sinusoid bergeser, bahwa koefisien bimodalitasnya hanya bergantung pada $\rho = R^2/C^2$ dengan minimum 7/15 di $\rho = 8/3$, dan bahwa perpindahan sayap yang langka sekalipun menurunkan lanskap di $W = P$ melalui pencilan yang menguasai momen keempat. Rumus campuran di $W = P$ cocok dengan simulasi sampai 0,0004 untuk perpindahan ala Lorenz. Rangkaian model bertingkat menunjukkan bahwa harmonik kedua bentuk gelombang dan variasi antarputaran menutup 77 sampai 84 persen selisih $c_m$ antara model sinusoid dan Lorenz untuk $m$ = 3 sampai 6, sedangkan hipotesis pertumbuhan amplitudo spiral tereliminasi oleh kontrol acak.
:::

::: {custom-style="KataKunci"}
Kata kunci: *embedding* tunda; *projection pursuit*; koefisien bimodalitas; atraktor Lorenz; data *surrogate*; sistem dinamik
:::

::: {custom-style="JudulAbstrak"}
Abstract
:::

::: {custom-style="Abstrak"}
Principal component analysis and projection pursuit look for a single best viewing direction in high-dimensional data. This note reverses that direction: a scalar time series is lifted to a higher dimension by delay embedding, and the whole landscape of a pattern index over the sphere of directions is studied as an object. The index is Sarle's bimodality coefficient, and the landscape is summarized by its mean over directions as a function of the window length $W = (m-1)\tau$. For the Lorenz system with $\rho$ = 28, 35 and 45, the window that maximizes bimodality obeys $W^{∗} \approx c_m P$, where $P$ is the rotation period. The constant $c_m$ rises from 0.48 ($m = 2$) to 0.75 ($m = 6$) and differs by at most 0.014 across $\rho$, while $P$ changes from 0.72 s to 0.54 s. The pattern is absent from phase-randomized surrogates, does not hold as a single constant for the Rössler system, and is absent for the Hénon map. A two-wing oscillator model explains the mechanism. Three analytic propositions show that the projection of a same-wing delay vector is always a single shifted sinusoid, that its bimodality coefficient depends only on $\rho = R^2/C^2$ with a minimum of 7/15 at $\rho = 8/3$, and that even rare wing switches depress the landscape at $W = P$ through outliers that dominate the fourth moment. The mixture formula at $W = P$ matches simulation to within 0.0004 for Lorenz-like switching. A ladder of models shows that the second harmonic of the in-wing waveform and rotation-to-rotation variability close 77 to 84 percent of the gap in $c_m$ between a sinusoidal model and Lorenz for $m$ = 3 to 6, while a growing-spiral hypothesis is ruled out by a shuffled control.
:::

::: {custom-style="KataKunci"}
Keywords: delay embedding; projection pursuit; bimodality coefficient; Lorenz attractor; surrogate data; dynamical systems
:::

# 1 Pendahuluan

Analisis komponen utama (PCA) memproyeksikan data berdimensi tinggi ke beberapa arah dengan variansi terbesar (Pearson, 1901). *Projection pursuit* memperumum ide ini. Setiap arah satuan $u$ diberi skor $I(u)$ yang mengukur seberapa "menarik" proyeksi data ke arah itu, biasanya seberapa jauh distribusinya dari Gaussian, lalu dicari arah dengan skor tertinggi (Friedman & Tukey, 1974; Huber, 1985; Jones & Sibson, 1987). Gaussian dipakai sebagai acuan "tidak menarik" karena sebagian besar proyeksi berdimensi rendah dari data berdimensi tinggi memang mendekati Gaussian (Diaconis & Freedman, 1984). Dalam kedua kerangka itu, yang dicari adalah satu puncak dari fungsi $I(u)$.

Naskah ini berangkat dari pertanyaan sebaliknya. Alih-alih menurunkan dimensi, deret waktu skalar $x(t)$ diangkat ke $\mathbb{R}^m$ melalui *embedding* tunda, lalu seluruh lanskap $I(u)$ di atas bola arah $S^{m-1}$ dipelajari, bukan hanya puncaknya. Pengangkatan linear tidak menambah informasi, karena data tetap berada di subruang berdimensi asal. *Embedding* tunda berbeda karena setiap koordinat membawa nilai deret pada waktu yang berbeda. Teorema Takens (1981) dan perluasannya (Sauer dkk., 1991) menjamin bahwa untuk $m$ yang cukup besar bentuk atraktor dapat direkonstruksi dari satu deret skalar. Pertanyaan yang muncul adalah fitur lanskap mana yang merupakan sifat atraktor, dan fitur mana yang hanya artefak dari pilihan koordinat $\tau$ dan $m$.

Bahan-bahan ide ini sudah dikenal secara terpisah. Gabungan *embedding* tunda dan PCA adalah analisis spektrum singular (Broomhead & King, 1986), dan PCA yang didahului pengangkatan nonlinear adalah PCA kernel (Schölkopf dkk., 1998). Karena setiap proyeksi vektor tunda adalah keluaran sebuah filter FIR, memaksimalkan non-Gaussianitas atas arah setara dengan ICA kanal tunggal (Davies & James, 2007). Fungsi indeks sebagai objek telah dipelajari untuk keperluan optimisasi, misalnya melalui kehalusan dan *squintability* (Laa & Cook, 2020; Zhang dkk., 2026). Fungsi atas semua arah di $S^{m-1}$ juga merupakan inti *Persistent Homology Transform* (Turner dkk., 2014), meski yang diukur di sana adalah topologi, bukan bentuk distribusi. Dari sisi rekonstruksi, Kugiumtzis (1996) berargumen bahwa parameter yang bermakna adalah panjang jendela $\tau_w = (m-1)\tau$, bukan $\tau$ dan $m$ secara terpisah, dan menyarankan jendela minimal satu periode orbit rata-rata untuk estimasi dimensi korelasi.

Sepanjang penelusuran penulis, belum ada studi yang memperlakukan lanskap indeks bentuk distribusi pada atraktor chaos sebagai objek yang bergantung pada parameter *embedding*. Naskah ini memberi tiga kontribusi.

1. Hukum skala empiris $W^{∗} \approx c_m P$ pada keluarga Lorenz: lokasi puncak lanskap rata-rata ditentukan oleh periode putaran $P$, dengan konstanta $c_m$ yang bergantung pada $m$ tetapi hampir tidak bergantung pada $\rho$ (Bagian 3).
2. Bentuk tertutup lanskap untuk model osilator dua sayap, termasuk efek perpindahan sayap (Bagian 4.1 sampai 4.4).
3. Rangkaian model bertingkat yang memisahkan kontribusi bentuk gelombang, variasi antarputaran, dan pertumbuhan amplitudo terhadap kenaikan $c_m$ (Bagian 4.5).

Bagian 2 menjelaskan metode. Bagian 5 membahas posisi hasil terhadap literatur, reproduksibilitas, dan keterbatasan, sedangkan Bagian 6 menutup dengan pertanyaan terbuka. Bukti lengkap ada di Lampiran A dan detail numerik di Lampiran B.

# 2 Kerangka dan metode

## 2.1 *Embedding* tunda dan bingkai filter

Untuk deret skalar $x(t)$, tunda $\tau$ dan dimensi $m$, vektor tunda dan panjang jendelanya adalah

::: {custom-style="Persamaan"}
$$v(t) = \left(x(t),\, x(t-\tau),\, \dots,\, x(t-(m-1)\tau)\right)^\top \in \mathbb{R}^m, \qquad W = (m-1)\tau.$$
:::

Setiap koordinat dinormalisasi ke rerata nol dan simpangan baku satu. Untuk arah satuan $u = (u_0, \dots, u_{m-1})^\top \in S^{m-1}$, proyeksinya adalah

::: {custom-style="Persamaan"}
$$y_u(t) = u^\top v(t) = \sum_{k=0}^{m-1} u_k\, x(t-k\tau).$$
:::

Persamaan (2) menunjukkan bahwa setiap arah adalah filter FIR dengan norma satu dan panjang $W$. Karena $y_{-u} = -y_u$ dan indeks yang dipakai tidak berubah oleh pencerminan, arah cukup dipandang modulo tanda.

## 2.2 Indeks bimodalitas dan lanskap rata-rata

Indeks utama adalah koefisien bimodalitas Sarle dengan koreksi sampel berhingga (Pfister dkk., 2013),

::: {custom-style="Persamaan"}
$$I_3(u) = \frac{\gamma_u^2 + 1}{\kappa_u + \dfrac{3(N-1)^2}{(N-2)(N-3)}},$$
:::

dengan $\gamma_u$ kemencengan dan $\kappa_u$ kurtosis berlebih dari proyeksi $y_u$, serta $N$ banyak titik. Nilai acuannya adalah 1/3 untuk Gaussian, 5/9 untuk distribusi seragam (ambang konvensional bimodalitas), 2/3 untuk distribusi arcsinus (proyeksi sinusoid murni), dan 1 untuk distribusi dua titik. Untuk distribusi simetris dan $N$ besar berlaku $I_3 \to \mu_2^2/\mu_4$, dengan $\mu_k$ momen pusat ke-$k$. Sebagai pemeriksa kedua dipakai aproksimasi negentropi $I_1(u) = \left(\mathbb{E}[\log\cosh y_u] - \mathbb{E}[\log\cosh \nu]\right)^2$ dengan $\nu \sim N(0,1)$ (Hyvärinen, 1998).

Lanskap diringkas sebagai rata-rata atas arah,

::: {custom-style="Persamaan"}
$$\overline{I}_3(W) = \mathbb{E}_{u \sim \operatorname{Unif}(S^{m-1})}\left[I_3(u)\right], \qquad W^{∗} = \arg\max_W \overline{I}_3(W), \qquad c_m = \frac{W^{∗}}{P}.$$
:::

Rata-rata ditaksir dengan 500 sampai 600 arah acak yang tersebar seragam. Lokasi puncak $W^{∗}$ diambil dari maksimum pada grid, lalu dihaluskan dengan interpolasi parabola melalui tiga titik di sekitarnya. Puncak yang jatuh di tepi grid ditandai dan tidak ditafsirkan. Ringkasan kontinu ini dipilih setelah dua ringkasan lain, yaitu fraksi arah dengan $I_3 > 5/9$ dan jumlah "pulau" arah bernilai tinggi, terbukti bergantung pada ambang dan kurang peka.

## 2.3 Sistem, periode, dan pembanding nol

Tiga sistem dipakai: atraktor Lorenz (1963) dengan $\sigma = 10$, $\beta = 8/3$ dan $\rho \in \{28, 35, 45\}$; atraktor Rössler (1976) dengan *a* = *b* = 0,2 dan *c* = 5,7; serta peta Hénon (1976) dengan *a* = 1,4 dan *b* = 0,3. Persamaan diferensial diintegrasikan dengan metode Runge-Kutta orde 4(5) dengan toleransi relatif dan absolut $10^{-8}$ sampai $10^{-9}$ (Lampiran B). Waktu pada kedua sistem kontinu ini tak berdimensi, dan satuannya ditulis s hanya untuk singkatnya. Hanya komponen $x(t)$ yang dipakai untuk lanskap.

Periode $P$ didefinisikan sebagai median selang antara dua maksimum lokal berturut-turut dari $z(t)$ untuk Lorenz dan dari $x(t)$ untuk Rössler. Untuk Hénon, $P$ adalah satu langkah. $P$ sengaja diukur dari keadaan penuh karena penaksir yang hanya memakai $x(t)$ tidak andal: maksimum $x$ melewatkan sekitar seperempat putaran Lorenz, sedangkan maksimum $|x|$ menghitung setiap putaran Rössler dua kali.

Pembanding nol adalah data *surrogate* acak-fase (Theiler dkk., 1992). Spektrum daya dipertahankan dan fase Fourier diacak, sehingga autokorelasi linear tetap sama tetapi struktur nonlinear hilang. Untuk setiap realisasi dibuat satu *surrogate*, dan lanskapnya dihitung dengan cara yang sama.

## 2.4 Rangkaian eksperimen

Eksperimen disusun bertahap. Fase 1 menguji kestabilan lokasi puncak pada Lorenz $\rho = 28$ terhadap kondisi awal dan panjang data. Fase 2 menguji apakah puncak mengikuti periode pada tiga nilai $\rho$, serta pada Rössler dan Hénon. Fase 3 sampai 3c memakai model mainan untuk menjelaskan mekanismenya. Pengaturan lengkap setiap fase ada di Lampiran B.

# 3 Hasil empiris

## 3.1 Lanskap pada satu *embedding*

Gambar 1 memperlihatkan atraktor Lorenz ($\rho = 28$), deret $x(t)$ yang diamati, rekonstruksinya untuk $m = 3$ dan *τ* = 0,12 s, serta lanskap $I_3(u)$ di $S^2$. Lanskap ini jauh dari datar. Nilainya berkisar dari 0,19 sampai 0,75, dan sekitar 25 persen luas bola memiliki $I_3 > 5/9$. Daerah bernilai tinggi muncul berpasangan karena $I(u) = I(-u)$.

Di arah maksimum, proyeksinya bimodal dengan dua lonjakan tajam di tepi (Gambar 1e). Dari arah itu atraktor terlihat dari samping sebagai kurva berbentuk S, dan lonjakan muncul di titik lipatan, tempat banyak titik menumpuk pada satu nilai proyeksi. Di arah minimum, proyeksinya berpuncak tajam di tengah dengan ekor tebal (Gambar 1f).

Indeks yang berbeda sepakat tentang di mana tidak ada struktur, tetapi tidak tentang di mana struktur paling kuat. Korelasi peringkat Spearman antara lanskap $I_1$ dan $I_3$ adalah 0,84, tetapi arah maksimum keduanya berjarak 41 derajat. Karena itu, analisis berikutnya memakai ringkasan rata-rata atas arah dan selalu memeriksa dua indeks.

![Gambar 1. (a) Atraktor Lorenz dengan $\rho = 28$. (b) Deret $x(t)$, satu-satunya data yang dipakai. (c) Rekonstruksi tunda dengan $m = 3$ dan *τ* = 0,12 s. (d) Lanskap $I_3(u)$ di $S^2$ dalam koordinat bola $(\theta, \varphi)$; garis hitam adalah kontur $I_3 = 5/9$, bintang menandai maksimum dan lingkaran menandai minimum. (e, f) Distribusi proyeksi di arah maksimum dan minimum, dibandingkan dengan Gaussian.](gambar/gambar1_lanskap.png){width=16cm}

## 3.2 Skala periode pada keluarga Lorenz

Pada $\rho = 28$, lokasi puncak $\overline{I}_3(W)$ stabil terhadap panjang data dan kondisi awal. Dengan 10 kondisi awal dan dua panjang deret (75 s dan 150 s), $W^{∗}$ berbeda paling banyak 0,008 s antarpanjang deret, pada grid berlangkah 0,08 s. Puncak $\overline{I}_1$ jatuh di $W$ yang sama, dalam satu langkah grid, untuk semua $m$ dari 2 sampai 7. *Surrogate* memberi $\overline{I}_3$ antara 0,31 dan 0,37 (rata-rata 0,345), dekat nilai Gaussian 1/3, tanpa puncak.

Untuk menguji apakah puncak mengikuti periode, eksperimen diulang pada $\rho$ = 28, 35 dan 45, dengan lima realisasi per sistem dan sumbu $W$ yang dinormalisasi oleh periode masing-masing. Periodenya adalah *P* = 0,724 ± 0,016 s, 0,627 ± 0,005 s dan 0,540 ± 0,007 s. Gambar 2 menunjukkan bahwa kurva untuk ketiga $\rho$ hampir berimpit setelah normalisasi, dan Tabel 1 merangkum lokasi puncaknya.

Untuk $m = 2$ dan $m = 3$, nilai $c_m$ berbeda paling banyak 0,002 antarnilai $\rho$. Untuk $m = 4$ sampai 6, perbedaannya 0,008 sampai 0,014 dan konsisten turun ketika $\rho$ naik, sedikit melebihi simpangan baku antarrealisasi (0,002 sampai 0,006). Jadi hukum skala

$$W^{∗} \approx c_m P$$

berlaku hampir tepat, bukan tepat. Konstanta $c_m$ sendiri naik dengan $m$, dari 0,48 pada $m = 2$ ke 0,75 pada $m = 6$.

: Tabel 1. Nilai $c_m = W^{∗}/P$ (rata-rata ± simpangan baku dari lima realisasi) pada keluarga Lorenz. Selisih yang dikutip di teks dihitung sebelum pembulatan.

| $\rho$ | *P* (s) | $m = 2$ | $m = 3$ | $m = 4$ | $m = 5$ | $m = 6$ |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 28 | 0,724 ± 0,016 | 0,478 ± 0,004 | 0,626 ± 0,005 | 0,660 ± 0,003 | 0,686 ± 0,004 | 0,748 ± 0,006 |
| 35 | 0,627 ± 0,005 | 0,479 ± 0,003 | 0,627 ± 0,004 | 0,655 ± 0,003 | 0,681 ± 0,003 | 0,744 ± 0,003 |
| 45 | 0,540 ± 0,007 | 0,479 ± 0,002 | 0,625 ± 0,003 | 0,647 ± 0,002 | 0,676 ± 0,003 | 0,741 ± 0,004 |

Untuk $m = 2$, kurva juga memiliki lembah tajam di *W*/*P* = 1,03 ± 0,02, 1,07 ± 0,01 dan 1,05 ± 0,02 untuk ketiga $\rho$. Pada $\rho = 28$ nilainya turun sampai 0,344, praktis sama dengan *surrogate*. Lembah ini sudah tidak tampak pada $m = 4$ (Gambar 2b).

Lokasi puncak tidak ditemukan oleh aturan baku pemilihan $\tau$. Minimum pertama informasi mutual (Fraser & Swinney, 1986) jatuh di *τ* = 0,164 ± 0,005 s, yang untuk $m = 3$ memberi *W* = 0,33 s, sedangkan puncak lanskap ada di $W^{∗} = c_3 P$ = 0,45 s. Nol pertama autokorelasi (2,4 ± 0,8 s) bahkan tidak stabil antarrealisasi.

![Gambar 2. Rata-rata $I_3$ atas arah sebagai fungsi $W/P$ untuk (a) $m = 2$ dan (b) $m = 4$ pada tiga nilai $\rho$. Pita menunjukkan ±1 simpangan baku antarrealisasi, garis putus-putus abu-abu adalah *surrogate* ($\rho = 28$), dan garis horizontal adalah nilai Gaussian 1/3. (c) $c_m = W^{∗}/P$ terhadap $m$, dengan batang galat ±1 simpangan baku dari lima realisasi.](gambar/gambar2_skala.png){width=16cm}

## 3.3 Batas cakupan: Rössler dan Hénon

Pada sistem Rössler (*P* = 5,96 s), kurva $\overline{I}_3(W)$ bergelombang dengan beberapa puncak yang hampir sama tinggi. Akibatnya lokasi maksimum global melompat antarnilai $m$: $W^{∗}/P$ bernilai 0,80; 1,67; 1,19 dan 1,63 untuk $m = 2$ sampai 5, dan jatuh di tepi grid untuk $m = 6$. Spektrum Rössler didominasi satu puncak sempit, sehingga *surrogate*-nya mempertahankan osilasi yang hampir sinusoidal dan memberi $\overline{I}_3$ = 0,45, jauh di atas 1/3. Untuk sistem seperti ini, perbandingan harus dilakukan terhadap *surrogate*-nya sendiri.

Pada peta Hénon, maksimum selalu berada di jendela terkecil, yaitu satu langkah per koordinat, tanpa puncak di dalam grid. Hal ini sesuai harapan untuk sistem tanpa putaran. Hukum skala di atas karena itu hanya diklaim untuk keluarga atraktor dengan geometri dua sayap seperti Lorenz.

# 4 Mekanisme: model osilator dua sayap

## 4.1 Model

Dua ciri Lorenz yang tampak dari $x(t)$ adalah osilasi di dalam sayap dan perpindahan tak teratur antarsayap. Model berikut mempertahankan keduanya tanpa chaos:

::: {custom-style="Persamaan"}
$$x(t) = s(t)\left[c + a\cos\theta(t)\right], \qquad \theta(t) = \frac{2\pi t}{P}, \qquad s(t) \in \{-1, +1\}.$$
:::

Tanda $s(t)$ hanya berganti di fase $\theta = \pi$, yaitu saat $|x|$ minimum, dan lama tinggal di satu sayap berdistribusi geometrik dengan rata-rata $D$ putaran. Syarat $a/c < 1$ menjamin bahwa $x$ tidak memotong nol di dalam sayap. Pada Lorenz $\rho = 28$, lintasan berpindah sayap kira-kira sekali setiap 2,2 putaran, sehingga *D* = 2,2 dipakai sebagai nilai "ala Lorenz". Dalam model ini tunda antarkoordinat menjadi beda fase $\psi = 2\pi\tau/P$. Pasangan atau vektor yang semua koordinatnya berada di sayap yang sama disebut sesayap.

## 4.2 Proyeksi pasangan sesayap

::: {custom-style="JudulProposisi"}
Proposisi 1 (bentuk proyeksi)
:::

::: {custom-style="Proposisi"}
Jika semua koordinat $v(t)$ berada di sayap yang sama, maka untuk setiap arah $u \in S^{m-1}$ berlaku
:::

::: {custom-style="Persamaan"}
$$y_u = s\left[C + R\cos(\theta - \delta)\right], \qquad C = c\sum_{k=0}^{m-1} u_k, \qquad R = a\,|Z|, \qquad \delta = -\arg Z, \qquad Z = \sum_{k=0}^{m-1} u_k\, e^{-ik\psi}.$$
:::

Buktinya (Lampiran A.2) hanya memakai $\cos(\theta - k\psi) = \operatorname{Re}\, e^{i(\theta - k\psi)}$: jumlah sinusoid berfrekuensi sama selalu berupa satu sinusoid. Proposisi ini berlaku untuk sembarang $m$, dan seluruh ketergantungan pada arah $u$ dan tunda $\tau$ terkumpul dalam dua bilangan, $C$ dan $R$.

## 4.3 Koefisien bimodalitas sebagai fungsi satu variabel

::: {custom-style="JudulProposisi"}
Proposisi 2 (fungsi bimodalitas)
:::

::: {custom-style="Proposisi"}
Misalkan $\theta$ seragam pada [0, 2π), $s = \pm 1$ dengan peluang sama dan bebas dari $\theta$, serta $C \neq 0$. Dengan $\rho = R^2/C^2$, proyeksi pada Proposisi 1 memenuhi $\mathbb{E}[y] = 0$, $\mu_2 = C^2 + R^2/2$ dan $\mu_4 = C^4 + 3C^2R^2 + 3R^4/8$, sehingga
:::

::: {custom-style="Persamaan"}
$$b(\rho) = \frac{\mu_2^2}{\mu_4} = \frac{2(2+\rho)^2}{3\rho^2 + 24\rho + 8}.$$
:::

::: {custom-style="Proposisi"}
Fungsi $b$ turun dari 1 (untuk $\rho \to 0$) ke minimum 7/15 di $\rho = 8/3$, lalu naik menuju 2/3 (untuk $\rho \to \infty$), karena
:::

::: {custom-style="Persamaan"}
$$b'(\rho) = \frac{8(\rho+2)(3\rho-8)}{(3\rho^2 + 24\rho + 8)^2}.$$
:::

Kedua ujungnya bermakna. Limit $\rho \to 0$ adalah distribusi dua titik $\pm C$ yang bimodal sempurna, sedangkan $\rho \to \infty$ adalah sinusoid murni dengan distribusi arcsinus. Gambar 3a menggambar $b(\rho)$ bersama nilai-nilai acuannya.

::: {custom-style="JudulProposisi"}
Akibat 1 (lanskap *m* = 2 tanpa perpindahan sayap)
:::

::: {custom-style="Proposisi"}
Untuk $m = 2$ dan $u = (\cos\alpha, \sin\alpha)$, berlaku $C^2 = c^2(1 + \sin 2\alpha)$ dan $R^2 = a^2(1 + \sin 2\alpha \cos\psi)$ dengan $\psi = 2\pi W/P$, sehingga
:::

::: {custom-style="Persamaan"}
$$\rho(\alpha, W) = \left(\frac{a}{c}\right)^2 \frac{1 + \sin 2\alpha\, \cos(2\pi W/P)}{1 + \sin 2\alpha}.$$
:::

::: {custom-style="Proposisi"}
Akibatnya lanskap rata-rata $L(W) = \mathbb{E}_\alpha\left[b(\rho(\alpha, W))\right]$ bergantung pada $W$ hanya melalui $\cos(2\pi W/P)$. Lanskap ini periodik dengan periode $P$, simetris terhadap $W = P/2$, dan stasioner di $W = P/2$.
:::

Di $W = P$, dan juga untuk $W \to 0$, berlaku $\rho = (a/c)^2$ untuk semua arah sehingga $L = b\left((a/c)^2\right)$. Di $W = P/2$, nilai $\rho$ menyapu seluruh [0, ∞) ketika $\alpha$ berubah. Secara numerik, untuk 0 < *a*/*c* < 1, titik $W = P/2$ adalah maksimum global lanskap pada $(0, P)$ jika dan hanya jika *a*/*c* > 0,4573. Di bawah ambang itu, maksimum berpindah ke ujung $W \to 0$ dan $W \to P$ (Gambar 3b).

Penjelasan geometrisnya sederhana. Di $W = P$, dan juga untuk $W \to 0$, kedua koordinat pasangan sesayap bernilai sama, sehingga titik-titiknya terletak pada diagonal dan setiap proyeksi hanyalah salinan berskala dari $x$ sendiri. Di $W = P/2$, pasangan sesayap berbentuk $(c + a\cos\theta,\; c - a\cos\theta)$, yaitu titik-titik pada segmen garis yang tegak lurus diagonal, satu segmen untuk setiap sayap. Arah di dekat diagonal meratakan setiap segmen menjadi hampir satu titik dan melihat dua gugus rapat, arah di dekat antidiagonal melihat satu sebaran arcsinus, dan arah di antaranya melewati minimum 7/15. Jika $a/c$ kecil, $x$ sendiri sudah hampir berupa dua titik, sehingga salinan berskala di $W = P$ sulit dikalahkan. Jika $a/c$ besar, setiap sayap $x$ melebar dan $b\left((a/c)^2\right)$ turun menuju 18/35 ≈ 0,51 saat *a*/*c* → 1, sehingga pemutaran segmen di $W = P/2$ memberi banyak arah pandangan yang lebih bimodal.

Simulasi mendukung ketiga hasil ini. Model dengan perpindahan langka (*D* = 40, *a*/*c* = 0,6) cocok dengan rumus untuk $W \le P/2$ dengan selisih maksimum 0,008 (Gambar 3c). Sinusoid murni memberi 0,6657 untuk semua $W$, sama dengan nilai harapan 2/3 setelah koreksi sampel berhingga Sarle untuk sekitar 4000 titik. Ambang $a/c$ juga konsisten dengan simulasi yang memuat perpindahan: pada model dengan *D* = 2,2, nilai *a*/*c* = 0,30 dan 0,45 tidak menghasilkan puncak di dalam grid untuk semua $m$, sedangkan *a*/*c* ≥ 0,60 menghasilkan puncak di dekat $P/2$ untuk $m = 2$, dengan *c*~2~ antara 0,45 dan 0,49.

## 4.4 Perpindahan sayap dan lembah di *W* = *P*

Rumus tanpa perpindahan tidak memiliki lembah di $W = P$ (garis hitam pada Gambar 3c), padahal simulasi dan Lorenz memilikinya. Lembah itu berasal dari pasangan yang melintasi perpindahan sayap.

::: {custom-style="JudulProposisi"}
Proposisi 3 (campuran di *W* = *P*)
:::

::: {custom-style="Proposisi"}
Untuk $m = 2$, pasangan beda sayap memiliki bentuk Proposisi 1 dengan $u_1$ diganti $-u_1$, yaitu $C' = c(u_0 - u_1)$ dan $R'^2 = a^2(1 - \sin 2\alpha\cos\psi)$. Di $W = P$, setiap pasangan melintasi tepat satu batas putaran, sehingga proporsi pasangan sesayap adalah $p = 1 - 1/D$ dan tidak bergantung pada $\theta$. Dengan $q = \sin 2\alpha$, momen kedua jenis pasangan bercampur secara linear dan
:::

::: {custom-style="Persamaan"}
$$b_{\text{mix}}(\alpha) = b\left((a/c)^2\right)\, \frac{\left[1 + (2p-1)q\right]^2}{1 + 2(2p-1)q + q^2}.$$
:::

::: {custom-style="JudulProposisi"}
Akibat 2 (arah tanpa sinyal sesayap)
:::

::: {custom-style="Proposisi"}
Pada dua arah diagonal berlaku
:::

::: {custom-style="Persamaan"}
$$b_{\text{mix}}(3\pi/4) = (1-p)\, b\left((a/c)^2\right), \qquad b_{\text{mix}}(\pi/4) = p\, b\left((a/c)^2\right).$$
:::

Akibat 2 menjelaskan mengapa perpindahan yang langka pun berpengaruh besar. Di arah $\alpha = 3\pi/4$, semua pasangan sesayap terproyeksi tepat ke nol karena $C = R = 0$. Satu-satunya variansi berasal dari pasangan beda sayap yang langka, yang menjadi pencilan dan menguasai momen keempat. Makin langka perpindahan, makin kecil $b$ di arah itu. Rata-ratanya atas $\alpha$ pun pulih ke nilai tanpa perpindahan dengan sangat lambat. Untuk *a*/*c* = 0,6, nilai tanpa perpindahan adalah 0,654, sedangkan rata-rata campuran adalah 0,483, 0,512 dan 0,578 untuk *p* = 0,9; 0,975 dan 0,999. Satu perpindahan per seribu putaran masih menurunkan lanskap rata-rata lebih dari 10 persen.

Dengan proporsi perpindahan yang benar-benar terjadi di simulasi (0,4268 untuk *D* = 2,2 dan 0,0235 untuk *D* = 40), Persamaan (10) memberi 0,4629 dan 0,5131. Perhitungan simulasi dengan 4000 arah seragam dan semua titik memberi 0,4633 dan 0,5187. Sisa selisih untuk *D* = 40 wajar, karena simulasi itu hanya memuat sekitar 47 perpindahan dalam 2000 putaran. Lembah Lorenz di dekat $W = P$ (Bagian 3.2) dengan demikian adalah jejak perpindahan sayap pada lanskap.

![Gambar 3. (a) Fungsi $b(\rho)$ dari Proposisi 2, dengan nilai acuan dua titik (1), sinusoid (2/3), seragam (5/9) dan Gaussian (1/3). (b) Lanskap $m = 2$ tanpa perpindahan sayap: rata-rata $b$ di $W = P/2$ dan di $W \to 0$ atau $W = P$ sebagai fungsi $a/c$; $P/2$ menjadi maksimum global jika *a*/*c* > 0,457. (c) Rumus tanpa perpindahan (hitam) dibandingkan dengan simulasi model dua sayap (*a*/*c* = 0,6) untuk *D* = 40 dan *D* = 2,2. Belah ketupat adalah prediksi Proposisi 3 di $W = P$ dengan proporsi perpindahan yang terealisasi.](gambar/gambar3_analitik.png){width=16cm}

## 4.5 Kenaikan *c*~*m*~ terhadap *m*

Untuk $m \ge 3$, Proposisi 1 dan 2 tetap berlaku, dan argumen simetri pada Akibat 1 juga tetap berjalan. Karena $u$ real, $|Z|$ tidak berubah oleh $\psi \to -\psi$, sehingga lanskap tanpa perpindahan simetris terhadap $\tau = P/2$ dan stasioner di sana, yaitu di $W = (m-1)P/2$ (Lampiran A.4). Secara numerik, titik itu adalah maksimum global untuk *a*/*c* antara 0,6 dan 1,0 pada $m = 3$ dan $m = 4$, sehingga rumus tanpa perpindahan meramalkan *c*~3~ = 1 dan *c*~4~ = 1,5. Model sinusoid dengan perpindahan ala Lorenz justru memuncak di 0,43*P* sampai 0,54*P* untuk $m = 3$ dan *a*/*c* antara 0,6 dan 1,0 (Tabel 2 dan model A pada Tabel 3), sedangkan Lorenz memberi *c*~3~ = 0,626. Jadi untuk $m \ge 3$ perpindahan sayap menentukan lokasi puncak, bukan sekadar mengoreksinya, dan $c_m$ belum dapat diturunkan secara tertutup. Karena itu kenaikan $c_m$ dipelajari dengan rangkaian model bertingkat, yang masing-masing menambahkan satu ciri Lorenz yang diukur langsung.

Tingkat amplitudo memengaruhi kemiringan. Pada model sinusoid dengan amplitudo tetap, kemiringan $c_6 - c_3$ naik dari 0,035 ke 0,069 ketika *a*/*c* naik dari 0,6 ke 0,9, tetapi $c_3$ tetap di bawah 0,53, jauh dari 0,626 pada Lorenz (Tabel 2).

Hipotesis pertumbuhan amplitudo spiral tidak terdukung. Di Lorenz, amplitudo osilasi tumbuh di dalam satu sayap sebelum lintasan terlempar ke sayap lain. Model spiral meniru hal ini dengan $A_{j+1} = \lambda A_j$, dimulai dari amplitudo acak antara 0,3 dan 0,5, sampai mencapai amplitudo kritis 0,9 yang memicu perpindahan. Sebagai kontrol, urutan amplitudo di dalam setiap sayap diacak, sementara barisan sayap dan himpunan amplitudonya dipertahankan. Dibanding kontrol, spiral menurunkan $c_3$ sebesar 0,009 sampai 0,015, yaitu menjauh dari nilai Lorenz, dan tidak mengubah $c_6$ lebih dari 0,004. Kemiringan spiral yang sedikit lebih besar berasal dari $c_3$ yang lebih rendah, bukan dari $c_6$ yang lebih tinggi.

: Tabel 2. Lokasi puncak pada model amplitudo tetap dan model spiral beserta kontrol urutan acaknya (rata-rata tiga realisasi). Kolom *D* adalah rata-rata lama tinggal terukur, dalam putaran. Kolom $c_6 - c_3$ dan selisih yang dikutip di teks dihitung sebelum pembulatan.

| Model | Parameter | *D* | $c_3$ | $c_6$ | $c_6 - c_3$ |
|:--|:--:|:--:|:--:|:--:|:--:|
| Amplitudo tetap | *a*/*c* = 0,60 | 2,21 | 0,426 | 0,461 | +0,035 |
| Amplitudo tetap | *a*/*c* = 0,75 | 2,19 | 0,488 | 0,539 | +0,051 |
| Amplitudo tetap | *a*/*c* = 0,90 | 2,29 | 0,528 | 0,598 | +0,069 |
| Spiral | *λ* = 1,25 | 5,17 | 0,475 | 0,537 | +0,063 |
| Urutan acak | *λ* = 1,25 | 5,17 | 0,489 | 0,536 | +0,047 |
| Spiral | *λ* = 1,6 | 3,26 | 0,479 | 0,538 | +0,059 |
| Urutan acak | *λ* = 1,6 | 3,26 | 0,490 | 0,534 | +0,044 |
| Spiral | *λ* = 2,0 | 2,76 | 0,484 | 0,541 | +0,056 |
| Urutan acak | *λ* = 2,0 | 2,76 | 0,494 | 0,539 | +0,045 |
| Lorenz | *ρ* = 28 | sekitar 2,2 | 0,626 | 0,748 | +0,122 |

Bentuk gelombang di dalam sayap diukur langsung dari lintasan Lorenz. Setiap putaran yang tidak berpindah sayap dipotong dari satu minimum $|x|$ ke minimum berikutnya, lalu dinormalisasi menjadi $g = (|x| - c)/a$, dengan $c$ rerata dan $a$ setengah jangkauan $|x|$ pada putaran itu (Gambar 4a). Dari lima lintasan dengan kondisi awal berbeda, yang masing-masing menyumbang sekitar 813 putaran sesayap, harmonik kedua bentuk rata-ratanya berukuran 0,171 ± 0,003 kali harmonik dasar, harmonik ketiga 0,021 ± 0,001, dan median $a/c$ adalah 1,00 ± 0,01, dengan rentang antarkuartil sekitar 0,75 sampai 1,21. Asimetrinya vertikal, bukan temporal: puncak $|x|$ lebih tinggi dan tajam, lembahnya lebar dan dangkal. Asimetri temporal memang tidak dapat memengaruhi lanskap rata-rata, karena membalik waktu hanya menukar urutan koordinat *embedding*, sedangkan rata-rata atas arah seragam tidak berubah oleh pertukaran itu.

Tangga model (Tabel 3 dan Gambar 4b) memakai perpindahan yang sama (*D* = 2,2) dan $a/c$ median Lorenz, dan setiap anak tangga menambahkan satu ciri. Model A memakai cosinus, model B memakai bentuk Lorenz yang dipotong ke harmonik 1 dan 2, model C memakai bentuk rata-rata penuh, model D menambahkan $a/c$ acak per putaran, dan model E memakai bentuk dan $a/c$ dari putaran Lorenz nyata yang diundi.

Harmonik kedua adalah pembawa utama. Model B menutup 61 persen selisih $c_3$ antara model A dan Lorenz, serta 43 sampai 52 persen untuk $m = 4$ sampai 6. Model C, yang memuat semua harmonik, tidak lebih baik dari B. Variasi antarputaran menambah efek di $m$ besar: model E menutup 77 sampai 84 persen selisih untuk semua $m$ dari 3 sampai 6.

Model D tidak dapat dipakai. Di Lorenz, putaran dengan $a/c$ besar memiliki lembah ternormalisasi yang lebih dangkal (korelasi 0,62 ± 0,02 antara $a/c$ dan nilai minimum $g$). Memasangkan bentuk rata-rata dengan $a/c$ acak membuat sekitar 41 persen putaran memotong nol, sehingga nilai-nilai itu harus dibuang dan model menjadi bias. Bentuk dan amplitudo terikat satu sama lain, dan hanya model E yang mempertahankan ikatan itu.

Sisa selisih antara model E dan Lorenz adalah 0,018 sampai 0,025. Sisa ini memperkirakan gabungan pengaruh dua ciri yang tidak ada di model mana pun: lintasan perpindahan sayap yang sebenarnya (di model berupa lompatan seketika) dan variasi periode antarputaran. Dua run independen dengan lintasan Lorenz berbeda memberi $c_6$ = 0,728 dan 0,714 untuk model E, sehingga sisa itu hanya diketahui secara kasar.

: Tabel 3. Tangga model: $c_m = W^{∗}/P$ (rata-rata tiga realisasi; simpangan baku antarrealisasi paling besar 0,008). Model D tidak sah karena mengabaikan ikatan bentuk dan amplitudo (lihat teks). Kolom $c_6 - c_3$ serta persentase dan sisa yang dikutip di teks dihitung sebelum pembulatan.

| Model | $c_2$ | $c_3$ | $c_4$ | $c_5$ | $c_6$ | $c_6 - c_3$ |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| A, cosinus | 0,495 | 0,543 | 0,549 | 0,585 | 0,619 | +0,076 |
| B, harmonik 1 dan 2 | 0,491 | 0,594 | 0,606 | 0,629 | 0,676 | +0,083 |
| C, bentuk rata-rata | 0,488 | 0,587 | 0,601 | 0,626 | 0,673 | +0,086 |
| D, ditambah *a*/*c* acak | 0,468 | 0,566 | 0,567 | 0,603 | 0,638 | +0,071 |
| E, putaran nyata | 0,461 | 0,608 | 0,635 | 0,669 | 0,728 | +0,120 |
| Lorenz *ρ* = 28 | 0,478 | 0,626 | 0,660 | 0,686 | 0,748 | +0,122 |

![Gambar 4. (a) Bentuk gelombang di dalam sayap: 40 putaran Lorenz nyata (abu-abu), bentuk rata-rata (hitam, model C), potongannya ke harmonik 1 dan 2 (oranye, model B), dan cosinus (biru, model A). Loncatan kecil di fase 0,5 adalah batas antarputaran, tempat dua minimum berurutan tidak sama dalam. (b) $c_m$ terhadap $m$ untuk model A, B dan E (rata-rata ± simpangan baku dari tiga realisasi) dan untuk Lorenz $\rho = 28$ (lima realisasi). Model C hampir berimpit dengan B, dan model D tidak sah; keduanya tercantum di Tabel 3.](gambar/gambar4_tangga.png){width=16cm}

# 5 Pembahasan

## 5.1 Posisi terhadap literatur

Hasil ini sejalan dengan argumen Kugiumtzis (1996) bahwa panjang jendela, bukan $\tau$ dan $m$ secara terpisah, adalah parameter yang bermakna, dan bahwa jendela berkaitan dengan periode orbit rata-rata. Namun arahnya berbeda. Kugiumtzis menyarankan $\tau_w \ge \tau_p$ sebagai batas bawah untuk estimasi dimensi korelasi, sedangkan lanskap bimodalitas di sini memuncak di $W^{∗} = c_m P$ dengan $c_m$ antara 0,48 dan 0,75, yaitu sebelum satu periode penuh. Jendela yang memaksimalkan struktur distribusi proyeksi lebih pendek daripada jendela yang disarankan untuk mengukur dimensi, karena kedua kriteria menjawab pertanyaan yang berbeda.

Bingkai filter pada Persamaan (2) menempatkan hasil ini dekat dengan ICA kanal tunggal (Davies & James, 2007), yang mencari satu filter $u$ dengan non-Gaussianitas maksimum. Objek di sini berbeda, yaitu rata-rata atas semua filter sebagai fungsi panjang jendela, bukan filter optimalnya. Demikian pula, studi lanskap indeks dalam *projection pursuit* (Laa & Cook, 2020; Zhang dkk., 2026) memakai kehalusan dan *squintability* untuk menilai kemudahan optimisasi, bukan ketergantungan lanskap pada parameter *embedding* suatu sistem dinamik. Peringatan Bickel dkk. (2018), bahwa *projection pursuit* dapat menemukan struktur palsu pada data Gaussian berdimensi tinggi, juga relevan. Rata-rata atas arah dan pembanding *surrogate* dipakai justru untuk menghindari jebakan itu.

## 5.2 Reproduksibilitas

Hasil yang hanya melibatkan model mainan dapat direproduksi hampir persis, karena semua bilangan acak berasal dari benih tetap. Run Fase 3 dan 3b di dua mesin dengan sistem operasi berbeda (Windows dan Linux) memberi nilai indeks yang berbeda paling banyak $3 \times 10^{-16}$, setingkat galat pembulatan mesin. Hasil yang mengintegrasikan persamaan Lorenz hanya dapat direproduksi secara statistik. Perbedaan pembulatan mesin sudah cukup untuk memisahkan lintasan: menggeser kondisi awal sebesar $10^{-14}$ mengubah jumlah putaran sesayap dari 824 menjadi 830, dan run di mesin lain memberi 817. Karena itu semua besaran yang melibatkan Lorenz dilaporkan sebagai rata-rata dan simpangan baku dari beberapa lintasan, dan ketidakpastian model E diperkirakan dari dua run dengan lintasan berbeda.

## 5.3 Keterbatasan

Beberapa keterbatasan perlu dicatat.

- Hukum skala hanya diuji pada satu keluarga atraktor, dan hasil Rössler menunjukkan bahwa hukum itu tidak universal.
- Hanya satu keluarga indeks yang dipakai, yaitu bentuk distribusi proyeksi satu dimensi, dengan negentropi sebagai pemeriksa.
- Resolusi grid $W$ membatasi presisi $W^{∗}$ sebelum interpolasi. Grid Fase 1 kasar (0,08 s), sehingga simpangan nol antarrealisasi di sana sebagian adalah artefak grid.
- Proposisi 3 eksak hanya di $W = P$. Di jendela lain, peluang sebuah pasangan melintasi perpindahan bergantung pada fase $\theta$.
- Nilai $c_m$ untuk $m \ge 3$ belum diturunkan secara analitik, dan pergeseran kecil $c_m$ terhadap $\rho$ untuk $m \ge 4$ belum dijelaskan.
- Model memperlakukan perpindahan sayap sebagai lompatan seketika, padahal di Lorenz perpindahan adalah satu ekskursi panjang yang melewati dekat nol.

# 6 Kesimpulan

Lanskap indeks bimodalitas pada atraktor Lorenz hasil *embedding* tunda memiliki struktur yang diatur oleh dinamika, bukan oleh pilihan koordinat. Lokasi puncak rata-ratanya mengikuti $W^{∗} \approx c_m P$ pada tiga nilai $\rho$. Model osilator dua sayap menjelaskan tiga pola utamanya: puncak di $P/2$ untuk $m = 2$ berasal dari geometri di dalam sayap dan hanya menjadi maksimum global jika *a*/*c* > 0,457; lembah di $W = P$ berasal dari perpindahan sayap, dan rumus campurannya cocok dengan simulasi; serta kenaikan $c_m$ terhadap $m$ terutama berasal dari harmonik kedua bentuk gelombang dan variasi antarputaran.

Empat pertanyaan tetap terbuka.

1. Rumus tertutup lanskap dengan perpindahan sayap untuk $W$ sembarang dan $m \ge 3$, yang akan menurunkan $c_m$ secara analitik.
2. Perubahan kualitatif lanskap ketika $m$ melewati ambang Takens. Argumen proyeksi berdimensi rendah (Diaconis & Freedman, 1984) mengisyaratkan bahwa perubahannya bertahap.
3. Peta ketidaksepakatan antarindeks sebagai objek tersendiri.
4. Keluarga atraktor lain dan data nyata.

# Pernyataan kontribusi dan penggunaan AI

Penulis merumuskan ide dan pertanyaan riset, menurunkan Proposisi 1 sampai 3, mengajukan hipotesis harmonik dan spiral, serta menjalankan dan mereproduksi semua simulasi. Kode simulasi, pemeriksaan numerik dan simbolik, perluasan argumen simetri untuk sembarang $m$ di Lampiran A.4, penelusuran literatur, gambar, dan draf naskah disusun dengan bantuan asisten AI (Claude, Anthropic). Penulis bertanggung jawab penuh atas isi naskah.

# Ketersediaan kode dan data

Skrip Python untuk semua eksperimen (invarians_lanskap.py, fase1_W.py, fase2_sistem.py, fase3_model_mainan.py, fase3b_spiral.py dan fase3c_gelombang.py) beserta skrip gambar dan keluaran CSV tersedia dari penulis dan akan diunggah ke repositori publik sebelum naskah dikirim.

# Daftar pustaka

::: {custom-style="Pustaka"}
Bickel, P. J., Kur, G., & Nadler, B. (2018). Projection pursuit in high dimensions. *Proceedings of the National Academy of Sciences*, *115*(37), 9151-9156. https://doi.org/10.1073/pnas.1801177115

Broomhead, D. S., & King, G. P. (1986). Extracting qualitative dynamics from experimental data. *Physica D*, *20*(2-3), 217-236.

Davies, M. E., & James, C. J. (2007). Source separation using single channel ICA. *Signal Processing*, *87*(8), 1819-1832. https://doi.org/10.1016/j.sigpro.2007.01.011

Diaconis, P., & Freedman, D. (1984). Asymptotics of graphical projection pursuit. *The Annals of Statistics*, *12*(3), 793-815.

Fraser, A. M., & Swinney, H. L. (1986). Independent coordinates for strange attractors from mutual information. *Physical Review A*, *33*(2), 1134-1140.

Friedman, J. H., & Tukey, J. W. (1974). A projection pursuit algorithm for exploratory data analysis. *IEEE Transactions on Computers*, *C-23*(9), 881-890.

Hénon, M. (1976). A two-dimensional mapping with a strange attractor. *Communications in Mathematical Physics*, *50*(1), 69-77.

Huber, P. J. (1985). Projection pursuit. *The Annals of Statistics*, *13*(2), 435-475.

Hyvärinen, A. (1998). New approximations of differential entropy for independent component analysis and projection pursuit. Dalam *Advances in Neural Information Processing Systems 10* (hlm. 273-279). MIT Press.

Jones, M. C., & Sibson, R. (1987). What is projection pursuit? *Journal of the Royal Statistical Society: Series A*, *150*(1), 1-36.

Kugiumtzis, D. (1996). State space reconstruction parameters in the analysis of chaotic time series: The role of the time window length. *Physica D*, *95*(1), 13-28. https://doi.org/10.1016/0167-2789(96)00054-1

Laa, U., & Cook, D. (2020). Using tours to visually investigate properties of new projection pursuit indexes with application to problems in physics. *Computational Statistics*, *35*, 1171-1205. https://doi.org/10.1007/s00180-020-00954-8

Lorenz, E. N. (1963). Deterministic nonperiodic flow. *Journal of the Atmospheric Sciences*, *20*(2), 130-141.

Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. *Philosophical Magazine*, *2*(11), 559-572.

Pfister, R., Schwarz, K. A., Janczyk, M., Dale, R., & Freeman, J. B. (2013). Good things peak in pairs: A note on the bimodality coefficient. *Frontiers in Psychology*, *4*, 700. https://doi.org/10.3389/fpsyg.2013.00700

Rössler, O. E. (1976). An equation for continuous chaos. *Physics Letters A*, *57*(5), 397-398.

Sauer, T., Yorke, J. A., & Casdagli, M. (1991). Embedology. *Journal of Statistical Physics*, *65*(3-4), 579-616.

Schölkopf, B., Smola, A., & Müller, K.-R. (1998). Nonlinear component analysis as a kernel eigenvalue problem. *Neural Computation*, *10*(5), 1299-1319.

Takens, F. (1981). Detecting strange attractors in turbulence. Dalam D. Rand & L.-S. Young (Ed.), *Dynamical Systems and Turbulence, Warwick 1980* (Lecture Notes in Mathematics, Vol. 898, hlm. 366-381). Springer.

Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D. (1992). Testing for nonlinearity in time series: The method of surrogate data. *Physica D*, *58*(1-4), 77-94.

Turner, K., Mukherjee, S., & Boyer, D. M. (2014). Persistent homology transform for modeling shapes and surfaces. *Information and Inference*, *3*(4), 310-344. https://doi.org/10.1093/imaiai/iau011

Zhang, H. S., Cook, D., Langrené, N., & Leung, J. W. Y. (2026). Squintability and other metrics for assessing projection pursuit indexes, and guiding optimization choices. *Journal of Computational and Graphical Statistics*. https://doi.org/10.1080/10618600.2026.2614076
:::

# Lampiran A. Bukti

Lampiran ini merapikan turunan penulis untuk Proposisi 1 sampai 3. Indeks koordinat ditulis mulai dari nol, sehingga $u = (u_0, \dots, u_{m-1})$ dan koordinat ke-$k$ adalah $x(t - k\tau)$.

## A.1 Momen fase seragam

Untuk $\theta$ seragam pada [0, 2π), identitas $\cos^2\theta = (1 + \cos 2\theta)/2$, $\cos^3\theta = (3\cos\theta + \cos 3\theta)/4$ dan $\cos^4\theta = (3 + 4\cos 2\theta + \cos 4\theta)/8$ memberi

::: {custom-style="PersamaanA"}
$$\mathbb{E}[\cos\theta] = \mathbb{E}[\cos^3\theta] = 0, \qquad \mathbb{E}[\cos^2\theta] = \frac{1}{2}, \qquad \mathbb{E}[\cos^4\theta] = \frac{3}{8}.$$
:::

Karena $\theta - \delta$ (modulo $2\pi$) juga seragam, momen yang sama berlaku untuk $X = \cos(\theta - \delta)$ berapa pun nilai $\delta$.

## A.2 Bukti Proposisi 1

Jika semua koordinat berada di sayap yang sama, koordinat ke-$k$ adalah $s\,[c + a\cos(\theta - k\psi)]$. Maka

::: {custom-style="PersamaanA"}
$$y_u = \sum_{k=0}^{m-1} u_k\, s\left[c + a\cos(\theta - k\psi)\right] = s\left[c\sum_{k=0}^{m-1} u_k + a\sum_{k=0}^{m-1} u_k \cos(\theta - k\psi)\right].$$
:::

Dengan $\cos x = \operatorname{Re}\, e^{ix}$,

::: {custom-style="PersamaanA"}
$$\sum_{k=0}^{m-1} u_k \cos(\theta - k\psi) = \operatorname{Re}\left[e^{i\theta} \sum_{k=0}^{m-1} u_k e^{-ik\psi}\right] = \operatorname{Re}\left[e^{i\theta} Z\right].$$
:::

Tulis $Z = |Z|\, e^{i\phi}$. Maka $\operatorname{Re}[e^{i\theta} Z] = |Z| \cos(\theta + \phi) = |Z| \cos(\theta - \delta)$ dengan $\delta = -\phi = -\arg Z$. Dengan $C = c\sum_k u_k$ dan $R = a|Z|$, diperoleh $y_u = s[C + R\cos(\theta - \delta)]$. ∎

Untuk $m = 2$ dan $u = (\cos\alpha, \sin\alpha)$ berlaku $C^2 = c^2(\cos\alpha + \sin\alpha)^2 = c^2(1 + \sin 2\alpha)$ dan

::: {custom-style="PersamaanA"}
$$|Z|^2 = (\cos\alpha + \sin\alpha\cos\psi)^2 + (\sin\alpha\sin\psi)^2 = 1 + \sin 2\alpha\cos\psi.$$
:::

## A.3 Bukti Proposisi 2

Tulis $X = \cos(\theta - \delta)$. Karena $s$ bebas dari $\theta$, $\mathbb{E}[s] = 0$ dan $s^2 = s^4 = 1$, diperoleh $\mathbb{E}[y] = \mathbb{E}[s]\,\mathbb{E}[C + RX] = 0$. Dengan momen pada A.1, suku berpangkat ganjil dari $X$ hilang, sehingga $\mu_2 = \mathbb{E}[(C + RX)^2] = C^2 + R^2/2$ dan

::: {custom-style="PersamaanA"}
$$\mu_4 = \mathbb{E}\left[(C + RX)^4\right] = C^4 + 6C^2R^2\cdot\frac{1}{2} + R^4\cdot\frac{3}{8} = C^4 + 3C^2R^2 + \frac{3R^4}{8}.$$
:::

Distribusi $y$ simetris, sehingga kemencengannya nol, dan untuk $N$ besar koreksi Sarle menuju 3, sehingga $I_3 \to 1/(\kappa + 3) = \mu_2^2/\mu_4$. Dengan $R^2 = \rho C^2$,

::: {custom-style="PersamaanA"}
$$b(\rho) = \frac{C^4(1 + \rho/2)^2}{C^4(1 + 3\rho + 3\rho^2/8)} = \frac{2(2+\rho)^2}{3\rho^2 + 24\rho + 8}.$$
:::

Dengan pembilang $N(\rho) = (1 + \rho/2)^2$ dan penyebut $M(\rho) = 1 + 3\rho + 3\rho^2/8$,

::: {custom-style="PersamaanA"}
$$N'M - NM' = \left(1 + \frac{\rho}{2}\right)\left[M - \left(1 + \frac{\rho}{2}\right)\left(3 + \frac{3\rho}{4}\right)\right] = \left(1 + \frac{\rho}{2}\right)\left(\frac{3\rho}{4} - 2\right),$$
:::

yang setara dengan Persamaan (8). Untuk $\rho \ge 0$ faktor $\rho + 2$ positif, sehingga tanda $b'$ sama dengan tanda $3\rho - 8$: $b$ turun pada [0, 8/3) dan naik sesudahnya. Jadi $\rho = 8/3$ adalah minimum global pada [0, ∞), dengan $b(8/3) = 2(14/3)^2/(64/3 + 72) = 7/15$. Limitnya adalah $b(0) = 1$ dan $b(\rho) \to 2/3$ untuk $\rho \to \infty$. ∎

## A.4 Bukti Akibat 1

Dari A.2, $\rho = R^2/C^2 = (a/c)^2 (1 + \sin 2\alpha\cos\psi)/(1 + \sin 2\alpha)$ dengan $\psi = 2\pi W/P$ untuk $m = 2$, sehingga $L(W)$ bergantung pada $W$ hanya melalui $\cos\psi$. Turunannya terhadap $W$ memuat faktor $\sin(2\pi W/P)$, yang nol di $W = P/2$, sehingga titik itu stasioner. Di $W = P$, $\cos\psi = 1$ dan $\rho = (a/c)^2$ untuk semua $\alpha$. Di $W = P/2$, $\cos\psi = -1$ dan

::: {custom-style="PersamaanA"}
$$\rho = \left(\frac{a}{c}\right)^2 \frac{1 - q}{1 + q}, \qquad q = \sin 2\alpha.$$
:::

Nilai $q$ berada pada selang (−1, 1]. Untuk setiap $k \ge 0$, nilai $q = (1 - k)/(1 + k)$ berada di selang itu dan memberi $\rho = k\,(a/c)^2$, sehingga $\rho$ menyapu seluruh [0, ∞). Ambang 0,4573 diperoleh dengan menyelesaikan $L(P/2) = b\left((a/c)^2\right)$ secara numerik (metode Brent, rata-rata atas $2 \times 10^5$ nilai $\alpha$). Penyisiran *a*/*c* dari 0,30 sampai 1,00 tidak menemukan maksimum lokal lain di dalam $(0, P/2)$. Untuk *a*/*c* jauh di atas 1, di luar wilayah sah model, urutan itu berubah lagi.

Argumen simetri ini berlaku untuk sembarang $m$. Menurut Proposisi 1, $C$ tidak bergantung pada $\tau$, dan karena $u$ real, $Z(-\psi) = \overline{Z(\psi)}$ sehingga $|Z|$ tidak berubah oleh $\psi \to -\psi$. Bersama periodisitas $2\pi$ dalam $\psi$, lanskap tanpa perpindahan memenuhi $L(\tau) = L(P - \tau)$ dan karena itu stasioner di $\tau = P/2$, yaitu di $W = (m-1)P/2$. Penyisiran numerik dengan 40.000 arah acak menunjukkan bahwa titik ini adalah maksimum global untuk *a*/*c* antara 0,6 dan 1,0 pada $m = 3$ dan $m = 4$.

## A.5 Bukti Proposisi 3

Untuk pasangan beda sayap, $v = s\left(c + a\cos\theta,\; -[c + a\cos(\theta - \psi)]\right)^\top$, sehingga

::: {custom-style="PersamaanA"}
$$y' = s\left[c(u_0 - u_1) + a\left(u_0\cos\theta - u_1\cos(\theta - \psi)\right)\right].$$
:::

Ini adalah bentuk Proposisi 1 dengan $u_1$ diganti $-u_1$. Langkah A.2 kemudian memberi $C'^2 = c^2(1 - \sin 2\alpha)$ dan $R'^2 = a^2(1 - \sin 2\alpha\cos\psi)$. Di $W = P$ berlaku $\psi = 2\pi$ dan $\cos\psi = 1$. Dengan $q = \sin 2\alpha$, $K_2 = c^2 + a^2/2$ dan $K_4 = c^4 + 3c^2a^2 + 3a^4/8$, Proposisi 2 memberi

::: {custom-style="PersamaanA"}
$$\mu_2^{\text{sama}} = K_2(1+q), \qquad \mu_4^{\text{sama}} = K_4(1+q)^2, \qquad \mu_2^{\text{beda}} = K_2(1-q), \qquad \mu_4^{\text{beda}} = K_4(1-q)^2,$$
:::

dengan label sama untuk pasangan sesayap dan label beda untuk pasangan beda sayap. Pasangan $(x(t), x(t-P))$ selalu melintasi tepat satu batas putaran, dan perpindahan di batas itu terjadi dengan peluang $1/D$ secara bebas dari $\theta$. Kedua komponen bererata nol, sehingga momen pusat campurannya adalah rata-rata berbobot. Untuk momen kedua, $\mu_2 = p\,\mu_2^{\text{sama}} + (1-p)\,\mu_2^{\text{beda}} = K_2\left[1 + (2p-1)q\right]$, dan untuk momen keempat

::: {custom-style="PersamaanA"}
$$\mu_4 = K_4\left[p(1+q)^2 + (1-p)(1-q)^2\right] = K_4\left[1 + 2(2p-1)q + q^2\right].$$
:::

Rasio $\mu_2^2/\mu_4$ memberi Persamaan (10), karena $K_2^2/K_4 = b\left((a/c)^2\right)$ menurut Proposisi 2. Untuk $q = -1$, pembilangnya $[1 - (2p-1)]^2 = 4(1-p)^2$ dan penyebutnya $4(1-p)$, sehingga $b_{\text{mix}} = (1-p)\,b\left((a/c)^2\right)$. Untuk $q = 1$, dengan cara yang sama, $b_{\text{mix}} = p\,b\left((a/c)^2\right)$. ∎

Sebagai contoh, untuk *D* = 40 berlaku *p* = 39/40, sehingga $b_{\text{mix}}(3\pi/4) = b\left((a/c)^2\right)/40$. Untuk *a*/*c* = 0,6 nilainya 0,654/40 ≈ 0,0164.

# Lampiran B. Detail numerik

: Tabel B1. Pengaturan setiap fase. Banyak arah acak adalah 600 pada Fase 1 dan 500 pada fase lain.

| Fase | Sistem atau model | Data per realisasi | Grid jendela | Realisasi | Benih |
|:--|:--|:--|:--|:--:|:--:|
| 1 | Lorenz *ρ* = 28 | 75 s dan 150 s, langkah 0,01 | *W* = 0,08 sampai 1,28 s, langkah 0,08 | 10 | 2026 |
| 2 | Lorenz (tiga *ρ*), Rössler, Hénon | 20.000 sampel; langkah 0,005 (Lorenz), 0,05 (Rössler) | *W*/*P* = 0,06 sampai 2,0, langkah 0,06; Hénon 1 sampai 24 langkah | 5 | 2027 |
| 3 | dua sayap, *a*/*c* = 0,6 | 200.000 (*D* = 40) dan 60.000 (*D* = 2,2) sampel, 100 sampel per putaran | *W*/*P* = 0,02 sampai 2,0, langkah 0,02 | 1 | 31 |
| 3b | amplitudo tetap, spiral, urutan acak | 60.000 sampel | *W*/*P* = 0,03 sampai 1,6, langkah 0,03 | 3 | 2718 |
| 3c | tangga model A sampai E | 60.000 sampel | *W*/*P* = 0,03 sampai 1,6, langkah 0,03 | 3 | 1618 |

Indeks dihitung pada sekitar 4000 titik per proyeksi (subsampel teratur) mulai Fase 2, dan pada setiap titik ketiga di Fase 1. Arah acak dibangkitkan sebagai vektor Gaussian baku yang dinormalisasi, sehingga tersebar seragam di bola, dan himpunan arah yang sama dipakai untuk semua $\tau$ dalam satu eksperimen. Transien yang dibuang adalah 20 s (Fase 1), 30 s (Lorenz Fase 2), 300 s (Rössler) dan 1000 iterasi (Hénon). Kondisi awal Lorenz diundi dari distribusi normal dengan rerata (1, 1, 20) dan simpangan baku 5, sedangkan kondisi awal Rössler memakai rerata (1, 1, 0) dan simpangan baku 0,5. Integrasi memakai `solve_ivp` (RK45) dengan toleransi $10^{-8}$ pada Fase 1 dan $10^{-9}$ pada fase lain. Kemencengan dan kurtosis berlebih dihitung dari momen sampel tanpa koreksi bias; untuk sekitar 4000 titik, selisihnya terhadap versi terkoreksi pada Persamaan (3) dapat diabaikan.

Lokasi puncak di Fase 1 diambil sebagai maksimum setelah perataan bergerak tiga titik. Mulai Fase 2, lokasi puncak diambil dengan interpolasi parabola tiga titik, dan puncak di tepi grid ditandai. Bentuk gelombang di Fase 3c diukur dari satu lintasan Lorenz *ρ* = 28 sepanjang 1000 s (langkah 0,005) yang berawal di (1, 1, 20), dengan 200 titik fase per putaran. Statistik bentuk gelombang di Bagian 4.5 dihitung dari lima lintasan sepanjang itu, yaitu lintasan tersebut dan empat lintasan dengan kondisi awal acak. Semua perhitungan memakai Python dengan NumPy, SciPy dan Matplotlib, dan semua berkas CSV dibaca dan ditulis dengan pengodean UTF-8 secara eksplisit.
