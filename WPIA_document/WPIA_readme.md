<!-- omit in toc -->
# 波動粒子相互作用解析(WPIA)コード開発メモ

<!-- omit in toc -->
## 概要
本書は、非線形波動粒子相互作用の直接観測を解析するためのコードである Wave–Particle Interaction Analyzer（WPIA）コードの実装方法をまとめたものである。本書の目的は、Shoji et al. (2021)[^1] において示された、Electromagnetic ion cyclotron wave（EMIC）と陽子の間で生じた非線形波動粒子相互作用について、Arase衛星による直接観測結果（Figures 3 and 4）を再現することである。


<!-- omit in toc -->
## 目次
- [I. WPIAの仕組み](#i-wpiaの仕組み)
  - [1. 波動粒子相互作用とエネルギー輸送](#1-波動粒子相互作用とエネルギー輸送)
  - [2. 共鳴条件付近の粒子分布](#2-共鳴条件付近の粒子分布)
- [II. WPIA手順](#ii-wpia手順)
  - [1. 環境構築](#1-環境構築)
  - [2. 粒子データの時系列再構築と衛星座標系における観測機器の視線方向ベクトルの導出](#2-粒子データの時系列再構築と衛星座標系における観測機器の視線方向ベクトルの導出)
  - [3. 座標変換](#3-座標変換)
  - [4. 磁場解析 (MGF)](#4-磁場解析-mgf)
  - [5. 波動解析](#5-波動解析)
  - [6. 粒子位相角の算出](#6-粒子位相角の算出)
  - [7. フラックスからカウント数への変換](#7-フラックスからカウント数への変換)
  - [8. 電場解析 (PWE-EFD)](#8-電場解析-pwe-efd)
  - [9. WPIA物理量](#9-wpia物理量)
  - [10. 可視化](#10-可視化)
- [文責](#文責)
- [謝辞](#謝辞)

## I. WPIAの仕組み

### 1. 波動粒子相互作用とエネルギー輸送
波動粒子相互作用の本質は、電磁場エネルギーの輸送とプラズマ粒子の運動エネルギー変化の間のエネルギー交換過程にある。

電場$\mathbf{E}$と磁場$\mathbf{B}$の関係式であるFaraday's law of induction $\nabla \times \mathbf{E} = - \partial\mathbf{B} / \partial t$とAmpère's circuital law $\nabla \times \mathbf{B} = \mu_{0} \mathbf{J} + \varepsilon_{0} \mu_{0} \, \partial\mathbf{E} / \partial t$から、以下の電磁場エネルギーに関する輸送方程式が導かれる。
$$
\frac{\partial u}{\partial t} + \nabla \cdot \mathbf{S} = - \mathbf{E} \cdot \mathbf{J}. \tag{1}
$$
ここで、$\varepsilon_{0}$と$\mu_{0}$、$\mathbf{J}$はそれぞれ真空の誘電率および透磁率、電流密度である。$u$は電磁場エネルギー密度、$\mathbf{S}$はPoynting fluxであり、それぞれ以下のように定義される。
$$
\begin{align*}
  u &:= \frac{1}{2} \left( \varepsilon_{0} \left| \mathbf{E} \right|^{2} + \frac{1}{\mu_{0}} \left| \mathbf{B} \right|^{2} \right),   \cr
  \mathbf{S} &:= \frac{1}{\mu_{0}} \left( \mathbf{E} \times \mathbf{B} \right)
\end{align*}
$$

一方、プラズマ側のエネルギー変化を考える。単一粒子の運動エネルギー $K$の時間変化は、ローレンツ力のうち電場がなす仕事 $q \mathbf{v} \cdot \mathbf{E}$で与えられる。系全体の運動エネルギー変化率は、速度分布関数 $f_{s} \left(\mathbf{r}, \mathbf{v}, t \right)$を用いてそのモーメントをとることで、以下のように記述できる。
$$
\begin{align*}
  \sum_{s} \iiint \frac{\mathrm{d} K_{\mathrm{s}}}{\mathrm{d} t} \, f_{s} \left( \mathbf{r}, \mathbf{v}, t \right) \, \mathrm{d}^{3} \mathbf{v} & = \mathbf{E} \cdot \sum_{s} q_{s} \iiint \mathbf{v} \, f_{s} \left( \mathbf{r}, \mathbf{v}, t \right) \, \mathrm{d}^{3} \mathbf{v} \cr
  & = \mathbf{E} \cdot \sum_{s} q_{s} n_{s} \mathbf{V}_{s} \cr
  & = + \mathbf{E} \cdot \mathbf{J} . \tag{2}
\end{align*}
$$
ここで、$n_{s}$と $\mathbf{V}_{s}$はプラズマ種$s$の数密度および平均流速である。

<figure align="center">
  <img src="./images/energy_transfer.png" alt="Figure 1" width='500'>
  <figcaption>図1. 電磁場エネルギーの輸送と粒子の運動エネルギー変化の概念図</figcaption>
</figure>

式(1)と式(2)を比較すると、項 $\mathbf{E} \cdot \mathbf{J}$ (単位時間、単位体積あたりのジュール熱)が電磁場とプラズマの間のエネルギー交換を媒介していることがわかる。図1に示すように、このエネルギー輸送は以下のように解釈できる。

- エネルギーの連続性:
もし電流が存在しない($\mathbf{J} = \mathbf{0}$)場合、式(1)は $\partial u / \partial t + \nabla \cdot \mathbf{S} = 0$となり、エネルギー流束の流入($\mathbf{S}_{\mathrm{in}}$)・流出($\mathbf{S}_{\mathrm{out}}$)のみで閉じられた単純な連続の式となる。
- 相互作用の鍵としての $\mathbf{E} \cdot \mathbf{J}$
電流が生じている領域では、この項が場と粒子のエネルギーを結合させる。
  - $\mathbf{E} \cdot \mathbf{J} > 0$の場合: 電磁場が粒子に仕事をなし、エネルギーは**電磁場からプラズマへ**と受け渡される (粒子の加速)。
  - $\mathbf{E} \cdot \mathbf{J} < 0$の場合: 粒子が電磁場に対して仕事をなし、エネルギーは**プラズマから電磁場へ**と受け渡される (波動の励起)。

巨視的な視点からは、波動粒子相互作用は式(1)と式(2)に現れる結合項 $\mathbf{E} \cdot \mathbf{J}$を通じた、場と粒子の間のエネルギー遷移として記述される。一方、次節で詳述するように、このプロセスは微視的には、特定の共鳴条件を満たす粒子群が波動の位相に同期し、電磁場から正味の仕事を受ける(あるいは場に仕事をなす)過程として解釈される。

### 2. 共鳴条件付近の粒子分布


## II. WPIA手順

波動粒子相互作用における場と粒子のエネルギー交換、および波動の周波数変化を記述する指標として、物理量$W_{\mathrm{Eint}}$[^2][^3]と $W_{\mathrm{Bint}}$[^4]を、以下のように定義する。
$$
\begin{align*}
  W_{\mathrm{Eint}} &:= \sum_{i=1}^{N} \mathbf{v}_{\perp i} \cdot \mathbf{E}_{\mathrm{w}i} \tag{3}, \\
  W_{\mathrm{Bint}} &:= \sum_{i=1}^{N} \mathbf{v}_{\perp i} \cdot \mathbf{B}_{\mathrm{w}i}. \tag{4}.
\end{align*}
$$
WPIAでは、これらの指標の算出に加え、位相空間 $\left( K, \alpha, \zeta \right)$ (ここで $K$: 運動エネルギー、$\alpha$: ピッチ角、$\zeta := \phi - \psi$、$\phi$: ジャイロ位相、$\psi$: 波の位相) 上における粒子分布の可視化を行う。
以下に、Arase衛星のLEP-iデータを用いた具体的な解析フローを示す。

### 1. 環境構築
解析には、SPEDASのPython実装である`pyspedas`およびArase衛星用プラグイン`ergpyspedas`を基盤とし、データ構造の操作に`xarray`、`numpy`、`scipy`を用いる。
```bash
# 依存ライブラリのインストール
python -m pip install pyspedas numpy xarray scipy matplotlib
# ERG-SC公式プラグイン（最新開発版）のインストール
python -m pip install --upgrade --force-reinstall git+https://github.com/ergsc-devel/pyspedas_plugin.git
```

### 2. 粒子データの時系列再構築と衛星座標系における観測機器の視線方向ベクトルの導出
LEP-i[^5]から陽子の3D Fluxデータをロードする。(ここではShoji et al. (2021)[^1]で使用した版のデータを使用する。) WPIAの実行に際し、スピン周期(8 sec)でまとめられたデータを、個々の粒子カウントが実際に発生した観測時刻へと展開する必要がある。
```python
import pyspedas as psp
import ergpyspedas as ergpy
import xarray as xr

time_range = ['2017-11-15/16:10:00', '2017-11-15/16:27:00']
ergpy.lepi(time_range, datatype='3dflux', get_support_data=True, version='v03_00')

# flux3d本体 (dims = time, v1(energy), v2(channel number), v3(spin phase))
flux3d = psp.get_data('erg_lepi_l2_3dflux_FPDU', xarray=True).sel(
    time=slice(*time_range)
)

# 各軸の座標
time_ax     = flux3d.time.values        # shape = (T,)
energy_ax   = flux3d.v1.values          # (E,)
channel_ax  = flux3d.v2.values          # (C,)
spin_ax     = flux3d.v3.values          # (S,)

time_num, energy_num, channel_num, spin_num = len(time_ax), len(energy_ax), len(channel_ax), len(spin_ax)   # T, E, C, S
```

<figure align="center">
  <img src="./images/Asamura_2018_Fig_1.png" alt="Figure 2" height='200'>
  <figcaption>図2. Location of LEPi on ERG satellite. (According to Asamura et al. (2018).)</figcaption>
</figure>

<figure align="center">
  <img src="./images/Asamura_2018_Fig_2.png" alt="Figure 3" height='400'>
  <figcaption>図3. Channel definition of LEPi. LEPi has a planar FOV in the Y_sc‒Z_sc plane. Illustrated direction corresponds to the velocity direction of the incoming particles. Channels 0-8, and e are wide channels, while channels 9-d are narrow channels. The positions of the center of the anodes in the satellite frame of reference are shown in parentheses. 本WPIAでは、Channel 0‒8を使用する。 (According to Asamura et al. (2018).)</figcaption>
</figure>

<figure align="center">
  <img src="./images/Asamura_2018_Fig_9.png" alt="Figure 4" height='300'>
  <figcaption>図4. Measurement timing of LEPi with the index pulse. Spin phase is determined by the reception timing of the index pulse. (According to Asamura et al. (2018).)</figcaption>
</figure>

<figure align="center">
  <img src="./images/Katoh_2018_Fig_2.png" alt="Figure 5" height='700'>
  <figcaption>図5. 観測の参考イメージとして、MEP-eでの観測ステップを引用。 (a) Variation of the pitch angle at the center of the field-of-view for sensor channels of the MEP-e during one spin period of 8 sec under the assumed condition. Schematics show in the upper panel represent the FOV of each sensor channel every 2 sec, where the color of rectangles corresponds to those of plotted lines. (b) Energy range measured by MEP-e during one spin period, where 16 energy steps are swept every 0.25 sec. (According to Katoh et al. (2018).)</figcaption>
</figure>

Asamura et al. (2018) [^6] によれば、LEP-iは1スピン（8秒）の間に16のスピン位相（0.5秒/ステップ）をもち、各ステップ内でエネルギーを32段階（15.625ミリ秒/ステップ）で掃引することで、速度空間の全容を捉える(図2-4[^6]と図5[^7]参照)。WPIAでは、この観測サイクルを考慮し、3D Fluxデータの各ビンを実際の観測時刻へ再割り当てする。併せて、衛星座標系(Spinning satellite Geometry Axis; SGA coordinate system)[^8]における各チャンネルの視線方向(単位速度ベクトル)を算出し、$\mathbf{E} \cdot \mathbf{J}$の計算に備える。
```python
# 1 spin time = 8 sec
# 1 spin phase time = 0.5 sec
# 1 energy time step = 15625 μsec
spin_offset_ns      = np.arange(spin_num, dtype='timedelta64[ns]') * 500_000_000        # 0.5 sec = 500,000,000 nsec
energy_offset_ns    = (np.arange(energy_num, dtype='int64') * 15_625_000 + 7_812_500).astype('timedelta64[ns]')

offset_ns           = energy_offset_ns[:, None] + spin_offset_ns    # (E, 1) + (S) -> (E, S)

flux_E_TS_C = flux3d.transpose('v1_dim', 'time', 'v3_dim', 'v2_dim').values.reshape(energy_num, time_num*spin_num, channel_num)

# xarray.DataArrayをエネルギーごとに生成
flux_data_arrays = {}
for energy_i in range(energy_num):
    time_flat = (time_ax[:, None] + offset_ns[energy_i][None, :]).reshape(-1) # ((T, 1) + (1, S) -> (T, S)).reshape(-1) -> (TxS,)

    flux_data_arrays[energy_i] = xr.DataArray(
        flux_E_TS_C[energy_i],
        dims=['time', 'channel'],
        coords={
            'time':         time_flat,
            'channel':      channel_ax,
            'energy_keV':   energy_ax[energy_i]
        },
        attrs=flux3d.attrs,
        name=f'flux_energy_{energy_i}'
    )

print(flux_data_arrays)

fidu_angle_dict     = psp.get_data('erg_lepi_l2_3dflux_FIDU_Angle_sga', xarray=True)
fidu_angle  = fidu_angle_dict.astype(float)     # shape (2, 3, 16)
AZ_deg_mid = fidu_angle[0, 1, :channel_num]      # shape (C,)
# AZ_deg_mid =  [ 78.75  56.25  33.75  11.25 -11.25 -33.75 -56.25 -78.75]

theta_sga = AZ_deg_mid                  # (C,)
varphi_sga = -90. * np.ones(spin_num)   # (S,)

# (TxS, C)の2次元配列を生成
theta_sga_time = np.tile(theta_sga, (time_num*spin_num, 1))   # (TxS, C)
varphi_sga_times = np.tile(varphi_sga, (time_num, 1)).reshape(-1, 1)   # (TxS, 1)
varphi_sga_time = np.tile(varphi_sga_times, (1, channel_num))   # (TxS, C)

angle_sga_time = np.stack((theta_sga_time, varphi_sga_time), axis=2)   # (TxS, C, 2)

# theta_sga, varphi_sga -> Vx_sga, Vy_sga, Vz_sga
vector_sga_time = np.zeros((time_num*spin_num, channel_num, 3))   # (TxS, C, 3)
vector_sga_time[:, :, 0] = np.cos(np.radians(angle_sga_time[:, :, 0])) * np.cos(np.radians(angle_sga_time[:, :, 1]))
vector_sga_time[:, :, 1] = np.cos(np.radians(angle_sga_time[:, :, 0])) * np.sin(np.radians(angle_sga_time[:, :, 1]))
vector_sga_time[:, :, 2] = np.sin(np.radians(angle_sga_time[:, :, 0]))

v_unit_vector_sga_energy_channel_list = {}
for energy_i in range(energy_num):
    time_flat = (time_ax[:, None] + offset_ns[energy_i][None, :]).reshape(-1)

    for channel_i in range(channel_num):
        v_unit_vector_sga_energy_channel_list[energy_i, channel_i] = xr.DataArray(
            vector_sga_time[:, channel_i, :],
            dims=['time', 'xyz'],
            coords={
                'time': time_flat,
                'xyz': ['x', 'y', 'z'],
                'channel': channel_ax[channel_i],
                'energy_keV':   energy_ax[energy_i]
            },
            name=f'v_unit_vector_sga_{energy_i}_{channel_i}'
        )
        dot_ = (v_unit_vector_sga_energy_channel_list[energy_i, channel_i] * v_unit_vector_sga_energy_channel_list[energy_i, channel_i]).sum(dim='xyz')
        print(f'v_unit_vector_sga_energy_channel_list[{energy_i}, {channel_i}] = ', np.nanmin(dot_), np.nanmax(dot_), np.nanmean(dot_))
```


### 3. 座標変換
粒子検出器のユニットベクトルを以下の座標系間で変換：
- **SGA (Spin-linked Geocentric Attitude)**
- **SGI (Spin-linked Geocentric Inertial)**
- **DSI (Despun Geocentric Inertial)**: 解析のメインとなる座標系。

### 4. 磁場解析 (MGF)
- DSI座標系の MGF `256hz` および `8sec` データをロード。
- 移動平均（デフォルト窓幅：100秒）を用いて背景磁場 ($\mathbf{B}_0$) を算出。
- 背景磁場からの変動成分（垂直成分 $\mathbf{B}_{\perp}$）を抽出。

### 5. 波動解析
- **スペクトログラム**: $B_z$ 成分のパワースペクトル密度（PSD）を計算。
- **フィルタリング**: バターワース・バンドパスフィルタ（例：0.45 - 0.75 Hz）を適用し、特定の波動モード（EMIC波など）を抽出。

### 6. 粒子位相角の算出
- **ピッチ角** ($\alpha$): 粒子速度ベクトルと背景磁場 $\mathbf{B}_0$ とのなす角。
- **ゼータ角** ($\zeta$): 垂直波動磁場に対する粒子の位相角。

### 7. フラックスからカウント数への変換
- 微分数フラックス ($J$) を以下の式でカウント数 ($C$) に変換：
  $$C = \tau \cdot \epsilon \cdot G \cdot E \cdot J$$
  ここで $\tau$ はサンプリング時間、$\epsilon$ は検出効率、$G$ は幾何学的因子、$E$ はエネルギー。

### 8. 電場解析 (PWE-EFD)
- PWE-EFD `256hz` データをロード。
- $\mathbf{E} \cdot \mathbf{B} = 0$ を仮定して $E_z$ 成分を推定。
- 垂直電場変動成分 ($\mathbf{E}_{\perp}$) を抽出。

### 9. WPIA物理量
エネルギー交換を評価するための $W_{Eint}$ および $W_{Bint}$ を算出：
- $W_{Eint} \propto \frac{q E_w C}{\epsilon G} \sin^2 \alpha$
- $W_{Bint} \propto \frac{q B_w C}{\epsilon G} \sin^2 \alpha$

### 10. 可視化
- 磁場のスペクトログラムのプロット。
- $(\alpha, \zeta)$ 空間における粒子分布の可視化。
- 時刻・エネルギーごとの $W_{Eint}$ および $W_{Bint}$ のサマリープロット。

## 文責
- **作成者**: 齋藤 幸碩 ([Researchmap](https://researchmap.jp/koseki_saito))
- **最終更新日**: 2026/03/17
- **所属・役職**: 東北大学大学院理学研究科地球物理学専攻　特任研究員

## 謝辞
本書は、科学研究費助成事業 基盤研究(S) 23H05429『惑星放射線帯消失モデルの実証と能動的制御方法の開拓』の一環として作成された。作成にあたっては、東北大学大学院理学研究科地球物理学専攻の加藤雄人教授、ならびに名古屋大学大学院工学研究科電気工学専攻の竹内亘平様より、多くの議論と貴重なコメントをいただいた。ここに記して感謝申し上げる。


[^1]: Shoji, M., Y. Miyoshi, L. M. Kistler, K. Asamura, A. Matsuoka, Y. Kasaba, S. Matsuda, Y. Kasahara, and I. Shinohara (2021), ``Discovery of proton hill in the phase space duiring interactions between ions and electromagnetic ion cyclotron waves'', *Scientific Reports **11***, 13480, doi:[10.1038/s41598-021-92541-0](https://doi.org/10.1038/s41598-021-92541-0).
[^2]: Fukuhara, H., H. Kojima, Y. Ueda, Y. Omura, Y. Katoh, and H. Yamakawa (2009), ``A new instrument for th study of wave-particle interactions in space: One-chip Wave-Particle Interaction Analyzer'', *Earth, Planets and Space **61***(6), 765‒778, doi:[10.1186/BF03353183](https://doi.org/10.1186/BF03353183).
[^3]: Katoh, Y., M. Kitahara, H. Kojima, Y. Omura, S. Kasahara, M. Hirahara, Y. Miyoshi, K. Seki, K. Asamura, T. Takashima, and T. Ono (2013), ``Significance of Wave-Particle Interaction Analyzer for direct measurements of nonlinear wave-particle interactions'', *Annales Geophysicae **31***(3), 503‒512, doi:[10.5194/angeo-31-503-2013](https://doi.org/10.5194/angeo-31-503-2013).
[^4]: Shoji, M., Y. Miyoshi, Y. Katoh, K. Keika, V. Angelopoulos, S. Kasahara, K. Asamura, S. Nakamura, and Y. Omura (2017), ``Ion hole formation and nonlinear generation of electromagnetic ion cyclotron waves: THEMIS observations'', *Geophysical Research Letters **44***(17), 8,730‒8,738, doi:[10.1002/2017GL074254](https://doi.org/10.1002/2017GL074254).
[^5]: ERG Science Center (ERG-SC) wiki/ErgSat/Lepi https://ergsc.isee.nagoya-u.ac.jp/mw/index.php/ErgSat/Lepi, Accessed 2026-03-18.
[^6]: Asamura, K., Y. Kazama, S. Yokota, S. Kasahara, and Y. Miyoshi (2018), ``Low-energy particle experiments‒ion mass analyzer (LEPi) onboard the ERG (Arase) satellite'', *Earth, Planets and Space **70***, 70, doi:[10.1186/s40623-018-0846-0](https://doi.org/10.1186/s40623-018-0846-0).
[^7]: Katoh, Y., H. Kojima, M. Hikishima, T. Takashima, K. Asamura, Y. Miyoshi, Y. Kasahara, T. Mitani, N. Higashio, A. Matsuoka, M. Ozaki, S. Yagitani, S. Yokota, S. Matsuda, M. Kitahara, and I. Shinohara (2018), ``Software-type Wave‒Particle Interaction Analyzer on board the Arase satellite'', *Earth, Planets and Space **70***, 4, doi:[10.1186/s40623-017-0771-7](https://doi.org/10.1186/s40623-017-0771-7).
[^8]: ERG Science Center (2020), ``Definition of science coordinate systems for the Arase satellite'', https://ergsc.isee.nagoya-u.ac.jp/assets/howto/ERG_Coordinate_System_202004.pdf, Accessed 2026-03-18.

<!--
- Takeuchi, K., Y. Miyoshi, K. Asamura, K. Terasawa, C.-W. Jun, Y. Kasahara, Y. Kasaba, S. Matsuda, F. Tsuchiya, A. Kumamoto, T. Hori, A. Shinbori, A. Matsuoka, M. Teramoto, K. Yamamoto, I. Shinohara, and N. Kitamura (2025), ``Direct Measurement of Pitch Angle Scattering in EMIC-Proton Interactions by the WPIA method'', *Japan Geoscience Union Meeting 2025*, PEM13-P16.
-->