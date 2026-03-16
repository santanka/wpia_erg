# ドキュメント: Shoji et al. 2021 図の再現 (Takeuchi手法)

このドキュメントは、Jupyter Notebook `Shoji_et_al_2021_fig_3_reproduce_Takeuchi.ipynb` および `Shoji_et_al_2021_fig_4_reproduce_Takeuchi.ipynb` で実装されているデータ処理手順の概要と解説をまとめたものです。

## 概要
本解析の目的は、Shoji et al. (2021) で提示された波動粒子相互作用解析（WPIA）の結果（特に図3および図4）を再現することです。座標系の取り扱いやフラックス処理には「Takeuchi手法」を採用しています。

## データ処理手順

### 1. 環境構築
- `ERG_DATA_DIR`: Arase（あらせ）衛星データのローカル保存パスを設定。
- 使用ライブラリ: `pyspedas`, `ergpyspedas`, `numpy`, `xarray`, `scipy`。

### 2. 3Dフラックスデータ処理 (LEP-i)
- LEP-iインストゥルメントから `3dflux` データをロード。
- 各エネルギー・チャネルごとにフラックスを時系列データとして再構成（フラット化）。
- スピン位相やエネルギー・ステップのオフセットを考慮し、正確な観測時刻を決定。

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
- **作成者**: [名前を入力してください]
- **最終更新日**: 2026-03-16
- **所属/プロジェクト**: WPIA ERG プロジェクト

## 参考文献
- Shoji, M., et al. (2021). [論文タイトル/DOI等を入力]
