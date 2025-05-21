# 都市モデル最適化ツール

本リポジトリは、都市計画・立地モデルに基づく空間経済の解析を目的とした最適化ツールです。C言語で実装されたメインのロジックに加え、Pythonスクリプトによる補助的な処理・可視化が含まれています。

## 📁 ディレクトリ構成

```
C/
├── .vscode/ # VSCode用の設定ファイル
├── build/Debug/ # ビルド済みファイルの出力先
├── data/ # 入力データ（CSVやパラメータファイルなど）
├── src/ # ソースコード本体
│ ├── core/ # 数値計算に用いるファイル。使用する関数が全て1つのファイルに入っている。
│ ├── function/ # 各関数を、それぞれのファイルに保存したもの。関数のテスト用。
│ └── test_code/ # 検証用のコード・ユニットテスト等
├── run_with_log.bat # 実行時ログを出力するバッチファイル
```

python/ # Pythonスクリプト（前処理・後処理・可視化など）

.gitignore # Git無視ファイル設定
README.md # このファイル

## 🛠 ビルド方法（C言語）

本ツールはCMake等は未使用で、手動またはバッチビルドを想定しています。

### Windows + VSCode環境

1. `.vscode/` の `tasks.json` に基づいて `Ctrl + Shift + B` でビルド可能。
2. 出力ファイルは `build/Debug/` に生成されます。

## 🚀 実行方法

### バッチファイル経由（Windows）

```bash
run_with_log.bat
```

## 🧠 解いている最適化問題の概要

本ツールは、都市空間における人口分布・地価・賃金が内生的に決定される **短期均衡モデル** を対象とし、以下の目的関数を最小化することで均衡解を求めます。

### 🎯 目的関数（ZSD）

![image](https://github.com/user-attachments/assets/2f14e05d-b1c0-4a09-ab15-8ffa95b386ad)

- \( R_i \)：地点 \(i\) における地価  
- \( W_j \)：地点 \(j\) における賃金  
- \( m_i \)：企業の立地分布  
- \( n_{ij} \)：住民の通勤パターン  
- \( \pi_i, v_{ij} \)：それぞれ企業利潤・家計効用  
- \( \bar{S}_i \)：土地供給量

土地市場・労働市場の均衡条件を制約として課しています。

### ⚖ 効用関数の例（べき乗型）

\[
f(S_{ij}^H, L_{ij}^H) = \frac{(S_{ij}^H)^{\alpha_1} (L_{ij}^H)^{\alpha_2}}{t T_{ij}}, \quad
g(S_i^F, L_i^F) = (S_i^F)^{\beta_1} (L_i^F)^{\beta_2}
\]

---

## ⚙ アルゴリズム（FISTA）

本モデルの最適化は、**FISTA（Fast Iterative Shrinkage-Thresholding Algorithm）** をベースにした反復法により解かれます。

### 🔁 アルゴリズムの概要

1. 初期化：\( L_0, \epsilon, \eta \) などのパラメータ設定  
2. Backtracking によりステップ長を調整  
3. 勾配降下による更新：

\[
p_k^C = \bar{p}_k^C - \frac{1}{L_k} \nabla ZSD(\bar{p}_k^C)
\]

4. 収束判定：

\[
\left| \frac{p_k^C(i) - p_{k-1}^C(i)}{p_{k-1}^C(i)} \right| < \epsilon \quad \forall i
\]

5. モメンタム項による次ステップの生成：

\[
t_{k+1} = \frac{1 + \sqrt{1 + 4 t_k^2}}{2}, \quad
\bar{p}_{k+1}^C = p_k^C + \frac{t_k - 1}{t_{k+1}} (p_k^C - p_{k-1}^C)
\]

6. 勾配方向が一致した場合の Restart 処理

### 💡 特徴・工夫点

- **凸性の保証**：効用関数と制約を工夫し、全体を凸計画問題として定式化  
- **Backtracking による収束の安定性**  
- **射影操作により負の値を回避**  
- **モメンタムと Restart を活用した高速収束**

---

