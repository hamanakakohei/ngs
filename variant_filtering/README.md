# バリアントフィルタリングツール

`filter_annovar_result.R`：**Annovar** および他のバリアントコーラー（例：**XHMM**）の結果を統合し、遺伝子に対するアノテーションを付けて、
アレル頻度・遺伝形式などを元にバリアントと遺伝子をフィルターする。

---

## 🔧 使い方

```bash
# AD用
Rscript filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out filtered.AD.txt \
  --af_threshold_tommo 0.0005 \
  --af_threshold_exac_all 0.00005 \
  --af_threshold_exac_eas 0.0005 \
  --inheritance AD \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz,cleaned_G2P.tsv.gz,cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz

# AR用
Rscript filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out filtered.AR.txt \
  --af_threshold_tommo 0.01 \
  --af_threshold_exac_all 0.01 \
  --af_threshold_exac_eas 0.01 \
  --inheritance AR \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz,cleaned_G2P.tsv.gz,cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz

# XL用
Rscript filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out filtered.XL.txt \
  --af_threshold_exac_all 0.00005 \
  --af_threshold_exac_eas 0.0005 \
  --inheritance XL \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz,cleaned_G2P.tsv.gz,cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz
```

### 必須引数

| オプション | 説明 |
|------------|------|
| `--annovar_result` | Annovarによるアノテーション結果（TSV） |
| `--out` | 出力ファイル名（TSV） |

### フィルター条件

| オプション | 説明 |
|------------|------|
| `--af_threshold_hgvd` | HGVDのアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_tommo` | ToMMoのアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_exac_all` | ExAC（全集団）のアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_exac_eas` | ExAC（東アジア集団）のアレル頻度の上限値。指定しなければフィルターしない。 |

### 遺伝形式に基づくフィルター

| オプション | 説明 |
|------------|------|
| `--inheritance` | 遺伝形式を指定：`AD`、`AR`、`XL` |
| `--sample_filter` | 指定したサンプルで観測されるバリアントのみ出力。`--inheritance`でARを指定時は、当該サンプルに2つ以上のバリアントがある遺伝子のみ抽出される。|
| `--gene_mode_of_inheritance_filter` | GenCCやG2Pに基づき、`--inheritance`で指定された遺伝形式に一致する遺伝子のバリアントのみを出力。トリオ解析など探索的解析では無効化を推奨。 |

### 遺伝子に対するアノテーションファイルの統合

| オプション | 説明 |
|------------|------|
| `--gene_annotations` | 遺伝子に関するアノテーションファイルを複数指定可（GenCC, G2P, PanelAppなど）。各ファイルには 遺伝子名を入れた`Gene.refGene` 列が必要。 |

### 他のバリアントコーラーの結果の統合

| オプション | 説明 |
|------------|------|
| `--other_caller_results` | XHMMなど他のコーラーの出力ファイル。複数指定可。事前に整形が必要（例：XHMMは `preprocess_XHMM.py` を使用）。 |

---

## 🔨 補助スクリプト（前処理）

入力に使うファイルを準備するためのスクリプトも含まれています。

### `--gene_annotations` 用

| スクリプト | 説明 |
|------------|------|
| `preprocess_PanelApp.py` | PanelAppの遺伝子リストを整形 |
| `preprocess_GenCC.py` | GenCCアノテーションを整形 |
| `preprocess_G2P.py` | G2Pアノテーションを整形 |

### `--other_caller_results` 用

| スクリプト | 説明 |
|------------|------|
| `preprocess_XHMM.py` | XHMMの出力を整形 |
| （その他） | 他のコーラー結果を入力したい場合は、列名の要件など `cat_another_caller_variants` 関数を参照してください |


---

```markdown
## 🔧 セットアップ方法

このスクリプトは、**RとPythonの両方のパッケージ**を使用します。  
以下の手順で、必要な依存関係をまとめてインストールできます。

### ✅ 1. Conda 環境の構築（推奨：mamba 使用）

[mamba](https://github.com/mamba-org/mamba) は `conda` の高速代替です。未導入の場合は先に `mamba` をインストールしてください。

以下のコマンドを実行して、必要なライブラリを一括インストールします：

```bash
mamba install -c conda-forge \
    pandas numpy requests \
    r-tidyverse r-argparse r-glue
```

> 📌 補足：`conda-forge` チャネルを利用することで、安定かつ最新のパッケージを取得できます。  
> 初回のみ、以下のコマンドで設定しておくと便利です：
> ```bash
> conda config --add channels conda-forge
> conda config --set channel_priority strict
> ```

---

### ✅ 2. スクリプトの取得

必要なリポジトリを `git clone` で取得してください。

```bash
git clone https://github.com/hamanakakohei/ngs
git clone https://github.com/hamanakakohei/misc
```

---

### ✅ 3. 実行例や使い方

詳しい使い方や実行例は、`ngs/filter_annovar_result_wrapper.sh` を参照してください。

```

---

## 🔎 補足アドバイス

- Python のみ使用する場合は `r-〜` パッケージは不要です。
- R のみ使用する場合は `pandas numpy requests` を省略できます。
- 複数の環境を管理したい場合は、`environment.yml` を使った方法もおすすめです（必要であれば生成します）。

---

このまま **README.md に貼り付けて使えます**。必要があれば、次に「使い方」セクションも一緒に整備しましょうか？
