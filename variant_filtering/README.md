# バリアントフィルタリングツール

`filter_annovar_result.R`：**Annovar** および他のバリアントコーラー（例：**XHMM**）の結果を統合し、遺伝子に対するアノテーションを付けて、
アレル頻度・遺伝形式などでバリアントと遺伝子をフィルターする。



---



## 🔧 ０．インストール方法

```bash
# condaの仮想環境内で
mamba install -c conda-forge \
    pandas numpy requests \
    r-tidyverse r-argparse r-glue

git clone https://github.com/hamanakakohei/ngs
git clone https://github.com/hamanakakohei/misc
```



---



## 🔧 １．使い方

```bash
# AD用
filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out exome_summary.filtered.AD.txt \
  --af_threshold_tommo 0.0005 \
  --af_threshold_exac_all 0.00005 \
  --af_threshold_exac_eas 0.0005 \
  --inheritance AD \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz cleaned_G2P.tsv.gz cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz

# AR用
filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out exome_summary.filtered.AR.txt \
  --af_threshold_tommo 0.01 \
  --af_threshold_exac_all 0.01 \
  --af_threshold_exac_eas 0.01 \
  --inheritance AR \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz cleaned_G2P.tsv.gz cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz

# XL用
filter_annovar_result.R \
  --annovar_result exome_summary.txt \
  --out exome_summary.filtered.XL.txt \
  --af_threshold_exac_all 0.00005 \
  --af_threshold_exac_eas 0.0005 \
  --inheritance XL \
  --sample_filter Sample_10000 \
  --gene_mode_of_inheritance_filter \
  --gene_annotations cleaned_GenCC.tsv.gz cleaned_G2P.tsv.gz cleaned_PanelApp.tsv.gz \
  --other_caller_results cleaned_xhmm_Sample_10000.tsv.gz
```

### １－１．入力と出力

| オプション | 説明 |
|------------|------|
| `--annovar_result` | Annovarによるアノテーション結果（TSV） |
| `--out` | 出力ファイル名（TSV） |

### １－２．アレル頻度でフィルター

| オプション | 説明 |
|------------|------|
| `--af_threshold_hgvd` | HGVDのアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_tommo` | ToMMoのアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_exac_all` | ExAC（全集団）のアレル頻度の上限値。指定しなければフィルターしない。 |
| `--af_threshold_exac_eas` | ExAC（東アジア集団）のアレル頻度の上限値。指定しなければフィルターしない。 |

### １－３．遺伝形式でフィルター

| オプション | 説明 |
|------------|------|
| `--inheritance` | 遺伝形式を指定：`AD`、`AR`、`XL` |
| `--sample_filter` | 指定したサンプルで見られるバリアントのみ出力。`--inheritance`でARを指定時は、当該サンプルに2つ以上のバリアントがある遺伝子のみ抽出される。|
| `--gene_mode_of_inheritance_filter` | GenCCやG2Pに基づき、`--inheritance`で指定された遺伝形式で疾患を伝える遺伝子のバリアントのみを出力。トリオ解析など探索的解析ではこのオプションは外す。 |

### １－４．バリアント機能でフィルター

```bash
# 以下のパターンに該当するバリアントは除かれる
`ExonicFunc.refGene == "synonymous SNV" & SpliceAI_max_score < 0.1`
`ExonicFunc.refGene == "." & !str_detect(GeneDetail.refGene, "UTR") &　SpliceAI_max_score < 0.1`
```

### １－５．遺伝子に対するアノテーションファイルの統合

| オプション | 説明 |
|------------|------|
| `--gene_annotations` | 遺伝子に関するアノテーションファイルを、スペース区切りで複数指定可（GenCC, G2P, PanelAppなど）。各ファイルには 遺伝子名を入れた`Gene.refGene` 列が必要。 |

### １－６．他のバリアントコーラー結果の統合

| オプション | 説明 |
|------------|------|
| `--other_caller_results` | XHMMなど他のコーラーの出力ファイル。スペース区切りで複数指定可。事前に整形が必要（例：XHMMは `preprocess_XHMM.py` を使用）。 |



---



## 🔨 ２．前処理用の補助スクリプト

オプションに与えるファイルを準備するためのスクリプトも含まれています。

### ２－１．`--gene_annotations` 用

| スクリプト | 説明 |
|------------|------|
| `preprocess_PanelApp.py` | PanelAppの遺伝子リストを整形 |
| `preprocess_GenCC.py` | GenCCアノテーションを整形 |
| `preprocess_G2P.py` | G2Pアノテーションを整形 |

処理済みのファイルはngs/variant_filtering/result以下にあります：
- cleaned_G2P_2025-04-28.tsv.gz
- cleaned_GenCC_2025_05_20.tsv.gz
- cleaned_PanelApp_ataxia_2025_05_20.tsv.gz

### ２－２．`--other_caller_results` 用

| スクリプト | 説明 |
|------------|------|
| `preprocess_XHMM.py` | XHMMの出力を整形 |
| （その他） | 他のコーラー結果を入力したい場合は、列名の要件など`misc/utils/annovar.R` 内の`cat_another_caller_variants` 関数を参照してください |

#### ２－２－１．使い方

```bash
preprocess_XHMM.py \
  --xhmm_result data.segdup.strvar.haplo.deciph.omim.xcnv.gene \
  --sample Sample_10000 \
  --depth_threshold_for_homo_del 10 \
  --out cleaned_XHMM_Sample_10000.tsv.gz
```

| 引数 | 説明 |
|------|------|
| `--xhmm_result` | XHMM結果ファイル（遺伝子情報付き） |
| `--sample` | 対象とするサンプル名（XHMM結果ファイルのSAMPLE列からこの指定したサンプルのみを取り出す） |
| `--depth_threshold_for_homo_del` | ホモ欠失（1/1）とみなすMEAN_ORIG_RDの閾値（デフォルト: `10`） |
| `--out` | 出力ファイル名 |
