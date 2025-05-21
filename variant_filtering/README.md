# バリアントフィルタリングツール

このリポジトリには、**Annovar** および他のバリアントコーラー（例：**XHMM**）の結果を統合し、  
**アレル頻度・遺伝形式・対象サンプル**・**遺伝子アノテーション**などを元に柔軟にバリアントをフィルターするRスクリプト `filter_annovar_result.R` が含まれています。

---

## 📌 主な機能

- Annovarと他のバリアントコーラーの出力を統合
- 以下の条件に基づくフィルタリングに対応：
  - アレル頻度（HGVD, ToMMo, ExAC）
  - 遺伝形式（AD, AR, XL）
  - 特定サンプルのバリアントのみ抽出
  - 遺伝子アノテーション情報（GenCC、G2P等）

---

## 🔧 使い方

```bash
Rscript filter_annovar_result.R [オプション]

必須引数
オプション	説明
--annovar_result	Annovarによるアノテーション結果（TSVファイル）
--out	出力ファイル名（TSV形式）
フィルター条件
オプション	説明
--af_threshold_hgvd	HGVDのアレル頻度の上限値。指定しなければフィルターしない（NA扱い）
--af_threshold_tommo	ToMMoのアレル頻度の上限値
--af_threshold_exac_all	ExAC（全集団）のアレル頻度の上限値
--af_threshold_exac_eas	ExAC（東アジア集団）のアレル頻度の上限値
遺伝形式に基づくフィルター
オプション	説明
--inheritance	遺伝形式を指定：AD（常染色体優性）、AR（常染色体劣性）、XL（X連鎖）
--sample_filter	指定したサンプルで観測されるバリアントのみ出力。AR指定時は、当該サンプルに2つ以上のバリアントがある遺伝子のみ抽出される
--gene_mode_of_inheritance_filter	GenCCやG2Pに基づき、指定された遺伝形式に一致する遺伝子のバリアントのみを出力。トリオ解析など探索的解析では無効化を推奨
追加のアノテーションファイル・外部ツール出力の統合
オプション	説明
--gene_annotations	遺伝子に関するアノテーションファイルを複数指定可（GenCC, G2P, PanelAppなど）。各ファイルには Gene.refGene 列が必要
--other_caller_results	XHMMなど他のコーラーの出力ファイル。複数指定可。事前に整形が必要（例：XHMMは preprocess_XHMM.py を使用）
