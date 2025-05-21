# variant filtering

filter_annovar_result.R：
Annovarと他コーラー（例：XHMM）のバリアント表を統合しつつ、アレル頻度、遺伝形式、対象サンプルなどの条件でフィルターする。

options:
  --annovar_result ANNOVAR_RESULT
                        Annovar結果のテーブルファイル
  --out OUT             出力ファイル名
  --af_threshold_hgvd AF_THRESHOLD_HGVD
                        HGVDのアレル頻度の上限値。指定しない場合はフィルターしない
  --af_threshold_tommo AF_THRESHOLD_TOMMO
                        ToMMoのアレル頻度の上限値。指定しない場合はフィルターしない
  --af_threshold_exac_all AF_THRESHOLD_EXAC_ALL
                        ExAC ALLのアレル頻度の上限値。指定しない場合はフィルターしない
  --af_threshold_exac_eas AF_THRESHOLD_EXAC_EAS
                        ExAC EASのアレル頻度の上限値。指定しない場合はフィルターしない
  --inheritance INHERITANCE
                        出力したい遺伝形式（AD,AR,XLのどれか）を指定する
                        AD,ARなら常染色体のみ、XLならX染色体のバリアントに絞る
                        この指定に基づき、--sample_filter and/or
                        --gene_mode_of_inheritance_filterが実行される
  --sample_filter SAMPLE_FILTER
                        指定したサンプルが保持するバリアントのみを出力
                        --inheritanceでARを指定した場合、そのサンプルで2つ以上のバリアントがある遺伝子に絞り込まれる
  --gene_mode_of_inheritance_filter
                        指定された遺伝形式（--inheritance）と一致する遺伝子のバリアントのみを残す GenCC や
                        G2P 由来の MOI (mode of inheritance) 情報を使用
                        MOI列の名前は内部で想定されている（詳細は
                        filter_by_gene_mode_of_inheritance 関数を参照）
                        トリオ解析時など新規の遺伝形式を探索する時は、このフィルターを外す
  --gene_annotations GENE_ANNOTATIONS [GENE_ANNOTATIONS ...]
                        GenCC、G2P、候補遺伝子リストなどの遺伝子に対するアノテーションファイルを指定
                        複数ファイルの指定はカンマかスペース区切り 各ファイルは 'Gene.refGene'
                        という列に遺伝子名を記載する
  --other_caller_results OTHER_CALLER_RESULTS [OTHER_CALLER_RESULTS ...]
                        別のコーラー（例：XHMM）の出力ファイルを指定する、カンマかスペース区切りで複数指定可
                        事前に整形が必要で、XHMMの場合はpreprocess_XHMM.pyを使う
                        その他のコーラーの場合、フォーマットはcat_another_caller_variants関数を参照

使用例は、filter_annovar_result_wrapper.shを参照。

オプションに与えるファイルを用意するためのヘルパースクリプト：
--gene_annotations用:
- preprocess_PanelApp.py
- preprocess_GenCC.py
- preprocess_G2P.py

--other_caller_results用:
- preprocess_XHMM.py

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
🧪 実行例
bash
コピーする
Rscript filter_annovar_result.R \
  --annovar_result example.tsv \
  --out filtered.tsv \
  --inheritance AR \
  --sample_filter Proband001 \
  --af_threshold_exac_all 0.01 \
  --gene_annotations gencc.tsv g2p.tsv \
  --other_caller_results xhmm_results.tsv
また、実行用のラッパースクリプト filter_annovar_result_wrapper.sh も参考にしてください。

🔨 補助スクリプト（前処理）
入力に使うファイルを準備するためのスクリプトも同梱されています。

--gene_annotations 用
スクリプト	説明
preprocess_PanelApp.py	PanelAppの遺伝子リストを整形
preprocess_GenCC.py	GenCCアノテーションを整形
preprocess_G2P.py	G2Pアノテーションを整形
--other_caller_results 用
スクリプト	説明
preprocess_XHMM.py	XHMMの出力を整形
（その他）	その他のツールについては cat_another_caller_variants 関数を参照してください
📄 ライセンス
MIT License（LICENSE を参照）

🙏 謝辞
本ツールは Annovar の出力を元に、外部アノテーションやコーラー結果との統合を通じて、臨床・研究向けの高度なバリアント絞り込みを実現することを目的としています。
