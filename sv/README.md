# Smoove スクリプト

このディレクトリには、[smoove](https://github.com/brentp/smoove) を用いた構造変異検出パイプラインのスクリプトを配置しています。

## ディレクトリ構成

- `01_smoove_call.sh`：サンプルごとの構造変異検出
- `02_smoove_merge.sh`：検出結果のマージ
- `03_smoove_genotype.sh`：共通サイトにおけるジェノタイプ推定
- `04_smoove_paste.sh`：ジェノタイプ結果の統合

## 使用方法

各スクリプトには対応する `qsub` スクリプトも用意されています。
実行例や依存環境については、上位ディレクトリの `README.md` も参照してください。

---

この `README.md` は、Git が空のディレクトリを無視しないようにするためにも含められています。
