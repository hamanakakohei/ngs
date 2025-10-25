最新のEHdnのスクリプトだとoutlier解析の出力が変わっている、、、
最新のでEHdn profileして古いのでEHdn merge, outlier.py locus, outlier.py motifをしたら、期待通りの出力になる
top_case_zscore列は、high_case_counts列やcounts列と全く別物なので注意。
基本的に後者の2列の値を使う。
