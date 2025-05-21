# variant filtering

Annovarが出力するバリアント表を


Annovarと他コーラーのバリアント表を統合しつつ、アレル頻度、遺伝形式、対象サンプルなどの条件でフィルターする。

options:
  -h, --help            show this help message and exit
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




アレル頻度でフィルターする


- GenCC
- OMIM
- G2P

疾患の候補遺伝子リスト

他のバリアントコーラー（例：XHMM）の結果を整形してAnnovarの結果に加える
