#!/usr/bin/env python3


import argparse
import pandas as pd

# 与えるファイルの復習
# outlier_locus.tsvの中身：
# contig  start  end   motif  gene  region    top_case_zscore  high_case_counts       counts
# chr10   1000   1001  CCG    CHD3  intronic  12.0             DA001:12.0,DA002:1.14  DA003:1.14,DA004:1.14
#
# outlier_motif.tsvの中身：
# motif  top_case_zscore  high_case_counts       counts
# CCG    1.44             DA001:2.32,DA002:2.13  DA003:1.13,DA004:1.13
#
# 各列の内容（参照：https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/04_Outlier_quickstart.md）
# - top_case_zscore:    Top z-score of a case sample
# - high_case_counts:   Counts of case samples corresponding to z-score greater than 1.0
# - counts:             Nonzero counts for all samples


def select_one_family(df, proband, others, controls):
    # 指定した患者が含まれる行に限定
    df = df[df['counts'].str.contains(proband)].copy()

    members = [proband] + others
    read_counts_by_member = {member: [] for member in members}
    more_expanded_all = []
    more_expanded_ctrl = []

    for _, row in df.iterrows():
        # 各行の'counts'列（「サンプル：カウント」が「,」でつながっている）を「サンプル：カウント」の辞書にする
        count_dict = {
            sample: float(count)
            for sample, count in (
                item.split(':') for item in row['counts'].split(',') 
            )
        }

        # 家族の各メンバーのカウント（なければ 0.0）を得て辞書内のリストに追加していく
        member_counts = {member: count_dict.get(member, 0.0) for member in members}
        for member in members:
            read_counts_by_member[member].append(member_counts[member])

        # proband のリード数を得て、
        # 全サンプルで患者より大きいカウントを持つ人の数と
        # コントロールで患者より大きいカウントを持つ人の数を調べて、リストに追加する
        proband_count = member_counts[proband]

        more_expanded_all.append(
            sum(c > proband_count for    c in count_dict.values())
        )

	more_expanded_ctrl.append(
    	    sum(c > proband_count for s, c in count_dict.items() if s in controls
	)

    # 作ったリストを基の表に列として追加する
    # サンプル名を列名にして、motif/locus read countの値を入れる
    df['more_expanded_samples_all'] = more_expanded_all
    df['more_expanded_samples_ctrl'] = more_expanded_ctrl

    for member in members:
        df[member] = read_counts_by_member[member]

    return df.drop(columns=['top_case_zscore'])


def rename_motif_proband( motif_proband, members ):
  # motifファイルをlocusファイと結合する時に列名が被るので変える
  motif_proband = motif_proband.\
    rename({
      'high_case_counts':           'high_case_counts__motif',
      'counts':                     'counts__motif',
      'more_expanded_samples_all':  'more_expanded_samples_all__motif',
      'more_expanded_samples_ctrl': 'more_expanded_samples_ctrl__motif'
    }, axis=1)

  # 各サンプルのmotif read countは列名がサンプル名なので、それも変える
  for member in members:
    motif_proband = motif_proband.\
      rename( { member: f'{member}_motif' }, axis=1 )
  
  return motif_proband


def run():
  parser = argparse.ArgumentParser( description='' )
  parser.add_argument( "--proband",                 type=str, required=True, help="発端者のID" ) 
  parser.add_argument( "--other_family_members",    type=str, required=True, help="他の家族のIDリストのファイル" )
  parser.add_argument( "--controls",                type=str, required=True, help="健常者のIDリストのファイル" )
  parser.add_argument( "--outlier_locus",           type=str, required=True, help="outlier locus解析の結果" )
  parser.add_argument( "--outlier_motif",           type=str, required=True, help="outlier motif解析の結果" )
  parser.add_argument( "--out",                     type=str, required=True, help="" )
  args = parser.parse_args()
  proband = args.proband
  
  with open( args.other_family_members, 'r' ) as f:
    others   = [i.strip() for i in f.readlines()]    

  with open( args.controls, 'r' ) as f:
    controls = [i.strip() for i in f.readlines()]    

  locus_df = pd.read_table( args.outlier_locus )
  motif_df = pd.read_table( args.outlier_motif )

  locus_proband = select_one_family( locus_df, proband, others, controls )
  motif_proband = select_one_family( motif_df, proband, others, controls )
  motif_proband = rename_motif_proband( motif_proband, [ proband ] + others )

  pd.merge( locus_proband, motif_proband, how='outer' ).\
    to_csv( args.out, sep='\t', index=False )


if __name__ == "__main__":
  run()
