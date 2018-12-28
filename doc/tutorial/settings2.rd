=begin JA

= 設定ファイルを用いた実験設定の変更

# * 杉山耕一朗 (sugiyama)
#   * $Id: settings2.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $


本文書では設定ファイル(NAMELIST ファイル) を用いた実験設定の変更方法に
ついて記す. 

設定ファイルを変更した後の実際の計算実行の方法については「ごくらくdcpam5」
http://www.gfd-dennou.org/library/dcpam/dcpam5/dcpam5_latest/doc/gokuraku/
を参照されたい.


== 時間刻み・積分期間を変更するには

時間刻みや積分期間は, 設定ファイル(NAMELIST ファイル) に &timeset nml
を用いて設定する. deepconv では時間に関する情報は(({秒単位}))で設定する. 

例えば, 設定ファイルに以下のように設定されている場合, 長い時間ステップ
は 0.4 秒, 短い時間ステップは 0.01 秒, 積分時間は 400 秒である. 

    &timeset_nml
      DelTimeLong  = 0.4d0       !長いタイムステップ
      DelTimeShort = 1.0d-2      !短いタイムステップ(音波関連項)
      IntegPeriod  = 400.0d0     !積分時間 
      ... 
    /

これを積分時間 1000 秒にするためには以下のように変更する. 

    &timeset_nml
      DelTimeLong  = 0.4d0       !長いタイムステップ
      DelTimeShort = 1.0d-2      !短いタイムステップ(音波関連項)
      IntegPeriod  = 1000.0d0    !積分時間 
      ...
    /

deepconv では時間刻みを 2 つ設定する必要がある. 短い時間ステップは音波
に関連する項を解くために用いるものであり, それ以外の項は長い時間ステッ
プを用いる. そのため, 短い時間ステップは音波に対する CFL 条件を満たす必
要がある. deepconv 実行時の標準出力に CFL 条件に関する出力が含まれてい
るので, それをチェックすると良い.

    $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf
                         :        
    *** MESSAGE [cflcheck] ***  Sound Wave Velocity = 346.9604303613247
    *** MESSAGE [cflcheck] ***  DelTimeShort = 0.10000000000000000E-01
    *** MESSAGE [cflcheck] ***  Courant number for DelTimeSort = 0.2453380731118955
                         :        
    ########## PREDICTION OF CALCULATION ###########
    Start Date             2014-02-27T12:22:28+09:00
    Current Date           2014-02-27T12:23:05+09:00
    Progress     50.00%  [************             ]
    Remaining CPU TIME      0.310000E+02
    Completion Date        2014-02-27T12:23:36+09:00
    *** MESSAGE [cflcheck] ***  Courant number of VelX for DelTimeLong = 0.93223408165580368
    *** MESSAGE [cflcheck] ***  Courant number of VelY for DelTimeLong = 0.
    *** MESSAGE [cflcheck] ***  Courant number of VelZ for DelTimeLong = 0.44468541992527427
                         :        
                         :        




=end JA


=begin HTML
<hr />
<small>
  $Id: settings2.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $
</small>
=end HTML

