=begin JA

= 設定ファイルを用いた実験設定の変更

# * 杉山耕一朗 (sugiyama)
#   * $Id: settings5.rd,v 1.2 2014/03/04 04:44:05 sugiyama Exp $


本文書では設定ファイル(NAMELIST ファイル) を用いた実験設定の変更方法に
ついて記す. 

設定ファイルを変更した後の実際の計算実行の方法については「ごくらくdcpam5」
http://www.gfd-dennou.org/library/dcpam/dcpam5/dcpam5_latest/doc/gokuraku/
を参照されたい.


== リスタート計算を行うには

この節では, deepconv/arare5 でのリスタート方法について述べる. 
ここで言うリスタートとは, ある期間積分した後で, その最後の状態から
計算を再開することを指す. 

deepcnov/arare5 のリスタート計算は, 以下の手順により行う.

  * リスタートファイルの指定,
  * 計算再開時刻の指定.
  * 積分時間 / 積分終了時刻の指定.

つまり, 再計算のためには, それ以前の計算においてリスタートファイルを
出力しておく必要がある

=== [前回の計算で行うこと] リスタートファイルの出力のための設定 

リスタートファイルは以下のように指定することで出力される. この例では
200 秒毎にリスタートに必要な情報をリスタートファイルに格納する. 

    !!!
    !!!入出力ファイルに関する設定
    !!!
    &restartfileio_nml      
      ...
      OutputFile = "arare-catseys_restart.nc",  !生成するリスタートファイルの名称
    /

    !!!
    !!!積分時間に関する設定
    !!!
    &timeset_nml
      ... 
      DelTimeOutput = 200.0d0    !出力時間間隔
    /


=== [次回の計算のために行うこと] リスタート計算を行うための設定

リスタート計算を行うためには, リスタートする時刻を決める必要がある. 
例えば, リスタートファイル名が input.nc であり, 

    % ncdump -v time input.nc
      netcdf input {
          ...
          double time(time) ;
                  time:long_name = "time" ;
                  time:units = "sec" ;
          ...
        time = 0, 199, 200, 399, 400 ;
    }

ならば, t = 200 sec もしくは t = 400 sec から計算をリスタートさせること
ができる.  (deepconv は時間積分に 2 level scheme を使っているので, リス
タートするためには時刻が 2 つ必要となる. 上記の例では, t = 199 sec と
t = 200 sec, t = 399 sec と t = 400 sec は組である). 


リスタート開始時刻が t = 400 sec, リスタートファイルが
arare-catseys_init.nc, 積分時間が 1000 sec の場合, 以下のように設定する. 
以下の例で FilePrefix を変更しているのは, 前回の計算の結果を上書きしないためである. 

    !!!入出力ファイルに関する設定
    !!!
    &restartfileio_nml      
      InputFile  = "arare-catseys_init.nc", ! 入力に使うリスタートファイルの名称
      ...
    /

    !!!積分時間に関する設定
    !!!
    &timeset_nml
      ... 
      RestartTime   = 400.0d0
      IntegPeriod   = 1000.0d0     !積分時間 
      ... 
    /
      
    !!! データ出力の全体設定
    !!! 
    &gtool_historyauto_nml
      ... 
      FilePrefix = "venus.02._"
      ...
    /


=end JA


=begin HTML
<hr />
<small>
  $Id: settings5.rd,v 1.2 2014/03/04 04:44:05 sugiyama Exp $
</small>
=end HTML

