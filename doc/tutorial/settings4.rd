=begin JA

= 設定ファイルを用いた実験設定の変更

# * 杉山耕一朗 (sugiyama)
#   * $Id: settings4.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $


本文書では設定ファイル(NAMELIST ファイル) を用いた実験設定の変更方法に
ついて記す. 

設定ファイルを変更した後の実際の計算実行の方法については「ごくらくdcpam5」
http://www.gfd-dennou.org/library/dcpam/dcpam5/dcpam5_latest/doc/gokuraku/
を参照されたい.


== 出力設定を変更するには

解析用のヒストリデータの出力に関する設定は &gtool_historyauto_nml を編
集することで変更する. 以下では最低限の説明を行います. 詳細は ((<gtool5
チュートリアル集
|URL:http://www.gfd-dennou.org/library/gtool/gt5tutorial/>)) をご覧下さ
い.

例えば, ヒストリーデータの設定として以下が設定されているとします. 

    ! データ出力の全体設定
    !
    &gtool_historyauto_nml
       FilePrefix = "arare-catseye_" ! ファイルの接頭詞
       IntValue = 10.0,              ! 出力間隔の数値
       IntUnit = 'sec',              ! 出力間隔の単位
       AllOutput = .false.
    /
    !
    ! データ出力の個別設定
    !
    &gtool_historyauto_nml
      Name = 'PTemp, VelX, VelZ, Exner, Km',     ! 出力変数
      IntUnit = 'sec',                           ! 出力間隔の単位
    /
    &gtool_historyauto_nml
      Name = 'PTempAdv,PTempTurb',               ! 出力変数
      SpaceAverage = .true.,.true.,.false.
      IntUnit = 'sec',                           ! 出力間隔の単位
    /


それぞれの意味は, 以下の通りです. 出力する変数を増やす場合には, 変数名
を (({Name})) に追加して下さい. また, (({AllOutput})) を (({.true.}))
にすると, 全ての変数が出力されます.


    FilePrex
       (文字型) データのファイル名の接頭詞. 例えば"exp1-" と指定すれば, 変数
       "VelX" の出力ファイル名は"exp1-VelX.nc" となる. また, "data01/" のようにス
       ラッシュを含む文字列を指定することで, カレントディレクトリ以外の場所
       に出力するよう設定することも可能である.

    IntValue
       (実数型) 出力間隔の数値

    IntUnit
       (文字型) 出力間隔の単位. "sec", "min", "hour", "day", "month", "year" な
       どが使用可能である. 使用可能な単位の詳細については, gtool5 ライブラリ: 
       dc date types モジュール (http://www.gfd-dennou.org/library/gtool/gtool5/
       gtool5_current/doc/code_reference/classes/dc_date_types.html) の"Char-
       acters list for unit" を参照されたい.

    AllOutput
       (論理型) HistoryAutoAddVariable によってプログラム内で登録された
       変数を全て出力するためのフラグ. デフォルトではNAMELIST による変
       数ごとの個別出力設定に示すように, 変数は明示的に指定しない限り出
       力されませんが, この項目を ".true." とすることで, 全ての変数が出
       力されます. 

    Precision
       (文字型) データの精度. "oat" (単精度実数型), "double" (倍精度実数型), "int" (整数型) を指定可能

    Name
       (文字型) 出力する変数名. 複数指定する場合はカンマでつなげる. 

    SpaceAverage
       (論理型配列) 空間平均のフラグ. 配列の1 番目, 2 番目, 3 番目が, 経度, 緯
       度, 高度に対応する.




=end JA


=begin HTML
<hr />
<small>
  $Id: settings4.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $
</small>
=end HTML

