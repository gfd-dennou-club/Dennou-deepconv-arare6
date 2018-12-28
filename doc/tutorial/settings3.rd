=begin JA

= 設定ファイルを用いた実験設定の変更

# * 杉山耕一朗 (sugiyama)
#   * $Id: settings3.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $


本文書では設定ファイル(NAMELIST ファイル) を用いた実験設定の変更方法に
ついて記す. 

設定ファイルを変更した後の実際の計算実行の方法については「ごくらくdcpam5」
http://www.gfd-dennou.org/library/dcpam/dcpam5/dcpam5_latest/doc/gokuraku/
を参照されたい.


== 物理定数・惑星に関する定数を変更するには

惑星に関する定数は, 設定ファイル(NAMELIST ファイル) に &constants nmlを
用いて設定する. 比熱と乾燥成分の分子量を決める方法には以下の 3 通りがあ
ります.


===比熱と分子量を陽に与える場合

重力加速度, 温位の基準圧力, 地表面温度・圧力, 乾燥成分の比熱と分子量を
以下のように与えます. この場合, (({定積比熱 (CvDry) や気体定数
(GasRDry) はモデル内部で計算されます})).

    &constants_nml
      Grav         = 9.8d0,        !重力             [m/s]
      PressBasis   = 1000.0d2,     !(温位の)基準圧力 [Pa]
      TempSfc      = 300.0d0,      !地表面温度       [K]
      PressSfc     = 1000.0d2,     !地表面圧力       [Pa]
      CpDry        = 1004.0d0,     !乾燥成分の比熱
      MolWtDry     = 0.029d0,      !乾燥成分の分子量
    /

モデルで実際に使われる値をチェックするためには deepconv の標準出力を
チェックするのが良いでしょう. CvDry や GasRDry の値も確認できます. 

    $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf
                         :        
    *** MESSAGE [constants_init] ***  Grav = 9.800000000000000
    *** MESSAGE [constants_init] ***  PressBasis = 100000.
    *** MESSAGE [constants_init] ***  TempSfc = 300.
    *** MESSAGE [constants_init] ***  PressSfc = 100000.
    *** MESSAGE [constants_init] ***  CpDry    = 1004.
    *** MESSAGE [constants_init] ***  CpDryMol = 29.11600000000000
    *** MESSAGE [constants_init] ***  CvDry    = 717.3103448275862
    *** MESSAGE [constants_init] ***  GasRDry  = 286.6896551724137
    *** MESSAGE [constants_init] ***  MolWtDry = 0.2900000000000000
    *** MESSAGE [constants_init] ***  DayTime  = 86400.
                         :        


===比熱と気体定数を与える場合

重力加速度, 温位の基準圧力, 地表面温度・圧力, 乾燥成分の比熱と気体定数
を以下のように与えます. この場合, (({定積比熱や乾燥成分の分子量はモデル
内部で計算されます})). 

    &constants_nml
      Grav         = 9.8d0,        !重力             [m/s]
      PressBasis   = 1000.0d2,     !(温位の)基準圧力 [Pa]
      TempSfc      = 300.0d0,      !地表面温度       [K]
      PressSfc     = 1000.0d2,     !地表面圧力       [Pa]
      CpDry        = 1004.0d0,     !乾燥成分の比熱
      GasRDry      = 286.7d0,      !乾燥成分の気体定数
    /


===物質名とそのモル比を与える場合

重力加速度, 温位の基準圧力, 地表面温度・圧力, 乾燥成分に相当する物質と
モル比を以下のように与えます. 窒素 80 %, 酸素 20 % としています. (({定
積比熱, 定圧比熱, 乾燥成分の分子量はモデル内部で計算されます})).

    &constants_nml
      Grav         = 9.8d0,        !重力             [m/s]
      PressBasis   = 1000.0d2,     !(温位の)基準圧力 [Pa]
      TempSfc      = 300.0d0,      !地表面温度       [K]
      PressSfc     = 1000.0d2,     !地表面圧力       [Pa]
      SpcDrySymbol(1) = 'N2-g',    !乾燥成分の化学種名
      SpcDrySymbol(2) = 'O2-g',    !乾燥成分の化学種名
      SpcDryMolFr(1)  = 0.8d0,     !乾燥成分の存在度
      SpcDryMolFr(2)  = 0.2d0,     !乾燥成分の存在度
    /





=end JA


=begin HTML
<hr />
<small>
  $Id: settings3.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $
</small>
=end HTML

