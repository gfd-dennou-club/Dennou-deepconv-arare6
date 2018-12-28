=begin JA

= 凝結物の扱いに関する解説

# * 杉山耕一朗 (sugiyama)
#   * $Id: condense.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $


deepconv/arare5 では複数凝結物を扱えるよう設計しているが, その扱いは少々面倒です. 
本文書では凝結物の内部的扱いや, 凝結物を追加する方法を説明します. 

== 凝結物の設定

設定ファイル (Namelist) に凝結成分の名前とモル比を指定します. 

    &composition_nml
      SpcWetSymbol(1) = 'H2O-g'        , !湿潤成分
      SpcWetSymbol(2) = 'NH3-g'        , !湿潤成分
      SpcWetSymbol(3) = 'H2S-g'        , !湿潤成分
      SpcWetSymbol(4) = 'H2O-s-Cloud'  , !湿潤成分
      SpcWetSymbol(5) = 'H2O-s-Rain'   , !湿潤成分
      SpcWetSymbol(6) = 'NH4SH-s-Cloud', !湿潤成分
      SpcWetSymbol(7) = 'NH4SH-s-Rain' , !湿潤成分
      SpcWetMolFr(1)  = 1.0d-3 ,         !湿潤成分の存在度
      SpcWetMolFr(2)  = 4.0d-4 ,         !湿潤成分の存在度
      SpcWetMolFr(3)  = 6.0d-5 ,         !湿潤成分の存在度
    /

指定可能な凝結物名は src/chem/chemdata.f90 で定義されている. 
雲と雨は "-Cloud" と "-Rain" によって区別する. 


== 凝結物に関する変数・配列

src/setup/composition.f90 で, 凝結物関連の配列の初期化を行っている. 

1) src/chem/chemdata.f90 を参照し, それぞれの凝結物の ID を調べる. 

     SpcSymbol:         SpcID: 
       H2O-g,               5
       NH3-g,               8
       H2S-g,              10
       H2O-l-Cloud,         7
       H2O-l-Rain,          7
       NH4SH-s-Cloud,      11
       NH4SH-s-Rain        11

2) 各カテゴリーに含まれる物質の数を決める. do ループを回す回数に利用する. 

      気体:        GasNum = 3
      凝結物(雲):  CloudNum = 2
      凝結物(雨):  RainNum = 2

3) 各カテゴリーの配列添え字を決める. 気体だけに操作したい場合等々で利用する.
    
      IdxG = 1, 2, 3, 0, 0, 0, ...
      IdxC = 4, 6, 0, 0, 0, 0, ...
      IdxR = 5, 7, 0, 0, 0, 0, ...
    
4) 凝結(Condensation)を生じる物質の数と, それらの配列添え字を決める. 上記の例では H2O の凝結のみが生じる. CondNum は do ループを回す回数として利用する. 
    
      CondNum = 1
      IdxCG = 1, 0, 0, 0, 0, 0, ...
      IdxCC = 4, 0, 0, 0, 0, 0, ...
      IdxCR = 5, 0, 0, 0, 0, 0, ...

飽和蒸気圧を計算する場合には, 以下のように指定している. 

    do i = 1, CondNum
      n = IdxCC(i)                              ! n => 4 ( i => 1 )
      call xyz_SvapPress( SpcID(n), xyz_Temp )  ! SpcID => 7 
    end
    
5) NH4SH の生成反応に関与する物質の配列添え字
    
      IdxNH3    = 2
      IdxH2S    = 3
      IdxNH4SHc = 6
      IdxNH4SHr = 7


== 凝結物のデータの追加方法

src/chemdata/chemdata.f90 にデータを追加します. 

1) 凝結物名を ChemData_SpcSymbol に追加します. 命名則は, "化学種名-相" です. 配列添字が物質に対する ID となります. 

    ChemData_SpcSymbol = (/ &
      & "N2-g   ",  &              !#01
      & "H2-g   ",  &              !#02
      & "He-g   ",  &              !#03
      & "He-g   ",  &              !#04
      & "H2O-g  ",  &              !#05
      & "H2O-l  ",  &              !#06
      & "H2O-s  ",  &              !#07
      & "NH3-g  ",  &              !#08
      & "NH3-s  ",  &              !#09
      & "H2S-g  ",  &              !#10
      & "NH4SH-s",  &              !#11
      & "CO2-g  ",  &              !#12
      & "CO2-s  ",  &              !#13
      & "CH4-g  ",  &              !#14
      & "CH4-l  ",  &              !#15
      & "CH4-s  "   &              !#16
      & /)

2-1) 気相に対して必ず与える必要のある物性値

   * 相 ('Gas' => 気相, 'Liq' => 液相, 'Sol' => 固相)
   * 分子量
   * 分子を構成する元素とその個数
   * [option] 基準状態のエントロピー
   * [option] 基準状態のエンタルピー
   * [option] 比熱の観測値 (温度依存性)

エントロピー, エンタルピー, 比熱の温度依存性は deepconv/arare5 の数値計
算では使っていないが, 物性値のチェックのために与えている. 

2-2) 液相・固相に対して必ず与える必要のある物性値

   * 相 ('Gas' => 気相, 'Liq' => 液相, 'Sol' => 固相)
   * 分子量
   * 分子を構成する元素とその個数
   * 飽和蒸気圧の式 (Antoine の式) の係数
   * [option] 飽和蒸気圧の観測値 (温度依存性)

飽和蒸気圧の観測値は, 与えた飽和蒸気圧の式のチェックのために用いている. 


=end JA


=begin HTML
<hr />
<small>
  $Id: condense.rd,v 1.1 2014/03/01 20:11:22 sugiyama Exp $
</small>
=end HTML

