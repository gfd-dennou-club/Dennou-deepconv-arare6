=begin TOPLINK
[((<English|URL:mars-plumetest.htm.en>)) |
((<Japanese|URL:mars-plumetest.htm>))]
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
=end TOPLINK


=begin JA

= 地球大気条件下でのサーマル上昇実験

#* 履歴
#  * 2011/03/04(杉山耕一朗) 新規作成


# * 山下 達也 (yamasita), 杉山 耕一朗 (sugiyama), 小高 正嗣 (odakker)
#   * $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $

=end JA
=begin EN

= A Numerical experiment of ascending hot plume

#* History
#  * 2011/03/01 (Tatsuya Yamashita) Initial release


# * Tatsuya YAMASHITA (yamasita), Ko-ichiro SUGIYAMA (sugiyama), Masatsugu ODAKA (odakker)
#   * $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $

=end EN


=begin JA

ここでは, 単独のホットプリュームを上昇させる実験を実行する方法を説明します. 

この計算には以下の物理過程を用いています. 

  * 乱流拡散
    * 1.5 次のクロージャーモデル (Klemp and Wilhelmson, 1978)
  * 凝結・蒸発
    * Kessler のパラメタリゼーション

=end JA
=begin EN

A method to perform an test experiment
of an ascending hot plume is described. 

Following physical processes are used in this experiment.

  * Turbulent mixing
    * 1.5 order closure (Klemp and Wilhelmson, 1978)
  * Condensation and evaporation
    * Kessler parametarization

=end EN


=begin JA
== 概要

本実験はデフォルトの設定で行なうことが出来ます. 
デフォルト設定での計算をする場合はこの「概要」のみを読めば十分です. 
デフォルトとは異なる設定で計算をする場合は, 「詳細」をご覧ください. 

例えばカレントディレクトリがソースツリーディレクトリである場合, 
以下のコマンドを打つと計算が実行されます. 
計算には 10 分程度かかると思います. 


  $ ./bin/arare -N=arare-earth.conf


計算終了後, gpview などで描画を行ないます. 
例えば gpview で描画する場合, 以下のようにコマンドを打ちます. 

  $ gpview --wsn 4 arare-earth_PotTemp.nc@PotTemp --anim t

gpview の使用方法は
((<GPhys コマンドチュートリアル|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>))
などを参照してください.


=end JA
=begin EN
== Overview

If you perform a calculation under default configuration, 
you have only to see this section "Overview". 
If not so, see the section "Details". 


If you finished compling the source for moist convection with major component 
condensation, you can execute the experiment. 
In order to execute the experiment, command as follows. 
It will take ten minutes or so. 

  $ ./bin/arare -N=arare-earth.conf


After finishing the calculation, next step is visualization. 
For example, if you use gpview, enter a command as follows.

  $ gpview --wsn 4 arare-earth_PotTemp.nc@PotTemp --anim t

If you want to know How to use gpview, for example, 
See
((<this page|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>)). 

=end EN

#=begin HTML
#<tr>
#<img src="mars-moistconv-denscloud.png" width="30%">
#<img src="mars-moistconv-velz.png" width="30%">
#<img src="mars-moistconv-pottemp.png" width="30%">
#</tr>
#=end HTML


=begin JA
== 参考文献
=end JA

=begin EN
== References
=end EN

=begin

* Klemp, J. B., Wilhelmson, R. B., 1978: 
  The simulation of three-dimensional convective storm dynamics", 
  J. Atmos. Sci., 35, pp. 1070 -- 1096.


=end

=begin HTML
<hr />
<small>
  $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $
</small>
=end HTML

