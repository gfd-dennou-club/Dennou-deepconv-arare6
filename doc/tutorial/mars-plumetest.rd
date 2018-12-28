=begin TOPLINK
[((<English|URL:mars-plumetest.htm.en>)) |
((<Japanese|URL:mars-plumetest.htm>))]
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
=end TOPLINK


=begin JA

= 火星サーマル上昇実験

#* 履歴
#  * 2011/03/01(山下達也) 新規作成


# * 山下 達也 (yamasita), 杉山 耕一朗 (sugiyama), 小高 正嗣 (odakker)
#   * $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $

=end JA
=begin EN

= A Numerical experiment of ascending hot plume

#* History
#  * 2011/03/01 (Tatsuya Yamashita) Initial release


# * Tatsuya YAMASHITA (yamasita), Ko-ichiro SUGIYAMA (sugiyama), Masatsugu ODAKA (odakker)
#   * $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $

=end EN


=begin JA

ここでは, 小高他(2006)で行なわれた, 
単独のホットプリュームを上昇させる実験を実行する方法を説明します. 

この計算には以下の物理過程を用いています. 

  * 乱流拡散
    * 1.5 次のクロージャーモデル (Klemp and Wilhelmson, 1978)
  * 凝結・蒸発
    * 拡散成長 (e.g., Tobie et al., 2003)

=end JA
=begin EN

A method to perform an test experiment
of an ascending hot plume is described. 

Following physical processes are used in this experiment.

  * Turbulent mixing
    * 1.5 order closure (Klemp and Wilhelmson, 1978)
  * Condensation and evaporation
    * Diffusional growth (e.g., Tobie et al., 2003)

=end EN


=begin JA
== 概要

本実験はデフォルトの設定で行なうことが出来ます. 
デフォルト設定での計算をする場合はこの「概要」のみを読めば十分です. 
デフォルトとは異なる設定で計算をする場合は, 「詳細」をご覧ください. 


火星主成分凝結対流用のソースのコンパイルが完了すれば, 
計算を実行することができます. 
コンパイル時には make ではなく, make mmc と打つことに注意して下さい. 
例えばカレントディレクトリがソースツリーディレクトリである場合, 
以下のコマンドを打つと計算が実行されます. 
計算には 10 分程度かかると思います. 


  $ ./bin/arare -N=arare-mmc.conf


計算終了後, gpview などで描画を行ないます. 
例えば gpview で描画する場合, 以下のようにコマンドを打ちます. 

  $ gpview --wsn 4 --range -4.0:4.0 MarsCond_PotTemp.nc@PotTempDist,t=1000

gpview の使用方法は
((<GPhys コマンドチュートリアル|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>))
などを参照してください.

MarsCond_VelZ.nc の変数 VelZ や 
MarsCond_DensCloud.nc の変数 DensCloud 
MarsCond_PotTemp.nc の変数 PotTempDist 
に関して描画を行なうと以下のような図が得られます. 


=end JA
=begin EN
== Overview

If you perform a calculation under default configuration, 
you have only to see this section "Overview". 
If not so, see the section "Details". 


If you finished compling the source for moist convection with major component 
condensation, you can execute the experiment. 
Note that command "make mmc" when compiling.
In order to execute the experiment, command as follows. 
It will take ten minutes or so. 

  $ ./bin/arare -N=arare-mmc.conf


After finishing the calculation, next step is visualization. 
For example, if you use gpview, enter a command as follows.

  $ gpview --wsn 4 --range -4.0:4.0 MarsCond_PotTemp.nc@PotTempDist,t=1000

If you want to know How to use gpview, for example, 
See
((<this page|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>)). 

After drawing, you can get figures as follows. 

=end
=begin HTML
<tr>
<img src="mars-moistconv-denscloud.png" width="30%">
<img src="mars-moistconv-velz.png" width="30%">
<img src="mars-moistconv-pottemp.png" width="30%">
</tr>
=end HTML
=begin


=end EN

=begin JA

== 詳細

デフォルトとは異なる設定で実験を実行する場合, 
以下の 4 つのステップで計算を行います.

  * 使用するプログラムの選択・コピー
  * 設定パラメータの変更
  * ソースのコンパイル
  * 計算の実行

=== 使用するプログラムの選択・コピー

サーマル上昇実験を行なう為には, 
次に示すように src-mmc 以下にあるプログラムを手動でコピーして下さい. 

  $ cp main/arare_091012_masstest6-8_advlong-filter.f90 main/arare.f90
  $ cp io/historyfileio_mmconv2_masstest7.f90 io/historyfileio.f90
  $ cp physocs/densitycloud_turb-masstest2.f90 physics/densitycloud.f90
  $ cp physics/latentheat_all.f90 physics/latentheat.f90
  $ cp physics/masscondense_threshold-masstest4.f90 physics/masscondense.f90
  $ cp physics/adiation_balance4.f90 physics/radiation.f90
  $ cp setup/storepottemp_mmconv2.f90 setup/storepottemp.f90
  $ cp setup/storedenscloud_2d-5.f90 setup/storedenscloud.f90
  $ cp moist/moistbuoyancy_org.f90 moist/moistbuoyancy.f90 
  $ cp env/basicenv_smooth.f90 env/basicenv.f90
  $ cp numdiffusion_smooth_alpha1e-4.f90 util/numdiffusion.f90


=== 設定パラメータの変更

ソースツリーディレクトリ内の火星主成分凝結対流用の namelist ファイル 
arare-mmc.conf を編集します. 

namelist ファイルの各項目の詳細については
((<こちら|URL:../../doc/tutorial/namelist.htm>)) をご覧ください. 

編集項目は多数存在しますが, 
以下では(編集しても問題が起こらないと期待される)
そのいくつかをピックアップします. 

* 積分時間

  TimeInt      = 1800.0d0

デフォルトでは 1800 sec となっていますが, 
計算機の性能や気分に合わせて長くしたり短くしたりしてみて下さい. 

* データ出力時間間隔

  TimeDisp     = 100.0d0

デフォルトでは 100 sec となっていますが, 
計算機の性能や気分に合わせて長くしたり短くしたりしてみて下さい. 

* サーマルの振幅

  DelMax       = 4.0d0

値を小さくし過ぎると, 雲が出来る前にサーマルがつぶれてしまいます. 
また値を大きくし過ぎると, 物理的に意味のある解が得られなくなります. 
基本場を考慮すると, 10 K を超えないようにするのが良いでしょう. 
また DelMax を負値にすると冷たい気塊が下降する実験を行なうことができます. 

* 初期におけるサーマルの中心位置(鉛直方向)の領域に対する割合

  ZcRate       = 0.0d0

0 (下端)から 1 (上端)の間で指定します. 
デフォルトではモデル下端にサーマルを配置しています. 
例えばちょうど中心に配置する場合は, 0.5d0 とします. 

* 臨界飽和比

  SatRatioCr = 1.0d0

凝結が始まるまでにどれほど過冷却にならなければならないかを示すパラメータです. 
物理的に意味がある計算を行なう為には 1 以上に設定する必要があります. 
ちなみに火星大気では 1.0 -- 1.35 の間の値をとると言われています. 


=== ソースのコンパイル

ソースツリーディレクトリで, make mmc コマンドを実行してください. 
既にコンパイルを実行している場合には, 
make mmc を行なう前に make clean を実行して下さい. 

  $ make clean
  $ make mmc


=== 計算の実行

最後に計算を実行します. 
リスタートデータ, 地表面リスタートデータといくつかのヒストリデータ
ファイルが出力されます. 

実行プログラム bin/arare は引数に何も指定しないと, 作業ディレクトリ内の
arare.conf を namelist ファイルとして読み込みます. 
火星主成分凝結計算を行なう場合, arare-mmc.conf を arare.conf としてコピー
してから, 実行します. 

  $ cp arare-mmc.conf arare.conf
  $ ./bin/arare

または namelist ファイルを指定して実行することもできます. 

  $ ./bin/arare -N=arare-mmc.conf

実行時には以下のオプションを与えることができます.

  -N=(NAMELIST ファイルのパス) または --namelist=(NAMELIST ファイルのパス)
     NAMELIST ファイルを陽に指定する.

  -D または --debug
     デバッグメッセージを出力する.

  -H または --help
     ヘルプメッセージを表示して終了する.

また, 標準出力および標準エラー出力を保存することもできます. 
例えばシェルが bash の場合,

  $ /work/deepconv/bin/arare > arare.log 2> arare_error.log

とします. 
これにより, 標準出力ファイル arare.log, 標準エラー出力ファイル 
arare_error.log が生成されます. 


計算が完了すると, 以下のようなデータが出力されます. 


  MarsCond_restart2.nc     リスタートファイル
  MarsCond_BasicZ.nc       基本場変数
  MarsCond_DensCloud.nc    雲密度
  MarsCond_Exner.nc        無次元圧力関数
  MarsCond_H2O-g.nc        水蒸気混合比(本計算ではゼロが出力されるだけ)
  MarsCond_H2O-s-Cloud.nc  雲水混合比(本計算ではゼロが出力されるだけ)
  MarsCond_H2O-s-Rain.nc   雨水混合比(本計算ではゼロが出力されるだけ)
  MarsCond_Kh.nc           乱流拡散係数(運動量)
  MarsCond_Km.nc           乱流拡散係数(スカラー量)
  MarsCond_PotTemp.nc      温位
  MarsCond_SatRatio.nc     飽和比
  MarsCond_VelX.nc         速度(X 成分)
  MarsCond_VelZ.nc         速度(Z 成分)
  MarsCond_Zprof.nc        診断量


得られた netcdf ファイルに関して, 描画を行なうと, 
「概要」で示したような図が得られます. 


=== リスタートデータからの実行

途中から続けて計算を行ないたい場合には, arare-mmc.conf を編集します. 

例えば MarsCond_restart2.nc を基にリスタート計算を行ないたい場合は, 
以下のように記述します. 

  (前)  InitFile    = ""   !初期値ファイル
  (後)  InitFile    = "MarsCond_restart2.nc"

リスタート計算後に出力されるリスタートファイルの名前を編集するには, 
以下の項目を編集します. 
特に, リスタートを開始するファイルとリスタート計算後に出力されるファイル
の名前が同一である場合, 上書きされてリスタートを開始するファイルが無くな
ってしまうので注意して下さい. 

  ReStartFile = "MarsCond_restart3.nc"

arare-mmc.conf の編集が終わったら, 以下のコマンドで実行を行ってください.
なお, 出力されるデータファイルの名前はデフォルトでは決まっているので, 
以前計算したデータはリネームするか, 別のディレクトリに移しておくことを
お薦めします. 

  $ ./bin/arare -N=arare-mmc.conf


=end JA

=begin EN

== Details

If you perform a calculation under non-default configuration, 
the experiment is performed with the following 4 steps:

  * Choosing and copying the programs
  * Changing the configuration
  * Compling
  * Execution of experiments

=== Choosing and copying the programs

For a Numerical experiment of ascending hot plume, 
you have to copy following programs in src-mmc directory.

  $ cp main/arare_091012_masstest6-8_advlong-filter.f90 main/arare.f90
  $ cp io/historyfileio_mmconv2_masstest7.f90 io/historyfileio.f90
  $ cp physocs/densitycloud_turb-masstest2.f90 physics/densitycloud.f90
  $ cp physics/latentheat_all.f90 physics/latentheat.f90
  $ cp physics/masscondense_threshold-masstest4.f90 physics/masscondense.f90
  $ cp physics/adiation_balance4.f90 physics/radiation.f90
  $ cp setup/storepottemp_mmconv2.f90 setup/storepottemp.f90
  $ cp setup/storedenscloud_2d-5.f90 setup/storedenscloud.f90
  $ cp moist/moistbuoyancy_org.f90 moist/moistbuoyancy.f90 
  $ cp env/basicenv_smooth.f90 env/basicenv.f90
  $ cp numdiffusion_smooth_alpha1e-4.f90 util/numdiffusion.f90


=== Changing the configuration

You may edit a namelist file for moist convection with major component 
condensation, "arare-mmc.conf". 

For details of namelist files, see
((<this page|URL:../../doc/tutorial/namelist.htm>)).

For example, you can edit the following entries easily. 

* Integration time

  TimeInt      = 1800.0d0

* Output time interval

  TimeDisp     = 100.0d0

* Amplitude of potential temperature of a plume

  DelMax       = 4.0d0

* Ratio of the altitude of a plume to vertical domain size

  ZcRate       = 0.0d0

* Critical saturation ratio

  SatRatioCr = 1.0d0


=== Compling

Before "make mmc", you should command "make clean".

  $ make clean
  $ make mmc



=== Execution of experiments

In order to execute the experiment, command as follows. 


If you give no argument, the execution program bin/arare read arare.conf
as a namelist file. 
When you perform experiments for Martian atmospheric convection with
major component condensation, you have to copy arare-mmc.conf as arare.conf, 
and execute bin/arare.

  $ cp arare-mmc.conf arare.conf
  $ ./bin/arare

Alternatively, you can also execute bon/arare by specifying a namelist file.

  $ ./bin/arare -N=arare-mmc.conf

When executing bin/arare, you can give following arguments.

  -N=(Path to namelist file) or --namelist=(Path to namelist file)
     Specifying a namelist file explicitly.

  -D or --debug
     Output the debug message.

  -H or --help
     Output the help message. 


Furthermore, you can save standard output and standard error output. 
If your shell is bash, command as follows. 

  $ ./bin/arare > arare.log 2> arare_error.log

Then you can get a standard output file "arare.log" and
a standard error output file "arare_error.log".


If you have finished calculation, you can obtain following data files. 

  MarsCond_restart2.nc     Restart file
  MarsCond_BasicZ.nc       Basic states
  MarsCond_DensCloud.nc    Cloud of density
  MarsCond_Exner.nc        The Exner function
  MarsCond_H2O-g.nc        Mixing ration for H2O vapor
  MarsCond_H2O-s-Cloud.nc  Mixing ration for H2O cloud
  MarsCond_H2O-s-Rain.nc   Mixing ration for H2O rain
  MarsCond_Kh.nc           Turbulent diffusion coefficient(momentum)
  MarsCond_Km.nc           Turbulent diffusion coefficient(scalar)
  MarsCond_PotTemp.nc      Potential Temperature
  MarsCond_SatRatio.nc     Saturation ratio
  MarsCond_VelX.nc         Horizontal velocity
  MarsCond_VelZ.nc         Vertical velocity
  MarsCond_Zprof.nc        Diagnosis

You can get figures shown above by vizualisation for data files. 

=end EN


=begin JA
== 参考文献
=end JA

=begin EN
== References
=end EN

=begin

* M. Odaka, T. Kitamori, K. Sugiyama, K. Nakajima and Y.-Y. Hayashi, 2006:
  Numerical simulation of Martian atmospheric convection, 
  Proc. of the 20th ISAS Atmospheric Science Symposium, JAXA/ISAS, 
  103 -- 106

* Klemp, J. B., Wilhelmson, R. B., 1978: 
  The simulation of three-dimensional convective storm dynamics", 
  J. Atmos. Sci., 35, pp. 1070 -- 1096.

* Tobie, G., Forget, F., Lott, F., 2003: 
  Numerical simulation of winter polar wave clouds observed by
  Mars Global Surveyor Mars Orbiter Laser Altimeter", 
  Icarus, 35, 33 -- 49.

=end

=begin HTML
<hr />
<small>
  $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $
</small>
=end HTML

