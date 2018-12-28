=begin JA

= ケルビンヘルムホルツ不安定の実験 

# * 杉山耕一朗 (sugiyama)
#   * $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $

=end JA
=begin EN

= KelvinHelmholtz Instability Experiment 

# * Ko-ichiro Sugiyama
#   * $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $

=end EN

=begin JA
ケルビンヘルムホルツ不安定の実験を実行する方法を説明します. 
=end JA
=begin EN
A method to perform a KelvinHelmholtz instability experiment is described. 
=end EN


=begin JA
== 概要
本実験は以下の 3 つのステップで行います.

  * 実験ディレクトリの準備
  * 初期値の準備
  * 実験の実行
=end JA
=begin EN
== Overview
This experiment is performed with the following 3 steps:

  * Preparation of directory for experiments
  * Preparation of initial condition
  * Execution of experiments
=end EN


=begin JA

== 実験用ディレクトリ準備

deepconv の数値実験はソースツリー内部では行わず
ソースツリーとは別の外部ディレクトリにて行うことを推奨いたします. 

まず deepconv/arare5 ソースのトップディレクトリ(以下の例では arare5-YYYYMMDD とする)に移動してください. 
以下では deepconv/arare5 ソースディレクトリの隣に ../deepconv-exp/catseye-exp ディレクトリを作成し, そこで実験を行うことにします. 
次のように  ../deepconv-exp/catseye-exp ディレクトリを作成し, 移動してください. 

  $ mkdir -p ../deepconv-exp/catseye-exp
  $ cd ../deepconv-exp/catseye-exp

次に, このディレクトリに実行ファイルと設定ファイル置き場を作成します. 

  $ mkdir bin
  $ mkdir conf

最後に作成したディレクトリに実行ファイルと設定ファイルをコピーします. 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/*catseye*.conf conf

なお, 実行ファイルと設定ファイル (NAMELIST ファイル) があれば, どのディレクトリにおいても計算を行うことができます. 


=end JA
=begin EN

== Preparation of a directory for experiments

Let us move to the top directory of deepconv/arare5 src tree
(assuming arare5-YYYYMMDD in the following example). 
Here, we perform an experiment in ../arare5-exp/catseye-exp directory. 
Please create the directory and enter there as follows: 

  $ mkdir -p ../deepconv-exp/catseye-exp
  $ cd ../deepconv-exp/catseye-exp

Then, pleaase create the directories for executable files and configuration
files as follows:

  $ mkdir bin
  $ mkdir conf

Finally, executable files and configuration files are copied as follows: 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/*catseye*.conf conf

Note that you can perform an experiment in any directory by using executable files and configuration (NAMELIST) files.

=end EN


=begin JA

== 初期値データファイルの作成

arare_init-data と 
((<arare-catseye_init-data.conf|URL:../../exp_setup_files/arare-catseye_init-data.conf>))
を用いて初期値ファイル arare-catseye_init.nc を作成します. 

  $ ./bin/arare_init-data -N=conf/arare-catseye_init-data.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye_init-data.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye_init-data.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named IsoThermal...
   *** MESSAGE [main] ***  Making Initial data (disturb) named tanh...
   *** MESSAGE [main] ***  Making Initial wind data (disturb) named tanh...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = arare-catseys_init.nc


== 実験の実行

実行ファイル arare と NAMELIST ファイル
((<arare-catseye.conf|URL:../../exp_setup_files/arare-catseye.conf>))
を用いて, 以下のように arare を実行してください. 
プログラム終了には数分かかります. 

(なお, クロスコンパイル環境では以下の方法でプログラムを
実行することはできないので注意してください. その場合の実行方法
に関しては, その環境でのプログラム実行マニュアルなどを参照ください. )

  $ ./bin/arare -N=conf/arare-catseye.conf | tee catseye.log

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Exner.nc" is closed
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Km.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.320010E-01
   time-integration       0.610998E+02  (1.02 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.611318E+02  (1.02 minutes)

この場合, 約 1 分の時間積分に要しました. 
計算結果は VelX.nc や PTemp.nc として出力されます. 
また, リスタートファイルが arare-catseye_restart.nc として出力されます. 

=end JA

=begin EN

== Create initial data file

Create initial data file "denscurrent-dry_restart.nc"
using "bin/arare_init-data" and ((<"arare-catseye_init-data.conf"|URL:../../exp_setup_files/arare-catseye_init-data.conf>)).
  
  $ ./bin/arare_init-data -N=conf/arare-catseye_init-data.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye_init-data.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye_init-data.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named IsoThermal...
   *** MESSAGE [main] ***  Making Initial data (disturb) named tanh...
   *** MESSAGE [main] ***  Making Initial wind data (disturb) named tanh...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = arare-catseys_init.nc

== Run the experiment

Using an executable files 'arare' and a NAMELIST file
((<arare-catseye.conf|URL:../../exp_setup_files/arare-catseye.conf>)), 
execute 'arare' as follows. 
This program will be finished in few minutes - tens of minutes. 

  $ ./bin/arare -N=conf/arare-catseye.conf | tee catseye.log

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Exner.nc" is closed
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Km.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.320010E-01
   time-integration       0.610998E+02  (1.02 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.611318E+02  (1.02 minutes)

In this case, about 1 minites integration is performed. 
History data are output to 'VelX.nc' and 'PTemp.nc' etc., 
and a restart data is output to 'arare-catseye_restart.nc'. 

=end EN

=begin JA
== 結果の可視化
((<簡単な解析・可視化|URL:./visualization.htm>)) を参照してください. 

((<"IMG:catseye_img01.png">))

=end JA
=begin EN
== Visualization

Please see ((<First step analysis and visualization|URL:./visualization.htm.en>)).

((<"IMG:catseye_img01.png">))

=end EN


=begin JA
== 参考文献
* 坪木と榊原, CReSS マニュアル第 2 版.
=end JA

=begin EN
== References
* Tsuboki and Sakakibara, CReSS User Manual (in Japanese)
=end EN




=begin HTML
<hr />
<small>
  $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
</small>
=end HTML

