=begin JA

= 重力流の実験 by Straka et al. (1993)

# * 杉山耕一朗 (sugiyama)
#   * $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $

=end JA
=begin EN

= Density Current Experiment by Straka et al. (1993)

# * Yoshiyuki O. Takahashi (yot), Shin-ichi Takehiro (takepiro)
#   * $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $

=end EN

=begin JA
Straka et al. (1993) で行われた, 重力流の実験を実行する方法を説明します. 
=end JA
=begin EN
A method to perform a density current experiment by Straka et al. (1993) is described. 
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
以下では deepconv/arare5 ソースディレクトリの隣に ../deepconv-exp/s93-exp ディレクトリを作成し, そこで実験を行うことにします. 
次のように  ../deepconv-exp/s93-exp ディレクトリを作成し, 移動してください. 

  $ mkdir -p ../deepconv-exp/s93-exp
  $ cd ../deepconv-exp/s93-exp

次に, このディレクトリに実行ファイルと設定ファイル置き場を作成します. 

  $ mkdir bin
  $ mkdir conf

最後に作成したディレクトリに実行ファイルと設定ファイルをコピーします. 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin/
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/arare-DensCurrent-dry*.conf conf/

なお, 実行ファイルと設定ファイル (NAMELIST ファイル) があれば, どのディレクトリにおいても計算を行うことができます. 


=end JA
=begin EN

== Preparation of a directory for experiments

Let us move to the top directory of deepconv/arare5 src tree
(assuming arare5-YYYYMMDD in the following example). 
Here, we perform an experiment in ../arare5-exp/s93-exp directory. 
Please create the directory and enter there as follows: 

  $ mkdir -p ../deepconv-exp/s93-exp
  $ cd ../deepconv-exp/s93-exp

Then, pleaase create the directories for executable files and configuration
files as follows:

  $ mkdir bin
  $ mkdir conf

Finally, executable files and configuration files are copied as follows: 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin/
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/arare-DensCurrent-dry*.conf conf/

Note that you can perform an experiment in any directory by using executable files and configuration (NAMELIST) files.

=end EN


=begin JA

== 初期値データファイルの作成

arare_init-data と 
((<arare-DensCurrent-dry_init.conf|URL:../../exp_setup_files/arare-DensCurrent-dry_init.conf>))
を用いて初期値ファイル denscurrent-dry_restart.nc を作成します. 

  $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry_init.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry_init.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named DRY...
   *** MESSAGE [main] ***  Making Initial data (disturb) named CosXZ...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = denscurrent-dry_restart.nc

== 実験の実行

実行ファイル arare と NAMELIST ファイル
((<arare-DensCurrent-dry.conf|URL:../../exp_setup_files/arare-DensCurrent-dry.conf>))
を用いて, 以下のように arare を実行してください. 
プログラム終了には数分〜数十分かかります. 

(なお, クロスコンパイル環境では以下の方法でプログラムを
実行することはできないので注意してください. その場合の実行方法
に関しては, その環境でのプログラム実行マニュアルなどを参照ください. )

  $ ./bin/arare -N=conf/arare-DensCurrent-dry.conf | tee s93.log

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_ElstEnrgy.nc" is closed
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_PotEnrgy.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.720040E-01
   time-integration       0.234947E+03  (3.92 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.235019E+03  (3.92 minutes)

この場合, 約 4 分の時間積分が行われます. 
計算結果は VelX.nc や PTemp.nc として出力されます. 
また, リスタートファイルが denscurrent-dry_restart2.nc として出力されます. 

=end JA

=begin EN

== Create initial data file

Create initial data file "denscurrent-dry_restart.nc"
using "bin/arare_init-data" and ((<"arare-DensCurrent-dry_init.conf"|URL:../../exp_setup_files/arare-DensCurrent-dry_init.conf>)).
  
  $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry_init.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry_init.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named DRY...
   *** MESSAGE [main] ***  Making Initial data (disturb) named CosXZ...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = denscurrent-dry_restart.nc

== Run the experiment

Using an executable files 'arare' and a NAMELIST file
((<arare-DensCurrent-dry.conf|URL:../../exp_setup_files/arare-DensCurrent-dry.conf>)), 
execute 'arare' as follows. 
This program will be finished in few minutes - tens of minutes. 

  $ ./bin/arare -N=conf/arare-DensCurrent-dry.conf | tee s93.log

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_ElstEnrgy.nc" is closed
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_PotEnrgy.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.720040E-01
   time-integration       0.234947E+03  (3.92 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.235019E+03  (3.92 minutes)

In this case, about 4 minites integration is performed. 
History data are output to 'VelX.nc' and 'PTemp.nc' etc., 
and a restart data is output to 'denscurrent-dry_restart2.nc'. 

=end EN

=begin JA
== 結果の可視化
((<簡単な解析・可視化|URL:./visualization.htm>)) を参照してください. 

((<"IMG:exp-s93_img01.png">))

=end JA
=begin EN
== Visualization

Please see ((<First step analysis and visualization|URL:./visualization.htm.en>)).

((<"IMG:exp-s93_img01.png">))

=end EN


=begin JA
== 参考文献
=end JA

=begin EN
== References
=end EN

=begin
* Straka, J. M., Wilhelmson, R. B., Wicker, A. L., Anderson, J. R., Droegemeier, K. K., 1993 : Numerical solutions of a non-linear density current : A benchmark solution and comparisons. International Journal for Numerical Methods in Fluids , 17, 1--22

=end


=begin HTML
<hr />
<small>
  $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
</small>
=end HTML

