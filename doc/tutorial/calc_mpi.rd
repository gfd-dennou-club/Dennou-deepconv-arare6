=begin JA

= MPI を用いた並列計算を行うには

# * 杉山耕一朗 (sugiyama)
#   * $Id: calc_mpi.rd,v 1.2 2014/03/04 14:31:21 sugiyama Exp $


deepconv/arare5 は MPI (Message Passing Interface) を用いて並列化されている. 
本文書では MPI を用いた並列計算の方法を示す. 


== 並列化の概要

deepconv/arare5 では, 水平方向 (X 方向 & Y 方向) に領域分割することで並
列化を行っている. X 方向と Y 方向は等間隔に分割する. 鉛直方向 (Z 方向)
には分割しない. 

== コンパイル

deepconv/arare5 の並列計算を行うためには, 以下のライブラリが必要です. 

* MPI ライブラリ
* MPI コンパイラでコンパイルした gtool ライブラリ
  * ((<gtool5 インストールガイド|URL:http://www.gfd-dennou.org/library/gtool/gtool5/gtool5_current/INSTALL.htm>)) の「MPI 用にビルドする場合には」を参照下さい. 

deepconv/arare5 をコンパイルする際は, (({Config.mk})) において, 

    FC = mpifrt          # MPI Fortran90 コンパイラ
    CPPFLAGS = LIB_MPI   # MPI 

を指定する必要があります. 詳細は((<インストールガイド|URL:../../INSTALL.htm>))を参照下さい. 


== 設定ファイルの書き方

設定ファイル(NAMELIST ファイル) の &gridset_nml に並列数を追加する必要
があります.

X 方向の並列数を 5, Y 方向共に並列数を 10 とする場合, 以下のように
(({xsub, ysub})) を指定します.

    &gridset_nml
      xsub  = 5                  ! X 方向の並列数
      ysub  = 10                 ! Y 方向の並列数
      xdim  = 50                 ! X 方向刻み点数
      ydim  = 50                 ! Y 方向刻み点数
      zdim  = 50                 ! Z 方向刻み点数
    /

水平鉛直 2 次元の計算を行う場合には, X 方向の並列数を指定します.  
Y 方向には並列できません (ysub を記載しなくても良いです; 
デフォルト値 ysub=1 が使われます). 

    &gridset_nml
      xsub  = 5                  ! X 方向の並列数
      xdim  = 50                 ! X 方向刻み点数
      ydim  = 1                  ! Y 方向刻み点数
      zdim  = 50                 ! Z 方向刻み点数
    /


== 並列計算の実行

計算の実行の手順は, 逐次版と同じく, 

* 初期値の準備
* 実験の実行

である. 

初期値の準備のためには, プロセス数 N ( N = xsub x ysub ) の場合, 
以下のように用意する. 

    $ mpiexec -n N ./arare_init-data -N=arare_init.conf

プロセス毎に初期値ファイルが出力される. 上記を実行することで, N 個の
ファイル (init_rank000000.nc, init_rank000001.nc, init_rank000002.nc, ...)
が出力される.

実験は以下のように実行する. 

    $ mpiexec -n N ./arare -N=arare.conf

上記によって出力される変数もプロセス毎に出力される
(VelX_rank000000.nc, VelX_rank000001.nc, VelX_rank000002.nc, ...,
PTemp_rank000000.nc, PTemp_rank000001.nc, PTemp_rank000002.nc, ...). 


== 出力データの統合

=== データの統合 

並列計算を実行した場合, 計算に使用したプロセス毎にファイルが出力されます.
例えば 128 並列した場合, 1 つの変数に対して 128 個のファイルが生成され
ます. このままでは解析を行うのが面倒ですので, 各変数を 1 つのファイルに
まとめるためのスクリプトを用意しています.

arare_unite.rb は, 設定ファイルより MPI の並列数とファイル名の接頭詞を
読み込むことで, 変数毎に結合したファイルを作成します. 但し, netCDF の
2 GB の壁問題を念頭に, 1 ファイルの大きさが 2 GB を超えないよう調整して
います.

   USAGE:
     $ ruby script/arare_unite.rb <設定ファイル>

   想定するディレクトリ構造
     ./arare.conf         設定ファイル
     ./data/              生データ置き場
            VelX_rank000000.nc
            VelX_rank000001.nc
            VelX_rank000002.nc
            ....
     ./TIME_XXXX-XXXXX/   結合したファイル置き場. (スクリプトにより自動生成される)
            VelX.nc

#=== 平均
#
#時間方向にファイルが分割されていると, 時系列を表示することが面倒です.
#そのため, 予め X 方向と Y 方向に平均した結果を netCDF ファイルとして保
#管するためのスクリプトを用意しています.
#
#arare_meanXY は, カレントディレクトリ内の TIME_XXXX ディレクトリの中の
#ファイルをオープンして, XY 平均・結合して, MEAN_XY ディレクトリ内にファ
#イル出力します. 速度 (VelX, VelY, VelZ) については, 自乗平均 (RMS) も同
#時に計算して出力します.
#
#   USAGE:
#     $ ruby script/arare_meanXY.rb
#
#   想定するディレクトリ構造
#     ./TIME_0000-1000/   結合したファイル置き場.
#            VelX.nc
#     ./TIME_1100-2000/   結合したファイル置き場.
#            VelX.nc
#     ./MEAN_XY/          XY 平均して時間方向に結合したファイル置き場 (スクリプトにより自動生成される)
#            VelX.nc
#            VelXRMS.nc
#

=end JA


=begin HTML
<hr />
<small>
  $Id: calc_mpi.rd,v 1.2 2014/03/04 14:31:21 sugiyama Exp $
</small>
=end HTML

