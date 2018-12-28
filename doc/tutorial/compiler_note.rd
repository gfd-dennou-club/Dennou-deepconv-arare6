=begin JA

= いくつかのコンパイラに関する注意書き

# * 森川 靖大, 高橋 芳幸
# * 杉山 耕一朗, dcpam からコピー
#   * $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $

=end JA
=begin EN

= Notes about some compilers

# * Yasuhiro Morikawa, Yoshiyuki O. Takahashi
# * Ko-ichiro Sugiyama, copy from dcpam
#   * $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $

=end EN

=begin HTML
<a name="ifort">
<a name="g95">
=end HTML

=begin JA


== PGI Fortran

=== コンパイルオプション

浮動小数点演算の精度を保証するためのオプションを必ず付けて下さい. 
deepconv においても, コンパイラーの最適化を強くかけた場合に, 
計算結果がおかしくなる例が報告されています. 

   -Kieee


== Intel Fortran

=== コンパイルオプション

浮動小数点演算の精度を保証するためのオプションを必ず付けて下さい. 
deepconv においても, コンパイラーの最適化を強くかけた場合に, 
計算結果がおかしくなる例が報告されています. 

   -fp-model strict -prec-div
   (もしくは, -fp-model precise -prec-div)



== Intel Fortran, G95 Fortran

=== セグメンテーションエラー

OS およびカーネルのバージョンによっては, 解像度を上げた計算を
行った際にセグメンテーションエラーが生じる場合があります.

この場合には, メモリのスタック領域の使用容量の最大値が比較的小さく
設定されている可能性があります. 

多くの場合にはこの使用容量の最大値は
ユーザでも設定可能なソフトリミットと呼ばれるものであるため, 
以下のように最大値を増やすよう設定可能です. 
ただし, システム側で決められている最大値であるハードリミット以上に
大きくソフトリミットを設定することはできません. 


==== スタック領域のソフトリミットをハードリミットまで増やす

: csh, tcsh でスタック領域のソフトリミットの設定値を変更する

  スタック領域のハードリミットを調べます
  
    > limit -h | grep stack
    stacksize    XXXXXXXX

  ここで XXXXXXXX として表示されるのが
  スタック領域のハードリミットです.
  (数値または "unlimited" が表示されます). 
  この XXXXXXXX をスタック領域のソフトリミットに設定します.

    > limit stacksize XXXXXXXX

: sh, bash でスタック領域のソフトリミットの設定値を変更する

  スタック領域のハードリミットを調べます
  
    $ ulimit -Ha |grep stack
    stack size            (kbytes, -s) XXXXXXXX

  ここで XXXXXXXX として表示されるのが
  スタック領域のハードリミットです.
  (数値または "unlimited" が表示されます). 
  この XXXXXXXX をスタック領域のソフトリミットに設定します.

    $ ulimit -s XXXXXXXX

##
## 安藤さんから Mac の場合について情報をもらうこと.
##

# 以下はメモ書き (森川が大阪府大の重氏宛に書いたメールより抜粋)
#
#== プログラムが動かなかった原因
#
#プログラムが動かなかった原因はスタック領域に割り当てられるメモリの不足
#です. 実メモリを潤沢に積んでいても, 各ユーザ毎にメモリの使用容量に限界
#が設定されており, limit や ulimit -a コマンドでその限界が表示できます. 
#どうやら Linux Kernel 2.6 系ではスタックに割り当てられるソフトリミット
#が約 8 MB になっているため, メモリ不足でセグメンテーションエラーが発生
#していたようです.
#
#ちなみに, メモリにはスタック領域以外にヒープ領域もあり, プログラムがス
#タック領域を使用するかどうかはプログラムに書き方 (特に各変数の宣言文を
#どのように書くか) に拠ります. dcpam5 では, 再コンパイルせずに解像度を
#変えるため, また私的な可読性を優先するなどの理由から配列の動的割付に自
#動配列を用いています. この配列に関してはメモリのスタック領域が使用され
#ることになります. 例えば T42L20 の場合, 使用メモリサイズは全体で約 72
#MB です. このうちスタック領域を利用する分がソフトリミットの 8 MB を超
#えるためにエラーが生じています.
#
#なお, 自動配列などに関しては, 気象庁 Fortran 標準コーディングルール
#http://www.mri-jma.go.jp/Project/mrinpd/coderule.html の「3.2 動的割り
#付け」が参考になるかと思います. (dcpam でも基本的にはこのコーディング
#ルールを準拠しています).
#
#


また, Intel Fortran の ver.10.0 以上ではコンパイルオプション (FFLAGS) に

  -heap-arrays

を追加してください. これは次のような事情によります.

Intel Fortran の ver.10.0 以上では一時的なメモリ割り当てを
スタックへの割り付けることが増えました.
このため, 上記のように limit または ulimit により制限を緩めても
セグメンテーションエラーが出ることがあります.
-heap-arrays を指定すると一時的なメモリ割り当てが
ヒープ領域上におこなわれるようになるので, セグメンテーションエラー
が回避できます. 

参考：http://software.intel.com/en-us/articles/intel-fortran-compiler-increased-stack-usage-of-80-or-higher-compilers-causes-segmentation-fault/


=end JA

=begin EN
== Intel Fortran, G95 Fortran

=== Segmentation fault

When a high resolution calculation is performed, a segmentation fault
might be caused according to version of OS and kernel.

In this case, there is a possibility that the maximum value of the use
capacity of a stack area of the memory is set comparatively small. 

In many cases, the maximum value of this use capacity is the one that
is called "soft limit" that users can set. Therefore, it can be set
that the maximum value is increased as follows.
However, a soft limit cannot be greatly set more than a "hard limit"
that is the maximum value that has been decided on the system side.

==== A soft limit of a stack area is increased to a hard limit

: With csh or tcsh, A setting value of a soft limit of a stack area is changed

  A hard limit of a stack area is examined as follows.
  
    > limit -h | grep stack
    stacksize    XXXXXXXX

  "XXXXXXXX" is a value of a hard limit.
  (Numerical value or "unlimited" is displayed).
  This "XXXXXXXX" is set to a soft limit of a stack area as follows. 

    > limit stacksize XXXXXXXX

: With sh or bash, A setting value of a soft limit of a stack area is changed

  A hard limit of a stack area is examined as follows.
  
    $ ulimit -Ha |grep stack
    stack size            (kbytes, -s) XXXXXXXX

  "XXXXXXXX" is a value of a hard limit.
  (Numerical value or "unlimited" is displayed).
  This "XXXXXXXX" is set to a soft limit of a stack area as follows. 

    $ ulimit -s XXXXXXXX


=end EN

=begin HTML
<a name="frt">
=end HTML

=begin JA
== Fujitsu Fortran

=== NAMELIST 読み込みエラー

ファイルに記述される NAMELIST が多い場合に, まれに一部の
NAMELIST が読み込まれない場合があるようです. 

ファイルにはちゃんと書き込まれているはずなのに, NAMELIST の情報が実行
プログラムに反映されない場合には, 実行プログラムによって出力されるメッ
セージをチェックしてください. もしも以下のようなメッセージが表示されて
いれば, ここで述べるような症状が現れていることになります.

 !*** WARNING [XxxxNmlRead] ***  NAMELIST group "xxxx_nml" is not found
    in "xxxx_xxxx.nml" (iostat=190).
                        ^^^^^^^^^^ これが特徴.

この場合には, 読み込まれない NAMELIST の前に空行やコメント行を挿入する
などし, 再度プログラムの実行と動作のチェックを行ってください. 
空行やコメント行を何行か足すことで上記の症状は回避されるようです.
(原因は不明です). 

=end JA

=begin EN
== Fujitsu Fortran

=== NAMELIST loading error

A part of NAMELIST group name is not loaded in rare cases 
when a lot of NAMELIST group names is described in a file. 

Please check messages output by an execution program when information
of NAMELIST is not reflected in the execution program though it is
sure to be written in the file correctly.
The symptom described in the above-mentioned will appear if  
following messages are displayed.

 !*** WARNING [XxxxNmlRead] ***  NAMELIST group "xxxx_nml" is not found
    in "xxxx_xxxx.nml" (iostat=190).
                        ^^^^^^^^^^ This is the feature. 

In this case, please insert some null lines or comment lines in
front of the NAMELIST, and check operation with repeated
execution of the program.
The above-mentioned symptom seems to be evaded by some null lines or
comment lines. (The cause is uncertain).

=end EN


=begin HTML
<hr />
<small>
  $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $
</small>
=end HTML
