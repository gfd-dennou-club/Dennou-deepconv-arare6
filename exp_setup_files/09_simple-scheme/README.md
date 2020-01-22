
= 拡散方程式

== コンパイル

$ gt5frt diffuse.f90 -o diffuse  -I ../../include -L ../../lib -larare

== 実行

$./diffuse -N=diffuse.conf

== 可視化

* 2D アニメーションの例.
  $ gpview diffuse.nc@zeta,y=0,t=0:1e4:50 --anim t --nocont

* z=500 の断面. 異なる時間の値の重ね合わせ.
  $ gpview --overplot 10 diffuse.nc@zeta,y=0,t=0:1e4:10,z=500 --anim t

== 参考

https://www.cosmo.sci.hokudai.ac.jp/~gfdlab/comptech/resume/200_diveq/2012_0202-ogihara.pdf

-----------------------------------------------------

= 移流方程式

== コンパイル

$ gt5frt advect.f90 -o advect  -I ../../include -L ../../lib -larare

== 実行

$ ./advect -N=advect.conf

== 可視化

* 2D アニメーションの例.
  $ gpview advect.nc@zeta,y=0,t=0:1e4:50 --anim t --nocont

* z=500 の断面. 異なる時間の値の重ね合わせ.
  $ gpview --overplot 10 advect.nc@zeta,y=0,t=0:1e4:10,z=500 --anim t

= 音波

== コンパイル

$ gt5frt soundwave.f90 -o soundwave  -I ../../include -L ../../lib -larare

