!= Module DataSet
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: dataset.f90,v 1.2 2011/06/17 19:04:00 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
! モデルで利用する物理・化学的な定数を決めるための変数参照型モジュール
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!== Future Plans
!

module dataset
  !
  ! モデルで利用する物理・化学的な定数を決めるための変数参照型モジュール
  !
  
  !モジュール読み込み
  
  !暗黙の型宣言禁止
  implicit none
  
  !save 属性
  save
  
  !公開変数
  
contains
  
  subroutine dataset_init(namelist_filename)
    !
    !=概要
    !
    !NameList ファイルから情報を取得する.
    !このサブルーチン内で, 化学情報の初期化を行っている
    !
    !=凝縮成分の取り扱いについて
    !
    !計算に利用する凝縮成分は, NAMELIST ファイル内の SpcWetSymbol に
    !書かれている. 書き方は以下の通り. 
    !  
    !  H2O-g, H2O-l-Cloud, H2O-l-Rain, NH3-g, H2S-g, NH4SH-s-Cloud
    !
    !ポイントは, 液相と固相には, 雲 or 雨を表す文字列を付加することである.
    !
    !それぞれの化学種には, ChemData.f90 で定義した ID 番号(ChemData_SpcID)が
    !振られている. その ID を, SpcWetID の第 1 次元に保管する. 
    !気相と凝縮相の対応は, ChemData.f90 で定義した ChemData_SameSpc を
    !利用する. 上記の例では, SpcWetID は以下のようになる. 
    !
    ! Symbol:   H2O-g, H2O-l-Cloud, H2O-l-Rain, NH3-g, H2S-g, NH4SH-s-Cloud
    ! ID(:,1):      5,           6,          6,     8,     10,           11
    ! ID(:,2):      0,           5,          5,     0,      0,            8
    ! ID(:,3):      0,           0,          0,     0,      0,           10
    !
    !利用しない部分にはゼロを代入しておく. 
    !
    
    
    !暗黙の型宣言禁止
    implicit none

    

  end subroutine dataset_init
  
end module dataset
