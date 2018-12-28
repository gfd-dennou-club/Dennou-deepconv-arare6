!= Module TimeFilter
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: timefilter.f90,v 1.5 2007-04-19 14:29:33 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview
!
!アッセリンの時間フィルタ. leap-frog の後に適用. 
! 
!== Error Handling
!
!== Bugs
!
!== Note
!
!== Future Plans
!
!

module TimeFilter

  !モジュール読み込み
  use GridSet,  only: DimXMin,     & !x 方向の配列の下限
    &                 DimXMax,     & !x 方向の配列の上限
    &                 DimZMin,     & !z 方向の配列の下限
    &                 DimZMax,     & !z 方向の配列の上限
    &                 SpcNum         !化学種の数

  !暗黙の型宣言禁止
  implicit none

  !属性の設定
  private

  !関数を public 属性に設定
  public AsselinFilter_xz
  public AsselinFilter_pz
  public AsselinFilter_xr
  public AsselinFilter_xza

  interface AsselinFilter_xza
    module procedure AsselinFilter_aaa
  end interface 
  interface AsselinFilter_xz
    module procedure AsselinFilter_aa
  end interface 
  interface AsselinFilter_pz
    module procedure AsselinFilter_aa
  end interface
  interface AsselinFilter_xr
    module procedure AsselinFilter_aa
  end interface

  !変数定義
  real(8) :: tfil = 1.0d-1  !アッセリンの時間フィルタの係数

  !値を save する
  save tfil
  
contains
  
  subroutine AsselinFilter_aa(aa_VarA, aa_VarN, aa_VarB)
    !
    ! 時間フィルター; Asselin のタイムフィルターを利用
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in)     :: aa_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8), intent(inout)  :: aa_VarN(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8), intent(in)     :: aa_VarB(DimXMin:DimXMax, DimZMin:DimZMax)  
    real(8)                 :: aa_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !時間フィルタ
    aa_Var  = aa_VarN + tfil * ( aa_VarB  - 2.0d0 * aa_VarN + aa_VarA ) 
    aa_VarN = aa_Var
    
  end subroutine AsselinFilter_aa
  

  subroutine AsselinFilter_aaa( aaa_VarA, aaa_VarN, aaa_VarB )
    !
    ! 時間フィルター; Asselin のタイムフィルターを利用
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in)     :: aaa_VarA(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8), intent(inout)  :: aaa_VarN(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8), intent(in)     :: aaa_VarB(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8)                 :: aaa_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)

    !時間フィルタ
    aaa_Var  = aaa_VarN + tfil * ( aaa_VarB  - 2.0d0 * aaa_VarN + aaa_VarA ) 
    aaa_VarN = aaa_Var
    
  end subroutine AsselinFilter_aaa
  
end module TimeFilter
