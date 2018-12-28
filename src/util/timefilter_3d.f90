!= Module TimeFilter_3D
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: timefilter_3d.f90,v 1.4 2007-08-20 07:33:02 odakker Exp $
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

module TimeFilter_3d

  !モジュール読み込み
  use dc_types, only : DP

  use GridSet_3d, only: DimXMin,     & !x 方向の配列の下限
    &                 DimXMax,     & !x 方向の配列の上限
    &                 DimYMin,     & !y 方向の配列の下限
    &                 DimYMax,     & !y 方向の配列の上限
    &                 DimZMin,     & !z 方向の配列の下限
    &                 DimZMax,     & !z 方向の配列の上限
    &                 SpcNum         !化学種の数

  !暗黙の型宣言禁止
  implicit none

  !属性の設定
  private

  !関数を public 属性に設定
  public AsselinFilter

  !変数定義
  real(DP) :: tfil = 1.0d-1  !アッセリンの時間フィルタの係数

  !値を save する
  save tfil
  
contains
  
  subroutine AsselinFilter(VarA, VarN, VarB)
    !
    ! 時間フィルター; Asselin のタイムフィルターを利用
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)     :: VarA(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
    real(DP), intent(inout)  :: VarN(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)  
    real(DP), intent(in)     :: VarB(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
    real(DP)                 :: Var(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)


    !時間フィルタ
    Var  = VarN + tfil * ( VarB  - 2.0d0 * VarN + VarA ) 
    VarN = Var

    
  end subroutine AsselinFilter
  
end module TimeFilter_3d
