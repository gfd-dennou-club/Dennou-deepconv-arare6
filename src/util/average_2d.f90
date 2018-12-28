!= 平均演算モジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: average.f90,v 1.5 2007-04-11 11:59:58 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module average
  !
  !平均操作を行うための関数を束ねたパッケージ型モジュール. 
  !水平 Arakawa-C, 鉛直 Lorentz グリッドとする. 
  !

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !関数を public にする
  public xz_pz, xr_pr
  public pz_xz, pr_xr
  public xz_xr, pz_pr
  public xr_xz, pr_pz
  public pz_xr
  public xr_pz
  public pr_xz
  public xz_pr

  !--------------------------------------------
  interface xz_pz
    module procedure xa_pa
  end interface xz_pz

  interface xr_pr
    module procedure xa_pa
  end interface xr_pr
  !--------------------------------------------
  interface pz_xz
    module procedure pa_xa
  end interface pz_xz

  interface pr_xr
    module procedure pa_xa
  end interface pr_xr
  !--------------------------------------------  
  interface xz_xr
    module procedure az_ar
  end interface xz_xr

  interface pz_pr
    module procedure az_ar
  end interface pz_pr
  !--------------------------------------------
  interface xr_xz
    module procedure ar_az
  end interface xr_xz

  interface pr_pz
    module procedure ar_az
  end interface pr_pz
  !--------------------------------------------

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xa_pa( pz_var ) 
    !
    ! x 方向: 半整数格子から整数格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: pz_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: xa_pa(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値

    xa_pa(DimXMin+1 : DimXMax, DimZMin : DimZMax) =              &
      &  (                                                       &
      &     pz_var(DimXMin+1 : DimXMax,   DimZMin : DimZMax)     &
      &   + pz_var(DimXMin   : DimXMax-1, DimZMin : DimZMax)     &
      &   ) * 5.0d-1
    xa_pa(DimXMin,:) = xa_pa(DimXMin+1,:) 
    
  end function xa_pa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pa_xa( xz_var ) 
    !
    ! x 方向: 整数格子から半整数格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xz_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: pa_xa(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値
    
    pa_xa(DimXMin : DimXMax-1, DimZMin : DimZMax) =             &
      &  (                                                      &
      &     xz_var(DimXMin+1 : DimXMax   , DimZMin : DimZMax)   &
      &   + xz_var(DimXMin   : DimXMax-1 , DimZMin : DimZMax)   &
      &   ) * 5.0d-1
    pa_xa(DimXMax,:) = pa_xa(DimXMax-1,:)

  end function pa_xa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function az_ar( xr_var )
    !
    ! z 方向: 半整数格子から整数格子への変換
    !
  
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xr_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: az_ar(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値
    
    az_ar(DimXMin : DimXMax, DimZMin+1 : DimZMax) =              &
      &  (                                                       &
      &     xr_var(DimXMin : DimXMax, DimZMin+1 : DimZMax   )    &
      &   + xr_var(DimXMin : DimXMax, DimZMin   : DimZMax-1 )    &
      &   ) * 5.0d-1
    az_ar(:,DimZMin) = az_ar(:,DimZMin+1)

  end function az_ar
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function ar_az( xz_var ) 
    !
    ! z 方向: 整数格子点から半整数格子点への変換
    !

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xz_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: ar_az(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値    

    ar_az(DimXMin : DimXMax, DimZMin : DimZMax-1) =              &
      &  (                                                       &
      &     xz_var( DimXMin : DimXMax, DimZMin+1 :DimZMax   )    &
      &   + xz_var( DimXMin : DimXMax, DimZMin   :DimZMax-1 )    &
      &   ) * 5.0d-1

    ar_az(:,DimZMax) = ar_az(:,DimZMax-1)

  end function ar_az

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pz_xr( xr_var ) 
    !
    ! xz 方向: xr から pz な格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xr_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: pz_xr(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値    

    pz_xr(DimXMin : DimXMax-1, DimZMin+1 : DimZMax) =             &
      &  (                                                        &
      &     xr_var(DimXMin+1 : DimXMax,   DimZMin+1 : DimZMax  )  &
      &   + xr_var(DimXMin+1 : DimXMax,   DimZMin   : DimZMax-1)  &
      &   + xr_var(DimXMin   : DimXMax-1, DimZMin+1 : DimZMax  )  &
      &   + xr_var(DimXMin   : DimXMax-1, DimZMin   : DimZMax-1)  &
      &   ) * 2.5d-1
    pz_xr(DimXMax,:) = 0.0d0
    pz_xr(:,DimZMin) = 0.0d0

  end function pz_xr
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xr_pz( pz_var ) 
    !
    ! xz 方向: pz から xr な格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: pz_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: xr_pz(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値
    
    xr_pz(DimXMin+1 : DimXMax, DimZMin : DimZMax-1 ) =           &
      &  (                                                       &
      &     pz_var(DimXMin+1 :DimXMax,   DimZMin+1 : DimZMax  )  &
      &   + pz_var(DimXMin+1 :DimXMax,   DimZMin   : DimZMax-1)  &
      &   + pz_var(DimXMin   :DimXMax-1, DimZMin+1 : DimZMax  )  &
      &   + pz_var(DimXMin   :DimXMax-1, DimZMin   : DimZMax-1)  &
      &  ) * 2.5d-1
    xr_pz(DimXMin, :) = 0.0d0
    xr_pz(:, DimZMax) = 0.0d0

  end function xr_pz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pr_xz( xz_var ) 
    !
    ! xz 方向: 整数格子から半整数格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xz_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: pr_xz(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値
    
    pr_xz(DimXMin : DimXMax-1, DimZMin : DimZMax-1) =              &
      &  (                                                         &
      &     xz_var(DimXMin+1 : DimXMax,   DimZMin+1 : DimZMax   )  &
      &   + xz_var(DimXMin+1 : DimXMax,   DimZMin   : DimZMax-1 )  &
      &   + xz_var(DimXMin   : DimXMax-1, DimZMin+1 : DimZMax   )  &
      &   + xz_var(DimXMin   : DimXMax-1, DimZMin   : DimZMax-1 )  &
      &  ) * 2.5d-1

    pr_xz(DimXMax, :) = 0.0d0
    pr_xz(:, DimZMax) = 0.0d0

  end function pr_xz
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xz_pr( pr_var )
    !
    ! xz 方向: 半整数格子から整数格子への変換
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! 配列の X 方向の下限
      &                   DimXMax,  &! 配列の X 方向の上限
      &                   DimZMin,  &! 配列の Z 方向の下限
      &                   DimZMax    ! 配列の Z 方向の上限
  
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: pr_var(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均演算の対象となる変数
    real(DP)              :: xz_pr(DimXMin:DimXMax, DimZMin:DimZMax)
                                     !平均値

    xz_pr(DimXMin+1 : DimXMax, DimZMin+1 : DimZMax) =               &
      &  (                                                          & 
      &     pr_var(DimXMin+1 : DimXMax,   DimZMin+1 : DimZMax   )   &
      &   + pr_var(DimXMin+1 : DimXMax,   DimZMin   : DimZMax-1 )   &
      &   + pr_var(DimXMin   : DimXMax-1, DimZMin+1 : DimZMax   )   &
      &   + pr_var(DimXMin   : DimXMax-1, DimZMin   : DimZMax-1 )   &
      &   ) * 2.5d-1
    xz_pr(DimXMin, :) = 0.0d0
    xz_pr(:, DimZMin) = 0.0d0

  end function xz_pr
  
end module average
