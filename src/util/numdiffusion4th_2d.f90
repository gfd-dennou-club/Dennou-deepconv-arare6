!= Module NumDiffusion
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: numdiffusion.f90,v 1.19 2011-02-28 12:28:38 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module NumDiffusion4th
  !
  !数値拡散項の計算モジュール
  ! 4 次精度中心差分を利用
  !

  !モジュール呼び出し
  use dc_types,  only : DP

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !関数を public にする
  public NumDiffusion_Init
  public xz_NumDiffScalar
  public xz_NumDiffKm
  public xza_NumDiffScalar
  public pz_NumDiffVelX
  public xr_NumDiffVelZ

  !変数の定義
  real(DP), save  :: NuH   = 0.0d0   !数値粘性の係数 (水平方向)
  real(DP), save  :: NuV   = 0.0d0   !数値粘性の係数 (鉛直方向)
! real(DP), save  :: Alpha = 1.0d-4
  real(DP), save  :: Alpha = 0.0d0
  real(DP), save  :: Alpha_Velocity = 1.0d0

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine NumDiffusion_init(cfgfile)
    ! 
    ! NumDiffusion モジュールの初期化ルーチン
    !

    !モジュール読み込み
    use dc_message,  only: MessageNotify
    use mpi_wrapper, only: myrank
    use axesset,     only: dx, dy, dz     !格子点間隔
    use timeset,     only: DelTimeLong    !長い時間ステップ

    !暗黙の型宣言禁止
    implicit none

    ! 変数の設定
    character(*), intent(in) :: cfgfile

    ! Namelist から Alpha の値を決め直す.
    NAMELIST /numdiffusion/ Alpha, Alpha_Velocity
    open (10, FILE=cfgfile)
    read(10, NML=numdiffusion)
    close(10)
    
    ! CReSS マニュアルでは, 2 次精度中心差分の場合, Alpha < 1/8 くらいが適当と述べている.
    NuH = Alpha * ( dx ** 4.0d0 ) / ( 2.0d0 * DelTimeLong )
    NuV = Alpha * ( dz ** 4.0d0 ) / ( 2.0d0 * DelTimeLong )
    
    !確認
    if (myrank == 0) then
      call MessageNotify( "M", "NumDiffusion_init", "Alpha = %f", d=(/Alpha/) )
      call MessageNotify( "M", "NumDiffusion_init", "NuH = %f", d=(/NuH/) )
      call MessageNotify( "M", "NumDiffusion_init", "NuV = %f", d=(/NuV/) )
      call MessageNotify( "M", "NumDiffusion_init", "Alpha_Velocity = %f", d=(/Alpha_Velocity/) )
      if (Alpha > 0.125) then
        call MessageNotify( "E", "NumDiffusion_init", "Alpha = %f > 1/8", d=(/Alpha/) )
        stop
      end if
    end if

  end subroutine NumDiffusion_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xz_NumDiffScalar(xz_Scalar)
    !
    ! x, z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: xz_Scalar(imin:imax, kmin:kmax)
                                                    !スカラー量
    real(DP)             :: xz_NumDiffScalar(imin:imax, kmin:kmax)
                                                    !水平方向の数値拡散
    
    xz_NumDiffScalar =   &
      &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xz_Scalar ))))) &
      &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xz_Scalar ))))) 
    
  end function xz_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xz_NumDiffKm(xz_Scalar)
    !
    ! x, z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: xz_Scalar(imin:imax, kmin:kmax)
                                                    !スカラー量
    real(DP)             :: xz_NumDiffKm(imin:imax, kmin:kmax)
                                                    !水平方向の数値拡散
    
    xz_NumDiffKm =   &
      &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xz_Scalar ))))) &
      &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xz_Scalar ))))) 
    
  end function xz_NumDiffKm

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xza_NumDiffScalar(xza_Scalar)
    !
    ! x, z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax, ncmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: xza_Scalar(imin:imax, kmin:kmax, ncmax)
                                                    !スカラー量
    real(DP)             :: xza_NumDiffScalar(imin:imax, kmin:kmax, ncmax)
                                                    !水平方向の数値拡散
    integer             :: s
    

    do s = 1, ncmax
      xza_NumDiffScalar(:,:,s) =   &
        &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xza_Scalar(:,:,s) ))))) &
        &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xza_Scalar(:,:,s) ))))) 
    end do

  end function xza_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function pz_NumDiffVelX(pz_VarX)
    !
    ! z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: pz_dx_xz, xz_dx_pz, &
      &                   pz_dz_pr, pr_dz_pz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pz_VarX(imin:imax, kmin:kmax)
                                                    !物理量
    real(DP)             :: pz_NumDiffVelX(imin:imax, kmin:kmax)
                                                    !数値拡散
    
    !速度の拡散係数だけ変更する場合にAlpha_Velocityの値を1以外の値にする
    pz_NumDiffVelX = &
      & -  Alpha_Velocity * NuH * ( pz_dx_xz( xz_dx_pz( pz_dx_xz( xz_dx_pz( pz_VarX ) ) ) ) ) &
      & - Alpha_Velocity * NuV * ( pz_dz_pr( pr_dz_pz( pz_dz_pr( pr_dz_pz( pz_VarX ) ) ) ) )
    
  end function pz_NumDiffVelX
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xr_NumDiffVelZ(xr_VarZ)
    !
    ! x 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xr_dx_pr, pr_dx_xr, &
      &                   xr_dz_xz, xz_dz_xr

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xr_VarZ(imin:imax, kmin:kmax)
                                                    !物理量
    real(DP)             :: xr_NumDiffVelZ(imin:imax, kmin:kmax)
                                                    !数値拡散

    !速度の拡散係数だけ変更する場合にAlpha_Velocityの値を1以外の値にする
    xr_NumDiffVelZ = &
      & - Alpha_Velocity * NuH * ( xr_dx_pr( pr_dx_xr( xr_dx_pr( pr_dx_xr( xr_VarZ ) ) ) ) )&
      & - Alpha_Velocity * NuV * ( xr_dz_xz( xz_dz_xr( xr_dz_xz( xz_dz_xr( xr_VarZ ) ) ) ) )
    
  end function xr_NumDiffVelZ

end module NumDiffusion4th
