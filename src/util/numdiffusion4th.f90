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
  public xyz_NumDiffScalar
  public xyz_NumDiffKm
  public xyza_NumDiffScalar
  public pyz_NumDiffVelX
  public xqz_NumDiffVelY
  public xyr_NumDiffVelZ

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
    NuH = Alpha * ( SQRT( dx * dy ) ** 4.0d0 ) / ( 2.0d0 * DelTimeLong )
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

  function xyz_NumDiffScalar(xyz_Scalar)
    !
    ! x, y, z 方向に半格子ずれた点での数値拡散項を評価
    !

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax
    use differentiate_center2,                    &
      &             only: xyz_dx_pyz, pyz_dx_xyz, &
      &                   xyz_dy_xqz, xqz_dy_xyz, &
      &                   xyz_dz_xyr, xyr_dz_xyz
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(8), intent(in) :: xyz_Scalar(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量
    real(8)             :: xyz_NumDiffScalar(imin:imax,jmin:jmax,kmin:kmax)
                                                    !水平方向の数値拡散
    
    xyz_NumDiffScalar =   &
      &  - NuH * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_Scalar ))))) &
      &  - NuH * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_Scalar ))))) &
      &  - NuV * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_Scalar ))))) 
    
  end function xyz_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyz_NumDiffKm(xyz_Scalar)
    !
    ! x, z 方向に半格子ずれた点での数値拡散項を評価
    !

    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax
    use differentiate_center2,                    &
      &             only: xyz_dx_pyz, pyz_dx_xyz, &
      &                   xyz_dy_xqz, xqz_dy_xyz, &
      &                   xyz_dz_xyr, xyr_dz_xyz
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(8), intent(in) :: xyz_Scalar(imin:imax,jmin:jmax,kmin:kmax)
                                                    !スカラー量
    real(8)             :: xyz_NumDiffKm(imin:imax,jmin:jmax,kmin:kmax)
                                                    !水平方向の数値拡散
    
    xyz_NumDiffKm =   &
      &  - NuH * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_Scalar ))))) &
      &  - NuH * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_Scalar ))))) &
      &  - NuV * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_Scalar ))))) 
    
  end function xyz_NumDiffKm

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyza_NumDiffScalar(xyza_Scalar)
    !
    ! x, z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax, ncmax
    use differentiate_center2,                    &
      &             only: xyz_dx_pyz, pyz_dx_xyz, &
      &                   xyz_dy_xqz, xqz_dy_xyz, &
      &                   xyz_dz_xyr, xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(8), intent(in) :: xyza_Scalar(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                    !スカラー量
    real(8)             :: xyza_NumDiffScalar(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                    !水平方向の数値拡散
    integer             :: s
    
    do s = 1, ncmax
      xyza_NumDiffScalar(:,:,:,s) =   &
        &  - NuH * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyza_Scalar(:,:,:,s) ))))) &
        &  - NuH * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyza_Scalar(:,:,:,s) ))))) &
        &  - NuV * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyza_Scalar(:,:,:,s) ))))) 
    end do

  end function xyza_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function pyz_NumDiffVelX(pyz_VarX)
    !
    ! z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax
    use differentiate_center2,                    &
      &             only: pyz_dx_xyz, xyz_dx_pyz, &
      &                   pyz_dy_pqz, pqz_dy_pyz, &
      &                   pyz_dz_pyr, pyr_dz_pyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pyz_VarX(imin:imax,jmin:jmax,kmin:kmax)
                                                    !物理量
    real(8)             :: pyz_NumDiffVelX(imin:imax,jmin:jmax,kmin:kmax)
                                                    !数値拡散
    
    !速度の拡散係数だけ変更する場合にAlpha_Velocityの値を1以外の値にする
    pyz_NumDiffVelX = &
      & - Alpha_Velocity * NuH * ( pyz_dx_xyz( xyz_dx_pyz( pyz_dx_xyz( xyz_dx_pyz( pyz_VarX ) ) ) ) ) &
      & - Alpha_Velocity * NuH * ( pyz_dy_pqz( pqz_dy_pyz( pyz_dy_pqz( pqz_dy_pyz( pyz_VarX ) ) ) ) ) &
      & - Alpha_Velocity * NuV * ( pyz_dz_pyr( pyr_dz_pyz( pyz_dz_pyr( pyr_dz_pyz( pyz_VarX ) ) ) ) )
    
  end function pyz_NumDiffVelX
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xqz_NumDiffVelY(xqz_VarY)
    !
    ! z 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax
    use differentiate_center2,                    &
      &             only: xqz_dx_pqz, pqz_dx_xqz, &
      &                   xqz_dy_xyz, xyz_dy_xqz, &
      &                   xqz_dz_xqr, xqr_dz_xqz
      
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: xqz_VarY(imin:imax,jmin:jmax,kmin:kmax)
                                                    !物理量
    real(8)             :: xqz_NumDiffVelY(imin:imax,jmin:jmax,kmin:kmax)
                                                    !数値拡散
    
    !速度の拡散係数だけ変更する場合にAlpha_Velocityの値を1以外の値にする
    xqz_NumDiffVelY = &
      & - Alpha_Velocity * NuH * ( xqz_dx_pqz(pqz_dx_xqz(xqz_dx_pqz(pqz_dx_xqz( xqz_VarY ))))) &
      & - Alpha_Velocity * NuH * ( xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz( xqz_VarY ))))) &
      & - Alpha_Velocity * NuV * ( xqz_dz_xqr(xqr_dz_xqz(xqz_dz_xqr(xqr_dz_xqz( xqz_VarY )))))
    
  end function xqz_NumDiffVelY
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyr_NumDiffVelZ(xyr_VarZ)
    !
    ! x 方向に半格子ずれた点での数値拡散項を評価
    !
    
    !モジュール読み込み
    use dc_types,   only: DP
    use gridset,    only: imin, imax, jmin, jmax, kmin, kmax
    use differentiate_center2,                    &
      &             only: xyr_dx_pyr, pyr_dx_xyr, &
      &                   xyr_dy_xqr, xqr_dy_xyr, &
      &                   xyr_dz_xyz, xyz_dz_xyr

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: xyr_VarZ(imin:imax,jmin:jmax,kmin:kmax)
                                                    !物理量
    real(8)             :: xyr_NumDiffVelZ(imin:imax,jmin:jmax,kmin:kmax)
                                                    !数値拡散

    !速度の拡散係数だけ変更する場合にAlpha_Velocityの値を1以外の値にする
    xyr_NumDiffVelZ = &
      & - Alpha_Velocity * NuH * ( xyr_dx_pyr(pyr_dx_xyr(xyr_dx_pyr(pyr_dx_xyr( xyr_VarZ ))))) &
      & - Alpha_Velocity * NuH * ( xyr_dy_xqr(xqr_dy_xyr(xyr_dy_xqr(xqr_dy_xyr( xyr_VarZ ))))) & 
      & - Alpha_Velocity * NuV * ( xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr( xyr_VarZ )))))
    
  end function xyr_NumDiffVelZ

end module NumDiffusion4th
