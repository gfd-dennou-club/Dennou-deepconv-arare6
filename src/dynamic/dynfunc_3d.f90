!= Module DynFunc_3D
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu 
! Version::   $Id: dynfunc_3d.f90,v 1.5 2008-06-26 09:24:27 odakker2 Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
!モデルの力学過程を計算するために必要となる関数群を束ねたモジュール
!具体的には以下の項を計算するための関数を格納する.  
!  * 移流項
!  * 浮力項
!  * 気圧傾度力項
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!  * エクスナー関数の空間方向の離散化において, 2 次精度の離散化を陽に利用しているため,
!    気圧傾度力項の計算プログラムにおいて differentiate_center4 モジュールを指定する
!    ことはできないので注意.
!
!== Future Plans
!

module DynFunc_3d
  !
  !陽開放を用いた力学過程の各項の計算モジュール. 
  !具体的には以下の項を計算するための関数を格納する.  
  !  * 移流項
  !  * 浮力項
  !  * 気圧傾度力項
  !

  !モジュール読み込み
  use dc_types, only: DP
  
  use gridset, only:  DimXMin,           &! x 方向の配列の下限
    &                 DimXMax,           &! x 方向の配列の上限
    &                 DimYMin,           &! y 方向の配列の下限
    &                 DimYMax,           &! y 方向の配列の上限
    &                 DimZMin,           &! z 方向の配列の下限
    &                 DimZMax,           &! z 方向の配列の上限
    &                 SpcNum              !
  use damping_3d_v2,  only: DampSound           !音波の減衰係数
  use basicset_3d, only: xyz_EffMolWtBasicZ, &
    &                    xyz_PotTempBasicZ   !基本場の温位
  use constants, only: CpDry,             &!乾燥成分の比熱
    &                  Grav                !重力加速度
  use xyz_base_module, only: &
    &                xyz_avr_pyz, xyr_avr_pyr, xqz_avr_pqz, &
    &                pyz_avr_xyz, pyr_avr_xyr, pqz_avr_xqz, &
    &                xyz_avr_xqz, pyz_avr_pqz, xyr_avr_xqr, &
    &                xqz_avr_xyz, pqz_avr_pyz, xqr_avr_xyr, &
    &                xyz_avr_xyr, pyz_avr_pyr, xqz_avr_xqr, &
    &                xyr_avr_xyz, pyr_avr_pyz, xqr_avr_xqz 
!  use StorePotTemp_3d,  only: StorePotTempAdv
!  use StoreMixRt_3d,    only: StoreMixRtAdv

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !移流項計算のための関数を public にする
  public xyz_AdvScalar
  public xyz_AdvKm
  public xyza_AdvScalar
  public pyz_AdvVelX
  public xqz_AdvVelY
  public xyr_AdvVelZ

  !浮力項計算のための関数を public にする
  public xyr_Buoy

  !気圧傾度力の計算のための関数を public にする
  public pyz_GradPi
  public xqz_GradPi


contains


!!!------------------------------------------------------------------------!!!
  function xyz_AdvScalar(xyz_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP), intent(in) :: xyz_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !スカラー量
    real(DP)             :: xyz_AdvScalar &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !スカラー量の水平移流
    
    xyz_AdvScalar =                                    &
      & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyz_Var))  &
      & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyz_Var))  &
      & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyz_Var))    

!    call StorePotTempAdv( xyz_AdvScalar )   

  end function xyz_AdvScalar


!!!------------------------------------------------------------------------!!!
  function xyz_AdvKm(xyz_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP), intent(in) :: xyz_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !スカラー量
    real(DP)             :: xyz_AdvKm &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !スカラー量の水平移流
    
    xyz_AdvKm =                               &
      & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyz_Var))  &
      & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyz_Var))  &
      & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyz_Var))    

  end function xyz_AdvKm

!!!------------------------------------------------------------------------!!!
  function xyza_AdvScalar(xyza_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz
!    use xyz_deriv_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP), intent(in) :: xyza_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax,SpcNum)
                                                        !スカラー量
    real(DP)             :: xyza_AdvScalar &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax,SpcNum)
                                                        !スカラー量の水平移流
    integer :: s ! ループ変数

    do s = 1, SpcNum        
      xyza_AdvScalar(:,:,:,s) =                                  &
        & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyza_Var(:,:,:,s))) &
        & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyza_Var(:,:,:,s))) &
        & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyza_Var(:,:,:,s)))    
    end do

!    call StoreMixRtAdv( xyza_AdvScalar )   

  end function xyza_AdvScalar
!!!------------------------------------------------------------------------!!!
  function pyz_AdvVelX(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: xyz_dx_pyz, pqz_dy_pyz, pyr_dz_pyz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP)             :: pyz_AdvVelX(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風の移流
    
!    pz_AdvVelX = 0.0d0  !初期化

    pyz_AdvVelX =                                                 &
      & - pyz_VelX * pyz_avr_xyz( xyz_dx_pyz( pyz_VelX ) )        &
      & - pyz_avr_pqz( pqz_avr_xqz( xqz_VelY ) * pqz_dy_pyz( pyz_VelX ) ) &
      & - pyz_avr_pyr( pyr_avr_xyr( xyr_VelZ ) * pyr_dz_pyz( pyz_VelX ) )

    
  end function pyz_AdvVelX


!!!------------------------------------------------------------------------!!!
  function xyr_AdvVelZ(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: pyr_dx_xyr, xqr_dy_xyr, xyz_dz_xyr
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP)             :: xyr_AdvVelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風の移流
  
!    xr_AdvVelZ = 0.0d0  !初期化
    xyr_AdvVelZ =                                                 &
      & - xyr_avr_pyr( pyr_avr_pyz( pyz_VelX ) * pyr_dx_xyr( xyr_VelZ ) ) &
      & - xyr_avr_xqr( xqr_avr_xqz( xqz_VelY ) * xqr_dy_xyr( xyr_VelZ ) ) &
      & - xyr_VelZ * xyr_avr_xyz( xyz_dz_xyr( xyr_VelZ ) )
    
  end function xyr_AdvVelZ

!!!------------------------------------------------------------------------!!!
  function xqz_AdvVelY(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use xyz_deriv_c4_module, only: pqz_dx_xqz, xyz_dy_xqz, xqr_dz_xqz
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !水平風速
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風速
    real(DP)             :: xqz_AdvVelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !鉛直風の移流
  
!    xr_AdvVelZ = 0.0d0  !初期化
    xqz_AdvVelY =                                                 &
      & - xqz_avr_pqz( pqz_avr_pyz( pyz_VelX ) * pqz_dx_xqz( xqz_VelY ) ) &
      & - xqz_VelY * xqz_avr_xyz( xyz_dy_xqz( xqz_VelY ) ) &
      & - xqz_avr_xqr( xqr_avr_xyr( xyr_VelZ ) * xqr_dz_xqz( xqz_VelY ) )
    
  end function xqz_AdvVelY

 
!!!------------------------------------------------------------------------!!!
  function xyr_Buoy(xyz_PotTemp)
    !
    ! 鉛直方向の運動方程式に現れる浮力項を計算

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xyz_PotTemp(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !温位擾乱
    real(DP)              :: xyr_Buoy &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !浮力項

!    !初期化
!    xr_Buoy = 0.0d0

    !浮力項の計算
    xyr_Buoy = Grav * xyr_avr_xyz(xyz_PotTemp / xyz_PotTempBasicZ)

  end function xyr_Buoy


!!!------------------------------------------------------------------------!!!
  function pyz_GradPi(xyz_Exner, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! z 方向に半格子ずれた点での圧力傾度力項の計算. 
    ! 音波減衰項を含めた形式で定式化してあることに注意.
    !
    
    !モジュール読み込み

    use xyz_deriv_module, only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, pyz_dx_xyz
        
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xyz_Exner &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !エクスナー関数の擾乱
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !水平速度
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !水平速度
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !鉛直速度
    real(DP)             :: pyz_GradPi &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !圧力傾度力
    real(DP)             :: xyz_DivVel &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !速度の収束

    !速度の収束
    xyz_DivVel = xyz_dx_pyz( pyz_VelX ) + xyz_dy_xqz( xqz_VelY ) &
      &          + xyz_dz_xyr( xyr_VelZ )
    
    !圧力傾度
!    pyz_GradPi = 0.0d0
    pyz_GradPi =  &
      & + pyz_avr_xyz( CpDry * xyz_PotTempBasicZ / xyz_EffMolWtBasicZ )    &
      &   * pyz_dx_xyz( xyz_Exner ) &
      & - pyz_dx_xyz( DampSound * xyz_DivVel )

  end function pyz_GradPi

!!!------------------------------------------------------------------------!!!
  function xqz_GradPi(xyz_Exner, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! y 方向に半格子ずれた点での圧力傾度力項の計算. 
    ! 音波減衰項を含めた形式で定式化してあることに注意.
    !
    
    !モジュール読み込み

    use xyz_deriv_module, only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, xqz_dy_xyz
        
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: xyz_Exner &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !エクスナー関数の擾乱
    real(DP), intent(in)  :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !水平速度
    real(DP), intent(in)  :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !水平速度
    real(DP), intent(in)  :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !鉛直速度
    real(DP)              :: xqz_GradPi &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !圧力傾度力
    real(DP)              :: xyz_DivVel &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !速度の収束

    !速度の収束
    xyz_DivVel = xyz_dx_pyz( pyz_VelX ) + xyz_dy_xqz( xqz_VelY ) &
      &          + xyz_dz_xyr( xyr_VelZ )
    
    !圧力傾度
!    xqz_GradPi = 0.0d0
    xqz_GradPi =  &
      & + xqz_avr_xyz( CpDry * xyz_PotTempBasicZ / xyz_EffMolWtBasicZ )    &
      &   * xqz_dy_xyz( xyz_Exner )               &
      & - xqz_dy_xyz( DampSound * xyz_DivVel )
    
  end function xqz_GradPi
  
end module DynFunc_3d
