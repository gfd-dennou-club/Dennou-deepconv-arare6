!= Module DynFunc
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: dynfunc.f90,v 1.12 2011-02-28 12:00:23 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
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
!    気圧傾度力項の計算プログラムにおいて differentiate_center4 モジュールを指定することはできないので注意.
!
!== Future Plans
!

module DynFunc
  !
  !陽開放を用いた力学過程の各項の計算モジュール. 
  !具体的には以下の項を計算するための関数を格納する.  
  !  * 移流項
  !  * 浮力項
  !  * 気圧傾度力項
  !

  !モジュール読み込み
  use gridset, only:  DimXMin,           &! x 方向の配列の下限
    &                 DimXMax,           &! x 方向の配列の上限
    &                 DimZMin,           &! z 方向の配列の下限
    &                 DimZMax,           &! z 方向の配列の上限
    &                 SpcNum              !
  use damping,  only: DampSound           !音波の減衰係数
  use basicset, only: xz_PotTempBasicZ,  &!基本場の温位
    &                 xz_EffMolWtBasicZ   !基本場の分子量効果
  use constants,only: CpDry,             &!乾燥成分の比熱
    &                 Grav                !重力加速度
  use average,  only: xz_avr_pz, xz_avr_xr, &
    &                 pz_avr_xz, pz_avr_pr, &
    &                 pr_avr_xr, pr_avr_pz, &
    &                 xr_avr_pr, xr_avr_xz
!  use StorePotTemp,  only: StorePotTempAdv
!  use StoreMixRt,    only: StoreMixRtAdv
!  use StoreMom,    only: StoreMomAdv
!  use StoreBuoy,    only: StoreBuoyTemp

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !移流項計算のための関数を public にする
  public xz_AdvScalar
  public xz_AdvKm
  public xza_AdvScalar
  public pz_AdvVelX
  public xr_AdvVelZ

  !浮力項計算のための関数を public にする
  public xr_Buoy

  !気圧傾度力の計算のための関数を public にする
  public pz_GradPi


contains


!!!------------------------------------------------------------------------!!!
  function xz_AdvScalar(xz_Var, pz_VelX, xr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !水平風速
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !鉛直風速
    real(8), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量
    real(8)             :: xz_AdvScalar(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量の水平移流
    
    xz_AdvScalar =                               &
      & - xz_avr_pz(pz_VelX * pz_dx_xz(xz_Var))  &
      & - xz_avr_xr(xr_VelZ * xr_dz_xz(xz_Var))    

!    call StorePotTempAdv( xz_AdvScalar )   

  end function xz_AdvScalar



!!!------------------------------------------------------------------------!!!
  function xz_AdvKm(xz_Var, pz_VelX, xr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !水平風速
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !鉛直風速
    real(8), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量
    real(8)             :: xz_AdvKm(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量の水平移流
    
    xz_AdvKm =                                   &
      & - xz_avr_pz(pz_VelX * pz_dx_xz(xz_Var))  &
      & - xz_avr_xr(xr_VelZ * xr_dz_xz(xz_Var))    
    
  end function xz_AdvKm



  function xza_AdvScalar(xza_Var, pz_VelX, xr_VelZ)
    !
    ! x, z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !水平風速
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !鉛直風速
    real(8), intent(in) :: xza_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
                                                        !スカラー量
    real(8)             :: xza_AdvScalar(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
                                                        !スカラー量の水平移流
    integer             :: s

    do s = 1, SpcNum
      xza_AdvScalar(:,:,s) =                               &
        & - xz_avr_pz(pz_VelX * pz_dx_xz(xza_Var(:,:,s)))  &
        & - xz_avr_xr(xr_VelZ * xr_dz_xz(xza_Var(:,:,s)))    
    end do

!    call StoreMixRtAdv( xza_AdvScalar )   
    
  end function xza_AdvScalar
  

!!!------------------------------------------------------------------------!!!
  function pz_AdvVelX(pz_VelX, xr_VelZ)
    !
    ! z 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use differentiate_center4, only: xz_dx_pz, pr_dz_pz
!    use differentiate_center2, only: xz_dx_pz, pr_dz_pz

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !水平風速
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !鉛直風速
    real(8)             :: pz_AdvVelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量の水平移流
    
!    pz_AdvVelX = 0.0d0  !初期化
    pz_AdvVelX =                                                 &
      & - pz_VelX * pz_avr_xz( xz_dx_pz( pz_VelX ) )             &
      & - pz_avr_pr( pr_avr_xr( xr_VelZ ) * pr_dz_pz( pz_VelX ) )
    
!    call StoreMomAdv( pz_AdvVelX )   

  end function pz_AdvVelX


!!!------------------------------------------------------------------------!!!
  function xr_AdvVelZ(xr_VelZ, pz_VelX)
    !
    ! x 方向に半格子ずれた点における移流を計算
    !
    
    !モジュール読み込み
    use differentiate_center4, only: pr_dx_xr, xz_dz_xr
!    use differentiate_center2, only: pr_dx_xr, xz_dz_xr
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !水平風速
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !鉛直風速
    real(8)             :: xr_AdvVelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !スカラー量の水平移流
  
!    xr_AdvVelZ = 0.0d0  !初期化
    xr_AdvVelZ =                                                 &
      & - xr_avr_pr( pr_avr_pz( pz_VelX ) * pr_dx_xr( xr_VelZ ) ) &
      & - xr_VelZ * xr_avr_xz( xz_dz_xr( xr_VelZ ) )
    
  end function xr_AdvVelZ
  

!!!------------------------------------------------------------------------!!!
  function xr_Buoy(xz_PotTemp)
    !
    ! 鉛直方向の運動方程式に現れる浮力項を計算
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in)  :: xz_PotTemp(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !温位擾乱
    real(8)              :: xr_Buoy(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !浮力項

!    !初期化
!    xr_Buoy = 0.0d0

    !浮力項の計算
    xr_Buoy = Grav * xr_avr_xz(xz_PotTemp / xz_PotTempBasicZ)

!    call StoreBuoyTemp(xz_avr_xr(xr_Buoy))

  end function xr_Buoy


!!!------------------------------------------------------------------------!!!
  function pz_GradPi(xz_Exner, pz_VelX, xr_VelZ)
    !
    ! z 方向に半格子ずれた点での圧力傾度力項の計算. 
    ! 音波減衰項を含めた形式で定式化してあることに注意.
    !
    
    !モジュール読み込み
    use differentiate_center2,  only: pz_dx_xz, xz_dx_pz, xz_dz_xr
        
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in)  :: xz_Exner(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !エクスナー関数の擾乱
    real(8), intent(in)  :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !水平速度
    real(8), intent(in)  :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !鉛直速度
    real(8)              :: pz_GradPi(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !圧力傾度力
    real(8)              :: xz_DivVel(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !速度の収束

    !速度の収束
    xz_DivVel =  xz_dx_pz( pz_VelX ) + xz_dz_xr( xr_VelZ )
    
    !圧力傾度
!    pz_GradPi = 0.0d0
!    pz_GradPi =  &
!      & pz_avr_xz( CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )    &
!      &   * ( pz_dx_xz( xz_Exner ) - pz_dx_xz( DampSound * xz_DivVel ) )  
    pz_GradPi =  &
      & pz_avr_xz( CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )    &
      &   * pz_dx_xz( xz_Exner ) &
      & - pz_dx_xz( DampSound * xz_DivVel ) 
    
  end function pz_GradPi

!!!--------------------------------------------------------------------!!!
  function xr_GradPi(xz_Exner, pz_VelX, xr_VelZ)
    !
    ! x 方向に半格子ずれた点での圧力傾度力項の計算. 
    ! 音波減衰項を含めた形式で定式化してあることに注意.
    !
    
    !モジュール読み込み
    use differentiate_center2,  only: xr_dz_xz, xz_dx_pz, xz_dz_xr

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(in)  :: xz_Exner(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !エクスナー関数の擾乱
    real(8), intent(in)  :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !水平速度
    real(8), intent(in)  :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !鉛直速度
    real(8)              :: xr_GradPi(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !圧力傾度力
    real(8)              :: xz_DivVel(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !速度の収束

    !速度の収束
    xz_DivVel =  xz_dx_pz( pz_VelX ) + xz_dz_xr( xr_VelZ )
    
    !速度 w の圧力勾配
!    xr_GradPi = 0.0d0
!    xr_GradPi =   &
!      & xr_avr_xz(CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )   &
!      &   * ( xr_dz_xz( xz_Exner ) - xr_dz_xz( DampSound * xz_DivVel ) )
    xr_GradPi =   &
      & xr_avr_xz(CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )   &
      &   * xr_dz_xz( xz_Exner ) &
      & - xr_dz_xz( DampSound * xz_DivVel ) 
    
  end function xr_GradPi
  
end module DynFunc
