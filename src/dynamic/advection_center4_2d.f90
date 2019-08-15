!= Module advection_center4_2d
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro)
! Version::   $Id: advection_center4_2d.f90,v 1.3 2014/07/08 00:55:25 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module advection_center4_2d
  !
  ! 2D 計算用の移流計算モジュール. 
  !
  !   移流: 4 次中央差分
  !   数値拡散 (4 階): 2 次中央差分
  !
  ! リープフロッグで, 移流を中央差分で計算するために, 
  ! 数値拡散項を追加している. 
  ! 
  ! エクスナー関数の基本場の移流は圧力方程式で考慮しているので, ここには含めない.

  !モジュール読み込み
  !
  use dc_types,   only : DP

  !暗黙の型宣言禁止
  !
  implicit none

  ! 変数の設定
  !
  real(DP), save, private :: NuHh  = 0.0d0         !熱に対する数値粘性の係数 (水平方向)
  real(DP), save, private :: NuVh  = 0.0d0         !熱に対する数値粘性の係数 (鉛直方向)
  real(DP), save, private :: NuHm  = 0.0d0         !運動量に対する数値粘性の係数 (水平方向)
  real(DP), save, private :: NuVm  = 0.0d0         !運動量に対する数値粘性の係数 (鉛直方向)
  character(*), parameter :: module_name = 'advection_center4_2d'
                                          ! モジュールの名称.
                                          ! Module name
  real(DP), allocatable, save, private :: xyr_DPTempBZDz(:,:,:)   !基本場の鉛直微分
  real(DP), allocatable, save, private :: xyrf_DQMixBZDz(:,:,:,:) !基本場の鉛直微分

  !public
  !
  public advection_center4_2d_init
  public advection_center4_2d_dry
  public advection_center4_2d_tracer

contains

  subroutine advection_center4_2d_init( AlphaNDiff, NDiffRatio )
    !
    ! 初期化ルーチン
    !

    ! モジュール読み込み
    !
    use timeset,    only : DelTimeLong
    use axesset,    only : dx, dz             ! 格子間隔
    use dc_message, only : MessageNotify
    use gridset,     only : imin, imax,      &
      &                     jmin, jmax,      &
      &                     kmin, kmax,      &
      &                     ncmax
    use basicset,    only : xyz_PTempBZ,     &!基本場の温位
      &                     xyzf_QMixBZ
    use differentiate_center4, &
      &              only : xyr_dz_xyz

    !暗黙の型宣言禁止
    !
    implicit none
    
    ! 変数の定義
    !
    real(DP), intent(in) :: AlphaNDiff  !数値拡散の係数. 
    real(DP), intent(in) :: NDiffRatio
    integer              :: f

    !-------------------------------------------------------------------
    ! 数値拡散係数を決める
    !
    ! CReSS マニュアルの記述に従って NuH, NuV を決める.
    ! 運動量と熱に対する数値拡散の大きさを変えられるように NDiffRatio を乗じている.
    ! 
    NuHh = AlphaNDiff * ( dx ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    NuVh = AlphaNDiff * ( dz ** 4.0d0 ) / (2.0d0 * DelTimeLong)

    NuHm = NuHh * NDiffRatio
    NuVm = NuVh * NDiffRatio

    !-------------------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", module_name, "NuHh = %f", d=(/NuHh/) )
    call MessageNotify( "M", module_name, "NuVh = %f", d=(/NuVh/) )
    call MessageNotify( "M", module_name, "NuHm = %f", d=(/NuHm/) )
    call MessageNotify( "M", module_name, "NuVm = %f", d=(/NuVm/) )

    !-------------------------------------------------------------------
    ! 配列の用意
    !
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

  end subroutine advection_center4_2d_init
  
!!!------------------------------------------------------------------------!!!
  
  subroutine advection_center4_2d_dry( &
    & pyz_VelXB,    pyz_VelXN,          & ! (in)
    & xyr_VelZB,    xyr_VelZN,          & ! (in)
    & xyz_PTempB,   xyz_PTempN,         & ! (in)
    & xyz_ExnerB,   xyz_ExnerN,         & ! (in)
    & xyz_KmB,      xyz_KmN,            & ! (in)
    & pyz_DVelXDtAdv,  pyz_VelXnDiff,      & !(out)
    & xyr_DVelZDtAdv,  xyr_VelZnDiff,      & !(out)
    & xyz_DPTempDtAdv, xyz_PTempNDiff,     & !(out)
    & xyz_DExnerDtAdv, xyz_ExnerNDiff,     & !(out)
    & xyz_DKmDtAdv,    xyz_KmNDiff         & !(out)
    & )
    ! 
    ! 移流計算 (乾燥)
    !
    !   移流: 4 次中央差分
    !   数値拡散 (4 階): 2 次中央差分
    !
    ! リープフロッグで, 移流を中央差分で計算するために, 
    ! 数値拡散項を追加している. 
    !

    !モジュール読み込み
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x 方向の配列の下限
      &                  imax,            &! x 方向の配列の上限
      &                  jmin,            &! y 方向の配列の下限
      &                  jmax,            &! y 方向の配列の上限
      &                  kmin,            &! z 方向の配列の下限
      &                  kmax              ! z 方向の配列の上限
    use axesset,  only : dx, dz            ! 格子間隔

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 配列の定義
    !
    real(DP), intent(in)    :: pyz_VelXB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmN(imin:imax,jmin:jmax,kmin:kmax)

    real(DP), intent(out)   :: pyz_DVelXDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_DVelZDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DKmDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DExnerDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DPTempDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: pyz_VelXnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_KmNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_PTempNDiff(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)                :: fct1, fct2
    integer                 :: i, k
    integer, parameter      :: j = 1
    
      
    ! 微分に用いる係数を予め計算
    !
    fct1 = 9.0d0 / 8.0d0
    fct2 = 1.0d0 / 24.0d0

    !---------------------------------------------------------------------
    ! 速度 U
    ! 

    ! 移流
    !
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
          
          pyz_DVelXDtAdv(i,j,k) =                                             &
            & - pyz_VelXN(i,j,k)                                           &
            &   * (                                                        &
            &       fct1 * (   pyz_VelXN(i+1,j,k) - pyz_VelXN(i-1,j,k) )   &
            &     - fct2 * (   pyz_VelXN(i+2,j,k) + pyz_VelXN(i+1,j,k)     &
            &                - pyz_VelXN(i-1,j,k) - pyz_VelXN(i-2,j,k) )   &
            &     ) * 5.0d-1 / dx                                          &
            & - (                                                          &
            &   + ( xyr_VelZN(i+1,j,k) + xyr_VelZN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( pyz_VelXN(i,j,k+1) - pyz_VelXN(i,j,k)   )   &
            &       - fct2 * ( pyz_VelXN(i,j,k+2) - pyz_VelXN(i,j,k-1) )   &
            &       )                                                      &
            &   + ( xyr_VelZN(i+1,j,k-1) + xyr_VelZN(i,j,k-1) )            & 
            &     * (                                                      &
            &         fct1 * ( pyz_VelXN(i,j,k)   - pyz_VelXN(i,j,k-1) )   &
            &       - fct2 * ( pyz_VelXN(i,j,k+1) - pyz_VelXN(i,j,k-2) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dz
        
      end do
    end do

    ! 数値拡散
    !
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
          
          pyz_VelXnDiff(i,j,k) =                        &
            & - (                                       &
            &     + pyz_VelXB(i+2,j,k)                  &
            &     + pyz_VelXB(i-2,j,k)                  &
            &     - pyz_VelXB(i+1,j,k) * 4.0d0          &
            &     - pyz_VelXB(i-1,j,k) * 4.0d0          &
            &     + pyz_VelXB(i  ,j,k) * 6.0d0          &
            &   ) * NuHm / ( dx ** 4.0d0 )              &
            & - (                                       & 
            &     + pyz_VelXB(i,j,k+2)                  &
            &     + pyz_VelXB(i,j,k-2)                  &
            &     - pyz_VelXB(i,j,k+1) * 4.0d0          &
            &     - pyz_VelXB(i,j,k-1) * 4.0d0          &
            &     + pyz_VelXB(i,j,k  ) * 6.0d0          &
            &   ) * NuVm / ( dz ** 4.0d0 )
        
      end do
    end do

    ! 値の確定
    !
    pyz_DVelXDtAdv(imin:imin+1,:,:) = 1.0d10
    pyz_DVelXDtAdv(imax-1:imax,:,:) = 1.0d10
    pyz_DVelXDtAdv(:,:,kmin:kmin+1) = 1.0d10
    pyz_DVelXDtAdv(:,:,kmax-1:kmax) = 1.0d10

    pyz_VelXnDiff(imin:imin+1,:,:) = 1.0d10
    pyz_VelXnDiff(imax-1:imax,:,:) = 1.0d10
    pyz_VelXnDiff(:,:,kmin:kmin+1) = 1.0d10
    pyz_VelXnDiff(:,:,kmax-1:kmax) = 1.0d10
    
    !---------------------------------------------------------------------
    ! 速度 W
    ! 

    ! 移流
    ! 
    do k = kmin+2, kmax-2
      do i = imin+2, imax-2
        
          xyr_DVelZDtAdv(i,j,k) =                                          &
            & - (                                                          &
            &   + ( pyz_VelXN(i,j,k+1) + pyz_VelXN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i+1,j,k) - xyr_VelZN(i,j,k)   )   &
            &       - fct2 * ( xyr_VelZN(i+2,j,k) - xyr_VelZN(i-1,j,k) )   &
            &       )                                                      &
            &   + ( pyz_VelXN(i-1,j,k+1) + pyz_VelXN(i-1,j,k) )            & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i,j,k)   - xyr_VelZN(i-1,j,k) )   &
            &       - fct2 * ( xyr_VelZN(i+1,j,k) - xyr_VelZN(i-2,j,k) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dx                                            &
            & - xyr_VelZN(i,j,k)                                           &
            &   * (                                                        &
            &       fct1 * (   xyr_VelZN(i,j,k+1) - xyr_VelZN(i,j,k-1) )   &
            &     - fct2 * (   xyr_VelZN(i,j,k+2) + xyr_VelZN(i,j,k+1)     &
            &                - xyr_VelZN(i,j,k-1) - xyr_VelZN(i,j,k-2) )   &
            &     ) * 5.0d-1 / dz
          
      end do
    end do
     
    ! 数値拡散
    !
    do k = kmin+2, kmax-2
      do i = imin+2, imax-2
        
          xyr_VelZnDiff(i,j,k) =                        &
            & - (                                       &
            &     + xyr_VelZB(i+2,j,k)                  &
            &     + xyr_VelZB(i-2,j,k)                  &
            &     - xyr_VelZB(i+1,j,k) * 4.0d0          &
            &     - xyr_VelZB(i-1,j,k) * 4.0d0          &
            &     + xyr_VelZB(i  ,j,k) * 6.0d0          &
            &   ) * NuHm / ( dx ** 4.0d0 )              &
            & - (                                       &
            &     + xyr_VelZB(i,j,k+2)                  &
            &     + xyr_VelZB(i,j,k-2)                  &
            &     - xyr_VelZB(i,j,k+1) * 4.0d0          &
            &     - xyr_VelZB(i,j,k-1) * 4.0d0          &
            &     + xyr_VelZB(i,j,k  ) * 6.0d0          &
            &   ) * NuVm / ( dz ** 4.0d0 )
        
      end do
    end do

    ! 値の確定
    !
    xyr_DVelZDtAdv(imin:imin+1,:,:) = 1.0d10
    xyr_DVelZDtAdv(imax-1:imax,:,:) = 1.0d10
    xyr_DVelZDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyr_DVelZDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyr_VelZnDiff(imin:imin+1,:,:) = 1.0d10
    xyr_VelZnDiff(imax-1:imax,:,:) = 1.0d10
    xyr_VelZnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyr_VelZnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! 温位
    !     

    ! 移流
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
          
          xyz_DPTempDtAdv(i,j,k) =                                            &
            & - (                                                             &
            &      pyz_VelXN(i,j,k)                                           &
            &        * (                                                      &
            &            fct1 * ( xyz_PTempN(i+1,j,k) - xyz_PTempN(i,j,k)   ) &
            &          - fct2 * ( xyz_PTempN(i+2,j,k) - xyz_PTempN(i-1,j,k) ) &
            &          )                                                      &
            &    + pyz_VelXN(i-1,j,k)                                         &
            &        * (                                                      &
            &            fct1 * ( xyz_PTempN(i,j,k)   - xyz_PTempN(i-1,j,k) ) &
            &          - fct2 * ( xyz_PTempN(i+1,j,k) - xyz_PTempN(i-2,j,k) ) &
            &          )                                                      &
            &   ) * 5.0d-1 / dx                                               &
            & - (                                                             &
            &      xyr_VelZN(i,j,k)                                           &
            &        * (                                                      &
            &            fct1 * ( xyz_PTempN(i,j,k+1) - xyz_PTempN(i,j,k)   ) &
            &          - fct2 * ( xyz_PTempN(i,j,k+2) - xyz_PTempN(i,j,k-1) ) &
            &          )                                                      &
            &    + xyr_VelZN(i,j,k-1)                                         &
            &        * (                                                      &
            &            fct1 * ( xyz_PTempN(i,j,k)   - xyz_PTempN(i,j,k-1) ) &
            &          - fct2 * ( xyz_PTempN(i,j,k+1) - xyz_PTempN(i,j,k-2) ) &
            &          )                                                      &
            &   ) * 5.0d-1 / dz                                               &
            & - (                                                             &
            &      xyr_VelZN(i,j,k)   * xyr_DPTempBZDz(i,j,k)                 &
            &    + xyr_VelZN(i,j,k-1) * xyr_DPTempBZDz(i,j,k-1)               & 
            &   ) * 5.0d-1 
        
      end do
    end do

    ! 数値拡散
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
        
          xyz_PTempNDiff(i,j,k) =                        &
            & - (                                        &
            &     + xyz_PTempB(i+2,j,k)                  &
            &     + xyz_PTempB(i-2,j,k)                  &
            &     - xyz_PTempB(i+1,j,k) * 4.0d0          &
            &     - xyz_PTempB(i-1,j,k) * 4.0d0          &
            &     + xyz_PTempB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )               &
            & - (                                        &
            &       xyz_PTempB(i,j,k+2)                  &
            &     + xyz_PTempB(i,j,k-2)                  &
            &     - xyz_PTempB(i,j,k+1) * 4.0d0          &
            &     - xyz_PTempB(i,j,k-1) * 4.0d0          &
            &     + xyz_PTempB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )

      end do
    end do

    ! 値の確定
    !
    xyz_DPTempDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DPTempDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DPTempDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DPTempDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_PTempnDiff(imin:imin+1,:,:) = 1.0d10
    xyz_PTempnDiff(imax-1:imax,:,:) = 1.0d10
    xyz_PTempnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_PTempnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! エクスナー関数
    !     

    ! 移流
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
          
          xyz_DExnerDtAdv(i,j,k) =                                               &
            & - (                                                             &
            &      pyz_VelXN(i,j,k)                                           &
            &        * (                                                      &
            &            fct1 * ( xyz_ExnerN(i+1,j,k) - xyz_ExnerN(i,j,k)   ) &
            &          - fct2 * ( xyz_ExnerN(i+2,j,k) - xyz_ExnerN(i-1,j,k) ) &
            &          )                                                      &
            &    + pyz_VelXN(i-1,j,k)                                         &
            &        * (                                                      &
            &            fct1 * ( xyz_ExnerN(i,j,k)   - xyz_ExnerN(i-1,j,k) ) &
            &          - fct2 * ( xyz_ExnerN(i+1,j,k) - xyz_ExnerN(i-2,j,k) ) &
            &          )                                                      &
            &   ) * 5.0d-1 / dx                                               &
            & - (                                                             &
            &      xyr_VelZN(i,j,k)                                           &
            &        * (                                                      &
            &            fct1 * ( xyz_ExnerN(i,j,k+1) - xyz_ExnerN(i,j,k)   ) &
            &          - fct2 * ( xyz_ExnerN(i,j,k+2) - xyz_ExnerN(i,j,k-1) ) &
            &          )                                                      &
            &    + xyr_VelZN(i,j,k-1)                                         &
            &        * (                                                      &
            &            fct1 * ( xyz_ExnerN(i,j,k)   - xyz_ExnerN(i,j,k-1) ) &
            &          - fct2 * ( xyz_ExnerN(i,j,k+1) - xyz_ExnerN(i,j,k-2) ) &
            &          )                                                      &
            &   ) * 5.0d-1 / dz 
        
      end do
    end do

    ! 数値拡散
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
        
          xyz_ExnerNDiff(i,j,k) =                        &
            & - (                                        &
            &     + xyz_ExnerB(i+2,j,k)                  &
            &     + xyz_ExnerB(i-2,j,k)                  &
            &     - xyz_ExnerB(i+1,j,k) * 4.0d0          &
            &     - xyz_ExnerB(i-1,j,k) * 4.0d0          &
            &     + xyz_ExnerB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )               &
            & - (                                        &
            &       xyz_ExnerB(i,j,k+2)                  &
            &     + xyz_ExnerB(i,j,k-2)                  &
            &     - xyz_ExnerB(i,j,k+1) * 4.0d0          &
            &     - xyz_ExnerB(i,j,k-1) * 4.0d0          &
            &     + xyz_ExnerB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )

      end do
    end do

    ! 値の確定
    !
    xyz_DExnerDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DExnerDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DExnerDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DExnerDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_ExnernDiff(imin:imin+1,:,:) = 1.0d10
    xyz_ExnernDiff(imax-1:imax,:,:) = 1.0d10
    xyz_ExnernDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_ExnernDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! 乱流拡散係数
    !     

    ! 移流
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
          
          xyz_DKmDtAdv(i,j,k) =                                            &
            & - (                                                       &
            &      pyz_VelXN(i,j,k)                                     &
            &        * (                                                &
            &            fct1 * ( xyz_KmN(i+1,j,k) - xyz_KmN(i,j,k)   ) &
            &          - fct2 * ( xyz_KmN(i+2,j,k) - xyz_KmN(i-1,j,k) ) &
            &          )                                                &
            &    + pyz_VelXN(i-1,j,k)                                   &
            &        * (                                                &
            &            fct1 * ( xyz_KmN(i,j,k)   - xyz_KmN(i-1,j,k) ) &
            &          - fct2 * ( xyz_KmN(i+1,j,k) - xyz_KmN(i-2,j,k) ) &
            &          )                                                &
            &   ) * 5.0d-1 / dx                                         &
            & - (                                                       &
            &      xyr_VelZN(i,j,k)                                     &
            &        * (                                                &
            &            fct1 * ( xyz_KmN(i,j,k+1) - xyz_KmN(i,j,k)   ) &
            &          - fct2 * ( xyz_KmN(i,j,k+2) - xyz_KmN(i,j,k-1) ) &
            &          )                                                &
            &    + xyr_VelZN(i,j,k-1)                                   &
            &        * (                                                &
            &            fct1 * ( xyz_KmN(i,j,k)   - xyz_KmN(i,j,k-1) ) &
            &          - fct2 * ( xyz_KmN(i,j,k+1) - xyz_KmN(i,j,k-2) ) &
            &          )                                                &
            &   ) * 5.0d-1 / dz
        
      end do
    end do

    ! 数値拡散
    ! 
    do k = kmin + 2, kmax - 2
      do i = imin + 2, imax - 2
        
          xyz_KmNDiff(i,j,k) =                        &
            & - (                                     &
            &     + xyz_KmB(i+2,j,k)                  &
            &     + xyz_KmB(i-2,j,k)                  &
            &     - xyz_KmB(i+1,j,k) * 4.0d0          &
            &     - xyz_KmB(i-1,j,k) * 4.0d0          &
            &     + xyz_KmB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )            &
            & - (                                     &
            &       xyz_KmB(i,j,k+2)                  &
            &     + xyz_KmB(i,j,k-2)                  &
            &     - xyz_KmB(i,j,k+1) * 4.0d0          &
            &     - xyz_KmB(i,j,k-1) * 4.0d0          &
            &     + xyz_KmB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )
          
      end do
    end do

    ! 値の確定
    !
    xyz_DKmDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DKmDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DKmDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DKmDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_KmnDiff(imin:imin+1,:,:) = 1.0d10
    xyz_KmnDiff(imax-1:imax,:,:) = 1.0d10
    xyz_KmnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_KmnDiff(:,:,kmax-1:kmax) = 1.0d10
       
  end subroutine advection_center4_2d_dry

!!!------------------------------------------------------------------------!!!
  
  subroutine advection_center4_2d_tracer(&
    & pyz_VelXN,   xyr_VelZN,          & ! (in)
    & xyzf_QMixB,  xyzf_QMixN,         & ! (in)
    & xyzf_QMixAdv, xyzf_QMixNDiff     & !(out)
    & )
    ! 
    ! 移流計算
    !
    !   移流: 4 次中央差分
    !   数値拡散 (4 階): 2 次中央差分
    !
    ! リープフロッグで, 移流を中央差分で計算するために, 
    ! 数値拡散項を追加している. 
    !

    !モジュール読み込み
    !
    use dc_types, only: DP
    use gridset, only :  imin,            &! x 方向の配列の下限
      &                  imax,            &! x 方向の配列の上限
      &                  jmin,            &! y 方向の配列の下限
      &                  jmax,            &! y 方向の配列の上限
      &                  kmin,            &! z 方向の配列の下限
      &                  kmax,            &! z 方向の配列の上限
      &                  ncmax           
    use axesset,  only : dx, dz            ! 格子間隔

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 配列の定義
    !
    real(DP), intent(in)    :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyzf_QMixB(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixN(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(out)   :: xyzf_QMixAdv(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out)   :: xyzf_QMixNDiff(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)

    real(DP)                :: fct1, fct2
    integer                 :: i, j, k, s

      
    ! 微分に用いる係数を予め計算
    !
    fct1 = 9.0d0 / 8.0d0
    fct2 = 1.0d0 / 24.0d0
    
    j = 1

    !---------------------------------------------------------------------
    ! 混合比
    !     
   
    ! 移流
    ! 
    do s = 1, ncmax
      do k = kmin + 2, kmax - 2
        do i = imin + 2, imax - 2
          
            xyzf_QMixAdv(i,j,k,s) =                                                 &
              & - (                                                                 &
              &      pyz_VelXN(i,j,k)                                               &
              &        * (                                                          &
              &            fct1 * ( xyzf_QMixN(i+1,j,k,s) - xyzf_QMixN(i,j,k,s)   ) &
              &          - fct2 * ( xyzf_QMixN(i+2,j,k,s) - xyzf_QMixN(i-1,j,k,s) ) &
              &          )                                                          &
              &    + pyz_VelXN(i-1,j,k)                                             &
              &        * (                                                          &
              &            fct1 * ( xyzf_QMixN(i,j,k,s)   - xyzf_QMixN(i-1,j,k,s) ) &
              &          - fct2 * ( xyzf_QMixN(i+1,j,k,s) - xyzf_QMixN(i-2,j,k,s) ) &
              &          )                                                          &
              &   ) * 5.0d-1 / dx                                                   &
              & - (                                                                 &
              &      xyr_VelZN(i,j,k)                                               &
              &        * (                                                          &
              &            fct1 * ( xyzf_QMixN(i,j,k+1,s) - xyzf_QMixN(i,j,k,s)   ) &
              &          - fct2 * ( xyzf_QMixN(i,j,k+2,s) - xyzf_QMixN(i,j,k-1,s) ) &
              &          )                                                          &
              &    + xyr_VelZN(i,j,k-1)                                             &
              &        * (                                                          &
              &            fct1 * ( xyzf_QMixN(i,j,k,s)   - xyzf_QMixN(i,j,k-1,s) ) &
              &          - fct2 * ( xyzf_QMixN(i,j,k+1,s) - xyzf_QMixN(i,j,k-2,s) ) &
              &          )                                                          &
              &   ) * 5.0d-1 / dz                                                   &
              & - (                                                                 &
              &      xyr_VelZN(i,j,k)   * xyrf_DQMixBZDz(i,j,k,s)                   &
              &    + xyr_VelZN(i,j,k-1) * xyrf_DQMixBZDz(i,j,k-1,s)                 &
              &   ) * 5.0d-1 
            
        end do
      end do
    end do
    
    ! 数値拡散
    ! 
    do s = 1, ncmax
      do k = kmin + 2, kmax - 2
        do i = imin + 2, imax - 2
          
            xyzf_QMixNDiff(i,j,k,s) =                        &
              & - (                                          &
              &       xyzf_QMixB(i+2,j,k,s)                  &
              &     + xyzf_QMixB(i-2,j,k,s)                  &
              &     - xyzf_QMixB(i+1,j,k,s) * 4.0d0          &
              &     - xyzf_QMixB(i-1,j,k,s) * 4.0d0          &
              &     + xyzf_QMixB(i  ,j,k,s) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )                 &
              & - (                                          &
              &       xyzf_QMixB(i,j,k+2,s)                  &
              &     + xyzf_QMixB(i,j,k-2,s)                  &
              &     - xyzf_QMixB(i,j,k+1,s) * 4.0d0          &
              &     - xyzf_QMixB(i,j,k-1,s) * 4.0d0          &
              &     + xyzf_QMixB(i,j,k  ,s) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
          
          end do
      end do
    end do

    ! 値の確定
    !
    xyzf_QMixAdv(imin:imin+1,:,:,:) = 1.0d10
    xyzf_QMixAdv(imax-1:imax,:,:,:) = 1.0d10
    xyzf_QMixAdv(:,:,kmin:kmin+1,:) = 1.0d10
    xyzf_QMixAdv(:,:,kmax-1:kmax,:) = 1.0d10

    xyzf_QMixnDiff(imin:imin+1,:,:,:) = 1.0d10
    xyzf_QMixnDiff(imax-1:imax,:,:,:) = 1.0d10
    xyzf_QMixnDiff(:,:,kmin:kmin+1,:) = 1.0d10
    xyzf_QMixnDiff(:,:,kmax-1:kmax,:) = 1.0d10
    
  end subroutine advection_center4_2d_tracer
  
end module advection_center4_2d
