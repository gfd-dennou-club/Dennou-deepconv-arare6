!= Module advection_center4_std
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro)
! Version::   $Id: advection_center4_std.f90,v 1.3 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module advection_center4_std
  !
  ! 移流計算モジュール. 微分平均演算モジュールを利用. 
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
  character(*), parameter :: module_name = 'advection_center4_std'
                                          ! モジュールの名称.
                                          ! Module name
  real(DP), allocatable, save, private :: xyr_DPTempBZDz(:,:,:)   !基本場の鉛直微分
  real(DP), allocatable, save, private :: xyrf_DQMixBZDz(:,:,:,:) !基本場の鉛直微分

  !public
  !
  public advection_center4_std_init
  public advection_center4_std_main

contains

  subroutine advection_center4_std_init( AlphaNDiff, NDiffRatio )
    !
    ! 初期化ルーチン
    !

    ! モジュール読み込み
    !
    use timeset,     only : DelTimeLong
    use axesset,     only : dx, dy, dz       ! 格子間隔
    use dc_message,  only : MessageNotify
    use gridset,     only : imin, imax,      &
      &                     jmin, jmax,      &
      &                     kmin, kmax,      &
      &                     ncmax,           &
      &                     FlagCalc3D
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
!    NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    if ( FlagCalc3D ) then 
      NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    else
      NuHh = AlphaNDiff * ( dx ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    end if
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


  end subroutine advection_center4_std_init

!!!------------------------------------------------------------------------!!!
  
  subroutine advection_center4_std_main(    &
    & pyz_VelXB,    pyz_VelXN,          & ! (in)
    & xqz_VelYB,    xqz_VelYN,          & ! (in)
    & xyr_VelZB,    xyr_VelZN,          & ! (in)
    & xyz_PTempB,   xyz_PTempN,         & ! (in)
    & xyz_ExnerB,   xyz_ExnerN,         & ! (in)
    & xyzf_QMixB,   xyzf_QMixN,         & ! (in)
    & xyz_KmB,      xyz_KmN,            & ! (in)
    & pyz_DVelXDtAdv,  pyz_VelXnDiff,      & !(out)
    & xqz_DVelYDtAdv,  xqz_VelYnDiff,      & !(out)
    & xyr_DVelZDtAdv,  xyr_VelZnDiff,      & !(out)
    & xyz_DPTempDtAdv, xyz_PTempNDiff,     & !(out)
    & xyz_DExnerDtAdv, xyz_ExnerNDiff,     & !(out)
    & xyzf_QMixAdv, xyzf_QMixNDiff,     & !(out)
    & xyz_DKmDtAdv,    xyz_KmNDiff         & !(out)
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
    use dc_types, only : DP
    use gridset,  only : imin,            &! x 方向の配列の下限
      &                  imax,            &! x 方向の配列の上限
      &                  jmin,            &! y 方向の配列の下限
      &                  jmax,            &! y 方向の配列の上限
      &                  kmin,            &! z 方向の配列の下限
      &                  kmax,            &! z 方向の配列の上限
      &                  ncmax
    use average,  only : pyz_xyz, pyz_pqz, pqz_xqz, pyz_pyr, pyr_xyr, &
      &                  xqz_pqz, pqz_pyz, xqz_xyz, xqz_xqr, xqr_xyr, &
      &                  xyr_pyr, pyr_pyz, xyr_xqr, xqr_xqz, xyr_xyz, &
      &                  xyz_pyz, xyz_xqz, xyz_xyr
    use differentiate_center4,                                                                              &
      &          only : pyz_c4dx_xyz => pyz_dx_xyz, xqz_c4dy_xyz => xqz_dy_xyz, xyr_c4dz_xyz => xyr_dz_xyz, &
      &                 xyz_c4dx_pyz => xyz_dx_pyz, pqz_c4dy_pyz => pqz_dy_pyz, pyr_c4dz_pyz => pyr_dz_pyz, &
      &                 pqz_c4dx_xqz => pqz_dx_xqz, xyz_c4dy_xqz => xyz_dy_xqz, xqr_c4dz_xqz => xqr_dz_xqz, &
      &                 pyr_c4dx_xyr => pyr_dx_xyr, xqr_c4dy_xyr => xqr_dy_xyr, xyz_c4dz_xyr => xyz_dz_xyr
!    use differentiate_center2,                              &
!      &          only : pyz_dx_xyz, xyz_dx_pyz, pyz_dy_pqz, &
!      &                 pqz_dy_pyz, pyz_dz_pyr, pyr_dz_pyz, &
!      &                 xqz_dx_pqz, pqz_dx_xqz, xqz_dy_xyz, &
!      &                 xyz_dy_xqz, xqz_dz_xqr, xqr_dz_xqz, &
!      &                 xyr_dx_pyr, pyr_dx_xyr, xyr_dy_xqr, &
!      &                 xqr_dy_xyr, xyr_dz_xyz, xyz_dz_xyr
    use differentiate_center2,                              &
      &          only : xyz_dx4_xyz, pyz_dx4_pyz, xqz_dx4_xqz, xyr_dx4_xyr, &
      &                 xyz_dy4_xyz, pyz_dy4_pyz, xqz_dy4_xqz, xyr_dy4_xyr, &
      &                 xyz_dz4_xyz, pyz_dz4_pyz, xqz_dz4_xqz, xyr_dz4_xyr

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 配列の定義
    !
    real(DP), intent(in)    :: pyz_VelXB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixB(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixN(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmN(imin:imax,jmin:jmax,kmin:kmax)

    real(DP), intent(out)   :: pyz_DVelXDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_DVelYDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_DVelZDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DKmDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DExnerDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DPTempDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyzf_QMixAdv(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out)   :: pyz_VelXnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_KmNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_PTempNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyzf_QMixNDiff(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)

    integer                 :: f


    !------------------------------------------------------------------------
    ! 乱流拡散係数
    !

    ! Advection term
    !
    xyz_DKmDtAdv =                                           &
      & - xyz_pyz( pyz_VelXN * pyz_c4dx_xyz( xyz_KmN ) )  &
      & - xyz_xqz( xqz_VelYN * xqz_c4dy_xyz( xyz_KmN ) )  &
      & - xyz_xyr( xyr_VelZN * xyr_c4dz_xyz( xyz_KmN ) )  

    ! numerical diffusion term
    !
!    xyz_KmNDiff =                                                             &
!      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_KmB ))))) &
!      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_KmB ))))) &
!      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_KmB ))))) 
    xyz_KmNDiff =                                                             &
      &  - NuHh * xyz_dx4_xyz( xyz_KmB )  &
      &  - NuHh * xyz_dy4_xyz( xyz_KmB )  &
      &  - NuVh * xyz_dz4_xyz( xyz_KmB )
        
    !------------------------------------------------------------------------
    ! 温位
    !

    ! Advection term
    !
    xyz_DPTempDtAdv =                                           &
      & - xyz_pyz( pyz_VelXN * pyz_c4dx_xyz( xyz_PTempN ) )  &
      & - xyz_xqz( xqz_VelYN * xqz_c4dy_xyz( xyz_PTempN ) )  &
      & - xyz_xyr( xyr_VelZN * xyr_c4dz_xyz( xyz_PTempN ) )  &
      & - xyz_xyr( xyr_VelZN * xyr_DPTempBZDz ) 
    
    ! numerical diffusion term
    !
!    xyz_PTempNDiff =                                                             &
!      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_PTempB ))))) &
!      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_PTempB ))))) &
!      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_PTempB ))))) 
    xyz_PTempNDiff =                                                             &
      &  - NuHh * xyz_dx4_xyz( xyz_PTempB )  &
      &  - NuHh * xyz_dy4_xyz( xyz_PTempB )  &
      &  - NuVh * xyz_dz4_xyz( xyz_PTempB )

    !------------------------------------------------------------------------
    ! エクスナー関数
    !

    ! Advection term
    !
    xyz_DExnerDtAdv =                                           &
      & - xyz_pyz( pyz_VelXN * pyz_c4dx_xyz( xyz_ExnerN ) )  &
      & - xyz_xqz( xqz_VelYN * xqz_c4dy_xyz( xyz_ExnerN ) )  &
      & - xyz_xyr( xyr_VelZN * xyr_c4dz_xyz( xyz_ExnerN ) ) 

    ! numerical diffusion term
    !
!    xyz_ExnerNDiff =                                                             &
!      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_ExnerB ))))) &
!      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_ExnerB ))))) &
!      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_ExnerB )))))     
    xyz_ExnerNDiff =                                                             &
      &  - NuHh * xyz_dx4_xyz( xyz_ExnerB )  &
      &  - NuHh * xyz_dy4_xyz( xyz_ExnerB )  &
      &  - NuVh * xyz_dz4_xyz( xyz_ExnerB )
    
    
    !---------------------------------------------------------------------
    ! 混合比
    ! 
    
    do f = 1, ncmax

      ! Advection term
      !
      xyzf_QMixAdv(:,:,:,f) =                                          &
        & - xyz_pyz( pyz_VelXN * pyz_c4dx_xyz( xyzf_QMixN(:,:,:,f) ) ) &
        & - xyz_xqz( xqz_VelYN * xqz_c4dy_xyz( xyzf_QMixN(:,:,:,f) ) ) &
        & - xyz_xyr( xyr_VelZN * xyr_c4dz_xyz( xyzf_QMixN(:,:,:,f) ) ) &
        & - xyz_xyr( xyr_VelZN * xyrf_DQMixBZDz(:,:,:,f) )
      
      ! numerical diffusion term
      !
!      xyzf_QMixNDiff(:,:,:,f) =                                                             &
!        &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyzf_QMixB(:,:,:,f) ))))) &
!        &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyzf_QMixB(:,:,:,f) ))))) &
!        &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyzf_QMixB(:,:,:,f) ))))) 
      xyzf_QMixNDiff(:,:,:,f) =                         &
        &  - NuHh * xyz_dx4_xyz( xyzf_QMixB(:,:,:,f) )  &
        &  - NuHh * xyz_dy4_xyz( xyzf_QMixB(:,:,:,f) )  &
        &  - NuVh * xyz_dz4_xyz( xyzf_QMixB(:,:,:,f) )
  
    end do

    !------------------------------------------------------------------------
    ! X 方向の速度
    !

    ! advection
    !
    pyz_DVelXDtAdv =                                                      &
      & - pyz_VelXN * pyz_xyz( xyz_c4dx_pyz( pyz_VelXN ) )             &
      & - pyz_pqz( pqz_xqz( xqz_VelYN ) * pqz_c4dy_pyz( pyz_VelXN ) )  &
      & - pyz_pyr( pyr_xyr( xyr_VelZN ) * pyr_c4dz_pyz( pyz_VelXN ) )

    ! Numerical diffusion term 
    !
!    pyz_VelXnDiff =                                                             &
!      & - NuHm * ( pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz( pyz_VelXB ))))) &
!      & - NuHm * ( pyz_dy_pqz(pqz_dy_pyz(pyz_dy_pqz(pqz_dy_pyz( pyz_VelXB ))))) &
!      & - NuVm * ( pyz_dz_pyr(pyr_dz_pyz(pyz_dz_pyr(pyr_dz_pyz( pyz_VelXB )))))    
    pyz_VelXnDiff =                        &
      & - NuHm * pyz_dx4_pyz( pyz_VelXB )  &
      & - NuHm * pyz_dy4_pyz( pyz_VelXB )  &
      & - NuVm * pyz_dz4_pyz( pyz_VelXB )

    !------------------------------------------------------------------------
    ! Y 方向の速度
    !

    ! advection
    !
    xqz_DVelYDtAdv =                                                      &
      & - xqz_pqz( pqz_pyz( pyz_VelXN ) * pqz_c4dx_xqz( xqz_VelYN ) )  &
      & - xqz_VelYN * xqz_xyz( xyz_c4dy_xqz( xqz_VelYN ) )             &
      & - xqz_xqr( xqr_xyr( xyr_VelZN ) * xqr_c4dz_xqz( xqz_VelYN ) )

    ! Numerical diffusion term
    !
!    xqz_VelYnDiff =                                                             &
!      & - NuHm * ( xqz_dx_pqz(pqz_dx_xqz(xqz_dx_pqz(pqz_dx_xqz( xqz_VelYB ))))) &
!      & - NuHm * ( xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz( xqz_VelYB ))))) &
!      & - NuVm * ( xqz_dz_xqr(xqr_dz_xqz(xqz_dz_xqr(xqr_dz_xqz( xqz_VelYB )))))
    xqz_VelYnDiff =                        &
      & - NuHm * xqz_dx4_xqz( xqz_VelYB )  &
      & - NuHm * xqz_dy4_xqz( xqz_VelYB )  &
      & - NuVm * xqz_dz4_xqz( xqz_VelYB )

    !------------------------------------------------------------------------
    ! Z 方向の速度
    !

    ! Advection term
    !
    xyr_DVelZDtAdv =                                                       &
      & - xyr_pyr( pyr_pyz( pyz_VelXN ) * pyr_c4dx_xyr( xyr_VelZN ) )   &
      & - xyr_xqr( xqr_xqz( xqz_VelYN ) * xqr_c4dy_xyr( xyr_VelZN ) )   &
      & - xyr_VelZN * xyr_xyz( xyz_c4dz_xyr( xyr_VelZN ) )
       
    ! Numerical diffusion term
    !
!    xyr_VelZnDiff =                                                             &
!      & - NuHm * ( xyr_dx_pyr(pyr_dx_xyr(xyr_dx_pyr(pyr_dx_xyr( xyr_VelZB ))))) &
!      & - NuHm * ( xyr_dy_xqr(xqr_dy_xyr(xyr_dy_xqr(xqr_dy_xyr( xyr_VelZB ))))) & 
!      & - NuVm * ( xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr( xyr_VelZB )))))
    xyr_VelZnDiff =                        &
      & - NuHm * xyr_dx4_xyr( xyr_VelZB )  &
      & - NuHm * xyr_dy4_xyr( xyr_VelZB )  & 
      & - NuVm * xyr_dz4_xyr( xyr_VelZB )
    
  end subroutine advection_center4_std_main
  
end module advection_center4_std
