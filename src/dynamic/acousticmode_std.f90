!= Module acousticmode_std
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro)
! Version::   $Id: acousticmode_std.f90,v 1.3 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module acousticmode_std
  !
  ! 音波モードに関する計算ルーチンを束ねたモジュール
  !
  !   水平 explicit
  !   鉛直 implicit
  !
  ! 微分平均モジュールを用いて計算を行う

  !モジュール読み込み
  !
  use dc_types,   only : DP

  !暗黙の型宣言禁止
  !
  implicit none

  !属性の指定
  !
  private

  ! 変数の定義
  !
  real(DP), save :: beta  = 0.5d0         !クランクニコルソン法なら 0.5
                                          !完全陰解法なら 1

  real(DP), save :: AlphaH = 0.0d0        !音波減衰項の減衰係数 (水平方向)
  real(DP), save :: AlphaV = 0.0d0        !音波減衰項の減衰係数 (鉛直方向)

  real(DP), allocatable, save :: A(:)     !係数行列の対角成分
  real(DP), allocatable, save :: B(:)     !係数行列の上三角部分
  real(DP), allocatable, save :: C(:)     !係数行列の下三角部分
  real(DP), allocatable, save :: AL1(:)   !LU 分解の結果 L (1 次元配列)
  integer,  allocatable, save :: IP(:)    !部分ピボット交換の情報を格納

  real(DP), allocatable, save :: xyr_CpVPTempBZ(:,:,:)       !係数行列の計算に用いる配列
  real(DP), allocatable, save :: xyr_CpDensVPTempSQBZ(:,:,:) !係数行列の計算に用いる配列
  real(DP), allocatable, save :: xyr_DensVPTempBZ(:,:,:)     !係数行列の計算に用いる配列
  real(DP), allocatable, save :: xyz_VelSoundSQBZ(:,:,:)     !係数行列の計算に用いる配列
  real(DP), allocatable, save :: xyz_CpDensVPTempSQBZ(:,:,:) !係数行列の計算に用いる配列

  character(*), parameter:: module_name = 'acousticmode_std'
                                                  ! モジュールの名称.
                                                  ! Module name
  ! public 
  !
  public acousticmode_std_init
  public acousticmode_std_exp
  public acousticmode_std_imp
  
contains

  subroutine acousticmode_std_exp(                   &
    & pyz_VelXN, xqz_VelYN, xyr_VelZN, xyz_ExnerN,   & !(IN)
    & xyz_VelDivN,                                   & !(OUT)
    & pyz_PGrad, pyz_SWF,                            & !(OUT)
    & xqz_PGrad, xqz_SWF                             & !(OUT)
    & )
    
    ! モジュール呼び出し
    !
    use dc_types,  only : DP
    use gridset,   only : imin,            &! x 方向の配列の下限
      &                   imax,            &! x 方向の配列の上限
      &                   jmin,            &! y 方向の配列の下限
      &                   jmax,            &! y 方向の配列の上限
      &                   kmin,            &! z 方向の配列の下限
      &                   kmax,            &! z 方向の配列の上限
      &                   nx, ny, nz
    use constants, only : CpDry             ! 乾燥成分の比熱
    use basicset,  only : pyz_VPTempBZ,    &! 基本場の温位
                          xqz_VPTempBZ      ! 基本場の温位
    use differentiate_center2,             &
      &            only : xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                   pyz_dx_xyz, xqz_dy_xyz

    ! 暗黙の型宣言禁止
    !
    implicit none
    
    ! 変数の宣言
    !
    real(DP), intent(in)  :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xqz_VelYN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: xqz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: pyz_SWF(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: xqz_SWF(1:nx,1:ny,1:nz)

    real(DP)              :: aaa_tmp(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------------
    ! 速度の収束を計算
    !
    xyz_VelDivN =                 &
      &   xyz_dx_pyz( pyz_VelXN ) &
      & + xyz_dy_xqz( xqz_VelYN ) &
      & + xyz_dz_xyr( xyr_VelZN )

    !------------------------------------------------------------------
    ! X 方向

!    aaa_tmp   = + AlphaH * pyz_dx_xyz( xyz_VelDivN )     
    aaa_tmp   =   pyz_dx_xyz( AlphaH * xyz_VelDivN )     
    pyz_SWF   =   aaa_tmp(1:nx,1:ny,1:nz)
    
    aaa_tmp   = - CpDry * pyz_VPTempBZ * pyz_dx_xyz( xyz_ExnerN ) 
    pyz_PGrad =   aaa_tmp(1:nx,1:ny,1:nz)

    !------------------------------------------------------------------
    ! Y 方向
    
    aaa_tmp   =   AlphaH * xqz_dy_xyz( xyz_VelDivN ) 
!    aaa_tmp   =   xqz_dy_xyz( AlphaH * xyz_VelDivN ) 
    xqz_SWF   =   aaa_tmp(1:nx,1:ny,1:nz)

    aaa_tmp   = - CpDry * xqz_VPTempBZ * xqz_dy_xyz( xyz_ExnerN ) 
    xqz_PGrad =   aaa_tmp(1:nx,1:ny,1:nz)

  end subroutine Acousticmode_std_exp

!!!------------------------------------------------------------!!!

  subroutine acousticmode_std_imp(                   &
    & pyz_VelXA, xqz_VelYA, xyr_VelZN, xyz_VelDivN,  & !(IN)
    & xyz_ExnerN,                                    & !(IN)
    & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs, & !(IN)
    & xyz_ExnerA,                                    & !(OUT)
    & xyr_PGrad, xyr_SWF                             & !(OUT)
    & )
    !
    ! 陰解法を用いたエクスナー関数/鉛直速度の計算. 
    !
    
    ! モジュールの読み込み
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x 方向の配列の下限
      &                  imax,            &! x 方向の配列の上限
      &                  jmin,            &! y 方向の配列の下限
      &                  jmax,            &! y 方向の配列の上限
      &                  kmin,            &! z 方向の配列の下限
      &                  kmax,            &! z 方向の配列の上限
      &                  nx, ny, nz,      &! 物理領域の大きさ
      &                  nxny              ! 物理領域の大きさ (nx * ny)
    use constants,only : CpDry             ! 乾燥成分の比熱
!    use timeset,  only : DelTimeShort, TimeN
    use timeset,  only : DelTimeShort
    use basicset, only : xyz_VPTempBZ,    &!基本場の仮温位
                         xyr_VPTempBZ      !基本場の仮温位
    use axesset,  only : dz
    use differentiate_center2, &
      &           only : xyr_dz_xyz, xyz_dz_xyr, xyz_dx_pyz, xyz_dy_xqz
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    !
    real(DP), intent(in)   :: pyz_VelXA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 u [τ+Δτ]
    real(DP), intent(in)   :: xqz_VelYA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 v [τ+Δτ]
    real(DP), intent(in)   :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 w [τ]
    real(DP), intent(in)   :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
                                                           !\Div \Dvect{u}
    real(DP), intent(in)   :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !無次元圧力
    real(DP), intent(in)   :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z 方向の外力項[t]
    real(DP), intent(in)   :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !エクスナー関数の外力項[t]
    real(DP), intent(in)   :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !エクスナー関数の外力項[t]
    real(DP), intent(out)  :: xyz_ExnerA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !無次元圧力[τ+Δτ]
    real(DP), intent(out)  :: xyr_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out)  :: xyr_SWF(1:nx,1:ny,1:nz)
    
    ! 作業変数定義
    !
    real(DP)               :: D(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: E(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: F(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: F0(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: xyr_DExnerDz(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: xyr_DVelDivDz(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: aaa_tmp(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: dt            ! 短い時間格子間隔
    integer                :: INFO          ! 解のコンディション
    integer                :: i, j, k
      
    real(DP)               :: TX(nz,nxny)    !解行列を転置したもの
    character(1),parameter :: TRANS = 'N'
    

    !---------------------------------------------------------------
    ! Initialize
    !
    xyz_ExnerA = 0.0d0
    dt = DelTimeShort

    !---------------------------------------------------------------
    !行列計算のための係数
    
    ! 共通して現れる微分を先に計算 
    !
    xyr_DExnerDZ =  xyr_dz_xyz( xyz_ExnerN ) 
    
    xyr_DVelDivDZ = xyr_dz_xyz( xyz_VelDivN )
    
    E =                                      &
      & - ( 1.0d0 - beta ) * xyr_DExnerDZ    &
      & + (                                  &
      &     AlphaV * xyr_DVelDivDZ           &
      &   + xyr_DVelZDtNl                    &
      &   )                                  &
      &   / xyr_CpVPTempBZ
    
    F0  =                                                            &
      & + xyr_DensVPTempBZ                                           &
      &   * (                                                        &
      &     + xyr_VelZN                                              &  
      &     - xyr_CpVPTempBZ * ( 1.0d0 - beta) * xyr_DExnerDZ * dt   &
      &     + AlphaV * xyr_DVelDivDZ  * dt                           &
      &     + xyr_DVelZDtNl  * dt                                    &
      &    )
    
    F =                                            &
      & - beta * dt                                & 
      &   * xyz_VelSoundSQBZ                       &  
      &   / xyz_CpDensVPTempSQBZ                   &  
      &   * xyz_dz_xyr( F0 )                       &
      & + ( xyz_DExnerDtNl + xyz_DExnerDtNs ) * dt
    
    D =                                                           &
      & + xyz_ExnerN                                              &
      & - (1.0d0 - beta) * dt                                     &
      &   * xyz_VelSoundSQBZ                                      &  
      &   / xyz_CpDensVPTempSQBZ                                  &  
      &   * xyz_dz_xyr( xyr_DensVPTempBZ * xyr_VelZN )            &
      & - xyz_VelSoundSQBZ * dt                                   &
      &   / (CpDry * xyz_VPTempBZ )                               &
      &   * ( xyz_dx_pyz( pyz_VelXA ) + xyz_dy_xqz( xqz_VelYA ) ) &
      & + F
    
    D(:,:,1) =                                     &
      & + D(:,:,1)                                 &
      & - beta * (dt * dt)                         &
      &   * xyz_VelSoundSQBZ(:,:,1)                &  
      &   / xyz_CpDensVPTempSQBZ(:,:,1)            &  
      &   * xyr_CpDensVPTempSQBZ(:,:,0)            &
      &   * E(:,:,0)                               &
      &   / dz
    
    D(:,:,nz) =                                    &
      & + D(:,:,nz)                                &
      & + beta * (dt * dt)                         &
      &   * xyz_VelSoundSQBZ(:,:,nz)               &  
      &   / xyz_CpDensVPTempSQBZ(:,:,nz)           &  
      &   * xyr_CpDensVPTempSQBZ(:,:,nz)           &
      &   * E(:,:,nz)                              &
      &   / dz

!    call HistoryAutoPut(TimeN, 'D', D(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'E', E(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'F', F(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENs', xyz_DExnerDtNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENl', xyz_DExnerDtNl(1:nx,1:ny,1:nz))
    
    !-----------------------------------------------------------
    !連立一次方程式の解を求める
    
    ! LAPACK の仕様に合わせて変形 
    !
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          TX(k, i + nx * (j - 1)) = D(i,j,k)
        end do
      end do
    end do
    
    !解行列の計算. LAPACK を使用. 
    !
    call DGTTRS(TRANS, nz, nxny, C, A, B, AL1, IP, TX, nz, INFO)
    
!    !解のコンディションをチェック. 
!    !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyz_ExnerA(i,j,k) = TX(k, i + nx * (j - 1 ))
        end do
      end do
    end do
    

    !------------------------------------------------------------
    ! 鉛直速度
    !
!    aaa_tmp =  AlphaV * xyr_dz_xyz( xyz_VelDivN ) 
    aaa_tmp =  xyr_dz_xyz( AlphaV * xyz_VelDivN ) 
    xyr_SWF =  aaa_tmp(1:nx,1:ny,1:nz)
    
    aaa_tmp =                                               &
      & - CpDry * xyr_VPTempBZ                              &
      &   * (                                               &
      &         beta           * xyr_dz_xyz( xyz_ExnerA )   &
      &       + (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerN )   &
      &     )                                                
    xyr_PGrad =  aaa_tmp(1:nx,1:ny,1:nz)
    
  end subroutine Acousticmode_std_imp
  

!!!--------------------------------------------------------------------!!!
  subroutine acousticmode_std_init( AlphaSound )
    !
    !エクスナー関数を陰解法で解く際に必要となる, 係数行列の要素を決め, 
    !LU 分解を行う. 
    !
    
    ! モジュール読み込み
    !
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : FlagCalc3D,      &
      &                    imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x 方向の物理領域の上限
      &                    ny,              &! x 方向の物理領域の上限
      &                    nz                ! y 方向の物理領域の上限
    use constants,  only : CpDry           ! 乾燥成分の比熱
    use timeset,    only : DelTimeShort
    use axesset,    only : dx, dy, dz        ! 格子間隔
    use basicset,   only : xyz_VelSoundBZ,  &!基本場の音速 
      &                    xyz_DensBZ,      &!基本場の密度
      &                    xyz_VPTempBZ      !基本場の温位
    use average,    only : xyr_xyz
    
    !暗黙の型宣言禁止
    !
    implicit none

    !変数の定義
    !
    real(DP), intent(in) :: AlphaSound
    real(DP)             :: r_CpDensVPTempSQBZ(kmin:kmax)
    real(DP)             :: z_VelSoundSQBZ(kmin:kmax) 
    real(DP)             :: z_CpDensVPTempSQBZ(kmin:kmax) 
    real(DP)             :: dt      ! 短い時間格子
    integer              :: INFO    !解のコンディションチェック

    !----------------------------------------------------------------
    ! 初期化

    ! 音波減衰項の減衰係数を決める
    ! 
    ! 気象庁予報課報告別冊 49 p53 に従い, 水平と鉛直とを分けて考える. 
    !
!    AlphaH = AlphaSound * ( Min( dx * dx, dy * dy ) ) / DelTimeShort
!!    AlphaH = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz ) ) / DelTimeShort
!    AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz ) ) / DelTimeShort
    if ( FlagCalc3D ) then 
      AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
    else
      AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
    end if


    !-------------------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", module_name, "AlphaH = %f", d=(/AlphaH/) )
    call MessageNotify( "M", module_name, "AlphaV = %f", d=(/AlphaV/) )

    ! 変数名が長すぎたので, 名前を置き換える
    !
    dt = DelTimeShort

    ! 配列の割り付け
    !
    allocate( A(1:nz) )
    allocate( B(2:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_VelSoundSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CpDensVPTempSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_CpDensVPTempSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DensVPTempBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_CpVPTempBZ(imin:imax,jmin:jmax,kmin:kmax) )
   
    !----------------------------------------------------------------
    ! 係数行列および共通して利用される配列の値を決める
    !
    xyz_VelSoundSQBZ     = xyz_VelSoundBZ * xyz_VelSoundBZ 
    xyz_CpDensVPTempSQBZ = CpDry * xyz_DensBZ * xyz_VPTempBZ * xyz_VPTempBZ
    xyr_CpDensVPTempSQBZ = xyr_xyz( xyz_CpDensVPTempSQBZ )
    xyr_DensVPTempBZ     = xyr_xyz( xyz_DensBZ * xyz_VPTempBZ )
    xyr_CpVPTempBZ       = xyr_xyz( CpDry * xyz_VPTempBZ )

    ! 鉛直方向のみの変数については, 基本場を使うので, 
    ! nx, ny の値で代表させることとした. 
    !
    z_VelSoundSQBZ        = xyz_VelSoundSQBZ(nx,ny,:)
    z_CpDensVPTempSQBZ(:) = xyz_CpDensVPTempSQBZ(nx,ny,:)
    r_CpDensVPTempSQBZ(:) = xyr_CpDensVPTempSQBZ(nx,ny,:)

    !
    !
    A(2:nz-1) =                               &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(2:nz-1)           &
      &    / z_CpDensVPTempSQBZ(2:nz-1)       &
      &    * (                                &
      &         r_CpDensVPTempSQBZ(2:nz-1)    &
      &       + r_CpDensVPTempSQBZ(1:nz-2)    &
      &       )                               &
      &    * ( dt * dt )                      &
      &    / ( dz * dz )

    A(1) =                                    &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(1)                &
      &    / z_CpDensVPTempSQBZ(1)            &
      &   * r_CpDensVPTempSQBZ(1)             &
      &   * ( dt * dt )                       &
      &   / ( dz * dz ) 

    A(nz) =                                   &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(nz)               &
      &    / z_CpDensVPTempSQBZ(nz)           &
      &   * r_CpDensVPTempSQBZ(nz-1)          &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )  

    B(2:nz) =                                 &
      & - ( beta * beta )                     &
      &    * z_VelSoundSQBZ(1:nz-1)           &
      &    / z_CpDensVPTempSQBZ(1:nz-1)       &
      &   * r_CpDensVPTempSQBZ(1:nz-1)        &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )
    
    C(1:nz-1) =                               &
      & - ( beta * beta )                     &
      &    * z_VelSoundSQBZ(2:nz)             &
      &    / z_CpDensVPTempSQBZ(2:nz)         &
      &   * r_CpDensVPTempSQBZ(1:nz-1)        &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )

    !----------------------------------------------------------------
    ! 係数行列を LU 分解
    !

    ! 配列の割り当て
    !
    allocate( AL1(1:nz-2), IP(1:nz) )

    ! 解行列の計算. LAPACK を使用. 
    !
    call DGTTRF(nz, C, A, B, AL1, IP, INFO)
    
    ! 解のコンディションをチェック. 
    !
    if (INFO /= 0) then
      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
      stop
    end if

  end subroutine Acousticmode_std_init

end module Acousticmode_std
