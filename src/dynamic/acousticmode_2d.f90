!= Module acousticmode_2d
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro)
! Version::   $Id: acousticmode_2d.f90,v 1.3 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module acousticmode_2d
  !
  ! 音波モードに関する計算ルーチンを束ねたモジュール
  !
  !   水平 explicit
  !   鉛直 implicit
  !

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

  character(*), parameter:: module_name = 'acousticmode_2d'
                                                  ! モジュールの名称.
                                                  ! Module name
  ! public 
  !
  public acousticmode_2d_init
  public acousticmode_2d_exp
  public acousticmode_2d_imp

contains

  subroutine acousticmode_2d_exp(                    &
    & pyz_VelXN, xyr_VelZN,                          & !(IN)
    & xyz_ExnerN,                                    & !(IN)
    & xyz_VelDivN,                                   & !(OUT)
    & pyz_PGrad, pyz_SWF                             & !(OUT)
    & )
    
    ! モジュール呼び出し
    !
    use gridset,   only : imin,            &! x 方向の配列の下限
      &                   imax,            &! x 方向の配列の上限
      &                   jmin,            &! y 方向の配列の下限
      &                   jmax,            &! y 方向の配列の上限
      &                   kmin,            &! z 方向の配列の下限
      &                   kmax,            &! z 方向の配列の上限
      &                   nx, ny, nz
    use axesset,   only : dx, dz            ! 格子間隔
    use constants, only : CpDry             ! 乾燥成分の比熱
    use basicset,  only : pyz_VPTempBZ 

    ! 暗黙の型宣言禁止
    !
    implicit none
    
    ! 変数の宣言
    !
    real(DP), intent(in)  :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: pyz_SWF(1:nx,1:ny,1:nz)
    
    integer               :: i, k
    integer, parameter    :: j = 1
    
    !------------------------------------------------------------------
    ! 速度の収束を計算
    !
    do k = kmin + 1, kmax 
        do i = imin + 1, imax 

          xyz_VelDivN(i,j,k) =         &
            & + (                      &
            &     pyz_VelXN(i,   j, k) &
            &   - pyz_VelXN(i-1, j, k) &
            &   ) / dx                 &
            & + (                      &
            &     xyr_VelZN(i, j, k)   &
            &   - xyr_VelZN(i, j, k-1) &
            &   ) / dz
        
        end do
    end do

    ! 値を確定させる
    !
    xyz_VelDivN(imin,:,:) = 1.0d10
    xyz_VelDivN(:,:,kmin) = 1.0d10

    !------------------------------------------------------------------
    ! X 方向
    
    do k = 1, nz
        do i = 1, nx
        
          ! 音波減衰項
          !            
          pyz_SWF(i,j,k) =                 &
            &   AlphaH                     &
            &   * (                        &
            &       xyz_VelDivN(i+1, j, k) &
            &     - xyz_VelDivN(i,   j, k) &
            &     ) / dx
          
          ! 圧力傾度力
          !
          pyz_PGrad(i,j,k) =                &
            & - CpDry * pyz_VPTempBZ(i,j,k) &
            &   * (                         &
            &       xyz_ExnerN(i+1, j, k)   &
            &     - xyz_ExnerN(i,   j, k)   &
            &     ) / dx
          
        end do
    end do
    
  end subroutine Acousticmode_2d_exp

!!!------------------------------------------------------------!!!

  subroutine acousticmode_2d_imp(                    &
    & pyz_VelXA, xyr_VelZN, xyz_VelDivN,             & !(IN)
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
    use axesset,  only : dx, dz            ! 格子間隔
    use constants,only : CpDry             ! 乾燥成分の比熱
    use timeset,  only : DelTimeShort
    use basicset, only : xyz_VPTempBZ,    &!基本場の仮温位
                            xyr_VPTempBZ      !基本場の仮温位

    !暗黙の型宣言禁止
    !
    implicit none
    
    !入出力変数
    real(DP), intent(in)   :: pyz_VelXA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 u [τ+Δτ]
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
    
    !作業変数定義
    real(DP)               :: D(1:nx,1:ny,1:nz)  
    real(DP)               :: E(1:nx,1:ny,0:nz)
    real(DP)               :: F(1:nx,1:ny,1:nz)
    real(DP)               :: F0(1:nx,1:ny,0:nz)  
    real(DP)               :: xyr_DExnerDz(1:nx,1:ny,0:nz)  
    real(DP)               :: xyr_DVelDivDz(1:nx,1:ny,0:nz)  
    real(DP)               :: dt            ! 短い時間格子間隔
    integer                :: INFO          ! 解のコンディション
    integer                :: i, k
    integer, parameter     :: j = 1
      
    real(DP)               :: TX(nz,nxny)    !解行列を転置したもの
    character(1),parameter :: TRANS = 'N'
    
    
    !---------------------------------------------------------------
    ! Initialize
    !
    dt = DelTimeShort
    
    !---------------------------------------------------------------
    !行列計算のための係数
    
    ! 共通して現れる微分を先に計算 
    !
    do k = 0, nz
        do i = 1, nx
        
          xyr_DExnerDz(i,j,k) =              &
            &  (                             &
            &   + xyz_ExnerN(i, j, k+1)      &
            &   - xyz_ExnerN(i, j, k)        &
            &   ) / dz                        
          
          xyr_DVelDivDz(i,j,k) =             &
            &  (                             &
            &  + xyz_VelDivN(i, j, k+1)      &
            &  - xyz_VelDivN(i, j, k)        &
            &  ) / dz  

      end do
    end do
    
    !  添字の範囲は, 1:nx, 1:ny, 0:nz
    !
    do k = 0, nz
        do i = 1, nx
        
          E(i,j,k) =                                    &
            & - ( 1.0d0 - beta ) * xyr_DExnerDZ(i,j,k)  &
            & + (                                       &
            &     AlphaV * xyr_DVelDivDz(i,j,k)         &
            &   + xyr_DVelZDtNl(i,j,k)                  &
            &   )                                       &
            &   / xyr_CpVPTempBZ(i,j,k)  
        
        end do
    end do
    
    ! 被微分関数
    !
    do k = 0, nz
        do i = 1, nx
        
          F0(i,j,k)  =                                  &
            & + xyr_DensVPTempBZ(i,j,k)                 &
            &   * (                                     &
            &     + xyr_VelZN(i,j,k)                    &
            &     - xyr_CpVPTempBZ(i,j,k)               &
            &       * ( 1.0d0 - beta )                  &
            &       * xyr_DExnerDZ(i,j,k) * dt          &
            &     + AlphaV * xyr_DVelDivDz(i,j,k) * dt  &
            &     + xyr_DVelZDtNl(i,j,k)  * dt          &
            &    )
        
        end do
    end do
    
    !行列計算のための係数
    !  添字の範囲は, 1:nx, 1:ny, 1:nz
    !
    do k = 1, nz
        do i = 1, nx
        
          F(i,j,k) =                                  &
            & - beta * dt                             &
            &   * xyz_VelSoundSQBZ(i,j,k)             &
            &   / xyz_CpDensVPTempSQBZ(i,j,k)         &
            &   * (                                   &
            &       F0(i, j, k)                       &
            &     - F0(i, j, k-1)                     &
            &     ) / dz                              &
            & + xyz_DExnerDtNl(i,j,k) * dt            &
            & + xyz_DExnerDtNs(i,j,k) * dt
        
        end do
    end do
    
    !行列計算のための係数
    !  添字の範囲は, 1:nx, 1:ny, 1:nz
    !
    do k = 1, nz 
        do i = 1, nx
        
          D(i,j,k) =                                               &
            & + xyz_ExnerN(i,j,k)                                  &
            & - (1.0d0 - beta) * dt                                &
            &   * xyz_VelSoundSQBZ(i,j,k)                          &
            &   / xyz_CpDensVPTempSQBZ(i,j,k)                      &
            &   * (                                                &
            &       xyr_DensVPTempBZ(i,j,k)   * xyr_VelZN(i,j,k)   &
            &     - xyr_DensVPTempBZ(i,j,k-1) * xyr_VelZN(i,j,k-1) &
            &     ) / dz                                           &
            & - xyz_VelSoundSQBZ(i,j,k) * dt                       &
            &   / (CpDry * xyz_VPTempBZ(i,j,k))                    &
            &   * (                                                &
            &     + (                                              &
            &         pyz_VelXA(i,   j, k)                         &
            &       - pyz_VelXA(i-1, j, k)                         &
            &       ) / dx                                         &
            &     )                                                &
            & + F(i,j,k)
        
        end do
    end do
    
    ! 行列計算のための係数
    !
    do i = 1, nx
        
        D(i,j,1) =                                     &
          & + D(i,j,1)                                 &
          & - beta * (dt * dt)                         &
          &   * xyz_VelSoundSQBZ(i,j,1)                &
          &   / xyz_CpDensVPTempSQBZ(i,j,1)            &
          &   * xyr_CpDensVPTempSQBZ(i,j,0)            &
          &   * E(i,j,0)                               &
          &   / dz
      
        D(i,j,nz) =                                    &
          & + D(i,j,nz)                                &
          & + beta * (dt * dt)                         &
          &   * xyz_VelSoundSQBZ(i,j,nz)               &
          &   / xyz_CpDensVPTempSQBZ(i,j,nz)           &
          &   * xyr_CpDensVPTempSQBZ(i,j,nz)           &
          &   * E(i,j,nz)                              &
          &   / dz
      
      end do
    
    !-----------------------------------------------------------
    !連立一次方程式の解を求める
    
    ! LAPACK の仕様に合わせて変形 
    !
    do k = 1, nz
        do i = 1, nx
          TX(k,i) = D(i,j,k)
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
        do i = 1, nx
          xyz_ExnerA(i,j,k) = TX(k,i)
        end do
    end do
    
    !------------------------------------------------------------
    ! 鉛直速度
    !
    ! 鉛直速度の k = nz の値は境界条件によって決まるので, k 方向の
    ! ループは 1 ~ nz-1 とする.
    ! 
    
    do k = 1, nz-1
        do i = 1, nx
          
          ! 音波減衰項
          !            
          xyr_SWF(i,j,k) =                 &
            & + AlphaV                     &
            &   * (                        &
            &       xyz_VelDivN(i,j,k+1)   &
            &     - xyz_VelDivN(i,j,k)     &
            &     ) / dz
          
          ! 圧力傾度力
          !
          xyr_PGrad(i,j,k) =                &
            & - CpDry * xyr_VPTempBZ(i,j,k) &
            &   * (                         &
            &       beta                    &
            &       * (                     &
            &           xyz_ExnerA(i,j,k+1) &
            &         - xyz_ExnerA(i,j,k)   &
            &         )                     &
            &     + (1.0d0 - beta)          &
            &       * (                     &
            &           xyz_ExnerN(i,j,k+1) &
            &         - xyz_ExnerN(i,j,k)   &
            &         )                     &
            &     ) / dz

        end do
    end do

    ! 穴埋め
    !   サブルーチンの戻り値とするためには, 配列の値を確定させる必要があるため. 
    ! 
    xyz_ExnerA(imin:0,:,:) = 0.0d0
    xyz_ExnerA(nx+1:imax,:,:) = 0.0d0
    xyz_ExnerA(:,:,kmin:0) = 0.0d0
    xyz_ExnerA(:,:,nz+1:kmax) = 0.0d0
    xyr_SWF(:,:,nz)   = 0.0d0
    xyr_PGrad(:,:,nz) = 0.0d0

  end subroutine Acousticmode_2d_imp
  
!!!--------------------------------------------------------------------!!!

  subroutine acousticmode_2d_init( AlphaSound )
    !
    ! 音波モードの計算モジュールの初期化
    ! エクスナー関数を陰解法で解く際に必要となる, 係数行列の要素を決め, LU 分解を行う. 
    !
    
    ! モジュール読み込み
    !
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x 方向の物理領域の上限
      &                    ny,              &! x 方向の物理領域の上限
      &                    nz                ! y 方向の物理領域の上限
    use constants,  only : CpDry             ! 乾燥成分の比熱
    use timeset,    only : DelTimeShort
    use axesset,    only : dx, dz            ! 格子間隔
    use basicset,   only : xyz_VelSoundBZ,  &!基本場の音速 
      &                    xyz_DensBZ,      &!基本場の密度
      &                    xyz_VPTempBZ      !基本場の温位
    use average,    only : xyr_xyz 

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in) :: AlphaSound
    real(DP)             :: r_CpDensVPTempSQBZ(0:nz)
    real(DP)             :: z_VelSoundSQBZ(1:nz) 
    real(DP)             :: z_CpDensVPTempSQBZ(1:nz) 
    real(DP)             :: dt      ! 短い時間格子
    integer              :: INFO    !解のコンディションチェック
    integer              :: k

    !----------------------------------------------------------------
    ! 初期化

    ! 音波減衰項の減衰係数を決める
    ! 
    ! 気象庁予報課報告別冊 49 p53 に従い, 水平と鉛直とを分けて考える. 
    !
!    AlphaH = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort
!    AlphaV = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort
    AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
    AlphaV = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort

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
    z_VelSoundSQBZ(1:nz)     = xyz_VelSoundSQBZ(nx,ny,1:nz)
    z_CpDensVPTempSQBZ(1:nz) = xyz_CpDensVPTempSQBZ(nx,ny,1:nz)
    r_CpDensVPTempSQBZ(0:nz) = xyr_CpDensVPTempSQBZ(nx,ny,0:nz)

    do k = 2, nz-1
      A(k) =                                &
        & + 1.0d0                           &
        & + ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k)            &
        &    / z_CpDensVPTempSQBZ(k)        &
        &    * (                            &
        &         r_CpDensVPTempSQBZ(k)     &
        &       + r_CpDensVPTempSQBZ(k-1)   &
        &       )                           &
        &    * ( dt * dt )                  &
        &    / ( dz * dz )
    end do

    A(1) =                                  &
      & + 1.0d0                             &
      & + ( beta * beta )                   &
      &   * z_VelSoundSQBZ(1)               &
      &   / z_CpDensVPTempSQBZ(1)           &
      &   * r_CpDensVPTempSQBZ(1)           &
      &   * ( dt * dt )                     &
      &   / ( dz * dz ) 

    A(nz) =                                 &
      & + 1.0d0                             &
      & + ( beta * beta )                   &
      &   * z_VelSoundSQBZ(nz)              &
      &   / z_CpDensVPTempSQBZ(nz)          &
      &   * r_CpDensVPTempSQBZ(nz-1)        &
      &   * ( dt * dt )                     &
      &   / ( dz * dz )  

    do k = 2, nz
      B(k) =                                &
        & - ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k-1)          &
        &    / z_CpDensVPTempSQBZ(k-1)      &
        &   * r_CpDensVPTempSQBZ(k-1)       &
        &   * ( dt * dt )                   &
        &   / ( dz * dz )
    end do
    
    do k = 1, nz-1
      C(k) =                                &
        & - ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k+1)          &
        &    / z_CpDensVPTempSQBZ(k+1)      &
        &   * r_CpDensVPTempSQBZ(k)         &
        &   * (dt * dt )                    &
        &   / ( dz * dz )
    end do

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

  end subroutine Acousticmode_2d_init
  
end module Acousticmode_2d
