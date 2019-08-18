!---------------------------------------------------------------------
!     Copyright (C) GFD Dennou Club, 2004, 2005. All rights reserved.
!---------------------------------------------------------------------
!= Module DisturbEnv
!
!   * Developer: SUGIYAMA Ko-ichiro, ODAKA Masatsugu
!   * Version: $Id: initialdata_disturb.f90,v 1.17 2014/07/08 00:59:08 sugiyama Exp $ 
!   * Tag Name: $Name:  $
!   * Change History: 


module initialdata_disturb
  !
  !擾乱のデフォルト値を与えるためのルーチン. 
  !
  
  !暗黙の型宣言禁止
  implicit none

  !公開要素
  public initialdata_disturb_random
  public initialdata_disturb_gaussXZ
  public initialdata_disturb_gaussXY
  public initialdata_disturb_gaussYZ
  public initialdata_disturb_gaussXYZ
  public initialdata_disturb_cosXZ
  public initialdata_disturb_cosXY
  public initialdata_disturb_cosYZ
  public initialdata_disturb_cosXYZ
  public initialdata_disturb_coneXZ
  public initialdata_disturb_coneXY
  public initialdata_disturb_coneYZ
  public initialdata_disturb_dryreg
  public initialdata_disturb_square
  public initialdata_disturb_moist
  public initialdata_disturb_tanh
  public initialdata_disturb_tanh_sin
  public initialdata_disturb_circleXZ

contains

  subroutine initialdata_disturb_random( DelMax, Zpos, xyz_Var )
    !
    ! 初期擾乱を乱数で与える
    !

    !モジュール読み込み
    use mpi_wrapper,only: myrank
    use dc_types,   only: DP
    use axesset,    only: z_Z                   ! Z 座標軸(スカラー格子点)
    use gridset,    only: nx, imin, imax,      &! 配列サイズ (X 方向)
      &                   ny, jmin, jmax,      &! 配列サイズ (Y 方向)
      &                   kmin, kmax            ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Zpos
    real(DP), intent(out) :: xyz_Var(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: Random(imin:imax, jmin:jmax)
    real(DP)              :: mean
    integer               :: i, j, k, kpos
    integer, allocatable  :: seed(:)
    integer               :: seedsize
    
    ! 初期化
    xyz_Var = 0.0d0
    kpos    = 1

    ! 0.0--1.0 の擬似乱数発生
    !  mpi の場合に, 各 CPU の持つ乱数が異なるよう調整している.  
    !
    call random_seed(size=seedsize)  !初期値のサイズを取得
    allocate(seed(seedsize))         !配列割り当て

    ! システムクロックを乱数の種として利用
    do i = 1, seedsize
      call system_clock(count=seed(i)) !時間取得
    end do
    
    ! ノード毎にシードを変える.
    seed = seed * ( myrank + 2) * 10 
    
    ! 乱数の種を与える
    call random_seed(put=seed(:))
    deallocate(seed)

    ! 乱数生成. 値の範囲は 0-1 となる. 
    do j = 1, ny
      do i = 1, nx
        call random_number(random(i,j))
      end do
    end do
    
    ! 指定された高度の配列添字を用意   
    do k = kmin, kmax
      if ( z_Z(k) >= Zpos ) then 
        kpos = k
        exit
      end if
    end do

    ! 擾乱が全体としてはゼロとなるように調整. 平均からの差にする. 
    mean = sum( Random(1:nx, 1:ny) ) / real((nx * ny),8)
    do j = 1, ny
      do i = 1, nx
        xyz_Var(i, j, kpos) = DelMax * ( Random(i,j) - mean )
      end do
    end do
    
  end subroutine initialdata_disturb_random

!!!
!!! Klemp and Wilhelmson (1978) の初期値
!!! 

!  subroutine initialdata_disturb_kw1978(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
!    !
!    ! Klemp and Wilhelmson (1978) の初期値
!
!    implicit none
!
!    real(DP), intent(in)   :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
!    real(DP), intent(out)  :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
!    real(DP)               :: xyz_beta(imin:imax, jmin:jmax, kmin:kmax)
!    integer                :: i, j, k
!
!    do k = kmin, kmax
!      do j = jmin, jmax
!        do i = imin, imax
!          beta(i,j,k) =                               &
!            & (                                       &
!            &      ( ( x_X(i) - Xc ) / Xr ) ** 2.0d0  &
!            &    + ( ( y_Y(j) - Yc ) / Yr ) ** 2.0d0  &
!            &    + ( ( z_Z(k) - Zc ) / Zr ) ** 2.0d0  &
!            &  ) ** 5.0d-1
!        end do
!      end do
!    end do
!    
!    where ( beta < 1.0d0 )
!      xyz_Var = DelMax * ( dcos( Pi * 5.0d-1 * beta ) ** 2.0d0 )
!    end where
!
!  end subroutine initialdata_disturb_kw1978



!!!
!!! tanh
!!! 

  subroutine initialdata_disturb_tanh(VarMean, VarDel, Zc, Zr, aaz_Var, aaz_VarBZ)
    !
    ! tanh 型のシア
    !   A(z) = A0 + A1 \tanh( (z - Zc) / h )

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)            :: VarMean, VarDel, Zc, Zr
    real(DP), intent(in), optional  :: aaz_VarBZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(out)           :: aaz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer                         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          aaz_Var(i,j,k) = VarMean + VarDel * tanh( (z_Z(k) - Zc) / Zr ) 
        end do
      end do
    end do

    if ( present( aaz_VarBZ) ) then 
      aaz_Var = aaz_Var - aaz_VarBZ
    end if

  end subroutine initialdata_disturb_tanh


  subroutine initialdata_disturb_tanh_sin(VarMean, VarDel, Zc, Zr, aaz_Var, aaz_VarBZ)
    !
    ! tanh 型のシア
    !   A(z) = A0 + A1 \tanh( (z - Zc) / h )

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,  XMax,    &! X 座標軸(スカラー格子点)
      &                  z_Z,  ZMax      ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)            :: VarMean, VarDel, Zc, Zr
    real(DP), intent(in), optional  :: aaz_VarBZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(out)           :: aaz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter             :: Pi = 3.1415926535897932385d0 
    integer                         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          aaz_Var(i,j,k) =                                  &
            &   VarMean                                     &
            & + VarDel * tanh( (z_Z(k) - Zc) / Zr )         &
            & + VarDel * sin( x_X(i) * 2.0d0 * Pi / XMax )  &
            &          * sin( z_Z(k) * 2.0d0 * Pi / ZMax )
        end do
      end do
    end do

    if ( present( aaz_VarBZ) ) then 
      aaz_Var = aaz_Var - aaz_VarBZ
    end if

  end subroutine initialdata_disturb_tanh_sin


!!!
!!! 円錐形
!!! 

  subroutine initialdata_disturb_ConeXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! 円錐型の初期値

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  y_Y             ! Y 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((x_X(i) - Xc) / Xr )**2.0d0   &
            &                    +  ((y_Y(j) - Yc) / Yr )**2.0d0) )
          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeXY


  subroutine initialdata_disturb_ConeXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! 円錐型の初期値
    
    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((x_X(i) - Xc) / Xr )**2.0d0   &
            &                    +  ((z_Z(k) - Zc) / Zr )**2.0d0) )

          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeXZ 


  subroutine initialdata_disturb_ConeYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! 円錐型の初期値

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: y_Y,           &! Y 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((y_Y(i) - Yc) / Yr )**2.0d0   &
            &                    +  ((z_Z(k) - Zc) / Zr )**2.0d0) )

          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeYZ

  
!!!
!!! 一定
!!!  
  
  subroutine initialdata_disturb_circleXZ(DelMax, Xc, Xr, Zc, xyz_Var)
    !
    ! 円柱型の初期値
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Zc, Xr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k

    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          if ( ( (x_X(i) - Xc) ** 2.0d0 + (z_Z(k) - Zc) ** 2.0d0 ) < Xr*Xr ) then 
            xyz_Var(i,j,k) = DelMax
          else
            xyz_Var(i,j,k) = 0.0d0
          end if
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_circleXZ


!!!
!!! ガウシアン
!!!  
  
  subroutine initialdata_disturb_gaussXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! ガウシアンな初期値
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k

    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1 ) 
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussXZ  

  subroutine initialdata_disturb_gaussXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! ガウシアンな初期値
    !

    !モジュール読み込み
    use dc_types, only: DP
    use axesset,  only: x_X,           &! X 座標軸(スカラー格子点)
      &                 y_Y             ! Y 座標軸(スカラー格子点)
    use gridset,  only: imin, imax,    &! 配列サイズ (X 方向)
      &                 jmin, jmax,    &! 配列サイズ (Y 方向)
      &                 kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1 )
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussXY


  subroutine initialdata_disturb_gaussYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ガウシアンな初期値
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: y_Y,           &! Y 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1 )
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussYZ


  subroutine initialdata_disturb_gaussXYZ(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ガウシアンな初期値
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  y_Y,           &! Y 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1   &
            &                - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1 ) 
        end do
      end do
    end do
    
!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where

  end subroutine initialdata_disturb_gaussXYZ

!!!
!!! 
!!!
  subroutine initialdata_disturb_cosXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! 重力流の計算に利用する擾乱
    ! A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: xyz_X,         &! X 座標軸(スカラー格子点)
      &                  xyz_Z           ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 

    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXZ
  
  
  subroutine initialdata_disturb_cosXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! 重力流の計算に利用する擾乱
    ! A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: xyz_X,         &! X 座標軸(スカラー格子点)
      &                  xyz_Y           ! Y 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列の X 方向の上限
      &                  jmin, jmax,    &! 配列の Y 方向の上限
      &                  kmin, kmax      ! 配列の Z 方向の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 
    
    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXY
  

  subroutine initialdata_disturb_cosYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! 重力流の計算に利用する擾乱
    ! A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: xyz_Y,        &! X 座標軸(スカラー格子点)
      &                  xyz_Z          ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,   &! 配列サイズ (X 方向)
      &                  jmin, jmax,   &! 配列サイズ (Y 方向)
      &                  kmin, kmax     ! 配列サイズ (Z 方向)
    
    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 

    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )

    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosYZ
  
  
  subroutine initialdata_disturb_cosXYZ(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! 重力流の計算に利用する擾乱
    ! A [cos(πL) + 1]*0.5 ( L < 1.0 ) or 0.0

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: xyz_X,        &! X 座標軸(スカラー格子点)
      &                  xyz_Y,        &! X 座標軸(スカラー格子点)
      &                  xyz_Z          ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,   &! 配列サイズ (X 方向)
      &                  jmin, jmax,   &! 配列サイズ (Y 方向)
      &                  kmin, kmax     ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 
    
    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXYZ

!!!
!!!---------------------------------------------------------------------
!!!

  subroutine initialdata_disturb_dryreg( &
    & XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, &
    & xyzf_QMix)
    !
    ! 矩形な乾燥領域を作る
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,             &! X 座標軸(スカラー格子点)
      &                  y_Y,             &! X 座標軸(スカラー格子点)
      &                  z_Z               ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,      &! 配列サイズ (X 方向)
      &                  jmin, jmax,      &! 配列サイズ (Y 方向)
      &                  kmin, kmax,      &! 配列サイズ (Z 方向)
      &                  ncmax             ! 計算領域のマージン
    use basicset, only: xyzf_QMixBZ
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax
    real(DP), intent(out) :: xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, 1:ncmax)
    integer               :: i, j, k, s
    
    ! XposMin:XposMax,ZposMin:ZposMax で囲まれた領域の初期の湿度をゼロにするために
    ! 基本場と逆符号の水蒸気擾乱を与える
    do s = 1, ncmax
      do k = kmin,kmax  
        do j = jmin, jmax
          do i = imin,imax
            if (z_Z(k) >= ZposMin .AND. z_Z(k) < ZposMax &
              & .AND. y_Y(j) >= YposMin .AND. y_Y(j) < YposMax &
              & .AND. x_X(i) >= XposMin .AND. x_X(i) < XposMax) then
              xyzf_QMix(i,j,k,s) = - xyzf_QMixBZ(i,j,k,s)
            end if
          end do
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_dryreg
  
  
  subroutine initialdata_disturb_square( &
    & DelMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, &
    & xyz_Var)
    !
    ! 立方体
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X 座標軸(スカラー格子点)
      &                  y_Y,           &! Y 座標軸(スカラー格子点)
      &                  z_Z             ! Z 座標軸(スカラー格子点)
    use gridset,   only: imin, imax,    &! 配列サイズ (X 方向)
      &                  jmin, jmax,    &! 配列サイズ (Y 方向)
      &                  kmin, kmax      ! 配列サイズ (Z 方向)

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: DelMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k

  
    do k = kmin, kmax  
      do j = jmin, jmax
        do i = imin, imax
          if (z_Z(k) >= ZposMin .AND. z_Z(k) <= ZposMax &
            & .AND. y_Y(j) >= YposMin .AND. y_Y(j) <= YposMax &
            & .AND. x_X(i) >= XposMin .AND. x_X(i) <= XposMax) then
            xyz_Var(i,j,k) = DelMax
!            write(*,*) x_X(i), y_Y(j), z_Z(k), xyz_Var(i,j,k)
          end if
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_square
  
  
  subroutine initialdata_disturb_moist(Hum, xyzf_QMix)
    !
    ! 湿度 Hum な混合比の鉛直分布を作成
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use gridset,   only: nx, imin, imax,   &! 配列サイズ (X 方向)
      &                  ny, jmin, jmax,   &! 配列サイズ (Y 方向)
      &                  nz, kmin, kmax,   &! 配列サイズ (Z 方向)
      &                  ncmax              ! 計算領域のマージン
    use basicset,  only: xyz_TempBZ,       &! 基本場の温度
      &                  xyz_PressBZ,      &! 基本場の圧力
      &                  xyzf_QMixBZ        ! 基本場の混合比
    use composition, only:                 &
      &                  MolWtWet,         &!凝縮成分の分子量
      &                  SpcWetMolFr        !凝縮成分の初期モル比
    use constants, only: MolWtDry           !乾燥成分の分子量
    use eccm,      only: eccm_molfr

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: Hum
    real(DP), intent(out) :: xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, 1:ncmax)
    real(DP)              :: zf_MolFr(kmin:kmax, 1:ncmax)
    integer               :: i, j, k, s
  
    ! 湿度ゼロなら何もしない
    if ( Hum == 0.0d0 ) return

    ! 水平一様なので, i=0 だけ計算. 
    i = 1
    j = 1
    call eccm_molfr( SpcWetMolFr(1:ncmax), Hum, xyz_TempBZ(i,j,:), &
      &              xyz_PressBZ(i,j,:), zf_MolFr )
    
    !気相のモル比を混合比に変換
    do s = 1, ncmax
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyzf_QMix(i,j,k,s) = zf_MolFr(k,s) * MolWtWet(s) / MolWtDry - xyzf_QMixBZ(i,j,k,s)
          end do
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_moist
  
end module initialdata_disturb
