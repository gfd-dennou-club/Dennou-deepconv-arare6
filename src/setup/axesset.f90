!= 3 次元 (xyz 方向) 等間隔交互格子 格子点設定モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: axesset.f90,v 1.16 2014/07/08 01:05:32 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module axesset
  != 3 次元 (xyz 方向) 等間隔交互格子 格子点設定モジュール
  !
  !== 概要
  !
  ! axesset は, 3 次元 (xyz 方向) 等間隔交互格子を用いた有限差分法に基づく
  ! 軸情報を提供する. 
  ! 

  !モジュール読み込み
  !
  use dc_types,      only: DP

  !暗黙の型宣言禁止
  !
  implicit none

  !デフォルトは非公開
  !
  private

  !公開要素
  ! 初期化手続き
  !
  public :: axesset_init

  ! 変数定義
  real(DP), public, save :: Xmin = 0.0d0                ! x 座標の始点・終点
  real(DP), public, save :: Xmax = 1.0d4                ! x 座標の始点・終点
  real(DP), public, save :: Ymin = 0.0d0                ! x 座標の始点・終点
  real(DP), public, save :: Ymax = 1.0d4                ! x 座標の始点・終点
  real(DP), public, save :: Zmin = 0.0d0                ! z 座標の始点・終点
  real(DP), public, save :: Zmax = 1.0d4                ! z 座標の始点・終点
  real(DP), public, save :: DX  
  real(DP), public, save :: DY           
  real(DP), public, save :: DZ           
  real(DP), allocatable, public, save :: x_X(:)         ! 半整数格子点座標
  real(DP), allocatable, public, save :: p_X(:)         ! 整数格子点座標
  real(DP), allocatable, public, save :: x_dx(:)        ! 半整数格子点間隔
  real(DP), allocatable, public, save :: p_dx(:)        ! 整数格子点間隔
  real(DP), allocatable, public, save :: y_Y(:)         ! 半整数格子点座標
  real(DP), allocatable, public, save :: q_Y(:)         ! 整数格子点座標
  real(DP), allocatable, public, save :: y_dy(:)        ! 半整数格子点間隔
  real(DP), allocatable, public, save :: q_dy(:)        ! 整数格子点間隔
  real(DP), allocatable, public, save :: z_Z(:)         ! 半整数格子点座標
  real(DP), allocatable, public, save :: r_Z(:)         ! 整数格子点座標
  real(DP), allocatable, public, save :: z_dz(:)        ! 半整数格子点間隔
  real(DP), allocatable, public, save :: r_dz(:)        ! 整数格子点間隔
  real(DP), allocatable, public, save :: xyz_X(:,:,:)   ! x 座標(半整数格子)
  real(DP), allocatable, public, save :: xyz_Y(:,:,:)   ! y 座標(半整数格子)
  real(DP), allocatable, public, save :: xyz_Z(:,:,:)   ! z 座標(半整数格子)
  real(DP), allocatable, public, save :: xyz_dX(:,:,:)  ! x 格子間隔(半整数格子)
  real(DP), allocatable, public, save :: xyz_dY(:,:,:)  ! y 格子間隔(半整数格子)
  real(DP), allocatable, public, save :: xyz_dZ(:,:,:)  ! z 格子間隔(半整数格子)

  real(DP), public, save :: DelX                        !格子点間隔
  real(DP), public, save :: DelY                        !格子点間隔
  real(DP), public, save :: DelZ                        !格子点間隔
  real(DP), allocatable, public, save :: s_X(:)         !X 座標軸(スカラー格子点)
  real(DP), allocatable, public, save :: f_X(:)         !X 座標軸(ベクトル格子点)
  real(DP), allocatable, public, save :: s_Y(:)         !Y 座標軸(スカラー格子点)
  real(DP), allocatable, public, save :: f_Y(:)         !Y 座標軸(ベクトル格子点)
  real(DP), allocatable, public, save :: s_Z(:)         !Z 座標軸(スカラー格子点)
  real(DP), allocatable, public, save :: f_Z(:)         !Z 座標軸(ベクトル格子点)

contains
  !--------------------------------------------------------------------
  subroutine axesset_init
    !
    ! 格子点座標配列と格子点間隔配列の初期化    

    ! モジュール呼び出し
    !
    use dc_types,      only: DP
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use mpi_wrapper,   only: MPIWrapperCartCreate, &
      &                      MPIWrapperCartShift,  &
      &                      MPIWrapperCommFree
    use gridset,       only: FlagCalc3D,    &
      &                      xsub, ysub,    & ! 
      &                      nx, ny, nz,    & ! 格子点数
      &                      imin, imax, jmin, jmax, kmin, kmax ! 配列の上限値と下限値
    use namelist_util, only: namelist_filename

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 変数定義
    !
    real(DP),allocatable :: xy_X(:,:)! x 座標(半整数格子, 作業配列)
    real(DP),allocatable :: xy_Y(:,:)! y 座標(半整数格子, 作業配列)
    real(DP),allocatable :: yz_Z(:,:)! z 座標(半整数格子, 作業配列)
    integer              :: unit     ! 設定ファイル用装置番号
    integer              :: comm_cart
    logical, parameter   :: periodic = .false.

    !設定ファイルから読み込む出力ファイル情報
    NAMELIST /axesset_nml/ xmin, xmax, ymin, ymax, zmin, zmax

    !設定ファイルから出力ファイルに記載する情報を読み込む
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=axesset_nml)
    close(unit)
    
    ! 配列の上下限の値, 座標値と格子点間隔を設定
    ! * 1 次元用のサブルーチンを用いる
    !
    call MPIWrapperCartCreate(xsub, ysub, periodic, comm_cart)
    call x_axis_init(comm_cart)
    call y_axis_init(comm_cart)
    call z_axis_init
    call MPIWrapperCommFree(comm_cart)
    
    ! 3 次元格子点座標配列の設定
    ! * 組み込み関数 spread を用いる. 
    ! * 中間配列として 2 次元格子点座標配列を作り, それを 3 次元に拡張する.
    ! 
    allocate(xy_X(imin:imax,jmin:jmax))
    allocate(xy_Y(imin:imax,jmin:jmax))
    allocate(yz_Z(jmin:jmax,kmin:kmax))
    
    allocate(xyz_X(imin:imax,jmin:jmax,kmin:kmax))
    allocate(xyz_Y(imin:imax,jmin:jmax,kmin:kmax))
    allocate(xyz_Z(imin:imax,jmin:jmax,kmin:kmax))
    allocate(xyz_dX(imin:imax,jmin:jmax,kmin:kmax))
    allocate(xyz_dY(imin:imax,jmin:jmax,kmin:kmax))
    allocate(xyz_dZ(imin:imax,jmin:jmax,kmin:kmax))
    
    xy_X  = spread(x_X, 2,size(y_Y))
    xyz_X = spread(xy_X,3,size(z_Z))
    
    xy_X   = spread(x_dX, 2,size(y_dY))
    xyz_dX = spread(xy_X,3,size(z_dZ))
    
    xy_Y  = spread(y_Y, 1,size(x_X))
    xyz_Y = spread(xy_Y,3,size(z_Z))
    
    xy_Y   = spread(y_dY, 1,size(x_dX))
    xyz_dY = spread(xy_Y,3,size(z_dZ))
    
    yz_Z  = spread(z_Z, 1,size(y_Y))
    xyz_Z = spread(yz_Z,1,size(x_X))
    
    yz_Z   = spread(z_dZ, 1,size(y_dY))
    xyz_dZ = spread(yz_Z,1,size(x_dX))
    
    deallocate(xy_X)
    deallocate(xy_Y)
    deallocate(yz_Z)

    ! 読み込んだ情報を出力
    call MessageNotify( "M", "axesset_init", "XMin = %f", d=(/XMin/)    )
    call MessageNotify( "M", "axesset_init", "XMax = %f", d=(/XMax/)    )
    call MessageNotify( "M", "axesset_init", "YMin = %f", d=(/YMin/)    )
    call MessageNotify( "M", "axesset_init", "YMax = %f", d=(/YMax/)    )
    call MessageNotify( "M", "axesset_init", "ZMin = %f", d=(/ZMin/)    )
    call MessageNotify( "M", "axesset_init", "ZMax = %f", d=(/ZMax/)    )
    call MessageNotify( "M", "axesset_init", "dx = %f", d=(/dx/)    )
    call MessageNotify( "M", "axesset_init", "dy = %f", d=(/dy/)    )
    call MessageNotify( "M", "axesset_init", "dz = %f", d=(/dz/)    )

    !----------------------------------------------------------------------
    ! arare4 用の設定
    !
    DelX = DX
    DelY = DY
    DelZ = DZ

    allocate(s_X(imin:imax), f_X(imin:imax))
    allocate(s_Y(jmin:jmax), f_Y(jmin:jmax))
    allocate(s_Z(kmin:kmax), f_Z(kmin:kmax))

    s_X = x_X
    f_X = p_X
    s_Y = y_Y
    f_Y = q_Y
    s_Z = z_Z
    f_Z = r_Z

  contains

    subroutine x_axis_init(comm_cart)
      !
      != 格子点座標配列と格子点間隔配列の初期化    
      !
      
      ! 暗黙の型宣言禁止
      implicit none
      
      integer, intent(in) :: comm_cart
      
      ! 作業変数
      integer             :: ix   ! do ループ添字
      integer             :: rank, tmp
      integer             :: lp
      integer, parameter  :: direction = 1
      
      allocate(x_X(imin:imax))
      allocate(p_X(imin:imax))
      allocate(x_dx(imin:imax))
      allocate(p_dx(imin:imax))
      
      ! 初期化
      lp = -999
      
      ! 等間隔格子
      dx = (xmax - xmin) / ( nx * xsub )
      
      ! 隣り合う CPU の個数を得る
      sx: do ix = 1, xsub
        call MPIWrapperCartShift( comm_cart, direction, ix, rank, tmp )
        if ( rank < 0 ) then 
          lp = ix - 1
          exit sx
        end if
      end do sx
      
      !確認
      if ( lp == -999 ) then 
        call MessageNotify( "E", "axesset_init", "lp (x) is undifined" )
      end if
      
      ! 配列の作成
      do ix = imin, imax
        p_X(ix) = XMin + dx * nx * lp + dx * ix 
        x_X(ix) = XMin + dx * nx * lp + dx * ix - dx * 0.5d0
        x_dx(ix) = dx
        p_dx(ix) = dx
      end do
      
    end subroutine x_axis_init
    
    !--------------------------------------------------------------------
    subroutine y_axis_init(comm_cart)
      !
      != y 方向の座標値と格子点間隔を設定する
      !
      
      ! 暗黙の型宣言禁止
      implicit none
      
      integer, intent(in) :: comm_cart
      
      ! 作業変数
      integer             :: jy   ! ループ添字
      integer             :: lp
      integer             :: rank, tmp
      integer, parameter  :: direction = 0
      
      allocate(y_Y(jmin:jmax))
      allocate(q_Y(jmin:jmax))
      allocate(y_dy(jmin:jmax))
      allocate(q_dy(jmin:jmax))
      
      ! 初期化
      lp = -999
      
      ! 2 次元の場合は ymax = dx とする. 
      if (.NOT. FlagCalc3D ) then 
        ymax = dx
      end if
      
      ! 等間隔格子
      dy = (ymax - ymin) / ( ny * ysub )
      
      ! 隣り合う CPU の個数を得る
      sy: do jy = 1, ysub
        call MPIWrapperCartShift( comm_cart, direction, jy, tmp, rank )
        if ( rank < 0 ) then 
          lp = jy - 1
          exit sy
        end if
      end do sy
      
      !確認
      if ( lp == -999 ) then 
        call MessageNotify( "E", "axesset_init", "lp (y) is undifined" )
      end if
      
      ! 配列の作成
      do jy = jmin, jmax
        q_Y(jy) = YMin + dy * ny * lp  + dy * jy
        y_Y(jy) = YMin + dy * ny * lp  + dy * jy - dy * 0.5d0
        y_dy(jy) = dy
        q_dy(jy) = dy
      end do
      
    end subroutine y_axis_init


    subroutine z_axis_init
      !
      != z 方向の座標値と格子点間隔を設定する
      !
      
      ! 暗黙の型宣言禁止
      implicit none
      
      ! 作業変数
      integer                 :: kz   ! ループ添字
      
      allocate(z_Z(kmin:kmax))
      allocate(r_Z(kmin:kmax))
      allocate(z_dz(kmin:kmax))
      allocate(r_dz(kmin:kmax))
      
      ! 等間隔格子
      dz = (zmax - zmin) / nz
      
      ! 配列の作成
      do kz = kmin, kmax
        r_Z(kz) = ZMin + dz * kz
        z_Z(kz) = ZMin + dz * (kz - 0.5)
        r_dz(kz) = dz
        z_dz(kz) = dz
      end do
    
    end subroutine z_axis_init
  
  end subroutine axesset_init

end module axesset

