!= 格子設定用モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: gridset.f90,v 1.14 2014/07/08 01:05:33 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module gridset
  !
  != 格子設定用モジュール
  !
  !引数に与えられた NAMELIST ファイルから, 格子点情報を取得し, 
  !保管するための変数参照型モジュール
  !
  ! 2D 版と 3D 版で RegXMin, RegYMin, RegZMin の定義が異なることに注意
  !
  !== Procedures List
  ! gridset_init   :: 初期化ルーチン
  ! gridset_check  :: 設定内容の整合性をチェックするためのルーチン
  !

  !暗黙の型宣言禁止
  implicit none
  
  !save 属性
  private
  
  !公開変数
  integer, public, save :: NX = 10 ! x 方向格子点数
  integer, public, save :: NY = 10 ! y 方向格子点数
  integer, public, save :: NZ = 10 ! z 方向格子点数
  integer, public, save :: NCMAX = 1  ! 組成配列要素数
  integer, public, save :: Xmg = 2 ! x 方向糊代格子点数
  integer, public, save :: Ymg = 2 ! y 方向糊代格子点数
  integer, public, save :: Zmg = 2 ! z 方向糊代格子点数
  integer, public, save :: imin    ! x 方向の配列の下限 
  integer, public, save :: imax    ! x 方向の配列の上限
  integer, public, save :: jmin    ! y 方向の配列の下限 
  integer, public, save :: jmax    ! y 方向の配列の下限 
  integer, public, save :: kmin    ! z 方向の配列の下限 
  integer, public, save :: kmax    ! z 方向の配列の下限 
  integer, public, save :: xsub = 1
                                   ! Number of MPI processes (X direction)
  integer, public, save :: ysub = 1
                                   ! Number of MPI processes (Y direction)
  logical, public, save :: FlagCalc3D = .true.
  logical, public, save :: FlagCalcMoist = .true.

  integer, private, save :: xdim = 1 
  integer, private, save :: ydim = 1
  integer, private, save :: zdim = 1
  integer, public, save  :: nxny = 100 ! x 方向と y 方向の総格子点数

  integer, public, save :: MarginX       !境界のグリッド数
  integer, public, save :: MarginY       !境界のグリッド数
  integer, public, save :: MarginZ       !境界のグリッド数
  integer, public, save :: SpcNum        !化学種の数
  integer, public, save :: DimXMin       ! x 方向の配列の下限
  integer, public, save :: DimXMax       ! x 方向の配列の上限
  integer, public, save :: DimYMin       ! y 方向の配列の下限
  integer, public, save :: DimYMax       ! y 方向の配列の上限
  integer, public, save :: DimZMin       ! z 方向の配列の下限
  integer, public, save :: DimZMax       ! z 方向の配列の上限
  integer, public, save :: RegXMin       ! x 方向の物理領域の下限
  integer, public, save :: RegXMax       ! x 方向の物理領域の上限
  integer, public, save :: RegYMin       ! y 方向の物理領域の下限
  integer, public, save :: RegYMax       ! y 方向の物理領域の上限
  integer, public, save :: RegZMin       ! z 方向の物理領域の下限
  integer, public, save :: RegZMax       ! z 方向の物理領域の上限
  integer, public, save :: FileNX        !ファイル出力用
  integer, public, save :: FileNY        !ファイル出力用
  integer, public, save :: FileNZ        !ファイル出力用
  integer, public, save :: FileXMin      !ファイル出力用
  integer, public, save :: FileXMax      !ファイル出力用
  integer, public, save :: FileYMin      !ファイル出力用
  integer, public, save :: FileYMax      !ファイル出力用
  integer, public, save :: FileZMin      !ファイル出力用
  integer, public, save :: FileZMax      !ファイル出力用

  integer, public, save :: im = 10 ! x 方向格子点数
  integer, public, save :: jm = 10 ! y 方向格子点数
  integer, public, save :: km = 10 ! z 方向格子点数

  ! サブルーチンの公開
  public gridset_init

contains

  subroutine gridset_init
    !
    != 初期化ルーチン
    !
    ! 設定ファイルから情報を読み込み格子点数を計算する
    !

    !モジュール読み込み
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use mpi_wrapper,   only: nprocs
    use namelist_util, only: namelist_filename
    
    !暗黙の型宣言禁止
    implicit none

    !内部変数
    integer            :: unit                !設定ファイル用装置番号

    !-----------------------------------------------------------------
    ! 設定ファイルから情報を読み込み
    !
    NAMELIST /gridset_nml/ xdim, ydim, zdim, NCMAX, Xmg, Ymg, Zmg, xsub, ysub

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=gridset_nml)
    close(unit)

    !-----------------------------------------------------------------    
    ! NX, NY, NZ を決める
    !
    nx = xdim / xsub 
    imin = 0  - xmg   ! 3D 版と異なる
    imax = nx + xmg
   
    if (ydim == 1) then
      ! 二次元
      ymg = 0
      ny = 1
      FlagCalc3D = .false. 
    else
      ! 三次元
      ny = ydim / ysub
      FlagCalc3D = .true. 
    end if
    jmin = 0  - ymg   ! 3D 版と異なる
    jmax = ny + ymg

    nz = zdim
    kmin = 0  - zmg   ! 3D 版と異なる
    kmax = nz + zmg    

    nxny = nx * ny

    ! トレーサー計算のスイッチ
    !
    if ( ncmax == 0 ) then 
!      FlagCalcTracer = .false.
      FlagCalcMoist = .false.
    else
!      FlagCalcTracer = .true.
      FlagCalcMoist = .true.
    end if

    !-----------------------------------------------------------------    
    ! 読み込んだ情報を出力
    !
    call MessageNotify( "M", "gridset_init", "xsub = %d",  i=(/xsub/) )
    call MessageNotify( "M", "gridset_init", "ysub = %d",  i=(/ysub/) )
    call MessageNotify( "M", "gridset_init", "xdim = %d",   i=(/xdim/) )
    call MessageNotify( "M", "gridset_init", "ydim = %d",   i=(/ydim/) )
    call MessageNotify( "M", "gridset_init", "zdim = %d",   i=(/zdim/) )
    call MessageNotify( "M", "gridset_init", "[1 node] NX = %d",   i=(/NX/) )
    call MessageNotify( "M", "gridset_init", "[1 node] NY = %d",   i=(/NY/) )
    call MessageNotify( "M", "gridset_init", "[1 node] NZ = %d",   i=(/NZ/) )
    call MessageNotify( "M", "gridset_init", "[1 node] NCMAX = %d",   i=(/NCMAX/) )
    call MessageNotify( "M", "gridset_init", "[1 node] xmg  = %d", i=(/Xmg/) )
    call MessageNotify( "M", "gridset_init", "[1 node] ymg  = %d", i=(/Ymg/) )
    call MessageNotify( "M", "gridset_init", "[1 node] zmg  = %d", i=(/Zmg/) )
    call MessageNotify( "M", "gridset_init", "[1 node] imin = %d", i=(/imin/) )
    call MessageNotify( "M", "gridset_init", "[1 node] imax = %d", i=(/imax/) )
    call MessageNotify( "M", "gridset_init", "[1 node] jmin = %d", i=(/jmin/) )
    call MessageNotify( "M", "gridset_init", "[1 node] jmax = %d", i=(/jmax/) )
    call MessageNotify( "M", "gridset_init", "[1 node] kmin = %d", i=(/kmin/) )
    call MessageNotify( "M", "gridset_init", "[1 node] kmax = %d", i=(/kmax/) )

    !-----------------------------------------------------------------
    ! 値のチェック
    !    
    call gridset_check

    !-----------------------------------------------------------------
    ! arare4 用の変数を設定
    !
    im = nx
    jm = ny
    km = nz

    MarginX = Xmg
    MarginY = Ymg
    MarginZ = Zmg
    SpcNum  = ncmax

    RegXMin = 0    !3D 版と異なる
    RegXMax = NX
    RegYMin = 0    !3D 版と異なる
    RegYMax = NY
    RegZMin = 0    !3D 版と異なる
    RegZMax = NZ

    DimXMin = RegXMin - MarginX
    DimXMax = RegXMax + MarginX
    DimYMin = RegYMin - MarginY
    DimYMax = RegYMax + MarginY
    DimZMin = RegZMin - MarginZ
    DimZMax = RegZMax + MarginZ

    FileNX = NX
    FileNY = NY
    FileNZ = NZ

    FileXMin = RegXMin + 1    !3D 版と異なる
    FileXMax = RegXMax
    FileYMin = RegYMin + 1    !3D 版と異なる
    FileYMax = RegYMax
    FileZMin = RegZMin + 1    !3D 版と異なる
    FileZMax = RegZMax

  contains

    subroutine gridset_check
      !
      != 設定内容の整合性をチェックするためのルーチン
      !
      ! 設定内容に矛盾がないか確認をする. 
      ! 問題が有る場合はエラーメッセージを出して終了する.
      ! 

      !暗黙の型宣言禁止
      implicit none
      
      ! nprocs が xsub or ysub で割り切れるかチェック. 
      !
      if ( mod( nprocs, xsub ) /= 0) then
        call MessageNotify( "E", "gridset_init", "mod( nprocs / xsub ) /= 0" )
      end if
      
      if ( mod( nprocs, ysub ) /= 0) then
        call MessageNotify( "E", "gridset_init", "mod( nprocs / ysub ) /= 0" )
      end if
      
      if (xsub * ysub /= nprocs) then
        call MessageNotify( "E", "gridset_init", "xsub * ysub /= nprocs" )
      end if
      
      ! X 方向のマージンの大きさのチェック
      ! 
      if ( xdim < xmg ) then
        call MessageNotify( "E", "gridset_init", "xdim < Xmg" )
      end if
      
      ! xdim が xsub で割り切れるか確認
      ! 
      if ( mod(xdim, xsub) /= 0  ) then
        call MessageNotify( "E", "gridset_init", "mod(xdim, xsub) /= 0" )
      end if
      
      ! Y 方向のマージンの大きさのチェック
      ! 
      if ( ydim < ymg ) then
        call MessageNotify( "E", "gridset_init", "ydim < Ymg" )
      end if
      
      ! ydim が xsub で割り切れるか確認
      ! 
      if ( mod(ydim, ysub) /= 0  ) then
        call MessageNotify( "E", "gridset_init", "mod(ydim, ysub) /= 0" )
      end if
      
      ! Z 方向のマージンの大きさのチェック
      ! 
      if ( zdim < zmg ) then
        call MessageNotify( "E", "gridset_init", "zdim < Zmg" )
      end if
      
    end subroutine gridset_check

  end subroutine gridset_init
  
end module gridset
