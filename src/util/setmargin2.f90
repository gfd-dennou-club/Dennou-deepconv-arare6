!= Module Damping
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: setmargin2.f90,v 1.4 2015/02/20 08:03:32 sugiyama Exp $
! Tag Name::  $Name: arare5-20150220 $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module SetMargin
  !
  ! 3 次元 (xyz 方向) 不等間隔交互格子 有限差分モデル用 境界条件モジュール
  ! X 方向・Y 方向には MPI を用いた並列化を行えるようにしている
  !

  use dc_types, only : DP
  use gridset,  only : imin, imax, jmin, jmax, kmin, kmax,    &
    &                  xmg, ymg, zmg, nx, ny, nz, ncmax, xsub, ysub
  use namelist_util,only: namelist_filename
  use mpi_wrapper, only : myrank, nprocs, &
    &                     MPIWrapperISend, MPIWrapperIRecv, MPIWrapperWait, &
    &                     MPIWrapperCartCreate, MPIWrapperCartShift, MPIWrapperCommFree

  implicit none

  private

  integer, save :: Rrank = 1
                           ! 右隣の CPU の rank
  integer, save :: Lrank = 1
                           ! 左隣の CPU の rank
  integer, save :: Urank = 1
                           ! 上側の CPU の rank
  integer, save :: Drank = 1
                           ! 下側の CPU の rank

  public :: SetMargin_init
  public :: SetMargin_xyzf 
  public :: SetMargin_xyz
  public :: SetMargin_pyz
  public :: SetMargin_xqz
  public :: SetMargin_xyr

  interface SetMargin_xyz
    module procedure SetMargin_aaz
  end interface

  interface SetMargin_pyz
    module procedure SetMargin_aaz
  end interface

  interface SetMargin_xqz
    module procedure SetMargin_aaz
  end interface

  interface SetMargin_xyr
    module procedure SetMargin_aar
  end interface
  
contains

!!!-------------------------------------------------------------
  subroutine XcyclYcycl(aaa_Var)

    implicit none
    real(DP),intent(inout) :: aaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(8)              :: Xsbuf_a( Xmg, jmin:jmax, kmin:kmax )
    real(8)              :: Xrbuf_a( Xmg, jmin:jmax, kmin:kmax )
    real(8)              :: Xsbuf_b( Xmg, jmin:jmax, kmin:kmax )
    real(8)              :: Xrbuf_b( Xmg, jmin:jmax, kmin:kmax )
    real(8), allocatable :: Ysbuf_a( :,:,: )
    real(8), allocatable :: Yrbuf_a( :,:,: )
    real(8), allocatable :: Ysbuf_b( :,:,: )
    real(8), allocatable :: Yrbuf_b( :,:,: )
    integer              :: ireqs_a, ireqr_a, ireqs_b, ireqr_b
    integer              :: ny2, nz2

!    real(8) :: Ysbuf_a( nx, 0:ymg, nz ) 
!    real(8) :: Yrbuf_a( nx, 0:ymg, nz ) 
!    real(8) :: Ysbuf_b( nx, 0:ymg, nz ) 
!    real(8) :: Yrbuf_b( nx, 0:ymg, nz ) 


    !!!
    !!! Y 方向に通信する
    !!!   single CPU で動かす場合は, MPIWrapper* サブルーチンは何もしないことに注意
    !!!

    !-------------------------------
    ! 送信する部分配列(sbuf_a, sbuf_b)を用意する
    if ( ymg > 0 ) then 
      
      allocate( Ysbuf_a( nx, ymg, nz ) )
      allocate( Yrbuf_a( nx, ymg, nz ) )
      allocate( Ysbuf_b( nx, ymg, nz ) )
      allocate( Yrbuf_b( nx, ymg, nz ) )
      
      Ysbuf_a( 1:nx, 1:ymg, 1:nz ) = aaa_Var( 1:nx, ny-ymg+1:ny, 1:nz )
      Ysbuf_b( 1:nx, 1:ymg, 1:nz ) = aaa_Var( 1:nx, 1:ymg, 1:nz )
      
      !-------------------------------
      ! 送受信開始
      !
      ! サブルーチンを call する順番は 2CPU の場合を考慮して決める必要がある
      ! rank0 にて上側・下側の順に送信する
      !  => rank0 は rank1 から下側・上側の順にデータを受診する
      
      Yrbuf_a = Ysbuf_b  ! single CPU で用いるときのための処理
      Yrbuf_b = Ysbuf_a  ! single CPU で用いるときのための処理
      
      call MPIWrapperISend( Urank, nx, ymg, nz, Ysbuf_a, ireqs_a ) !送信
      call MPIWrapperISend( Drank, nx, ymg, nz, Ysbuf_b, ireqs_b ) !送信

      call MPIWrapperIRecv( Drank, nx, ymg, nz, Yrbuf_b, ireqr_b ) !受信      
      call MPIWrapperIRecv( Urank, nx, ymg, nz, Yrbuf_a, ireqr_a ) !受信
      
      !-------------------------------
      ! 送受信の終了を待ってから代入
      
      call MPIWrapperWait( ireqs_a )       
      call MPIWrapperWait( ireqr_a )      
      call MPIWrapperWait( ireqs_b )
      call MPIWrapperWait( ireqr_b )
      aaa_var( 1:nx, ny+1:ny+ymg, 1:nz ) = Yrbuf_a( 1:nx, 1:ymg, 1:nz ) 
      aaa_var( 1:nx, -ymg+1:0, 1:nz ) = Yrbuf_b( 1:nx, 1:ymg, 1:nz ) 
      
    end if

    !!!
    !!! X 方向に周期境界条件を適用
    !!!   single CPU で動かす場合は, MPIWeapper* サブルーチンは何もしないことに注意
    !!!
    !-------------------------------
    ! 送信する部分配列(sbuf_a, sbuf_b)を用意する
    Xsbuf_a( 1:xmg, jmin:jmax, kmin:kmax ) = aaa_Var( nx-xmg+1:nx, jmin:jmax, kmin:kmax )  ! nx - xmg + 1 ~ nx を保管
    Xsbuf_b( 1:xmg, jmin:jmax, kmin:kmax ) = aaa_Var( 1:xmg, jmin:jmax, kmin:kmax )  ! 1 ~ xmg を保管

    !-------------------------------
    ! 送受信開始
    ! サブルーチンを call する順番は 2CPU の場合を考慮して決める必要がある
    ! rank0 にて右側・左側の順に送信する
    !  => rank0 は rank1 より左側・右側の順にデータを受診
    
    Xrbuf_a = Xsbuf_b  ! single CPU で用いるときのための処理
    Xrbuf_b = Xsbuf_a  ! single CPU で用いるときのための処理

    ny2 = ny + ymg * 2
    nz2 = nz + zmg * 2

    call MPIWrapperISend( Rrank, xmg, ny2, nz2, Xsbuf_a, ireqs_a ) !送信
    call MPIWrapperISend( Lrank, xmg, ny2, nz2, Xsbuf_b, ireqs_b ) !送信

    call MPIWrapperIRecv( Lrank, xmg, ny2, nz2, Xrbuf_b, ireqr_b ) !受信    
    call MPIWrapperIRecv( Rrank, xmg, ny2, nz2, Xrbuf_a, ireqr_a ) !受信
    
    !-------------------------------
    ! 送受信の終了を待ってから代入

    call MPIWrapperWait( ireqs_a )    
    call MPIWrapperWait( ireqr_a )
    call MPIWrapperWait( ireqs_b )       
    call MPIWrapperWait( ireqr_b )

    aaa_var( nx+1:nx+xmg, jmin:jmax, kmin:kmax ) = Xrbuf_a( 1:xmg, jmin:jmax, kmin:kmax ) ! nx + 1 ~ nx + xmg に代入
    aaa_var( -xmg+1:0, jmin:jmax, kmin:kmax ) = Xrbuf_b( 1:xmg, jmin:jmax, kmin:kmax ) ! - xmg + 1 ~ 0 に代入

  end subroutine XcyclYcycl

!!!-------------------------------------------------------------
  subroutine ZSym(aaz_Var)
    ! z 方向に対称境界条件を適用する
    !
    implicit none

    real(DP),intent(inout) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    integer :: k

    do k = 1, zmg
      aaz_Var(:,:,1-k)  = aaz_Var(:,:,k)
      aaz_Var(:,:,nz+k) = aaz_Var(:,:,nz+1-k)
    end do

  end subroutine ZSym

!!!-------------------------------------------------------------
  subroutine ZAntSym(aar_Var)
    ! z 方向に対称境界条件を適用する
    !
    implicit none

    real(DP),intent(inout) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax) 
    integer :: k

    aar_Var(:,:,0) = 0.0d0
    aar_Var(:,:,nz) = 0.0d0
    
    do k = 1, zmg-1
      aar_Var(:,:,-k) = - aar_Var(:,:,k)
    end do
    
    do k = 1, zmg
      aar_Var(:,:,nz+k) = - aar_Var(:,:,nz-k)
    end do
    
  end subroutine ZAntSym

!!!-------------------------------------------------------------
  subroutine SetMargin_xyzf(xyzf_Var)

    implicit none
    real(DP),intent(inout) :: xyzf_Var(imin:imax,jmin:jmax,kmin:kmax,1:ncmax) 
    integer                :: s


    do s = 1, ncmax

      ! x, y 方向に周期境界条件を適用する
      !
      call XCyclYCycl( xyzf_Var(:,:,:,s) ) !inout

      ! z 方向に対象境界条件を適用する
      !
      call ZSym( xyzf_Var(:,:,:,s) ) !inout
      
    end do

  end subroutine SetMargin_xyzf

!!!-------------------------------------------------------------
  subroutine SetMargin_aaz(aaz_Var)

    implicit none
    real(DP),intent(inout) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 


    ! x, y 方向に周期境界条件を適用する
    !
    call XCyclYCycl( aaz_Var(:,:,:) ) !inout
    
    ! z 方向に対象境界条件を適用する
    !
    call ZSym( aaz_Var(:,:,:) ) !inout
       
  end subroutine SetMargin_aaz

 
!!!-------------------------------------------------------------
  subroutine SetMargin_aar(aar_Var)

    implicit none
    real(DP),intent(inout) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax) 


    ! x, y 方向に周期境界条件を適用する
    !
    call XCyclYCycl( aar_Var(:,:,:) ) !inout
    
    ! z 方向に対象境界条件を適用する
    !
    call ZAntSym( aar_Var(:,:,:) ) !inout
    
    
  end subroutine SetMargin_aar

!!!-------------------------------------------------------------
  subroutine SetMargin_init
    !
    ! 初期化
    !
    ! Initialization of SetMargin 
    !

    ! モジュール引用 ; USE statements
    !
    use dc_message,  only: MessageNotify

    ! 作業変数
    ! Work variables
    !
    integer, parameter :: disp = 1
    integer            :: direction
    integer            :: comm_cart
    logical, parameter :: periodic = .true.

    ! 通信先の取得   
    !
    call MPIWrapperCartCreate(xsub, ysub, periodic, comm_cart)
    direction = 0
    call MPIWrapperCartShift(comm_cart, direction, disp, Urank, Drank)
    direction = 1
    call MPIWrapperCartShift(comm_cart, direction, disp, Lrank, Rrank)
    call MPIWrapperCommFree(comm_cart)

    if (myrank == 0) then 
      call MessageNotify( "M", "SetMargin_init", "Rrank = %d",   i=(/Rrank/) )
      call MessageNotify( "M", "SetMargin_init", "Lrank = %d",   i=(/Lrank/) )
      call MessageNotify( "M", "SetMargin_init", "Urank = %d",   i=(/Urank/) )
      call MessageNotify( "M", "SetMargin_init", "Drank = %d",   i=(/Drank/) )
    end if
    
  end subroutine SetMargin_init
  
end module SetMargin
