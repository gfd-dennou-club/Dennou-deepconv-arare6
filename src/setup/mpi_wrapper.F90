!= MPI 関連ルーチン
!
!= MPI related routines
!
! Authors::   Yoshiyuki O. Takahashi, SUGIYAMA Ko-ichiro
! Version::   $Id: mpi_wrapper.F90,v 1.10 2014/03/04 04:49:44 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module mpi_wrapper
  !
  != MPI 関連ルーチン
  !
  != MPI related routines
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! MPI 関係の変数の管理と MPI 関係ラッパールーチンのモジュール. 
  !
  ! This is a module containing MPI-related variables and wrapper routines. 


  ! モジュール引用 ; USE statements
  !

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types, only: DP    ! 倍精度実数型. Double precision.

#ifdef LIB_MPI
  ! MPI
  !
  use mpi
#endif

  ! 宣言文 ; Declaration statements
  !
  implicit none

  ! 公開手続き
  ! Public procedure
  !
  public :: MPIWrapperInit
  public :: MPIWrapperCartCreate
  public :: MPIWrapperCartShift
  public :: MPIWrapperCommFree
  public :: MPIWrapperFinalize
  public :: MPIWrapperISend
  public :: MPIWrapperIRecv
  public :: MPIWrapperWait
  public :: MPIWrapperGather
  public :: MPIWrapperAllreduce

  ! 非公開変数
  !
  private
!  integer, save :: mpi_comm_cart_save
!  logical, save :: FLAG_MPI_Cart = .false. 
!                           ! FLAG

  ! 公開変数
  ! Public variables
  !
  integer, save, public :: nprocs
                           ! Number of MPI processes
  integer, save, public :: myrank
                           ! My rank
  logical, save, public :: FLAG_LIB_MPI
                           ! FLAG

  ! 非公開変数
  ! Private variables
  !
  interface MPIWrapperISend
    module procedure &
      MPIWrapperISend_dble_1d, &
      MPIWrapperISend_dble_2d, &
      MPIWrapperISend_dble_3d, &
      MPIWrapperISend_dble_4d
  end interface

  interface MPIWrapperIRecv
    module procedure &
      MPIWrapperIRecv_dble_1d, &
      MPIWrapperIRecv_dble_2d, &
      MPIWrapperIRecv_dble_3d, &
      MPIWrapperIRecv_dble_4d
  end interface

  interface MPIWrapperAbort
    module procedure &
      MPIWrapperStop
  end interface

contains

  !------------------------------------------------------------------------------------

  subroutine MPIWrapperInit
    !
    ! MPI の初期化
    !
    ! Initialization of MPI
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr

#endif

    ! set default values
    ! 
    FLAG_LIB_MPI = .false.
    nprocs = 1
    myrank = 0

#ifdef LIB_MPI

    ! call MPI initial subroutine
    ! 
    FLAG_LIB_MPI = .true.
    call mpi_init( ierr )
    call mpi_comm_size( mpi_comm_world, nprocs, ierr )
    call mpi_comm_rank( mpi_comm_world, myrank, ierr )

#endif

  end subroutine MPIWrapperInit

  !------------------------------------------------------------------------------------

  subroutine MPIWrapperCartCreate( xsub, ysub, flag_periodic, mpi_comm_cart )
    !
    ! MPI cart の初期化
    !
    ! Initialization of MPI cart
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    !
    integer, intent(in)    :: xsub, ysub  ! CPU nums of X and Y direction
    logical, intent(in)    :: flag_periodic
    integer, intent(inout) :: mpi_comm_cart

    ! 作業変数
    !
    integer, parameter :: disp = 1
    integer, parameter :: ndim = 2
    logical, parameter :: reorder = .false.
    logical            :: periodic(ndim)
    integer            :: idivid(ndim)
    integer            :: direction
    integer            :: ierr

#ifdef LIB_MPI
      
    idivid(1) = ysub
    idivid(2) = xsub
    periodic(1) = flag_periodic
    periodic(2) = flag_periodic
    
    call mpi_cart_create(                                 &
      &  mpi_comm_world, ndim, idivid, periodic, reorder, & ! (IN)
      &  mpi_comm_cart, ierr                              & ! (OUT)
      & )    

!    if (.not. Flag_MPI_Cart) then 
!      Flag_MPI_Cart = .true.
!      mpi_comm_cart_save = mpi_comm_cart
!    end if

#endif
  end subroutine MPIWrapperCartCreate

  !------------------------------------------------------------------------------------

  subroutine MPIWrapperCommFree(mpi_comm_cart)
    !
    ! MPI の終了処理
    !
    ! Finalization of MPI
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数 
    !
    integer, intent(in) :: mpi_comm_cart

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: cart
    integer :: ierr

    cart = mpi_comm_cart
    call mpi_comm_free( cart, ierr )

#endif

  end subroutine MPIWrapperCommFree

  !-------------------------------------------------------------------------------------
  subroutine MPIWrapperFinalize
    !
    ! MPI の終了処理
    !
    ! Finalization of MPI
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr

!    if (FLAG_MPI_Cart) then 
!      call mpi_comm_free( mpi_comm_cart_save, ierr )
!    end if
    call mpi_finalize( ierr )

#endif

  end subroutine MPIWrapperFinalize

  !-------------------------------------------------------------------------------------
  subroutine MPIWrapperCartShift(mpi_comm_cart, direction, disp, rank1, rank2)

    integer, intent(in)   :: mpi_comm_cart
    integer, intent(in)   :: direction 
    integer, intent(in)   :: disp
    integer, intent(out)  :: rank1
    integer, intent(out)  :: rank2
    integer               :: ierr

    ! デフォルト値. 
    !
    rank1 = -1
    rank2 = -1

#ifdef LIB_MPI

    call mpi_cart_shift(                  &
      &  mpi_comm_cart, direction, disp,  & !(IN)
      &  rank1, rank2, ierr               & !(OUT)
      & )

#endif
    
  end subroutine MPIWrapperCartShift

  !------------------------------------------------------------------------------------

  subroutine MPIWrapperStop
    !
    ! MPI の異常終了処理
    !
    ! Abort of MPI
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: errorcode = 9
    integer :: ierr

    call mpi_abort( mpi_comm_world, errorcode, ierr )
    call MPIWrapperFinalize
    stop

#endif 

  end subroutine MPIWrapperstop

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperWait( ireq )
    !
    ! MPI 通信終了まで待機
    !
    ! Wait finishing MPI transfer
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 作業変数
    ! Work variables
    !
    integer, intent(inout) :: ireq ! request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: istatus( MPI_STATUS_SIZE )

    call mpi_wait( ireq, istatus, ierr )

#endif 

  end subroutine MPIWrapperWait

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperGather( &
    & nz, z_Var,               & ! (in)
    & a_Var                    & ! (out)
    & )
    !
    ! 1D 倍精度配列の非ブロッキング通信(送信)
    !
    ! Non-blocking transfer (send) of real(8) 1D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: nz
    real(8) , intent(in ) :: z_Var(nz)
    real(8) , intent(out) :: a_Var(nz * nprocs)


#ifdef LIB_MPI
    ! 作業変数
    ! Work variables
    !
    integer , parameter   :: root = 0
    integer               :: i, ierr

    i = myrank * nz + 1

    call mpi_gather( z_Var, nz, mpi_double_precision,    &
      &              a_Var(i), nz, mpi_double_precision, &
      &              root, mpi_comm_world, ierr)

#endif 

  end subroutine MPIWrapperGather


  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperAllreduce( &
    &   nz,  sbuf,                & ! (in)
    &   rbuf                      & ! (out)
    & )
    !
    ! 集団通信サブルーチン
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: nz
    real(8) , intent(in ) :: sbuf(nz)
    real(8) , intent(out) :: rbuf(nz)

    ! 作業変数
    ! Work variables
    !
    integer               :: ierr

    ! シングル版対応.
    !
    rbuf = sbuf   

#ifdef LIB_MPI

    call mpi_allreduce(                   &
      &    sbuf, rbuf, nz,                &
      &    mpi_double_precision, MPI_SUM, &
      &    mpi_comm_world, ierr  )

#endif 

  end subroutine MPIWrapperAllreduce
  

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperISend_dble_1d( &
    & idest, im,                      & ! (in)
    & buf,                            & ! (in)
    & ireq                            & ! (out)
    & )
    !
    ! 1D 倍精度配列の非ブロッキング通信(送信)
    !
    ! Non-blocking transfer (send) of real(8) 1D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idest
                              ! Process number of destination
    integer , intent(in ) :: im
                              ! Size of 1st dimension of sent data
    real(DP), intent(in ) :: buf( im )
                              ! Array to be sent
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI
    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_isend( buf, isize, &
      mpi_double_precision, idest, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperISend_dble_1d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperIRecv_dble_1d( &
    & idep, im,                       & ! (in)
    & buf,                            & ! (out)
    & ireq                            & ! (out)
    & )
    !
    ! 1D 倍精度配列の非ブロッキング通信(受信)
    !
    ! Non-blocking transfer (receive) of real(8) 1D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idep
                              ! Process number of departure
    integer , intent(in ) :: im
                              ! Size of 1st dimension of received data
    real(DP), intent(out) :: buf( im )
                              ! Array to be received
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_irecv( buf, isize, &
      mpi_double_precision, idep, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperIRecv_dble_1d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperISend_dble_2d( &
    & idest, im, jm,                  & ! (in)
    & buf,                            & ! (in)
    & ireq                            & ! (out)
    & )
    !
    ! 2D 倍精度配列の非ブロッキング通信(送信)
    !
    ! Non-blocking transfer (send) of real(8) 2D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idest
                              ! Process number of destination
    integer , intent(in ) :: im
                              ! Size of 1st dimension of sent data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of sent data
    real(DP), intent(in ) :: buf( im, jm )
                              ! Array to be sent
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_isend( buf, isize, &
      mpi_double_precision, idest, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperISend_dble_2d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperIRecv_dble_2d( &
    & idep, im, jm,                   & ! (in)
    & buf,                            & ! (out)
    & ireq                            & ! (out)
    & )
    !
    ! 2D 倍精度配列の非ブロッキング通信(受信)
    !
    ! Non-blocking transfer (receive) of real(8) 2D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idep
                              ! Process number of destination
    integer , intent(in ) :: im
                              ! Size of 1st dimension of received data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of received data
    real(DP), intent(out) :: buf( im, jm )
                              ! Array to be received
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_irecv( buf, isize, &
      mpi_double_precision, idep, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperIRecv_dble_2d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperISend_dble_3d( &
    & idest, im, jm, km,              & ! (in)
    & buf,                            & ! (in)
    & ireq                            & ! (out)
    & )
    !
    ! 3D 倍精度配列の非ブロッキング通信(送信)
    !
    ! Non-blocking transfer (send) of real(8) 3D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idest
                              ! Process number of destination
    integer , intent(in ) :: im
                              ! Size of 1st dimension of sent data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of sent data
    integer , intent(in ) :: km
                              ! Size of 3rd dimension of sent data
    real(DP), intent(in ) :: buf( im, jm, km )
                              ! Array to be sent
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

!    write(*,*) "*** iSEND *** ", myrank, " => ", idest

    isize = size( buf )

    call mpi_isend( buf, isize, &
      mpi_double_precision, idest, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperISend_dble_3d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperIRecv_dble_3d( &
    & idep, im, jm, km,               & ! (in)
    & buf,                            & ! (out)
    & ireq                            & ! (out)
    & )
    !
    ! 3D 倍精度配列の非ブロッキング通信(受信)
    !
    ! Non-blocking transfer (receive) of real(8) 3D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idep
                              ! Process number of departure
    integer , intent(in ) :: im
                              ! Size of 1st dimension of received data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of received data
    integer , intent(in ) :: km
                              ! Size of 3rd dimension of received data
    real(DP), intent(inout) :: buf( im, jm, km )
                              ! Array to be received
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

!    write(*,*) "*** iRECV *** ", myrank, " <= ", idep

    call mpi_irecv( buf, isize, &
      mpi_double_precision, idep, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperIRecv_dble_3d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperISend_dble_4d( &
    & idest, im, jm, km, lm,          & ! (in)
    & buf,                            & ! (in)
    & ireq                            & ! (out)
    & )
    !
    ! 4D 倍精度配列の非ブロッキング通信(送信)
    !
    ! Non-blocking transfer (send) of real(8) 4D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idest
                              ! Process number of destination
    integer , intent(in ) :: im
                              ! Size of 1st dimension of sent data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of sent data
    integer , intent(in ) :: km
                              ! Size of 3rd dimension of sent data
    integer , intent(in ) :: lm
                              ! Size of 4th dimension of sent data
    real(DP), intent(in ) :: buf( im, jm, km, lm )
                              ! Array to be sent
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_isend( buf, isize, &
      mpi_double_precision, idest, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperISend_dble_4d

  !--------------------------------------------------------------------------------------

  subroutine MPIWrapperIRecv_dble_4d( &
    & idep, im, jm, km, lm,           & ! (in)
    & buf,                            & ! (out)
    & ireq                            & ! (out)
    & )
    !
    ! 4D 倍精度配列の非ブロッキング通信(受信)
    !
    ! Non-blocking transfer (receive) of real(8) 4D array
    !

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    ! input/output variables
    !
    integer , intent(in ) :: idep
                              ! Process number of departure
    integer , intent(in ) :: im
                              ! Size of 1st dimension of received data
    integer , intent(in ) :: jm
                              ! Size of 2nd dimension of received data
    integer , intent(in ) :: km
                              ! Size of 3rd dimension of received data
    integer , intent(in ) :: lm
                              ! Size of 4th dimension of received data
    real(DP), intent(out) :: buf( im, jm, km, lm )
                              ! Array to be received
    integer , intent(out) :: ireq
                              ! Request number

#ifdef LIB_MPI

    ! 作業変数
    ! Work variables
    !
    integer :: ierr
    integer :: isize

    isize = size( buf )

    call mpi_irecv( buf, isize, &
      mpi_double_precision, idep, 1, mpi_comm_world, &
      ireq, ierr )

#endif 

  end subroutine MPIWrapperIRecv_dble_4d

  !--------------------------------------------------------------------------------------

end module mpi_wrapper
