!= �֤����ΰ�פ��ͤ��������뤿��Υ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: setmargin2.f90,v 1.3 2014/03/04 05:55:07 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module SetMargin
  !
  ! ͭ�º�ʬ��ǥ��� �������⥸�塼��
  ! X ������Y �����ˤ� MPI ���Ѥ������󲽤���ǽ
  !

  !���ۤη�����ػ�
  implicit none

  !°��
  private

  integer, save, private :: Rrank = 1   ! ���٤� CPU �� rank
  integer, save, private :: Lrank = 1   ! ���٤� CPU �� rank
  integer, save, private :: Urank = 1   ! ��¦�� CPU �� rank
  integer, save, private :: Drank = 1   ! ��¦�� CPU �� rank

  !�ؿ��θ���
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

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine XcyclYcycl(aaa_Var)
    !
    ! X ������ Y �������Ф��Ƽ�����������Ŭ��
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only : DP
    use gridset,     only : imin, imax, jmin, jmax, kmin, kmax,    &
      &                     xmg, ymg, zmg, nx, ny, nz, FlagCalc3D
    use mpi_wrapper, only : MPIWrapperISend, MPIWrapperIRecv, MPIWrapperWait

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP),intent(inout) :: aaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: Xsbuf_a( Xmg, jmin:jmax, kmin:kmax )
    real(DP)               :: Xrbuf_a( Xmg, jmin:jmax, kmin:kmax )
    real(DP)               :: Xsbuf_b( Xmg, jmin:jmax, kmin:kmax )
    real(DP)               :: Xrbuf_b( Xmg, jmin:jmax, kmin:kmax )
    real(8), allocatable   :: Ysbuf_a( :,:,: )
    real(8), allocatable   :: Yrbuf_a( :,:,: )
    real(8), allocatable   :: Ysbuf_b( :,:,: )
    real(8), allocatable   :: Yrbuf_b( :,:,: )
!    real(DP)               :: Ysbuf_a( nx, ymg, nz )
!    real(DP)               :: Yrbuf_a( nx, ymg, nz )
!    real(DP)               :: Ysbuf_b( nx, ymg, nz )
!    real(DP)               :: Yrbuf_b( nx, ymg, nz )
    integer                :: ireqs_a, ireqr_a, ireqs_b, ireqr_b
    integer                :: ny2, nz2

    !!!
    !!! Y �������̿�����
    !!!   single CPU ��ư��������, MPIWeapper* ���֥롼����ϲ��⤷�ʤ����Ȥ����
    !!!
    if ( FlagCalc3D ) then 
      
      allocate( Ysbuf_a( nx, ymg, nz ) )
      allocate( Yrbuf_a( nx, ymg, nz ) )
      allocate( Ysbuf_b( nx, ymg, nz ) )
      allocate( Yrbuf_b( nx, ymg, nz ) )

      !-------------------------------
      ! ����������ʬ����(sbuf_a, sbuf_b)���Ѱդ���
      !
      Ysbuf_a( 1:nx, 1:ymg, 1:nz ) = aaa_Var( 1:nx, ny-ymg+1:ny, 1:nz )
      Ysbuf_b( 1:nx, 1:ymg, 1:nz ) = aaa_Var( 1:nx, 1:ymg, 1:nz )
      
      !-------------------------------
      ! ����������
      !
      ! ���֥롼����� call ������֤� 2CPU �ξ����θ���Ʒ���ɬ�פ�����
      ! rank0 �ˤƾ�¦����¦�ν����������
      !  => rank0 �� rank1 ���鲼¦����¦�ν�˥ǡ�������Ǥ���

      Yrbuf_a = Ysbuf_b  ! single CPU ���Ѥ���Ȥ��Τ���ν���
      Yrbuf_b = Ysbuf_a  ! single CPU ���Ѥ���Ȥ��Τ���ν���
      
      call MPIWrapperISend( Urank, nx, ymg, nz, Ysbuf_a, ireqs_a ) !����
      call MPIWrapperISend( Drank, nx, ymg, nz, Ysbuf_b, ireqs_b ) !����

      call MPIWrapperIRecv( Drank, nx, ymg, nz, Yrbuf_b, ireqr_b ) !����      
      call MPIWrapperIRecv( Urank, nx, ymg, nz, Yrbuf_a, ireqr_a ) !����
      
      !-------------------------------
      ! �������ν�λ���ԤäƤ�������
      
      call MPIWrapperWait( ireqs_a )       
      call MPIWrapperWait( ireqr_a )
      call MPIWrapperWait( ireqs_b )
      call MPIWrapperWait( ireqr_b )

      aaa_var( 1:nx, ny+1:ny+ymg, 1:nz ) = Yrbuf_a( 1:nx, 1:ymg, 1:nz ) 
      aaa_var( 1:nx, -ymg+1:0, 1:nz ) = Yrbuf_b( 1:nx, 1:ymg, 1:nz ) 

    end if
    
    !!!
    !!! X �����˼�����������Ŭ��
    !!!   single CPU ��ư��������, MPIWeapper* ���֥롼����ϲ��⤷�ʤ����Ȥ����
    !!!
    !-------------------------------
    ! ����������ʬ����(sbuf_a, sbuf_b)���Ѱդ���
    Xsbuf_a( 1:xmg, jmin:jmax, kmin:kmax ) = aaa_Var( nx-xmg+1:nx, jmin:jmax, kmin:kmax )  ! nx - xmg + 1 ~ nx ���ݴ�
    Xsbuf_b( 1:xmg, jmin:jmax, kmin:kmax ) = aaa_Var( 1:xmg, jmin:jmax, kmin:kmax )  ! 1 ~ xmg ���ݴ�

    !-------------------------------
    ! ����������
    !
    ! ���֥롼����� call ������֤� 2CPU �ξ����θ���Ʒ���ɬ�פ�����
    ! rank0 �ˤƱ�¦����¦�ν����������
    !  => rank0 �� rank1 ��꺸¦����¦�ν�˥ǡ��������

    Xrbuf_a = Xsbuf_b  ! single CPU ���Ѥ���Ȥ��Τ���ν���
    Xrbuf_b = Xsbuf_a  ! single CPU ���Ѥ���Ȥ��Τ���ν���

    ny2 = ny + ymg * 2
    nz2 = nz + zmg * 2

    call MPIWrapperISend( Rrank, xmg, ny2, nz2, Xsbuf_a, ireqs_a ) !����
    call MPIWrapperISend( Lrank, xmg, ny2, nz2, Xsbuf_b, ireqs_b ) !����
    
    call MPIWrapperIRecv( Lrank, xmg, ny2, nz2, Xrbuf_b, ireqr_b ) !����    
    call MPIWrapperIRecv( Rrank, xmg, ny2, nz2, Xrbuf_a, ireqr_a ) !����
   
    !-------------------------------
    ! �������ν�λ���ԤäƤ�������

    call MPIWrapperWait( ireqs_a )    
    call MPIWrapperWait( ireqr_a )
    call MPIWrapperWait( ireqs_b )       
    call MPIWrapperWait( ireqr_b )

    aaa_var( nx+1:nx+xmg, jmin:jmax, kmin:kmax ) = Xrbuf_a( 1:xmg, jmin:jmax, kmin:kmax ) ! nx + 1 ~ nx + xmg ������
    aaa_var( -xmg+1:0, jmin:jmax, kmin:kmax ) = Xrbuf_b( 1:xmg, jmin:jmax, kmin:kmax ) ! - xmg + 1 ~ 0 ������
    
  end subroutine XcyclYcycl

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ZSym(aaz_Var)
    !
    ! z �������оζ�������Ŭ�Ѥ���
    !

    use dc_types, only : DP
    use gridset,  only : imin, imax, jmin, jmax, kmin, kmax, zmg, nz
    
    implicit none
    
    real(DP),intent(inout) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    integer                :: k

    do k = 1, zmg
      aaz_Var(:,:,1-k)  = aaz_Var(:,:,k)
      aaz_Var(:,:,nz+k) = aaz_Var(:,:,nz+1-k)
    end do
    
  end subroutine ZSym

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ZAntSym(aar_Var)
    !
    ! z �������оζ�������Ŭ�Ѥ���
    !

    use dc_types, only : DP
    use gridset,  only : imin, imax, jmin, jmax, kmin, kmax, zmg, nz

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

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine SetMargin_xyzf(xyzf_Var)

    use dc_types, only : DP
    use gridset,  only : imin, imax, jmin, jmax, kmin, kmax, ncmax

    implicit none

    real(DP),intent(inout) :: xyzf_Var(imin:imax,jmin:jmax,kmin:kmax,1:ncmax) 
    integer                :: s

    do s = 1, ncmax

      ! x, y �����˼�����������Ŭ�Ѥ���
      !
      call XCyclYCycl( xyzf_Var(:,:,:,s) ) !inout

      ! z �������оݶ�������Ŭ�Ѥ���
      !
      call ZSym( xyzf_Var(:,:,:,s) ) !inout
      
    end do

  end subroutine SetMargin_xyzf

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine SetMargin_aaz(aaz_Var)

    use dc_types, only : DP
    use gridset,  only : imin, imax, jmin, jmax, kmin, kmax

    implicit none

    real(DP),intent(inout) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 


    ! x, y �����˼�����������Ŭ�Ѥ���
    !
    call XCyclYCycl( aaz_Var ) !inout
    
    ! z �������оݶ�������Ŭ�Ѥ���
    !
    call ZSym( aaz_Var ) !inout
       
  end subroutine SetMargin_aaz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  subroutine SetMargin_aar(aar_Var)


    use dc_types, only : DP
    use gridset,  only : imin, imax, jmin, jmax, kmin, kmax

    implicit none

    real(DP),intent(inout) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax) 

    ! x, y �����˼�����������Ŭ�Ѥ���
    !
    call XCyclYCycl( aar_Var ) !inout
    
    ! z �������оݶ�������Ŭ�Ѥ���
    !
    call ZAntSym( aar_Var ) !inout
    
    
  end subroutine SetMargin_aar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  subroutine SetMargin_init
    !
    ! �����
    !
    ! Initialization of SetMargin 
    !

    ! �⥸�塼����� ; USE statements
    !
    use dc_message,  only : MessageNotify
    use gridset,     only : xsub, ysub
    use mpi_wrapper, only : MPIWrapperCartCreate, &
      &                     MPIWrapperCartShift,  &
      &                     MPIWrapperCommFree
    
    ! ����ѿ�
    ! Work variables
    !
    integer, parameter :: disp = 1
    integer            :: direction
    integer            :: comm_cart
    logical, parameter :: periodic = .true.

    ! �̿���μ���   
    !
    call MPIWrapperCartCreate(xsub, ysub, periodic, comm_cart)
    direction = 0
    call MPIWrapperCartShift(comm_cart, direction, disp, Urank, Drank)
    direction = 1
    call MPIWrapperCartShift(comm_cart, direction, disp, Lrank, Rrank)
    call MPIWrapperCommFree(comm_cart)

    call MessageNotify( "M", "SetMargin_init", "Rrank = %d",   i=(/Rrank/) )
    call MessageNotify( "M", "SetMargin_init", "Lrank = %d",   i=(/Lrank/) )
    call MessageNotify( "M", "SetMargin_init", "Urank = %d",   i=(/Urank/) )
    call MessageNotify( "M", "SetMargin_init", "Drank = %d",   i=(/Drank/) )
    
  end subroutine SetMargin_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

end module SetMargin
