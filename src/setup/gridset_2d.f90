!= �ʻ������ѥ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: gridset.f90,v 1.14 2014/07/08 01:05:33 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module gridset
  !
  != �ʻ������ѥ⥸�塼��
  !
  !������Ϳ����줿 NAMELIST �ե����뤫��, �ʻ�������������, 
  !�ݴɤ��뤿����ѿ����ȷ��⥸�塼��
  !
  ! 2D �Ǥ� 3D �Ǥ� RegXMin, RegYMin, RegZMin ��������ۤʤ뤳�Ȥ����
  !
  !== Procedures List
  ! gridset_init   :: ������롼����
  ! gridset_check  :: �������Ƥ�������������å����뤿��Υ롼����
  !

  !���ۤη�����ػ�
  implicit none
  
  !save °��
  private
  
  !�����ѿ�
  integer, public, save :: NX = 10 ! x �����ʻ�����
  integer, public, save :: NY = 10 ! y �����ʻ�����
  integer, public, save :: NZ = 10 ! z �����ʻ�����
  integer, public, save :: NCMAX = 1  ! �����������ǿ�
  integer, public, save :: Xmg = 2 ! x ��������ʻ�����
  integer, public, save :: Ymg = 2 ! y ��������ʻ�����
  integer, public, save :: Zmg = 2 ! z ��������ʻ�����
  integer, public, save :: imin    ! x ����������β��� 
  integer, public, save :: imax    ! x ����������ξ��
  integer, public, save :: jmin    ! y ����������β��� 
  integer, public, save :: jmax    ! y ����������β��� 
  integer, public, save :: kmin    ! z ����������β��� 
  integer, public, save :: kmax    ! z ����������β��� 
  integer, public, save :: xsub = 1
                                   ! Number of MPI processes (X direction)
  integer, public, save :: ysub = 1
                                   ! Number of MPI processes (Y direction)
  logical, public, save :: FlagCalc3D = .true.
  logical, public, save :: FlagCalcMoist = .true.

  integer, private, save :: xdim = 1 
  integer, private, save :: ydim = 1
  integer, private, save :: zdim = 1
  integer, public, save  :: nxny = 100 ! x ������ y ��������ʻ�����

  integer, public, save :: MarginX       !�����Υ���åɿ�
  integer, public, save :: MarginY       !�����Υ���åɿ�
  integer, public, save :: MarginZ       !�����Υ���åɿ�
  integer, public, save :: SpcNum        !���ؼ�ο�
  integer, public, save :: DimXMin       ! x ����������β���
  integer, public, save :: DimXMax       ! x ����������ξ��
  integer, public, save :: DimYMin       ! y ����������β���
  integer, public, save :: DimYMax       ! y ����������ξ��
  integer, public, save :: DimZMin       ! z ����������β���
  integer, public, save :: DimZMax       ! z ����������ξ��
  integer, public, save :: RegXMin       ! x ������ʪ���ΰ�β���
  integer, public, save :: RegXMax       ! x ������ʪ���ΰ�ξ��
  integer, public, save :: RegYMin       ! y ������ʪ���ΰ�β���
  integer, public, save :: RegYMax       ! y ������ʪ���ΰ�ξ��
  integer, public, save :: RegZMin       ! z ������ʪ���ΰ�β���
  integer, public, save :: RegZMax       ! z ������ʪ���ΰ�ξ��
  integer, public, save :: FileNX        !�ե����������
  integer, public, save :: FileNY        !�ե����������
  integer, public, save :: FileNZ        !�ե����������
  integer, public, save :: FileXMin      !�ե����������
  integer, public, save :: FileXMax      !�ե����������
  integer, public, save :: FileYMin      !�ե����������
  integer, public, save :: FileYMax      !�ե����������
  integer, public, save :: FileZMin      !�ե����������
  integer, public, save :: FileZMax      !�ե����������

  integer, public, save :: im = 10 ! x �����ʻ�����
  integer, public, save :: jm = 10 ! y �����ʻ�����
  integer, public, save :: km = 10 ! z �����ʻ�����

  ! ���֥롼����θ���
  public gridset_init

contains

  subroutine gridset_init
    !
    != ������롼����
    !
    ! ����ե����뤫�������ɤ߹��߳ʻ�������׻�����
    !

    !�⥸�塼���ɤ߹���
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use mpi_wrapper,   only: nprocs
    use namelist_util, only: namelist_filename
    
    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    integer            :: unit                !����ե������������ֹ�

    !-----------------------------------------------------------------
    ! ����ե����뤫�������ɤ߹���
    !
    NAMELIST /gridset_nml/ xdim, ydim, zdim, NCMAX, Xmg, Ymg, Zmg, xsub, ysub

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=gridset_nml)
    close(unit)

    !-----------------------------------------------------------------    
    ! NX, NY, NZ �����
    !
    nx = xdim / xsub 
    imin = 0  - xmg   ! 3D �ǤȰۤʤ�
    imax = nx + xmg
   
    if (ydim == 1) then
      ! �󼡸�
      ymg = 0
      ny = 1
      FlagCalc3D = .false. 
    else
      ! ������
      ny = ydim / ysub
      FlagCalc3D = .true. 
    end if
    jmin = 0  - ymg   ! 3D �ǤȰۤʤ�
    jmax = ny + ymg

    nz = zdim
    kmin = 0  - zmg   ! 3D �ǤȰۤʤ�
    kmax = nz + zmg    

    nxny = nx * ny

    ! �ȥ졼�����׻��Υ����å�
    !
    if ( ncmax == 0 ) then 
!      FlagCalcTracer = .false.
      FlagCalcMoist = .false.
    else
!      FlagCalcTracer = .true.
      FlagCalcMoist = .true.
    end if

    !-----------------------------------------------------------------    
    ! �ɤ߹������������
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
    ! �ͤΥ����å�
    !    
    call gridset_check

    !-----------------------------------------------------------------
    ! arare4 �Ѥ��ѿ�������
    !
    im = nx
    jm = ny
    km = nz

    MarginX = Xmg
    MarginY = Ymg
    MarginZ = Zmg
    SpcNum  = ncmax

    RegXMin = 0    !3D �ǤȰۤʤ�
    RegXMax = NX
    RegYMin = 0    !3D �ǤȰۤʤ�
    RegYMax = NY
    RegZMin = 0    !3D �ǤȰۤʤ�
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

    FileXMin = RegXMin + 1    !3D �ǤȰۤʤ�
    FileXMax = RegXMax
    FileYMin = RegYMin + 1    !3D �ǤȰۤʤ�
    FileYMax = RegYMax
    FileZMin = RegZMin + 1    !3D �ǤȰۤʤ�
    FileZMax = RegZMax

  contains

    subroutine gridset_check
      !
      != �������Ƥ�������������å����뤿��Υ롼����
      !
      ! �������Ƥ�̷�⤬�ʤ�����ǧ�򤹤�. 
      ! ���꤬ͭ����ϥ��顼��å�������Ф��ƽ�λ����.
      ! 

      !���ۤη�����ػ�
      implicit none
      
      ! nprocs �� xsub or ysub �ǳ���ڤ�뤫�����å�. 
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
      
      ! X �����Υޡ�������礭���Υ����å�
      ! 
      if ( xdim < xmg ) then
        call MessageNotify( "E", "gridset_init", "xdim < Xmg" )
      end if
      
      ! xdim �� xsub �ǳ���ڤ�뤫��ǧ
      ! 
      if ( mod(xdim, xsub) /= 0  ) then
        call MessageNotify( "E", "gridset_init", "mod(xdim, xsub) /= 0" )
      end if
      
      ! Y �����Υޡ�������礭���Υ����å�
      ! 
      if ( ydim < ymg ) then
        call MessageNotify( "E", "gridset_init", "ydim < Ymg" )
      end if
      
      ! ydim �� xsub �ǳ���ڤ�뤫��ǧ
      ! 
      if ( mod(ydim, ysub) /= 0  ) then
        call MessageNotify( "E", "gridset_init", "mod(ydim, ysub) /= 0" )
      end if
      
      ! Z �����Υޡ�������礭���Υ����å�
      ! 
      if ( zdim < zmg ) then
        call MessageNotify( "E", "gridset_init", "zdim < Zmg" )
      end if
      
    end subroutine gridset_check

  end subroutine gridset_init
  
end module gridset
