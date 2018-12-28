!= 3 ���� (xyz ����) ���ֳָ�߳ʻ� �ʻ�������⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: axesset.f90,v 1.16 2014/07/08 01:05:32 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module axesset
  != 3 ���� (xyz ����) ���ֳָ�߳ʻ� �ʻ�������⥸�塼��
  !
  !== ����
  !
  ! axesset ��, 3 ���� (xyz ����) ���ֳָ�߳ʻҤ��Ѥ���ͭ�º�ʬˡ�˴�Ť�
  ! ��������󶡤���. 
  ! 

  !�⥸�塼���ɤ߹���
  !
  use dc_types,      only: DP

  !���ۤη�����ػ�
  !
  implicit none

  !�ǥե���Ȥ������
  !
  private

  !��������
  ! �������³��
  !
  public :: axesset_init

  ! �ѿ����
  real(DP), public, save :: Xmin = 0.0d0                ! x ��ɸ�λ���������
  real(DP), public, save :: Xmax = 1.0d4                ! x ��ɸ�λ���������
  real(DP), public, save :: Ymin = 0.0d0                ! x ��ɸ�λ���������
  real(DP), public, save :: Ymax = 1.0d4                ! x ��ɸ�λ���������
  real(DP), public, save :: Zmin = 0.0d0                ! z ��ɸ�λ���������
  real(DP), public, save :: Zmax = 1.0d4                ! z ��ɸ�λ���������
  real(DP), public, save :: DX  
  real(DP), public, save :: DY           
  real(DP), public, save :: DZ           
  real(DP), allocatable, public, save :: x_X(:)         ! Ⱦ�����ʻ�����ɸ
  real(DP), allocatable, public, save :: p_X(:)         ! �����ʻ�����ɸ
  real(DP), allocatable, public, save :: x_dx(:)        ! Ⱦ�����ʻ����ֳ�
  real(DP), allocatable, public, save :: p_dx(:)        ! �����ʻ����ֳ�
  real(DP), allocatable, public, save :: y_Y(:)         ! Ⱦ�����ʻ�����ɸ
  real(DP), allocatable, public, save :: q_Y(:)         ! �����ʻ�����ɸ
  real(DP), allocatable, public, save :: y_dy(:)        ! Ⱦ�����ʻ����ֳ�
  real(DP), allocatable, public, save :: q_dy(:)        ! �����ʻ����ֳ�
  real(DP), allocatable, public, save :: z_Z(:)         ! Ⱦ�����ʻ�����ɸ
  real(DP), allocatable, public, save :: r_Z(:)         ! �����ʻ�����ɸ
  real(DP), allocatable, public, save :: z_dz(:)        ! Ⱦ�����ʻ����ֳ�
  real(DP), allocatable, public, save :: r_dz(:)        ! �����ʻ����ֳ�
  real(DP), allocatable, public, save :: xyz_X(:,:,:)   ! x ��ɸ(Ⱦ�����ʻ�)
  real(DP), allocatable, public, save :: xyz_Y(:,:,:)   ! y ��ɸ(Ⱦ�����ʻ�)
  real(DP), allocatable, public, save :: xyz_Z(:,:,:)   ! z ��ɸ(Ⱦ�����ʻ�)
  real(DP), allocatable, public, save :: xyz_dX(:,:,:)  ! x �ʻҴֳ�(Ⱦ�����ʻ�)
  real(DP), allocatable, public, save :: xyz_dY(:,:,:)  ! y �ʻҴֳ�(Ⱦ�����ʻ�)
  real(DP), allocatable, public, save :: xyz_dZ(:,:,:)  ! z �ʻҴֳ�(Ⱦ�����ʻ�)

  real(DP), public, save :: DelX                        !�ʻ����ֳ�
  real(DP), public, save :: DelY                        !�ʻ����ֳ�
  real(DP), public, save :: DelZ                        !�ʻ����ֳ�
  real(DP), allocatable, public, save :: s_X(:)         !X ��ɸ��(�����顼�ʻ���)
  real(DP), allocatable, public, save :: f_X(:)         !X ��ɸ��(�٥��ȥ�ʻ���)
  real(DP), allocatable, public, save :: s_Y(:)         !Y ��ɸ��(�����顼�ʻ���)
  real(DP), allocatable, public, save :: f_Y(:)         !Y ��ɸ��(�٥��ȥ�ʻ���)
  real(DP), allocatable, public, save :: s_Z(:)         !Z ��ɸ��(�����顼�ʻ���)
  real(DP), allocatable, public, save :: f_Z(:)         !Z ��ɸ��(�٥��ȥ�ʻ���)

contains
  !--------------------------------------------------------------------
  subroutine axesset_init
    !
    ! �ʻ�����ɸ����ȳʻ����ֳ�����ν����    

    ! �⥸�塼��ƤӽФ�
    !
    use dc_types,      only: DP
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use mpi_wrapper,   only: MPIWrapperCartCreate, &
      &                      MPIWrapperCartShift,  &
      &                      MPIWrapperCommFree
    use gridset,       only: FlagCalc3D,    &
      &                      xsub, ysub,    & ! 
      &                      nx, ny, nz,    & ! �ʻ�����
      &                      imin, imax, jmin, jmax, kmin, kmax ! ����ξ���ͤȲ�����
    use namelist_util, only: namelist_filename

    ! ���ۤη�����ػ�
    !
    implicit none

    ! �ѿ����
    !
    real(DP),allocatable :: xy_X(:,:)! x ��ɸ(Ⱦ�����ʻ�, �������)
    real(DP),allocatable :: xy_Y(:,:)! y ��ɸ(Ⱦ�����ʻ�, �������)
    real(DP),allocatable :: yz_Z(:,:)! z ��ɸ(Ⱦ�����ʻ�, �������)
    integer              :: unit     ! ����ե������������ֹ�
    integer              :: comm_cart
    logical, parameter   :: periodic = .false.

    !����ե����뤫���ɤ߹�����ϥե��������
    NAMELIST /axesset_nml/ xmin, xmax, ymin, ymax, zmin, zmax

    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=axesset_nml)
    close(unit)
    
    ! ����ξ岼�¤���, ��ɸ�ͤȳʻ����ֳ֤�����
    ! * 1 �����ѤΥ��֥롼������Ѥ���
    !
    call MPIWrapperCartCreate(xsub, ysub, periodic, comm_cart)
    call x_axis_init(comm_cart)
    call y_axis_init(comm_cart)
    call z_axis_init
    call MPIWrapperCommFree(comm_cart)
    
    ! 3 �����ʻ�����ɸ���������
    ! * �Ȥ߹��ߴؿ� spread ���Ѥ���. 
    ! * �������Ȥ��� 2 �����ʻ�����ɸ�������, ����� 3 �����˳�ĥ����.
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

    ! �ɤ߹������������
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
    ! arare4 �Ѥ�����
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
      != �ʻ�����ɸ����ȳʻ����ֳ�����ν����    
      !
      
      ! ���ۤη�����ػ�
      implicit none
      
      integer, intent(in) :: comm_cart
      
      ! ����ѿ�
      integer             :: ix   ! do �롼��ź��
      integer             :: rank, tmp
      integer             :: lp
      integer, parameter  :: direction = 1
      
      allocate(x_X(imin:imax))
      allocate(p_X(imin:imax))
      allocate(x_dx(imin:imax))
      allocate(p_dx(imin:imax))
      
      ! �����
      lp = -999
      
      ! ���ֳֳʻ�
      dx = (xmax - xmin) / ( nx * xsub )
      
      ! �٤�礦 CPU �θĿ�������
      sx: do ix = 1, xsub
        call MPIWrapperCartShift( comm_cart, direction, ix, rank, tmp )
        if ( rank < 0 ) then 
          lp = ix - 1
          exit sx
        end if
      end do sx
      
      !��ǧ
      if ( lp == -999 ) then 
        call MessageNotify( "E", "axesset_init", "lp (x) is undifined" )
      end if
      
      ! ����κ���
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
      != y �����κ�ɸ�ͤȳʻ����ֳ֤����ꤹ��
      !
      
      ! ���ۤη�����ػ�
      implicit none
      
      integer, intent(in) :: comm_cart
      
      ! ����ѿ�
      integer             :: jy   ! �롼��ź��
      integer             :: lp
      integer             :: rank, tmp
      integer, parameter  :: direction = 0
      
      allocate(y_Y(jmin:jmax))
      allocate(q_Y(jmin:jmax))
      allocate(y_dy(jmin:jmax))
      allocate(q_dy(jmin:jmax))
      
      ! �����
      lp = -999
      
      ! 2 �����ξ��� ymax = dx �Ȥ���. 
      if (.NOT. FlagCalc3D ) then 
        ymax = dx
      end if
      
      ! ���ֳֳʻ�
      dy = (ymax - ymin) / ( ny * ysub )
      
      ! �٤�礦 CPU �θĿ�������
      sy: do jy = 1, ysub
        call MPIWrapperCartShift( comm_cart, direction, jy, tmp, rank )
        if ( rank < 0 ) then 
          lp = jy - 1
          exit sy
        end if
      end do sy
      
      !��ǧ
      if ( lp == -999 ) then 
        call MessageNotify( "E", "axesset_init", "lp (y) is undifined" )
      end if
      
      ! ����κ���
      do jy = jmin, jmax
        q_Y(jy) = YMin + dy * ny * lp  + dy * jy
        y_Y(jy) = YMin + dy * ny * lp  + dy * jy - dy * 0.5d0
        y_dy(jy) = dy
        q_dy(jy) = dy
      end do
      
    end subroutine y_axis_init


    subroutine z_axis_init
      !
      != z �����κ�ɸ�ͤȳʻ����ֳ֤����ꤹ��
      !
      
      ! ���ۤη�����ػ�
      implicit none
      
      ! ����ѿ�
      integer                 :: kz   ! �롼��ź��
      
      allocate(z_Z(kmin:kmax))
      allocate(r_Z(kmin:kmax))
      allocate(z_dz(kmin:kmax))
      allocate(r_dz(kmin:kmax))
      
      ! ���ֳֳʻ�
      dz = (zmax - zmin) / nz
      
      ! ����κ���
      do kz = kmin, kmax
        r_Z(kz) = ZMin + dz * kz
        z_Z(kz) = ZMin + dz * (kz - 0.5)
        r_dz(kz) = dz
        z_dz(kz) = dz
      end do
    
    end subroutine z_axis_init
  
  end subroutine axesset_init

end module axesset

