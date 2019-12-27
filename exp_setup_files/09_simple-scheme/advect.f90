! Sample program 
! 
!
program advect

!!!
!!! �⥸�塼�������. 
!!!

  ! gtool5 �⥸�塼�� (I/O) �����饤�֥��.
  ! gtool5 : http://www.gfd-dennou.org/arch/gtool/
  !
  use dc_types,          only: STRING, DP   ! DP �� double precision �ΰ�̣
  use gtool_history,     only: HistoryCreate, HistoryAddVariable, HistoryPut, HistoryClose

  ! �������⥸�塼�� (src/setup/ �ʲ��Υե�����)
  !
  use mpi_wrapper,     only : MPIWrapperInit,                          &
    &                         MPIWrapperFinalize
  use argset,          only : argset_init
  use gridset,         only : gridset_init,                            &
    &                         imin, imax, jmin, jmax, kmin, kmax,      &
    &                         nx, ny, nz, ncmax
  use timeset,         only : timeset_init,                            &
    &                         TimesetDelTimeHalf, TimesetProgress,     &
    &                         TimeA, TimeN, TimeB,                     &
    &                         Nstep, NstepShort, NstepLong,            &
    &                         NstepOutput, FlagInitialRun,             &
    &                         DelTimeLong
  use axesset,         only : axesset_init, x_X, y_Y, z_Z, XMax, ZMax
  use namelist_util,   only : NmlutilInit
  
  ! �������⥸�塼��
  !
  use setmargin,       only : SetMargin_init
  use setmargin,       only : SetMargin_xyz
  use differentiate_center4
  use average
  
!!!
!!! ���ۤη�����ػ�
!!! 
  implicit none

!!!
!!! �ѿ���� 
!!!  
  ! ͽ���ѿ�
  !
  real(DP), allocatable :: xyz_ZetaAl(:,:,:)    ! ���λ������
  real(DP), allocatable :: xyz_ZetaNl(:,:,:)    ! ���ߤ���
  real(DP), allocatable :: xyz_ZetaBl(:,:,:)    ! ���λ������
  real(DP), allocatable :: xyz_DZetaDtNl(:,:,:) ! tendency
  
  ! ʪ���ѥ�᥿
  ! ����ե����� (.conf) �Ǥϻ���Ǥ��ʤ����ܤ��ۤ˻���.  
  !
!  real(DP), parameter :: nu    = 1.0d-1  ! �Ȼ�����
  real(DP), parameter :: VelX0 = 1.0d0   ! ��ή®�� X ����
  real(DP), parameter :: VelZ0 = 0.0d0   ! ��ή®�� Z ����
  real(DP), parameter :: sigma = 0.1     ! ���ʬ�ۤ��礭��
  real(DP), parameter :: DelMax = 1.0d0  ! ���ι⤵
  real(DP), parameter :: Xc = 5.0d2      ! �����濴���� (���) (0 < Xc < 1)
  real(DP), parameter :: Zc = 5.0d2      ! �����濴���� (���) (0 < Zc < 1)
  real(DP), parameter :: Xr = 3.0d2      ! ������
  real(DP), parameter :: Zr = 3.0d2      ! ������

  ! I/O
  !
  character(STRING) :: cfgfile ! NAMELIST �ե�����̾ ; NAMELIST fine name

  ! do �롼���ѿ� ; do loop variable
  !
  integer :: it = 0
  integer :: i, k


!!!
!!! ����� (main/arare.f90 �� MainInit ���֥롼����δ�ñ��)
!!! 
  call MainInit
  
!!!
!!! ���������.
!!!

  ! �濴�˻���Ϳ����
  !
  do k = kmin, kmax
     do i = imin, imax
        xyz_ZetaNl(i,:,k) =                                    & 
             &    DelMax                                       &
             &  * dexp(                                        &
             &       - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1  &
             &       - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1  &
             &    )        
    end do
  end do

!  write(*,*) xyz_ZetaNl(95:105,1,100)

  ! �������
  !
  call SetMargin_xyz( xyz_ZetaNl )
  xyz_ZetaBl = xyz_ZetaNl

  
  !����ͤν���
  !
  call HistoryPut('zeta', xyz_ZetaNl(1:nx,1:ny,1:nz))

!!!
!!! ������ʬ
!!!
  do while ( Nstep <= NstepLong )

     ! tendency (����) �η׻�.
     ! �׻��ΰ������Τ���ˤϿ��ͳȻ��बɬ�פ�������̵��.
     ! Ĺ���ַ׻��Ǥ��ʤ������Ϥ褷�Ȥ���.
     ! �ǥե���ȤǤ� VelZ0 = 0 �Ȥ��Ƥ���.
     ! 
     xyz_DZetaDtNl =                                     &
          &  - VelX0 * xyz_pyz(pyz_dx_xyz(xyz_ZetaNl))   &
          &  - VelZ0 * xyz_xyr(xyr_dz_xyz(xyz_ZetaNl))   


     ! ���̤���ʬ (leap-frog)
     !
     xyz_ZetaAl = xyz_ZetaBl + 2.0d0 * DelTimeLong * xyz_DZetaDtNl

     write(*,*) xyz_ZetaAl(95:105,1,100)
     
     ! Asselin �Υ�����ե��륿��������
     !
     call AsselinTimeFilter
     
     ! ������� ; Boundary condition
     !
     call SetMargin_xyz( xyz_ZetaAl )  
          
     ! �ҥ��ȥ�ե��������. 
     !
     if ( it == NstepOutput ) then
        call HistoryPut('zeta', xyz_ZetaNl(1:nx,1:ny,1:nz))
        it = 0
     end if
     
     ! �롼�פ�󤹤���ν���
     !
     xyz_ZetaNl = xyz_ZetaAl
     xyz_ZetaBl = xyz_ZetaNl
     it = it + 1
     
     ! ����οʹ�
     !
     call TimesetProgress

  end do

!!!
!!! ��λ����
!!!

  ! ���ϥե�����Υ�����
  !
  call HistoryClose

  ! MPI END
  !
  call MPIWrapperFinalize
  
  stop

contains

!!!
!!! ������ץ��� (src/main/arare.f90 ���������֥롼����δ�ά��)  
!!!
  subroutine MainInit
  
    ! MPI
    !
    call MPIWrapperInit

    ! NAMELIST �ե�����̾���ɤ߹���
    !
    call argset_init( cfgfile ) !(out)
    
    ! NAMELIST �ե�����̾�Υ⥸�塼��ؤ���Ͽ
    !
    call NmlutilInit( cfgfile ) !(in)
    
    ! ������ʬ�ν����
    !
    call timeset_init
    
    ! �ʻ�������ν����
    !
    call gridset_init
    
    ! ���ξ���ν����                    
    !! ��) axesset �� gridset ��ޡ������Ƥ⹽��ʤ�.
    !
    call axesset_init
    
    ! I/O �ե�����̾�ν���� 
    !! ��) �������֥롼����.
    !
    call fileset_init
    
    ! �ޡ����������ν����
    !
    call SetMargin_Init
    
    ! �����ѿ��ν����
    !
    call VariableAllocate
    
    ! t=0 ���������ʬ��Ԥ�����, �ǽ�ΰ���ܤϥ����顼ˡ������. 
    !
    if ( FlagInitialRun ) call TimesetDelTimeHalf

  end subroutine MainInit
  
!!!
!!! ����� allocate
!!!
  subroutine VariableAllocate
    !
    !������Ȥ���, ����������, �ͤȤ��ƥ������������.
    !
    
    !���ۤη�����ػ�
    implicit none

    !����������
    allocate( xyz_ZetaBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ZetaNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ZetaAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DZetaDtNl(imin:imax,jmin:jmax,kmin:kmax) )

    !���������
    xyz_ZetaBl    = 0.0d0
    xyz_ZetaNl    = 0.0d0
    xyz_ZetaAl    = 0.0d0
    xyz_DZetaDtNl = 0.0d0
    
  end subroutine VariableAllocate

!!!
!!! ���֥ե��륿��. �꡼�ץե�å��ǲ򤯻���ɬ��. 
!!!
  
  subroutine AsselinTimeFilter
    !
    ! ���֥ե��륿��; Asselin �Υ�����ե��륿��������    
    !   t = 0.0 �ξ��ˤ� tfil = 0.0d0, ����ʳ��� tfil = 1.0d-1
    !   (t = 0 �λ��ϥ����顼ˡ��, ����ʳ��ϥ꡼�ץե�å�ˡ����ʬ���뤿��)
    !
    use TimeSet, only: tfil

    ! ���ۤη�����ػ�    
    !
    implicit none

    ! ����ѿ�
    !
    real(DP) :: tfil2

    tfil2 = 1.0d0 - 2.0d0 * tfil 
    xyz_ZetaNl = tfil * ( xyz_ZetaBl + xyz_ZetaAl ) + tfil2 * xyz_ZetaNl

  end subroutine AsselinTimeFilter
  

!!!
!!! src/io �ʲ��ǹԤäƤ��뤳�Ȥδ�ά��. gtool5 �����Ѥ������. 
!!! deepconv �ν���͡��ꥹ�����ȥե������ HistoryPut,
!!! deepconv �η�̤� HistoryAutoPut ��Ȥäƽ��Ϥ��Ƥ���.
!!!
!!! gtool5 : http://www.gfd-dennou.org/arch/gtool/
!!!  
  
  ! I/O �ν����. 
  ! HistoryPut, HistoryAddVarible ������. 
  !
  subroutine fileset_init
    
    call HistoryCreate( &                                  ! �ҥ��ȥ꡼����
      file='advect.nc', title='2D diffusion model',                             &
      source='Sample program of deepconv/arare6',                                &
      institution='GFD_Dennou Club deepconv project',                            &
      dims=(/'x','y','z','t'/), dimsizes=(/nx,ny,nz,0/),                         &
      longnames=(/'X-coordinate','Y-coordinate','Z-coordinate','time        '/), &
      units=(/'m','m','m','1'/),                                                 &
      origin=0.0, interval=real(NstepOutput * DelTimeLong) )
    
    call HistoryPut('x',x_X(1:nx))               ! �ѿ�����
    call HistoryPut('y',y_Y(1:ny))               ! �ѿ�����
    call HistoryPut('z',z_Z(1:nz))               ! �ѿ�����

    call HistoryAddVariable( &                   ! �ѿ����
         varname='zeta', dims=(/'x','y','z','t'/), & 
         longname='variable', units='1', xtype='float')       !���ϻ��� float ��.

  end subroutine fileset_init

end program advect
