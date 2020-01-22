! Sample program 
! 
!
program soundwave

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
    &                         nx, ny, nz, ncmax, FlagCalc3D
  use timeset,         only : timeset_init,                            &
    &                         TimesetDelTimeHalf, TimesetProgress,     &
    &                         TimeA, TimeN, TimeB,                     &
    &                         Nstep, NstepShort, NstepLong,            &
    &                         NstepOutput, FlagInitialRun,             &
    &                         DelTimeLong, DelTimeShort
  use axesset,         only : axesset_init, x_X, y_Y, z_Z, XMax, ZMax, &
    &                         dx, dy, dz
  use constants,       only : constants_init, CpDry
  use namelist_util,   only : NmlutilInit
  
  ! �������⥸�塼��
  !
  use setmargin,       only : SetMargin_init, & 
    &                         SetMargin_pyz, SetMargin_xqz, SetMargin_xyr, SetMargin_xyz
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
  real(DP), allocatable :: pyz_VelXAs(:,:,:)
  real(DP), allocatable :: pyz_VelXNs(:,:,:)
  real(DP), allocatable :: pyz_VelXAl(:,:,:)
  real(DP), allocatable :: pyz_VelXNl(:,:,:)
  real(DP), allocatable :: pyz_VelXBl(:,:,:)
  real(DP), allocatable :: xqz_VelYAs(:,:,:)
  real(DP), allocatable :: xqz_VelYNs(:,:,:)
  real(DP), allocatable :: xqz_VelYAl(:,:,:)
  real(DP), allocatable :: xqz_VelYNl(:,:,:)
  real(DP), allocatable :: xqz_VelYBl(:,:,:)
  real(DP), allocatable :: xyr_VelZAs(:,:,:)
  real(DP), allocatable :: xyr_VelZNs(:,:,:)
  real(DP), allocatable :: xyr_VelZAl(:,:,:)
  real(DP), allocatable :: xyr_VelZNl(:,:,:)
  real(DP), allocatable :: xyr_VelZBl(:,:,:)
  real(DP), allocatable :: xyz_ExnerAs(:,:,:)
  real(DP), allocatable :: xyz_ExnerNs(:,:,:)
  real(DP), allocatable :: xyz_ExnerAl(:,:,:)
  real(DP), allocatable :: xyz_ExnerNl(:,:,:)
  real(DP), allocatable :: xyz_ExnerBl(:,:,:)
  real(DP), allocatable :: xyz_PTempAs(:,:,:)
  real(DP), allocatable :: xyz_PTempNs(:,:,:)
  real(DP), allocatable :: xyz_PTempAl(:,:,:)
  real(DP), allocatable :: xyz_PTempNl(:,:,:)
  real(DP), allocatable :: xyz_PTempBl(:,:,:)

  real(DP), allocatable :: pyz_DVelXDtNs(:,:,:)
  real(DP), allocatable :: pyz_DVelXDtNl(:,:,:)
  real(DP), allocatable :: xqz_DVelYDtNs(:,:,:)
  real(DP), allocatable :: xqz_DVelYDtNl(:,:,:)
  real(DP), allocatable :: xyr_DVelZDtNs(:,:,:)
  real(DP), allocatable :: xyr_DVelZDtNl(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNs(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNl(:,:,:)
  real(DP), allocatable :: xyz_DPTempDtNs(:,:,:)
  real(DP), allocatable :: xyz_DPTempDtNl(:,:,:)

  real(DP), allocatable :: xyz_VelDivNs(:,:,:)

  real(DP), parameter :: AlphaSound = 5.0d-2
  real(DP)            :: AlphaH
  real(DP)            :: AlphaV  

  real(DP), parameter :: xyz_PTempBZ = 300.0d0      !����Ū�������
  real(DP), parameter :: xyz_DensBZ = 1.1627d0      !����Ū�������
  real(DP), parameter :: xyz_VelSoundBZ = 346.96d0  !����Ū�������

  real(DP), parameter :: xyr_PTempBZ = 300.0d0      !����Ū�������
  real(DP), parameter :: xyr_DensBZ = 1.1627d0      !����Ū�������  
  
  ! ʪ���ѥ�᥿
  ! ����ե����� (.conf) �Ǥϻ���Ǥ��ʤ����ܤ��ۤ˻���.  
  ! 
!  real(DP), parameter :: nu    = 1.0d-1  ! �Ȼ�����
  real(DP), parameter :: VelX0 = 1.0d0   ! ��ή®�� X ����
  real(DP), parameter :: VelZ0 = 0.0d0   ! ��ή®�� Z ����
  real(DP), parameter :: DelMax = 5.0d-2 ! ���ι⤵
  real(DP), parameter :: Xc = 5.0d3      ! �����濴����
  real(DP), parameter :: Zc = 5.0d3      ! �����濴����
  real(DP), parameter :: Xr = 4.0d2      ! ������
  real(DP), parameter :: Zr = 4.0d2      ! ������

  ! I/O
  !
  character(STRING) :: cfgfile ! NAMELIST �ե�����̾ ; NAMELIST fine name

  ! do �롼���ѿ� ; do loop variable
  !
  integer :: it = 0
  integer :: i, k, tau


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
        xyz_ExnerNl(i,:,k) =                                   & 
             &    DelMax                                       &
             &  * dexp(                                        &
             &       - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1  &
             &       - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1  &
             &    )        
    end do
  end do

  ! �������
  !
  call SetMargin_xyz( xyz_ExnerNl )

  ! �����
  !
  xyz_ExnerBl = xyz_ExnerNl

!  write(*,*) xyz_ExnerNl(1,1,:)
  

!!!
!!! ���ȸ����
!!!
  if ( FlagCalc3D ) then
     AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
     AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
  else
     AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
     AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
  end if

!  write(*,*) AlphaH
!  write(*,*) AlphaV
 
!!!
!!! ������ʬ
!!!
  do while ( Nstep <= NstepLong )

     ! tendency (����) �η׻�.
     ! �׻��ΰ������Τ���ˤϿ��ͳȻ��बɬ�פ�������̵��.
     ! Ĺ���ַ׻��Ǥ��ʤ������Ϥ褷�Ȥ���.
     ! �ǥե���ȤǤ� VelZ0 = 0 �Ȥ��Ƥ���.
     ! 
!     xyz_DZetaDtNl =                                     &
!          &  - VelX0 * xyz_pyz(pyz_dx_xyz(xyz_ZetaNl))   &
!          &  - VelZ0 * xyz_xyr(xyr_dz_xyz(xyz_ZetaNl))   
         
     ! ���̤���ʬ (leap-frog)
     !
!     xyz_ZetaAl = xyz_ZetaBl + 2.0d0 * DelTimeLong * xyz_DZetaDtNl
          
     ! û�����֥��ƥåפν���ͺ���.
     ! Initial values set up for time integration with short time step.
     !
     pyz_VelXNs  = pyz_VelXBl
     xqz_VelYNs  = xqz_VelYBl
     xyr_VelZNs  = xyr_VelZBl
     xyz_ExnerNs = xyz_ExnerBl
     xyz_PTempNs = xyz_PTempBl
     
     ! û�����֥��ƥåפλ�����ʬ. �����顼ˡ������.
     ! Time integration with short time step.
     !
     Euler: do tau = 1, NstepShort
        
        ! ®�٤�ȯ��
        ! divergence of velocity
        !
        xyz_VelDivNs =                       &
             &    xyz_dx_pyz( pyz_VelXNs )   &
             &  + xyz_dy_xqz( xqz_VelYNs )   &
             &  + xyz_dz_xyr( xyr_VelZNs )
       
        ! ���Ϲ�η׻�. �۲��� (HE-VE)
        ! pressure terms 
        !
        pyz_DVelXDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * pyz_dx_xyz( xyz_ExnerNs ) &
             &   + AlphaH * pyz_dx_xyz( xyz_VelDivNs )

        xqz_DVelYDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * xqz_dy_xyz( xyz_ExnerNs ) &
             &   + AlphaH * xqz_dy_xyz( xyz_VelDivNs )

        xyr_DVelZDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * xyr_dz_xyz( xyz_ExnerNs ) &
             &   + AlphaV * xyr_dz_xyz( xyz_VelDivNs )

        xyz_DExnerDtNs =                                                    &
             & - xyz_VelSoundBZ * xyz_VelSoundBZ / (CpDry * xyz_PTempBZ )   &
             &   * (  + xyz_dx_pyz( pyz_VelXNs )                            &
             &        + xyz_dy_xqz( xqz_VelYNs )                            &
             &        + xyz_dz_xyr( xyr_VelZNs )                            & 
             !             &        + xyz_dz_xyr( xyr_DensBZ * xyr_PTempBZ * xyr_VelZNs ) &
!             &           / (xyz_DensBZ * xyz_PTempBZ )                      &
             &     ) 
        
        ! ��ʬ
        ! time integration
        !
        pyz_VelXAs = pyz_VelXNs  + DelTimeShort * ( pyz_DVelXDtNl  + pyz_DVelXDtNs  )
        xqz_VelYAs = xqz_VelYNs  + DelTimeShort * ( xqz_DVelYDtNl  + xqz_DVelYDtNs  )
        xyr_VelZAs = xyr_VelZNs  + DelTimeShort * ( xyr_DVelZDtNl  + xyr_DVelZDtNs  )
        xyz_ExnerAs= xyz_ExnerNs + DelTimeShort * ( xyz_DExnerDtNl + xyz_DExnerDtNs )

        ! �������
        ! boundary condition
        !
        call SetMargin_pyz( pyz_VelXAs  )
        call SetMargin_xqz( xqz_VelYAs  )
        call SetMargin_xyr( xyr_VelZAs  )
        call SetMargin_xyz( xyz_ExnerAs )
        
        ! û�����֥��ƥåפΥ롼�פ�󤹤���ν���
        ! Renew prognostic variables for next short time step integration.
        !
        xyz_ExnerNs  = xyz_ExnerAs
        pyz_VelXNs   = pyz_VelXAs
        xqz_VelYNs   = xqz_VelYAs
        xyr_VelZNs   = xyr_VelZAs
        
     end do Euler
     
     ! �ǽ�Ū��û�����֥��ƥåפǤ��ͤ�Ĺ�����֥��ƥåפǤ��ͤȤߤʤ�
     ! Renew prognostic variables for next long time step integration.
     !
     xyz_ExnerAl  = xyz_ExnerAs
     pyz_VelXAl   = pyz_VelXAs
     xqz_VelYAl   = xqz_VelYAs
     xyr_VelZAl   = xyr_VelZAs
         
     ! Asselin �Υ�����ե��륿��������
     !
     call AsselinTimeFilter
     
     ! ������� ; Boundary condition
     !
     call SetMargin_xyz( xyz_PTempAl )

     ! �ҥ��ȥ�ե��������. 
     !
     if ( it == NstepOutput ) then
        call HistoryPut('VelX',  pyz_VelXNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('VelY',  xqz_VelYNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('VelZ',  xyr_VelZNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('Exner', xyz_ExnerNl(1:nx,1:ny,1:nz ))
        call HistoryPut('PTemp', xyz_PTempNl(1:nx,1:ny,1:nz ))
        it = 0
     end if
     
     ! �롼�פ�󤹤���ν���
     !
     pyz_VelXNl  = pyz_VelXAl
     pyz_VelXBl  = pyz_VelXNl
     xqz_VelYNl  = xqz_VelYAl
     xqz_VelYBl  = xqz_VelYNl
     xyr_VelZNl  = xyr_VelZAl
     xyr_VelZBl  = xyr_VelZNl
     xyz_ExnerNl = xyz_ExnerAl
     xyz_ExnerBl = xyz_ExnerNl
     xyz_PTempNl = xyz_PTempAl
     xyz_PTempBl = xyz_PTempNl
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

    ! ����ξ���ν����
    !
    call constants_init
    
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
    allocate( pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax) )
    
    allocate( pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DPTempDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax) )

    allocate( xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax) )
    
    !���������
    pyz_VelXAs = 0.0d0; pyz_VelXNs = 0.0d0; pyz_VelXAl = 0.0d0; pyz_VelXNl = 0.0d0; pyz_VelXBl = 0.0d0
    xqz_VelYAs = 0.0d0; xqz_VelYNs = 0.0d0; xqz_VelYAl = 0.0d0; xqz_VelYNl = 0.0d0; xqz_VelYBl = 0.0d0
    xyr_VelZAs = 0.0d0; xyr_VelZNs = 0.0d0; xyr_VelZAl = 0.0d0; xyr_VelZNl = 0.0d0; xyr_VelZBl = 0.0d0
    xyz_ExnerAs= 0.0d0; xyz_ExnerNs= 0.0d0; xyz_ExnerAl= 0.0d0; xyz_ExnerNl= 0.0d0; xyz_ExnerBl= 0.0d0
    xyz_PTempAs= 0.0d0; xyz_PTempNs= 0.0d0; xyz_PTempAl= 0.0d0; xyz_PTempNl= 0.0d0; xyz_PTempBl= 0.0d0
    
    pyz_DVelXDtNs = 0.0d0;    pyz_DVelXDtNl = 0.0d0
    xqz_DVelYDtNs = 0.0d0;    xqz_DVelYDtNl = 0.0d0
    xyr_DVelZDtNs = 0.0d0;    xyr_DVelZDtNl = 0.0d0
    xyz_DExnerDtNs = 0.0d0;   xyz_DExnerDtNl = 0.0d0
    xyz_DPTempDtNs = 0.0d0;   xyz_DPTempDtNl = 0.0d0

    xyz_VelDivNs = 0.0d0
     
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
    xyz_PTempNl = tfil * ( xyz_PTempBl + xyz_PTempAl ) + tfil2 * xyz_PTempNl

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
      file='soundwave.nc', title='2D diffusion model',                           &
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
         varname='VelX', dims=(/'x','y','z','t'/), & 
         longname='Velocity (x-direction)', units='m/s', xtype='float')       !���ϻ��� float ��.
    call HistoryAddVariable( &                   ! �ѿ����
         varname='VelY', dims=(/'x','y','z','t'/), & 
         longname='Velocity (y-direction)', units='m/s', xtype='float')       !���ϻ��� float ��.
    call HistoryAddVariable( &                   ! �ѿ����
         varname='VelZ', dims=(/'x','y','z','t'/), & 
         longname='Velocity (z-direction)', units='m/s', xtype='float')       !���ϻ��� float ��.
    call HistoryAddVariable( &                   ! �ѿ����
         varname='Exner', dims=(/'x','y','z','t'/), & 
         longname='Exner function', units='1', xtype='float')       !���ϻ��� float ��.
    call HistoryAddVariable( &                   ! �ѿ����
         varname='PTemp', dims=(/'x','y','z','t'/), & 
         longname='Potential temperature', units='K', xtype='float')       !���ϻ��� float ��.

  end subroutine fileset_init

end program soundwave
