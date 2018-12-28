!= Module acousticmode_2d
!
! Authors::   �����̰�ϯ(SUGIYAMA Ko-ichiro)
! Version::   $Id: acousticmode_2d.f90,v 1.3 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module acousticmode_2d
  !
  ! ���ȥ⡼�ɤ˴ؤ���׻��롼�����«�ͤ��⥸�塼��
  !
  !   ��ʿ explicit
  !   ��ľ implicit
  !

  !�⥸�塼���ɤ߹���
  !
  use dc_types,   only : DP

  !���ۤη�����ػ�
  !
  implicit none

  !°���λ���
  !
  private

  ! �ѿ������
  !
  real(DP), save :: beta  = 0.5d0         !����󥯥˥��륽��ˡ�ʤ� 0.5
                                          !��������ˡ�ʤ� 1

  real(DP), save :: AlphaH = 0.0d0        !���ȸ����θ��그�� (��ʿ����)
  real(DP), save :: AlphaV = 0.0d0        !���ȸ����θ��그�� (��ľ����)

  real(DP), allocatable, save :: A(:)     !����������г���ʬ
  real(DP), allocatable, save :: B(:)     !��������ξ廰����ʬ
  real(DP), allocatable, save :: C(:)     !��������β�������ʬ
  real(DP), allocatable, save :: AL1(:)   !LU ʬ��η�� L (1 ��������)
  integer,  allocatable, save :: IP(:)    !��ʬ�ԥܥåȸ򴹤ξ�����Ǽ

  real(DP), allocatable, save :: xyr_CpVPTempBZ(:,:,:)       !��������η׻����Ѥ�������
  real(DP), allocatable, save :: xyr_CpDensVPTempSQBZ(:,:,:) !��������η׻����Ѥ�������
  real(DP), allocatable, save :: xyr_DensVPTempBZ(:,:,:)     !��������η׻����Ѥ�������
  real(DP), allocatable, save :: xyz_VelSoundSQBZ(:,:,:)     !��������η׻����Ѥ�������
  real(DP), allocatable, save :: xyz_CpDensVPTempSQBZ(:,:,:) !��������η׻����Ѥ�������

  character(*), parameter:: module_name = 'acousticmode_2d'
                                                  ! �⥸�塼���̾��.
                                                  ! Module name
  ! public 
  !
  public acousticmode_2d_init
  public acousticmode_2d_exp
  public acousticmode_2d_imp

contains

  subroutine acousticmode_2d_exp(                    &
    & pyz_VelXN, xyr_VelZN,                          & !(IN)
    & xyz_ExnerN,                                    & !(IN)
    & xyz_VelDivN,                                   & !(OUT)
    & pyz_PGrad, pyz_SWF                             & !(OUT)
    & )
    
    ! �⥸�塼��ƤӽФ�
    !
    use gridset,   only : imin,            &! x ����������β���
      &                   imax,            &! x ����������ξ��
      &                   jmin,            &! y ����������β���
      &                   jmax,            &! y ����������ξ��
      &                   kmin,            &! z ����������β���
      &                   kmax,            &! z ����������ξ��
      &                   nx, ny, nz
    use axesset,   only : dx, dz            ! �ʻҴֳ�
    use constants, only : CpDry             ! ������ʬ����Ǯ
    use basicset,  only : pyz_VPTempBZ 

    ! ���ۤη�����ػ�
    !
    implicit none
    
    ! �ѿ������
    !
    real(DP), intent(in)  :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: pyz_SWF(1:nx,1:ny,1:nz)
    
    integer               :: i, k
    integer, parameter    :: j = 1
    
    !------------------------------------------------------------------
    ! ®�٤μ�«��׻�
    !
    do k = kmin + 1, kmax 
        do i = imin + 1, imax 

          xyz_VelDivN(i,j,k) =         &
            & + (                      &
            &     pyz_VelXN(i,   j, k) &
            &   - pyz_VelXN(i-1, j, k) &
            &   ) / dx                 &
            & + (                      &
            &     xyr_VelZN(i, j, k)   &
            &   - xyr_VelZN(i, j, k-1) &
            &   ) / dz
        
        end do
    end do

    ! �ͤ���ꤵ����
    !
    xyz_VelDivN(imin,:,:) = 1.0d10
    xyz_VelDivN(:,:,kmin) = 1.0d10

    !------------------------------------------------------------------
    ! X ����
    
    do k = 1, nz
        do i = 1, nx
        
          ! ���ȸ����
          !            
          pyz_SWF(i,j,k) =                 &
            &   AlphaH                     &
            &   * (                        &
            &       xyz_VelDivN(i+1, j, k) &
            &     - xyz_VelDivN(i,   j, k) &
            &     ) / dx
          
          ! ���Ϸ�����
          !
          pyz_PGrad(i,j,k) =                &
            & - CpDry * pyz_VPTempBZ(i,j,k) &
            &   * (                         &
            &       xyz_ExnerN(i+1, j, k)   &
            &     - xyz_ExnerN(i,   j, k)   &
            &     ) / dx
          
        end do
    end do
    
  end subroutine Acousticmode_2d_exp

!!!------------------------------------------------------------!!!

  subroutine acousticmode_2d_imp(                    &
    & pyz_VelXA, xyr_VelZN, xyz_VelDivN,             & !(IN)
    & xyz_ExnerN,                                    & !(IN)
    & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs, & !(IN)
    & xyz_ExnerA,                                    & !(OUT)
    & xyr_PGrad, xyr_SWF                             & !(OUT)
    & )
    !
    ! ����ˡ���Ѥ����������ʡ��ؿ�/��ľ®�٤η׻�. 
    !

    ! �⥸�塼����ɤ߹���
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x ����������β���
      &                  imax,            &! x ����������ξ��
      &                  jmin,            &! y ����������β���
      &                  jmax,            &! y ����������ξ��
      &                  kmin,            &! z ����������β���
      &                  kmax,            &! z ����������ξ��
      &                  nx, ny, nz,      &! ʪ���ΰ���礭��
      &                  nxny              ! ʪ���ΰ���礭�� (nx * ny)
    use axesset,  only : dx, dz            ! �ʻҴֳ�
    use constants,only : CpDry             ! ������ʬ����Ǯ
    use timeset,  only : DelTimeShort
    use basicset, only : xyz_VPTempBZ,    &!���ܾ�β�����
                            xyr_VPTempBZ      !���ܾ�β�����

    !���ۤη�����ػ�
    !
    implicit none
    
    !�������ѿ�
    real(DP), intent(in)   :: pyz_VelXA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !®�� u [��+����]
    real(DP), intent(in)   :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !®�� w [��]
    real(DP), intent(in)   :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
                                                           !\Div \Dvect{u}
    real(DP), intent(in)   :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !̵��������
    real(DP), intent(in)   :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z �����γ��Ϲ�[t]
    real(DP), intent(in)   :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !�������ʡ��ؿ��γ��Ϲ�[t]
    real(DP), intent(in)   :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !�������ʡ��ؿ��γ��Ϲ�[t]
    real(DP), intent(out)  :: xyz_ExnerA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !̵��������[��+����]
    real(DP), intent(out)  :: xyr_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out)  :: xyr_SWF(1:nx,1:ny,1:nz)
    
    !����ѿ����
    real(DP)               :: D(1:nx,1:ny,1:nz)  
    real(DP)               :: E(1:nx,1:ny,0:nz)
    real(DP)               :: F(1:nx,1:ny,1:nz)
    real(DP)               :: F0(1:nx,1:ny,0:nz)  
    real(DP)               :: xyr_DExnerDz(1:nx,1:ny,0:nz)  
    real(DP)               :: xyr_DVelDivDz(1:nx,1:ny,0:nz)  
    real(DP)               :: dt            ! û�����ֳʻҴֳ�
    integer                :: INFO          ! ��Υ���ǥ������
    integer                :: i, k
    integer, parameter     :: j = 1
      
    real(DP)               :: TX(nz,nxny)    !������ž�֤������
    character(1),parameter :: TRANS = 'N'
    
    
    !---------------------------------------------------------------
    ! Initialize
    !
    dt = DelTimeShort
    
    !---------------------------------------------------------------
    !����׻��Τ���η���
    
    ! ���̤��Ƹ������ʬ����˷׻� 
    !
    do k = 0, nz
        do i = 1, nx
        
          xyr_DExnerDz(i,j,k) =              &
            &  (                             &
            &   + xyz_ExnerN(i, j, k+1)      &
            &   - xyz_ExnerN(i, j, k)        &
            &   ) / dz                        
          
          xyr_DVelDivDz(i,j,k) =             &
            &  (                             &
            &  + xyz_VelDivN(i, j, k+1)      &
            &  - xyz_VelDivN(i, j, k)        &
            &  ) / dz  

      end do
    end do
    
    !  ź�����ϰϤ�, 1:nx, 1:ny, 0:nz
    !
    do k = 0, nz
        do i = 1, nx
        
          E(i,j,k) =                                    &
            & - ( 1.0d0 - beta ) * xyr_DExnerDZ(i,j,k)  &
            & + (                                       &
            &     AlphaV * xyr_DVelDivDz(i,j,k)         &
            &   + xyr_DVelZDtNl(i,j,k)                  &
            &   )                                       &
            &   / xyr_CpVPTempBZ(i,j,k)  
        
        end do
    end do
    
    ! ����ʬ�ؿ�
    !
    do k = 0, nz
        do i = 1, nx
        
          F0(i,j,k)  =                                  &
            & + xyr_DensVPTempBZ(i,j,k)                 &
            &   * (                                     &
            &     + xyr_VelZN(i,j,k)                    &
            &     - xyr_CpVPTempBZ(i,j,k)               &
            &       * ( 1.0d0 - beta )                  &
            &       * xyr_DExnerDZ(i,j,k) * dt          &
            &     + AlphaV * xyr_DVelDivDz(i,j,k) * dt  &
            &     + xyr_DVelZDtNl(i,j,k)  * dt          &
            &    )
        
        end do
    end do
    
    !����׻��Τ���η���
    !  ź�����ϰϤ�, 1:nx, 1:ny, 1:nz
    !
    do k = 1, nz
        do i = 1, nx
        
          F(i,j,k) =                                  &
            & - beta * dt                             &
            &   * xyz_VelSoundSQBZ(i,j,k)             &
            &   / xyz_CpDensVPTempSQBZ(i,j,k)         &
            &   * (                                   &
            &       F0(i, j, k)                       &
            &     - F0(i, j, k-1)                     &
            &     ) / dz                              &
            & + xyz_DExnerDtNl(i,j,k) * dt            &
            & + xyz_DExnerDtNs(i,j,k) * dt
        
        end do
    end do
    
    !����׻��Τ���η���
    !  ź�����ϰϤ�, 1:nx, 1:ny, 1:nz
    !
    do k = 1, nz 
        do i = 1, nx
        
          D(i,j,k) =                                               &
            & + xyz_ExnerN(i,j,k)                                  &
            & - (1.0d0 - beta) * dt                                &
            &   * xyz_VelSoundSQBZ(i,j,k)                          &
            &   / xyz_CpDensVPTempSQBZ(i,j,k)                      &
            &   * (                                                &
            &       xyr_DensVPTempBZ(i,j,k)   * xyr_VelZN(i,j,k)   &
            &     - xyr_DensVPTempBZ(i,j,k-1) * xyr_VelZN(i,j,k-1) &
            &     ) / dz                                           &
            & - xyz_VelSoundSQBZ(i,j,k) * dt                       &
            &   / (CpDry * xyz_VPTempBZ(i,j,k))                    &
            &   * (                                                &
            &     + (                                              &
            &         pyz_VelXA(i,   j, k)                         &
            &       - pyz_VelXA(i-1, j, k)                         &
            &       ) / dx                                         &
            &     )                                                &
            & + F(i,j,k)
        
        end do
    end do
    
    ! ����׻��Τ���η���
    !
    do i = 1, nx
        
        D(i,j,1) =                                     &
          & + D(i,j,1)                                 &
          & - beta * (dt * dt)                         &
          &   * xyz_VelSoundSQBZ(i,j,1)                &
          &   / xyz_CpDensVPTempSQBZ(i,j,1)            &
          &   * xyr_CpDensVPTempSQBZ(i,j,0)            &
          &   * E(i,j,0)                               &
          &   / dz
      
        D(i,j,nz) =                                    &
          & + D(i,j,nz)                                &
          & + beta * (dt * dt)                         &
          &   * xyz_VelSoundSQBZ(i,j,nz)               &
          &   / xyz_CpDensVPTempSQBZ(i,j,nz)           &
          &   * xyr_CpDensVPTempSQBZ(i,j,nz)           &
          &   * E(i,j,nz)                              &
          &   / dz
      
      end do
    
    !-----------------------------------------------------------
    !ϢΩ�켡�������β�����
    
    ! LAPACK �λ��ͤ˹�碌���ѷ� 
    !
    do k = 1, nz
        do i = 1, nx
          TX(k,i) = D(i,j,k)
        end do
    end do
    
    !�����η׻�. LAPACK �����. 
    !
    call DGTTRS(TRANS, nz, nxny, C, A, B, AL1, IP, TX, nz, INFO)
    
!    !��Υ���ǥ�����������å�. 
!    !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

    do k = 1, nz
        do i = 1, nx
          xyz_ExnerA(i,j,k) = TX(k,i)
        end do
    end do
    
    !------------------------------------------------------------
    ! ��ľ®��
    !
    ! ��ľ®�٤� k = nz ���ͤ϶������ˤ�äƷ�ޤ�Τ�, k ������
    ! �롼�פ� 1 ~ nz-1 �Ȥ���.
    ! 
    
    do k = 1, nz-1
        do i = 1, nx
          
          ! ���ȸ����
          !            
          xyr_SWF(i,j,k) =                 &
            & + AlphaV                     &
            &   * (                        &
            &       xyz_VelDivN(i,j,k+1)   &
            &     - xyz_VelDivN(i,j,k)     &
            &     ) / dz
          
          ! ���Ϸ�����
          !
          xyr_PGrad(i,j,k) =                &
            & - CpDry * xyr_VPTempBZ(i,j,k) &
            &   * (                         &
            &       beta                    &
            &       * (                     &
            &           xyz_ExnerA(i,j,k+1) &
            &         - xyz_ExnerA(i,j,k)   &
            &         )                     &
            &     + (1.0d0 - beta)          &
            &       * (                     &
            &           xyz_ExnerN(i,j,k+1) &
            &         - xyz_ExnerN(i,j,k)   &
            &         )                     &
            &     ) / dz

        end do
    end do

    ! �����
    !   ���֥롼���������ͤȤ��뤿��ˤ�, ������ͤ���ꤵ����ɬ�פ����뤿��. 
    ! 
    xyz_ExnerA(imin:0,:,:) = 0.0d0
    xyz_ExnerA(nx+1:imax,:,:) = 0.0d0
    xyz_ExnerA(:,:,kmin:0) = 0.0d0
    xyz_ExnerA(:,:,nz+1:kmax) = 0.0d0
    xyr_SWF(:,:,nz)   = 0.0d0
    xyr_PGrad(:,:,nz) = 0.0d0

  end subroutine Acousticmode_2d_imp
  
!!!--------------------------------------------------------------------!!!

  subroutine acousticmode_2d_init( AlphaSound )
    !
    ! ���ȥ⡼�ɤη׻��⥸�塼��ν����
    ! �������ʡ��ؿ��򱢲�ˡ�ǲ򤯺ݤ�ɬ�פȤʤ�, ������������Ǥ���, LU ʬ���Ԥ�. 
    !
    
    ! �⥸�塼���ɤ߹���
    !
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x ������ʪ���ΰ�ξ��
      &                    ny,              &! x ������ʪ���ΰ�ξ��
      &                    nz                ! y ������ʪ���ΰ�ξ��
    use constants,  only : CpDry             ! ������ʬ����Ǯ
    use timeset,    only : DelTimeShort
    use axesset,    only : dx, dz            ! �ʻҴֳ�
    use basicset,   only : xyz_VelSoundBZ,  &!���ܾ�β�® 
      &                    xyz_DensBZ,      &!���ܾ��̩��
      &                    xyz_VPTempBZ      !���ܾ�β���
    use average,    only : xyr_xyz 

    !���ۤη�����ػ�
    implicit none

    real(DP), intent(in) :: AlphaSound
    real(DP)             :: r_CpDensVPTempSQBZ(0:nz)
    real(DP)             :: z_VelSoundSQBZ(1:nz) 
    real(DP)             :: z_CpDensVPTempSQBZ(1:nz) 
    real(DP)             :: dt      ! û�����ֳʻ�
    integer              :: INFO    !��Υ���ǥ����������å�
    integer              :: k

    !----------------------------------------------------------------
    ! �����

    ! ���ȸ����θ��그�������
    ! 
    ! ����ģͽ�������̺� 49 p53 �˽���, ��ʿ�ȱ�ľ�Ȥ�ʬ���ƹͤ���. 
    !
!    AlphaH = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort
!    AlphaV = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort
    AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
    AlphaV = AlphaSound * ( Min( dx * dx, dz * dz ) ) / DelTimeShort

    !-------------------------------------------------------------------
    ! ����
    !
    call MessageNotify( "M", module_name, "AlphaH = %f", d=(/AlphaH/) )
    call MessageNotify( "M", module_name, "AlphaV = %f", d=(/AlphaV/) )

    ! �ѿ�̾��Ĺ�������Τ�, ̾�����֤�������
    !
    dt = DelTimeShort

    ! ����γ���դ�
    !
    allocate( A(1:nz) )
    allocate( B(2:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_VelSoundSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CpDensVPTempSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_CpDensVPTempSQBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DensVPTempBZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_CpVPTempBZ(imin:imax,jmin:jmax,kmin:kmax) )
   
    !----------------------------------------------------------------
    ! �������󤪤�Ӷ��̤������Ѥ����������ͤ����
    !
    xyz_VelSoundSQBZ     = xyz_VelSoundBZ * xyz_VelSoundBZ
    xyz_CpDensVPTempSQBZ = CpDry * xyz_DensBZ * xyz_VPTempBZ * xyz_VPTempBZ
    xyr_CpDensVPTempSQBZ = xyr_xyz( xyz_CpDensVPTempSQBZ )
    xyr_DensVPTempBZ     = xyr_xyz( xyz_DensBZ * xyz_VPTempBZ )
    xyr_CpVPTempBZ       = xyr_xyz( CpDry * xyz_VPTempBZ )

    ! ��ľ�����Τߤ��ѿ��ˤĤ��Ƥ�, ���ܾ��Ȥ��Τ�, 
    ! nx, ny ���ͤ���ɽ�����뤳�ȤȤ���. 
    !
    z_VelSoundSQBZ(1:nz)     = xyz_VelSoundSQBZ(nx,ny,1:nz)
    z_CpDensVPTempSQBZ(1:nz) = xyz_CpDensVPTempSQBZ(nx,ny,1:nz)
    r_CpDensVPTempSQBZ(0:nz) = xyr_CpDensVPTempSQBZ(nx,ny,0:nz)

    do k = 2, nz-1
      A(k) =                                &
        & + 1.0d0                           &
        & + ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k)            &
        &    / z_CpDensVPTempSQBZ(k)        &
        &    * (                            &
        &         r_CpDensVPTempSQBZ(k)     &
        &       + r_CpDensVPTempSQBZ(k-1)   &
        &       )                           &
        &    * ( dt * dt )                  &
        &    / ( dz * dz )
    end do

    A(1) =                                  &
      & + 1.0d0                             &
      & + ( beta * beta )                   &
      &   * z_VelSoundSQBZ(1)               &
      &   / z_CpDensVPTempSQBZ(1)           &
      &   * r_CpDensVPTempSQBZ(1)           &
      &   * ( dt * dt )                     &
      &   / ( dz * dz ) 

    A(nz) =                                 &
      & + 1.0d0                             &
      & + ( beta * beta )                   &
      &   * z_VelSoundSQBZ(nz)              &
      &   / z_CpDensVPTempSQBZ(nz)          &
      &   * r_CpDensVPTempSQBZ(nz-1)        &
      &   * ( dt * dt )                     &
      &   / ( dz * dz )  

    do k = 2, nz
      B(k) =                                &
        & - ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k-1)          &
        &    / z_CpDensVPTempSQBZ(k-1)      &
        &   * r_CpDensVPTempSQBZ(k-1)       &
        &   * ( dt * dt )                   &
        &   / ( dz * dz )
    end do
    
    do k = 1, nz-1
      C(k) =                                &
        & - ( beta * beta )                 &
        &    * z_VelSoundSQBZ(k+1)          &
        &    / z_CpDensVPTempSQBZ(k+1)      &
        &   * r_CpDensVPTempSQBZ(k)         &
        &   * (dt * dt )                    &
        &   / ( dz * dz )
    end do

    !----------------------------------------------------------------
    ! ��������� LU ʬ��
    !

    ! ����γ������
    !
    allocate( AL1(1:nz-2), IP(1:nz) )

    ! �����η׻�. LAPACK �����. 
    !
    call DGTTRF(nz, C, A, B, AL1, IP, INFO)
    
    ! ��Υ���ǥ�����������å�. 
    !
    if (INFO /= 0) then
      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
      stop
    end if

  end subroutine Acousticmode_2d_init
  
end module Acousticmode_2d
