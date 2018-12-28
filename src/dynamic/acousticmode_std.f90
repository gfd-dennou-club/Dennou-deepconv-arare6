!= Module acousticmode_std
!
! Authors::   �����̰�ϯ(SUGIYAMA Ko-ichiro)
! Version::   $Id: acousticmode_std.f90,v 1.3 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module acousticmode_std
  !
  ! ���ȥ⡼�ɤ˴ؤ���׻��롼�����«�ͤ��⥸�塼��
  !
  !   ��ʿ explicit
  !   ��ľ implicit
  !
  ! ��ʬʿ�ѥ⥸�塼����Ѥ��Ʒ׻���Ԥ�

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

  character(*), parameter:: module_name = 'acousticmode_std'
                                                  ! �⥸�塼���̾��.
                                                  ! Module name
  ! public 
  !
  public acousticmode_std_init
  public acousticmode_std_exp
  public acousticmode_std_imp
  
contains

  subroutine acousticmode_std_exp(                   &
    & pyz_VelXN, xqz_VelYN, xyr_VelZN, xyz_ExnerN,   & !(IN)
    & xyz_VelDivN,                                   & !(OUT)
    & pyz_PGrad, pyz_SWF,                            & !(OUT)
    & xqz_PGrad, xqz_SWF                             & !(OUT)
    & )
    
    ! �⥸�塼��ƤӽФ�
    !
    use dc_types,  only : DP
    use gridset,   only : imin,            &! x ����������β���
      &                   imax,            &! x ����������ξ��
      &                   jmin,            &! y ����������β���
      &                   jmax,            &! y ����������ξ��
      &                   kmin,            &! z ����������β���
      &                   kmax,            &! z ����������ξ��
      &                   nx, ny, nz
    use constants, only : CpDry             ! ������ʬ����Ǯ
    use basicset,  only : pyz_VPTempBZ,    &! ���ܾ�β���
                          xqz_VPTempBZ      ! ���ܾ�β���
    use differentiate_center2,             &
      &            only : xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                   pyz_dx_xyz, xqz_dy_xyz

    ! ���ۤη�����ػ�
    !
    implicit none
    
    ! �ѿ������
    !
    real(DP), intent(in)  :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xqz_VelYN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_VelDivN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: xqz_PGrad(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: pyz_SWF(1:nx,1:ny,1:nz)
    real(DP), intent(out) :: xqz_SWF(1:nx,1:ny,1:nz)

    real(DP)              :: aaa_tmp(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------------
    ! ®�٤μ�«��׻�
    !
    xyz_VelDivN =                 &
      &   xyz_dx_pyz( pyz_VelXN ) &
      & + xyz_dy_xqz( xqz_VelYN ) &
      & + xyz_dz_xyr( xyr_VelZN )

    !------------------------------------------------------------------
    ! X ����

!    aaa_tmp   = + AlphaH * pyz_dx_xyz( xyz_VelDivN )     
    aaa_tmp   =   pyz_dx_xyz( AlphaH * xyz_VelDivN )     
    pyz_SWF   =   aaa_tmp(1:nx,1:ny,1:nz)
    
    aaa_tmp   = - CpDry * pyz_VPTempBZ * pyz_dx_xyz( xyz_ExnerN ) 
    pyz_PGrad =   aaa_tmp(1:nx,1:ny,1:nz)

    !------------------------------------------------------------------
    ! Y ����
    
    aaa_tmp   =   AlphaH * xqz_dy_xyz( xyz_VelDivN ) 
!    aaa_tmp   =   xqz_dy_xyz( AlphaH * xyz_VelDivN ) 
    xqz_SWF   =   aaa_tmp(1:nx,1:ny,1:nz)

    aaa_tmp   = - CpDry * xqz_VPTempBZ * xqz_dy_xyz( xyz_ExnerN ) 
    xqz_PGrad =   aaa_tmp(1:nx,1:ny,1:nz)

  end subroutine Acousticmode_std_exp

!!!------------------------------------------------------------!!!

  subroutine acousticmode_std_imp(                   &
    & pyz_VelXA, xqz_VelYA, xyr_VelZN, xyz_VelDivN,  & !(IN)
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
    use constants,only : CpDry             ! ������ʬ����Ǯ
!    use timeset,  only : DelTimeShort, TimeN
    use timeset,  only : DelTimeShort
    use basicset, only : xyz_VPTempBZ,    &!���ܾ�β�����
                         xyr_VPTempBZ      !���ܾ�β�����
    use axesset,  only : dz
    use differentiate_center2, &
      &           only : xyr_dz_xyz, xyz_dz_xyr, xyz_dx_pyz, xyz_dy_xqz
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 

    ! ���ۤη�����ػ�
    !
    implicit none

    ! �������ѿ�
    !
    real(DP), intent(in)   :: pyz_VelXA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !®�� u [��+����]
    real(DP), intent(in)   :: xqz_VelYA(imin:imax,jmin:jmax,kmin:kmax) 
                                                           !®�� v [��+����]
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
    
    ! ����ѿ����
    !
    real(DP)               :: D(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: E(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: F(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: F0(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: xyr_DExnerDz(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: xyr_DVelDivDz(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: aaa_tmp(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)               :: dt            ! û�����ֳʻҴֳ�
    integer                :: INFO          ! ��Υ���ǥ������
    integer                :: i, j, k
      
    real(DP)               :: TX(nz,nxny)    !������ž�֤������
    character(1),parameter :: TRANS = 'N'
    

    !---------------------------------------------------------------
    ! Initialize
    !
    xyz_ExnerA = 0.0d0
    dt = DelTimeShort

    !---------------------------------------------------------------
    !����׻��Τ���η���
    
    ! ���̤��Ƹ������ʬ����˷׻� 
    !
    xyr_DExnerDZ =  xyr_dz_xyz( xyz_ExnerN ) 
    
    xyr_DVelDivDZ = xyr_dz_xyz( xyz_VelDivN )
    
    E =                                      &
      & - ( 1.0d0 - beta ) * xyr_DExnerDZ    &
      & + (                                  &
      &     AlphaV * xyr_DVelDivDZ           &
      &   + xyr_DVelZDtNl                    &
      &   )                                  &
      &   / xyr_CpVPTempBZ
    
    F0  =                                                            &
      & + xyr_DensVPTempBZ                                           &
      &   * (                                                        &
      &     + xyr_VelZN                                              &  
      &     - xyr_CpVPTempBZ * ( 1.0d0 - beta) * xyr_DExnerDZ * dt   &
      &     + AlphaV * xyr_DVelDivDZ  * dt                           &
      &     + xyr_DVelZDtNl  * dt                                    &
      &    )
    
    F =                                            &
      & - beta * dt                                & 
      &   * xyz_VelSoundSQBZ                       &  
      &   / xyz_CpDensVPTempSQBZ                   &  
      &   * xyz_dz_xyr( F0 )                       &
      & + ( xyz_DExnerDtNl + xyz_DExnerDtNs ) * dt
    
    D =                                                           &
      & + xyz_ExnerN                                              &
      & - (1.0d0 - beta) * dt                                     &
      &   * xyz_VelSoundSQBZ                                      &  
      &   / xyz_CpDensVPTempSQBZ                                  &  
      &   * xyz_dz_xyr( xyr_DensVPTempBZ * xyr_VelZN )            &
      & - xyz_VelSoundSQBZ * dt                                   &
      &   / (CpDry * xyz_VPTempBZ )                               &
      &   * ( xyz_dx_pyz( pyz_VelXA ) + xyz_dy_xqz( xqz_VelYA ) ) &
      & + F
    
    D(:,:,1) =                                     &
      & + D(:,:,1)                                 &
      & - beta * (dt * dt)                         &
      &   * xyz_VelSoundSQBZ(:,:,1)                &  
      &   / xyz_CpDensVPTempSQBZ(:,:,1)            &  
      &   * xyr_CpDensVPTempSQBZ(:,:,0)            &
      &   * E(:,:,0)                               &
      &   / dz
    
    D(:,:,nz) =                                    &
      & + D(:,:,nz)                                &
      & + beta * (dt * dt)                         &
      &   * xyz_VelSoundSQBZ(:,:,nz)               &  
      &   / xyz_CpDensVPTempSQBZ(:,:,nz)           &  
      &   * xyr_CpDensVPTempSQBZ(:,:,nz)           &
      &   * E(:,:,nz)                              &
      &   / dz

!    call HistoryAutoPut(TimeN, 'D', D(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'E', E(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'F', F(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENs', xyz_DExnerDtNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENl', xyz_DExnerDtNl(1:nx,1:ny,1:nz))
    
    !-----------------------------------------------------------
    !ϢΩ�켡�������β�����
    
    ! LAPACK �λ��ͤ˹�碌���ѷ� 
    !
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          TX(k, i + nx * (j - 1)) = D(i,j,k)
        end do
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
      do j = 1, ny
        do i = 1, nx
          xyz_ExnerA(i,j,k) = TX(k, i + nx * (j - 1 ))
        end do
      end do
    end do
    

    !------------------------------------------------------------
    ! ��ľ®��
    !
!    aaa_tmp =  AlphaV * xyr_dz_xyz( xyz_VelDivN ) 
    aaa_tmp =  xyr_dz_xyz( AlphaV * xyz_VelDivN ) 
    xyr_SWF =  aaa_tmp(1:nx,1:ny,1:nz)
    
    aaa_tmp =                                               &
      & - CpDry * xyr_VPTempBZ                              &
      &   * (                                               &
      &         beta           * xyr_dz_xyz( xyz_ExnerA )   &
      &       + (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerN )   &
      &     )                                                
    xyr_PGrad =  aaa_tmp(1:nx,1:ny,1:nz)
    
  end subroutine Acousticmode_std_imp
  

!!!--------------------------------------------------------------------!!!
  subroutine acousticmode_std_init( AlphaSound )
    !
    !�������ʡ��ؿ��򱢲�ˡ�ǲ򤯺ݤ�ɬ�פȤʤ�, ������������Ǥ���, 
    !LU ʬ���Ԥ�. 
    !
    
    ! �⥸�塼���ɤ߹���
    !
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : FlagCalc3D,      &
      &                    imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x ������ʪ���ΰ�ξ��
      &                    ny,              &! x ������ʪ���ΰ�ξ��
      &                    nz                ! y ������ʪ���ΰ�ξ��
    use constants,  only : CpDry           ! ������ʬ����Ǯ
    use timeset,    only : DelTimeShort
    use axesset,    only : dx, dy, dz        ! �ʻҴֳ�
    use basicset,   only : xyz_VelSoundBZ,  &!���ܾ�β�® 
      &                    xyz_DensBZ,      &!���ܾ��̩��
      &                    xyz_VPTempBZ      !���ܾ�β���
    use average,    only : xyr_xyz
    
    !���ۤη�����ػ�
    !
    implicit none

    !�ѿ������
    !
    real(DP), intent(in) :: AlphaSound
    real(DP)             :: r_CpDensVPTempSQBZ(kmin:kmax)
    real(DP)             :: z_VelSoundSQBZ(kmin:kmax) 
    real(DP)             :: z_CpDensVPTempSQBZ(kmin:kmax) 
    real(DP)             :: dt      ! û�����ֳʻ�
    integer              :: INFO    !��Υ���ǥ����������å�

    !----------------------------------------------------------------
    ! �����

    ! ���ȸ����θ��그�������
    ! 
    ! ����ģͽ�������̺� 49 p53 �˽���, ��ʿ�ȱ�ľ�Ȥ�ʬ���ƹͤ���. 
    !
!    AlphaH = AlphaSound * ( Min( dx * dx, dy * dy ) ) / DelTimeShort
!!    AlphaH = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz ) ) / DelTimeShort
!    AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz ) ) / DelTimeShort
    if ( FlagCalc3D ) then 
      AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
    else
      AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
    end if


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
    z_VelSoundSQBZ        = xyz_VelSoundSQBZ(nx,ny,:)
    z_CpDensVPTempSQBZ(:) = xyz_CpDensVPTempSQBZ(nx,ny,:)
    r_CpDensVPTempSQBZ(:) = xyr_CpDensVPTempSQBZ(nx,ny,:)

    !
    !
    A(2:nz-1) =                               &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(2:nz-1)           &
      &    / z_CpDensVPTempSQBZ(2:nz-1)       &
      &    * (                                &
      &         r_CpDensVPTempSQBZ(2:nz-1)    &
      &       + r_CpDensVPTempSQBZ(1:nz-2)    &
      &       )                               &
      &    * ( dt * dt )                      &
      &    / ( dz * dz )

    A(1) =                                    &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(1)                &
      &    / z_CpDensVPTempSQBZ(1)            &
      &   * r_CpDensVPTempSQBZ(1)             &
      &   * ( dt * dt )                       &
      &   / ( dz * dz ) 

    A(nz) =                                   &
      & + 1.0d0                               &
      & + ( beta * beta )                     &
      &    * z_VelSoundSQBZ(nz)               &
      &    / z_CpDensVPTempSQBZ(nz)           &
      &   * r_CpDensVPTempSQBZ(nz-1)          &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )  

    B(2:nz) =                                 &
      & - ( beta * beta )                     &
      &    * z_VelSoundSQBZ(1:nz-1)           &
      &    / z_CpDensVPTempSQBZ(1:nz-1)       &
      &   * r_CpDensVPTempSQBZ(1:nz-1)        &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )
    
    C(1:nz-1) =                               &
      & - ( beta * beta )                     &
      &    * z_VelSoundSQBZ(2:nz)             &
      &    / z_CpDensVPTempSQBZ(2:nz)         &
      &   * r_CpDensVPTempSQBZ(1:nz-1)        &
      &   * ( dt * dt )                       &
      &   / ( dz * dz )

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

  end subroutine Acousticmode_std_init

end module Acousticmode_std
