!= ��ñ��ή����: ��ή�Ȼ���������
!
! Baker et al. (1998) ���Ѥ���줿�Ȼ������������ήǴ�����Ѥ���.
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA, Masatsugu
! Version::   $Id: turbulence_constkm.f90,v 1.3 2014/06/07 17:34:27 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Turbulence_constKm
  !
  ! ��ñ��ή����: ��ή�Ȼ���������
  !

  !�⥸�塼���ɤ߹��� 
  !
  use dc_types, only: DP, STRING

  !���ۤη�����ػ�
  !
  implicit none

  !°���λ���
  !
  private
  
  !�ؿ��� public ������
  !
  public Turbulence_constKm_Init
  public Turbulence_constKm_Forcing

  !�ѿ����
  !
  real(DP), save :: Cm     = 2.0d-1          !��ή���ͥ륮�����Ǽ��η��� 
  real(DP), save :: MixLen = 0.0d0           !ʿ�Ѻ����Υ
  real(DP), save :: ConstKm = 0.0d0          !��ư�̤��Ф�����ή�Ȼ�����   
  real(DP), save :: ConstKh = 0.0d0          !Ǯ���Ф�����ή�Ȼ�����
  real(DP), save :: FactorDispHeat = 0.0d0   !�����Ǯ���θ���뤫�Υ����å�  
  real(DP), save :: FactorDExnerDtTurb =1.0d0!��������������ή�Ȼ�����θ���뤫�Υ����å�

  character(*), parameter:: module_name = 'turbulence_constKm'

contains

!!!------------------------------------------------------------------------!!!
  subroutine turbulence_constKm_init
    !
    ! Turbulence �⥸�塼��ν�����롼����
    ! 

    !�⥸�塼���ɤ߹��� 
    !
    use dc_types,      only : STRING
    use axesset,       only : dx,            &! x �����γʻ����ֳ�
      &                       dy,            &! y �����γʻ����ֳ�
      &                       dz              ! z �����γʻ����ֳ�
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use gridset,       only : FlagCalc3D 
    use namelist_util, only : namelist_filename

    ! ���ۤη�����ػ�
    !
    implicit none

    ! ����ѿ�
    !
    integer :: unit

    !-------------------------------------------------------------------
    ! NAMELIST �����������
    NAMELIST /turbulence_constKm_nml/ &
      & Cm, ConstKm, ConstKh, FactorDExnerDtTurb, FactorDispHeat

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=turbulence_constKm_nml)
    close(unit)

    !-------------------------------------------------------------------
    ! �����Υ
    ! 2 �����׻��ξ��ˤ� DelY �˰�¸���ʤ��褦�ˤ��뤿��� if ʸ������.
    ! 
    if ( FlagCalc3D ) then 
      MixLen = (dx * dy * dz ) ** (1.0d0 / 3.0d0)
    else
      MixLen = sqrt( dx * dz ) 
    end if

    !-------------------------------------------------------------------
    ! tendency �ν���
    !
    call turbulence_constkm_output
   
    !-------------------------------------------------------------------
    ! Output
    !
    call MessageNotify( "M", module_name, "Cm = %f", d=(/Cm/))
    call MessageNotify( "M", module_name, "MixLen= %f", d=(/MixLen/))
    call MessageNotify( "M", module_name, "ConstKm= %f", d=(/ConstKm/))
    call MessageNotify( "M", module_name, "ConstKh= %f", d=(/ConstKh/))
    call MessageNotify( "M", module_name, &
      &                 "FactorDispHeat= %f", d=(/ FactorDispHeat /))
    call MessageNotify( "M", module_name, &
      &                 "FactorDExnerDtTurb= %f", d=(/ FactorDExnerDtTurb /))

  end subroutine turbulence_constKm_init

!!!------------------------------------------------------------------------!!!

  subroutine turbulence_constKm_forcing(       &
    & pyz_VelXBl, xqz_VelYBl,   xyr_VelZBl,    &
    & xyz_PTempBl,                             &
    & pyz_DVelXDt, xqz_DVelYDt,  xyr_DVelZDt,  &
    & xyz_DPTempDt,xyz_DExnerDt,               &
    & xyz_KmAl, xyz_KhAl &
    )
    
    !�⥸�塼����ɤ߹���
    !
    use gtool_historyauto, &
      &            only : HistoryAutoPut
    use dc_types,  only : DP
    use constants, only : CpDry
    use timeset,   only : TimeN
    use gridset,   only : imin,             &! x ����������β���
      &                   imax,             &! x ����������ξ��
      &                   jmin,             &! y ����������β���
      &                   jmax,             &! y ����������ξ��
      &                   kmin,             &! z ����������β���
      &                   kmax,             &! z ����������ξ��
      &                   nx,ny,nz       
    use basicset,  only : xyz_ExnerBZ,      &!���ܾ�Υ������ʡ��ؿ�
      &                   xyz_DensBZ,       &!���ܾ��̩��
      &                   pyz_DensBZ,       &!���ܾ��̩��
      &                   xqz_DensBZ,       &!���ܾ��̩��
      &                   xyr_DensBZ         !���ܾ��̩��
    use average,   only : pqz_xyz, pyz_xyz, &
      &                   xqz_xyz, xqr_xyz, &
      &                   xyr_xyz, pyr_xyz
    use differentiate_center2,              &
      &            only : xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                   pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
      &                   pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
      &                   pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
      &                   xyr_dx_pyr, xyr_dy_xqr, xqr_dz_xqz, &
      &                   xqz_dx_pqz, xqr_dy_xyr, xqz_dz_xqr  
    use DExnerDt,  only : xyz_DExnerDt_xyz

    ! ���ۤη�����ػ�
    !
    implicit none

    ! �������ѿ� 
    !
    real(DP),intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ʿ®��
    real(DP),intent(in)    :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ʿ®��
    real(DP),intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ľ®��
    real(DP),intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !����
    real(DP),intent(inout) :: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !�����顼�̤ο�ʿ��ή�Ȼ�
    real(DP),intent(inout) :: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !�����顼�̤ο�ʿ��ή�Ȼ�
    real(DP),intent(inout) :: xyr_DVelZDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !�����顼�̤ο�ʿ��ή�Ȼ�
    real(DP),intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(out)   :: xyz_KmAl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή�Ȼ�����
    real(DP),intent(out)   :: xyz_KhAl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή�Ȼ�����

    ! ����ѿ� 
    !
    real(DP)            :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή�Ȼ�����
    real(DP)            :: xyz_KhBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή�Ȼ�����
    real(DP)            :: xyz_DispPI(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή���ͥ륮���ξû�
    real(DP)            :: xyz_DispHeat(imin:imax,jmin:jmax,kmin:kmax)
                                                    !��ή���ͥ륮���ξû�
    real(DP)            :: xyz_Turb(imin:imax,jmin:jmax,kmin:kmax)
                                                    !
    real(DP)            :: pyz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xqz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyr_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyz_DensBZKhBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)            :: xyz_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)


    !----------------------------------
    ! �Ȼ������λ���ȯŸ�ϲ򤫤ʤ�
    !
    xyz_KmBl = ConstKm
    xyz_KhBl = ConstKh
    xyz_KmAl = ConstKm
    xyz_KhAl = ConstKh
    
    !---------------------------------------------------------
    ! ����ѿ������
    !
    xyz_DensBZKhBl      = xyz_DensBZ * xyz_KhBl
    xyz_DensBZKmBl      = xyz_DensBZ * xyz_KmBl
    
    !--------------------------------
    ! ���̤� tendency
    !
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_DensBZKhBl )                          &
      &   * xyr_dz_xyz( xyz_PTempBl ) ) / xyz_DensBZ

    xyz_DispHeat = (xyz_KmBl ** 3.0d0)                              &
      & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0)) &
      & * FactorDispHeat 

    xyz_DPTempDt = xyz_DPTempDt + xyz_Turb + xyz_DispHeat

    call HistoryAutoPut(TimeN, 'DPTempDtDisp', xyz_DispHeat(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtTurb', xyz_Turb(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! Turb.u
    !
    pyz_Turb =                                                      &
      &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) ) &
      & + pyz_dy_pqz(                                               &
      &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )        &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )        &
      &   )                                                         &
      & + pyz_dz_pyr(                                               &
      &       pyr_xyz( xyz_DensBZKmBl ) * pyr_dx_xyr( xyr_VelZBl )  &
      &     + pyr_xyz( xyz_DensBZKmBl ) * pyr_dz_pyz( pyz_VelXBl )  &
      &   ) / pyz_DensBZ 
!      & - 2.0d0 * pyz_dx_xyz( ( xyz_KmBl ** 2.0d0 ) )              &
!      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )

    pyz_DVelXDt = pyz_DVelXDt + pyz_Turb

    call HistoryAutoPut(TimeN, 'DVelXDtTurb', pyz_Turb(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! Turb.v
    !
    xqz_Turb =                                                      &
      &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) ) &
      & + xqz_dx_pqz(                                               &
      &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )        &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )        &
      &   )                                                         &
      & + xqz_dz_xqr(                                               &
      &       xqr_xyz( xyz_DensBZKmBl ) * xqr_dy_xyr( xyr_VelZBl )  &
      &     + xqr_xyz( xyz_DensBZKmBl ) * xqr_dz_xqz( xqz_VelYBl )  &
      &   ) / xqz_DensBZ
!      & - 2.0d0 * xqz_dy_xyz( ( xyz_KmBl ** 2.0d0 ) )             &
!      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    
    xqz_DVelYDt = xqz_DVelYDt + xqz_Turb

    call HistoryAutoPut(TimeN, 'DVelYDtTurb', xqz_Turb(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! Turb.w
    !
    xyr_Turb = &
      & + 2.0d0 * xyr_dz_xyz( xyz_DensBZKmBl * xyz_dz_xyr( xyr_VelZBl ) ) &
      &    / xyr_DensBZ                                                   &
      & + xyr_dx_pyr(                                                     &
      &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )               &
      &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )               &
      &   )                                                               &
      & + xyr_dy_xqr(                                                     &
      &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )               &
      &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )               &
      &   )                                                          
!      & - 2.0d0 * xyr_dz_xyz( xyz_DensBZ * ( xyz_KmBl ** 2.0d0 ) )  &
!      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )       &
!      &   / xyr_xyz( xyz_DensBZ ) 

    xyr_DVelZDt = xyr_DVelZDt + xyr_Turb

    call HistoryAutoPut(TimeN, 'DVelZDtTurb', xyr_Turb(1:nx, 1:ny, 1:nz))

    !--------------------
    ! Exner function
    !
    xyz_DispPI = xyz_DExnerDt_xyz( xyz_DispHeat ) * FactorDExnerDtTurb
    xyz_DExnerDt = xyz_DExnerDt + xyz_DispPI
    
    call HistoryAutoPut(TimeN, 'DExnerDtDisp', xyz_DispPI(1:nx, 1:ny, 1:nz))
    
  end subroutine Turbulence_constKm_forcing

!!!------------------------------------------------------------------------!!!

  subroutine turbulence_constkm_output

    ! �⥸�塼��
    !
    use gtool_historyauto, only : HistoryAutoAddVariable

    !���ۤη�����ػ�
    !
    implicit none

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (x)', &
      & units='m.s-2',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (y)', &
      & units='m.s-2',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (z)', &
      & units='m.s-2',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of exner function', &
      & units='s-1',    &
      & xtype='float')

  end subroutine turbulence_constkm_output


end module Turbulence_constKm


