!= Module HeatFluxConst
!
! Authors::   ODAKA Masatsugu
! Version::   $Id: surfaceflux_const.f90,v 1.5 2014/03/04 04:49:42 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2012. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Surfaceflux_const
  !
  ! ����, ������������ΰ����Ǯ����ư�̡�ʪ���ե�å�����Ϳ�������ˤ��
  ! ��ư��, ����, �ŷ���ʬ���Ѳ�Ψ��׻�����⥸�塼��. 
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP, STRING

  !���ۤη�����ػ�
  implicit none

  !°���λ���
  private

  !�ؿ��� public ������
  public surfaceflux_const_init
  public surfaceflux_const_forcing

  !�ѿ����
  real(DP), save :: SfcXMomFluxBtm = 0.0d0   ! X �����α�ư�̥ե�å��� (��������) 
  real(DP), save :: SfcXMomFluxTop = 0.0d0   ! X �����α�ư�̥ե�å��� (��������) 
  real(DP), save :: SfcYMomFluxBtm = 0.0d0   ! Y �����α�ư�̥ե�å��� (��������) 
  real(DP), save :: SfcYMomFluxTop = 0.0d0   ! Y �����α�ư�̥ե�å��� (��������) 
  real(DP), save :: SfcHeatFluxBtm = 0.0d0   ! Ǯ�ե�å��� (��������) 
  real(DP), save :: SfcHeatFluxTop = 0.0d0   ! Ǯ�ե�å��� (��������) 
  real(DP), save :: SfcQmixFluxBtm = 0.0d0   ! ʪ���ե�å��� (��������) 
  real(DP), save :: SfcQmixFluxTop = 0.0d0   ! ʪ���ե�å��� (��������) 

  character(STRING), parameter:: module_name = 'surfaceflux_const'
                              ! �⥸�塼���̾��.
                              ! Module name
  real(DP), save :: FactorDExnerDtSurf = 1.0d0
 

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_const_init
    !
    !NAMELIST ����ɬ�פʾ�����ɤ߼��, �����ѿ��������Ԥ�. 
    !

    !�⥸�塼���ɤ߹���
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable    
    use gridset,           only : ncmax
    use composition,       only : SpcWetSymbol
    use namelist_util,     only : namelist_filename

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    integer    :: l, unit

    !---------------------------------------------------------------
    ! NAMELIST �����������
    !
    NAMELIST /surfaceflux_const_nml/ &
      & SfcXMomFluxBtm, SfcXMomFluxTop, SfcYMomFluxBtm, SfcYMomFluxTop, &
      & SfcHeatFluxBtm, SfcHeatFluxTop, SfcQmixFluxBtm, SfcQmixFluxTop, &
      & FactorDExnerDtSurf 

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_const_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "SfcXMomFluxBtm = %f", &
      &                  d=(/SfcXMomFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcXMomFluxTop = %f", &
      &                  d=(/SfcXMomFluxTop/) )
    call MessageNotify( "M", module_name, "SfcYMomFluxBtm = %f", &
      &                  d=(/SfcYMomFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcYMomFluxTop = %f", &
      &                  d=(/SfcYMomFluxTop/) )
    call MessageNotify( "M", module_name, "SfcHeatFluxBtm = %f", &
      &                  d=(/SfcHeatFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcHeatFluxTop = %f", &
      &                  d=(/SfcHeatFluxTop/) )
    call MessageNotify( "M", module_name, "SfcQmixFluxBtm = %f", &
      &                  d=(/SfcQmixFluxBtm/) )
    call MessageNotify( "M", module_name, "SfcQmixFluxTop = %f", &
      &                  d=(/SfcQmixFluxTop/) )
    call MessageNotify( 'M', module_name, "FactorDExnerDtSurf = %f", &
      &                  d=(/FactorDExnerDtSurf/) )

    !---------------------------------------------------------------
    ! �ե�������Ϥ����
    !
    call HistoryAutoAddVariable(      &
      & varname='PTempSfcFlux',       &
      & dims=(/'x','y','z','t'/),         &
      & longname='surface potential temperature flux (heat flux divided by density and specific heat)', &
      & units='K.m.s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(      &
      & varname='ExnerSfcFlux',       &
      & dims=(/'x','y','z','t'/),         &
      & longname='surface exner function flux (heat flux divided by density and specific heat)', &
      & units='s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfcFlux',    &
      & dims=(/'x','y','z','t'/),     &
      & longname='surface flux of x-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfcFlux',    &
      & dims=(/'x','y','z','t'/),     &
      & longname='surface flux of y-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname=trim(SpcWetSymbol(l))//'SfcFlux', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='surface flux of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio (mass flux divided by density)',  &
        & units='m.s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(  &
      & varname='DPTempDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='potential temperature tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='exner function tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='x-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='y-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfc', & 
        & dims=(/'x','y','z','t'/),     &
        & longname=trim(SpcWetSymbol(l))//' mixing ratio tendency by surface flux',  &
        & units='s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(       &
      & varname='SfcHeatFlux',         &
      & dims=(/'x','y','z','t'/),          &
      & longname='surface heat flux',  &
      & units='W.m-2',                 &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcXMomFlux',                        &
      & dims=(/'x','y','z','t'/),                         &
      & longname='surface x-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcYMomFlux',                        &
      & dims=(/'x','y','z','t'/),                         &
      & longname='surface y-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(                               &
        & varname=trim(SpcWetSymbol(l))//'_SfcMassFlux',         &
        & dims=(/'x','y','z','t'/),                              &
        & longname=trim(SpcWetSymbol(l))//' surface mass flux',  &
        & units='kg.m-2.s-1',                                    &
        & xtype='float')
    end do

  end subroutine Surfaceflux_Const_init


!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Const_forcing( &
    &   pyz_DVelXDt, xqz_DVelYDt, xyz_DPTempDt, xyz_DExnerDt, xyzf_DQMixDt &
    & )
    ! 

    !�⥸�塼���ɤ߹���
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use gridset,           only : imin,         & !x ����������β���
      &                           imax,         & !x ����������ξ��
      &                           jmin,         & !y ����������β���
      &                           jmax,         & !y ����������ξ��
      &                           kmin,         & !z ����������β���
      &                           kmax,         & !z ����������ξ��
      &                           nx, ny, nz, ncmax
    use axesset,           only : xyz_dz          !z �����γʻ����ֳ�
!    use xyz_base_module,   only : xyz_pyz,  &
!      &                           xyz_xqz,  &
!      &                           pyz_xyz,  &
!      &                           xqz_xyz
    use basicset,          only : xyz_DensBZ     !���ܾ��̩��
    use constants,         only : CpDry
    use composition,       only : SpcWetSymbol
    use timeset,           only : TimeN
    use DExnerDt,          only : xyz_DExnerDt_xyz_xyzf

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP) :: pyr_VelXflux (imin:imax,jmin:jmax,kmin:kmax) !®�٥ե�å���
    real(DP) :: xqr_VelYflux (imin:imax,jmin:jmax,kmin:kmax) !®�٥ե�å���
    real(DP) :: xyr_PTempFlux(imin:imax,jmin:jmax,kmin:kmax) !���٥ե�å���
    real(DP) :: xyr_ExnerFlux(imin:imax,jmin:jmax,kmin:kmax) !Exner �ؿ�������
    real(DP) :: xyrf_QMixFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                             !ʪ���ե�å���
    real(DP) :: pyr_SfcXMomFlux(imin:imax,jmin:jmax,kmin:kmax) !��ư�̥ե�å���
    real(DP) :: xqr_SfcYMomFlux(imin:imax,jmin:jmax,kmin:kmax) !��ư�̥ե�å���
    real(DP) :: xyr_SfcHeatFlux(imin:imax,jmin:jmax,kmin:kmax) !Ǯ�ե�å���

    real(DP) :: pyz_DVelXDtFlux (imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtFlux (imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DPTempDtFlux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtFlux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyzf_DQMixDtFlux(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    integer  :: s             !�롼���ѿ�

!-----

    ! �����
    xyr_SfcHeatFlux = 0.0d0
    xyr_SfcHeatFlux = 0.0d0
    
    pyr_SfcXMomFlux = 0.0d0
    pyr_SfcXMomFlux = 0.0d0

    xqr_SfcYMomFlux = 0.0d0
    xqr_SfcYMomFlux = 0.0d0

    xyrf_QMixFlux = 0.0d0
    xyrf_QMixFlux = 0.0d0
    

    ! NAMLIST �����ꤷ���ͤ�����
    !
    xyr_SfcHeatFlux(:,:,1)  = SfcHeatFluxBtm
    xyr_SfcHeatFlux(:,:,nz) = SfcHeatFluxTop
    
    pyr_SfcXMomFlux(:,:,1)  = SfcXMomFluxBtm
    pyr_SfcXMomFlux(:,:,nz) = SfcXMomFluxTop

    xqr_SfcYMomFlux(:,:,1)  = SfcYMomFluxBtm
    xqr_SfcYMomFlux(:,:,nz) = SfcYMomFluxTop

    xyrf_QMixFlux(:,:,1,:)  = SfcQmixFluxBtm
    xyrf_QMixFlux(:,:,nz,:) = SfcQmixFluxTop

    ! ®�٥ե�å���, ���٥ե�å����ؤ��Ѵ�
    ! 
    pyr_VelXFlux  = pyr_SfcXMomFlux / xyz_DensBz
    xqr_VelYFlux  = xqr_SfcYMomFlux / xyz_DensBz
    xyr_PtempFlux = xyr_SfcHeatFlux / xyz_DensBz / CpDry

    ! Exner �ؿ�ͽ�󼰤ζ�����η׻�
    !
    xyr_ExnerFlux = xyz_DExnerDt_xyz_xyzf(xyr_PTempFlux, xyrf_QMixFlux) * FactorDExnerDtSurf 

    ! �����Ѳ�Ψ�η׻�
    !
    xyz_DPTempDtFlux = xyr_PTempFlux / xyz_dz
    xyz_DExnerDtFlux = xyr_ExnerFlux / xyz_dz
    do s = 1, ncmax
      xyzf_DQMixDtFlux(:,:,:,s) = xyrf_QMixFlux(:,:,:,s) / xyz_dz
    end do
    pyz_DVelXDtFlux  = pyr_VelXFlux / xyz_dz
    xqz_DVelYDtFlux  = xqr_VelYFlux / xyz_dz

    ! ����������ؤβû�
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtFlux
    xyz_DExnerDt = xyz_DExnerDt + xyz_DExnerDtFlux
    do s = 1, ncmax
      xyzf_DQMixDt(:,:,:,s) = xyzf_DQMixDt(:,:,:,s) + xyzf_DQMixDtFlux(:,:,:,s)
    end do
    pyz_DVelXDt = pyz_DVelXDt + pyz_DVelXDtFlux
    xqz_DVelYDt = xqz_DVelYDt + xqz_DVelYDtFlux

    ! ����
    ! 
    call HistoryAutoPut(TimeN, 'DPTempDtSfc', xyz_DPTempDtFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtSfc', xyz_DExnerDtFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSfc',  pyz_DVelXDtFlux (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSfc',  xqz_DVelYDtFlux (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfc', &
        & xyzf_DQMixDtFlux(1:nx,1:ny,1:nz,s))
    end do

    call HistoryAutoPut(TimeN, 'PTempSfcFlux', xyr_PTempFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerSfcFlux', xyr_ExnerFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelXSfcFlux',  pyr_VelXFlux (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelYSfcFlux',  xqr_VelYFlux (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s))//'SfcFlux', &
        & xyrf_QMixFlux(1:nx,1:ny,1:nz,s))
    end do

    call HistoryAutoPut(TimeN, 'SfcHeatFlux', xyr_SfcHeatFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', pyr_SfcXMomFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', xqr_SfcYMomFlux(1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s))//'SfcMassFlux', &
        & xyz_DensBZ(1:nx,1:ny,1:nz) * xyrf_QMixFlux(1:nx,1:ny,1:nz,s))
    end do

  end subroutine Surfaceflux_Const_forcing
  
end module Surfaceflux_Const
