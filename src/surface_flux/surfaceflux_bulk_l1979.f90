!= Module HeatFlux
!
! Authors::   ODAKA Masatsugu
! Version::   $Id: surfaceflux_bulk_l1979.f90,v 1.1 2012/07/30 08:08:29 odakker Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module Surfaceflux_bulk_l1979
  !
  != $B2<It6-3&$G$N%U%i%C%/%9$N7W;;%b%8%e!<%k(B
  !
  != Surface flux
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! Louis et al. (1979) $B$NJ}K!$K4p$E$$$FCOI=LL%U%i%C%/%9$r7W;;(B. 
  !
  ! Surface fluxes are calculated by using the scheme by Louis et al. (1979). 
  !
  !== References
  !
  ! Louis, J-F., 1979:
  ! A parametric model of vertical eddy fluxes in the atmosphere. 
  ! Boundary Layer Meteorology, 17, 187-202.
  !== Variable List
  !
  !== Procedures List
  !

  ! $B%b%8%e!<%k0zMQ(B; USE statement
  !  

  ! GTOOL $BJQ?t$H<jB3$-(B
  ! GTOOL variables and procedures
  !
  use dc_types, only: DP, STRING
  use dc_iounit, only: FileOpen
  use dc_message, only: MessageNotify
  use gtool_historyauto, only: HistoryAutoAddVariable, HistoryAutoPut

  ! $BJBNs=hM}MQJQ?t(B
  ! Parallel processing variable 
  !       

  ! $B3J;RE@@_Dj(B
  ! Grid points settings
  !
  use gridset, only: imin, & ! x $BJ}8~G[Ns2<8B(B 
    &                      & ! Upper limit of array in x
    &                imax, & ! x $BJ}8~G[Ns>e8B(B 
    &                      & ! Lower limit of array in x
    &                jmin, & ! y $BJ}8~G[Ns2<8B(B 
    &                      & ! Upper limit of array in y
    &                jmax, & ! y $BJ}8~G[Ns>e8B(B 
    &                      & ! Lower limit of array in y
    &                kmin, & ! z $BJ}8~G[Ns2<8B(B 
    &                      & ! Upper limit of array in z
    &                kmax, & ! z $BJ}8~G[Ns>e8B(B 
    &                      & ! Lower limit of array in z
    &                nx,   & ! x $BJ}8~3J;RE@?t(B 
    &                      & ! Number of grid point in x
    &                ny,   & ! y $BJ}8~3J;RE@?t(B 
    &                      & ! Number of grid point in y
    &                nz,   & ! z $BJ}8~3J;RE@?t(B 
    &                      & ! Number of grid point in z
    &                ncmax   ! $BAH@.$N?t(B       
                             ! Number of spices

  ! $B:BI8<4$H1i;;;R(B
  ! Axes and operator settings
  !
  use axesset, only: z_dz,        & ! z $BJ}8~3J;RE@4V3V(B 
    &                             & ! Grid size in z
    &                xyz_avr_pyz, & ! $BJ?6QA`:n(B
    &                             & ! Average operator
    &                xyz_avr_xqz, & ! $BJ?6QA`:n(B
    &                             & ! Average operator
    &                pyz_avr_xyz, & ! $BJ?6QA`:n(B
    &                             & ! Average operator
    &                xqz_avr_xyz    ! $BJ?6QA`:n(B
                                    ! Average operator

  ! $B4pK\>lJQ?t(B
  ! Basic state variables
  !
  use basicset, only: xyz_ExnerBZ,  & ! $B05NO4X?t(B
    &                               & ! Exner function
    &                 xyz_PressBZ,  & ! $B05NO(B
    &                               & ! Pressure
    &                 xyz_PTempBZ,  & ! $B290L(B
    &                               & ! Potential temperature
    &                 xyz_TempBZ,   & ! $B29EY(B
    &                               & ! Temperature
    &                 xyzf_QMixBZ,  & ! $B:.9gHf(B
    &                               & ! Mixing ration
    &                 xyz_DensBZ      ! $BL)EY(B    
                                      ! Density

  ! $BDj?t(B
  ! Constatns
  !
  use constants, only: Grav,       & ! $B=ENO2CB.EY(B
    &                              & ! Gravity
    &                  MolWtDry,   & ! $B4%Ag6u5$J,;RNL(B 
    &                              & ! Molecular weight of dry air
    &                  PressBasis, & ! $B290L$N4p=`05NO(B
    &                              & ! Reference pressure
    &                  TempSfc,    & ! $BCOI=LL29EY(B
    &                              & ! Surface temperature
    &                  PressSfc,   & ! $BCOI=LL5$05(B
    &                              & ! Surface pressure
    &                  CpDry,      & ! $B4%Ag6u5$Dj05HfG.(B 
    &                              & ! Specific heat of dry air
    &                  GasRDry       ! $B4%Ag6u5$5$BNDj?t(B
                                     ! Gas constant of dry air

  use constants0, only: FKarm        ! $B%+%k%^%sDj?t(B
                                     ! Karmann constant

  ! $BAH@.$K4X$9$kJQ?t(B
  ! Setting for atmospheric composition
  ! 
  use composition, only: IdxCG,       & ! $B6E7k2aDx(B($B5$BN(B)$B$NG[NsE:;z(B 
    &                                 & ! Index of vapor
    &                    IdxCC,       & ! $B6E7k2aDx(B($B1@(B)$B$NG[NsE:;z(B 
    &                                 & ! Index of cloud
    &                    SpcWetID,    & !$B6E7k@.J,$N2=3X<o(BID
    &                                 & ! ID number of moist condensation spices
    &                    CondNum,     & ! $B1@$N?t(B
    &                                 & ! Number of cloud component
    &                    MolWtWet,    & ! $B6E7k@.J,J,;RNL(B
    &                                 & ! Molecular weight of moist air
    &                    SpcWetSymbol   ! $B6E7k@.J,$N2=3X<oL>(B
                                        ! Chemical formula of condensation spices

  ! $B2=3XNL7W;;$K4X$9$k@_Dj(B
  ! Settings for chemical calculation
  !
  use chemcalc, only: SvapPress ! $BK0OB>x5$05(B
                                ! Saturation vapor pressure

  ! NAMELIST $B$K4X$9$k@_Dj(B
  ! Settings for NAMELIST
  !
  use namelist_util, only: namelist_filename ! NAMELIST $B%U%!%$%kL>(B
                                             ! NAMELIST file name
  ! $B;~9o$K4X$9$k@_Dj(B
  ! Setting for time
  !
  use timeset, only:  TimeN ! $B;~9o(B t 
                            ! Time "t"

  ! $B;~4VJQ2=9`(B
  ! Tendency of variable
  use DExnerDt, only: xy_DExnerDt_xy_xyf ! D$\pi$/Dt


  ! $B0EL[$N7?@k8@6X;_(B
  ! Implicit none
  !
  implicit none

  ! $BB0@-$N;XDj(B
  !
  private

  ! $B8x3+<jB3$-(B
  ! Public procedure
  !
  public surfaceflux_bulk_init
  public surfaceflux_bulk_forcing

  ! $BHs8x3+JQ?t(B
  ! Privete variables
  ! 
  logical, save:: FlagConstBulkCoef
                            ! Flag for using constant bulk coefficient
  logical, save:: FlagUseOfBulkCoefInNeutralCond
                            ! Flag for using bulk coefficient in neutral condition
  real(DP), save:: ConstBulkCoef
                            ! $B%P%k%/78?t0lDjCM(B. 
                            ! Steady value of bulk coefficient
  real(DP), save :: VelMinForRi = 1.0d-8 
                            ! $B%j%A%c!<%I?t7W;;MQB.EY2<8BCM(B
                            ! Lower limit of velocity for Ri
  real(DP), save :: SfcRoughLength = 1.0d-2
                            ! $BADEYD9$5(B
                            ! Roughness length
  real(DP), save :: VelBulkCoefMin = 0.0d0
                            ! $ u $ $B%P%k%/78?t:G>.CM(B. 
                            ! Minimum value of $ u $ bulk coefficient
  real(DP), save :: TempBulkCoefMin = 0.0d0
                            ! $ T $ $B%P%k%/78?t:G>.CM(B. 
                            ! Minimum value of $ T $ bulk coefficient
  real(DP), save :: QmixBulkCoefMin = 0.0d0
                            ! $ q $ $B%P%k%/78?t:G>.CM(B. 
                            ! Minimum value of $ q $ bulk coefficient
  real(DP), save :: VelBulkCoefMax = 1.0d2
                            ! $ u $ $B%P%k%/78?t:GBgCM(B. 
                            ! Maximum value of $ u $ bulk coefficient
  real(DP), save :: TempBulkCoefMax = 1.0d2
                            ! $ T $ $B%P%k%/78?t:GBgCM(B. 
                            ! Maximum value of $ T $ bulk coefficient
  real(DP), save :: QmixBulkCoefMax = 1.0d2
                            ! $ q $ $B%P%k%/78?t:GBgCM(B. 
                            ! Maximum value of $ q $ bulk coefficient

  real(DP), save :: DPTempDtFluxMin = 0.0d0
                            ! $B290L%U%i%C%/%9:G>.CM(B. 
                            ! Minimum value of potential temp. flux
  real(DP), save :: ExnerFluxMin = 0.0d0
                            ! $B05NO4X?t%U%i%C%/%9:G>.CM(B. 
                            ! Minimum value of exner function flux
  real(DP), save :: QmixFluxMin = 0.0d0
                            ! $B:.9gHf%U%i%C%/%9:G>.CM(B. 
                            ! Minimum value of mixing ratio flux
  real(DP), save  :: Vel0 = 0.0d0  ! $B2<AX$G$N?eJ?B.EY?s>e$2CM(B
                                   ! 

  character(*), parameter:: module_name = 'surfaceflux_bulk_l1979'
                                   ! $B%b%8%e!<%k$NL>>N(B.
                                   ! Module name

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Bulk_init
    !
    ! NAMELIST $B$+$iI,MW$J>pJs$rFI$_<h$j(B, $B;~4V4XO"$NJQ?t$N@_Dj$r9T$&(B. 
    !

    ! $B0EL[$N7?@k8@6X;_(B
    ! Implicit none
    !
    implicit none

    ! $B:n6HJQ?t(B
    ! Work variables
    !
    integer    :: l,   & ! $BAH@.J}8~$K2s$k(B DO $B%k!<%WMQ:n6HJQ?t(B
      &                & ! Work variables for DO loop in dimension of constituents
      &           unit   ! $B=PNOAuCVHV9f(B
                         ! Device number 

    !---------------------------------------------------------------
    ! NAMELIST $B$+$i>pJs$r<hF@(B
    !
    NAMELIST /surfaceflux_bulk_nml/ &
      &  FlagConstBulkCoef,                                 &
      &  FlagUseOfBulkCoefInNeutralCond, ConstBulkCoef,     &
      ! 
      &  VelMinForRi, SfcRoughLength, Vel0,                 &
      !
      &  VelBulkCoefMin, TempBulkCoefMin, QmixBulkCoefMin,  &
      &  VelBulkCoefMax, TempBulkCoefMax, QmixBulkCoefMax             


    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_bulk_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "SfcRoughLength = %f",  &
      &  d=(/SfcRoughLength/))
    call MessageNotify( "M", module_name, "VelMinForRi = %f",     &
      &  d=(/VelMinForRi/) )
    call MessageNotify( "M", module_name, "Vel0 = %f", d=(/Vel0/))
    
    call MessageNotify( 'M', module_name, "FlagConstBulkCoef              = %b", l = (/ FlagConstBulkCoef /) )
    call MessageNotify( 'M', module_name, "FlagUseOfBulkCoefInNeutralCond = %b", l = (/ FlagUseOfBulkCoefInNeutralCond /) )
    call MessageNotify( 'M', module_name, "ConstBulkCoef   = %f", &
      &  d = (/ ConstBulkCoef   /) )
    call MessageNotify( 'M', module_name, "VelBulkCoefMin  = %f", &
      &  d = (/ VelBulkCoefMin  /) )
    call MessageNotify( 'M', module_name, "TempBulkCoefMin = %f", &
      &  d = (/ TempBulkCoefMin /) )
    call MessageNotify( 'M', module_name, "QmixBulkCoefMin = %f", &
      &  d = (/ QmixBulkCoefMin /) )
    call MessageNotify( "M", module_name, "VelBulkCoefMax = %f",  &
      &  d=(/VelBulkCoefMax/))
    call MessageNotify( "M", module_name, "TempBulkCoefMax = %f", &
      &  d=(/TempBulkCoefMax/))
    call MessageNotify( "M", module_name, "QmixBulkCoefMax = %f", &
      &  d=(/QmixBulkCoefMax/))

    call HistoryAutoAddVariable(      &
      & varname='PTempSfcFlux',       &
      & dims=(/'x','y','t'/),         &
      & longname='surface potential temperature flux (heat flux divided by density and specific heat)', &
      & units='K.m.s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(      &
      & varname='ExnerSfcFlux',       &
      & dims=(/'x','y','t'/),         &
      & longname='surface exner function flux (heat flux divided by density and specific heat)', &
      & units='s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfcFlux',    &
      & dims=(/'x','y','t'/),     &
      & longname='surface flux of x-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfcFlux',    &
      & dims=(/'x','y','t'/),     &
      & longname='surface flux of y-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfcFlux', & 
        & dims=(/'x','y','t'/),     &
        & longname='surface flux of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio (mass flux divided by density)',  &
        & units='m.s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(  &
      & varname='PTempSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='potential temperature tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='exner function tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='x-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfc',         &
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
      & dims=(/'x','y','t'/),          &
      & longname='surface heat flux',  &
      & units='W.m-2',                 &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcXMomFlux',                        &
      & dims=(/'x','y','t'/),                         &
      & longname='surface x-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcYMomFlux',                        &
      & dims=(/'x','y','t'/),                         &
      & longname='surface y-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(                               &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfcMassFlux',         &
        & dims=(/'x','y','t'/),                              &
        & longname=trim(SpcWetSymbol(l))//' surface mass flux',  &
        & units='kg.m-2.s-1',                                    &
        & xtype='float')
    end do

  end subroutine Surfaceflux_Bulk_init
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Bulk_forcing( &
    &   pyz_VelX, xqz_VelY, xyz_PTemp, xyz_Exner, xyzf_QMix, &
    &   pyz_DVelXDt, xqz_DVelYDt, xyz_DPTempDt, xyz_DExnerDt, xyzf_DQMixDt &
    & )
    ! 
    ! $B2<It6-3&$+$i$N%U%i%C%/%9$K$h$k29EY$NJQ2=N($r(B,
    ! $B%P%k%/J}K!$K4p$E$$$F7W;;$9$k(B.
    !

    ! $B0EL[$N7?@k8@6X;_(B
    ! Implicit none
    !
    implicit none

    ! $BJQ?t(B
    ! variables
    !
    real(DP), intent(in)   :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.(B
                              ! X-component velocity
    real(DP), intent(in)   :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.(B
                              ! Y-component velocity
    real(DP), intent(in)   :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B290L(B
                              ! Potential temperature
    real(DP), intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B05NO4X?t(B
                              ! Exner function
    real(DP), intent(in)   :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! $B:.9gHf(B
                              ! Mixing ration
    real(DP), intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! $BIwB.;~4VJQ2=N((B
                              ! X-component velocity tendency
    real(DP), intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! $BIwB.;~4VJQ2=N((B
                              ! Y-component velocity tendency
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B290L;~4VJQ2=N((B
                              ! Potential tempreture tendency
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B05NO4X?t;~4VJQ2=N((B
                              ! Exner function tendency
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! $B:.9gHf;~4VJQ2=N((B
                              ! Mixing ratio tendency

    ! $B:n6HJQ?t(B
    ! Work variables
    real(DP) :: xy_SurfBulkRiNum(imin:imax,jmin:jmax)
                              ! $B%P%k%/%j%A%c!<%I%=%s?t(B
                              ! Bulk Richardson number
    real(DP) :: xy_SurfRoughLength(imin:imax,jmin:jmax)
                              ! $BADEYD9$5(B
                              ! Roughness length
    real(DP) :: xy_SurfVelBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($B1?F0NL(B)
                              ! Bulk coefficient for momentum
    real(DP) :: xy_SurfTempBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($BG.(B)
                              ! Bulk coefficient for heat
    real(DP) :: xy_SurfQmixBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($B:.9gHf(B)
                              ! Bulk coefficient for mixing ratio
    real(DP) :: py_VelXflux (imin:imax,jmin:jmax)
                              ! x $BJ}8~B.EY%U%i%C%/%9(B
                              ! velocity flux in x
    real(DP) :: xq_VelYflux (imin:imax,jmin:jmax)
                              ! y $BJ}8~B.EY%U%i%C%/%9(B
                              ! celocity flux in y
    real(DP) :: xy_DPTempDtFlux(imin:imax,jmin:jmax)
                              ! $B290L%U%i%C%/%9(B
                              ! potential temperature flux
    real(DP) :: xy_ExnerFlux(imin:imax,jmin:jmax)
                              !
                              !
    real(DP) :: xyf_QMixFlux(imin:imax,jmin:jmax,ncmax)
                              ! $B6E7k@.J,:.9gHf%U%i%C%/%9(B
                              ! Mixing ratio flux
    real(DP) :: xyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.(B (xyz $B3J;R(B)
                              ! X-component velocity (xyz grid)
    real(DP) :: xyz_VelY(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.(B (xyz $B3J;R(B)
                              ! Y-component velocity (xyz grid)
    real(DP) :: xyz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.@dBPCM(B (xyz $B3J;R(B)
                              ! Absolute value of horizontal velocity (xyz grid)
    real(DP) :: pyz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.@dBPCM(B (pyz $B3J;R(B)
                              ! Absolute value of horizontal velocity (pyz grid)
    real(DP) :: xqz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B?eJ?IwB.@dBPCM(B (xqz $B3J;R(B)
                              ! Absolute value of horizontal velocity (xqz grid)
    real(DP) :: xyz_PTempAll(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B290L(B($B4pK\>l(B + $B>qMp(B)
                              ! Total value of potential temperature
    real(DP) :: xyz_ExnerAll (imin:imax,jmin:jmax,kmin:kmax)
                              ! $B05NO4X?t(B($B4pK\>l(B + $B>qMp(B)
                              ! Total value of exner function
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! $B:.9gHf(B($B4pK\>l(B + $B>qMp(B)
                              ! Total value of mixing ratios
    real(DP) :: xy_DPTempDtBulk(imin:imax,jmin:jmax)
                              ! $B;~4VJQ2=(B($B290L(B)
                              ! potential temperature tendency by surface flux
    real(DP) :: xy_DExnerDtBulk(imin:imax,jmin:jmax)
                              ! $B;~4VJQ2=(B($B05NO4X?t(B)
                              ! Exner function tendency by surface flux
    real(DP) :: xyf_DQMixDtBulk(imin:imax,jmin:jmax, ncmax)
                              ! $B;~4VJQ2=(B($B:.9gHf(B)
                              ! Mixing ratio tendency by surface flux
    real(DP) :: py_DVelXDtBulk (imin:imax,jmin:jmax)
                              ! $B;~4VJQ2=(B(U)
                              ! x-component velocity tendency by surface flux
    real(DP) :: xq_DVelYDtBulk (imin:imax,jmin:jmax)
                              ! $B;~4VJQ2=(B(V)
                              ! y-component velocity tendency by surface flux
    real(DP) :: xyz_DPTempDtBulk(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B;~4VJQ2=(B($B290L(B)
                              ! potential temperature tendency by surface flux
    real(DP) :: xyz_DExnerDtBulk(imin:imax,jmin:jmax,kmin:kmax)
                              ! $B;~4VJQ2=(B($B05NO4X?t(B)
                              ! Exner function tendency by surface flux
    real(DP) :: xyzf_DQMixDtBulk(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! $B;~4VJQ2=(B($B:.9gHf(B)
                              ! Mixing ratio tendency by surface flux
    real(DP) :: pyz_DVelXDtBulk (imin:imax,jmin:jmax,kmin:kmax)
                              ! $B;~4VJQ2=(B(U)
                              ! x-component velocity tendency by surface flux
    real(DP) :: xqz_DVelYDtBulk (imin:imax,jmin:jmax,kmin:kmax)
                              ! $B;~4VJQ2=(B(V)
                              ! y-component velocity tendency by surface flux
    real(DP) :: ExnerBZSfc    ! $BCOI=LL05NO4X?t(B 
                              ! Basic state Exner function at the surface
    real(DP) :: xy_PressSfc(imin:imax,jmin:jmax)
                              ! Total pressure at the surface
    integer  :: kz            ! $BG[NsE:;z(B
                              ! Arrzy index
    integer  :: s             ! $BAH@.J}8~$K2s$k(B DO $B%k!<%WMQ:n6HJQ?t(B
                              ! Work variables for DO loop in dimension of constituen                         

    ! $B=i4|2=(B
    ! Initialization
    ! 
    kz = 1

    ! $BADEYD9$5$N;XDj(B
    ! Specify surface length
    xy_SurfRoughLength = SfcRoughLength

    ! $BA4NL$N7W;;(B
    ! Calculate total value of thermodynamic variables
    ! 
    xyz_PTempAll  = xyz_PTemp + xyz_PTempBZ
    xyz_ExnerAll  = xyz_Exner + xyz_ExnerBZ
    xyzf_QMixAll  = xyzf_QMix + xyzf_QMixBZ

    ! Perturbation component of Exner function at the surface is assumed 
    ! to be same as that at the lowest layer. (YOT, 2011/09/03)
    ! 
    ExnerBZSfc    = (PressSfc / PressBasis) ** (GasRDry / CpDry)
    xy_PressSfc   = (PressBasis * xyz_ExnerAll(:,:,kz))**(CpDry / GasRDry)

    ! xyz $B3J;RE@$NB.EY$N7W;;(B
    ! Calculate velocities at xyz grid points
    !
    xyz_VelX = xyz_avr_pyz(pyz_VelX)
    xyz_VelY = xyz_avr_xqz(xqz_VelY)

    ! $B?eJ?IwB.$N@dBPCM$N7W;;(B
    ! Calculate of absoluto horizontal velocities
    ! 
    xyz_AbsVel = SQRT( xyz_VelX**2 + xyz_VelY**2 + Vel0**2 )
    pyz_AbsVel = pyz_avr_xyz(xyz_AbsVel)
    xqz_AbsVel = xqz_avr_xyz(xyz_AbsVel)

    ! $B%P%k%/(B $ R_i $ $B?t;;=P(B
    ! Calculate bulk $ R_i $
    !
    xy_SurfBulkRiNum =                                 &
      &   Grav / ( xyz_PTempAll(:,:,1) )               & 
      &   * ( xyz_PTempAll(:,:,1)                      &
      &     -TempSfc / (ExnerBZSfc + xyz_Exner(:,:,1)))&
      &   / max( xyz_AbsVel(:,:,1), VelMinForRi )**2   &
      &   * z_dz(1) * 0.5d0

    ! $B%P%k%/78?t$N7W;;(B
    ! Bulk coefficients are calculated.
    ! 
    call BulkCoef( &
      & xy_SurfBulkRiNum,    & ! (in)
      & xy_SurfRoughLength,  & ! (in)
      & xy_SurfVelBulkCoef,  & ! (out)
      & xy_SurfTempBulkCoef, & ! (out)
      & xy_SurfQmixBulkCoef  & ! (out)
      & )
    ! $B%U%i%C%/%9$N7W;;(B
    ! Surface fluxes are calculated. 
    !
    xy_DPTempDtFlux = - xy_SurfTempBulkCoef * xyz_AbsVel(:,:,kz)  &
      & * ( xyz_PTempAll(:,:,kz) - TempSfc / ( ExnerBZSfc + xyz_Exner(:,:,kz) ) )

    xyf_QMixFlux = 0.0d0
    do s = 1, CondNum
      xyf_QMixFlux(:,:,IdxCG(s)) =                                 &
        &     - xy_SurfQmixBulkCoef * xyz_AbsVel(:,:,kz)           &
        &       * (                                                &
        &            xyzf_QMixAll(:,:,kz,s)                        &
        &          - SvapPress( SpcWetID(IdxCC(s)), TempSfc )      &
        &             / ( xy_PressSfc )                            &
        &             * (MolWtWet(IdxCG(s)) / MolWtDry)            &
        &         )
    end do
    !
    xy_ExnerFlux = xy_DExnerDt_xy_xyf(xy_DPTempDtFlux, xyf_QMixFlux, kz)
    !
    py_VelXFlux = - xy_SurfVelBulkCoef * pyz_AbsVel(:,:,kz) * pyz_VelX(:,:,kz)
    !
    xq_VelYFlux = - xy_SurfVelBulkCoef * xqz_AbsVel(:,:,kz) * xqz_VelY(:,:,kz)


    ! $B%U%i%C%/%9$N2<8BCM$r@_Dj(B
    ! Set lower limit of surface fluxes
    ! 
    xy_DPTempDtFlux = max( DPTempDtFluxMin, xy_DPTempDtFlux )
    xy_ExnerFlux = max( ExnerFluxMin, xy_ExnerFlux )
    xyf_QMixFlux = max( QMixFluxMin, xyf_QMixFlux )


    ! $BCOI=%U%i%C%/%9$K$h$k;~4VJQ2=$r7W;;(B
    ! Tendencies by surface fluxes (convergences of fluxes) are calculated.
    !
    xy_DPTempDtBulk = - ( 0.0d0 - xy_DPTempDtFlux ) / z_dz(kz)
    !
    xy_DExnerDtBulk = - ( 0.0d0 - xy_ExnerFlux ) / z_dz(kz)
    !
    xyf_DQMixDtBulk = - ( 0.0d0 - xyf_QMixFlux ) / z_dz(kz)
    !
    py_DVelXDtBulk = - ( 0.0d0 - py_VelXFlux ) / z_dz(kz)
    !
    xq_DVelYDtBulk = - ( 0.0d0 - xq_VelYFlux ) / z_dz(kz)

    
    ! $B2>0z?tG[Ns$X3JG<(B
    ! Add tendency by surface flux convergence
    !
    xyz_DPTempDt(:,:,kz) = xyz_DPTempDt(:,:,kz) + xy_DPTempDtBulk
    xyz_DExnerDt(:,:,kz) = xyz_DExnerDt(:,:,kz) + xy_DExnerDtBulk
    do s = 1, ncmax
      xyzf_DQMixDt(:,:,kz,s) = xyzf_DQMixDt(:,:,kz,s) + xyf_DQMixDtBulk(:,:,s)
    end do
    pyz_DVelXDt (:,:,kz) = pyz_DVelXDt (:,:,kz) + py_DVelXDtBulk
    xqz_DVelYDt (:,:,kz) = xqz_DVelYDt (:,:,kz) + xq_DVelYDtBulk

    ! $B=PNO(B
    ! Output
    !
    xyz_DPTempDtBulk = 0.0d0
    xyz_DExnerDtBulk = 0.0d0
    xyzf_DQMixDtBulk = 0.0d0
    pyz_DVelXDtBulk  = 0.0d0
    xqz_DVelYDtBulk  = 0.0d0

    xyz_DPTempDtBulk(:,:,kz) = xy_DPTempDtBulk
    xyz_DExnerDtBulk(:,:,kz) = xy_DExnerDtBulk
    do s = 1, ncmax
      xyzf_DQMixDtBulk(:,:,kz,s) = xyf_DQMixDtBulk(:,:,s)
    end do
    pyz_DVelXDtBulk (:,:,kz) = py_DVelXDtBulk
    xqz_DVelYDtBulk (:,:,kz) = xq_DVelYDtBulk
    !
    call HistoryAutoPut(TimeN, 'PTempSfc', xyz_DPTempDtBulk(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerSfc', xyz_DExnerDtBulk(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelXSfc',  pyz_DVelXDtBulk (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelYSfc',  xqz_DVelYDtBulk (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfc', &
        & xyzf_DQMixDtBulk(1:nx,1:ny,1:nz,s))
    end do


    call HistoryAutoPut(TimeN, 'PTempSfcFlux', xy_DPTempDtFlux(1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'ExnerSfcFlux', xy_ExnerFlux(1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'VelXSfcFlux',  py_VelXFlux (1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'VelYSfcFlux',  xq_VelYFlux (1:nx,1:ny))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfcFlux', &
        & xyf_QMixFlux(1:nx,1:ny,s))
    end do

    call HistoryAutoPut(TimeN, 'SfcHeatFlux', &
      & CpDry * xyz_DensBZ(1:nx,1:ny,1) * xy_DPTempDtFlux(1:nx,1:ny) &
      &   * ( ExnerBZSfc + xyz_Exner(1:nx,1:ny,1) ) )
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', &
      & xyz_DensBZ(1:nx,1:ny,1) * py_VelXFlux (1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'SfcYMomFlux', &
      & xyz_DensBZ(1:nx,1:ny,1) * xq_VelYFlux (1:nx,1:ny))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfcMassFlux', &
        & xyz_DensBZ(1:nx,1:ny,1) * xyf_QMixFlux(1:nx,1:ny,s))
    end do

  end subroutine Surfaceflux_Bulk_forcing
!!!------------------------------------------------------------------------!!!
  subroutine BulkCoef(     & 
    & xy_SurfBulkRiNum,    & ! (in)
    & xy_SurfRoughLength,  & ! (in)
    & xy_SurfVelBulkCoef,  & ! (out)
    & xy_SurfTempBulkCoef, & ! (out)
    & xy_SurfQmixBulkCoef  & ! (out)
    & )


    implicit none

    real(DP),intent(in) :: xy_SurfBulkRiNum(imin:imax,jmin:jmax)
                              ! $B%P%k%/%j%A%c!<%I%=%s?t(B
                              ! Bulk Richardson number
    real(DP),intent(in) :: xy_SurfRoughLength(imin:imax,jmin:jmax)
                              ! $BADEYD9$5(B
                              ! Roughness length
    real(DP),intent(out) :: xy_SurfVelBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($B1?F0NL(B)
                              ! Bulk coefficient for momentum
    real(DP),intent(out) :: xy_SurfTempBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($BG.(B)
                              ! Bulk coefficient for heat
    real(DP),intent(out) :: xy_SurfQmixBulkCoef(imin:imax,jmin:jmax)
                              ! $B%P%k%/78?t(B($B:.9gHf(B)
                              ! Bulk coefficient for mixing ratio

    ! $B:n6HJQ?t(B
    ! Work variables
    !
    real(DP) :: xy_SurfBulkCoefInNeutCond(imin:imax, jmin:jmax)

    integer:: i               ! $B7PEYJ}8~$K2s$k(B DO $B%k!<%WMQ:n6HJQ?t(B
                              ! Work variables for DO loop in longitude
    integer:: j               ! $B0^EYJ}8~$K2s$k(B DO $B%k!<%WMQ:n6HJQ?t(B
                              ! Work variables for DO loop in latitude

    ! $B<B9TJ8(B ; Executable statement
    !
    if ( FlagConstBulkCoef ) then

      ! Use of constant bulk coefficient
      !

      xy_SurfVelBulkCoef  = ConstBulkCoef
      xy_SurfTempBulkCoef = ConstBulkCoef
      xy_SurfQmixBulkCoef = ConstBulkCoef

    else

      ! Parameterization by Louis et al. (1981)
      !
      ! $BCfN)%P%k%/78?t$N7W;;(B
      ! Calculate bulk coefficient in neutral condition
      !
      xy_SurfBulkCoefInNeutCond  = &
        & ( FKarm &
        & / log ( z_dz(1) / xy_SurfRoughLength ) )**2

      if ( FlagUseOfBulkCoefInNeutralCond ) then

        ! $BCfN)>r7o$G$N%P%k%/78?t$N@_Dj(B
        ! Set bulk coefficient in neutral condition
        !

        xy_SurfVelBulkCoef  = xy_SurfBulkCoefInNeutCond
        xy_SurfTempBulkCoef = xy_SurfBulkCoefInNeutCond

        xy_SurfQmixBulkCoef = xy_SurfTempBulkCoef

      else

        do j = jmin, jmax
          do i = imin, imax

            if ( xy_SurfBulkRiNum(i,j) > 0.0_DP ) then 

              xy_SurfVelBulkCoef(i,j) =                                        &
                &   xy_SurfBulkCoefInNeutCond(i,j)                             &
                &   / (   1.0_DP                                               &
                &       + 4.7_DP * xy_SurfBulkRiNum(i,j)                      &
                &     )**2

              xy_SurfTempBulkCoef(i,j) = xy_SurfVelBulkCoef(i,j) 

              xy_SurfQmixBulkCoef(i,j) = xy_SurfTempBulkCoef(i,j)

            else

              xy_SurfVelBulkCoef(i,j) =                                       &
                &   xy_SurfBulkCoefInNeutCond(i,j)                            &
                &   * (   1.0_DP                                              &
                &       - 9.4_DP * xy_SurfBulkRiNum(i,j)                     &
                &           / (   1.0_DP                                      &
                &             + 69.56_DP * xy_SurfBulkCoefInNeutCond(i,j)      &
                &               * sqrt( - z_dz(1)  &
                &                         / xy_SurfRoughLength(i,j)           &
                &                         * xy_SurfBulkRiNum(i,j)             &
                &                     )                                       &
                &           )                                                 &
                &     )

              xy_SurfTempBulkCoef(i,j) =                                      &
                &   xy_SurfBulkCoefInNeutCond(i,j)                            &
                &   * (   1.0_DP                                              &
                &       - 9.4_DP * xy_SurfBulkRiNum(i,j)                     &
                &         / (   1.0_DP                                        &
                &             + 49.82_DP * xy_SurfBulkCoefInNeutCond(i,j)      &
                &               * sqrt( - z_dz(1) &
                &                         / xy_SurfRoughLength(i,j)           &
                &                         * xy_SurfBulkRiNum(i,j)             &
                &                     )                                       &
                &           )                                                 &
                &     )

              xy_SurfQmixBulkCoef(i,j) = xy_SurfTempBulkCoef(i,j)
    
           end if
         end do
        end do
      end if
    end if

    ! $B:GBg(B/$B:G>.(B $BH=Dj(B
    ! Measure maximum/minimum
    !
    do i = imin, imax
      do j = jmin, jmax

        xy_SurfVelBulkCoef(i,j)  = &
          & max( min( xy_SurfVelBulkCoef(i,j), VelBulkCoefMax ), &
          &      VelBulkCoefMin )

        xy_SurfTempBulkCoef(i,j) = &
          & max( min( xy_SurfTempBulkCoef(i,j), TempBulkCoefMax ), &
          &      TempBulkCoefMin )

        xy_SurfQmixBulkCoef(i,j) = &
          & max( min( xy_SurfQmixBulkCoef(i,j), QmixBulkCoefMax ), &
          &      QmixBulkCoefMin )

      end do
    end do

  end subroutine BulkCoef
!!!------------------------------------------------------------------------!!!  
end module Surfaceflux_bulk_l1979
