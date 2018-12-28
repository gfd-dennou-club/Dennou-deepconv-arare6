!= Module Turbulence_kw1978
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA, Masatsugu
! Version::   $Id: turbulence_kw1978.f90,v 1.26 2014/06/07 17:34:27 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Turbulence_kw1978_v1
  !
  ! Klemp and Wilhelmson (1978) $B$NMpN.2aDx(B
  ! $BMpN.%(%M%k%.!<$N;~4VH/E8J}Dx<0$r2r$/$3$H$GMpN.3H;678?t$r7h$a$k(B
  !

  !$B%b%8%e!<%kFI$_9~$_(B 
  use dc_types, only: DP, STRING

  !$B0EL[$N7?@k8@6X;_(B
  implicit none

  !$B4X?t$r(B public $B$K@_Dj(B
  public Turbulence_KW1978_Init
  public Turbulence_KW1978_Forcing

  !$BJQ?tDj5A(B
  real(DP), save, private :: Cm     = 2.0d-1          !$BMpN.%(%M%k%.!<?GCG<0$N78?t(B 
  real(DP), save, private :: MixLen = 0.0d0           !$BJ?6Q:.9g5wN%(B
  real(DP), save, public  :: KmMax  = 0.0d0           !$BMpN.3H;678?t$N:GBgCM(B

  real(DP), save, private :: FactorDExnerDtTurb = 1.0d0 !$BMpN.3H;69`$r9MN8$9$k$+$N%9%$%C%A(B
  logical,  save, private :: FlagArare4  = .true.       !$B4pK\>l$NL)EY$r9MN8$9$k$+$N%9%$%C%A(B

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:) 
                                                  !$B4pK\>l$N1tD>HyJ,(B
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:)
                                                  !$B4pK\>l$N1tD>HyJ,(B

  character(*), parameter:: module_name = 'turbulence_kw1978'

contains

!!!------------------------------------------------------------------------!!!
  subroutine turbulence_kw1978_init
    !
    ! Turbulence $B%b%8%e!<%k$N=i4|2=%k!<%A%s(B
    ! 

    !$B%b%8%e!<%kFI$_9~$_(B 
    !    
    use gridset,       only : imin, imax, jmin, jmax, kmin, kmax, ncmax
    use axesset,       only : dx,            &! x $BJ}8~$N3J;RE@4V3V(B
      &                       dy,            &! y $BJ}8~$N3J;RE@4V3V(B
      &                       dz              ! z $BJ}8~$N3J;RE@4V3V(B
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use timeset,       only : DelTimeLong
    use gridset,       only : FlagCalc3D 
    use namelist_util, only : namelist_filename
    use differentiate_center2, &
      &                only : xyr_dz_xyz
    use basicset,      only : xyz_PTempBZ,   &! $B4pK\>l$N290L(B
      &                       xyzf_QMixBZ

    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    !$B:n6HJQ?t(B
    !
    integer :: unit, f

    !-------------------------------------------------------------------
    ! NAMELIST $B$+$i>pJs$r<hF@(B
    !
    NAMELIST /turbulence_kw1978_nml/     &
      & Cm, KmMax,                       &
      & FactorDExnerDtTurb, FlagArare4

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=turbulence_kw1978_nml)
    close(unit)

    !-------------------------------------------------------------------
    ! $B:.9g5wN%(B
    ! 2 $B<!857W;;$N>l9g$K$O(B DelY $B$K0MB8$7$J$$$h$&$K$9$k$?$a$K(B if $BJ8$rMxMQ(B.
    ! 
    if ( FlagCalc3D ) then 
      MixLen = ( dx * dy * dz ) ** (1.0d0 / 3.0d0)
    else
      MixLen = sqrt( dx * dz )
    end if
    
    !-------------------------------------------------------------------
    ! KmMax $B$,@_Dj$5$l$F$$$J$$>l9g(B. 
    ! $B0BDj@-2r@O$G$O(B, dt / l**2 < 0.5 $B$rK~$?$9I,MW$,$"$k(B. $B$3$3$G$O(B 0.1 $B$K$7$F$*$$$?(B. 
    !
    if (KmMax == 0.0d0) then 
       KmMax = 0.1 * (MixLen ** 2.0d0) / (DelTimeLong * 2.0d0) !LeapFrog
    end if

    !-------------------------------------------------------------------
    ! tendency $B$N=PNO(B
    !
    call turbulence_kw1978_output

    !-------------------------------------------------------------------
    ! $BG[Ns$NMQ0U(B
    !
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

    !-------------------------------------------------------------------
    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "Cm = %f", d=(/Cm/))
    call MessageNotify( "M", &
      & module_name, "KmMax = %f", d=(/KmMax/))
    call MessageNotify( "M", &
      & module_name, "MixLen = %f", d=(/MixLen/))
    call MessageNotify( "M", &
      & module_name, "FlagArare4 = %b", l=(/ FlagArare4 /))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtTurb= %f", d=(/ FactorDExnerDtTurb /))

  end subroutine turbulence_kw1978_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine turbulence_KW1978_forcing(       &
    & pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,   &
    & xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl,  &
    & xyz_KmBl,    xyz_KhBl,    xyz_CDensBl,  &
    & pyz_DVelXDt, xqz_DVelYDt, xyr_DVelZDt,  &
    & xyz_DPTempDt,xyz_DExnerDt,xyzf_DQMixDt, &
    & xyz_DKmDt,   xyz_DCDensDt               &
    & )

    !$B%b%8%e!<%kFI$_9~$_(B 
    !
    use dc_types, only: DP, STRING    
    use gtool_historyauto, only: HistoryAutoPut
    use composition, only: SpcWetSymbol
    use timeset, only:  TimeN
    use gridset, only:  imin,           &! x $BJ}8~$NG[Ns$N2<8B(B
      &                 imax,           &! x $BJ}8~$NG[Ns$N>e8B(B
      &                 jmin,           &! y $BJ}8~$NG[Ns$N2<8B(B
      &                 jmax,           &! y $BJ}8~$NG[Ns$N>e8B(B
      &                 kmin,           &! z $BJ}8~$NG[Ns$N2<8B(B
      &                 kmax,           &! z $BJ}8~$NG[Ns$N>e8B(B
      &                 nx,ny,nz,ncmax
    use basicset, only: xyz_PTempBZ,     &!$B4pK\>l$N290L(B
      &                 xyzf_QMixBZ,        &!$B4pK\>l$N:.9gHf(B 
      &                 xyz_ExnerBZ,        &!$B4pK\>l$N%(%/%9%J!<4X?t(B
      &                 xyz_DensBZ,         &!$B4pK\>l$NL)EY(B
      &                 pyz_DensBZ,         &!$B4pK\>l$NL)EY(B
      &                 xqz_DensBZ,         &!$B4pK\>l$NL)EY(B
      &                 xyr_DensBZ           !$B4pK\>l$NL)EY(B
    use constants,only: Grav,          &
      &                 MolWtDry,      & 
      &                 CpDry
    use average,  only: xyz_pyz, xyr_pyr, xqz_pqz, &
      &                 pyz_xyz, pyr_xyr, pqz_xqz, &
      &                 xyz_xqz, pyz_pqz, xyr_xqr, &
      &                 xqz_xyz, pqz_pyz, xqr_xyr, &
      &                 xyz_xyr, pyz_pyr, xqz_xqr, &
      &                 xyr_xyz, pyr_pyz, xqr_xqz, &
      &                 pqz_xyz, pyr_xyz, xqr_xyz, &
      &                 xyz_pqz, xyz_pyr, xyz_xqr
    use differentiate_center2,  &
      &           only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                 pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
      &                 pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
      &                 pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
      &                 pyr_dz_pyz, pyr_dx_xyr, xyr_dx_pyr, &
      &                 xqr_dz_xqz, xqr_dy_xyr, xyr_dy_xqr, &
      &                 pqz_dy_pyz, pqz_dx_xqz, xqz_dx_pqz, &
      &                 xqr_dy_xyr, xqr_dz_xqz, xqz_dz_xqr  
    use DExnerDt, only: xyz_DExnerDt_xyz

    ! $B0EL[$N7?@k8@6X;_(B
    ! 
    implicit none

    ! $BF~=PNOJQ?t(B
    !
    real(DP),intent(in) :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B?eJ?B.EY(B
    real(DP),intent(in) :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B?eJ?B.EY(B
    real(DP),intent(in) :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B1tD>B.EY(B
    real(DP),intent(in) :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B290L(B
    real(DP),intent(in) :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BL5<!8505NO(B
    real(DP),intent(in) :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                    !$B6E=L@.J,$N:.9gHf(B
    real(DP),intent(in) :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BMpN.3H;678?t(B
    real(DP),intent(in) :: xyz_KhBl(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BMpN.3H;678?t(B
    real(DP),intent(in) :: xyz_CDensBl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP),intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP),intent(inout):: xyr_DVelZDt(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP),intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP),intent(inout):: xyz_DKmDt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP),intent(inout) :: xyz_DCDensDt(imin:imax,jmin:jmax,kmin:kmax)

    ! $B:n6HJQ?t(B
    !
    real(DP)            :: xyz_Buoy(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B12G4@-78?t$N(B
    real(DP)            :: xyz_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B12G4@-78?t$N(B
    real(DP)            :: xyz_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B12G4@-78?t$N(B
    real(DP)            :: xyz_Shear(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B12G4@-78?t$N(B
    real(DP)            :: xyz_Diff(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B12G4@-78?t$N(B
    real(DP)            :: xyz_Disp(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BMpN.%(%M%k%.!<$N>C;6(B
    real(DP)            :: xyz_DispPI(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BMpN.%(%M%k%.!<$N>C;6(B
    real(DP)            :: xyz_DispHeat(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$BMpN.%(%M%k%.!<$N>C;6(B
    real(DP)            :: xyz_Turb(imin:imax,jmin:jmax,kmin:kmax)
                                                    !
    real(DP)            :: xyzf_Turb(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP)            :: pyz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xqz_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyr_Turb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: pyz_DVelXDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP)            :: xqz_DVelYDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP)            :: xyr_DVelZDt0(imin:imax,jmin:jmax,kmin:kmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP)            :: xyz_DPTempDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyz_DExnerDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyzf_DQMixDt0(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP)            :: xyz_DKmDt0(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyz_DCDensDt0(imin:imax,jmin:jmax,kmin:kmax)

    integer             :: f


    !----------------------------------
    ! $B=i4|2=(B
    !
    pyz_DVelXDt0 = pyz_DVelXDt
    xqz_DVelYDt0 = xqz_DVelYDt
    xyr_DVelZDt0 = xyr_DVelZDt
    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt
    xyzf_DQMixDt0 = xyzf_DQMixDt
    xyz_DKmDt0    = xyz_DKmDt
    xyz_DCDensDt0 = xyz_DCDensDt
    
    !----------------------------------
    ! $B3H;678?t$N;~4VH/E8(B ($B%(%M%k%.!<J}Dx<0$r(B Km $B$N<0$KJQ7A$7$?$b$N(B)
    !

    ! Buoyancy term
    !
    xyz_Buoy = xyz_BuoyMoistKm(xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl) 

    xyz_BuoyT =                                                  &
      &  - 3.0d0 * Grav * ( Cm * Cm  * MixLen * MixLen )         &
      &    * xyz_xyr(                                            &
      &        xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz        &
      &      ) / ( 2.0d0 * xyz_PTempBZ )
    
    xyz_BuoyM = xyz_Buoy - xyz_BuoyT

    xyz_Shear = &
      &   ( ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )                &
      & * (                                                        &
      &      ( xyz_dx_pyz( pyz_VelXBl ) ) ** 2.0d0                 &
      &    + ( xyz_dy_xqz( xqz_VelYBl ) ) ** 2.0d0                 &
      &    + ( xyz_dz_xyr( xyr_VelZBl ) ) ** 2.0d0                 &
      &    + 5.0d-1                                                &
      &      * (                                                   &
      &          (                                                 &
      &              xyz_pyr( pyr_dz_pyz( pyz_VelXBl ) )           &
      &            + xyz_pyr( pyr_dx_xyr( xyr_VelZBl ) )           &
      &           ) ** 2.0d0                                       &
      &        + (                                                 &
      &              xyz_xqr( xqr_dy_xyr( xyr_VelZBl ) )           &
      &            + xyz_xqr( xqr_dz_xqz( xqz_VelYBl ) )           &
      &           ) ** 2.0d0                                       &
      &        + (                                                 &
      &              xyz_pqz( pqz_dx_xqz( xqz_VelYBl ) )           &
      &            + xyz_pqz( pqz_dy_pyz( pyz_VelXBl ) )           &
      &           ) ** 2.0d0                                       &
      &        )                                                   &
      &   )                                                        &
      & - xyz_KmBl * (  xyz_dx_pyz( pyz_VelXBl )                   &
      &               + xyz_dy_xqz( xqz_VelYBl )                   &
      &               + xyz_dz_xyr( xyr_VelZBl ) ) / 3.0d0

    xyz_Diff =                                           &
      & + (                                              &
      &    + xyz_dx_pyz(pyz_dx_xyz(xyz_KmBl ** 2.0d0))   &
      &    + xyz_dy_xqz(xqz_dy_xyz(xyz_KmBl ** 2.0d0))   &
      &    + xyz_dz_xyr(xyr_dz_xyz(xyz_KmBl ** 2.0d0))   &
      &   ) * 5.0d-1                                     &
      & + (                                              &
      &    + (xyz_pyz(pyz_dx_xyz(xyz_KmBl))) ** 2.0d0    &
      &    + (xyz_xqz(xqz_dy_xyz(xyz_KmBl))) ** 2.0d0    &
      &    + (xyz_xyr(xyr_dz_xyz(xyz_KmBl))) ** 2.0d0    &
      &   )

    ! t - \Delta t $B$GI>2A(B
    !
    xyz_Disp = - (xyz_KmBl ** 2.0d0) * 5.0d-1 / (MixLen ** 2.0d0)

    ! tendency
    !
    xyz_DKmDt = ( xyz_Buoy + xyz_Shear + xyz_Diff + xyz_Disp ) + xyz_DKmDt0 

    call HistoryAutoPut(TimeN, 'DKmDtBuoy',  xyz_Buoy(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtBuoyT', xyz_BuoyT(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtBuoyM', xyz_BuoyM(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtShear', xyz_Shear(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DKmDtTurb',  xyz_Diff(1:nx,1:ny,1:nz)) 
    call HistoryAutoPut(TimeN, 'DKmDtDisp',  xyz_Disp(1:nx,1:ny,1:nz))    
    
    !--------------------------------
    ! $B1@L)EY$N(B tendency
    !
    if (FlagArare4) then 
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyz_CDensBl ) )
    else
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) )  &
      & + xyz_dz_xyr(                                                    &
      &     xyr_xyz( xyz_DensBZ * xyz_KhBl ) * xyr_dz_xyz( xyz_CDensBl ) &
      &   ) / xyz_DensBZ
    end if

    xyz_DCDensDt = xyz_DCDensDt0 + xyz_Turb

    call HistoryAutoPut(TimeN, 'DCDensDtTurb', xyz_Turb(1:nx,1:ny,1:nz))

    !--------------------------------
    ! $B290L$N(B tendency
    !
    if (FlagArare4) then 
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_DPTempBZDz )
    else
    xyz_Turb =                                                           &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
      & + xyz_dz_xyr(                                                    &
      &     xyr_xyz( xyz_DensBZ * xyz_KhBl )                             &
      &     * ( xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz )             &
      &   )                                                              &
      &   / xyz_DensBZ
    end if
    
    xyz_DispHeat = (xyz_KmBl ** 3.0d0)                                   &
      & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0))

    xyz_DPTempDt = ( xyz_Turb + xyz_DispHeat ) + xyz_DPTempDt0 

    call HistoryAutoPut(TimeN, 'DPTempDtDisp', xyz_DispHeat(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtTurb', xyz_Turb(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! $B:.9gHf$N(B tendency
    !
    do f = 1, ncmax    
      if (FlagArare4) then 
      xyzf_Turb(:,:,:,f) =                                                         &
        & + xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyrf_DQMixBZDz(:,:,:,f) )
      else
      xyzf_Turb(:,:,:,f) =                                                         &
        &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
        & + xyz_dz_xyr(                                                            &
        &     xyr_xyz( xyz_DensBZ * xyz_KhBl )                                     &
        &     * ( xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) + xyrf_DQMixBZDz(:,:,:,f) )   &
        &   ) / xyz_DensBZ 
      end if 
      
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtTurb', xyzf_Turb(1:nx,1:ny,1:nz,f))
    end do

    xyzf_DQMixDt = xyzf_DQMixDt0 + xyzf_Turb
    
    !--------------------------------
    ! VelX $B$N(B tendency
    !
    if (FlagArare4) then 
    pyz_Turb =                                                     &
      &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )&
      & + pyz_dy_pqz(                                              &
      &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )       &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )       &
      &   )                                                        &
      & + pyz_dz_pyr(                                              &
      &       pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )       &
      &     + pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )       &
      &   )                                                        &
      & - 2.0d0 * pyz_dx_xyz( ( xyz_KmBl ** 2.0d0 ) )              &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    pyz_Turb =                                                             &
      &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )        &
      & + pyz_dy_pqz(                                                      &
      &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )               &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )               &
      &   )                                                                &
      & + pyz_dz_pyr(                                                      &
      &       pyr_xyz( xyz_DensBZ * xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )  &
      &     + pyr_xyz( xyz_DensBZ * xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )  &
      &   ) / pyz_DensBZ                                                   &
      & - 2.0d0 * pyz_dx_xyz( ( xyz_KmBl ** 2.0d0 ) )                      &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    end if

    pyz_DVelXDt = pyz_DVelXDt0 + pyz_Turb

    call HistoryAutoPut(TimeN, 'DVelXDtTurb', pyz_Turb(1:nx, 1:ny, 1:nz))


    !--------------------------------
    ! VelY $B$N(B tendency
    !
    if (FlagArare4) then 
    xqz_Turb =                                                      &
      &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) ) &
      & + xqz_dx_pqz(                                               &
      &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )        &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )        &
      &   )                                                         &
      & + xqz_dz_xqr(                                               &
      &       xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )        &
      &     + xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )        &
      &   )                                                         &
      & - 2.0d0 * xqz_dy_xyz( ( xyz_KmBl ** 2.0d0 ) )               &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    xqz_Turb =                                                            &
      &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) )       &
      & + xqz_dx_pqz(                                                     &
      &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )              &
      &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )              &
      &   )                                                               &
      & + xqz_dz_xqr(                                                     &
      &       xqr_xyz( xyz_DensBZ * xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl ) &
      &     + xqr_xyz( xyz_DensBZ * xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl ) &
      &   ) / xqz_DensBZ                                                  &
      & - 2.0d0 * xqz_dy_xyz( ( xyz_KmBl ** 2.0d0 ) )                     &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    end if 
    
    xqz_DVelYDt = xqz_DVelYDt0 + xqz_Turb

    call HistoryAutoPut(TimeN, 'DVelYDtTurb', xqz_Turb(1:nx, 1:ny, 1:nz))


    !--------------------------------
    ! VelZ $B$N(B tendency
    !
    if (FlagArare4) then 
    xyr_Turb =                                                       &
      & + 2.0d0 * xyr_dz_xyz( xyz_KmBl * xyz_dz_xyr( xyr_VelZBl ) )  &
      & + xyr_dx_pyr(                                                &
      &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )          &
      &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )          &
      &   )                                                          &
      & + xyr_dy_xqr(                                                &
      &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )          &
      &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )          &
      &   )                                                          & 
      & - 2.0d0 * xyr_dz_xyz(  xyz_KmBl ** 2.0d0 )                   &
      &   / ( 3.0d0 * ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )
    else
    xyr_Turb =                                                       &
      & + 2.0d0                                                      &
      &   * xyr_dz_xyz(                                              &
      &       xyz_DensBZ * xyz_KmBl * xyz_dz_xyr( xyr_VelZBl )       &
      &     ) / xyr_DensBZ                                           &
      & + xyr_dx_pyr(                                                &
      &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )          &
      &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )          &
      &   )                                                          &
      & + xyr_dy_xqr(                                                &
      &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )          &
      &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )          &
      &   )                                                          & 
      & - 2.0d0 * xyr_dz_xyz( xyz_DensBZ * xyz_KmBl * xyz_KmBl )     &
      &   / ( 3.0d0 * ( Cm * Cm * MixLen * MixLen ) )                &
      &   / xyr_DensBZ 
    end if

    xyr_DVelZDt = xyr_DVelZDt0 + xyr_Turb

    call HistoryAutoPut(TimeN, 'DVelZDtTurb', xyr_Turb(1:nx, 1:ny, 1:nz))

    !--------------------
    ! Exner function
    !
    xyz_DispPI = xyz_DExnerDt_xyz( xyz_DispHeat ) * FactorDExnerDtTurb
    xyz_DExnerDt = xyz_DExnerDt0 + xyz_DispPI

    call HistoryAutoPut(TimeN, 'DExnerDtDisp', xyz_DispPI(1:nx, 1:ny, 1:nz))


  contains

  function xyz_BuoyMoistKm(xyz_PTemp, xyz_Exner, xyzf_QMix)
    !
    ! $BIbNO9`$N7W;;(B
    !

    !$B%b%8%e!<%k8F$S=P$7(B
    !
    use dc_types, only: DP, STRING    
    use composition, only: MolWtWet,        &
      &                 SpcWetID,           &
      &                 CondNum,            &!$B6E7k2aDx$N?t(B
      &                 IdxCG,              &!$B6E7k2aDx(B($B>x5$(B)$B$NG[NsE:$(;z(B
      &                 IdxCC,              &!$B6E7k2aDx(B($B1@(B)$B$NG[NsE:$(;z(B
      &                 GasNum,             &!$B5$BN$N?t(B
      &                 IdxG                 !$B5$BN$NG[NsE:$(;z(B
    use ChemCalc, only: xyz_LatentHeat       !$B@xG.(B
    !    &              ReactHeatNH4SH       !NH4SH $B$NH?1~G.(B
    use gtool_historyauto,                  &
      &           only: HistoryAutoPut
    use basicset, only: xyz_PTempBZ,        &!$B4pK\>l$N290L(B
      &                 xyz_QMixBZ,         &!$B4pK\>l$N:.9gHf(B
      &                 xyz_QMixBZPerMolWt, &!$B4pK\>l$N:.9gHf(B
      &                 xyz_EffMolWtBZ,     &!$B4pK\>l$N:.9gHf(B
      &                 xyz_ExnerBZ          !$B4pK\>l$N%(%/%9%J!<4X?t(B
    use average,  only: xyz_xyr
    use differentiate_center2,              &
      &           only: xyr_dz_xyz

    !$B0EL[$N7?@k8@6X;_(B
    implicit none
    
    !$BJQ?tDj5A(B
    real(DP), intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                                               !$B290L(B
    real(DP), intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                                               !$BL5<!8505NO(B
    real(DP), intent(in) :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !$B6E=L@.J,$N:.9gHf(B
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !$B6E=L@.J,$N:.9gHf(B
    real(DP) :: xyzf_QMixAll2(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !$B6E=L@.J,$N:.9gHf(B
    real(DP) :: xyz_BuoyMoistKm(imin:imax,jmin:jmax,kmin:kmax)
                                               !
    real(DP) :: xyzf_LatentHeat(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !$B@xG.(B
    real(DP) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                               !$B29EY(B
    real(DP) :: xyz_EffHeat(imin:imax,jmin:jmax,kmin:kmax)
                                               !
    real(DP) :: xyz_EffPTemp(imin:imax,jmin:jmax,kmin:kmax)    
                                               !$BIbNO$KBP$9$k29EY:9$N4sM?(B
    real(DP) :: xyz_EffMolWt(imin:imax,jmin:jmax,kmin:kmax)    
                                               !$BIbNO$KBP$9$kJ,;RNL:9$N4sM?(B
    real(DP) :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !$B:.9gHf(B/$BJ,;RNL(B
    integer              :: s
    
    
    !$B29EY(B, $B05NO(B, $B:.9gHf$NA4NL$r5a$a$k(B
    !$B>qMp@.J,$HJ?6Q@.J,$NB-$7;;(B
    xyz_TempAll     = ( xyz_PTemp + xyz_PTempBZ ) * ( xyz_Exner + xyz_ExnerBZ )
    xyzf_QMixAll    = xyzf_QMixBZ + xyzf_QMix
    xyzf_LatentHeat = 0.0d0
    
    !$B:n6HG[Ns$N=i4|2=(B. $B5$BN$N$_MxMQ(B
    do s = 1, GasNum
      xyzf_QMixPerMolWt(:,:,:,s) = xyzf_QMix(:,:,:,IdxG(s)) / MolWtWet(IdxG(s))
    end do
    
    !$B29EY$N8z2L(B
    xyz_EffPTemp = xyz_PTemp / xyz_PTempBZ 
    
    !$BJ,;RNL8z2L(B + $B0z$-$E$j$N8z2L(B
    xyz_EffMolWt =                                      &
      & + sum(xyzf_QMixPerMolWt, 4)                     &
      &    / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt )  &
      & - sum(xyzf_QMix, 4) / ( 1.0d0 + xyz_QMixBZ )
    
    !$B>x5$$,>xH/$9$k>l9g$N@xG.$r7W;;(B
    !  $BJ,;RNL$NItJ,$O$$$D$G$b8z$/$,@xG.$OK0OB$7$F$$$J$$$H8z$+$J$$$N$G(B, 
    !  $B1@$N:.9gHf$,%<%m$N;~$K$O(B, $B@xG.$N4sM?$O%<%m$H$J$k$h$&$KD4@a$7$F$$$k(B
    !
    xyzf_QMixAll2 = xyzf_QMixAll
!    xzya_QMixAll2(:,:,:,IdxNH3) = &
!      &  xyzf_QMixAll(:,:,:,IdxNH3) - xyzf_QMixAll(:,:,:,IdxH2S) 

    do s = 1, CondNum
      xyzf_LatentHeat(:,:,:,s) =                                                   &
        & xyz_LatentHeat( SpcWetID(IdxCC(s)), xyz_TempAll )                        &
        &  * xyzf_QMixAll2(:,:,:,IdxCG(s))                                         &
        &  * ( 5.0d-1 + sign( 5.0d-1, (xyzf_QMixAll2(:,:,:,IdxCC(s)) - 1.0d-4) ) )
    end do

    xyz_EffHeat = ( sum( xyzf_LatentHeat, 4 ) * xyz_EffMolWtBZ &
!      &            + ReactHeatNH4SH * &
!      &            (xyzf_QMixAll(:,:,:,IdxH2S) + xyzf_QMixAll(:,:,:,IdxH2S)) &
      &          ) / ( CpDry * xyz_ExnerBZ ) 

    
    !$BMpN.3H;678?t$N;~4VH/E8<0$NIbNO9`$r7h$a$k(B
    xyz_BuoyMoistKm = &
      &  - 3.0d0 * Grav * (( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ))   &
      &    * xyz_xyr(                                               &
      &        xyr_dz_xyz(                                          &
      &           xyz_EffHeat                                       &
      &         + xyz_PTempBZ / xyz_EffMolWtBZ                      &
      &           * ( 1.0d0 + xyz_EffPTemp + xyz_EffMolWt )         &
      &         )                                                   &
      &      )                                                      &
      &    / ( 2.0d0 * xyz_PTempBZ / xyz_EffMolWtBZ)                   

  end function xyz_BuoyMoistKm

  end subroutine turbulence_KW1978_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine turbulence_kw1978_output
    !
    ! tendency $B$N=PNO$N$?$a$N@_Dj(B
    !

    ! $B%b%8%e!<%k8F$S=P$7(B
    !
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : ncmax
    use composition,       only : SpcWetSymbol

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $B:n6HJQ?t(B
    !
    integer :: l

    !-----------------------------------------------------
    ! tendency $B=PNO$N$?$a$N@_Dj(B
    !
    call HistoryAutoAddVariable(  &
      & varname='DVelXDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (x)', &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelYDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (y)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (z)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of exner function', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DCDensDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of cloud density', &
      & units='kg.m-3.s-1',    &
      & xtype='double')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtTurb', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Turbulence term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
    end do

    call HistoryAutoAddVariable(  &
      & varname='DKmDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Diffusion term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoy',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoyT',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km (temperature)', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtBuoyM',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy term of Km (composition)', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtShear',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Shear term of Km', &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtDisp',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Dissipation term of Km', &
      & units='s-1',    &
      & xtype='double')

  end subroutine turbulence_kw1978_output

end module Turbulence_kw1978_v1
