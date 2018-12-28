!= Module Turbulence_kw1978
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA, Masatsugu
! Version::   $Id: turbulence_kw1978_v2.f90,v 1.3 2014/07/08 07:49:01 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Turbulence_kw1978_v2
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

  real(DP), save, private :: FactorDExnerDtTurb = 0.0d0 !$BMpN.3H;69`$r9MN8$9$k$+$N%9%$%C%A(B
  logical,  save, private :: FlagArare4  = .false.       !$B4pK\>l$NL)EY$r9MN8$9$k$+$N%9%$%C%A(B

  integer, save, private      :: IDTurbulence     = 0 !$BC;$$;~4V%9%F%C%W$N7W;;J}K!(B
  integer, parameter, private :: IDTurbulence_std = 1 ! 1: $BHyJ,J?6Q%b%8%e!<%k$NMxMQ(B.
  integer, parameter, private :: IDTurbulence_2D  = 2 ! 2: 2D $BHG7W;;%k!<%A%s(B
  integer, parameter, private :: IDTurbulence_3D  = 3 ! 3: 3D $BHG7W;;%k!<%A%s(B

  real(DP), allocatable, save, private :: xyr_DPTempBZDz(:,:,:)   !$B290L$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: xyr_DExnerBZDz(:,:,:)   !$B%(%/%9%J!<4X?t$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: xyz_DDensBZDz(:,:,:)    !$BL)EY$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: pyz_DDensBZDz(:,:,:)    !$BL)EY$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: xqz_DDensBZDz(:,:,:)    !$BL)EY$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: xyr_DDensBZDz(:,:,:)    !$BL)EY$N1tD>HyJ,(B
  real(DP), allocatable, save, private :: xyrf_DQMixBZDz(:,:,:,:) !$B:.9gHf$N1tD>HyJ,(B

  real(DP), save, private :: MixLen_MixLen = 0.0d0 
  real(DP), save, private :: Cm_Cm_MixLen_MixLen = 0.0d0 

  character(STRING), parameter:: module_name = 'turbulence_kw1978_v2'

contains

!!!------------------------------------------------------------------------!!!
  subroutine turbulence_kw1978_init
    !
    ! Turbulence $B%b%8%e!<%k$N=i4|2=%k!<%A%s(B
    ! 

    !$B%b%8%e!<%kFI$_9~$_(B 
    !
    use dc_types,      only : STRING
    use axesset,       only : dx,            &! x $BJ}8~$N3J;RE@4V3V(B
      &                       dy,            &! y $BJ}8~$N3J;RE@4V3V(B
      &                       dz              ! z $BJ}8~$N3J;RE@4V3V(B
    use gridset,       only : imin,           &! x $BJ}8~$NG[Ns$N2<8B(B
      &                       imax,           &! x $BJ}8~$NG[Ns$N>e8B(B
      &                       jmin,           &! y $BJ}8~$NG[Ns$N2<8B(B
      &                       jmax,           &! y $BJ}8~$NG[Ns$N>e8B(B
      &                       kmin,           &! z $BJ}8~$NG[Ns$N2<8B(B
      &                       kmax,           &! z $BJ}8~$NG[Ns$N>e8B(B
      &                       ncmax
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use timeset,       only : DelTimeLong
    use gridset,       only : FlagCalc3D 
    use namelist_util, only : namelist_filename
    use average,       only : xyz_xyr 
    use differentiate_center2,           &
      &                only : xyr_dz_xyz
    use basicset,      only : &
      &                 xyz_ExnerBZ,     &!$B4pK\>l$N%(%/%9%J!<4X?t(B
      &                 xyz_PTempBZ,     &!$B4pK\>l$N290L(B
      &                 xyzf_QMixBZ,     &!$B4pK\>l$N:.9gHf(B
      &                 xyz_DensBZ        !$B4pK\>l$NL)EY(B
    
    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $B:n6HJQ?t(B
    !
    integer :: unit
    integer :: s
    character(STRING) :: FlagTurbulence = ""     !$BMpN.7W;;J}K!$NA*Br(B

    !-------------------------------------------------------------------
    ! NAMELIST $B$+$i>pJs$r<hF@(B
    !
    NAMELIST /turbulence_kw1978_nml/    &
      & Cm, KmMax,                      &
      & FactorDExnerDtTurb, FlagArare4, &
      & FlagTurbulence

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=turbulence_kw1978_nml)
    close(unit)

    !-------------------------------------------------------------------
    ! $B:.9g5wN%(B
    ! 2 $B<!857W;;$N>l9g$K$O(B DelY $B$K0MB8$7$J$$$h$&$K$9$k$?$a$K(B if $BJ8$rMxMQ(B.
    ! 
    if ( FlagCalc3D ) then 
      MixLen = (dx * dy * dz ) ** (1.0d0 / 3.0d0)
    else
      MixLen = sqrt( dx * dz ) 
    end if

    ! $B$^$H$a$F$*$/(B
    !
    MixLen_MixLen = MixLen * MixLen
    Cm_Cm_MixLen_MixLen = Cm * Cm * MixLen * MixLen
    
    !-------------------------------------------------------------------
    ! KmMax $B$,@_Dj$5$l$F$$$J$$>l9g(B. 
    ! $B0BDj@-2r@O$G$O(B, dt / l**2 < 0.5 $B$rK~$?$9I,MW$,$"$k(B. $B$3$3$G$O(B 0.1 $B$K$7$F$*$$$?(B. 
    !
    if (KmMax == 0.0d0) then 
       KmMax = 0.1 * (MixLen ** 2.0d0) / (DelTimeLong * 2.0d0) !LeapFrog
    end if

    !-------------------------------------------------------------------
    ! $B7W;;J}K!$NA*Br(B
    !   std: $B=>Mh$NHyJ,J?6Q1i;;$rMxMQ(B
    !   3d, 2d: $BHyJ,J?6Q1i;;$r<jF0$G%$%s%i%$%sE83+$7$?HG(B
    !
    if ( FlagTurbulence == "std" ) then 
      IDTurbulence = IDTurbulence_std

    else
      if ( FlagCalc3D ) then 
        IDTurbulence = IDTurbulence_3d
        if ( FlagArare4 ) &
          & call MessageNotify( "E", module_name, "FlagArare4 should be .false.")
      else
        IDTurbulence = IDTurbulence_2d
        if ( FlagArare4 ) &
          & call MessageNotify( "E", module_name, "FlagArare4 should be .false.")
      end if

    end if

    !-------------------------------------------------------------------
    ! $B4pK\>l$NHyJ,(B
    !
    allocate( xyr_DPTempBZDz(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DExnerBZDz(imin:imax,jmin:jmax,kmin:kmax) )
    allocate(                                         &
      & xyz_DDensBZDz(imin:imax,jmin:jmax,kmin:kmax), &
      & pyz_DDensBZDz(imin:imax,jmin:jmax,kmin:kmax), &
      & xqz_DDensBZDz(imin:imax,jmin:jmax,kmin:kmax), &
      & xyr_DDensBZDz(imin:imax,jmin:jmax,kmin:kmax)  &
      & )
    allocate( xyrf_DQMixBZDz(imin:imax,jmin:jmax,kmin:kmax,ncmax) )
    
    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ ) 
    xyr_DExnerBZDz = xyr_dz_xyz( xyz_ExnerBZ ) 
    do s = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,s) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,s) )
    end do

    xyr_DDensBZDz  = xyr_dz_xyz( xyz_DensBZ ) 
    xyz_DDensBZDz  = xyz_xyr( xyr_DDensBZDz )

    ! $B?eJ?0lMM$J$N$G(B, $BJ?6QA`:n$OI,MW$J$$(B
    !
    pyz_DDensBZDz  = xyz_DDensBZDz
    xqz_DDensBZDz  = xyz_DDensBZDz

    !-------------------------------------------------------------------
    ! tendency $B$N=PNO(B
    !
    call turbulence_kw1978_output
   
    !-------------------------------------------------------------------
    ! Output
    !
    call MessageNotify( "M", module_name, "Cm = %f", d=(/Cm/))
    call MessageNotify( "M", module_name, "KmMax = %f", d=(/KmMax/))
    call MessageNotify( "M", module_name, "MixLen= %f", d=(/MixLen/))
    call MessageNotify( "M", module_name, &
      &                 "FlagArare4 = %b", l=(/ FlagArare4 /))
    call MessageNotify( "M", module_name, &
      &                 "FactorDExnerDtTurb = %f", d=(/FactorDExnerDtTurb/))
    call MessageNotify( "M", module_name, &
      &                 "IDTurbulence = %d", i=(/IDTurbulence/))

  end subroutine turbulence_kw1978_init
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine turbulence_KW1978_forcing(        &
    & pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,    &
    & xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl,   &
    & xyz_KmBl,    xyz_KhBl,    xyz_CDensBl,   &
    & pyz_DVelXDt, xqz_DVelYDt,  xyr_DVelZDt,  &
    & xyz_DPTempDt,xyz_DExnerDt, xyzf_DQMixDt, &
    & xyz_DKmDt,   xyz_DCDensDt                &
    & )
    
    use dc_types,  only : DP
    use gtool_historyauto, only : &
      &                   HistoryAutoPut
    use timeset,   only : TimeN
    use gridset,   only : imin,           &! x $BJ}8~$NG[Ns$N2<8B(B
      &                   imax,           &! x $BJ}8~$NG[Ns$N>e8B(B
      &                   jmin,           &! y $BJ}8~$NG[Ns$N2<8B(B
      &                   jmax,           &! y $BJ}8~$NG[Ns$N>e8B(B
      &                   kmin,           &! z $BJ}8~$NG[Ns$N2<8B(B
      &                   kmax,           &! z $BJ}8~$NG[Ns$N>e8B(B
      &                   ncmax,          &
      &                   nx,ny,nz       
    use DExnerDt,  only : xyz_DExnerDt_xyz
    use composition,only: SpcWetSymbol

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
    real(DP)            :: xyz_TurbD(imin:imax,jmin:jmax,kmin:kmax)


    real(DP)            :: xyzf_Turb(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP)            :: xyzf_TurbD(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                    !$B%9%+%i!<NL$N?eJ?MpN.3H;6(B
    real(DP)            :: xyz_DCDensDtTurb(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: xyz_DCDensDtTurbD(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: pyz_Turb(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: pyz_TurbD(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xqz_Turb(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: xqz_TurbD(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)            :: xyr_Turb(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: xyr_TurbD(imin:imax,jmin:jmax,kmin:kmax)

    integer             :: f

    
    !----------------------------------
    ! $B3FM=JsJQ?t$KBP$7$FMpN.3H;6$K$h$k(B tendency $B$r7W;;$9$k(B. 
    ! $BFbIt%5%V%k!<%A%s$rMxMQ(B.
    !
    select case ( IDTurbulence ) 
    case ( IDTurbulence_std )

      call turbulence_kw1978_std

    case ( IDTurbulence_3d )

      call turbulence_kw1978_center2_3d

    case ( IDTurbulence_2d )

      call turbulence_kw1978_center2_2d

    end select

    ! $BL)EY$N1tD>JQ2=$r9MN8$9$k$3$H$K$h$k4sM?$rF@$k(B
    !
    call turbulence_kw1978_densityeffect

    !-----------------------------------------
    ! Km $B$N(B tendency
    !
    xyz_DKmDt = ( xyz_Buoy + xyz_Shear + xyz_Diff + xyz_Disp ) + xyz_DKmDt

    ! tendency $B$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DKmDtBuoy',  xyz_Buoy(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtBuoyT', xyz_BuoyT(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtBuoyM', xyz_BuoyM(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtShear', xyz_Shear(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtTurb',  xyz_Diff(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DKmDtDisp',  xyz_Disp(1:nx,1:ny,1:nz) )
   
    !--------------------------------
    ! $B1@L)EY$N(B tendency
    !
    xyz_DCDensDt = xyz_DCDensDt + xyz_DCDensDtTurb

    call HistoryAutoPut(TimeN, 'DCDensDtTurb', xyz_DCDensDtTurb(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DCDensDtTurbD', xyz_DCDensDtTurbD(1:nx,1:ny,1:nz))
    
    !--------------------------------
    ! $B290L$N(B tendency
    !
    xyz_DPTempDt = ( xyz_Turb + xyz_DispHeat ) + xyz_DPTempDt
    
    call HistoryAutoPut(TimeN, 'DPTempDtDisp', xyz_DispHeat(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtTurb', xyz_Turb(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtTurbD', xyz_TurbD(1:nx, 1:ny, 1:nz))

    !--------------------------------
    ! $B:.9gHf$N(B tendency
    !
    xyzf_DQMixDt = xyzf_DQMixDt + xyzf_Turb

    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtTurb', xyzf_Turb(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtTurbD', xyzf_TurbD(1:nx,1:ny,1:nz,f))
    end do

    !--------------------------------
    ! $B3FB.EY@.J,$N(B tendency
    !
    pyz_DVelXDt = pyz_DVelXDt + pyz_Turb
    xqz_DVelYDt = xqz_DVelYDt + xqz_Turb
    xyr_DVelZDt = xyr_DVelZDt + xyr_Turb

    call HistoryAutoPut(TimeN, 'DVelXDtTurb', pyz_Turb(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtTurb', xqz_Turb(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtTurb', xyr_Turb(1:nx, 1:ny, 1:nz))

    call HistoryAutoPut(TimeN, 'DVelXDtTurbD', pyz_TurbD(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtTurbD', xqz_TurbD(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtTurbD', xyr_TurbD(1:nx, 1:ny, 1:nz))

    !--------------------
    ! tendency of Exner function 
    !
    xyz_DispPI = xyz_DExnerDt_xyz( xyz_DispHeat ) * FactorDExnerDtTurb
    xyz_DExnerDt = xyz_DExnerDt + xyz_DispPI
    
    call HistoryAutoPut(TimeN, 'DExnerDtDisp', xyz_DispPI(1:nx, 1:ny, 1:nz))

  contains

    subroutine turbulence_kw1978_densityeffect
      
      use average,  only: xyz_pyz, xyz_xyr, pyz_xyz, pyz_pyr, &
        &                 xqz_xyz, xqz_xqr, xyr_xyz
      use differentiate_center2, &
        &           only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
        &                 pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
        &                 pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
        &                 pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
        &                 xyr_dx_pyr, xyr_dy_xqr, xqr_dz_xqz, &
        &                 xqz_dx_pqz, xqr_dy_xyr, xqz_dz_xqr  
      use basicset, only: xyz_DensBZ

      ! $B1@L)EY(B
      !
      xyz_DCDensDtTurbD =                                             &
        &   xyz_KhBl                                               &
        &   * xyz_xyr( xyr_DDensBZDz * xyr_dz_xyz( xyz_CDensBl ) ) &
        &   / xyz_DensBZ
      
      ! $B290L(B
      !
      xyz_TurbD =                                               &
        &   xyz_KhBl                                            &
        & * xyz_xyr(                                            &
        &    xyr_DDensBZDz                                      &
        &    * ( xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz )   &
        &   ) / xyz_DensBZ

      ! $B:.9gHf(B
      !
      do f = 1, ncmax    
        xyzf_TurbD(:,:,:,f) =                                                      &
          &   xyz_KhBl                                                             &
          & * xyz_xyr(                                                             &
          &     xyz_DDensBZDz                                                      &
          &     * ( xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) + xyrf_DQMixBZDz(:,:,:,f) ) &
          &   ) / xyz_DensBZ
      end do
      
      ! u
      !
      pyz_TurbD =                                                          &
        &   pyz_xyz( xyz_KmBl * xyz_xyr( xyr_DDensBZDz ) / xyz_DensBZ )    &
        & * pyz_pyr( pyr_dx_xyr( xyr_VelZBl ) + pyr_dz_pyz( pyz_VelXBl ) )
      
      ! v
      !
      xqz_TurbD =                                                          &
        &   xqz_xyz( xyz_KmBl * xyz_xyr( xyr_DDensBZDz ) / xyz_DensBZ )    &
        & * xqz_xqr( xqr_dy_xyr( xyr_VelZBl ) + xqr_dz_xqz( xqz_VelYBl ) )

      ! w
      !
      xyr_TurbD =                                                                               &
        & + 2.0d0 * xyr_DDensBZDz                                                               &
        &   * (                                                                                 &
        &       xyr_xyz( xyz_KmBl * xyz_dz_xyr( xyr_VelZBl ) / xyz_DensBZ )                     &
        &     - xyr_xyz( xyz_KmBl * xyz_KmBl / xyz_DensBZ ) / ( 3.0d0 * Cm_Cm_MixLen_MixLen )   &
        &   ) 

    end subroutine turbulence_kw1978_densityeffect

    
    subroutine turbulence_kw1978_std

      !$B%b%8%e!<%k$NFI$_9~$_(B
      !
      use average,   only : xyz_pyz, xyz_xqz, xyz_xyr,  &
        &                   xyz_pqz, xyz_pyr, xyz_xqr,  &
        &                   pqz_xyz, pyz_xyz, &
        &                   xqz_xyz, xqr_xyz, &
        &                   xyr_xyz, pyr_xyz
      use differentiate_center2, only: &
        &                   xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
        &                   pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
        &                   pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
        &                   pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
        &                   xyr_dx_pyr, xyr_dy_xqr, xqr_dz_xqz, &
        &                   xqz_dx_pqz, xqr_dy_xyr, xqz_dz_xqr  
      use basicset,  only : xyz_ExnerBZ,     &!$B4pK\>l$N%(%/%9%J!<4X?t(B
        &                   xyz_DensBZ,      &!$B4pK\>l$NL)EY(B
        &                   pyz_DensBZ,      &!$B4pK\>l$NL)EY(B
        &                   xqz_DensBZ,      &!$B4pK\>l$NL)EY(B
        &                   xyr_DensBZ        !$B4pK\>l$NL)EY(B
      use constants, only : CpDry
      
      ! $B0EL[$N7?@k8@6X;_(B
      ! 
      implicit none

      ! $B:n6HJQ?t(B
      ! 
      real(DP) :: xyz_KmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKhBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)

!      real(DP) :: pyr_KmBl(imin:imax, jmin:jmax, kmin:kmax)
!      real(DP) :: pyr_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)

      !---------------------------------------------------------
      ! $B:n6HJQ?t$NDj5A(B
      !
      xyz_KmBlKmBl        = xyz_KmBl   * xyz_KmBl
      xyz_DensBZKhBl      = xyz_DensBZ * xyz_KhBl
      xyz_DensBZKmBl      = xyz_DensBZ * xyz_KmBl
      xyz_DensBZKmBlKmBl  = xyz_DensBZ * xyz_KmBl * xyz_KmBl

      !---------------------------------------------------------
      ! Km $B$N;~4VH/E8J}Dx<0$N3F9`(B

      ! $B;60o9`(B
      !
      xyz_Disp = - xyz_KmBlKmBl * 5.0d-1 / MixLen_MixLen

      ! $BIbNO9`(B
      !
      call turbulence_kw1978_BuoyancyKm
      
      ! $B%7%"!<9`(B
      !
      xyz_Shear =                                                    &
        &   Cm_Cm_MixLen_MixLen                                      &
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
      
      ! $B3H;69`(B
      !
      xyz_Diff =                                                     &
        &   5.0d-1                                                   &
        &   * (                                                      &
        &         xyz_dx_pyz(pyz_dx_xyz(xyz_KmBlKmBl))               &
        &       + xyz_dy_xqz(xqz_dy_xyz(xyz_KmBlKmBl))               &
        &       + xyz_dz_xyr(xyr_dz_xyz(xyz_KmBlKmBl))               &
        &     )                                                      &
        & + (                                                        &
        &      (xyz_pyz(pyz_dx_xyz(xyz_KmBl))) ** 2.0d0              &
        &    + (xyz_xqz(xqz_dy_xyz(xyz_KmBl))) ** 2.0d0              &
        &    + (xyz_xyr(xyr_dz_xyz(xyz_KmBl))) ** 2.0d0              &
        &   )
      
      !--------------------------------
      ! $B1@L)EY$N(B tendency
      !
      if (FlagArare4) then 
      xyz_DCDensDtTurb =                                                     &
        &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) ) &
        & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyz_CDensBl ) )
      else
      xyz_DCDensDtTurb =                                                     &
        &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_CDensBl ) ) &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_CDensBl ) ) &
        & + xyz_dz_xyr(                                                   &
        &     xyr_xyz( xyz_DensBZKhBl ) * xyr_dz_xyz( xyz_CDensBl )       &
        &   ) / xyz_DensBZ
      end if
      
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
      !
      ! \Dinv{\rho} \DP{\rho Kh \DP{\theta}{z}}{z} $B$r$=$N$^$^<B9T$9$k$+(B, 
      ! \DP{Kh \DP{\theta}{z}}{z} + \Dinv{\rho} Kh \DP{\theta}{z} \DP{\rho}{z} 
      ! $B$H$$$&7A$KJQ7A$7$F$+$i<B9T$9$k$+$G(B, $BCM$,$:$l$k(B. $B?tCMHyJ,$J$N$G(B. 
      ! $B$=$N$^$^<B9T$9$k$N$,:G$b8m:9$,>.$5$$$@$m$&$H$$$&9M$($N$b$H(B, $BL)EY$r(B
      ! $B9MN8$9$k$3$H$G@8$8$k$*$D$j$N9`$rJL$K7W;;$9$k$3$H$K$7$?(B. 
      !
      xyz_Turb =                                                           &
        &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyz_PTempBl ) )  &
        & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyz_PTempBl ) )  &
        & + xyz_dz_xyr(                                                    &
        &     xyr_xyz( xyz_DensBZKhBl )                                    &
        &     * ( xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz )             &
        &   ) / xyz_DensBZ
      end if

      xyz_DispHeat = (xyz_KmBl ** 3.0d0) &
        & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0))
      
      !--------------------------------
      ! $B:.9gHf$N(B tendency
      !
      do f = 1, ncmax    
        if (FlagArare4) then 
        xyzf_Turb(:,:,:,f) =                                                         &
          &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
          & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
          & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
          & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyrf_DQMixBZDz(:,:,:,f) )
        else
        xyzf_Turb(:,:,:,f) =                                                         &
          &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
          & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ) ) &
          & + xyz_dz_xyr(                                                            &
          &     xyr_xyz( xyz_DensBZKhBl )                                            &
          &     * ( xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) ) + xyrf_DQMixBZDz(:,:,:,f) )   &
          &   ) / xyz_DensBZ
        end if
      end do

      !--------------------------------
      ! Turb.u 
      !
      if (FlagArare4) then 
      pyz_Turb = &
        &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )  &
        & + pyz_dy_pqz(                                                &
        &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )         &
        &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )         & 
        &   )                                                          &
        & + pyz_dz_pyr(                                                &
        &       pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )         &
        &     + pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )         &
        &   )                                                          &
        & - 2.0d0 * pyz_dx_xyz( xyz_KmBlKmBl )                         &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
      else
      pyz_Turb =                                                       &
        &   2.0d0 * pyz_dx_xyz( xyz_KmBl * xyz_dx_pyz( pyz_VelXBl ) )  &
        & + pyz_dy_pqz(                                                &
        &       pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )         &
        &     + pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )         & 
        &   )                                                          &
        & + pyz_dz_pyr(                                                &
        &       pyr_xyz( xyz_DensBZKmBl ) * pyr_dx_xyr( xyr_VelZBl )   &
        &     + pyr_xyz( xyz_DensBZKmBl ) * pyr_dz_pyz( pyz_VelXBl )   &
        &   ) / pyz_DensBZ                                             &
        & - 2.0d0 * pyz_dx_xyz( xyz_KmBlKmBl )                         &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
      end if

      !-----------------------------------------------
      ! Turb.v
      !
      if (FlagArare4) then 
      xqz_Turb =                                                       &
        &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) )  &
        & + xqz_dx_pqz(                                                &
        &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )         &
        &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )         &
        &   )                                                          &
        & + xqz_dz_xqr(                                                &
        &       xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )         &
        &     + xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )         &
        &   )                                                          &
        & - 2.0d0 * xqz_dy_xyz( xyz_KmBlKmBl )                         &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
      else
      xqz_Turb =                                                       &
        &   2.0d0 * xqz_dy_xyz( xyz_KmBl * xyz_dy_xqz( xqz_VelYBl ) )  &
        & + xqz_dx_pqz(                                                &
        &       pqz_xyz( xyz_KmBl ) * pqz_dy_pyz( pyz_VelXBl )         &
        &     + pqz_xyz( xyz_KmBl ) * pqz_dx_xqz( xqz_VelYBl )         &
        &   )                                                          &
        & + xqz_dz_xqr(                                                &
        &       xqr_xyz( xyz_DensBZKmBl ) * xqr_dy_xyr( xyr_VelZBl )   &
        &     + xqr_xyz( xyz_DensBZKmBl ) * xqr_dz_xqz( xqz_VelYBl )   &
        &   ) / xqz_DensBZ                                             &
        & - 2.0d0 * xqz_dy_xyz( xyz_KmBlKmBl )                         &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
      end if

      !-----------------------------------------------
      ! Turb.w
      !
      if (FlagArare4) then 
      xyr_Turb =                                                            &
        & + 2.0d0 * xyr_dz_xyz( xyz_KmBl * xyz_dz_xyr( xyr_VelZBl ) )       &
        & + xyr_dx_pyr(                                                     &
        &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )               &
        &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )               &
        &   )                                                               &
        & + xyr_dy_xqr(                                                     &
        &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )               &
        &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )               &
        &   )                                                               &
        & - 2.0d0 * xyr_dz_xyz( xyz_KmBlKmBl )                              &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
      else
      xyr_Turb =                                                            &
        & + 2.0d0                                                           &
        &   * xyr_dz_xyz(                                                   &
        &       xyz_DensBZKmBl * xyz_dz_xyr( xyr_VelZBl )                   &
        &     ) / xyr_DensBZ                                                &
        & + xyr_dx_pyr(                                                     &
        &      pyr_xyz( xyz_KmBl ) * pyr_dz_pyz( pyz_VelXBl )               &
        &    + pyr_xyz( xyz_KmBl ) * pyr_dx_xyr( xyr_VelZBl )               &
        &   )                                                               &
        & + xyr_dy_xqr(                                                     &
        &      xqr_xyz( xyz_KmBl ) * xqr_dz_xqz( xqz_VelYBl )               &
        &    + xqr_xyz( xyz_KmBl ) * xqr_dy_xyr( xyr_VelZBl )               &
        &   )                                                               &
        & - 2.0d0 * xyr_dz_xyz( xyz_DensBZKmBlKmBl )                        &
        &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )                               &
        &   / xyr_DensBZ 
      end if

    end subroutine Turbulence_kw1978_std

    !!---------------------------------------------------------------------------!!

    subroutine turbulence_kw1978_center2_3d

      ! $B%b%8%e!<%k(B
      !
      use axesset,   only : dx, dy, dz
      use gridset,   only : imin,           &! x $BJ}8~$NG[Ns$N2<8B(B
        &                   imax,           &! x $BJ}8~$NG[Ns$N>e8B(B
        &                   jmin,           &! y $BJ}8~$NG[Ns$N2<8B(B
        &                   jmax,           &! y $BJ}8~$NG[Ns$N>e8B(B
        &                   kmin,           &! z $BJ}8~$NG[Ns$N2<8B(B
        &                   kmax,           &! z $BJ}8~$NG[Ns$N>e8B(B
        &                   ncmax,          &!
        &                   nx,ny,nz
      use basicset,  only : xyz_ExnerBZ,    &
        &                   xyz_DensBZ,     &
        &                   pyz_DensBZ,     &
        &                   xyr_DensBZ     
      use constants, only : CpDry

      ! $B0EL[$N7?@k8@6X;_(B
      !
      implicit none

      ! $B:n6HJQ?t(B
      ! 
      real(DP) :: xyz_KmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: pqz_KmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: pyr_KmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xqr_KmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKhBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: pyr_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xqr_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      integer  :: i, j, k, f


      ! Km $B$N<+>h(B
      !
      xyz_KmBlKmBl = xyz_KmBl * xyz_KmBl

      ! Km $B$N3J;RE@0LCV$NJQ99(B
      !
      do k = kmin, kmax
        do j = jmin, jmax-1
          do i = imin, imax-1          
            pqz_KmBl(i,j,k) =             &
              &  (                        &
              &     xyz_KmBl(i,   j,   k) &
              &   + xyz_KmBl(i+1, j,   k) &
              &   + xyz_KmBl(i,   j+1, k) &
              &   + xyz_KmBl(i+1, j+1, k) &
              &  ) * 0.25d0 
          end do
        end do
      end do
      
      ! Km $B$N3J;RE@0LCV$NJQ99(B 
      !
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax-1
            pyr_KmBl(i,j,k) =                &
              &  (                           &
              &    + xyz_KmBl(i,   j, k  )   &
              &    + xyz_KmBl(i+1, j, k  )   &
              &    + xyz_KmBl(i,   j, k+1)   &
              &    + xyz_KmBl(i+1, j, k+1)   &
              &  ) * 0.25d0 
          end do
        end do
      end do

      ! Km $B$N3J;RE@0LCV$NJQ99(B 
      !      
      do k = kmin, kmax-1
        do j = jmin, jmax-1
          do i = imin, imax
            xqr_KmBl(i,j,k) =              &
              &  (                         &
              &   + xyz_KmBl(i, j,   k  )  &
              &   + xyz_KmBl(i, j+1, k  )  &
              &   + xyz_KmBl(i, j,   k+1)  &
              &   + xyz_KmBl(i, j+1, k+1)  &
              &  ) * 0.25d0 
          end do
        end do
      end do

      ! Dens * Kh
      !
      xyz_DensBZKhBl = xyz_DensBZ * xyz_KhBl

      ! Dens * Km 
      !
      xyz_DensBZKmBl  = xyz_DensBZ * xyz_KmBl
      xyz_DensBZKmBlKmBl = xyz_DensBZ * xyz_KmBl * xyz_KmBl

      ! Dens * Km (pyr)
      !
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax-1            
            pyr_DensBZKmBl(i,j,k) =                &
              &  (                                 &
              &    + xyz_DensBZKmBl(i,   j, k  )   &
              &    + xyz_DensBZKmBl(i+1, j, k  )   &
              &    + xyz_DensBZKmBl(i,   j, k+1)   &
              &    + xyz_DensBZKmBl(i+1, j, k+1)   &
              &  ) * 0.25d0 
          end do
        end do
      end do      
      
      ! Dens * Km (xqr)
      !
      do k = kmin, kmax-1
        do j = jmin, jmax-1
          do i = imin, imax
            xqr_DensBZKmBl(i,j,k) =              &
              &  (                               &
              &   + xyz_DensBZKmBl(i, j,   k  )  &
              &   + xyz_DensBZKmBl(i, j+1, k  )  &
              &   + xyz_DensBZKmBl(i, j,   k+1)  &
              &   + xyz_DensBZKmBl(i, j+1, k+1)  &
              &  ) * 0.25d0 
          end do
        end do
      end do

      !---------------------------------------
      ! Km $B$N;~4VH/E8J}Dx<0$N3F9`(B
      !

      ! $B;60o9`(B
      !
      xyz_Disp = - xyz_KmBlKmBl * 5.0d-1 / ( MixLen ** 2.0d0 )

      ! $BIbNO9`(B
      !
      call turbulence_kw1978_BuoyancyKm
      
      ! $B%7%"!<9`(B
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            xyz_Shear(i,j,k) = &
              &   Cm_Cm_MixLen_MixLen                     &
              &   * (                                     &
              &     + (                                   &
              &         (                                 &
              &             pyz_VelXBl( i,  j, k )        &
              &           - pyz_VelXBl( i-1,j, k )        &
              &          ) / dx                           &
              &       ) ** 2.0d0                          &
              &     + (                                   &
              &         (                                 &
              &             xqz_VelYBl( i, j,   k )       &
              &           - xqz_VelYBl( i, j-1, k )       &
              &         ) / dy                            &
              &       ) ** 2.0d0                          &
              &     + (                                   &
              &         (                                 &
              &             xyr_VelZBl( i, j, k   )       &
              &           - xyr_VelZBl( i, j, k-1 )       &
              &         ) / dz                            &
              &       ) ** 2.0d0                          &
              &     + 5.0d-1                              &
              &       * (                                 &
              &           (                               &
              &           + (                             &
              &             + pyz_VelXBl( i,   j, k+1 )   &
              &             + pyz_VelXBl( i-1, j, k+1 )   &
              &             - pyz_VelXBl( i,   j, k-1 )   &
              &             - pyz_VelXBl( i-1, j, k-1 )   &
              &             ) * 0.25d0 / dz               &
              &           + (                             &
              &             + xyr_VelZBl( i+1, j, k   )   &
              &             + xyr_VelZBl( i+1, j, k-1 )   &
              &             - xyr_VelZBl( i-1, j, k   )   &
              &             - xyr_VelZBl( i-1, j, k-1 )   &
              &             ) * 0.25d0 / dx               &
              &           ) ** 2.0d0                      &
              &         + (                               &
              &           + (                             &
              &             + xyr_VelZBl( i, j+1, k   )   &
              &             + xyr_VelZBl( i, j+1, k-1 )   &
              &             - xyr_VelZBl( i, j-1, k   )   &
              &             - xyr_VelZBl( i, j-1, k-1 )   &
              &             ) * 0.25d0 / dy               &
              &           + (                             &
              &             + xqz_VelYBl( i, j,   k+1 )   &
              &             + xqz_VelYBl( i, j-1, k+1 )   &
              &             - xqz_VelYBl( i, j,   k-1 )   &
              &             - xqz_VelYBl( i, j-1, k-1 )   &
              &             ) * 0.25d0 / dz               &
              &           ) ** 2.0d0                      &
              &         + (                               &
              &           + (                             &
              &             + xqz_VelYBl( i+1, j,   k )   &
              &             + xqz_VelYBl( i+1, j-1, k )   &
              &             - xqz_VelYBl( i-1, j,   k )   &
              &             - xqz_VelYBl( i-1, j-1, k )   &
              &             ) * 0.25d0 / dx               &
              &           + (                             &
              &             + pyz_VelXBl( i,   j+1, k )   &
              &             + pyz_VelXBl( i-1, j+1, k )   &
              &             - pyz_VelXBl( i,   j-1, k )   &
              &             - pyz_VelXBl( i-1, j-1, k )   &
              &             ) * 0.25d0 / dy               &
              &           ) ** 2.0d0                      &
              &         )                                 &
              &     )                                     &
              & - xyz_KmBl( i, j, k )                     &
              &   * (                                     &
              &       (                                   &
              &         pyz_VelXBl( i,   j, k )           &
              &       - pyz_VelXBl( i-1, j, k )           &
              &       ) / dx                              &
              &     + (                                   &
              &         xqz_VelYBl( i, j,   k )           &
              &       - xqz_VelYBl( i, j-1, k )           &
              &       ) / dy                              &
              &     + (                                   &
              &         xyr_VelZBl( i, j, k   )           &
              &       - xyr_VelZBl( i, j, k-1 )           &
              &       ) / dz                              &
              &     ) / 3.0d0

          end do
        end do
      end do
            
      ! $B3H;69`(B
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            xyz_Diff(i,  j, k  ) =                              &
              & + 5.0d-1                                        &
              &   * (                                           &
              &      + (                                        &
              &         + xyz_KmBlKmBl( i+1, j,   k   )         &
              &         + xyz_KmBlKmBl( i-1, j,   k   )         &
              &         - xyz_KmBlKmBl( i,   j,   k   ) * 2.0d0 &
              &        ) / ( dx * dx )                          &
              &      + (                                        &
              &         + xyz_KmBlKmBl( i,   j+1, k   )         &
              &         + xyz_KmBlKmBl( i,   j-1, k   )         &
              &         - xyz_KmBlKmBl( i,   j,   k   ) * 2.0d0 &
              &        ) / ( dy * dy )                          &
              &      + (                                        &
              &         + xyz_KmBlKmBl( i,   j,   k+1 )         &
              &         + xyz_KmBlKmBl( i,   j,   k-1 )         &
              &         - xyz_KmBlKmBl( i,   j,   k   ) * 2.0d0 &
              &        ) / ( dz * dz )                          &
              &     )                                           &
              & + (                                             &
              &    + (                                          &
              &        (                                        &
              &         + xyz_KmBl( i+1, j,   k   )             &
              &         - xyz_KmBl( i-1, j,   k   )             &
              &        ) * 0.5d0 / dx                           &
              &      ) ** 2.0d0                                 & 
              &    + (                                          &
              &        (                                        &
              &         + xyz_KmBl( i,   j+1, k   )             &
              &         - xyz_KmBl( i,   j-1, k   )             &
              &        ) * 0.5d0 / dy                           &
              &      ) ** 2.0d0                                 & 
              &    + (                                          &
              &        (                                        &
              &         + xyz_KmBl( i,   j,   k+1 )             &
              &         - xyz_KmBl( i,   j,   k-1 )             &
              &        ) * 0.5d0 / dz                           &
              &      ) ** 2.0d0                                 &
              &   )
            
          end do
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.Dens
      !
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            xyz_DCDensDtTurb(i,  j,  k  ) =                  &
              & + (                                       &
              &    + (                                    &
              &         xyz_KhBl(i+1,j,  k  )             &
              &       + xyz_KhBl(i,  j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i+1,j,  k  )          &
              &       - xyz_CDensBl(i,  j,  k  )          &
              &      )                                    &
              &    - (                                    &
              &         xyz_KhBl(i,  j,  k  )             &
              &       + xyz_KhBl(i-1,j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i,  j,  k  )          &
              &       - xyz_CDensBl(i-1,j,  k  )          &
              &      )                                    &
              &   ) * 0.5d0 / ( dx * dx )                 &
              & + (                                       &
              &    + (                                    &
              &         xyz_KhBl(i,  j+1,k  )             &
              &       + xyz_KhBl(i,  j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i,  j+1,k  )          & 
              &       - xyz_CDensBl(i,  j,  k  )          &
              &      )                                    &
              &    - (                                    &
              &         xyz_KhBl(i,  j,  k  )             &
              &       + xyz_KhBl(i,  j-1,k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i,  j,  k  )          &
              &       - xyz_CDensBl(i,  j-1,k  )          &
              &      )                                    &
              &   ) * 0.5d0 / ( dy * dy )                 &
              & + (                                       &
              &    + (                                    &
              &         xyz_DensBZKhBl(i,  j,  k+1)       &
              &       + xyz_DensBZKhBl(i,  j,  k  )       &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i,  j,  k+1)          &
              &       - xyz_CDensBl(i,  j,  k  )          &
              &      )                                    &
              &    - (                                    &
              &         xyz_DensBZKhBl(i,  j,  k  )       &
              &       + xyz_DensBZKhBl(i,  j,  k-1)       &
              &      )                                    &
              &    * (                                    &
              &         xyz_CDensBl(i,  j,  k  )          &
              &       - xyz_CDensBl(i,  j,  k-1)          &
              &      )                                    &
              &   ) * 0.5d0 / ( dz  * dz )                &
              &     / xyz_DensBZ(i,  j,  k  )
            
          end do
        end do
      end do


      !-----------------------------------------------
      ! Turb.\theta
      !
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            xyz_Turb(i,  j,  k  ) =                       &
              & + (                                       &
              &    + (                                    &
              &         xyz_KhBl(i+1,j,  k  )             &
              &       + xyz_KhBl(i,  j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_PTempBl(i+1,j,  k  )          &
              &       - xyz_PTempBl(i,  j,  k  )          &
              &      )                                    &
              &    - (                                    &
              &         xyz_KhBl(i,  j,  k  )             &
              &       + xyz_KhBl(i-1,j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_PTempBl(i,  j,  k  )          &
              &       - xyz_PTempBl(i-1,j,  k  )          &
              &      )                                    &
              &   ) * 0.5d0 / ( dx * dx )                 &
              & + (                                       &
              &    + (                                    &
              &         xyz_KhBl(i,  j+1,k  )             &
              &       + xyz_KhBl(i,  j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_PTempBl(i,  j+1,k  )          & 
              &       - xyz_PTempBl(i,  j,  k  )          &
              &      )                                    &
              &    - (                                    &
              &         xyz_KhBl(i,  j,  k  )             &
              &       + xyz_KhBl(i,  j-1,k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyz_PTempBl(i,  j,  k  )          &
              &       - xyz_PTempBl(i,  j-1,k  )          &
              &      )                                    &
              &   ) * 0.5d0 / ( dy * dy )                 &
              & + (                                       &
              &    + (                                    &
              &         xyz_DensBZKhBl(i,  j,  k+1)       &
              &       + xyz_DensBZKhBl(i,  j,  k  )       &
              &      )                                    &
              &    * (                                    &
              &       + (                                 &
              &            xyz_PTempBl(i,  j,  k+1)       &
              &          - xyz_PTempBl(i,  j,  k  )       &
              &         ) / dz                            &
              &       + xyr_DPTempBZDz(i,  j,  k )        &
              &      )                                    &
              &    - (                                    &
              &         xyz_DensBZKhBl(i,  j,  k  )       &
              &       + xyz_DensBZKhBl(i,  j,  k-1)       &
              &      )                                    &
              &    * (                                    &
              &       + (                                 &
              &            xyz_PTempBl(i,  j,  k  )       &
              &          - xyz_PTempBl(i,  j,  k-1)       &
              &         ) / dz                            &
              &       + xyr_DPTempBZDz(i,  j,  k-1 )      &
              &      )                                    &
              &   ) * 0.5d0 / dz                          &
              &     / xyz_DensBZ(i,  j,  k  )
            
          end do
        end do
      end do
      
      xyz_DispHeat = (xyz_KmBl ** 3.0d0) &
        & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0))

      !-----------------------------------------------
      ! Turb.Q
      !
      
      do f = 1, ncmax
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              
              xyzf_Turb(i,  j,  k,  f  ) =                  &
                & + (                                       &
                &    + (                                    &
                &         xyz_KhBl(i+1,j,  k  )             &
                &       + xyz_KhBl(i,  j,  k  )             &
                &      )                                    &
                &    * (                                    &
                &         xyzf_QMixBl(i+1,j,  k,  f  )      &
                &       - xyzf_QMixBl(i,  j,  k,  f  )      &
                &      )                                    &
                &    - (                                    &
                &         xyz_KhBl(i,  j,  k  )             &
                &       + xyz_KhBl(i-1,j,  k  )             &
                &      )                                    &
                &    * (                                    &
                &         xyzf_QMixBl(i,  j,  k,  f  )      &
                &       - xyzf_QMixBl(i-1,j,  k,  f  )      &
                &      )                                    &
                &   ) * 0.5d0 / ( dx * dx )                 &
                & + (                                       &
                &    + (                                    &
                &         xyz_KhBl(i,  j+1,k  )             &
                &       + xyz_KhBl(i,  j,  k  )             &
                &      )                                    &
                &    * (                                    &
                &         xyzf_QMixBl(i,  j+1,k,  f  )      &
                &       - xyzf_QMixBl(i,  j,  k,  f  )      &
                &      )                                    &
                &    - (                                    &
                &         xyz_KhBl(i,  j,  k  )             &
                &       + xyz_KhBl(i,  j-1,k  )             &
                &      )                                    &
                &    * (                                    &
                &         xyzf_QMixBl(i,  j,  k,  f  )      &
                &       - xyzf_QMixBl(i,  j-1,k,  f  )      &
                &      )                                    &
                &   ) * 0.5d0 / ( dy * dy )                 &
                & + (                                       &
                &    + (                                    &
                &         xyz_DensBZKhBl(i,  j,  k+1)       &
                &       + xyz_DensBZKhBl(i,  j,  k  )       &
                &      )                                    &
                &    * (                                    &
                &       + (                                 & 
                &            xyzf_QMixBl(i,  j,  k+1,f  )   &
                &          - xyzf_QMixBl(i,  j,  k,  f  )   &
                &         ) / dz                            &
                &       + xyrf_DQMixBZDz(i,  j,  k,  f  )   &
                &      )                                    &
                &    - (                                    &
                &         xyz_DensBZKhBl(i,  j,  k  )       &
                &       + xyz_DensBZKhBl(i,  j,  k-1)       &
                &      )                                    &
                &    * (                                    &
                &       + (                                 & 
                &            xyzf_QMixBl(i,  j,  k,  f  )   &
                &          - xyzf_QMixBl(i,  j,  k-1,f  )   &
                &         ) / dz                            &
                &       + xyrf_DQMixBZDz(i,  j,  k-1,f  )   &
                &      )                                    &
                &   ) * 0.5d0 / dz                          &
                &     / xyz_DensBZ(i,  j,  k  )
              
            end do
          end do
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.u 
      !
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            pyz_Turb(i,j,k) =                                                          &
              & + 2.0d0                                                                &
              &   * (                                                                  &
              &      + xyz_KmBl(i+1,j,  k  )                                           &
              &        * ( pyz_VelXBl(i+1,j,  k  ) - pyz_VelXBl(i,  j,  k  ) )         &
              &      - xyz_KmBl(i,  j,  k  )                                           &
              &        * ( pyz_VelXBl(i,  j,  k  ) - pyz_VelXBl(i-1,j,  k  ) )         &
              &     ) / ( dx * dx )                                                    &
              & + (                                                                    &
              &      + pqz_KmBl( i,  j,  k  )                                          &
              &        * (                                                             &
              &             ( xqz_VelYBl(i+1,j,  k  ) - xqz_VelYBl(i,  j,  k  ) ) / dx &
              &           + ( pyz_VelXBl(i,  j+1,k  ) - pyz_VelXBl(i,  j,  k  ) ) / dy &
              &          )                                                             &
              &      - pqz_KmBl( i,  j-1,k  )                                          &
              &        * (                                                             &
              &             ( xqz_VelYBl(i+1,j-1,k  ) - xqz_VelYBl(i,  j-1,k  ) ) / dx &
              &           + ( pyz_VelXBl(i,  j,  k  ) - pyz_VelXBl(i,  j-1,k  ) ) / dy &
              &          )                                                             &
              &   ) / dy                                                               &
              & + (                                                                    &
              &      + pyr_DensBZKmBl(i,  j,  k  )                                     &
              &        * (                                                             &
              &             ( xyr_VelZBl(i+1,j,  k  ) - xyr_VelZBl(i,  j,  k  ) ) / dx &
              &           + ( pyz_VelXBl(i,  j,  k+1) - pyz_VelXBl(i,  j,  k  ) ) / dz &
              &          )                                                             &
              &      - pyr_DensBZKmBl(i,  j,  k-1)                                     &
              &        * (                                                             &
              &             ( xyr_VelZBl(i+1,j,  k-1) - xyr_VelZBl(i,  j,  k-1) ) / dx &
              &           + ( pyz_VelXBl(i,  j,  k  ) - pyz_VelXBl(i,  j,  k-1) ) / dz &
              &          )                                                             &
              &   ) / dz / pyz_DensBZ(i,  j,  k)                                       &
              & - 2.0d0                                                                &
              &   * ( xyz_KmBlKmBl(i+1,j,  k  ) - xyz_KmBlKmBl(i,  j,  k  ) ) / dx     &
              &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
            
          end do
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.v
      !
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            xqz_Turb(i,j,k) =                                                          &
              & + 2.0d0                                                                &
              &   * (                                                                  &
              &      + xyz_KmBl(i,  j+1,k  )                                           &
              &        * ( xqz_VelYBl(i,  j+1,k  ) - xqz_VelYBl(i,  j,  k  ) )         &
              &      - xyz_KmBl(i,  j,  k  )                                           &
              &        * ( xqz_VelYBl(i,  j,  k  ) - xqz_VelYBl(i,  j-1,k  ) )         &
              &     ) / ( dy * dy )                                                    &
              & + (                                                                    &
              &      + pqz_KmBl( i,  j,  k  )                                          &
              &        * (                                                             &
              &             ( pyz_VelXBl(i,  j+1,k  ) - pyz_VelXBl(i,  j,  k  ) ) / dy &
              &           + ( xqz_VelYBl(i+1,j,  k  ) - xqz_VelYBl(i,  j,  k  ) ) / dx &
              &          )                                                             &
              &      - pqz_KmBl( i-1,j,  k  )                                          &
              &        * (                                                             &
              &             ( pyz_VelXBl(i-1,j+1,k  ) - pyz_VelXBl(i-1,j,  k  ) ) / dy &
              &           + ( xqz_VelYBl(i,  j,  k  ) - xqz_VelYBl(i-1,j,  k  ) ) / dx &
              &          )                                                             &
              &   ) / dx                                                               &
              & + (                                                                    &
              &      + xqr_DensBZKmBl(i,  j,  k  )                                     &
              &        * (                                                             &
              &             ( xyr_VelZBl(i,  j+1,k  ) - xyr_VelZBl(i,  j,  k  ) ) / dy &
              &           + ( xqz_VelYBl(i,  j,  k+1) - xqz_VelYBl(i,  j,  k  ) ) / dz &
              &          )                                                             &
              &      - xqr_DensBZKmBl(i,  j,  k-1)                                     &
              &        * (                                                             &
              &             ( xyr_VelZBl(i,  j+1,k-1) - xyr_VelZBl(i,  j,  k-1) ) / dy &
              &           + ( xqz_VelYBl(i,  j  ,k  ) - xqz_VelYBl(i,  j,  k-1) ) / dz &
              &          )                                                             &
              &   ) / dz / xyr_DensBZ(i,  j,  k)                                       &
              & - 2.0d0                                                                &
              &   * ( xyz_KmBlKmBl(i,  j+1,k  ) - xyz_KmBlKmBl(i,  j,  k  ) ) / dy     &
              &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
            
          end do
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.w
      !
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            xyr_Turb(i,  j,  k  ) =                                                    &
              & + 2.0d0                                                                &
              &   * (                                                                  &
              &      + xyz_DensBZKmBl(i,  j,  k+1)                                     &
              &        * ( xyr_VelZBl(i,  j,  k+1) - xyr_VelZBl(i,  j,  k  ) )         &
              &      - xyz_DensBZKmBl(i,  j,  k  )                                     &
              &        * ( xyr_VelZBl(i,  j,  k  ) - xyr_VelZBl(i,  j,  k-1) )         &
              &     ) / ( dz * dz )                                                    &
              &   / xyr_DensBZ(i,  j,  k  )                                            &
              & + (                                                                    &
              &      + pyr_KmBl( i,  j,  k  )                                          &
              &        * (                                                             &
              &             ( pyz_VelXBl(i,  j,  k+1) - pyz_VelXBl(i,  j,  k  ) ) / dz &
              &           + ( xyr_VelZBl(i+1,j,  k  ) - xyr_VelZBl(i,  j,  k  ) ) / dx &
              &          )                                                             &
              &      - pyr_KmBl( i-1,j,  k  )                                          &
              &        * (                                                             &
              &             ( pyz_VelXBl(i-1,j,  k+1) - pyz_VelXBl(i-1,j,  k  ) ) / dz &
              &           + ( xyr_VelZBl(i,  j,  k  ) - xyr_VelZBl(i-1,j,  k  ) ) / dx &
              &          )                                                             &
              &   ) / dx                                                               &
              & + (                                                                    &
              &      + xqr_KmBl(i,  j,  k  )                                           &
              &        * (                                                             &
              &             ( xqz_VelYBl(i,  j,  k+1) - xqz_VelYBl(i,  j,  k  ) ) / dz &
              &           + ( xyr_VelZBl(i,  j+1,k  ) - xyr_VelZBl(i,  j,  k  ) ) / dy &
              &          )                                                             &
              &      - xqr_KmBl(i,  j-1,k  )                                           &
              &        * (                                                             &
              &             ( xqz_VelYBl(i,  j-1,k+1) - xqz_VelYBl(i,  j-1,k  ) ) / dz &
              &           + ( xyr_VelZBl(i,  j,  k  ) - xyr_VelZBl(i,  j-1,k  ) ) / dy &
              &          )                                                             &
              &   ) / dy                                                               &
              & - 2.0d0                                                                &
              &   * (                                                                  &
              &        xyz_DensBZKmBlKmBl(i,  j,  k+1)                                 &
              &      - xyz_DensBZKmBlKmBl(i,  j,  k  )                                 &
              &     ) / dz                                                             &
              &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )                                    &
              &   / xyr_DensBZ(i,  j,  k  )                                      
            
          end do
        end do
      end do
      
    end subroutine Turbulence_kw1978_center2_3d

    !!---------------------------------------------------------------------------!!

    subroutine turbulence_kw1978_center2_2d

      ! $B%b%8%e!<%k(B
      !
      use axesset,   only : dx, dz
      use gridset,   only : imin,           &! x $BJ}8~$NG[Ns$N2<8B(B
        &                   imax,           &! x $BJ}8~$NG[Ns$N>e8B(B
        &                   jmin,           &! y $BJ}8~$NG[Ns$N2<8B(B
        &                   jmax,           &! y $BJ}8~$NG[Ns$N>e8B(B
        &                   kmin,           &! z $BJ}8~$NG[Ns$N2<8B(B
        &                   kmax,           &! z $BJ}8~$NG[Ns$N>e8B(B
        &                   ncmax,          &!
        &                   nx, nz
      use basicset,  only : xyz_ExnerBZ,    &
        &                   xyz_DensBZ,     &
        &                   pyz_DensBZ,     &
        &                   xyr_DensBZ     
      use constants, only : CpDry

      ! $B0EL[$N7?@k8@6X;_(B
      !
      implicit none

      ! $B:n6HJQ?t(B
      ! 
      real(DP) :: xyz_KmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: pyr_KmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKhBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: xyz_DensBZKmBlKmBl(imin:imax, jmin:jmax, kmin:kmax)
      real(DP) :: pyr_DensBZKmBl(imin:imax, jmin:jmax, kmin:kmax)
      integer  :: i, k, f
      integer, parameter :: j = 1


      ! Km $B$N<+>h(B
      !
      xyz_KmBlKmBl = xyz_KmBl * xyz_KmBl

      ! Km $B$N3J;RE@0LCV$NJQ99(B 
      !
      do k = kmin, kmax-1
        do i = imin, imax-1
          pyr_KmBl(i,j,k) =                &
            &  (                           &
            &    + xyz_KmBl(i,   j, k  )   &
            &    + xyz_KmBl(i+1, j, k  )   &
            &    + xyz_KmBl(i,   j, k+1)   &
            &    + xyz_KmBl(i+1, j, k+1)   &
            &  ) * 0.25d0 
        end do
      end do

      ! Dens * Kh
      !
      xyz_DensBZKhBl = xyz_DensBZ * xyz_KhBl

      ! Dens * Km 
      !
      xyz_DensBZKmBl  = xyz_DensBZ * xyz_KmBl
      xyz_DensBZKmBlKmBl = xyz_DensBZ * xyz_KmBl * xyz_KmBl

      ! Dens * Km (pyr)
      !
      do k = kmin, kmax-1
        do i = imin, imax-1            
          pyr_DensBZKmBl(i,j,k) =                &
            &  (                                 &
            &    + xyz_DensBZKmBl(i,   j, k  )   &
            &    + xyz_DensBZKmBl(i+1, j, k  )   &
            &    + xyz_DensBZKmBl(i,   j, k+1)   &
            &    + xyz_DensBZKmBl(i+1, j, k+1)   &
            &  ) * 0.25d0 
        end do
      end do      
      
      !---------------------------------------
      ! Km $B$N;~4VH/E8J}Dx<0$N3F9`(B
      !

      ! $B;60o9`(B
      !
      xyz_Disp = - xyz_KmBlKmBl * 5.0d-1 / ( MixLen ** 2.0d0 )

      ! $BIbNO9`(B
      !
      call turbulence_kw1978_BuoyancyKm
      
      ! $B%7%"!<9`(B
      !
      do k = 1, nz
        do i = 1, nx
          
          xyz_Shear(i,j,k) = &
!            &   Cm_Cm_MixLen_MixLen                     &
            &   ( ( Cm ** 2.0d0 ) * ( MixLen ** 2.0d0 ) )  &
            &   * (                                     &
            &     + (                                   &
            &         (                                 &
            &             pyz_VelXBl( i,  j, k )        &
            &           - pyz_VelXBl( i-1,j, k )        &
            &          ) / dx                           &
            &       ) ** 2.0d0                          &
            &     + (                                   &
            &         (                                 &
            &             xyr_VelZBl( i, j, k   )       &
            &           - xyr_VelZBl( i, j, k-1 )       &
            &         ) / dz                            &
            &       ) ** 2.0d0                          &
            &     + 5.0d-1                              &
            &       * (                                 &
            &           (                               &
            &           + (                             &
            &             + pyz_VelXBl( i,   j, k+1 )   &
            &             + pyz_VelXBl( i-1, j, k+1 )   &
            &             - pyz_VelXBl( i,   j, k-1 )   &
            &             - pyz_VelXBl( i-1, j, k-1 )   &
            &             ) * 0.25d0 / dz               &
            &           + (                             &
            &             + xyr_VelZBl( i+1, j, k   )   &
            &             + xyr_VelZBl( i+1, j, k-1 )   &
            &             - xyr_VelZBl( i-1, j, k   )   &
            &             - xyr_VelZBl( i-1, j, k-1 )   &
            &             ) * 0.25d0 / dx               &
            &           ) ** 2.0d0                      &
            &         )                                 &
            &     )                                     &
            & - xyz_KmBl( i, j, k )                     &
            &   * (                                     &
            &       (                                   &
            &         pyz_VelXBl( i,   j, k )           &
            &       - pyz_VelXBl( i-1, j, k )           &
            &       ) / dx                              &
            &     + (                                   &
            &         xyr_VelZBl( i, j, k   )           &
            &       - xyr_VelZBl( i, j, k-1 )           &
            &       ) / dz                              &
            &     ) / 3.0d0
          
        end do
      end do
      
      ! $B3H;69`(B
      !
      do k = 1, nz
        do i = 1, nx

          xyz_Diff(i,  j, k  ) =                              &
            & + 5.0d-1                                        &
            &   * (                                           &
            &      + (                                        &
            &         + xyz_KmBlKmBl( i+1, j,   k   )         &
            &         + xyz_KmBlKmBl( i-1, j,   k   )         &
            &         - xyz_KmBlKmBl( i,   j,   k   ) * 2.0d0 &
            &        ) / ( dx * dx )                          &
            &      + (                                        &
            &         + xyz_KmBlKmBl( i,   j,   k+1 )         &
            &         + xyz_KmBlKmBl( i,   j,   k-1 )         &
            &         - xyz_KmBlKmBl( i,   j,   k   ) * 2.0d0 &
            &        ) / ( dz * dz )                          &
            &     )                                           &
            & + (                                             &
            &    + (                                          &
            &        (                                        &
            &         + xyz_KmBl( i+1, j,   k   )             &
            &         - xyz_KmBl( i-1, j,   k   )             &
            &        ) * 0.5d0 / dx                           &
            &      ) ** 2.0d0                                 & 
            &    + (                                          &
            &        (                                        &
            &         + xyz_KmBl( i,   j,   k+1 )             &
            &         - xyz_KmBl( i,   j,   k-1 )             &
            &        ) * 0.5d0 / dz                           &
            &      ) ** 2.0d0                                 &
            &   )
          
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.Dens
      !
      do k = 1, nz
        do i = 1, nx
          
          xyz_DCDensDtTurb(i,  j,  k  ) =                  &
            & + (                                       &
            &    + (                                    &
            &         xyz_KhBl(i+1,j,  k  )             &
            &       + xyz_KhBl(i,  j,  k  )             &
            &      )                                    &
            &    * (                                    &
            &         xyz_CDensBl(i+1,j,  k  )          &
            &       - xyz_CDensBl(i,  j,  k  )          &
            &      )                                    &
            &    - (                                    &
            &         xyz_KhBl(i,  j,  k  )             &
            &       + xyz_KhBl(i-1,j,  k  )             &
            &      )                                    &
            &    * (                                    &
            &         xyz_CDensBl(i,  j,  k  )          &
            &       - xyz_CDensBl(i-1,j,  k  )          &
            &      )                                    &
            &   ) * 0.5d0 / ( dx * dx )                 &
            & + (                                       &
            &    + (                                    &
            &         xyz_DensBZKhBl(i,  j,  k+1)       &
            &       + xyz_DensBZKhBl(i,  j,  k  )       &
            &      )                                    &
            &    * (                                    &
            &         xyz_CDensBl(i,  j,  k+1)          &
            &       - xyz_CDensBl(i,  j,  k  )          &
            &      )                                    &
            &    - (                                    &
            &         xyz_DensBZKhBl(i,  j,  k  )       &
            &       + xyz_DensBZKhBl(i,  j,  k-1)       &
            &      )                                    &
            &    * (                                    &
            &         xyz_CDensBl(i,  j,  k  )          &
            &       - xyz_CDensBl(i,  j,  k-1)          &
            &      )                                    &
            &   ) * 0.5d0 / ( dz  * dz )                &
            &     / xyz_DensBZ(i,  j,  k  )

        end do
      end do


      !-----------------------------------------------
      ! Turb.\theta
      !
      do k = 1, nz
        do i = 1, nx
          
          xyz_Turb(i,  j,  k  ) =                       &
            & + (                                       &
            &    + (                                    &
            &         xyz_KhBl(i+1,j,  k  )             &
            &       + xyz_KhBl(i,  j,  k  )             &
            &      )                                    &
            &    * (                                    &
            &         xyz_PTempBl(i+1,j,  k  )          &
            &       - xyz_PTempBl(i,  j,  k  )          &
            &      )                                    &
            &    - (                                    &
            &         xyz_KhBl(i,  j,  k  )             &
            &       + xyz_KhBl(i-1,j,  k  )             &
            &      )                                    &
            &    * (                                    &
            &         xyz_PTempBl(i,  j,  k  )          &
            &       - xyz_PTempBl(i-1,j,  k  )          &
            &      )                                    &
            &   ) * 0.5d0 / ( dx * dx )                 &
            & + (                                       &
            &    + (                                    &
            &         xyz_DensBZKhBl(i,  j,  k+1)       &
            &       + xyz_DensBZKhBl(i,  j,  k  )       &
            &      )                                    &
            &    * (                                    &
            &       + (                                 &
            &            xyz_PTempBl(i,  j,  k+1)       &
            &          - xyz_PTempBl(i,  j,  k  )       &
            &         ) / dz                            &
            &       + xyr_DPTempBZDz(i,  j,  k )        &
            &      )                                    &
            &    - (                                    &
            &         xyz_DensBZKhBl(i,  j,  k  )       &
            &       + xyz_DensBZKhBl(i,  j,  k-1)       &
            &      )                                    &
            &    * (                                    &
            &       + (                                 &
            &            xyz_PTempBl(i,  j,  k  )       &
            &          - xyz_PTempBl(i,  j,  k-1)       &
            &         ) / dz                            &
            &       + xyr_DPTempBZDz(i,  j,  k-1 )      &
            &      )                                    &
            &   ) * 0.5d0 / dz                          &
            &     / xyz_DensBZ(i,  j,  k  )

        end do
      end do
      
      xyz_DispHeat = (xyz_KmBl ** 3.0d0) &
        & / (xyz_ExnerBZ * CpDry * (Cm ** 2.0d0) * (MixLen ** 4.0d0))

      !-----------------------------------------------
      ! Turb.Q
      !
      
      do f = 1, ncmax
        do k = 1, nz
          do i = 1, nx
              
            xyzf_Turb(i,  j,  k,  f  ) =                  &
              & + (                                       &
              &    + (                                    &
              &         xyz_KhBl(i+1,j,  k  )             &
              &       + xyz_KhBl(i,  j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyzf_QMixBl(i+1,j,  k,  f  )      &
              &       - xyzf_QMixBl(i,  j,  k,  f  )      &
              &      )                                    &
              &    - (                                    &
              &         xyz_KhBl(i,  j,  k  )             &
              &       + xyz_KhBl(i-1,j,  k  )             &
              &      )                                    &
              &    * (                                    &
              &         xyzf_QMixBl(i,  j,  k,  f  )      &
              &       - xyzf_QMixBl(i-1,j,  k,  f  )      &
              &      )                                    &
              &   ) * 0.5d0 / ( dx * dx )                 &
              & + (                                       &
              &    + (                                    &
              &         xyz_DensBZKhBl(i,  j,  k+1)       &
              &       + xyz_DensBZKhBl(i,  j,  k  )       &
              &      )                                    &
              &    * (                                    &
              &       + (                                 & 
              &            xyzf_QMixBl(i,  j,  k+1,f  )   &
              &          - xyzf_QMixBl(i,  j,  k,  f  )   &
              &         ) / dz                            &
              &       + xyrf_DQMixBZDz(i,  j,  k,  f  )   &
              &      )                                    &
              &    - (                                    &
              &         xyz_DensBZKhBl(i,  j,  k  )       &
              &       + xyz_DensBZKhBl(i,  j,  k-1)       &
              &      )                                    &
              &    * (                                    &
              &       + (                                 & 
              &            xyzf_QMixBl(i,  j,  k,  f  )   &
              &          - xyzf_QMixBl(i,  j,  k-1,f  )   &
              &         ) / dz                            &
              &       + xyrf_DQMixBZDz(i,  j,  k-1,f  )   &
              &      )                                    &
              &   ) * 0.5d0 / dz                          &
              &     / xyz_DensBZ(i,  j,  k  )
            
          end do
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.u 
      !
      
      do k = 1, nz
        do i = 1, nx
          
          pyz_Turb(i,j,k) =                                                          &
            & + 2.0d0                                                                &
            &   * (                                                                  &
            &      + xyz_KmBl(i+1,j,  k  )                                           &
            &        * ( pyz_VelXBl(i+1,j,  k  ) - pyz_VelXBl(i,  j,  k  ) )         &
            &      - xyz_KmBl(i,  j,  k  )                                           &
            &        * ( pyz_VelXBl(i,  j,  k  ) - pyz_VelXBl(i-1,j,  k  ) )         &
            &     ) / ( dx * dx )                                                    &
            & + (                                                                    &
            &      + pyr_DensBZKmBl(i,  j,  k  )                                     &
            &        * (                                                             &
            &             ( xyr_VelZBl(i+1,j,  k  ) - xyr_VelZBl(i,  j,  k  ) ) / dx &
            &           + ( pyz_VelXBl(i,  j,  k+1) - pyz_VelXBl(i,  j,  k  ) ) / dz &
            &          )                                                             &
            &      - pyr_DensBZKmBl(i,  j,  k-1)                                     &
            &        * (                                                             &
            &             ( xyr_VelZBl(i+1,j,  k-1) - xyr_VelZBl(i,  j,  k-1) ) / dx &
            &           + ( pyz_VelXBl(i,  j,  k  ) - pyz_VelXBl(i,  j,  k-1) ) / dz &
            &          )                                                             &
            &   ) / dz / pyz_DensBZ(i,  j,  k)                                       &
            & - 2.0d0                                                                &
            &   * ( xyz_KmBlKmBl(i+1,j,  k  ) - xyz_KmBlKmBl(i,  j,  k  ) ) / dx     &
            &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )
          
        end do
      end do
      
      
      !-----------------------------------------------
      ! Turb.v
      !
      xqz_Turb = 0.0d0
      
      !-----------------------------------------------
      ! Turb.w
      !
      
      do k = 1, nz
        do i = 1, nx
            
          xyr_Turb(i,  j,  k  ) =                                                    &
            & + 2.0d0                                                                &
            &   * (                                                                  &
            &      + xyz_DensBZKmBl(i,  j,  k+1)                                     &
            &        * ( xyr_VelZBl(i,  j,  k+1) - xyr_VelZBl(i,  j,  k  ) )         &
            &      - xyz_DensBZKmBl(i,  j,  k  )                                     &
            &        * ( xyr_VelZBl(i,  j,  k  ) - xyr_VelZBl(i,  j,  k-1) )         &
            &     ) / ( dz * dz )                                                    &
            &   / xyr_DensBZ(i,  j,  k  )                                            &
            & + (                                                                    &
            &      + pyr_KmBl( i,  j,  k  )                                          &
            &        * (                                                             &
            &             ( pyz_VelXBl(i,  j,  k+1) - pyz_VelXBl(i,  j,  k  ) ) / dz &
            &           + ( xyr_VelZBl(i+1,j,  k  ) - xyr_VelZBl(i,  j,  k  ) ) / dx &
            &          )                                                             &
            &      - pyr_KmBl( i-1,j,  k  )                                          &
            &        * (                                                             &
            &             ( pyz_VelXBl(i-1,j,  k+1) - pyz_VelXBl(i-1,j,  k  ) ) / dz &
            &           + ( xyr_VelZBl(i,  j,  k  ) - xyr_VelZBl(i-1,j,  k  ) ) / dx &
            &          )                                                             &
            &   ) / dx                                                               &
            & - 2.0d0                                                                &
            &   * (                                                                  &
            &        xyz_DensBZKmBlKmBl(i,  j,  k+1)                                 &
            &      - xyz_DensBZKmBlKmBl(i,  j,  k  )                                 &
            &     ) / dz                                                             &
            &   / ( 3.0d0 * Cm_Cm_MixLen_MixLen )                                    &
            &   / xyr_DensBZ(i,  j,  k  )                                      
          
        end do
      end do
      
    end subroutine Turbulence_kw1978_center2_2d

!!!------------------------------------------------------------------------------!!!

    subroutine turbulence_kw1978_BuoyancyKm
      !
      
    use ChemCalc,    only: xyz_LatentHeat     !$B@xG.(B
    use basicset,    only: xyz_PTempBZ,         &!$B4pK\>l$N290L(B
      &                    xyz_VPTempBZ,        &! $B2>290L(B
      &                    xyzf_QMixBZ,         &!$B4pK\>l$N:.9gHf(B
      &                    xyz_QMixBZ,          &!$B4pK\>l$N:.9gHf(B
      &                    xyz_QMixBZPerMolWt,  &!$B4pK\>l$N:.9gHf(B
      &                    xyz_EffMolWtBZ,      &!$B4pK\>l$N:.9gHf(B
      &                    xyz_ExnerBZ           !$B4pK\>l$N%(%/%9%J!<4X?t(B
    use constants,   only: Grav,               &
      &                    MolWtDry,           & 
      &                    CpDry
    use composition, only: MolWtWet,           &
      &                    SpcWetID,           &
      &                    CondNum,            &!$B6E7k2aDx$N?t(B
      &                    GasNum,             &!$B5$BN$N?t(B
      &                    IdxCG,              &!$B6E7k2aDx(B($B>x5$(B)$B$NG[NsE:$(;z(B
      &                    IdxCC,              &!$B6E7k2aDx(B($B1@(B)$B$NG[NsE:$(;z(B
      &                    IdxG                 !$B5$BN$NG[NsE:$(;z(B
    use average,     only: xyz_xyr
    use differentiate_center2,                 &
      &              only: xyr_dz_xyz 
     
      
      !$B0EL[$N7?@k8@6X;_(B
      implicit none
      
      !$BJQ?tDj5A(B
      !
      real(DP) :: xyzf_LatentHeat(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !$B@xG.(B
      real(DP) :: xyz_TempBlAll(imin:imax,jmin:jmax,kmin:kmax)
                                               !$B29EY(B
      real(DP) :: xyz_EffHeat(imin:imax,jmin:jmax,kmin:kmax)
                                               !
      real(DP) :: xyz_EffPTemp(imin:imax,jmin:jmax,kmin:kmax)    
                                               !$BIbNO$KBP$9$k29EY:9$N4sM?(B
      real(DP) :: xyz_EffMolWt(imin:imax,jmin:jmax,kmin:kmax)    
                                               !$BIbNO$KBP$9$kJ,;RNL:9$N4sM?(B
      real(DP) :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, GasNum)
                                               !$B:.9gHf(B/$BJ,;RNL(B
      real(DP) :: xyzf_QMixBlAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !$B6E=L@.J,$N:.9gHf(B
      integer  :: s

      
      !$B29EY(B, $B05NO(B, $B:.9gHf$NA4NL$r5a$a$k(B
      !$B>qMp@.J,$HJ?6Q@.J,$NB-$7;;(B
      xyz_TempBlAll   = ( xyz_PTempBl + xyz_PTempBZ ) * ( xyz_ExnerBl + xyz_ExnerBZ )
      xyzf_QMixBlAll  = xyzf_QMixBZ + xyzf_QMixBl
      xyzf_LatentHeat = 0.0d0
      
      !$B:n6HG[Ns$N=i4|2=(B. $B5$BN$N$_MxMQ(B
      !
      do s = 1, GasNum
        xyzf_QMixPerMolWt(:,:,:,s) = xyzf_QMixBl(:,:,:,IdxG(s)) / MolWtWet(IdxG(s))
      end do
      
      !$B29EY$N8z2L(B
      !
      xyz_EffPTemp = xyz_PTempBl / xyz_PTempBZ 
    
      !$BJ,;RNL8z2L(B + $B0z$-$E$j$N8z2L(B
      !
      xyz_EffMolWt =                                        &
        & + sum( xyzf_QMixPerMolWt, 4 )                     &
        &    / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt )    &
        & - sum(xyzf_QMixBl, 4) / ( 1.0d0 + xyz_QMixBZ )
      
      !$B>x5$$,>xH/$9$k>l9g$N@xG.$r7W;;(B
      !  $BJ,;RNL$NItJ,$O$$$D$G$b8z$/$,@xG.$OK0OB$7$F$$$J$$$H8z$+$J$$$N$G(B, 
      !  $B1@$N:.9gHf$,%<%m$N;~$K$O(B, $B@xG.$N4sM?$O%<%m$H$J$k$h$&$KD4@a$7$F$$$k(B
      !
      do s = 1, CondNum
        xyzf_LatentHeat(:,:,:,s) =                                                   &
          & xyz_LatentHeat( SpcWetID(IdxCC(s)), xyz_TempBlAll )                      &
          &  * xyzf_QMixBlAll(:,:,:,IdxCG(s))                                        &
          &  * ( 5.0d-1 + sign( 5.0d-1, (xyzf_QMixBlAll(:,:,:,IdxCC(s)) - 1.0d-4) ) )
      end do
      
      xyz_EffHeat = ( sum( xyzf_LatentHeat, 4 ) * xyz_EffMolWtBZ ) &
        &           / ( CpDry * xyz_ExnerBZ ) 
         
      !$BMpN.3H;678?t$N;~4VH/E8<0$NIbNO9`$r7h$a$k(B
      !
      xyz_Buoy =                                                           &
        &  - 3.0d0 * Grav * Cm_Cm_MixLen_MixLen                            &
        &    * xyz_xyr(                                                    &
        &        xyr_dz_xyz(                                               &
        &           xyz_EffHeat                                            &
        &         + xyz_VPTempBZ * ( 1.0d0 + xyz_EffPTemp + xyz_EffMolWt ) &
        &         )                                                        &
        &      )                                                           &
        &    / ( 2.0d0 * xyz_VPTempBZ )
      
      xyz_BuoyT =                                           &
        &  - 3.0d0 * Grav * Cm_Cm_MixLen_MixLen             &
        &    * xyz_xyr(                                     &
        &        xyr_dz_xyz( xyz_PTempBl ) + xyr_DPTempBZDz &
        &      ) / ( 2.0d0 * xyz_PTempBZ )

      xyz_BuoyM = xyz_Buoy - xyz_BuoyT


    end subroutine Turbulence_kw1978_BuoyancyKm
    

  end subroutine Turbulence_KW1978_forcing
  
!!!------------------------------------------------------------------------!!!
  
  subroutine turbulence_kw1978_output
    
    ! $B%b%8%e!<%k(B
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

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (x)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtTurbD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (x) (Density)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (y)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtTurbD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (y) (density)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (z)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtTurbD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of velocity (z) (density)', &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DCDensDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of cloud density', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DCDensDtTurbD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of cloud density (density)', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtTurb',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtTurbD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Turbulence term of potential temperature (density)', &
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

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtTurb', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Turbulence term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtTurbD', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Turbulence term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio (density)',  &
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
      & longname='Buoyancy term of Km ', &
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

end module Turbulence_kw1978_v2
