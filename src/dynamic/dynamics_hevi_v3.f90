!= Module DynamicsHEVI
!
! Authors::   $B?y;39L0lO/(B(SUGIYAMA Ko-ichiro)
! Version::   $Id: dynamics_hevi_v3.f90,v 1.4 2014/07/08 00:55:26 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module DynamicsHEVI_v3
  !
  ! $BNO3X%3%"(B. 
  !   $B%?%$%`%9%W%j%C%HK!$rMxMQ(B. $B2;GH%b!<%I$H$=$l0J30$rJL!9$N;~4V9o$_$G2r$/(B. 
  !   $BC;$$;~4V%9%F%C%W$N7W;;$K$O(B HE-VI $BK!$rMxMQ(B.
  !
  ! v3: v1 $B$H(B v2 $B$r@0M}$9$k$3$H$G(B, $B7W;;J}K!$r@Z$jBX$($i$l$k$h$&$KD4@0$7$?HG(B
  !
  ! Note: 
  !  * $B%(%/%9%J!<4X?t$N6u4VJ}8~$NN%;62=$O(B 2 $B<!@:EY$G$"$k$?$a(B, $B5$0579EY(B
  !    $BNO9`$N7W;;%W%m%0%i%`$K$*$$$F(B differentiate_center4 $B%b%8%e!<%k$r(B
  !    $B;XDj$9$k$3$H$O$G$-$J$$$N$GCm0U(B.

  !$B%b%8%e!<%kFI$_9~$_(B
  use dc_types,   only : DP

  !$B0EL[$N7?@k8@6X;_(B
  implicit none

  !$BB0@-$N;XDj(B
  private

  !$BJQ?t$NDj5A(B
  character(*), parameter:: module_name = 'dynamics_hevi'
                                           ! $B%b%8%e!<%k$NL>>N(B.
                                           ! Module name
  real(DP), save, private :: FactorBuoyTemp    = 1.0d0
                                           !$BIbNO(B ($B29EY$N4sM?(B) $B$NM-L5(B
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorBuoyMolWt   = 1.0d0
                                           !$BIbNO(B ($BJ,;RNL8z2L(B) $B$NM-L5(B
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorBuoyLoading = 1.0d0
                                           !$BIbNO(B ($B2Y=E8z2L(B) $B$NM-L5(B
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorDExnerDtAdv    = 0.0d0
                                           !$B%(%/%9%J!<4X?t$N0\N.$NM-L5(B
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorDExnerDtExpnd  = 0.0d0
                                           !$B%(%/%9%J!<4X?t$NKDD%9`$NM-L5(B
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorExnerFall  = 0.0d0
                                           !$BMn2<$KH<$&2>290LJQ2=$KH<$&%(%/%9%J!<4X?t$N(Btendency
                                           !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  integer, save, private  :: IDAdvection               = 0
                                           ! $B0\N.7W;;J}K!(B
  integer, parameter      :: IDAdvectionCenter4_std    = 1
                                           ! 1: 4 $B<!Cf1{:9J,(B ($BHyJ,J?6Q%b%8%e!<%kMxMQ(B)
  integer, parameter      :: IDAdvectionCenter4_2D_dry = 2
                                           ! 2: 4 $B<!Cf1{:9J,(B (2D $BHG7W;;%k!<%A%s(B)
  integer, parameter      :: IDAdvectionCenter4_2D     = 3
                                           ! 2: 4 $B<!Cf1{:9J,(B (2D $BHG7W;;%k!<%A%s(B)
  integer, parameter      :: IDAdvectionCenter4_3D_dry = 4
                                           ! 3: 4 $B<!Cf1{:9J,(B (3D $BHG7W;;%k!<%A%s(B)
  integer, parameter      :: IDAdvectionCenter4_3D      = 5
                                           ! 3: 4 $B<!Cf1{:9J,(B (3D $BHG7W;;%k!<%A%s(B)
  integer, save, private  :: IDAcousticmode            = 0 
                                           ! $BC;$$;~4V%9%F%C%W$N7W;;J}K!(B
  integer, parameter      :: IDAcousticmode_std        = 1 
                                           ! 1: $BHyJ,J?6Q%b%8%e!<%k$NMxMQ(B. 
  integer, parameter      :: IDAcousticmode_2D         = 2
                                           ! 2: 2D $BHG7W;;%k!<%A%s(B
  integer, parameter      :: IDAcousticmode_3D         = 3
                                           ! 3: 3D $BHG7W;;%k!<%A%s(B
  
  ! public 
  !
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Dynamics_Init
    !
    ! $BNO3X%3%"(B 3D $BHG$N=i4|2=%k!<%A%s(B
    !

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,   only : DP, STRING
    use dc_iounit,  only : FileOpen
    use dc_message, only : MessageNotify
    use gridset,    only : FlagCalc3D, FlagCalcMoist
    use namelist_util, &
      &             only : namelist_filename
    use advection_center4_std, &
      &             only : advection_center4_std_init
    use advection_center4_2D, &
      &             only : advection_center4_2D_init
    use advection_center4_3D, &
      &             only : advection_center4_3D_init
    use acousticmode_std, &
      &             only : acousticmode_std_init
    use acousticmode_2D, &
      &             only : acousticmode_2D_init
    use acousticmode_3D, &
      &             only : acousticmode_3D_init

    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none
    
    ! $BJQ?t$NDj5A(B
    !
    real(DP)          :: AlphaSound = 5.0d-2    !$B2;GH8:?j9`$N78?t(B ($B5$>]D#?tCMM=Js2]Js9p!&JL:}(B49 $B$h$j(B)
    real(DP)          :: AlphaNDiff = 1.0d-3    !4$B<!$N?tCM3H;6$N78?t(B. CReSS $B%^%K%e%"%k$h$j(B
    real(DP)          :: NDiffRatio = 1.0d0     !$BB.EY$KBP$9$kG4@-$r>e$2$k>l9g$O?t;z$r(B 1 $B0J>e$K$9$k(B. 
    integer           :: unit                   !$BAuCVHV9f(B
    character(STRING) :: FlagAdvection = ""     !$B0\N.7W;;J}K!$NA*Br(B
    character(STRING) :: FlagAcousticmode = ""  !$B2;GH%b!<%I$N7W;;J}K!$NA*Br(B    

    
    !-------------------------------------------------------------------
    ! Namelist $B$+$i>pJs$r<hF@$9$k(B
    !
    NAMELIST /Dynamics_nml/                                 &
      & AlphaSound, AlphaNDiff, NDiffRatio,                 &
      & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
      & FactorDExnerDtAdv, FactorDExnerDtExpnd, FactorExnerFall,  &
      & FlagAdvection, FlagAcousticmode

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=dynamics_nml)
    close(unit)
    
    !-------------------------------------------------------------------
    ! $B%U%i%0=hM}(B
    !
    if ( FlagAdvection == "Center4" .OR. FlagAdvection == "" ) then 
      !
      ! 4 $B<!@:EYCf?4:9J,(B. $BHyJ,J?6Q%b%8%e!<%k$r;H$o$J$$>l9g(B. ($B%G%U%)%k%H(B)
      !
      if ( FlagCalc3D ) then 
!        if ( FlagCalcTracer .OR. FlagCalcMoist ) then 
        if ( FlagCalcMoist ) then 
          IDAdvection = IDAdvectionCenter4_3D
        else
          IDAdvection = IDAdvectionCenter4_3D_dry
        end if
      else
!        if ( FlagCalcTracer .OR. FlagCalcMoist ) then 
        if ( FlagCalcMoist ) then 
          IDAdvection = IDAdvectionCenter4_2D
        else
          IDAdvection = IDAdvectionCenter4_2D_dry
        end if
      end if
    else 
      !
      ! 4 $B<!@:EYCf?4:9J,(B with $BHyJ,J?6Q%b%8%e!<%k(B
      !
      IDAdvection = IDAdvectionCenter4_std
    end if
    
    if ( FlagAcousticmode == "Center2" .OR. FlagAcousticmode == "" ) then 
      !
      ! $BHyJ,J?6Q%b%8%e!<%k$r;H$o$J$$>l9g(B. ($B%G%U%)%k%H(B)
      !
      if ( FlagCalc3D ) then 
        IDAcousticmode = IDAcousticmode_3D
      else
        IDAcousticmode = IDAcousticmode_2D
      end if
    else
      !
      ! with $BHyJ,J?6Q%b%8%e!<%k(B
      !
      IDAcousticmode = IDAcousticmode_std
    end if

    !-------------------------------------------------------------------
    ! $B0\N.7W;;MQ$N7W;;%b%8%e!<%k$N=i4|2=(B:
    !
    select case ( IDAdvection )

    case ( IDAdvectionCenter4_std )
      call advection_center4_std_init( AlphaNDiff, NDiffRatio )

    case ( IDAdvectionCenter4_3D_dry, IDAdvectionCenter4_3D )
      call advection_center4_3D_init( AlphaNDiff, NDiffRatio )
      
    case ( IDAdvectionCenter4_2D_dry, IDAdvectionCenter4_2D )
      call advection_center4_2D_init( AlphaNDiff, NDiffRatio )
      
    end select

    !-------------------------------------------------------------------
    ! $B2;GH%b!<%I$N7W;;%b%8%e!<%k$N=i4|2=(B:
    !
    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      call acousticmode_std_init( AlphaSound )

    case ( IDAcousticmode_3D )
      call acousticmode_3D_init( AlphaSound )

    case ( IDAcousticmode_2D )
      call acousticmode_2D_init( AlphaSound )
      
    end select
    
    !-------------------------------------------------------------------
    ! tendency $B$N=PNO(B
    !
    call Dynamics_Tendency_output

    !-------------------------------------------------------------------
    ! $B=PNO(B
    !
    call MessageNotify( "M", module_name, "AlphaSound = %f", d=(/AlphaSound/) )
    call MessageNotify( "M", module_name, "AlphaNDiff = %f", d=(/AlphaNDiff/) )
    
    call MessageNotify( "M", module_name, "FactorBuoyTemp   = %f", d=(/FactorBuoyTemp/) )
    call MessageNotify( "M", module_name, "FactorBuoyMolWt  = %f", d=(/FactorBuoyMolWt/) )
    call MessageNotify( "M", module_name, "FactorBuoyLoading= %f", d=(/FactorBuoyLoading/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtAdv   = %f", d=(/FactorDExnerDtAdv/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtExpnd = %f", d=(/FactorDExnerDtExpnd/) )
    
    call MessageNotify( "M", module_name, "IDAdvection = %d", i=(/IDAdvection/) )
    call MessageNotify( "M", module_name, "IDAcousticmode = %d", i=(/IDAcousticmode/) )
 
    
  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Long_forcing(   &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,        & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
    & xyz_PTempBl, xyz_PTempNl,       & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
    & xyz_KmBl,    xyz_KmNl,          & ! (in)
    & pyz_DVelXDtNl,                  & ! (inout)
    & xqz_DVelYDtNl,                  & ! (inout)
    & xyr_DVelZDtNl,                  & ! (inout)
    & xyz_DPTempDtNl,                 & ! (inout)
    & xyz_DExnerDtNl,                 & ! (inout)
    & xyzf_DQMixDtNl,                 & ! (inout)
    & xyz_DKmDtNl                     & ! (inout)
    & )
    !
    ! $BNO3X%3%"(B: $BD9$$;~4V%9%F%C%W$GI>2A$9$k9`$N7W;;(B.
    !

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
                                                    !$BG[Ns%5%$%:(B
      &                   nx, ny, nz,              &!$BJ*M}NN0h$NG[Ns%5%$%:(B
      &                   FlagCalcMoist !$B%U%i%0(B
    use composition, &
      &             only : SpcWetSymbol
    use timeset,    only : TimeN
    use gtool_historyauto, &
      &             only : HistoryAutoPut
    use advection_center4_std, &
      &             only : advection_center4_std_main
    use advection_center4_2D, &
      &             only : advection_center4_2D_dry, advection_center4_2D_tracer
    use advection_center4_3D, &
      &             only : advection_center4_3D_dry, advection_center4_3D_tracer
    use DExnerDt,   only : xyz_DExnerDt_xyzf

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BJQ?t$NDj5A(B
    !
    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(inout) :: xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)                :: pyz_DVelXDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_DVelYDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_DVelZDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DKmDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DExnerDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_DPTempDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_QMixAdv(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: xyzf_QMixFall(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: pyz_VelXnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_VelYnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_VelZnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_KmNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_ExnerNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_ExnerFall(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_PTempNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_QMixNDiff(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    
    integer                 :: f
    
    
    !--------------------------------------------------------------------
    ! $B0\N.9`$N7W;;(B
    !
    select case ( IDAdvection )
      
    case ( IDAdvectionCenter4_std ) 

      call advection_center4_std_main(     &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyzf_QMixAdv, xyzf_QMixNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )
      
    case ( IDAdvectionCenter4_3D_dry ) 
      
      call advection_center4_3D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

    case ( IDAdvectionCenter4_3D ) 
      
      call advection_center4_3D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,   xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xqz_DVelYDtAdv,  xqz_VelYnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

      call advection_center4_3D_tracer(       &
        & pyz_VelXNl, xqz_VelYNl, xyr_VelZNl, & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,          & ! (in)
        & xyzf_QMixAdv, xyzf_QMixNDiff        & !(out)
        & )
      
    case ( IDAdvectionCenter4_2D_dry ) 
      
      ! Y $BJ}8~$O(B tendency $B$ONm(B. 
      !
      xqz_DVelYDtAdv   = 0.0d0
      xqz_VelYnDiff = 0.0d0

      call advection_center4_2D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

    case ( IDAdvectionCenter4_2D ) 
      
      ! Y $BJ}8~$O(B tendency $B$ONm(B. 
      !
      xqz_DVelYDtAdv   = 0.0d0
      xqz_VelYnDiff = 0.0d0
      
      call advection_center4_2D_dry(       &
        & pyz_VelXBl,   pyz_VelXNl,        & ! (in)
        & xyr_VelZBl,   xyr_VelZNl,        & ! (in)
        & xyz_PTempBl,  xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl,  xyz_ExnerNl,       & ! (in)
        & xyz_KmBl,     xyz_KmNl,          & ! (in)
        & pyz_DVelXDtAdv,  pyz_VelXnDiff,     & !(out)
        & xyr_DVelZDtAdv,  xyr_VelZnDiff,     & !(out)
        & xyz_DPTempDtAdv, xyz_PTempNDiff,    & !(out)
        & xyz_DExnerDtAdv, xyz_ExnerNDiff,    & !(out)
        & xyz_DKmDtAdv,    xyz_KmNDiff        & !(out)
        & )

      call advection_center4_2D_tracer(    &
        & pyz_VelXNl, xyr_VelZNl,          & ! (in)
        & xyzf_QMixBl,  xyzf_QMixNl,       & ! (in)
        & xyzf_QMixAdv, xyzf_QMixNDiff     & !(out)
        & )

    end select

    !--------------------------------------------------------------------
    ! tendency $B$N99?7(B
    !

    ! $B3H;678?t(B 
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_KmNDiff + xyz_DKmDtAdv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_DKmDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_KmNDiff(1:nx,1:ny,1:nz) )
       
    ! $B290L(B
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_PTempNDiff + xyz_DPTempDtAdv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',  xyz_DPTempDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DPTempDtDiff', xyz_PTempNDiff(1:nx,1:ny,1:nz) )

    ! $B:.9gHf(B
    !
    call QMixfall  ! $BMn2<9`(B

    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_QMixNDiff + xyzf_QMixAdv + xyzf_QMixFall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_QMixAdv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_QMixNDiff(1:nx,1:ny,1:nz,f))
    end do

    ! $B%(%/%9%J!<4X?t(B
    !
    xyz_ExnerFall   = xyz_DExnerDt_xyzf( xyzf_QMixFall ) * FactorExnerFall
    xyz_ExnerNDiff  = xyz_ExnerNDiff   * FactorDExnerDtAdv
    xyz_DExnerDtAdv = xyz_DExnerDtAdv  * FactorDExnerDtAdv
    xyz_DExnerDtNl  = xyz_DExnerDtNl + ( xyz_ExnerFall+ xyz_ExnerNDiff + xyz_DExnerDtAdv )
   
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_DExnerDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_ExnerNDiff(1:nx,1:ny,1:nz) )
    
    ! VelX 
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_VelXNDiff + pyz_DVelXDtAdv )
    
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_DVelXDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_VelXnDiff(1:nx,1:ny,1:nz) )
    
    ! VelY
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_VelYNDiff + xqz_DVelYDtAdv )
    
    call HistoryAutoPut(TimeN, 'DVelYDtAdv',  xqz_DVelYDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_VelYnDiff(1:nx,1:ny,1:nz) )
    
    ! VelZ
    !
    call BuoyancyLong_xyr       ! $BIbNO9`(B
    
    xyr_DVelZDtNl =                            &
      &   xyr_DVelZDtNl                        &
      & + (                                    &
      &      xyr_VelZnDiff                     & 
      &    + xyr_DVelZDtAdv                       &
      &    + (                                 &
      &        + xyr_BuoyM * FactorBuoyMolWt   &
      &        + xyr_BuoyD * FactorBuoyLoading &
      &        + xyr_BuoyT * FactorBuoyTemp    &
      &      )                                 &
      &   )
    
    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_DVelZDtAdv(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_VelZnDiff(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyT', xyr_BuoyT(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyM', xyr_BuoyM(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyD', xyr_BuoyD(1:nx,1:ny,1:nz) )
    
  contains
    
    subroutine BuoyancyLong_xyr
      !
      ! $BIbNO9`$N7W;;(B
      !
    
      use composition,  only: GasNum,       &! 
        &                     IdxG,         &!
        &                     MolWtWet       ! $B<>=a@.J,$NJ,;RNL(B
      use constants,   only : MolWtDry,     &! $B4%Ag@.J,$NJ,;RNL(B
        &                     Grav           ! $B=ENO2CB.EY(B
      use basicset,    only : xyz_PTempBZ,  &! $B4pK\>l$N290L(B
        &                     xyr_QMixBZ,   &! $B4pK\>l$N:.9gHf(B
        &                     xyr_QMixBZPerMolWt

      ! $B0EL[$N7?@k8@6X;_(B
      !
      implicit none
      
      ! $B:n6HJQ?t(B
      !
      real(DP)    :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)    :: xyz_QMixPerMolWtSum(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_QMixNlSumGas(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_QMixNlSumCnd(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)    :: xyz_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
      integer     :: i, j, k, f, n
      
      ! Buoyancy due to temperature disturbunce
      !
      xyz_BuoyT = Grav * xyz_PTempNl / xyz_PTempBZ
      
      ! $B3J;R0LCV$NJQ49(B
      !
      do k = kmin, kmax-1
        xyr_BuoyT(:,:,k) =                           &
          & (                                        &
          &   xyz_BuoyT(:,:,k+1) + xyz_BuoyT(:,:,k)  &
          & ) * 5.0d-1 
      end do

      !$B7jKd$a(B
      xyr_BuoyT(:,:,kmax) = 0.0d0

      ! $B4%Ag$N>l9g$O(B BuoyD, BuoyM $B$ONm(B. 
      !
      if ( .NOT. FlagCalcMoist ) then 

        xyr_BuoyM = 0.0d0
        xyr_BuoyD = 00d0

        ! $B%5%V%k!<%A%s$rH4$1$k(B
        ! 
        return
      end if
           
      ! Buoyancy due to molecular weight
      !
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do

      xyz_QMixPerMolWtSum = sum(xyzf_QMixPerMolWt, 4) 
      xyz_QMixNlSumGas    = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax

            xyr_BuoyM(i,j,k) =                                         &
              & + Grav                                                 &
              &   * (                                                  &
              &         xyz_QMixPerMolWtSum(i,j,k+1)                   &
              &       + xyz_QMixPerMolWtSum(i,j,k)                     &
              &     ) * 5.0d-1                                         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) )   &
              & - Grav                                                 &
              &   * (                                                  &
              &         xyz_QMixNlSumGas(i,j,k+1)                      &
              &       + xyz_QMixNlSumGas(i,j,k)                        &
              &     ) * 5.0d-1                                         &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) ) 
            
          end do
        end do
      end do

      ! $B7jKd$a(B
      !
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      xyz_QMixNlSumCnd = sum( xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4 ) 
      
      do k = kmin, kmax-1
        do j = jmin, jmax
          do i = imin, imax
            
            xyr_BuoyD(i,j,k) =                       &
              & - Grav                               &
              &   * (                                &
              &         xyz_QMixNlSumCnd(i,j,k+1)    &
              &       + xyz_QMixNlSumCnd(i,j,k)      &
              &     ) * 5.0d-1                       &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do

      ! $B7jKd$a(B
      !
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr


    subroutine QmixFall
      !
      ! $B1+N3$NMn2<$K$h$k0\N.$r5a$a$k(B. 
      ! 

      ! $B%b%8%e!<%k8F$S=P$7(B
      !
      use average,          only : xyr_xyz
      use differentiate_center4, &
        &                   only : xyz_dz_xyr
      use constants,        only : FactorJ
      use composition,      only : IdxR, RainNum
      use basicset,         only : xyzf_QMixBZ,  &!$B4pK\>l$N:.9gHf(B
        &                          xyz_DensBZ     !$B4pK\>l$NL)EY(B
      
      !$B0EL[$N7?@k8@6X;_(B
      !
      implicit none
      
      !$BJQ?tDj5A(B
      !
      real(DP)  :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !$B>x5$:.9gHf(B($B>qMp(B + $BJ?6Q>l(B)
      real(DP)  :: xyrf_QMixFallFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !$B1+N3$NMn2<8z2L(B
      real(DP)  :: xyz_VelZFall(imin:imax,jmin:jmax,kmin:kmax)
                                                 !$B1+N3Mn2<B.EY(B
      integer   :: s, iR

      ! $B=i4|2=(B
      !
      xyzf_QMixAll = max( 0.0d0, xyzf_QMixBl + xyzf_QMixBZ )
      xyrf_QMixFallFlux = 0.0d0
      xyzf_QMixFall = 0.0d0
      xyz_VelZFall = 0.0d0

      ! $BMn2<9`$N7W;;(B
      !
      do s = 1, RainNum

        iR = IdxR(s)

        !$B1+N3=*C<B.EY(B
        xyz_VelZFall = - 12.2d0 * FactorJ * ( xyzf_QMixAll(:,:,:,iR) ** 0.125d0 )
        
        ! $B%U%i%C%/%9$N7W;;(B
        !
        xyrf_QMixFallFlux(:,:,:,iR) =                                &
          &  xyr_xyz (                                               &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall    &
          &  )

        ! $B>eC<$N%U%i%C%/%9$O%<%m(B
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! $B1+N3Mn2<$K$h$k;~4VJQ2=(B 
        !        
        xyzf_QMixFall(:,:,:,iR) =                                      &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
        call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtFall', &
          & xyzf_QMixFall(1:nx,1:ny,1:nz,iR))
        
      end do
            
    end subroutine QmixFall

  end subroutine Dynamics_Long_forcing
  
  
!!!---------------------------------------------------------------!!!
  
  subroutine Dynamics_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xqz_VelYNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xqz_DVelYDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! $BC;$$;~4V%9%F%C%W$GI>2A$9$k9`$N7W;;(B.
    !
    

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                             !$BG[Ns%5%$%:(B
      &                  nx, ny, nz          !$BJ*M}NN0h$NBg$-$5(B
    use timeset,  only : DelTimeShort, TimeN
    use setmargin,only : SetMargin_xyzf, SetMargin_xyz, &
      &                  SetMargin_pyz, SetMargin_xqz, SetMargin_xyr
    use gtool_historyauto, &
      &           only : HistoryAutoPut
    use acousticmode_std, &
      &           only : acousticmode_std_exp, acousticmode_std_imp
    use acousticmode_2D, &
      &           only : acousticmode_2D_exp, acousticmode_2D_imp
    use acousticmode_3D, &
      &           only : acousticmode_3D_exp, acousticmode_3D_imp
    use constants,only : CvDry,    &! $B4%Ag@.J,$NHfG.(B
      &                  GasRDry
    
    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BJQ?t$NDj5A(B
    !
    real(DP), intent(in)    :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(out)   :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)
    
    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Expnd(imin:imax,jmin:jmax,kmin:kmax)

    ! subroutine $B$N0z?t$H$J$kG[Ns$O(B, $BJ*M}NN0h$NBg$-$5$K(B.
    !
    real(DP) :: pyz_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: xqz_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: xyr_PGrad(1:nx,1:ny,1:nz)
    real(DP) :: pyz_SWF(1:nx,1:ny,1:nz)
    real(DP) :: xqz_SWF(1:nx,1:ny,1:nz)
    real(DP) :: xyr_SWF(1:nx,1:ny,1:nz)

    !------------------------------------------------------------
    ! $B=i4|2=(B
    !
    pyz_DVelXDtNs = 0.0d0
    xqz_DVelYDtNs = 0.0d0
    xyr_DVelZDtNs = 0.0d0
    
    !------------------------------------------------------------
    ! $B?eJ?J}8~(B: $BM[2rK!(B
    !
    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      call acousticmode_std_exp(              &
        & pyz_VelXNs, xqz_VelYNs, xyr_VelZNs, & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(OUT)
        & pyz_PGrad, pyz_SWF,                 & !(OUT)
        & xqz_PGrad, xqz_SWF                  & !(OUT)
        & )
      
    case ( IDAcousticmode_3D )
      call acousticmode_3d_exp(               &
        & pyz_VelXNs, xqz_VelYNs, xyr_VelZNs, & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(IN)
        & pyz_PGrad, pyz_SWF,                 & !(OUT)
        & xqz_PGrad, xqz_SWF                  & !(OUT)
        & )
      
    case ( IDAcousticmode_2D )
      xqz_PGrad = 0.0d0
      xqz_SWF   = 0.0d0

      call acousticmode_2d_exp(               &
        & pyz_VelXNs, xyr_VelZNs,             & !(IN)
        & xyz_ExnerNs,                        & !(IN)
        & xyz_VelDivNs,                       & !(IN)
        & pyz_PGrad, pyz_SWF                  & !(OUT)
        & )
      
    end select

    ! tendency 
    !
    pyz_DVelXDtNs(1:nx,1:ny,1:nz) = pyz_PGrad + pyz_SWF
    xqz_DVelYDtNs(1:nx,1:ny,1:nz) = xqz_PGrad + xqz_SWF
    
    ! $BCM$NJ]4I(B
    !
!    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'VelDiv', pyz_VelXNs(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelDiv', xyr_VelZNs(1:nx,1:ny,1:nz))

    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad )
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF   )
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_PGrad )
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_SWF   )

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * ( pyz_DVelXDtNl + pyz_DVelXDtNs )
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * ( xqz_DVelYDtNl + xqz_DVelYDtNs )

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)
    

    !------------------------------------------------------------
    ! $B1tD>J}8~(B: $B1"2rK!(B
    !  $B6E7k9`$rC;$$;~4V%9%F%C%W$GI>2A$9$k>l9g$r9M$($F(B, $BD9$$;~4V%9%F%C%W$H(B
    !  $BC;$$;~4V%9%F%C%W$GI>2A$5$l$?(B tendency $B$rB-$79g$o$;$?$b$N$r0z?t$K(B. 
    !

    ! $BC;$$;~4V%9%F%C%W$GI>2A$9$k05NO$N<0$N(B tendency
    !
    xyz_Expnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &         * FactorDExnerDtExpnd

    ! tendency $B$N9g7W(B
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_Expnd

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DExnerDtExpnd', xyz_Expnd(1:nx,1:ny,1:nz))   

    select case ( IDAcousticmode )
      
    case ( IDAcousticmode_std )
      
      call acousticmode_std_imp(                             &
        & pyz_VelXAs, xqz_VelYAs, xyr_VelZNs, xyz_VelDivNs,  & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    case ( IDAcousticmode_3D )
      
      call acousticmode_3D_imp(                              &
        & pyz_VelXAs, xqz_VelYAs, xyr_VelZNs, xyz_VelDivNs,  & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    case ( IDAcousticmode_2D )

      call acousticmode_2D_imp(                              &
        & pyz_VelXAs, xyr_VelZNs, xyz_VelDivNs,              & !(IN)
        & xyz_ExnerNs,                                       & !(IN)
        & xyr_DVelZDtNl, xyz_DExnerDtNl, xyz_DExnerDtNs,     & !(IN)
        & xyz_ExnerAs,                                       & !(OUT)
        & xyr_PGrad, xyr_SWF                                 & !(OUT)
        & )
      
    end select
    
    ! tendency
    !
    xyr_DVelZDtNs(1:nx,1:ny,1:nz) = xyr_PGrad + xyr_SWF
    
    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * ( xyr_DVelZDtNl + xyr_DVelZDtNs )

    ! Set Margin
    !
    call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
    call SetMargin_xyr( xyr_VelZAs )  ! (inout)

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad',  xyr_PGrad )
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',    xyr_SWF   )
    
  end subroutine Dynamics_Short_forcing
  
!!!------------------------------------------------------------------------!!!

  subroutine Dynamics_Tendency_Output
    !
    ! $B=PNO$N@_Dj(B
    ! 

    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : ncmax             ! $BJ*<A?t(B
    use composition,       only : SpcWetSymbol, &
      &                           IdxR, RainNum
    
    implicit none
   
    integer :: f, iR

    call HistoryAutoAddVariable(                              &
      & varname='DPTempDtAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of potential temperature',   &
      & units='K.s-1',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                                      &
      & varname='DPTempDtDiff',                                          &
      & dims=(/'x','y','z','t'/),                                     &
      & longname='Numerical diffusion term of potential temperature', &
      & units='K.s-1',                                                &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DExnerDtAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of exner function',          &
      & units='s-1',                                          &
      & xtype='double')
    
    call HistoryAutoAddVariable(                               &
      & varname='DExnerDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                              &
      & longname='Numerical diffusion term of exner function', &
      & units='s-1',                                           & 
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DExnerDtExpnd',                                 &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Expanding term of exner function',          &
      & units='s-1',                                          &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='CDensAdv',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of cloud density',           &
      & units='kg.m-3.s-1',                                   & 
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='CDensDiff',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of cloud density', &
      & units='kg.m-3.s-1',                                   &
      & xtype='double')

    do f = 1, ncmax
      call HistoryAutoAddVariable(                            &
        & varname='D'//trim(SpcWetSymbol(f))//'DtAdv',              &
        & dims=(/'x','y','z','t'/),                           &
        & longname='Advection term of '                       &
        &           //trim(SpcWetSymbol(f))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                &
        & xtype='double')
      
      call HistoryAutoAddVariable(                            &
        & varname='D'//trim(SpcWetSymbol(f))//'DtDiff',             &
        & dims=(/'x','y','z','t'/),                           &
        & longname='Diffusion term of '                       &
        &           //trim(SpcWetSymbol(f))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                &
        & xtype='double')
    end do

    do f = 1, RainNum
      iR = IdxR(f)
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtFall',             &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Fall term of '                             &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio', &
        & units='kg.kg-1.s-1',                                 &
        & xtype='double')

    end do

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (x)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of velocity (x)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (x)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelXDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (x)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (y)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of velocity (y)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (y)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelYDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (y)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtAdv',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection term of velocity (z)',            &
      & units='m.s-2',                                        &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtDiff',                                   &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Numerical diffusion term of Velocity (z)',  &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyT',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (Temperature)',                    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyM',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (MolWt)',                          &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtBuoyD',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Buoyancy (Drag)',                           &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtPGrad',                                  &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Pressure gradient term of velocity (z)',    &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DVelZDtSWF',                                    &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Filter for acoustic mode (z)',              &
      & units='m.s-2',                                        &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='VelDiv',                                     &
      & dims=(/'x','y','z','t'/),                             &
      & longname='velocity divergence',                       &
      & units='s-2',                                          &
      & xtype='double')

    call HistoryAutoAddVariable(                              &
      & varname='DKmDtAdv',                                      &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Advection of Km',                           &
      & units='s-1',                                          &
      & xtype='double')
    
    call HistoryAutoAddVariable(                              &
      & varname='DKmDtDiff',                                     &
      & dims=(/'x','y','z','t'/),                             &
      & longname='Diffusion term of Km',                      &
      & units='s-1',                                          &
      & xtype='double')
    
!    call HistoryAutoAddVariable(  &
!      & varname='D', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix D',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='E', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix E',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='F', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix F',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNs',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENl', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNl',  &
!      & units='s-1',    &
!      & xtype='double')

  end subroutine Dynamics_Tendency_Output
  
  
end module DynamicsHEVI_v3

