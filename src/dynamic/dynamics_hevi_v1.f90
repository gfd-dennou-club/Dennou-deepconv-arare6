!= Module DynamicsHEVI
!
! Authors::   $B?y;39L0lO/(B(SUGIYAMA Ko-ichiro), ODAKA Masatsugu 
! Version::   $Id: dynamicshevi.f90,v 1.21 2014/05/28 15:27:42 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module DynamicsHEVI
  !
  ! $BNO3X%3%"(B. 
  !   $B%?%$%`%9%W%j%C%HK!$rMxMQ(B. $B2;GH%b!<%I$H$=$l0J30$rJL!9$N;~4V9o$_$G2r$/(B. 
  !   $BC;$$;~4V%9%F%C%W$N7W;;$K$O(B HE-VI $BK!$rMxMQ(B.
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
  
  real(DP), save, private  :: beta  = 5.0d-1     !$B%/%i%s%/%K%3%k%=%sK!$J$i(B 0.5
                                                 !$B40A41"2rK!$J$i(B 1
  integer, save, private   :: N = 10             !$B78?t9TNs(B/$B2~9TNs$N<!?t(B, $B@09g@#K!(B
  integer, save, private   :: M = 10             !$BJ}Dx<0$NAH?t(B
  integer, save, private   :: NUD = 1            !$B78?t9TNs$N>e;03QItJ,$NBSI}(B
  integer, save, private   :: NLD = 1            !$B78?t9TNs$N2<;03QItJ,$NBSI}(B
  integer, save, private   :: NAL = 1            !LU $BJ,2r$N7k2L(B L $B$N@09g@#K!(B
  integer, save, private   :: NA = 3             !NUD + NLD + 1

  real(DP), allocatable, save, private :: xyz_F1BZ(:,:,:)
                                                 !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B
  real(DP), allocatable, save, private :: xyr_F2BZ(:,:,:)
                                                 !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B

  real(DP), allocatable, save, private :: A(:)   !$B78?t9TNs$NBP3Q@.J,(B
  real(DP), allocatable, save, private :: B(:)   !$B78?t9TNs$N>e;03QItJ,(B
  real(DP), allocatable, save, private :: C(:)   !$B78?t9TNs$N2<;03QItJ,(B
  real(DP), allocatable, save, private :: AL1(:) !LU $BJ,2r$N7k2L(B L (1 $B<!85G[Ns(B)
  integer, allocatable, save, private  :: IP(:)  !$BItJ,%T%\%C%H8r49$N>pJs$r3JG<(B

  real(DP), save, private :: AlphaH = 0.0d0      !$B2;GH8:?j9`$N8:?j78?t(B
  real(DP), save, private :: AlphaV = 0.0d0      !$B2;GH8:?j9`$N8:?j78?t(B
  real(DP), save, private :: NuHh   = 0.0d0      !$B?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save, private :: NuVh   = 0.0d0      !$B?tCMG4@-$N78?t(B ($B1tD>J}8~(B)
  real(DP), save, private :: NuHm   = 0.0d0      !$B?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save, private :: NuVm   = 0.0d0      !$B?tCMG4@-$N78?t(B ($B1tD>J}8~(B)

  character(*), parameter :: module_name = 'DynamicHEVI'
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
  real(DP), save, private :: FactorDExnerDtAdv    = 1.0d0 
                                                 !$B%(%/%9%J!<4X?t$N0\N.$NM-L5(B
                                                 !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorDExnerDtExpnd  = 1.0d0
                                                 !$B%(%/%9%J!<4X?t$NKDD%9`$NM-L5(B
                                                 !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B. 

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:)   !$B4pK\>l$N1tD>HyJ,(B
!  real(DP), allocatable, save,private :: xyr_DExnerBZDz(:,:,:)   !$B4pK\>l$N1tD>HyJ,(B
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:) !$B4pK\>l$N1tD>HyJ,(B

  !public 
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Init

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,      only : DP
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use gridset,       only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
      &                       FlagCalc3D
    use timeset,       only : DelTimeShort, DelTimeLong
    use axesset,       only : dx, dy, dz        ! $B3J;R4V3V(B
    use namelist_util, only : namelist_filename
    use basicset,      only : xyz_PTempBZ,     &!$B4pK\>l$N290L(B
      &                       xyz_ExnerBZ,     &
      &                       xyzf_QMixBZ
    use differentiate_center4, &
      &                only : xyr_dz_xyz
    
    !$B0EL[$N7?@k8@6X;_(B
    implicit none
    
    real(DP)  :: AlphaSound = 5.0d-2  !$B2;GH8:?j9`$N78?t(B ($B5$>]D#?tCMM=Js2]Js9p!&JL:}(B49 $B$h$j(B)
    real(DP)  :: AlphaNDiff  = 1.0d-3 !4$B<!$N?tCM3H;6$N78?t(B. CReSS $B%^%K%e%"%k$h$j(B
    real(DP)  :: NDiffRatio = 1.0d0   !$BB.EY$KBP$9$kG4@-$r>e$2$k>l9g$O?t;z$r(B 1 $B0J>e$K$9$k(B. 
    integer   :: unit            !$BAuCVHV9f(B
    integer   :: f

    !-------------------------------------------------------------------
    ! Namelist $B$+$i>pJs$r<hF@$9$k(B
    !
    NAMELIST /Dynamics_nml/                                    &
         & AlphaSound, AlphaNDiff, NDiffRatio, beta,           &
         & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
         & FactorDExnerDtAdv, FactorDExnerDtExpnd

    !$B%U%!%$%k%*!<%W%s(B. $B>pJs<hF@(B. 
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=dynamics_nml)
    close(unit)

    !-------------------------------------------------------------------
    ! $B2;GH8:?j9`$N8:?j78?t$r7h$a$k(B
    ! 
    ! $B5$>]D#M=Js2]Js9pJL:}(B 49 p53 $B$K=>$$(B, $B?eJ?$H1tD>$H$rJ,$1$F9M$($k(B. 
    ! $B$J$*(B, 2 $B<!857W;;$N>l9g$K$O(B DelY $B$K0MB8$7$J$$$h$&$K$9$k(B. 
    !
    if ( FlagCalc3D ) then 
      AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
    else
      AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
    end if

    !-------------------------------------------------------------------
    ! $B?tCM3H;678?t$r7h$a$k(B
    !
    ! CReSS $B%^%K%e%"%k$N5-=R$K=>$C$F(B NuH, NuV $B$r7h$a$k(B.
    ! $B1?F0NL$HG.$KBP$9$k?tCM3H;6$NBg$-$5$rJQ$($i$l$k$h$&$K(B NDiffRatio $B$r>h$8$F$$$k(B.
    ! 
    if ( FlagCalc3D ) then     
      NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    else
      NuHh = AlphaNDiff * ( dx ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    end if
    NuVh = AlphaNDiff * ( dz ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    NuHm = NuHh * NDiffRatio
    NuVm = NuVh * NDiffRatio

    !-------------------------------------------------------------------
    ! $B=PNO(B
    !
    call MessageNotify( "M", module_name, "AlphaH = %f", d=(/AlphaH/) )
    call MessageNotify( "M", module_name, "AlphaV = %f", d=(/AlphaV/) )
    call MessageNotify( "M", module_name, "NuHh = %f",   d=(/NuHh/)   )
    call MessageNotify( "M", module_name, "NuVh = %f",   d=(/NuVh/)   )
    call MessageNotify( "M", module_name, "NuHm = %f",   d=(/NuHm/)   )
    call MessageNotify( "M", module_name, "NuVm = %f",   d=(/NuVm/)   )
    call MessageNotify( "M", module_name, "FactorBuoyTemp = %f",    d=(/FactorBuoyTemp/)    )
    call MessageNotify( "M", module_name, "FactorBuoyMolWt = %f",   d=(/FactorBuoyMolWt/)   )
    call MessageNotify( "M", module_name, "FactorBuoyLoading = %f", d=(/FactorBuoyLoading/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtAdv = %f",    d=(/FactorDExnerDtAdv/)    )
    call MessageNotify( "M", module_name, "FactorDExnerDtExpnd = %f",  d=(/FactorDExnerDtExpnd/)  )

    ! $B1"2rK!$N7W;;@_Dj$N=i4|2=(B
    !
    call DynamicsVI_init()

    ! tendency $B$N=PNO(B
    !
    call Dynamics_Tendency_output

    ! $BG[Ns$NMQ0U(B
    !
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
!    allocate( xyr_DExnerBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
!    xyr_DExnerBZDz = xyr_dz_xyz( xyz_ExnerBZ )
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Long_forcing(        &
    & pyz_VelXBl,  pyz_VelXNl,    & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,    & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,    & ! (in)
    & xyz_PTempBl, xyz_PTempNl,   & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,   & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,   & ! (in)
    & xyz_KmBl,    xyz_KmNl,      & ! (in)
    & pyz_DVelXDtNl,              & ! (inout)
    & xqz_DVelYDtNl,              & ! (inout)
    & xyr_DVelZDtNl,              & ! (inout)
    & xyz_DPTempDtNl,             & ! (inout)
    & xyz_DExnerDtNl,             & ! (inout)
    & xyzf_DQMixDtNl,             & ! (inout)
    & xyz_DKmDtNl                 & ! (inout)
    & )
    !
    ! $BNO3X%3%"(B: $BD9$$;~4V%9%F%C%W$GI>2A$9$k9`$N7W;;(B.
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,  only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x $BJ}8~$NG[Ns%5%$%:(B
      &                     jmin, jmax,       &! y $BJ}8~$NG[Ns%5%$%:(B
      &                     kmin, kmax,       &! z $BJ}8~$NG[Ns%5%$%:(B
      &                     nx, ny, nz,       &! $BJ*M}NN0h$N%5%$%:(B
      &                     ncmax              ! $BJ*<A?t(B
    use timeset,     only : TimeN
    use composition, only : SpcWetSymbol
    use basicset,    only : xyz_PTempBZ,       & ! $B290L$N4pK\>l(B
      &                     xyzf_QMixBZ,       & ! $B:.9gHf$N4pK\>l(B
      &                     xyr_QMixBZ,        &
      &                     xyr_QMixBZPerMolWt
    use average,     only : pyz_xyz, pyz_pqz, pqz_xqz, pyz_pyr, pyr_xyr, &
      &                     xqz_pqz, pqz_pyz, xqz_xyz, xqz_xqr, xqr_xyr, &
      &                     xyr_pyr, pyr_pyz, xyr_xqr, xqr_xqz, xyr_xyz, &
      &                     xyz_pyz, xyz_xqz, xyz_xyr
    use differentiate_center4,                                                                              &
      &              only : pyz_c4dx_xyz => pyz_dx_xyz, xqz_c4dy_xyz => xqz_dy_xyz, xyr_c4dz_xyz => xyr_dz_xyz, &
      &                     xyz_c4dx_pyz => xyz_dx_pyz, pqz_c4dy_pyz => pqz_dy_pyz, pyr_c4dz_pyz => pyr_dz_pyz, &
      &                     pqz_c4dx_xqz => pqz_dx_xqz, xyz_c4dy_xqz => xyz_dy_xqz, xqr_c4dz_xqz => xqr_dz_xqz, &
      &                     pyr_c4dx_xyr => pyr_dx_xyr, xqr_c4dy_xyr => xqr_dy_xyr, xyz_c4dz_xyr => xyz_dz_xyr
    use differentiate_center2,                              &
      &              only : pyz_dx_xyz, xyz_dx_pyz, pyz_dy_pqz, &
      &                     pqz_dy_pyz, pyz_dz_pyr, pyr_dz_pyz, &
      &                     xqz_dx_pqz, pqz_dx_xqz, xqz_dy_xyz, &
      &                     xyz_dy_xqz, xqz_dz_xqr, xqr_dz_xqz, &
      &                     xyr_dx_pyr, pyr_dx_xyr, xyr_dy_xqr, &
      &                     xqr_dy_xyr, xyr_dz_xyz, xyz_dz_xyr, &
      &                     xyz_dx_pyz, pyz_dx_xyz, xyz_dy_xqz, &
      &                     xqz_dy_xyz, xyz_dz_xyr, xyr_dz_xyz
    use composition,  only: GasNum,       &! 
      &                     IdxG,         &!
      &                     MolWtWet       ! $B<>=a@.J,$NJ,;RNL(B
    use constants,    only: MolWtDry,     &! $B4%Ag@.J,$NJ,;RNL(B
      &                     Grav           ! $B=ENO2CB.EY(B



    !$B0EL[$N7?@k8@6X;_(B
    implicit none

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

    real(DP)                :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: pyz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_NDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
    integer                 :: f

    !------------------------------
    ! tendency of Km
    ! 

    ! Advection term
    !
    xyz_Adv  =                                             &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_KmNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_KmNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_KmNl ))    

    ! Numerical diffusion term 
    !
    xyz_NDiff =                                                                &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_KmBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_KmBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_KmBl ))))) 

    ! tendency
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_NDiff + xyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DKmDtAdv',    xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff',  xyz_NDiff(1:nx,1:ny,1:nz))


    !------------------------------
    ! tendency of potential temperature
    ! 

    ! Advection term
    !
    xyz_Adv = &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_PTempNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_PTempNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_PTempNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_DPTempBZDz )

    ! numerical diffusion term
    !
    xyz_NDiff = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_PTempBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_PTempBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_PTempBl ))))) 

    ! sum
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_NDiff + xyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff', xyz_NDiff(1:nx,1:ny,1:nz))

    
    !--------------------------------------
    ! Exner function
    !

    ! $B%U%i%C%/%99`$N7W;;(B. 4 $B<!@:EYCf?4:9J,$rMQ$$$F7W;;(B
    !
    xyz_Adv = &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_ExnerNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_ExnerNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_ExnerNl ))
!      & - xyz_xyr( xyr_VelZNl * xyr_DExnerBZDz ) 
    
    ! numerical diffusion term
    !
    xyz_NDiff = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_ExnerBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_ExnerBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_ExnerBl ))))) 
    xyz_NDiff = 0.0d0

    ! sum
    !
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_NDiff + xyz_Adv ) * FactorDExnerDtAdv

    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_NDiff(1:nx,1:ny,1:nz))


    !------------------------------
    ! tendency of mixing ratio
    ! 

    do f = 1, ncmax
      ! Advection term
      !
      xyzf_Adv(:,:,:,f) = &
        & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xyr( xyr_VelZNl * xyrf_DQMixBZDz(:,:,:,f) )

      ! numerical diffusion term
      !
      xyzf_NDiff(:,:,:,f) = &
        &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ))))) &
        &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ))))) &
        &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) )))))

    end do
    
    ! $BMn2<9`(B
    !
    call QMixFall

    ! sum
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_NDiff + xyzf_Adv + xyzf_Fall )

    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',  xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', xyzf_NDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------
    ! tendency of VelX
    ! 

    ! Advection term
    !
    pyz_Adv  = &
      & - pyz_VelXNl * pyz_xyz( xyz_c4dx_pyz( pyz_VelXNl ) )            &
      & - pyz_pqz( pqz_xqz( xqz_VelYNl ) * pqz_c4dy_pyz( pyz_VelXNl ) ) &
      & - pyz_pyr( pyr_xyr( xyr_VelZNl ) * pyr_c4dz_pyz( pyz_VelXNl ) )

    ! Numerical diffusion term 
    !
    pyz_NDiff = &
      & - NuHm * ( pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz( pyz_VelXBl ))))) &
      & - NuHm * ( pyz_dy_pqz(pqz_dy_pyz(pyz_dy_pqz(pqz_dy_pyz( pyz_VelXBl ))))) &
      & - NuVm * ( pyz_dz_pyr(pyr_dz_pyz(pyz_dz_pyr(pyr_dz_pyz( pyz_VelXBl )))))

    ! sum
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',   pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_NDiff(1:nx,1:ny,1:nz))

    !------------------------------
    ! tendency of VelY
    !     

    ! Advection term
    !
    xqz_Adv  = &
      & - xqz_pqz( pqz_pyz( pyz_VelXNl ) * pqz_c4dx_xqz( xqz_VelYNl ) ) &
      & - xqz_VelYNl * xqz_xyz( xyz_c4dy_xqz( xqz_VelYNl ) ) &
      & - xqz_xqr( xqr_xyr( xyr_VelZNl ) * xqr_c4dz_xqz( xqz_VelYNl ) )

    ! Numerical diffusion term
    !
    xqz_NDiff = &
      & - NuHm * ( xqz_dx_pqz(pqz_dx_xqz(xqz_dx_pqz(pqz_dx_xqz( xqz_VelYBl ))))) &
      & - NuHm * ( xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz( xqz_VelYBl ))))) &
      & - NuVm * ( xqz_dz_xqr(xqr_dz_xqz(xqz_dz_xqr(xqr_dz_xqz( xqz_VelYBl )))))

    ! sum
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_NDiff + xqz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelYDtAdv',   xqz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_NDiff(1:nx,1:ny,1:nz))

    !------------------------------
    ! tendency of VelZ
    ! 

    ! $B:n6HJQ?t(B
    !
    do f = 1, GasNum
      xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,IdxG(f)) / MolWtWet(IdxG(f))
    end do

    ! Advection term
    !
    xyr_Adv  = &
      & - xyr_pyr( pyr_pyz( pyz_VelXNl ) * pyr_c4dx_xyr( xyr_VelZNl ) ) &
      & - xyr_xqr( xqr_xqz( xqz_VelYNl ) * xqr_c4dy_xyr( xyr_VelZNl ) ) &
      & - xyr_VelZNl * xyr_xyz( xyz_c4dz_xyr( xyr_VelZNl ) )

    ! Numerical diffusion term
    !
    xyr_NDiff = &
      & - NuHm * ( xyr_dx_pyr(pyr_dx_xyr(xyr_dx_pyr(pyr_dx_xyr( xyr_VelZBl ))))) &
      & - NuHm * ( xyr_dy_xqr(xqr_dy_xyr(xyr_dy_xqr(xqr_dy_xyr( xyr_VelZBl ))))) & 
      & - NuVm * ( xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr( xyr_VelZBl )))))

    ! Buoyancy due to temperature disturbunce
    !
    xyr_BuoyT = Grav * xyr_xyz( xyz_PTempNl / xyz_PTempBZ)

    ! Buoyancy due to molecular weight
    !
    xyr_BuoyM =                                                 &
      & + Grav * xyr_xyz( sum(xyzf_QMixPerMolWt, 4) )           &
      &    / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt )          &
      & - Grav * xyr_xyz( sum(xyzf_QMixNl(:,:,:,1:GasNum), 4) ) &
      &    / ( 1.0d0 + xyr_QmixBZ ) 

    ! Buoyancy due to loading
    !
    xyr_BuoyD =                                                       &
      & - Grav * xyr_xyz( sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) ) &
      &    / ( 1.0d0 + xyr_QMixBZ )

    ! sum
    !
    xyr_DVelZDtNl = xyr_DVelZDtNl                        &
      &             + (                                  &
      &                 xyr_NDiff                        &
      &               + xyr_Adv                          &
      &               + (                                &
      &                 + xyr_BuoyM * FactorBuoyMolWt    &
      &                 + xyr_BuoyD * FactorBuoyLoading  &
      &                 + xyr_BuoyT * FactorBuoyTemp     &
      &                 )                                &
      &               )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_NDiff(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyT', xyr_BuoyT(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyM', xyr_BuoyM(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyD', xyr_BuoyD(1:nx,1:ny,1:nz))    

  contains

    subroutine QmixFall
      !
      ! $B1+N3$NMn2<$K$h$k0\N.$r5a$a$k(B. 
      ! 

      ! $B%b%8%e!<%k8F$S=P$7(B
      !
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyz_DensBZ
      
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
      xyzf_Fall = 0.0d0
      xyz_VelZFall = 0.0d0

      ! $BMn2<9`$N7W;;(B
      !
      do s = 1, RainNum

        iR = IdxR(s)

        !$B1+N3=*C<B.EY(B
        xyz_VelZFall = - 12.2d0 * FactorJ * ( xyzf_QMixAll(:,:,:,iR) ** 0.125d0 )
        
        ! $B%U%i%C%/%9$N7W;;(B
        !
        xyrf_QMixFallFlux(:,:,:,iR) =                                  &
          &  xyr_xyz (                                             &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall      &
          &  )

        ! $B>eC<$N%U%i%C%/%9$O%<%m(B
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! $B1+N3Mn2<$K$h$k;~4VJQ2=(B 
        !        
        xyzf_Fall(:,:,:,iR) =                                      &
          &  - xyz_c4dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall

  end subroutine Dynamics_Long_forcing
 

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
    ! $BNO3X%3%"(B: $BC;$$;~4V%9%F%C%W$GI>2A$9$k9`$N7W;;(B.
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x $BJ}8~$NG[Ns%5%$%:(B
      &                     jmin, jmax,       &! y $BJ}8~$NG[Ns%5%$%:(B
      &                     kmin, kmax,       &! z $BJ}8~$NG[Ns%5%$%:(B
      &                     nx, ny, nz         ! $BJ*M}NN0h$N%5%$%:(B
    use timeset,     only : TimeN, DelTimeShort
    use constants,   only : CpDry, CvDry, GasRDry ! $B4%Ag@.J,$NHfG.(B
    use basicset,    only : pyz_VPTempBZ,     &! $B4pK\>l$N290L(B
      &                     xqz_VPTempBZ,     &! $B4pK\>l$N290L(B
      &                     xyr_VPTempBZ       ! $B4pK\>l$N290L(B
!    use basicset,    only : xyz_VPTempBZ       !$B4pK\>l$N2>290L(B
    use average,     only : xyr_xyz, xqz_xyz, pyz_xyz
    use differentiate_center2,                &
      &              only : xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                     xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                     pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz
    use setmargin,  only : SetMargin_xyzf, SetMargin_xyz, &
      &                    SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

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

    real(DP) :: pyz_DVelXDtPGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtPGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtPGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtSWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtSWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtSWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !--------------------------------------
    ! initialize: Divergence of velocity
    !
    xyz_VelDivNs =                 &
      &   xyz_dx_pyz( pyz_VelXNs ) &
      & + xyz_dy_xqz( xqz_VelYNs ) &
      & + xyz_dz_xyr( xyr_VelZNs )
!    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelDiv', xyr_VelZNs(1:nx,1:ny,1:nz))
    
    !--------------------------------------
    ! VelX
    !
    pyz_DVelXDtSWF   =   pyz_dx_xyz( AlphaH * xyz_VelDivNs ) 
!    pyz_DVelXDtPGrad = - pyz_xyz( CpDry * xyz_VPTempBZ ) * pyz_dx_xyz( xyz_ExnerNs )      
    pyz_DVelXDtPGrad = - CpDry * pyz_VPTempBZ * pyz_dx_xyz( xyz_ExnerNs )
    pyz_DVelXDtNs =   pyz_DVelXDtPGrad + pyz_DVelXDtSWF

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * ( pyz_DVelXDtNs + pyz_DVelXDtNl )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_DVelXDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_DVelXDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)

    !--------------------------------------
    ! VelY
    !    
    xqz_DVelYDtSWF   =   xqz_dy_xyz( AlphaH * xyz_VelDivNs ) 
!    xqz_DVelYDtPGrad = - xqz_xyz( CpDry * xyz_VPTempBZ ) * xqz_dy_xyz( xyz_ExnerNs )
    xqz_DVelYDtPGrad = - CpDry * xqz_VPTempBZ * xqz_dy_xyz( xyz_ExnerNs )
    xqz_DVelYDtNs =   xqz_DVelYDtPGrad + xqz_DVelYDtSWF
    
    ! Time integration
    !
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * (xqz_DVelYDtNs + xqz_DVelYDtNl)

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_DVelYDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_DVelYDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)

    !----------------------------------------
    ! $B<!$N;~9o$N(B Exner $B4X?t(B
    !

    ! $BC;$$;~4V%9%F%C%W$GI>2A$9$k05NO$N<0$N(B tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency $B$N9g7W(B
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    xyz_ExnerAs = xyz_Exner( &
      & pyz_VelXAs,          &
      & xqz_VelYAs,          &
      & xyr_VelZNs,          &
      & xyz_VelDivNs,        &
      & xyz_ExnerNs,         &
      & xyr_DVelZDtNl,       &
      & xyz_DExnerDtNl,      &
      & xyz_DExnerDtNs       &
      & )

    ! Set Margin
    !
    call SetMargin_xyz( xyz_ExnerAs ) ! (inout)

    !--------------------------------------
    ! VelZ
    !
    xyr_DVelZDtSWF =  xyr_dz_xyz( AlphaV * xyz_VelDivNs ) 
    xyr_DVelZDtPGrad =                                          &
      & - CpDry * xyr_VPTempBZ                               &
      &   * (                                                &
      &         beta           * xyr_dz_xyz( xyz_ExnerAs )   &
      &       + (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerNs )   &
      &     )
    xyr_DVelZDtNs = xyr_DVelZDtPGrad + xyr_DVelZDtSWF

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNs + xyr_DVelZDtNl)

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_DVelZDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_DVelZDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  end subroutine Dynamics_Short_forcing

  
!!!--------------------------------------------------------------------!!!
  subroutine DynamicsVI_init()
    !
    !$B%(%/%9%J!<4X?t$r1"2rK!$G2r$/:]$KI,MW$H$J$k(B, $B78?t9TNs$NMWAG$r7h$a(B, 
    !LU $BJ,2r$r9T$&(B. 
    !

    ! $B%b%8%e!<%kFI$_9~$_(B
    use dc_types,   only : DP
!    use dc_message, only : MessageNotify
    use gridset,    only : imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                    ny,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                    nz                ! y $BJ}8~$NJ*M}NN0h$N>e8B(B
    use constants,  only : CpDry           ! $B4%Ag@.J,$NHfG.(B
    use timeset,    only : DelTimeShort
    use axesset,    only : dz        ! $B3J;R4V3V(B
    use basicset,   only : xyz_VelSoundBZ,  &!$B4pK\>l$N2;B.(B 
      &                    xyz_DensBZ,      &!$B4pK\>l$NL)EY(B
      &                    xyz_VPTempBZ      !$B4pK\>l$N290L(B
    use average,    only : xyr_xyz

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    real(DP)  :: DTS ! $BC;$$;~4V3J;R(B

    DTS = DelTimeShort

    !$BG[Ns$N3d$jIU$1(B
    allocate( A(1:nz) )
    allocate( B(1+1:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_F1BZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_F2BZ(imin:imax,jmin:jmax,kmin:kmax) )

    !----------------------------------------------------------------
    ! $B78?t9TNs$H6&DL$7$FMxMQ$5$l$kG[Ns$NCM$r7h$a$k(B
    !----------------------------------------------------------------

    !$B78?t9TNs$N7W;;(B
    !  A, B, C $B$r5a$a$k:](B, F1BZ $B$H(B F2BZ $B$O(B X $BJ}8~$K0lMM$J$N$G(B. 
    !  nx, ny $B$NCM$GBeI=$5$;$k$3$H$H$7$?(B. 
    xyz_F1BZ =                                                &
      &  ( xyz_VelSoundBZ ** 2.0d0 )                          &
      &   / (CpDry * xyz_DensBZ * (xyz_VPTempBZ ** 2.0d0))

    xyr_F2BZ =                                                &
      &  xyr_xyz(                                             &
      &    CpDry * xyz_DensBZ * ( xyz_VPTempBZ ** 2.0d0 )     &
      &   )
        
    A(1+1: nz-1) =                                &
      & (beta ** 2.0d0)                           &
      &    * xyz_F1BZ(nx,ny,1+1: nz-1)            &  
      &    * (DTS ** 2.0d0)                       &
      &    * (                                    &
      &          xyr_F2BZ(nx,ny,1+1: nz-1)        &
      &            / dz                           &
      &        + xyr_F2BZ(nx,ny,1  : nz-2)        &
      &            / dz                           &
      &       )                                   &
      &    / dz                                   &
      & + 1.0d0

    A(1) =                                        &
      & (beta ** 2.0d0)                           &
      &   * xyz_F1BZ(nx,ny,1)                     &
      &   * xyr_F2BZ(nx,ny,1)                     &
      &     / dz                                  &
      &   * (DTS ** 2.0d0)                        &
      &   / dz                                    &
      & + 1.0d0                                         

    A(nz) =                                       &
      & (beta ** 2.0d0)                           &
      &   * xyz_F1BZ(nx,ny,nz)                    &
      &   * xyr_F2BZ(nx,ny,nz-1)                  &
      &     / dz                                  &
      &   * (DTS ** 2.0d0)                        &
      &   / dz                                    &
      & + 1.0d0                                         
    
    B(1+1:nz) =                                   &
      & - (beta ** 2.0d0)                         &
      &   * xyz_F1BZ(nx,ny,1:nz-1)                & 
      &   * xyr_F2BZ(nx,ny,1:nz-1)                &
      &   * (DTS ** 2.0d0)                        &
      &   / ( dz * dz ) 
    
    C(1: nz-1) =                                  &
      & - ( beta ** 2.0d0 )                       &
      &   * xyz_F1BZ(nx,ny,1+1:nz)                &
      &   * xyr_F2BZ(nx,ny,1  :nz-1)              & 
      &   * (DTS ** 2.0d0) &
      &   / ( dz * dz )


    !----------------------------------------------------------------
    ! $B78?t9TNs$r(B LU $BJ,2r(B
    !----------------------------------------------------------------
    !$BG[Ns$NBg$-$5$rJ]4I(B
    N   = nz  !$B78?t9TNs(B/$B2~9TNs$N<!?t(B, $B@09g@#K!(B
    M   = nx * ny 
                               !$BJ}Dx<0$NAH?t(B
    NUD = 1                    !$B78?t9TNs$N>e;03QItJ,$NBSI}(B
    NLD = 1                    !$B78?t9TNs$N2<;03QItJ,$NBSI}(B
    NAL = NLD                  !LU $BJ,2r$N7k2L(B L $B$N@09g@#K!(B
    NA  = NUD + NLD + 1

    !$BG[Ns$N3d$jEv$F(B
!    allocate( AL1(N), AL2(NAL, N), AU2(NA, N), IP(N) )
    allocate( AL1(N), IP(N) )

    !LU $BJ,2r$N<B9T(B
    !  LAPACK $B$NMxMQ(B
    call ResolvLU_Lapack( )

   
  end subroutine DynamicsVI_init
  

!!!--------------------------------------------------------------------!!!
  function xyz_Exner(      &
    & pyz_VelXAs,          &
    & xqz_VelYAs,          &
    & xyr_VelZNs,          &
    & xyz_VelDivNs,        &
    & xyz_ExnerNs,         &
    & xyr_DVelZDtNl,       &
    & xyz_DExnerDtNl,      &
    & xyz_DExnerDtNs       &
    & )
    !
    !$B1"2rK!$rMQ$$$?%(%/%9%J!<4X?t$N7W;;(B. 
    !
    
    ! $B%b%8%e!<%k$NFI$_9~$_(B
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x $BJ}8~$NG[Ns$N2<8B(B
      &                  imax,            &! x $BJ}8~$NG[Ns$N>e8B(B
      &                  jmin,            &! y $BJ}8~$NG[Ns$N2<8B(B
      &                  jmax,            &! y $BJ}8~$NG[Ns$N>e8B(B
      &                  kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                  kmax,            &! z $BJ}8~$NG[Ns$N>e8B(B
      &                  nx, ny, nz        ! $BJ*M}NN0h$NBg$-$5(B
    use constants,only : CpDry             ! $B4%Ag@.J,$NHfG.(B
    use timeset,  only : DelTimeShort, TimeN
    use basicset, only : xyz_VelSoundBZ,  &!$B4pK\>l$N2;B.(B
      &                  xyz_DensBZ,      &!$B4pK\>l$NL)EY(B
      &                  xyz_VPTempBZ,    &!$B4pK\>l$N2>290L(B
      &                  xyr_VPTempBZ      !$B4pK\>l$N2>290L(B
    use axesset,  only : dz
    use average,  only : xyr_xyz
    use differentiate_center2, &
      &           only : xyr_dz_xyz, xyz_dz_xyr, xyz_dx_pyz, xyz_dy_xqz
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 

    !$B0EL[$N7?@k8@6X;_(B
    implicit none
    
    !$BF~=PNOJQ?t(B
    real(DP), intent(in)   :: pyz_VelXAs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !$BB.EY(B u [$B&S(B+$B&$&S(B]
    real(DP), intent(in)   :: xqz_VelYAs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !$BB.EY(B v [$B&S(B+$B&$&S(B]
    real(DP), intent(in)   :: xyr_VelZNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !$BB.EY(B w [$B&S(B]
    real(DP), intent(in)   :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)   :: xyr_DVelZDtNl &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z $BJ}8~$N30NO9`(B[t]
    real(DP), intent(in)   :: xyz_DExnerDtNl &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z $BJ}8~$N30NO9`(B[t]
    real(DP), intent(in)   :: xyz_DExnerDtNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z $BJ}8~$N30NO9`(B[t]
    real(DP), intent(in)   :: xyz_ExnerNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !$BL5<!8505NO(B
    real(DP)               :: xyz_Exner &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !$BL5<!8505NO(B[$B&S(B+$B&$&S(B]

    !$BJQ?tDj5A(B
    real(DP)  :: D1(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: D2(1:nx,1:ny,1:nz)  
    real(DP)  :: D(nx*ny,nz)
    real(DP)  :: E(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: F(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: DTS ! $BC;$$;~4V3J;R4V3V(B
    integer   :: ix, jy, kz

    ! Initialize
    !
    DTS = DelTimeShort
    xyz_Exner = 0.0d0

    !$B9TNs7W;;$N$?$a$N78?t(B
    E =   &
      & - ( 1.0d0 - beta ) * xyr_dz_xyz( xyz_ExnerNs )              &
      & + ( AlphaV * xyr_dz_xyz( xyz_VelDivNs ) + xyr_DVelZDtNl )   &
      &    / ( CpDry * xyr_VPTempBZ ) 

     F = - beta * xyz_F1BZ * DTS                                     &
      & * xyz_dz_xyr(                                               &
      &    xyr_xyz( xyz_DensBZ * xyz_VPTempBZ)                      &
      &    * (                                                      &
      &         xyr_VelZNs                                          &
      &       - xyr_xyz(CpDry * xyz_VPTempBZ) * DTS                 &
      &         * (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerNs )        &
      &       + xyr_dz_xyz( AlphaV * xyz_VelDivNs ) * DTS           &
      &       + xyr_DVelZDtNl * DTS                                 &
      &      )                                                      &
      &   )                                                         &
      & + (xyz_DExnerDtNs + xyz_DExnerDtNl ) * DTS

    D1 = xyz_ExnerNs                                                &
      & - (1.0d0 - beta)                                            &
      &   * xyz_F1BZ * DTS                                          &
      &   * xyz_dz_xyr(                                             &
      &       xyr_xyz(xyz_DensBZ * xyz_VPTempBZ) * xyr_VelZNs       &
      &     )                                                       &
      & - (xyz_VelSoundBZ ** 2.0d0) * DTS                           &
      &   / ( CpDry * xyz_VPTempBZ )                                &
      &   * ( xyz_dx_pyz( pyz_VelXAs ) +  xyz_dy_xqz( xqz_VelYAs ) )&
      & + F

    D1(:,:,1) = D1(:,:,1)                           &
      & - beta * xyz_F1BZ(:,:,1) * (DTS ** 2.0d0)   &
      &   * xyr_F2BZ(:,:,1-1) * E(:,:,1-1)          &
      &   / dz
    
    D1(:,:,nz) = D1(:,:,nz)                         &
      & + beta * xyz_F1BZ(:,:,nz) * (DTS ** 2.0d0)  &
      &   * xyr_F2BZ(:,:,nz) * E(:,:,nz)            &
      &   / dz

!    call HistoryAutoPut(TimeN, 'D', D1(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'E', E(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'F', F(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENs', xyz_DExnerDtNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENl', xyz_DExnerDtNl(1:nx,1:ny,1:nz))

    D2 = D1(1:nx,1:ny,1:nz)

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          D(ix + nx * (jy - 1), kz) =  D2(ix,jy,kz)
        end do
      end do
    end do

    !-----------------------------------------------------------
    !$BO"N)0l<!J}Dx<0$N2r$r5a$a$k(B
    !------------------------------------------------------------

    !$B2r$N7W;;(B
    !  LAPACK $BMxMQ(B
    call LinSolv_Lapack( D )

    !$BLa$jCM$r=PNO(B
    do kz = 1, nz
      do jy = 1, ny 
        do ix = 1, nx 
          xyz_Exner(ix,jy,kz) = D(ix + nx * (jy - 1 ), kz)
        end do
      end do
    end do

  end function xyz_Exner

!!!--------------------------------------------------------------------!!!
  subroutine ResolvLU_Lapack(  )
    !
    !$B<B(B 3 $B9`9TNs$N(B LU $BJ,2r(B($BG\@:EY(B). LAPACK $BMxMQ(B
    !

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BJQ?tDj5A(B
    integer    :: INFO  !$B2r$N%3%s%G%#%7%g%s%A%'%C%/(B
    
    !$BJQ?t$N=i4|2=(B
    INFO = 0
    
    !$B2r9TNs$N7W;;(B. LAPACK $B$r;HMQ(B. 
    call DGTTRF(N, C, A, B, AL1, IP, INFO)
    
    !$B2r$N%3%s%G%#%7%g%s$r%A%'%C%/(B. 
!    if (INFO /= 0) then
!      call MessgaeNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if
    
  end subroutine ResolvLU_Lapack
  

!!!--------------------------------------------------------------------!!!
  subroutine LinSolv_Lapack( X )
    !
    !LU $BJ,2r$5$l$?<B(B 3 $B9`9TNs$NO"N)(B 1 $B<!J}Dx<0(B($BG\@:EY(B). LAPACK $BMxMQ(B
    !

    !$B0EL[$N7?@k8@6X;_(B
    implicit none
    
    !$BJQ?tDj5A(B
    real(DP), intent(inout) :: X(M, N)     !$BDj?t(B/$B2r9TNs(B
    real(DP)                :: TX(N, M)    !$B2r9TNs$rE>CV$7$?$b$N(B
    integer                :: NRHS         !
    integer                :: INFO
    integer                :: LDB
    character(1),parameter :: TRANS = 'N'

    !$BJQ?t$N=i4|2=(B
    NRHS = M
    INFO = 0
    LDB  = N
    TX = transpose( X )
    
    !$B2r9TNs$N7W;;(B. LAPACK $B$r;HMQ(B. 
    call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)

    !$B2r$N=PNO(B
    X = transpose( TX )
     
    !$B2r$N%3%s%G%#%7%g%s$r%A%'%C%/(B. 
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if
     
  end subroutine LinSolv_Lapack



  subroutine  Dynamics_Tendency_output
    !
    ! $B%U%!%$%k=PNO$NDj5A(B
    !
    
    !$B%b%8%e!<%k8F$S=P$7(B
    use gtool_historyauto, only: HistoryAutoAddVariable
    use composition,       only: SpcWetSymbol
    use gridset,           only: ncmax          ! $BJ*<A?t(B
    
    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BJQ?tDj5A(B
    integer :: l

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of potential temperature',  &
      & units='K.s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of potential temperature',&
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of exner function',  &
      & units='s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DExnerDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of exner function',&
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='CDensAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of cloud density',  &
      & units='kg.m-3.s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='CDensDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of cloud density',&
      & units='kg.m-3.s-1',    &
      & xtype='double')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtAdv', &
        & dims=(/'x','y','z','t'/),     &
        & longname='Advection term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
      
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtDiff', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Diffusion term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtFall', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Fall term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

    end do

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (x)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelXDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of velocity (x)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (x)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (x)',  &
      & units='m.s-2',    &
      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='VelXTndNs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='Velocity Tendency (x)',  &
!      & units='m.s-2',    &
!      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (y)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelYDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of velocity (y)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (y)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (y)',  &
      & units='m.s-2',    &
      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='VelYTndNs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='Velocity Tendency (y)',  &
!      & units='m.s-2',    &
!      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (z)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelZDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of Velocity (z)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyT',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (Temperature)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyM',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (MolWt)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (Drag)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (z)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (z)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='VelDiv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='velocity divergence',  &
      & units='s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection of Km',  &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtDiff', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Diffusion term of Km',  &
      & units='s-1',    &
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

  end subroutine Dynamics_Tendency_output

end module DynamicsHEVI

