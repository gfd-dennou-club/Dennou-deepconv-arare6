!= Module DynamicsHEVI
!
! Authors::   $B?y;39L0lO/(B(SUGIYAMA Ko-ichiro), ODAKA Masatsugu 
! Version::   $Id: dynamics_hevi_v2.f90,v 1.8 2014/07/08 00:58:06 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module DynamicsHEVI_v2
  !
  ! $BNO3X%3%"(B. 
  !   $B%?%$%`%9%W%j%C%HK!$rMxMQ(B. $B2;GH%b!<%I$H$=$l0J30$rJL!9$N;~4V9o$_$G2r$/(B. 
  !   $BC;$$;~4V%9%F%C%W$N7W;;$K$O(B HE-VI $BK!$rMxMQ(B.
  !
  !   $BHyJ,J?6Q1i;;%b%8%e!<%k$r;H$o$:$K=q$-2<$7$?HG(B. 
  !   $B7W;;B.EY$O(B v1 $B$h$j$b(B 2 $BG\0J>eAa$$$,(B, $B%G%P%C%0$O$7$K$/$$(B. 
  !
  ! Note: 
  !  * $B%(%/%9%J!<4X?t$N6u4VJ}8~$NN%;62=$O(B 2 $B<!@:EY$G$"$k$?$a(B, $B5$0579EY(B
  !    $BNO9`$N7W;;%W%m%0%i%`$K$*$$$F(B differentiate_center4 $B%b%8%e!<%k$r(B
  !    $B;XDj$9$k$3$H$O$G$-$J$$$N$GCm0U(B.
  !

  !$B%b%8%e!<%kFI$_9~$_(B
  use dc_types,   only : DP

  !$B0EL[$N7?@k8@6X;_(B
  implicit none

  !$BB0@-$N;XDj(B
  private
  
  real(DP), save,private :: beta  = 5.0d-1     !$B%/%i%s%/%K%3%k%=%sK!$J$i(B 0.5
                                               !$B40A41"2rK!$J$i(B 1 
  integer, save,private  :: N = 10             !$B78?t9TNs(B/$B2~9TNs$N<!?t(B, $B@09g@#K!(B
  integer, save,private  :: M = 10             !$BJ}Dx<0$NAH?t(B
  integer, save,private  :: NUD = 1            !$B78?t9TNs$N>e;03QItJ,$NBSI}(B
  integer, save,private  :: NLD = 1            !$B78?t9TNs$N2<;03QItJ,$NBSI}(B
  integer, save,private  :: NAL = 1            !LU $BJ,2r$N7k2L(B L $B$N@09g@#K!(B
  integer, save,private  :: NA = 3             !NUD + NLD + 1

  real(DP), allocatable, save,private :: xyz_F1BZ(:,:,:)
                                               !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B
  real(DP), allocatable, save,private :: xyr_CpVPTempBZ(:,:,:)
                                               !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B
  real(DP), allocatable, save,private :: xyr_CpDensVPTemp2BZ(:,:,:)
                                               !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B
  real(DP), allocatable, save,private :: xyr_DensVPTempBZ(:,:,:) 
                                               !$B78?t9TNs$N7W;;$KMQ$$$kG[Ns(B

  real(DP), allocatable, save,private :: A(:)  !$B78?t9TNs$NBP3Q@.J,(B
  real(DP), allocatable, save,private :: B(:)  !$B78?t9TNs$N>e;03QItJ,(B
  real(DP), allocatable, save,private :: C(:)  !$B78?t9TNs$N2<;03QItJ,(B
  real(DP), allocatable, save,private :: AL1(:)!LU $BJ,2r$N7k2L(B L (1 $B<!85G[Ns(B)
  integer,  allocatable, save,private :: IP(:) !$BItJ,%T%\%C%H8r49$N>pJs$r3JG<(B

  real(DP), save,private :: AlphaH = 0.0d0     !$B2;GH8:?j9`$N8:?j78?t(B ($B?eJ?J}8~(B)
  real(DP), save,private :: AlphaV = 0.0d0     !$B2;GH8:?j9`$N8:?j78?t(B ($B1tD>J}8~(B)
  real(DP), save,private :: NuHh   = 0.0d0     !$BG.$KBP$9$k?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save,private :: NuVh   = 0.0d0     !$BG.$KBP$9$k?tCMG4@-$N78?t(B ($B1tD>J}8~(B)
  real(DP), save,private :: NuHm   = 0.0d0     !$B1?F0NL$KBP$9$k?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save,private :: NuVm   = 0.0d0     !$B1?F0NL$KBP$9$k?tCMG4@-$N78?t(B ($B1tD>J}8~(B)

  character(*), parameter:: module_name = 'DynamicHEVI'
                                               ! $B%b%8%e!<%k$NL>>N(B.
                                               ! Module name
  real(DP), save,private :: FactorBuoyTemp    = 1.0d0
                                               !$BIbNO(B ($B29EY$N4sM?(B) $B$NM-L5(B
                                               !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save,private :: FactorBuoyMolWt   = 1.0d0
                                               !$BIbNO(B ($BJ,;RNL8z2L(B) $B$NM-L5(B
                                               !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save,private :: FactorBuoyLoading = 1.0d0
                                               !$BIbNO(B ($B2Y=E8z2L(B) $B$NM-L5(B
                                               !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorDExnerDtAdv    = 1.0d0     !$B%(%/%9%J!<4X?t$N0\N.$NM-L5(B
                                               !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.
  real(DP), save, private :: FactorDExnerDtExpnd  = 1.0d0     !$B%(%/%9%J!<4X?t$NKDD%9`$NM-L5(B
                                               !$B9MN8$7$J$$>l9g$OCM$r%<%m$K$9$k(B.

  real(DP), allocatable, save,private :: xyr_Dummy(:,:,:) !$B%@%_!<(B

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:)   !$B4pK\>l$N1tD>HyJ,(B
  real(DP), allocatable, save,private :: xyr_DExnerBZDz(:,:,:)   !$B4pK\>l$N1tD>HyJ,(B
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:) !$B4pK\>l$N1tD>HyJ,(B
  
  ! public 
  !
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Init
    !
    ! $BNO3X%3%"(B $B=i4|2=%k!<%A%s(B
    !

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
    real(DP)  :: AlphaNDiff = 1.0d-3  !4$B<!$N?tCM3H;6$N78?t(B. CReSS $B%^%K%e%"%k$h$j(B
    real(DP)  :: NDiffRatio = 1.0d0   !$BB.EY$KBP$9$kG4@-$r>e$2$k>l9g$O?t;z$r(B 1 $B0J>e$K$9$k(B. 
    integer   :: unit                 !$BAuCVHV9f(B
    integer   :: f

    !-------------------------------------------------------------------
    ! Namelist $B$+$i>pJs$r<hF@$9$k(B
    !
    NAMELIST /Dynamics_nml/                                    &
         & AlphaSound, AlphaNDiff, NDiffRatio, beta,           &
         & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
         & FactorDExnerDtAdv, FactorDExnerDtExpnd

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

    ! $B1tD>1"2rK!$rMQ$$$k$?$a$K(B, $B9TNs$N=`Hw$r9T$&(B. 
    !
    call Dynamics_VI_init

    ! tendency $B$N=PNO(B
    !
    call Dynamics_Tendency_output

    ! $BG[Ns$NMQ0U(B
    !
    allocate( xyr_Dummy(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DExnerBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_Dummy      = 0.0d0
    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
!    xyr_DExnerBZDz = xyr_dz_xyz( xyz_ExnerBZ )
    xyr_DExnerBZDz = 0.0d0
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
                                          !$BG[Ns%5%$%:(B
      &                   FlagCalc3D      !$B%U%i%0(B

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


    !------------------------------------------------------
    ! $B>r7oJ,4t(B
    !
    if ( FlagCalc3D ) then 
    
      call Dynamics3D_Long_forcing(       &
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
     
    else

      call Dynamics2D_Long_forcing(       &
        & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
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

    end if

  end subroutine Dynamics_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics3D_Long_forcing(   &
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
    ! $B0\N.7W;;(B (3D $BHG(B)
    !
    !   $B0\N.(B: 4 $B<!Cf1{:9J,(B
    !   $B?tCM3H;6(B (4 $B3,(B): 2 $B<!Cf1{:9J,(B
    !
    ! $B%j!<%W%U%m%C%0$G(B, $B0\N.$rCf1{:9J,$G7W;;$9$k$?$a$K(B, 
    ! $B?tCM3H;69`$rDI2C$7$F$$$k(B. 
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,    only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x $BJ}8~$NG[Ns%5%$%:(B
      &                     jmin, jmax,       &! y $BJ}8~$NG[Ns%5%$%:(B
      &                     kmin, kmax,       &! z $BJ}8~$NG[Ns%5%$%:(B
      &                     nx, ny, nz,       &! $BJ*M}NN0h$N%5%$%:(B
      &                     ncmax              ! $BJ*<A?t(B
    use timeset,     only : TimeN
    use axesset,     only : dx, dy, dz         ! $B3J;R4V3V(B
    use composition, only : SpcWetSymbol

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

    real(DP)             :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: pyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xqz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xqz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_nDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    integer              :: f


    !--------------------------------------------------------------------
    ! $BMpN.3H;678?t(B

    ! $B0\N.$*$h$S?tCM3H;6(B
    !    
    call AdvC4_nDiff_xyz( xyz_KmBl, xyz_KmNl, xyr_Dummy ) !(IN)
    
    ! tendency $B$N99?7(B
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_nDiff + xyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! $B290L(B

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyz( xyz_PTempBl, xyz_PTempNl, xyr_DPTempBZDz ) !(IN)
    
    ! tendency $B$N99?7(B
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_nDiff + xyz_Adv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',   xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff',  xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! $B%(%/%9%J!<4X?t(B

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyz( xyz_ExnerBl, xyz_ExnerNl, xyr_DExnerBZDz ) !(IN)
    
    ! tendency $B$N99?7(B
    !
    xyz_nDiff = 0.0d0
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_nDiff + xyz_Adv ) * FactorDExnerDtAdv
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))
    
    !--------------------------------------------------------------------
    ! $B:.9gHf(B
    ! 

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyzf( xyzf_QMixBl, xyzf_QMixNl, xyrf_DQMixBZDz ) !(IN)

    ! $BMn2<9`(B
    !
    call QMixFall

    ! tendency $B$N99?7(B    
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_nDiff + xyzf_Adv + xyzf_Fall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_nDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', &
        & xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------------------------------------------
    ! VelX, VelY, VelZ
    ! 

    ! $B0\N.9`!&?tCM3H;69`$r$^$H$a$F7W;;(B
    !
    call AdvC4_nDiff_pyz_xqz_xyr

    ! tendency of VelX
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )

    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_nDiff(1:nx,1:ny,1:nz))

    ! tendency of VelY
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_NDiff + xqz_Adv )

    call HistoryAutoPut(TimeN, 'DVelYDtAdv',  xqz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_nDiff(1:nx,1:ny,1:nz))

    ! Buoyancy 
    ! 
    call BuoyancyLong_xyr

    ! tendency of VelZ
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

    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_nDiff(1:nx,1:ny,1:nz))    
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
      use average,     only : xyr_xyz
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyzf_QMixBZ, &!$B4pK\>l$N:.9gHf(B
        &                     xyz_DensBZ
      use differentiate_center4, &
        &              only : xyz_dz_xyr
      
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
        xyzf_Fall(:,:,:,iR) =                                          &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall


    subroutine AdvC4_nDiff_pyz_xqz_xyr

      implicit none

      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            pyz_Adv(i,j,k) =                                                   &
              & - pyz_VelXNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   pyz_VelXNl(i+1,j,k) - pyz_VelXNl(i-1,j,k) )   &
              &     - fct2 * (   pyz_VelXNl(i+2,j,k) + pyz_VelXNl(i+1,j,k)     &
              &                - pyz_VelXNl(i-1,j,k) - pyz_VelXNl(i-2,j,k) )   &
              &     ) * 5.0d-1 / dx                                            &
              & - (                                                            &
              &   + ( xqz_VelYNl(i+1,j,k) + xqz_VelYNl(i,j,k) )                &  
              &     * (                                                        &  
              &         fct1 * ( pyz_VelXNl(i,j+1,k) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j+2,k) - pyz_VelXNl(i,j-1,k) )   &
              &       )                                                        &
              &   + ( xqz_VelYNl(i+1,j-1,k) + xqz_VelYNl(i,j-1,k) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j-1,k) )   &
              &       - fct2 * ( pyz_VelXNl(i,j+1,k) - pyz_VelXNl(i,j-2,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dy                                              &
              & - (                                                            &
              &   + ( xyr_VelZNl(i+1,j,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+2) - pyz_VelXNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i+1,j,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j,k-1) )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
            
          end do
        end do
      end do
      
      pyz_Adv(imin:imin+1,:,:) = 0.0d0
      pyz_Adv(imax-1:imax,:,:) = 0.0d0
      pyz_Adv(:,jmin:jmin+1,:) = 0.0d0
      pyz_Adv(:,jmax-1:jmax,:) = 0.0d0
      pyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      pyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xqz_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j+1,k) + pyz_VelXNl(i,j,k) )                &
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i+1,j,k) - xqz_VelYNl(i,j,k)   )   &
              &       - fct2 * ( xqz_VelYNl(i+2,j,k) - xqz_VelYNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j+1,k) + pyz_VelXNl(i-1,j,k) )            &
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k)   - xqz_VelYNl(i-1,j,k) )   &
              &       - fct2 * ( xqz_VelYNl(i+1,j,k) - xqz_VelYNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - xqz_VelYNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xqz_VelYNl(i,j+1,k) - xqz_VelYNl(i,j-1,k) )   &
              &     - fct2 * (   xqz_VelYNl(i,j+2,k) + xqz_VelYNl(i,j+1,k)     &
              &                - xqz_VelYNl(i,j-1,k) - xqz_VelYNl(i,j-2,k) )   &
              &     ) * 5.0d-1 / dy                                            &
              & - (                                                            &
              &   + ( xyr_VelZNl(i,j+1,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k+1) - xqz_VelYNl(i,j,k)   )   &
              &       - fct2 * ( xqz_VelYNl(i,j,k+2) - xqz_VelYNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i,j+1,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k)   - xqz_VelYNl(i,j,k-1) )   &
              &       - fct2 * ( xqz_VelYNl(i,j,k+1) - xqz_VelYNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
          end do
        end do
      end do
      
      xqz_Adv(imin:imin+1,:,:) = 0.0d0
      xqz_Adv(imax-1:imax,:,:) = 0.0d0
      xqz_Adv(:,jmin:jmin+1,:) = 0.0d0
      xqz_Adv(:,jmax-1:jmax,:) = 0.0d0
      xqz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xqz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyr_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j,k+1) + pyz_VelXNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i+2,j,k) - xyr_VelZNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j,k+1) + pyz_VelXNl(i-1,j,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i-1,j,k) )   &
              &       - fct2 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - (                                                            &
              &   + ( xqz_VelYNl(i,j,k+1) + xqz_VelYNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j+1,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i,j+2,k) - xyr_VelZNl(i,j-1,k) )   &
              &       )                                                        &
              &   + ( xqz_VelYNl(i,j-1,k+1) + xqz_VelYNl(i,j-1,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i,j-1,k) )   &
              &       - fct2 * ( xyr_VelZNl(i,j+1,k) - xyr_VelZNl(i,j-2,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dy                                              &
              & - xyr_VelZNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xyr_VelZNl(i,j,k+1) - xyr_VelZNl(i,j,k-1) )   &
              &     - fct2 * (   xyr_VelZNl(i,j,k+2) + xyr_VelZNl(i,j,k+1)     &
              &                - xyr_VelZNl(i,j,k-1) - xyr_VelZNl(i,j,k-2) )   &
              &     ) * 5.0d-1 / dz
          end do
        end do
      end do
      
      xyr_Adv(imin:imin+1,:,:) = 0.0d0
      xyr_Adv(imax-1:imax,:,:) = 0.0d0
      xyr_Adv(:,jmin:jmin+1,:) = 0.0d0
      xyr_Adv(:,jmax-1:jmax,:) = 0.0d0
      xyr_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyr_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            pyz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + pyz_VelXBl(i+2,j,k)                  &
              &     + pyz_VelXBl(i-2,j,k)                  &
              &     - pyz_VelXBl(i+1,j,k) * 4.0d0          &
              &     - pyz_VelXBl(i-1,j,k) * 4.0d0          &
              &     + pyz_VelXBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + pyz_VelXBl(i,j+2,k)                  &
              &     + pyz_VelXBl(i,j-2,k)                  &
              &     - pyz_VelXBl(i,j+1,k) * 4.0d0          &
              &     - pyz_VelXBl(i,j-1,k) * 4.0d0          &
              &     + pyz_VelXBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        & 
              &     + pyz_VelXBl(i,j,k+2)                  &
              &     + pyz_VelXBl(i,j,k-2)                  &
              &     - pyz_VelXBl(i,j,k+1) * 4.0d0          &
              &     - pyz_VelXBl(i,j,k-1) * 4.0d0          &
              &     + pyz_VelXBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
            
          end do
        end do
      end do
      
      pyz_nDiff(imin:imin+1,:,:) = 0.0d0
      pyz_nDiff(imax-1:imax,:,:) = 0.0d0
      pyz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      pyz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      pyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      pyz_nDiff(:,:,kmax-1:kmax) = 0.0d0

      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xqz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xqz_VelYBl(i+2,j,k)                  &
              &     + xqz_VelYBl(i-2,j,k)                  &
              &     - xqz_VelYBl(i+1,j,k) * 4.0d0          &
              &     - xqz_VelYBl(i-1,j,k) * 4.0d0          &
              &     + xqz_VelYBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xqz_VelYBl(i,j+2,k)                  &
              &     + xqz_VelYBl(i,j-2,k)                  &
              &     - xqz_VelYBl(i,j+1,k) * 4.0d0          &
              &     - xqz_VelYBl(i,j-1,k) * 4.0d0          &
              &     + xqz_VelYBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        & 
              &     + xqz_VelYBl(i,j,k+2)                  &
              &     + xqz_VelYBl(i,j,k-2)                  &
              &     - xqz_VelYBl(i,j,k+1) * 4.0d0          &
              &     - xqz_VelYBl(i,j,k-1) * 4.0d0          &
              &     + xqz_VelYBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xqz_nDiff(imin:imin+1,:,:) = 0.0d0
      xqz_nDiff(imax-1:imax,:,:) = 0.0d0
      xqz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xqz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xqz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xqz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyr_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xyr_VelZBl(i+2,j,k)                  &
              &     + xyr_VelZBl(i-2,j,k)                  &
              &     - xyr_VelZBl(i+1,j,k) * 4.0d0          &
              &     - xyr_VelZBl(i-1,j,k) * 4.0d0          &
              &     + xyr_VelZBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j+2,k)                  &
              &     + xyr_VelZBl(i,j-2,k)                  &
              &     - xyr_VelZBl(i,j+1,k) * 4.0d0          &
              &     - xyr_VelZBl(i,j-1,k) * 4.0d0          &
              &     + xyr_VelZBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j,k+2)                  &
              &     + xyr_VelZBl(i,j,k-2)                  &
              &     - xyr_VelZBl(i,j,k+1) * 4.0d0          &
              &     - xyr_VelZBl(i,j,k-1) * 4.0d0          &
              &     + xyr_VelZBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyr_nDiff(imin:imin+1,:,:) = 0.0d0
      xyr_nDiff(imax-1:imax,:,:) = 0.0d0
      xyr_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xyr_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xyr_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyr_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_pyz_xqz_xyr
    
    
    subroutine AdvC4_nDiff_xyz( xyz_VarBl, xyz_VarNl, xyr_DVarBZDz )
      
      implicit none
      
      real(DP), intent(in)  :: xyz_VarBl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyz_VarNl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyr_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax)
      
      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0

      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyz_Adv(i,j,k) =                                                  &
              & - (                                                           &
              &      pyz_VelXNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i+2,j,k) - xyz_VarNl(i-1,j,k) ) &
              &          )                                                    &
              &    + pyz_VelXNl(i-1,j,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i-1,j,k) ) &
              &          - fct2 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i-2,j,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dx                                             &
              & - (                                                           &
              &      xqz_VelYNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j+1,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j+2,k) - xyz_VarNl(i,j-1,k) ) &
              &          )                                                    &
              &    + xqz_VelYNl(i,j-1,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j-1,k) ) &
              &          - fct2 * ( xyz_VarNl(i,j+1,k) - xyz_VarNl(i,j-2,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dy                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+2) - xyz_VarNl(i,j,k-1) ) &
              &          )                                                    &
              &    + xyr_VelZNl(i,j,k-1)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j,k-1) ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k-2) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dz                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)   * xyr_DVarBZDz(i,j,k)                & 
              &    + xyr_VelZNl(i,j,k-1) * xyr_DVarBZDz(i,j,k-1)              &
              &   ) * 5.0d-1 
          end do
        end do
      end do
      
      xyz_Adv(imin:imin+1,:,:) = 0.0d0
      xyz_Adv(imax-1:imax,:,:) = 0.0d0
      xyz_Adv(:,jmin:jmin+1,:) = 0.0d0
      xyz_Adv(:,jmax-1:jmax,:) = 0.0d0
      xyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      ! 4 $B<!$N?tCM3H;6(B: 2 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyz_nDiff(i,j,k) =                            &
              & - (                                       &
              &     + xyz_VarBl(i+2,j,k)                  &
              &     + xyz_VarBl(i-2,j,k)                  &
              &     - xyz_VarBl(i+1,j,k) * 4.0d0          &
              &     - xyz_VarBl(i-1,j,k) * 4.0d0          &
              &     + xyz_VarBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j+2,k)                  &
              &     + xyz_VarBl(i,j-2,k)                  &
              &     - xyz_VarBl(i,j+1,k) * 4.0d0          &
              &     - xyz_VarBl(i,j-1,k) * 4.0d0          &
              &     + xyz_VarBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHh / ( dy ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j,k+2)                  &
              &     + xyz_VarBl(i,j,k-2)                  &
              &     - xyz_VarBl(i,j,k+1) * 4.0d0          &
              &     - xyz_VarBl(i,j,k-1) * 4.0d0          &
              &     + xyz_VarBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyz_nDiff(imin:imin+1,:,:) = 0.0d0
      xyz_nDiff(imax-1:imax,:,:) = 0.0d0
      xyz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xyz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyz
    
    
    subroutine AdvC4_nDiff_xyzf( xyzf_VarBl, xyzf_VarNl, xyrf_DVarBZDz ) 

      implicit none
      
      real(DP), intent(in)  :: xyzf_VarBl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyzf_VarNl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyrf_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      
      real(DP)              :: fct1, fct2
      integer               :: s, i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = jmin + 2, jmax - 2
            do i = imin + 2, imax - 2
              
              xyzf_Adv(i,j,k,s) =                                                     &
                & - (                                                                 &
                &      pyz_VelXNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i+2,j,k,s) - xyzf_VarNl(i-1,j,k,s) ) &
                &          )                                                          &
                &    + pyz_VelXNl(i-1,j,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i-1,j,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i-2,j,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dx                                                   &
                & - (                                                                 &
                &      xqz_VelYNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j+1,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j+2,k,s) - xyzf_VarNl(i,j-1,k,s) ) &
                &          )                                                          &
                &    + xqz_VelYNl(i,j-1,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j-1,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j+1,k,s) - xyzf_VarNl(i,j-2,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dy                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+2,s) - xyzf_VarNl(i,j,k-1,s) ) &
                &          )                                                          &
                &    + xyr_VelZNl(i,j,k-1)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j,k-1,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k-2,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dz                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)   * xyrf_DVarBZDz(i,j,k,s)                   &
                &    + xyr_VelZNl(i,j,k-1) * xyrf_DVarBZDz(i,j,k-1,s)                 & 
                &   ) * 5.0d-1 
            end do
          end do
        end do
      end do
      
      xyzf_Adv(imin:imin+1,:,:,:) = 0.0d0
      xyzf_Adv(imax-1:imax,:,:,:) = 0.0d0
      xyzf_Adv(:,jmin:jmin+1,:,:) = 0.0d0
      xyzf_Adv(:,jmax-1:jmax,:,:) = 0.0d0
      xyzf_Adv(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_Adv(:,:,kmax-1:kmax,:) = 0.0d0
      
      ! $B?tCM3H;6(B: 2 $B<!@:EYCf?4:9J,(B
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = jmin + 2, jmax - 2
            do i = imin + 2, imax - 2
              
              xyzf_nDiff(i,j,k,s) =                            &
                & - (                                          &
                &       xyzf_VarBl(i+2,j,k,s)                  &
                &     + xyzf_VarBl(i-2,j,k,s)                  &
                &     - xyzf_VarBl(i+1,j,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i-1,j,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i  ,j,k,s) * 6.0d0          &
                &   ) * NuHh / ( dx ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j+2,k,s)                  &
                &     + xyzf_VarBl(i,j-2,k,s)                  &
                &     - xyzf_VarBl(i,j+1,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j-1,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j  ,k,s) * 6.0d0          &
                &   ) * NuHh / ( dy ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j,k+2,s)                  &
                &     + xyzf_VarBl(i,j,k-2,s)                  &
                &     - xyzf_VarBl(i,j,k+1,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j,k-1,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j,k  ,s) * 6.0d0          &
                &   ) * NuVh / ( dz ** 4.0d0 )
            end do
          end do
        end do
      end do
      
      xyzf_nDiff(imin:imin+1,:,:,:) = 0.0d0
      xyzf_nDiff(imax-1:imax,:,:,:) = 0.0d0
      xyzf_nDiff(:,jmin:jmin+1,:,:) = 0.0d0
      xyzf_nDiff(:,jmax-1:jmax,:,:) = 0.0d0
      xyzf_nDiff(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_nDiff(:,:,kmax-1:kmax,:) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyzf
  

    subroutine BuoyancyLong_xyr
      
      use composition, only: GasNum,       &! 
        &                    IdxG,         &!
        &                    MolWtWet       ! $B<>=a@.J,$NJ,;RNL(B
      use constants,   only: MolWtDry,     &! $B4%Ag@.J,$NJ,;RNL(B
        &                    Grav           ! $B=ENO2CB.EY(B
      use basicset,    only: xyr_QMixBZ,         &
        &                    xyr_QMixBZPerMolWt, &
        &                    xyz_PTempBZ
      
      implicit none
      
      real(DP)              :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)              :: tmp1(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp2(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp3(imin:imax,jmin:jmax,kmin:kmax)
      integer               :: i, j, k, f, n
      
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do
      
      ! Buoyancy due to temperature disturbunce
      !
      do k = kmin, kmax - 1
        do j = jmin, jmax
          do i = imin, imax
            
            xyr_BuoyT(i,j,k) =                                  &
              & Grav                                            &
              & * (                                             &
              &     xyz_PTempNl(i,j,k+1) / xyz_PTempBZ(i,j,k+1) &
              &   + xyz_PTempNl(i,j,k)   / xyz_PTempBZ(i,j,k)   &
              &   ) * 5.0d-1
            
          end do
        end do
      end do
      
      xyr_BuoyT(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to molecular weight
      !
      tmp1 = sum(xyzf_QMixPerMolWt, 4) 
      tmp2 = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax - 1
        do j = jmin, jmax 
          do i = imin, imax 
            
            xyr_BuoyM(i,j,k) =                                       &
              & + Grav                                               &
              &   * ( tmp1(i,j,k+1) + tmp1(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) ) &
              & - Grav                                               &
              &   * ( tmp2(i,j,k+1) + tmp2(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 + xyr_QmixBZ(i,j,k) ) 
            
          end do
        end do
      end do
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      tmp3 = sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) 
      
      do k = kmin, kmax - 1
        do j = jmin, jmax 
          do i = imin, imax 
            
            xyr_BuoyD(i,j,k) =                                &
              & - Grav                                        &
              &   * ( tmp3(i,j,k+1) + tmp3(i,j,k) ) * 5.0d-1  &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do
      
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr
    
  end subroutine Dynamics3D_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics2D_Long_forcing( &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
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
    ! $B0\N.7W;;(B (2D $BHG(B)
    !
    !   $B0\N.(B: 4 $B<!Cf1{:9J,(B
    !   $B?tCM3H;6(B (4 $B3,(B): 2 $B<!Cf1{:9J,(B
    !
    ! $B%j!<%W%U%m%C%0$G(B, $B0\N.$rCf1{:9J,$G7W;;$9$k$?$a$K(B, 
    ! $B?tCM3H;69`$rDI2C$7$F$$$k(B. 
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,    only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x $BJ}8~$NG[Ns%5%$%:(B
      &                     jmin, jmax,       &! y $BJ}8~$NG[Ns%5%$%:(B
      &                     kmin, kmax,       &! z $BJ}8~$NG[Ns%5%$%:(B
      &                     nx, ny, nz,       &! $BJ*M}NN0h$N%5%$%:(B
      &                     ncmax              ! $BJ*<A?t(B
    use timeset,     only : TimeN
    use axesset,     only : dx, dz             ! $B3J;R4V3V(B
    use composition, only : SpcWetSymbol

    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
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

    real(DP)             :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: pyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_nDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    integer              :: f


    !--------------------------------------------------------------------
    ! $BMpN.3H;678?t(B

    ! $B0\N.$*$h$S?tCM3H;6(B
    !
    call AdvC4_nDiff_xyz( xyz_KmBl, xyz_KmNl, xyr_Dummy ) !(IN)
    
    ! tendency $B$N99?7(B
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_nDiff + xyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))
    
    !--------------------------------------------------------------------
    ! $B290L(B

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyz( xyz_PTempBl, xyz_PTempNl, xyr_DPTempBZDz ) !(IN)
 
    ! tendency $B$N99?7(B
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_nDiff + xyz_Adv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',   xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff',  xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! $B%(%/%9%J!<4X?t(B

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyz( xyz_ExnerBl, xyz_ExnerNl, xyr_DExnerBZDz ) !(IN)

    ! tendency $B$N99?7(B
    !
    xyz_nDiff = 0.0d0
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_nDiff + xyz_Adv ) * FactorDExnerDtAdv
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! $B:.9gHf(B
    ! 

    ! $B0\N.$K$D$$$F$O(B, $B4pK\>l$b9MN8$9$k(B.
    !
    call AdvC4_nDiff_xyzf( xyzf_QMixBl, xyzf_QMixNl, xyrf_DQMixBZDz ) !(IN)

    ! $BMn2<9`(B
    !
    call QMixFall

    ! tendency $B$N99?7(B    
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_nDiff + xyzf_Adv + xyzf_Fall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_nDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', &
        & xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------------------------------------------
    ! VelX, VelY, VelZ
    ! 

    ! $B0\N.9`!&?tCM3H;69`$r$^$H$a$F7W;;(B
    !
    call AdvC4_nDiff_pyz_xyr

    ! tendency of VelX
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_nDiff(1:nx,1:ny,1:nz))

    ! tendency of VelY
    !
    xqz_DVelYDtNl = 0.0d0

    ! Buoyancy 
    ! 
    call BuoyancyLong_xyr

    ! tendency of VelZ
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
      use average,     only : xyr_xyz
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyzf_QMixBZ, &!$B4pK\>l$N:.9gHf(B
        &                     xyz_DensBZ
      use differentiate_center4, &
        &              only : xyz_dz_xyr
      
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
          &  xyr_xyz (                                                 &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall      &
          &  )

        ! $B>eC<$N%U%i%C%/%9$O%<%m(B
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! $B1+N3Mn2<$K$h$k;~4VJQ2=(B (DelTime $B$r$+$1$F$"$k(B)
        !        
        xyzf_Fall(:,:,:,iR) =                                          &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall

    
    subroutine AdvC4_nDiff_pyz_xyr

      implicit none

      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            pyz_Adv(i,j,k) =                                                   &
              & - pyz_VelXNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   pyz_VelXNl(i+1,j,k) - pyz_VelXNl(i-1,j,k) )   &
              &     - fct2 * (   pyz_VelXNl(i+2,j,k) + pyz_VelXNl(i+1,j,k)     &
              &                - pyz_VelXNl(i-1,j,k) - pyz_VelXNl(i-2,j,k) )   &
              &     ) * 5.0d-1 / dx                                            &
              & - (                                                            &
              &   + ( xyr_VelZNl(i+1,j,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+2) - pyz_VelXNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i+1,j,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j,k-1) )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
            
          end do
        end do
      end do
      
      pyz_Adv(imin:imin+1,:,:) = 0.0d0
      pyz_Adv(imax-1:imax,:,:) = 0.0d0
      pyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      pyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyr_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j,k+1) + pyz_VelXNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i+2,j,k) - xyr_VelZNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j,k+1) + pyz_VelXNl(i-1,j,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i-1,j,k) )   &
              &       - fct2 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - xyr_VelZNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xyr_VelZNl(i,j,k+1) - xyr_VelZNl(i,j,k-1) )   &
              &     - fct2 * (   xyr_VelZNl(i,j,k+2) + xyr_VelZNl(i,j,k+1)     &
              &                - xyr_VelZNl(i,j,k-1) - xyr_VelZNl(i,j,k-2) )   &
              &     ) * 5.0d-1 / dz
          end do
        end do
      end do
      
      xyr_Adv(imin:imin+1,:,:) = 0.0d0
      xyr_Adv(imax-1:imax,:,:) = 0.0d0
      xyr_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyr_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            pyz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + pyz_VelXBl(i+2,j,k)                  &
              &     + pyz_VelXBl(i-2,j,k)                  &
              &     - pyz_VelXBl(i+1,j,k) * 4.0d0          &
              &     - pyz_VelXBl(i-1,j,k) * 4.0d0          &
              &     + pyz_VelXBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        & 
              &     + pyz_VelXBl(i,j,k+2)                  &
              &     + pyz_VelXBl(i,j,k-2)                  &
              &     - pyz_VelXBl(i,j,k+1) * 4.0d0          &
              &     - pyz_VelXBl(i,j,k-1) * 4.0d0          &
              &     + pyz_VelXBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
            
          end do
        end do
      end do
      
      pyz_nDiff(imin:imin+1,:,:) = 0.0d0
      pyz_nDiff(imax-1:imax,:,:) = 0.0d0
      pyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      pyz_nDiff(:,:,kmax-1:kmax) = 0.0d0

      do k = kmin + 2, kmax - 2
        do j = 1, ny          
          do i = imin + 2, imax - 2
            
            xyr_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xyr_VelZBl(i+2,j,k)                  &
              &     + xyr_VelZBl(i-2,j,k)                  &
              &     - xyr_VelZBl(i+1,j,k) * 4.0d0          &
              &     - xyr_VelZBl(i-1,j,k) * 4.0d0          &
              &     + xyr_VelZBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j,k+2)                  &
              &     + xyr_VelZBl(i,j,k-2)                  &
              &     - xyr_VelZBl(i,j,k+1) * 4.0d0          &
              &     - xyr_VelZBl(i,j,k-1) * 4.0d0          &
              &     + xyr_VelZBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyr_nDiff(imin:imin+1,:,:) = 0.0d0
      xyr_nDiff(imax-1:imax,:,:) = 0.0d0
      xyr_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyr_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_pyz_xyr
    
    
    subroutine AdvC4_nDiff_xyz( xyz_VarBl, xyz_VarNl, xyr_DVarBZDz )
      
      implicit none
      
      real(DP), intent(in)  :: xyz_VarBl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyz_VarNl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyr_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax)
      
      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0

      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyz_Adv(i,j,k) =                                                  &
              & - (                                                           &
              &      pyz_VelXNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i+2,j,k) - xyz_VarNl(i-1,j,k) ) &
              &          )                                                    &
              &    + pyz_VelXNl(i-1,j,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i-1,j,k) ) &
              &          - fct2 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i-2,j,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dx                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+2) - xyz_VarNl(i,j,k-1) ) &
              &          )                                                    &
              &    + xyr_VelZNl(i,j,k-1)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j,k-1) ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k-2) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dz                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)   * xyr_DVarBZDz(i,j,k)                &
              &    + xyr_VelZNl(i,j,k-1) * xyr_DVarBZDz(i,j,k-1)              &
              &   ) * 5.0d-1
          end do
        end do
      end do
      
      xyz_Adv(imin:imin+1,:,:) = 0.0d0
      xyz_Adv(imax-1:imax,:,:) = 0.0d0
      xyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      ! 4 $B<!$N?tCM3H;6(B: 2 $B<!@:EYCf?4:9J,(B
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyz_nDiff(i,j,k) =                            &
              & - (                                       &
              &     + xyz_VarBl(i+2,j,k)                  &
              &     + xyz_VarBl(i-2,j,k)                  &
              &     - xyz_VarBl(i+1,j,k) * 4.0d0          &
              &     - xyz_VarBl(i-1,j,k) * 4.0d0          &
              &     + xyz_VarBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j,k+2)                  &
              &     + xyz_VarBl(i,j,k-2)                  &
              &     - xyz_VarBl(i,j,k+1) * 4.0d0          &
              &     - xyz_VarBl(i,j,k-1) * 4.0d0          &
              &     + xyz_VarBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyz_nDiff(imin:imin+1,:,:) = 0.0d0
      xyz_nDiff(imax-1:imax,:,:) = 0.0d0
      xyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyz
    
    
    subroutine AdvC4_nDiff_xyzf( xyzf_VarBl, xyzf_VarNl, xyrf_DVarBZDz ) 

      implicit none
      
      real(DP), intent(in)  :: xyzf_VarBl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyzf_VarNl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyrf_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      
      real(DP)              :: fct1, fct2
      integer               :: s, i, j, k
      
      ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! $B0\N.9`$N7W;;(B. $B0\N.9`(B: 4 $B<!@:EYCf?4:9J,(B
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = 1, ny
            do i = imin + 2, imax - 2
              
              xyzf_Adv(i,j,k,s) =                                                     &
                & - (                                                                 &
                &      pyz_VelXNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i+2,j,k,s) - xyzf_VarNl(i-1,j,k,s) ) &
                &          )                                                          &
                &    + pyz_VelXNl(i-1,j,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i-1,j,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i-2,j,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dx                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+2,s) - xyzf_VarNl(i,j,k-1,s) ) &
                &          )                                                          &
                &    + xyr_VelZNl(i,j,k-1)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j,k-1,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k-2,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dz                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)   * xyrf_DVarBZDz(i,j,k,s)                   &
                &    + xyr_VelZNl(i,j,k-1) * xyrf_DVarBZDz(i,j,k-1,s)                 &
                &   ) * 5.0d-1 
            end do
          end do
        end do
      end do
      
      xyzf_Adv(imin:imin+1,:,:,:) = 0.0d0
      xyzf_Adv(imax-1:imax,:,:,:) = 0.0d0
      xyzf_Adv(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_Adv(:,:,kmax-1:kmax,:) = 0.0d0
      
      ! $B?tCM3H;6(B: 2 $B<!@:EYCf?4:9J,(B
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = 1, ny
            do i = imin + 2, imax - 2
              
              xyzf_nDiff(i,j,k,s) =                            &
                & - (                                          &
                &       xyzf_VarBl(i+2,j,k,s)                  &
                &     + xyzf_VarBl(i-2,j,k,s)                  &
                &     - xyzf_VarBl(i+1,j,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i-1,j,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i  ,j,k,s) * 6.0d0          &
                &   ) * NuHh / ( dx ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j,k+2,s)                  &
                &     + xyzf_VarBl(i,j,k-2,s)                  &
                &     - xyzf_VarBl(i,j,k+1,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j,k-1,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j,k  ,s) * 6.0d0          &
                &   ) * NuVh / ( dz ** 4.0d0 )
            end do
          end do
        end do
      end do
      
      xyzf_nDiff(imin:imin+1,:,:,:) = 0.0d0
      xyzf_nDiff(imax-1:imax,:,:,:) = 0.0d0
      xyzf_nDiff(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_nDiff(:,:,kmax-1:kmax,:) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyzf
  

    subroutine BuoyancyLong_xyr
      
      use composition, only: GasNum,       &! 
        &                    IdxG,         &!
        &                    MolWtWet       ! $B<>=a@.J,$NJ,;RNL(B
      use constants,   only: MolWtDry,     &! $B4%Ag@.J,$NJ,;RNL(B
        &                    Grav           ! $B=ENO2CB.EY(B
      use basicset,    only: xyr_QMixBZ,         &
        &                    xyr_QMixBZPerMolWt, &
        &                    xyz_PTempBZ
      
      implicit none
      
      real(DP)              :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)              :: tmp1(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp2(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp3(imin:imax,jmin:jmax,kmin:kmax)
      integer               :: i, j, k, f, n
      
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do
      
      ! Buoyancy due to temperature disturbunce
      !
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax
            
            xyr_BuoyT(i,j,k) =                                  &
              & Grav                                            &
              & * (                                             &
              &     xyz_PTempNl(i,j,k+1) / xyz_PTempBZ(i,j,k+1) &
              &   + xyz_PTempNl(i,j,k)   / xyz_PTempBZ(i,j,k)   &
              &   ) * 5.0d-1
            
          end do
        end do
      end do
      
      xyr_BuoyT(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to molecular weight
      !
      tmp1 = sum(xyzf_QMixPerMolWt, 4) 
      tmp2 = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax 
            
            xyr_BuoyM(i,j,k) =                                       &
              & + Grav                                               &
              &   * ( tmp1(i,j,k+1) + tmp1(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) ) &
              & - Grav                                               &
              &   * ( tmp2(i,j,k+1) + tmp2(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 + xyr_QmixBZ(i,j,k) ) 
            
          end do
        end do
      end do
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      tmp3 = sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) 
      
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax 
            
            xyr_BuoyD(i,j,k) =                                &
              & - Grav                                        &
              &   * ( tmp3(i,j,k+1) + tmp3(i,j,k) ) * 5.0d-1  &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do
      
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr
    
  end subroutine Dynamics2D_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                           !$BG[Ns%5%$%:(B
      &                    FlagCalc3D      !$B%U%i%0(B

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

    if ( FlagCalc3D ) then 
      call Dynamics3D_Short_forcing(  &
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
    else
      call Dynamics2D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    end if

  end subroutine Dynamics_Short_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics3D_Short_forcing(  &
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
    ! $BNO3X%3%"(B ($BC;$$;~4V%9%F%C%W(B) 3D $BHG(B
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gtool_historyauto, only: HistoryAutoPut 
    
    use gridset, only: &
      &                 imin,            &! x $BJ}8~$NG[Ns$N2<8B(B
      &                 imax,            &! x $BJ}8~$NG[Ns$N>e8B(B
      &                 jmin,            &! y $BJ}8~$NG[Ns$N2<8B(B
      &                 jmax,            &! y $BJ}8~$NG[Ns$N>e8B(B
      &                 kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                 kmax,            &! z $BJ}8~$NG[Ns$N>e8B(B
      &                 nx,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 ny,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 nz                ! y $BJ}8~$NJ*M}NN0h$N>e8B(B
    use constants,only: CpDry, CvDry, GasRDry ! $B4%Ag@.J,$NHfG.(B
    use timeset, only:  DelTimeShort,  TimeN
    use axesset, only:  dx, dy, dz        ! $B3J;R4V3V(B
    use basicset, only: xyz_VelSW,      &!$B4pK\>l$N2;B.(B 
      &                 xyz_VPTempBZ,      &!$B4pK\>l$N290L(B
      &                 pyz_VPTempBZ,      &!$B4pK\>l$N290L(B
      &                 xqz_VPTempBZ,      &!$B4pK\>l$N290L(B
      &                 xyr_VPTempBZ        !$B4pK\>l$N290L(B
    use setmargin,only: SetMargin_xyzf, SetMargin_xyz, &
      &                 SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

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

    real(DP) :: pyz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_PGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_SWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------     
    ! initialize: Divergence of velocity
    !
    call VelDivC2
    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))

    !------------------------------------------------------------
    ! VelX, VelY
    !  $B?eJ?J}8~$OM[2rK!$G2r$/(B. 
    !
    call PGrad_HE

    ! tendency 
    !
    pyz_DVelXDtNs = pyz_PGrad + pyz_SWF
    xqz_DVelYDtNs = xqz_PGrad + xqz_SWF

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * (pyz_DVelXDtNl + pyz_DVelXDtNs)
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * (xqz_DVelYDtNl + xqz_DVelYDtNs)

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)
    
    !------------------------------------------------------------
    ! Exner function
    !  $B@QJ,CM$rJV$9$3$H$KCm0U(B. 
    !

    ! $BC;$$;~4V%9%F%C%W$GI>2A$9$k05NO$N<0$N(B tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency $B$N9g7W(B
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    !$B%(%/%9%J!<4X?t$N7W;;(B
    ! 
    call Exner_HEVI

    !------------------------------------------------------------
    ! VelZ
    !
    call PGrad_VI

    ! tendency
    !
    xyr_DVelZDtNs = xyr_PGrad + xyr_SWF

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNl + xyr_DVelZDtNs)

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  contains

    subroutine VelDivC2
      
      implicit none
      integer              :: i, j, k
      
      do  k = kmin + 1, kmax 
        do j = jmin + 1, jmax 
          do i = imin + 1, imax 
            xyz_VelDivNs(i,j,k) =        &
              & + (                      &
              &     pyz_VelXNs(i,j,k)      &
              &   - pyz_VelXNs(i-1,j,k)    &
              &   ) / dx                 &
              & + (                      &
              &     xqz_VelYNs(i,j,k)      &
              &   - xqz_VelYNs(i,j-1,k)    &
              &   ) / dy                 &
              & + (                      &
              &     xyr_VelZNs(i,j,k)      &
              &   - xyr_VelZNs(i,j,k-1)    &
              &   ) / dz
          end do
        end do
      end do
      
      xyz_VelDivNs(imin,:,:) = 0.0d0 
      xyz_VelDivNs(:,jmin,:) = 0.0d0  
      xyz_VelDivNs(:,:,kmin) = 0.0d0 
      
    end subroutine VelDivC2


    subroutine PGrad_HE
    
      implicit none
      integer              :: i, j, k
      
      !------------------------------------------------------------------
      ! X $BJ}8~(B
      
      do k = kmin, kmax
        do j = jmin, jmax
          do i = imin, imax - 1

            ! $B2;GH8:?j9`(B
            !            
            pyz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i+1,j,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dx            
            
            ! $B05NO79EYNO(B
            !
            pyz_PGrad(i,j,k) =                &
              & - CpDry * pyz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i+1,j,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) / dx                                
          end do
        end do
      end do
      
      ! $B7jKd$a(B
      !
      pyz_SWF(imax,:,:)   = 0.0d0 
      pyz_PGrad(imax,:,:) = 0.0d0 
      

      !------------------------------------------------------------------
      ! Y $BJ}8~(B
      
      do k = kmin, kmax
        do j = jmin, jmax - 1
          do i = imin, imax

            ! $B2;GH8:?j9`(B
            !            
            xqz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j+1,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dy
               
            ! $B05NO79EYNO(B
            !             
            xqz_PGrad(i,j,k) =                &
              & - CpDry * xqz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i,j+1,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) /dy                     
            
          end do
        end do
      end do
      
      ! $B7jKd$a(B
      !
      xqz_SWF(:,jmax,:) = 0.0d0 
      xqz_PGrad(:,jmax,:) = 0.0d0 
      
    end subroutine PGrad_HE


    subroutine PGrad_VI
    
      implicit none
      integer               :: i, j, k

      do k = kmin, kmax - 1
        do j = jmin, jmax
          do i = imin, imax

            ! $B2;GH8:?j9`(B
            !            
            xyr_SWF(i,j,k) =                  &
              & + AlphaV                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j,k+1)   & 
              &     - xyz_VelDivNs(i,j,k)     & 
              &     ) / dz
            
            ! $B05NO79EYNO(B
            !
            xyr_PGrad(i,j,k) =                 &
              & - CpDry * xyr_VPTempBZ(i,j,k)  &
              &   * (                          &
              &       beta                     &
              &       * (                      &
              &           xyz_ExnerAs(i,j,k+1) &
              &         - xyz_ExnerAs(i,j,k)   &
              &         )                      &
              &     + (1.0d0 - beta)           &
              &       * (                      &
              &           xyz_ExnerNs(i,j,k+1) &
              &         - xyz_ExnerNs(i,j,k)   &
              &         )                      &
              &     ) / dz
            
          end do
        end do
      end do
      
      xyr_PGrad(:,:,kmax) = 0.0d0
      xyr_SWF(:,:,kmax)   = 0.0d0
      
    end subroutine PGrad_VI


    subroutine Exner_HEVI
      !
      !$B1"2rK!$rMQ$$$?%(%/%9%J!<4X?t$N7W;;(B. 
      !

      !$B0EL[$N7?@k8@6X;_(B
      implicit none

      !$B:n6HJQ?tDj5A(B
      real(DP)               :: D1(1:nx,1:ny,1:nz)  
      real(DP)               :: D(nx*ny,nz)
      real(DP)               :: E(1:nx,1:ny,0:nz)
      real(DP)               :: F(1:nx,1:ny,1:nz)
      real(DP)               :: F0(1:nx,1:ny,kmin:kmax-1)  
      real(DP)               :: dt ! $BC;$$;~4V3J;R4V3V(B
      integer                :: i, j, k
      
      real(DP)               :: X(M, N)     !$BDj?t(B/$B2r9TNs(B
      real(DP)               :: TX(N, M)    !$B2r9TNs$rE>CV$7$?$b$N(B
      integer                :: NRHS        
      integer                :: INFO
      integer                :: LDB
      character(1),parameter :: TRANS = 'N'
            
      ! Initialize
      !
      dt = DelTimeShort

      !---------------------------------------------------------------
      !$B9TNs7W;;$N$?$a$N78?t(B

      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 0:nz
      !  D(:,:,1) $B$r5a$a$k;~$K(B D(:,:,0) $B$NCM$,I,MW$K$J$k$?$a(B. 
      !
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            
            E(i,j,k) =                            &
              & - ( 1.0d0 - beta )                &
              &   * (                             &
              &     + xyz_ExnerNs(i,j,k+1)        & !
              &     - xyz_ExnerNs(i,j,k)          & ! xyz => xyr
              &     ) / dz                        &
              & + (                               &
              &   + AlphaV                        &
              &     * (                           &
              &       + xyz_VelDivNs(i,j,k+1)     & !
              &       - xyz_VelDivNs(i,j,k)       & ! xyz => xyr
              &       ) / dz                      &
              &   + xyr_DVelZDtNl(i,j,k)          &
              &   ) / xyr_CpVPTempBZ(i,j,k)  
            
          end do
        end do
      end do
      
      ! $BHoHyJ,4X?t(B
      !   $BG[Ns(B F0 $B$NE:;z$NHO0O$O(B, 1:nx, 1:ny, kmin:kmax
      !   $BG[Ns(B F $B$r5a$a$k:]$K(B F0 $B$r(B z $BJ}8~$KHyJ,$9$k$?$a(B. 
      !
      do k = kmin, kmax-1
        do j = 1, ny
          do i = 1, nx
            
            F0(i,j,k)  =                            &
              & + xyr_DensVPTempBZ(i,j,k)           &
              &   * (                               &
              &     + xyr_VelZNs(i,j,k)             & 
              &     - xyr_CpVPTempBZ(i,j,k)         &
              &       * (1.0d0 - beta)              &
              &       * (                           &
              &           xyz_ExnerNs(i,j,k+1)      &
              &         - xyz_ExnerNs(i,j,k)        &
              &         ) / dz * dt                 &
              &     + AlphaV                        &
              &       * (                           &
              &           xyz_VelDivNs(i,j,k+1)     &
              &         - xyz_VelDivNs(i,j,k)       &
              &         ) / dz * dt                 &
              &     + xyr_DVelZDtNl(i,j,k) * dt     &
              &     )
          end do
        end do
      end do
      
      !$B9TNs7W;;$N$?$a$N78?t(B
      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            F(i,j,k) = &
              & - beta * xyz_F1BZ(i,j,k) * dt   &
              &   * (                           &
              &       F0(i,j,k)                 &
              &     - F0(i,j,k-1)               &
              &     ) / dz                      &
              & + xyz_DExnerDtNl(i,j,k) * dt 
            
          end do
        end do
      end do
    
      !$B9TNs7W;;$N$?$a$N78?t(B
      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz 
        do j = 1, ny
          do i = 1, nx
            
            D1(i,j,k) =                                 &
              & + xyz_ExnerNs(i,j,k)                    &
              & - (1.0d0 - beta)                        &
              &   * xyz_F1BZ(i,j,k) * dt                &
              &   * (                                   &
              &       xyr_DensVPTempBZ(i,j,k)           &
              &       * xyr_VelZNs(i,j,k)               &
              &     - xyr_DensVPTempBZ(i,j,k-1)         &
              &       * xyr_VelZNs(i,j,k-1)             &
              &     ) / dz                              &
              & - (xyz_VelSW(i,j,k) ** 2.0d0) * dt      &
              &   / (CpDry * xyz_VPTempBZ(i,j,k))       &
              &   * (                                   &
              &     + (                                 &
              &         pyz_VelXAs(i,j,k)               &
              &       - pyz_VelXAs(i-1,j,k)             &
              &       ) / dx                            &
              &     + (                                 &
              &         xqz_VelYAs(i,j,k)               &
              &       - xqz_VelYAs(i,j-1,k)             &
              &       ) / dy                            &
              &     )                                   &
              & + F(i,j,k)
          
          end do
        end do
      end do
      
      ! $B9TNs7W;;$N$?$a$N78?t(B
      !
      do j = 1, ny
        do i = 1, nx

          D1(i,j,1) =                                    &
            & + D1(i,j,1)                                &
            & - beta * xyz_F1BZ(i,j,1) * (dt ** 2.0d0)   &
            &   * xyr_CpDensVPTemp2BZ(i,j,0)             &
            &   * E(i,j,0)                               &
            &   / dz
          
          D1(i,j,nz) =                                   &
            & + D1(i,j,nz)                               &
            & + beta * xyz_F1BZ(i,j,nz) * (dt ** 2.0d0)  &
            &   * xyr_CpDensVPTemp2BZ(i,j,nz)            &
            &   * E(i,j,nz)                              &
            &   / dz
        end do
      end do
      
      !-----------------------------------------------------------
      !$BO"N)0l<!J}Dx<0$N2r$r5a$a$k(B
      
      !$BJQ?t$N=i4|2=(B
      !
      NRHS = M
      INFO = 0
      LDB  = N
      
      ! LAPACK $B$N;EMM$K9g$o$;$FJQ7A(B 
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            D(i + nx * (j - 1), k) =  D1(i,j,k)
          end do
        end do
      end do
      
      TX = transpose( D )
      
      !$B2r9TNs$N7W;;(B. LAPACK $B$r;HMQ(B. 
      !
      call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)
      
      !$B2r$N%3%s%G%#%7%g%s$r%A%'%C%/(B. 
      !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

      !$BLa$jCM$r=PNO(B
      !
      X = transpose( TX )
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyz_ExnerAs(i,j,k) = X(i + nx * (j - 1 ), k)
          end do
        end do
      end do

      xyz_ExnerAs(imin:0,:,:) = 0.0d0
      xyz_ExnerAs(:,jmin:0,:) = 0.0d0
      xyz_ExnerAs(:,:,kmin:0) = 0.0d0
      xyz_ExnerAs(nx+1:imax,:,:) = 0.0d0
      xyz_ExnerAs(:,ny+1:jmax,:) = 0.0d0
      xyz_ExnerAs(:,:,nz+1:kmax) = 0.0d0
      
      ! $B$N$jBe$KCM$rF~$l$k(B
      !
      call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
      
    end subroutine Exner_HEVI
    
  end subroutine Dynamics3D_Short_forcing
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  subroutine Dynamics2D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! $BNO3X%3%"(B ($BC;$$;~4V%9%F%C%W(B) 2D $BHG(B
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gtool_historyauto, only: HistoryAutoPut 
    use gridset, only: &
      &                 imin,            &! x $BJ}8~$NG[Ns$N2<8B(B
      &                 imax,            &! x $BJ}8~$NG[Ns$N>e8B(B
      &                 jmin,            &! y $BJ}8~$NG[Ns$N2<8B(B
      &                 jmax,            &! y $BJ}8~$NG[Ns$N>e8B(B
      &                 kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                 kmax,            &! z $BJ}8~$NG[Ns$N>e8B(B
      &                 nx,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 ny,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 nz              ! y $BJ}8~$NJ*M}NN0h$N>e8B(B
    use constants,only: CpDry, CvDry, GasRDry ! $B4%Ag@.J,$NHfG.(B
    use timeset, only:  DelTimeShort, TimeN
    use axesset, only:  dx,  dz        ! $B3J;R4V3V(B
    use basicset, only: xyz_VelSW,      &!$B4pK\>l$N2;B.(B 
      &                 xyz_VPTempBZ,      &!$B4pK\>l$N290L(B
      &                 pyz_VPTempBZ,      &!$B4pK\>l$N290L(B
      &                 xyr_VPTempBZ      !$B4pK\>l$N290L(B
    use setmargin,only: SetMargin_xyzf, SetMargin_xyz, &
      &                 SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

    implicit none

    real(DP), intent(in)     :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout)  :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_PGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_SWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------
    ! initialize: Divergence of velocity
    !
    call VelDivC2

    !------------------------------------------------------------
    ! VelX, VelY
    !  $B?eJ?J}8~$OM[2rK!$G2r$/(B. 
    !
    call PGrad_HE

    ! tendency 
    !
    pyz_DVelXDtNs = pyz_PGrad + pyz_SWF

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * (pyz_DVelXDtNl + pyz_DVelXDtNs)

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)

    ! Y $BJ}8~$K$OCM$,%<%m(B
    xqz_VelYAs = 0.0d0

    
    !------------------------------------------------------------
    ! Exner function
    !  $B@QJ,CM$rJV$9$3$H$KCm0U(B. 
    !

    ! $BC;$$;~4V%9%F%C%W$GI>2A$9$k05NO$N<0$N(B tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency $B$N9g7W(B
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    ! $B%(%/%9%J!<4X?t$N7W;;(B
    !
    call Exner_HEVI

    !------------------------------------------------------------
    ! VelZ
    !
    call PGrad_VI

    ! tendency
    !
    xyr_DVelZDtNs = xyr_PGrad + xyr_SWF

    ! $BCM$NJ]4I(B
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNl + xyr_DVelZDtNs)

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  contains

    subroutine VelDivC2
      
      implicit none
      integer              :: i, j, k
      
      do  k = kmin + 1, kmax 
        do j = 1, ny
          do i = imin + 1, imax 
            xyz_VelDivNs(i,j,k) =          &
              & + (                        &
              &     pyz_VelXNs(i,j,k)      &
              &   - pyz_VelXNs(i-1,j,k)    &
              &   ) / dx                   &
              & + (                        &
              &     xyr_VelZNs(i,j,k)      &
              &   - xyr_VelZNs(i,j,k-1)    &
              &   ) / dz
          end do
        end do
      end do
      
      xyz_VelDivNs(imin,:,:) = 0.0d0 
      xyz_VelDivNs(:,:,kmin) = 0.0d0 
      
    end subroutine VelDivC2


    subroutine PGrad_HE
      
      implicit none
      integer              :: i, j, k
      
      !------------------------------------------------------------------
      ! X $BJ}8~(B
      
      do k = kmin, kmax
        do j = 1, ny          
          do i = imin, imax - 1

            ! $B2;GH8:?j9`(B
            !            
            pyz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i+1,j,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dx            
            
            ! $B05NO79EYNO(B
            !
            pyz_PGrad(i,j,k) =                &
              & - CpDry * pyz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i+1,j,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) / dx                                
            
          end do
        end do
      end do
      
      ! $B7jKd$a(B
      !
      pyz_SWF(imax,:,:)   = 0.0d0 
      pyz_PGrad(imax,:,:) = 0.0d0 
      
    end subroutine PGrad_HE


    subroutine PGrad_VI
    
      implicit none
      integer               :: i, j, k

      do k = kmin, kmax - 1
        do j = 1, ny          
          do i = imin, imax

            ! $B2;GH8:?j9`(B
            !            
            xyr_SWF(i,j,k) =                  &
              & + AlphaV                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j,k+1)   & 
              &     - xyz_VelDivNs(i,j,k)     & 
              &     ) / dz
            
            ! $B05NO79EYNO(B
            !
            xyr_PGrad(i,j,k) =                 &
              & - CpDry * xyr_VPTempBZ(i,j,k)  &
              &   * (                          &
              &       beta                     &
              &       * (                      &
              &           xyz_ExnerAs(i,j,k+1) &
              &         - xyz_ExnerAs(i,j,k)   &
              &         )                      &
              &     + (1.0d0 - beta)           &
              &       * (                      &
              &           xyz_ExnerNs(i,j,k+1) &
              &         - xyz_ExnerNs(i,j,k)   &
              &         )                      &
              &     ) / dz
            
          end do
        end do
      end do
      
      xyr_PGrad(:,:,kmax) = 0.0d0
      xyr_SWF(:,:,kmax)   = 0.0d0
      
    end subroutine PGrad_VI


    subroutine Exner_HEVI
      !
      !$B1"2rK!$rMQ$$$?%(%/%9%J!<4X?t$N7W;;(B. 
      !

      !$B0EL[$N7?@k8@6X;_(B
      implicit none

      !$B:n6HJQ?tDj5A(B
      real(DP)               :: D1(1:nx,1:ny,1:nz)  
      real(DP)               :: D(nx*ny,nz)
      real(DP)               :: E(1:nx,1:ny,0:nz)
      real(DP)               :: F(1:nx,1:ny,1:nz)
      real(DP)               :: F0(1:nx,1:ny,kmin:kmax-1)  
      real(DP)               :: dt ! $BC;$$;~4V3J;R4V3V(B
      integer                :: i, j, k
      
      real(DP)               :: X(M, N)     !$BDj?t(B/$B2r9TNs(B
      real(DP)               :: TX(N, M)    !$B2r9TNs$rE>CV$7$?$b$N(B
      integer                :: NRHS        
      integer                :: INFO
      integer                :: LDB
      character(1),parameter :: TRANS = 'N'
      
      
      ! Initialize
      !
      dt = DelTimeShort

      !---------------------------------------------------------------
      !$B9TNs7W;;$N$?$a$N78?t(B

      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 0:nz
      !  D(:,:,1) $B$r5a$a$k;~$K(B D(:,:,0) $B$NCM$,I,MW$K$J$k$?$a(B. 
      !
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            
            E(i,j,k) =                            &
              & - ( 1.0d0 - beta )                &
              &   * (                             &
              &     + xyz_ExnerNs(i,j,k+1)        & !
              &     - xyz_ExnerNs(i,j,k)          & ! xyz => xyr
              &     ) / dz                        &
              & + (                               &
              &   + AlphaV                        &
              &     * (                           &
              &       + xyz_VelDivNs(i,j,k+1)     & !
              &       - xyz_VelDivNs(i,j,k)       & ! xyz => xyr
              &       ) / dz                      &
              &   + xyr_DVelZDtNl(i,j,k)          &
              &   ) / xyr_CpVPTempBZ(i,j,k)  
            
          end do
        end do
      end do
      
      ! $BHoHyJ,4X?t(B
      !   $BG[Ns(B F0 $B$NE:;z$NHO0O$O(B, 1:nx, 1:ny, kmin:kmax
      !   $BG[Ns(B F $B$r5a$a$k:]$K(B F0 $B$r(B z $BJ}8~$KHyJ,$9$k$?$a(B. 
      !
      do k = kmin, kmax-1
        do j = 1, ny
          do i = 1, nx
            
            F0(i,j,k)  =                            &
              & + xyr_DensVPTempBZ(i,j,k)           &
              &   * (                               &
              &     + xyr_VelZNs(i,j,k)             & 
              &     - xyr_CpVPTempBZ(i,j,k)         &
              &       * (1.0d0 - beta)              &
              &       * (                           &
              &           xyz_ExnerNs(i,j,k+1)      & !
              &         - xyz_ExnerNs(i,j,k)        & ! xyz => xyr
              &         ) / dz * dt                 &
              &     + AlphaV                        &
              &       * (                           &
              &           xyz_VelDivNs(i,j,k+1)     &
              &         - xyz_VelDivNs(i,j,k)       &
              &         ) / dz * dt                 &
              &     + xyr_DVelZDtNl(i,j,k) * dt     &
              &     )
          end do
        end do
      end do
      
      !$B9TNs7W;;$N$?$a$N78?t(B
      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            F(i,j,k) = &
              & - beta * xyz_F1BZ(i,j,k) * dt   &
              &   * (                           &
              &       F0(i,j,k)                 &
              &     - F0(i,j,k-1)               &
              &     ) / dz                      &
              & + xyz_DExnerDtNl(i,j,k) * dt 
            
          end do
        end do
      end do
    
      !$B9TNs7W;;$N$?$a$N78?t(B
      !  $BE:;z$NHO0O$O(B, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz 
        do j = 1, ny
          do i = 1, nx
            
            D1(i,j,k) =                                 &
              & + xyz_ExnerNs(i,j,k)                    &
              & - (1.0d0 - beta)                        &
              &   * xyz_F1BZ(i,j,k) * dt                &
              &   * (                                   &
              &       xyr_DensVPTempBZ(i,j,k)           &
              &       * xyr_VelZNs(i,j,k)               &
              &     - xyr_DensVPTempBZ(i,j,k-1)         &
              &       * xyr_VelZNs(i,j,k-1)             &
              &     ) / dz                              &
              & - (xyz_VelSW(i,j,k) ** 2.0d0) * dt      &
              &   / (CpDry * xyz_VPTempBZ(i,j,k))       &
              &   * (                                   &
              &     + (                                 &
              &         pyz_VelXAs(i,j,k)               &
              &       - pyz_VelXAs(i-1,j,k)             &
              &       ) / dx                            &
              &     )                                   &
              & + F(i,j,k)
          
          end do
        end do
      end do
      
      ! $B9TNs7W;;$N$?$a$N78?t(B
      !
      do j = 1, ny
        do i = 1, nx

          D1(i,j,1) =                                    &
            & + D1(i,j,1)                                &
            & - beta * xyz_F1BZ(i,j,1) * (dt ** 2.0d0)   &
            &   * xyr_CpDensVPTemp2BZ(i,j,0)             &
            &   * E(i,j,0)                               &
            &   / dz
          
          D1(i,j,nz) =                                   &
            & + D1(i,j,nz)                               &
            & + beta * xyz_F1BZ(i,j,nz) * (dt ** 2.0d0)  &
            &   * xyr_CpDensVPTemp2BZ(i,j,nz)            &
            &   * E(i,j,nz)                              &
            &   / dz
        end do
      end do
      
      !-----------------------------------------------------------
      !$BO"N)0l<!J}Dx<0$N2r$r5a$a$k(B
      
      !$BJQ?t$N=i4|2=(B
      !
      NRHS = M
      INFO = 0
      LDB  = N
      
      ! LAPACK $B$N;EMM$K9g$o$;$FJQ7A(B 
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            D(i + nx * (j - 1), k) =  D1(i,j,k)
          end do
        end do
      end do
      
      TX = transpose( D )
      
      !$B2r9TNs$N7W;;(B. LAPACK $B$r;HMQ(B. 
      !
      call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)
      
      !$B2r$N%3%s%G%#%7%g%s$r%A%'%C%/(B. 
      !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

      !$BLa$jCM$r=PNO(B
      !
      X = transpose( TX )
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyz_ExnerAs(i,j,k) = X(i + nx * (j - 1 ), k)
          end do
        end do
      end do
      xyz_ExnerAs(imin:0,:,:) = 0.0d0
      xyz_ExnerAs(:,:,kmin:0) = 0.0d0
      xyz_ExnerAs(nx+1:imax,:,:) = 0.0d0
      xyz_ExnerAs(:,:,nz+1:kmax) = 0.0d0
      
      ! $B$N$jBe$KCM$rF~$l$k(B
      !
      call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
      
    end subroutine Exner_HEVI
    
  end subroutine Dynamics2D_Short_forcing

  
!!!--------------------------------------------------------------------!!!
  subroutine Dynamics_VI_init()
    !
    ! $BNO3X%3%"(B $B1"2rK!ItJ,$N=i4|2=%k!<%A%s(B
    ! * $B%(%/%9%J!<4X?t$r1"2rK!$G2r$/:]$KI,MW$H$J$k(B, $B78?t9TNs$NMWAG$r7h$a(B, 
    !   LU $BJ,2r$r9T$&(B. 
    !

    !$B%b%8%e!<%k8F$S=P$7(B
    use dc_message, only : MessageNotify
    use gridset, only:  kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                 kmax,            &! z $BJ}8~$NG[Ns$N>e8B(B
      &                 nx,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 ny,              &! x $BJ}8~$NJ*M}NN0h$N>e8B(B
      &                 nz                ! y $BJ}8~$NJ*M}NN0h$N>e8B(B
    use axesset, only:  dz
    use constants,only: CpDry             ! $B4%Ag@.J,$NHfG.(B
    use timeset, only : DelTimeShort
    use basicset, only: xyz_VelSW,       &!$B4pK\>l$N2;B.(B 
      &                 xyz_DensBZ,      &!$B4pK\>l$NL)EY(B
      &                 xyz_VPTempBZ      !$B4pK\>l$N290L(B

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    real(DP)  :: dt      ! $BC;$$;~4V3J;R(B
    integer   :: INFO    !$B2r$N%3%s%G%#%7%g%s%A%'%C%/(B
    integer   :: i, j, k

    !----------------------------------------------------------------
    ! $B=i4|2=(B

    ! $BJQ?tL>$,D9$9$.$?$N$G(B, $BL>A0$rCV$-49$($k(B
    !
    dt = DelTimeShort

    ! $BG[Ns$N3d$jIU$1(B
    !
    allocate( A(1:nz) )
    allocate( B(2:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_F1BZ(1:nx,1:ny,1:nz) )
    allocate( xyr_CpDensVPTemp2BZ(1:nx,1:ny,kmin:kmax) )
    allocate( xyr_DensVPTempBZ(1:nx,1:ny,kmin:kmax) )
    allocate( xyr_CpVPTempBZ(1:nx,1:ny,kmin:kmax) )

    !----------------------------------------------------------------
    ! $B78?t9TNs$*$h$S6&DL$7$FMxMQ$5$l$kG[Ns$NCM$r7h$a$k(B
    !   A, B, C $B$r5a$a$k:](B, $B4pK\>l(B (BZ) $B$NNL$O(B X $BJ}8~$K0lMM$J$N$G(B. 
    !   nx, ny $B$NCM$GBeI=$5$;$k$3$H$H$7$?(B. 
    !
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyz_F1BZ(i,j,k) =                                                    &
            &  ( xyz_VelSW(i,j,k) ** 2.0d0 )                                   &
            &  / ( CpDry * xyz_DensBZ(i,j,k) * (xyz_VPTempBZ(i,j,k) ** 2.0d0) )
        end do
      end do
    end do

    do k = kmin, kmax - 1
      do j = 1, ny        
        do i = 1, nx
          xyr_CpDensVPTemp2BZ(i,j,k)=                    &
            &  CpDry                                     &
            &  * (                                       &
            &    + xyz_DensBZ(i,j,k+1)                   &
            &      * ( xyz_VPTempBZ(i,j,k+1) ** 2.0d0 )  &
            &    + xyz_DensBZ(i,j,k)                     &
            &      * ( xyz_VPTempBZ(i,j,k) ** 2.0d0 )    &
            &    ) * 5.0d-1
        end do
      end do
    end do
    xyr_CpDensVPTemp2BZ(:,:,kmax) = 0.0d0  !$B7jKd$a(B
    
    do k = kmin, kmax-1
      do j = 1, ny
        do i = 1, nx
          xyr_DensVPTempBZ(i,j,k) =                             &
            & + (                                               &
            &   + xyz_DensBZ(i,j,k+1) * xyz_VPTempBZ(i,j,k+1)   &
            &   + xyz_DensBZ(i,j,k)   * xyz_VPTempBZ(i,j,k)     &
            &   ) * 5.0d-1    
        end do
      end do
    end do
    xyr_DensVPTempBZ(:,:,kmax) = 0.0d0  !$B7jKd$a(B

    do k = kmin, kmax-1
      do j = 1, ny
        do i = 1, nx
          xyr_CpVPTempBZ(i,j,k) =          &
            &   CpDry                      &
            &   * (                        &
            &     + xyz_VPTempBZ(i,j,k+1)  &
            &     + xyz_VPTempBZ(i,j,k)    &
            &     ) * 5.0d-1
        end do
      end do
    end do
    xyr_CpVPTempBZ(:,:,kmax) = 0.0d0
          
    do k = 2, nz-1
      A(k) =                                        &
        & + 1.0d0                                   &
        & + ( beta ** 2.0d0 )                       &
        &    * xyz_F1BZ(nx,ny,k) * ( dt * dt )      &
        &    * (                                    &
        &         xyr_CpDensVPTemp2BZ(nx,ny,k)      &
        &       + xyr_CpDensVPTemp2BZ(nx,ny,k-1)    &
        &       )                                   &
        &    / ( dz * dz )
    end do

    A(1) =                                   &
      & + 1.0d0                              &
      & + ( beta ** 2.0d0 )                  &
      &   * xyz_F1BZ(nx,ny,1) * ( dt * dt )  &
      &   * xyr_CpDensVPTemp2BZ(nx,ny,1)     &
      &   / ( dz * dz ) 

    A(nz) =                                  &
      & + 1.0d0                              &
      & + ( beta ** 2.0d0 )                  &
      &   * xyz_F1BZ(nx,ny,nz) * ( dt * dt ) &
      &   * xyr_CpDensVPTemp2BZ(nx,ny,nz-1)  &
      &   / ( dz * dz )  

    do k = 2, nz
      B(k) =                                     &
        & - ( beta ** 2.0d0 )                    &
        &   * xyz_F1BZ(nx,ny,k-1) * ( dt * dt )  &
        &   * xyr_CpDensVPTemp2BZ(nx,ny,k-1)     &
        &   / ( dz * dz )
    end do
    
    do k = 1, nz-1
      C(k) =                                     &
        & - ( beta ** 2.0d0 )                    &
        &   * xyz_F1BZ(nx,ny,k+1) * (dt * dt )   &
        &   * xyr_CpDensVPTemp2BZ(nx,ny,k)       &
        &   / ( dz * dz )
    end do

    !----------------------------------------------------------------
    ! $B78?t9TNs$r(B LU $BJ,2r(B
    !
    !$BG[Ns$NBg$-$5$rJ]4I(B
    N    = nz                   !$B78?t9TNs(B/$B2~9TNs$N<!?t(B, $B@09g@#K!(B
    M    = nx * ny              !$BJ}Dx<0$NAH?t(B
    NUD  = 1                    !$B78?t9TNs$N>e;03QItJ,$NBSI}(B
    NLD  = 1                    !$B78?t9TNs$N2<;03QItJ,$NBSI}(B
    NAL  = NLD                  !LU $BJ,2r$N7k2L(B L $B$N@09g@#K!(B
    NA   = NUD + NLD + 1
    INFO = 0

    ! $BG[Ns$N3d$jEv$F(B
    !
    allocate( AL1(N), IP(N) )

    ! $B2r9TNs$N7W;;(B. LAPACK $B$r;HMQ(B. 
    !
    call DGTTRF(N, C, A, B, AL1, IP, INFO)
    
    ! $B2r$N%3%s%G%#%7%g%s$r%A%'%C%/(B. 
    !
    if (INFO /= 0) then
      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
      stop
    end if
    
  end subroutine Dynamics_VI_init
  


  subroutine Dynamics_Tendency_Output
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
    integer :: f

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

    do f = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtAdv', &
        & dims=(/'x','y','z','t'/),     &
        & longname='Advection term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
      
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtDiff', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Diffusion term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtFall', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Fall term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
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

  end subroutine Dynamics_Tendency_Output

  
end module DynamicsHEVI_v2

