!= Module cloudphys_k1969
!
! Authors::   $B?y;39L0lO/(B(SUGIYAMA Ko-ichiro), $B>.9b@5;L(B (ODAKA Masatsugu), $B9b66K'9,(B (YOSHIYUKI Takahashi)
! Version::   $Id: cloudphys_k1969.f90,v 1.28 2014/03/04 05:55:05 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module cloudphys_k1969
  !
  !$BCH$+$$1+$N%P%k%/K!$rMQ$$$?(B, $B?e>x5$$H1+(B, $B1@$H1+$N:.9gHf$NJQ4978?t$r5a$a$k(B.
  !   * $BCfEg7r2p(B (1994) $B$GMxMQ$7$?Dj<0$r$=$N$^$^MxMQ(B. 
  ! 
  
  !$B%b%8%e!<%kFI$_9~$_(B
  use dc_types,   only : DP
  
  !$B0EL[$N7?@k8@6X;_(B
  implicit none
  
  !$BB0@-$N;XDj(B
  private

  !$B4X?t$r(B public $B$K$9$k(B
  public Cloudphys_K1969_Init
  public Cloudphys_K1969_forcing

  real(DP), save :: FactorCloud2Rain = 1.0d0 !$B1@$+$i1+$X$NJQ49$NM-L5(B 
                                             !$B1+$XJQ49$5$;$J$$>l9g$OCM$r%<%m$K$9$k(B. 
  real(DP), save :: FactorRain2Gas = 1.0d0   !$B1+$+$i>x5$$X$NJQ49$NM-L5(B 
                                             !$B>x5$$XJQ49$5$;$J$$>l9g$OCM$r%<%m$K$9$k(B. 
  real(DP), save, public :: FactorCloud2Gas = 1.0d0 
                                             !$B1@$+$i>x5$$X$NJQ49$NM-L5(B 
                                             !$B>x5$$XJQ49$5$;$J$$>l9g$OCM$r%<%m$K$9$k(B. 

!  real(DP), save :: FactorJ      = 1.0d0 !$B1@J*M}2aDx$N%Q%i%a!<%?(B
!                                         !$BLZ@1$G$O(B 3.0d0
!                                         !$BCO5e$G$O(B 1.0d0 $B$H$9$k(B
  real(DP), save :: AutoConvTime = 1.0d3 !$BJ;9g@.D9$N;~Dj?t(B [sec]
  real(DP), save :: QMixCr       = 1.0d-3 
                                         !$BJ;9g@.D9$r@8$8$kNW3&:.9gHf(B [kg/kg]
  real(DP), save :: FactorDExnerDtCloud = 1.0d0

contains  

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_Init
    !
    ! $B=i4|2=%k!<%A%s(B
    !

    !$B%b%8%e!<%k8F$S=P$7(B
    use dc_types,      only : DP, STRING
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use namelist_util, only : namelist_filename

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BFbItJQ?t(B
    integer  :: unit    !$BAuCVHV9f(B
    character(*), parameter:: module_name = 'Cloudphys_K1969_Init'

    !-----------------------------------------------------------
    ! NAMELIST $B$+$i>pJs$r<hF@(B
    !
    NAMELIST /cloudphys_k1969_nml/                &
      & AutoConvTime, QMixCr,                     &
      & FactorDExnerDtCloud,                      &
      & FactorCloud2Rain, FactorRain2Gas, FactorCloud2Gas

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=cloudphys_k1969_nml)
    close(unit)

    !-----------------------------------------------------------
    ! $B=PNO(B
    !
    call MessageNotify( "M", &
      &  module_name, "AutoConvTime = %f",  d=(/AutoConvTime/) )
    call MessageNotify( "M", &
      &  module_name, "QMixCr = %f",  d=(/QMixCr/) )
    call MessageNotify( "M", &
      &  module_name, "FactorCloud2Rain = %f",  d=(/FactorCloud2Rain/) )
    call MessageNotify( "M", &
      &  module_name, "FactorRain2Gas = %f",  d=(/FactorRain2Gas/) )
    call MessageNotify( "M", &
      &  module_name, "FactorCloud2Gas = %f",  d=(/FactorCloud2Gas/) )
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtCloud= %f", d=(/ FactorDExnerDtCloud /))

     !
     ! HistoryAuto
     !
     call Cloudphys_K1969_HistoryAuto

    
   end subroutine Cloudphys_K1969_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_forcing(         &
    & xyz_ExnerNl,                            &!(in)
    & xyz_DExnerDt, xyz_PTempAl, xyzf_QMixAl  &!(inout)
    & )

    !$B%b%8%e!<%k8F$S=P$7(B
    use dc_types,only : DP, STRING
    use gtool_historyauto,                 &
      &          only : HistoryAutoPut
    use timeset, only : DelTimeLong, TimeN
    use gridset, only : imin,              &!x $BJ}8~$NG[Ns$N2<8B(B
      &                 imax,              &!x $BJ}8~$NG[Ns$N>e8B(B
      &                 jmin,              &!y $BJ}8~$NG[Ns$N>e8B(B
      &                 jmax,              &!y $BJ}8~$NG[Ns$N>e8B(B
      &                 kmin,              &!z $BJ}8~$NG[Ns$N2<8B(B
      &                 kmax,              &!z $BJ}8~$NG[Ns$N>e8B(B
      &                 nx, ny, nz, ncmax   !$BJ*M}NN0h$NBg$-$5(B
    use constants,only: FactorJ,           &!
      &                 PressBasis,        &!$B290L$N4p=`05NO(B 
      &                 CpDry,             &!$B4%Ag@.J,$NHfG.(B
      &                 MolWtDry,          &!
      &                 GasRDry             !$B4%Ag@.J,$N5$BNDj?t(B 
    use basicset, only: xyz_DensBZ,        &!$B4pK\>l$NL)EY(B
      &                 xyz_PTempBZ,       &!$B4pK\>l$N290L(B
      &                 xyz_ExnerBZ,       &!$B4pK\>l$NL5<!8505NO(B
      &                 xyzf_QMixBZ         !$B4pK\>l$N:.9gHf(B
    use composition,                       &
      &           only:  MolWtWet,         &!
      &                 SpcWetID,          &!
      &                 SpcWetSymbol,      &!
      &                 CondNum,           &!$B6E7k2aDx$N?t(B
      &                 IdxCG,             &!$B6E7k2aDx(B($B>x5$(B)$B$NG[NsE:$(;z(B
      &                 IdxCC,             &!$B6E7k2aDx(B($B1@(B)$B$NG[NsE:$(;z(B
      &                 IdxCR,             &!$B6E7k2aDx(B($B1@(B)$B$NG[NsE:$(;z(B
      &                 GasNum,            &!$B>x5$$N?t(B
      &                 CloudNum,          &!$B1@$N?t(B
      &                 RainNum,           &!$B1+$N?t(B
      &                 IdxG,              &!$B>x5$$NG[NsE:$(;z(B
      &                 IdxC,              &!$B1@$NG[NsE:$(;z(B
      &                 IdxR,              &!$B1+$NG[NsE:$(;z(B
      &                 IdxNH3,            &!NH3($B>x5$(B)$B$NG[NsE:$(;z(B
      &                 IdxH2S,            &!H2S($B>x5$(B)$B$NG[NsE:$(;z(B
      &                 IdxNH4SHr           !NH4SH($B1+(B)$B$NG[NsE:$(;z(B
    use average, only : xyr_xyz
    use differentiate_center2,             &
      &          only : xyz_dz_xyr
    use ChemCalc,only : xyz_SvapPress, xyz_LatentHeat, ReactHeatNH4SH, xyz_DelQMixNH4SH    
    use MoistAdjust,                       &
      &          only : MoistAdjustSvapPress, MoistAdjustNH4SH
    use DExnerDt,only : xyz_DExnerDt_xyzf, xyz_DExnerDt_xyz
    use SetMargin,only: SetMargin_xyz, SetMargin_xyzf

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_PTempAl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyzf_QMixAl(imin:imax, jmin:jmax, kmin:kmax, ncmax)

    real(DP) :: xyz_PTempOrig(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_PTempWork(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_DelPTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_PTempCond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyzf_QMixOrig(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_QMixWork(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_DelQMix(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_QMixCond(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: DelTime
    integer  :: s
    integer  :: iG, iC, iR

    real(DP) :: xyzf_Cloud2Rain(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                          !$B1@$+$i1+$X$NJQ49NL(B
    real(DP) :: xyz_AutoConv(imin:imax,jmin:jmax,kmin:kmax)
                                          !$BK0OB:.9gHf(B
    real(DP) :: xyz_Collect(imin:imax,jmin:jmax,kmin:kmax)
                                          !$B5,3J2=$5$l$?@xG.(B
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                          !$B:.9gHf$N>qMp@.J,(B + $BJ?6Q@.J,(B
    real(DP) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                          !$B29EY$N>qMp@.J,(B + $BJ?6Q@.J,(B
    real(DP) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
                                          !$BA405(B
    real(DP) :: xyz_ExnerAll(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_NonSaturate(imin:imax,jmin:jmax,kmin:kmax)
                                          !$BL$K0OBEY(B($BK0OB:.9gHf$H>x5$$N:.9gHf$N:9(B)
    real(DP) :: xyzf_Rain2Gas(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyzf_Rain2GasNH4SH(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyzf_DelPTemp(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyz_DelPTempNH4SH(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtCondTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtCondQMix(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_QMixSat(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP) :: xyz_QMixHum(imin:imax,jmin:jmax,kmin:kmax) 


    !-------------------------------------------------------------
    ! $B=i4|CM$rJ]4I(B Store Initial Value
    !
    xyz_PTempOrig = xyz_PTempAl
    xyzf_QMixOrig = xyzf_QMixAl

    !-------------------------------------------------------------
    ! $B;~4V9o$_I}(B. Leap-frog $B$J$N$G(B, 2 \del t
    !
    DelTime = 2.0d0 * DelTimeLong

    !------------------------------------------
    ! $BA4%(%/%9%J!<4X?t!&A405$r7W;;(B. $B%5%V%k!<%A%sFb$G$OJQ2=$7$J$$(B.
    !
    xyz_ExnerAll = xyz_ExnerNl + xyz_ExnerBZ
    xyz_PressAll = PressBasis * (xyz_ExnerAll ** (CpDry / GasRDry))

    !------------------------------------------    
    ! $BCH$+$$1+$N%Q%i%a%?%j%<!<%7%g%s(B.
    ! * $B1@(B<-->$B1+(B $B$NJQ49$r9T$&(B.
    !
    ! Warm rain parameterization.
    ! * Conversion from cloud to rain.
    
    !$B$3$l$^$G$NCM$r:n6HG[Ns$KJ]4I(B
    ! Previous values are stored to work area.
    !
    xyzf_QMixWork = xyzf_QMixAl
    
    !$B1+$X$NJQ2=NL$r7W;;(B
    ! Conversion values are calculated.
    !    
    xyzf_QMixAll = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyzf_Cloud2Rain = 0.0d0

    do s = 1, CloudNum

      ! $BCM$rJ]4I(B
      !
      iC = IdxC(s)
      iR = IdxR(s)

      !$BJ;9g@.D9(B
      !
      xyz_AutoConv =                                           &
        & DelTime / AutoConvTime                               &
        & * max( 0.0d0, ( xyzf_QMixAll(:,:,:,iC) - QMixCr) )

      !$B>WFM9gBN@.D9(B
      !
      xyz_Collect =                                            &
        &  DelTime                                             &
        &  * 2.2d0 * FactorJ * xyzf_QMixAll(:,:,:,iC)          &
        &  * (xyzf_QMixAll(:,:,:,iR) * xyz_DensBZ) ** 0.875d0  

      !$B1@$NJQ49NL(B: $BJ;9g@.D9$H9gBN>WFM$NOB(B
      !  $B85!9$NJQ2=NL$r>e8BCM$H$7$F@_Dj$9$k(B. $BIi$NCM$H$J$k(B.
      !
      xyzf_Cloud2Rain(:,:,:,iC) =                                        &
        & - min( xyzf_QMixAll(:,:,:,iC), ( xyz_AutoConv + xyz_Collect ) )
      
      !$B1+$NJQ49NL(B. $BId9f$O1@$NJQ49NL$H$OH?BP(B. 
      xyzf_Cloud2Rain(:,:,:,iR) = - xyzf_Cloud2Rain(:,:,:,iC) 

      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iC))//'DtC2R',   &
        & xyzf_Cloud2Rain(1:nx,1:ny,1:nz,iC) / DelTime)

      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtC2R',   &
        & xyzf_Cloud2Rain(1:nx,1:ny,1:nz,iR) / DelTime)

    end do

    ! $BJQ2=NL$rB-$79~$`(B
    ! $B1@$+$i1+$XJQ49$5$;$J$$>l9g$O(B FactorCloud2Rain = 0.0 $B$H$9$k(B. 
    !
    xyzf_QMixAl = xyzf_QMixWork + xyzf_Cloud2Rain * FactorCloud2Rain


    !-------------------------------------------    
    ! $BCH$+$$1+$N%Q%i%a%?%j%<!<%7%g%s(B.
    ! * $B>x5$(B<-->$B1+(B $B$NJQ49$r9T$&(B
    !
    ! Warm rain parameterization.
    ! * Conversion from rain to vapor.
    
    !$B$3$l$^$G$NCM$r:n6HG[Ns$KJ]4I(B
    ! Previous values are stored to work area.
    !
    xyz_PTempWork = xyz_PTempAl
    xyzf_QMixWork = xyzf_QMixAl
    
    ! $B1+$+$i>x5$$X$N:.9gHfJQ2=$r5a$a$k(B
    ! * $B290L$N7W;;$K$*$$$F(B, $B:.9gHfJQ2=$,I,MW$H$J$k$?$a(B, 
    !   $B:.9gHfJQ2=$r(B 1 $B$D$NG[Ns$H$7$FMQ0U$9$k(B.
    !
    ! Conversion values are calculated.
    !

    !$B29EY(B, $B05NO(B, $B:.9gHf$NA4NL$r5a$a$k(B
    !$B>qMp@.J,$HJ?6Q@.J,$NB-$7;;(B
    !
    xyz_TempAll   = ( xyz_PTempAl + xyz_PTempBZ ) * ( xyz_ExnerNl + xyz_ExnerBZ )
    xyzf_QMixAll  = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry  = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )

    xyzf_Rain2Gas = 0.0d0
    xyzf_DelPTemp = 0.0d0
    xyzf_Rain2GasNH4SH = 0.0d0
    xyz_DelPTempNH4SH  = 0.0d0

    do s = 1, CondNum

       ! $BCM$rJ]4I(B
       !
       iG = IdxCG(s)
       iC = IdxCC(s)
       iR = IdxCR(s)

      !$BK0OB>x5$05$H:.9gHf$N:9(B($BK0OBEY(B)$B$r7W;;(B. 
      !  $B1+$+$i>x5$$X$NJQ49NL$OK0OBEY$KHfNc$9$k(B.
      !
      xyz_NonSaturate =                                   &
        & max(                                            &
        &   0.0d0,                                        &
        &   xyz_SvapPress(SpcWetID(iC), xyz_TempAll)      &
        &     * MolWtWet(iG) / ( MolWtDry * xyz_PressDry) &
        &     - xyzf_QMixAll(:,:,:,iG)                    &
        &    )

      !$B1+$NJQ49NL(B
      !  $B85!9$N1+N3$N:.9gHf0J>e$K>xH/$,@8$8$J$$$h$&$K>e8BCM$r@_Dj(B
      !
      xyzf_Rain2Gas(:,:,:,iR) =                                    &
        & - min(                                                   &
        &    DelTime * 4.85d-2 * FactorJ * xyz_NonSaturate         &
        &     * ( xyzf_QMixAll(:,:,:,iR) * xyz_DensBZ )** 0.65d0,  &
        &    xyzf_QMixAll(:,:,:,iR)                                &
        &   ) 

      !$B>x5$$NJQ49NL(B
      !  $B1+N3$NJQ49NL$H$OId9f$,5U$H$J$k(B
      !
      xyzf_Rain2Gas(:,:,:,iG) = - xyzf_Rain2Gas(:,:,:,iR) 
    
      ! xyzf_DelQMix $B$r85$K@xG.$r7W;;(B
      !
      xyzf_DelPTemp(:,:,:,s) =                          &
        & xyz_LatentHeat( SpcWetID(iR), xyz_TempAll )   &
        &  * xyzf_Rain2Gas(:,:,:,iR)                    &
        &  / (xyz_ExnerAll * CpDry) 

    end do

    !$BK0OB>x5$05$H:.9gHf$N:9(B($BK0OBEY(B)$B$r7W;;(B. 
    !  $B1+$+$i>x5$$X$NJQ49NL$OK0OBEY$KHfNc$9$k(B.
    !  $BL$K0OBEY$r5a$a$?$$$N$G(B, $B%^%$%J%9$r$+$1;;$7$F$$$k(B
    !  (DelQMixNH4SH $B$O(B, NH4SH $B$,A}2C$9$kJ}8~(B, $B$9$J$o$AK0OBEY$r@5$H$7$F$$$k(B)
    !
    if (IdxNH4SHr /= 0) then 
      xyz_NonSaturate =                                                 &
        & max(                                                          &
        &  0.0d0,                                                       &
        &   - xyz_DelQMixNH4SH(                                         &  
        &       xyz_TempAll, xyz_PressAll, xyz_PressDry,                & 
!        &       xyz_TempAll, xyz_PressAll,                              & 
        &       xyzf_QMixAll(:,:,:,IdxNH3), xyzf_QMixAll(:,:,:,IdxH2S), &
        &       MolWtWet(IdxNH3), MolWtWet(IdxH2S)                      &
        &     )                                                         &
        &  )

      !$B1+$NJQ49NL(B
      !  $B85!9$N1+N3$N:.9gHf0J>e$K>xH/$,@8$8$J$$$h$&$K>e8BCM$r@_Dj(B
      !
      xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) =                              &
        & - min(                                                         &
        &     DelTime * 4.85d-2 * FactorJ * xyz_NonSaturate              &
        &      * (xyzf_QMixAll(:,:,:,IdxNH4SHr) * xyz_DensBZ) ** 0.65d0, &
        &     xyzf_QMixAll(:,:,:,IdxNH4SHr)                              &
        &    ) 
     
      !$B>x5$$NJQ49NL(B
      !  $B1+N3$NJQ49NL$H$OId9f$,5U$H$J$k(B
      !
      xyzf_Rain2GasNH4SH(:,:,:,IdxNH3) =                           &
        & - xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) * MolWtWet(IdxNH3) &
        &   / MolWtWet(IdxNH4SHr)
      xyzf_Rain2GasNH4SH(:,:,:,IdxH2S) =                           &
        & - xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) * MolWtWet(IdxH2S) &
        &   / MolWtWet(IdxNH4SHr)

      xyz_DelPTempNH4SH                                          &
        & = ReactHeatNH4SH * xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) &
        &    / (xyz_ExnerAll * CpDry)

    end if

    !$BJQ2=NL$rB-$7;;(B
    !$B1+$+$i>x5$$X$NJQ49$r@Z$k>l9g$O(B FactorRain2Gas = 0.0 $B$H$9$k(B. 
    !
    xyzf_DelQMix = xyzf_Rain2Gas + xyzf_Rain2GasNH4SH 
    xyz_DelPTemp = sum(xyzf_DelPTemp, 4) + xyz_DelPTempNH4SH 

    ! $B290L$H:.9gHf$N7W;;(B. $B1+$+$i>x5$$X$NJQ49J,$rDI2C(B
    !
    xyz_PTempAl = xyz_PTempWork + xyz_DelPTemp * FactorRain2Gas
    xyzf_QMixAl = xyzf_QMixWork + xyzf_DelQMix * FactorRain2Gas

    call HistoryAutoPut(TimeN, 'PTempR2G',  &
      &  xyz_DelPTemp(1:nx,1:ny,1:nz) / DelTime)
    do s = 1, GasNum
      iG = IdxG(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtR2G',  &
        & xyzf_DelQMix(1:nx,1:ny,1:nz, iG) / DelTime)
    end do
    do s = 1, RainNum
      iR = IdxR(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtR2G',  &
        & xyzf_DelQMix(1:nx,1:ny,1:nz, iR) / DelTime)
    end do

    !-------------------------------------------
    ! $B<>=aK0OBD4@a(B
    ! * $B>x5$(B<-->$B1@$NJQ49$r9T$&(B.
    !
    ! Moist adjustment.
    ! * Conversion from vapor to cloud.
    !

    !$B$3$l$^$G$NCM$r:n6HG[Ns$KJ]4I(B
    ! Previous values are stored to work area.
    !
    xyz_PTempWork = xyz_PTempAl
    xyzf_QMixWork = xyzf_QMixAl
    
    ! $B4%Ag@.J,$N05NO$H:.9gHf$NA4NL$r5a$a$k(B
    !
    xyzf_QMixAll = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )
    
    ! $B<>=aK0OBD4@a(B
    !
    call MoistAdjustSvapPress(   &
      & xyz_PressDry,            & ! (in) 
      & xyz_ExnerNl,             & ! (in)
      & xyz_PTempAl,             & ! (inout)
      & xyzf_QMixAl,             & ! (inout)
      & FactorCloud2Gas          & ! (in)
      & )
    if (IdxNH4SHr /= 0) then 
      call MoistAdjustNH4SH(     &
        & xyz_PressDry,          & !(in)
        & xyz_ExnerNl,           & !(in)
        & xyz_PTempAl,           & !(inout)
        & xyzf_QMixAl,           & !(inout)
        & FactorCloud2Gas        & !(in)
        & )
    end if

    ! Output
    !
    xyz_PTempCond = (xyz_PTempAl - xyz_PTempWork) / DelTime
    xyzf_QMixCond = (xyzf_QMixAl - xyzf_QMixWork) / DelTime

    call HistoryAutoPut(TimeN, 'PTempG2C', xyz_PTempCond(1:nx,1:ny,1:nz))
    do s = 1, GasNum
      iG = idxG(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtG2C',  &
        & xyzf_QMixCond(1:nx,1:ny,1:nz,iG) )
    end do
    do s = 1, CloudNum
      iC = idxC(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iC))//'DtG2C',  &
        & xyzf_QMixCond(1:nx,1:ny,1:nz,iC) )
    end do


    !-----------------------------------
    ! $B%(%/%9%J!<4X?t$N(B tendency
    !

    ! $B=i4|CM$+$i$N:9$r<h$k(B
    !
    xyz_PTempCond = (xyz_PTempAl - xyz_PTempOrig) / DelTime
    xyzf_QMixCond = (xyzf_QMixAl - xyzf_QMixOrig) / DelTime

    xyz_DExnerDtCondTemp = xyz_DExnerDt_xyz( xyz_PTempCond )  * FactorDExnerDtCloud
    xyz_DExnerDtCondQMix = xyz_DExnerDt_xyzf( xyzf_QMixCond ) * FactorDExnerDtCloud

    xyz_DExnerDt  = xyz_DExnerDt + xyz_DExnerDtCondTemp + xyz_DExnerDtCondQMix

    call HistoryAutoPut(TimeN, 'DExnerDtCondTemp', xyz_DExnerDtCondTemp(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtCondQMix', xyz_DExnerDtCondQMix(1:nx,1:ny,1:nz))


    !----------------------------------------------------------------
    ! $BK0OB>x5$05$HJ?9UDj?t(B
    !
    xyz_TempAll   = ( xyz_PTempAl + xyz_PTempBZ ) * ( xyz_ExnerNl + xyz_ExnerBZ )
    xyzf_QMixAll  = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry  = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )

    do s = 1, CondNum
      iG = IdxCG(s)
      iC = IdxCC(s)

      !$BK0OB>x5$05(B
      xyz_QMixSat =                                  &
        & xyz_SvapPress(SpcWetID(iC), xyz_TempAll)    &
        &  * MolWtWet(iC) / MolWtDry / xyz_PressDry

      !$BAjBP<>EY(B
      xyz_QMixHum = &
        & xyzf_QMixAll(:,:,:,iG) / xyz_QMixSat * 100.0d0

      !$B=PNO(B
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtHum', xyz_QMixHum(1:nx, 1:ny, 1:nz))
    end do

    ! Set Margin
    !
    call SetMargin_xyz(xyz_PTempAl)
    call SetMargin_xyzf(xyzf_QMixAl)

  end subroutine Cloudphys_K1969_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )
    !
    ! $B4%Ag@.J,$NJ,05$N7W;;(B
    ! 

    !$B%b%8%e!<%k8F$S=P$7(B
    use dc_types,    only : DP
    use gridset,     only : imin, imax, &
      &                     jmin, jmax, &
      &                     kmin, kmax, &
      &                     ncmax
    use composition, only : MolWtWet,   &!
      &                     IdxG,       &!$B>x5$$NG[NsE:$(;z(B
      &                     GasNum
    use constants,   only : MolWtDry

    !$B0EL[$N7?@k8@6X;_(B
    implicit none
    
    !$BJQ?tDj5A(B
    real(DP), intent(in) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(in) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_PressDry_xyzf_xyz(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_QMixAllPerMolWt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    integer              :: f, iG

    ! $B:.9gHf(B/$BJ,;RNL$r7W;;(B
    ! 
    do f = 1, GasNum
      iG = IdxG(f)
      xyzf_QMixAllPerMolWt(:,:,:,f) = xyzf_QMixAll(:,:,:,iG) / MolWtWet(iG)
    end do

    ! $B4%Ag@.J,$NJ,05(B. 
    ! 
    xyz_PressDry_xyzf_xyz = xyz_PressAll / (1.0d0 + MolWtDry * sum( xyzf_QMixAllPerMolWt(:,:,:,1:GasNum), 4))

!!DBG
!!    xyz_PressDry_xyzf_xyz = xyz_PressAll

  end function xyz_PressDry_xyzf_xyz


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_HistoryAuto
    !
    ! tendency $B=PNO$N$?$a$N@_Dj(B
    !

    !$B%b%8%e!<%k8F$S=P$7(B
    use gtool_historyauto, only : HistoryAutoAddVariable
    use composition,       only : SpcWetSymbol,  & 
      &                           GasNum,        &!$B>x5$$N?t(B
      &                           CloudNum,      &!$B1@$N?t(B
      &                           RainNum,       &!$B1+$N?t(B  
      &                           IdxG,          &!$B>x5$$NG[NsE:$(;z(B
      &                           IdxC,          &!$B1@$NG[NsE:$(;z(B
      &                           IdxR            !$B1+$NG[NsE:$(;z(B

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BJQ?tDj5A(B
    integer  :: l, iG, iC, iR

    call HistoryAutoAddVariable(                                     &
      & varname='PTempR2G',                                          &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Condensation term of potential temperature (R2G)', &
      & units='K.s-1',                                               &
      & xtype='double' )

    call HistoryAutoAddVariable(                                     &
      & varname='PTempG2C',                                          &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Condensation term of potential temperature (G2C)', &
      & units='K.s-1',                                               &
      & xtype='double' )

    call HistoryAutoAddVariable(                                     &
      & varname='DExnerDtCondTemp',                                     &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Latent heat term of exner function (Temp)',        &
      & units='K.s-1',                                               &
      & xtype='double')

    call HistoryAutoAddVariable(                                     &
      & varname='DExnerDtCondQMix',                                     &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Latent heat term of exner function (QMix)',        &
      & units='K.s-1',                                               &
      & xtype='double')

    do l = 1, GasNum
      iG = IdxG(l)

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtR2G',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iG))//' mixing ratio (R2G)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtG2C',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iG))//' mixing ratio (G2C)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtHum',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Humidity of ' //trim(SpcWetSymbol(iG)),          &
        & units='1',                                                 &
        & xtype='double')
      
    end do
    
    do l = 1, CloudNum
      iC = IdxC(l)
      
      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iC))//'DtC2R',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iC))//' mixing ratio (C2R)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )
      
      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iC))//'DtG2C',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iC))//' mixing ratio (G2C)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

    end do


    do l = 1, RainNum
      iR = IdxR(l)

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtR2G',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio (R2G)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

     call HistoryAutoAddVariable(                                    &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtC2R',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio (C2R)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

    end do

  end subroutine Cloudphys_K1969_HistoryAuto
  
  
end module Cloudphys_k1969
