!= Module HeatFlux
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: surfaceflux_diff.f90,v 1.16 2014/11/07 06:46:44 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Surfaceflux_diff
  !
  !$B2<It6-3&$G$N%U%i%C%/%9$N7W;;%b%8%e!<%k(B. $B2<It6-3&$N3H;678?t$r7h$a$&$A$9$k(B. 
  !
  
  !$B%b%8%e!<%kFI$_9~$_(B
  use dc_types, only: DP, STRING

  !$B0EL[$N7?@k8@6X;_(B
  implicit none

  !$BB0@-$N;XDj(B
  private

  !$B4X?t$r(B public $B$K@_Dj(B
  public surfaceflux_diff_init
  public surfaceflux_diff_forcing

  !$BJQ?tDj5A(B
  real(DP), save  :: Kappa = 800.0d0               ! $B2<It6-3&$G$NMpN.3H;678?t(B
  real(DP), save  :: FactorDExnerDtSurf = 1.0d0    ! Flag for diabatice heating term in pressure equation 
  character(STRING), parameter:: module_name = 'surfaceflux_diff'

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Diff_init
    !
    !NAMELIST $B$+$iI,MW$J>pJs$rFI$_<h$j(B, $B;~4V4XO"$NJQ?t$N@_Dj$r9T$&(B. 
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use namelist_util,     only : namelist_filename
    use composition,       only : SpcWetSymbol
    use gridset,           only : ncmax

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BFbItJQ?t(B
    integer    :: l, unit

    !---------------------------------------------------------------    
    ! NAMELIST $B$+$i>pJs$r<hF@(B
    !
    NAMELIST /surfaceflux_diff_nml/ &
      & Kappa,                      &
      & FactorDExnerDtSurf 

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_diff_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "Kappa = %f", d=(/Kappa/))
    call MessageNotify( 'M', module_name, "FactorDExnerDtSurf = %f", &
      &                  d=(/FactorDExnerDtSurf/) )

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtFlux',         &
      & dims=(/'x','y','z','t'/), &
      & longname='surface flux of potential temperature', &
      & units='kg.kg-1.s-1',            &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='ExnerFlux',      &
      & dims=(/'x','y','z','t'/), &
      & longname='surface flux of Exner function', &
      & units='kg.kg-1.s-1',            &
      & xtype='double')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtFlux', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Surface Flux term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
    end do

  end subroutine Surfaceflux_Diff_init


!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Diff_forcing( &
    &   xyz_PTemp, xyzf_QMix,          &
    &   xyz_DPTempDt, xyz_DExnerDt,    &
    &   xyzf_DQMixDt                   &
    & )
    ! 
    ! $B2<It6-3&$+$i$N%U%i%C%/%9$K$h$k29EY$NJQ2=N($r(B,
    ! $B%P%k%/J}K!$K4p$E$$$F7W;;$9$k(B.
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut    
    use gridset,           only : imin,         & !x $BJ}8~$NG[Ns$N2<8B(B
      &                           imax,         & !x $BJ}8~$NG[Ns$N>e8B(B
      &                           jmin,         & !y $BJ}8~$NG[Ns$N2<8B(B
      &                           jmax,         & !y $BJ}8~$NG[Ns$N>e8B(B
      &                           kmin,         & !z $BJ}8~$NG[Ns$N2<8B(B
      &                           kmax,         & !z $BJ}8~$NG[Ns$N>e8B(B
      &                           nx, ny, nz, ncmax
    use axesset,           only : z_dz            !z $BJ}8~$N3J;RE@4V3V(B
    use timeset,           only : TimeN
    use composition,       only : SpcWetSymbol, GasNum
    use DExnerDt,          only : xyz_DExnerDt_xyz_xyzf

    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none
    
    !$BJQ?tDj5A(B
    !
    real(DP), intent(in)   :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                                           !$B290L$N>qMp@.J,(B    
    real(DP), intent(in)   :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                           !$B290L$N>qMp@.J,(B    
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP)               :: xyz_DPTempDt0(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyz_DExnerDt0(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyzf_DQMixDt0(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP)               :: xyz_Heatflux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyz_Exnerflux(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)               :: xyzf_QMixflux(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                            !$BCOI=LLG.%U%i%C%/%9(B
    integer, parameter     :: kz = 1        !$BG[NsE:;z(B
    integer                :: l             !$B%k!<%WJQ?t(B

    ! $B=i4|2=(B
    !
    xyz_HeatFlux  = 0.0d0
    xyz_ExnerFlux = 0.0d0
    xyzf_QMixFlux = 0.0d0

    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt
    xyzf_DQMixDt0 = xyzf_DQMixDt

    ! $BCOI=LLG.%U%i%C%/%9$K$h$k2CG.N($r7W;;(B
    !   * $BC10L$O(B K/s
    !   * $B%(%/%9%J!<4X?t$O4pK\>l$NCM$GBeI=$5$;$k(B.     
    !   * $B3J;RE@(B xz $B$G$O(B, $BJ*M}NN0h$N:G2<C<$NE:$(;z$O(B kz = 1

    xyz_HeatFlux(:,:,kz) =                                  &
!     &  - Kappa * xyz_PTemp(:,:,kz) * xyz_ExnerBZ(:,:,kz)  &  !check
      &  - Kappa * xyz_PTemp(:,:,kz)                        &
      &    / ( ( z_dz(kz) * 5.0d-1 ) ** 2.0d0 )  

    do l = 1, GasNum
      xyzf_QMixFlux(:,:,kz,l) =                        &
        &     - Kappa * xyzf_QMix(:,:,kz,l)            &
        &        / ( ( z_dz(kz) * 5.0d-1 ) ** 2.0d0 )  
    end do
    
    xyz_DPTempDt = xyz_DPTempDt0 + xyz_Heatflux
    xyzf_DQMixDt = xyzf_DQMixDt0 + xyzf_QMixflux
    
    xyz_ExnerFlux = xyz_DExnerDt_xyz_xyzf( xyz_HeatFlux, xyzf_QMixFlux ) * FactorDExnerDtSurf
    xyz_DExnerDt = xyz_DExnerDt0 + xyz_ExnerFlux
    
    call HistoryAutoPut(TimeN, 'DPTempDtFlux', xyz_HeatFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerFlux', xyz_ExnerFlux(1:nx,1:ny,1:nz))
    do l = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(l))//'DtFlux', xyzf_Qmixflux(1:nx,1:ny,1:nz,l))
    end do    

  end subroutine Surfaceflux_Diff_forcing
  
end module Surfaceflux_diff
