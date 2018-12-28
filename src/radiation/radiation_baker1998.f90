!= Module Radiation
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_baker1998.f90,v 1.3 2014/03/04 05:55:05 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_baker1998
  !
  ! Baker et al (1998) $B$rLO$7$?J|<M6/@)$rM?$($k$?$a$N%b%8%e!<%k(B
  !

  ! $B%b%8%e!<%kFI$_9~$_(B
  !
  use dc_types, only: DP, STRING
  
  ! $B0EL[$N7?@k8@6X;_(B
  !
  implicit none

  ! private $BB0@-(B
  !
  private

  !$BJQ?tDj5A(B
  !
  real(DP), save, allocatable, public :: xyz_DPTempDtRad(:,:,:)  !$BJ|<M2CG.$,B8:_$9$kNN0h(B
  real(DP), save, allocatable, public :: xyz_ExnerRad(:,:,:)  !$BJ|<M2CG.$,B8:_$9$kNN0h(B

  character(STRING), parameter :: module_name = 'radiation_baker1998'
                               ! $B%b%8%e!<%k$NL>>N(B.
                               ! Module name

  public Radiation_baker1998_init
  public Radiation_baker1998_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_baker1998_init
    !
    ! $B=i4|2=%k!<%A%s(B: BAKER1998 $B$GMQ$$$i$l$F$$$kG.6/@)$rMQ$$$k(B
    !

    ! $B%b%8%e!<%k8F$S=P$7(B
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x $BJ}8~$NG[Ns$N2<8B(B
      &                           imax,     & !x $BJ}8~$NG[Ns$N>e8B(B
      &                           jmin,     & !y $BJ}8~$NG[Ns$N2<8B(B
      &                           jmax,     & !y $BJ}8~$NG[Ns$N>e8B(B
      &                           kmin,     & !z $BJ}8~$NG[Ns$N2<8B(B
      &                           kmax        !z $BJ}8~$NG[Ns$N>e8B(B
    use axesset,           only : z_Z         !Z $B:BI8<4(B($B%9%+%i!<3J;RE@(B)
    use constants,         only : CpDry       !$BHfG.(B
    use basicset,          only : xyz_ExnerBZ, & !$B%(%/%9%J!<4X?t$N4pK\>l(B
      &                           xyz_DensBZ     !$B4pK\>l$NL)EY(B
    use DExnerDt,          only : xyz_DExnerDt_xyz
    
    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    !$BFbItJQ?t(B
    !
    integer  :: k             !$B%k!<%WJQ?t(B

    allocate( xyz_DPTempDtRad(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRad(imin:imax, jmin:jmax, kmin:kmax) )
    
    ! $B290L$NJ|<M6/@)9`(B.
    ! BAKER1998 $B$GMQ$$$i$l$F$$$kG.6/@)$N%W%m%U%!%$%k$rMQ$$$k(B
    !
    do k = kmin, kmax
      xyz_DPTempDtRad(:,:,k) = cal_Qsub( z_Z(k) ) / CpDry / xyz_DensBZ(:,:,k) / xyz_ExnerBZ(:,:,k)
    end do
    xyz_ExnerRad = xyz_DExnerDt_xyz( xyz_DPTempDtRad ) 

    ! Output
    !
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='double')

  end subroutine Radiation_baker1998_init


  subroutine Radiation_baker1998_forcing(  &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x $BJ}8~$NG[Ns$N2<8B(B
      &                           imax,     & !x $BJ}8~$NG[Ns$N>e8B(B
      &                           jmin,     & !y $BJ}8~$NG[Ns$N2<8B(B
      &                           jmax,     & !y $BJ}8~$NG[Ns$N>e8B(B
      &                           kmin,     & !z $BJ}8~$NG[Ns$N2<8B(B
      &                           kmax,     & !z $BJ}8~$NG[Ns$N>e8B(B
      &                           nx, ny, nz

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BF~=PNOJQ?t(B
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency $B$N99?7(B
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRad
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRad
    
    ! $B%U%!%$%k=PNO(B
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRad(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_baker1998_forcing
  

  function cal_Qsub( z0 )
    ! 
    ! BS1989 $B$GMQ$$$i$l$F$$$?G.6/@)$N<0(B
    ! $BF1$8%k!<%A%s$,(B surfaceflux_baker1998.f90 $B$K$b4^$^$l$F$$$k$N$GCm0U$9$k$3$H(B.

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types,          only : DP

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BJQ?tDj5A(B
    !
    real(DP), intent(in) :: z0
    real(DP)             :: cal_Qsub
    real(DP), parameter  :: z_Upper = 6.70d4      ! $B>eItD:E@9bEY(B
    real(DP), parameter  :: z_Lower = 2.70d4      ! $B2<ItD:E@9bEY(B
    real(DP), parameter  :: c_Upper = 2.70d-2     ! $BBh(B 1 $B78?t(B
    real(DP), parameter  :: c_Lower = 3.60d-3     ! $BBh(B 2 $B78?t(B
    real(DP), parameter  :: s_Upper = 7.50d3      ! $BBh(B 1 $BI8=`JP:9(B
    real(DP), parameter  :: s_Lower = 1.30d4      ! $BBh(B 2 $BI8=`JP:9(B

    cal_Qsub =                                                                      &
      & (                                                                           &
      &     c_Lower * exp( - ( z0 - z_Lower )**2.0d0 / ( 2.0d0 * s_Lower **2.0d0) ) &
      &   + c_Upper * exp( - ( z0 - z_Upper )**2.0d0 / ( 2.0d0 * s_Upper **2.0d0) ) &
      &  ) 
    
    
  end function cal_Qsub

end module Radiation_baker1998
