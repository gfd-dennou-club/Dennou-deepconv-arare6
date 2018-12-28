!= Module advection_center4_mpdata_3d
!
! Authors::   $B?y;39L0lO/(B(SUGIYAMA Ko-ichiro)
! Version::   $Id: advection_center4_mpdata_3d.f90,v 1.3 2014/07/08 00:55:25 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module advection_center4_mpdata_3d
  !
  ! 3D $B7W;;MQ$N0\N.7W;;%b%8%e!<%k(B. 
  !
  !   $B0\N.(B: 4 $B<!Cf1{:9J,(B
  !   $B?tCM3H;6(B (4 $B3,(B): 2 $B<!Cf1{:9J,(B
  !
  ! $B%j!<%W%U%m%C%0$G(B, $B0\N.$rCf1{:9J,$G7W;;$9$k$?$a$K(B, 
  ! $B?tCM3H;69`$rDI2C$7$F$$$k(B. 
  ! 

  !$B%b%8%e!<%kFI$_9~$_(B
  !
  use dc_types,   only : DP

  !$B0EL[$N7?@k8@6X;_(B
  !
  implicit none

  !$BB0@-$N;XDj(B
  !
  private

  ! $BJQ?t$N@_Dj(B
  !
  real(DP), save :: NuHh  = 0.0d0         !$BG.$KBP$9$k?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save :: NuVh  = 0.0d0         !$BG.$KBP$9$k?tCMG4@-$N78?t(B ($B1tD>J}8~(B)
  real(DP), save :: NuHm  = 0.0d0         !$B1?F0NL$KBP$9$k?tCMG4@-$N78?t(B ($B?eJ?J}8~(B)
  real(DP), save :: NuVm  = 0.0d0         !$B1?F0NL$KBP$9$k?tCMG4@-$N78?t(B ($B1tD>J}8~(B)
  character(*), parameter:: module_name = 'advection_center4_mpdata_3d'
                                          ! $B%b%8%e!<%k$NL>>N(B.
                                          ! Module name
  !public
  !
  public advection_center4_mpdata_3d_init
  public advection_center4_mpdata_3d_dry
  public advection_center4_mpdata_3d_tracer

contains

  subroutine advection_center4_mpdata_3d_init( AlphaNDiff, NDiffRatio )
    !
    ! $B=i4|2=%k!<%A%s(B
    !

    ! $B%b%8%e!<%kFI$_9~$_(B
    !
    use timeset,    only : DelTimeLong
    use axesset,    only : dx, dy, dz       ! $B3J;R4V3V(B
    use dc_message, only : MessageNotify

    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none
    
    ! $BJQ?t$NDj5A(B
    !
    real(DP), intent(in) :: AlphaNDiff  !$B?tCM3H;6$N78?t(B. 
    real(DP), intent(in) :: NDiffRatio
    
    !-------------------------------------------------------------------
    ! $B?tCM3H;678?t$r7h$a$k(B
    !
    ! CReSS $B%^%K%e%"%k$N5-=R$K=>$C$F(B NuH, NuV $B$r7h$a$k(B.
    ! $B1?F0NL$HG.$KBP$9$k?tCM3H;6$NBg$-$5$rJQ$($i$l$k$h$&$K(B NDiffRatio $B$r>h$8$F$$$k(B.
    ! 
    NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    NuVh = AlphaNDiff * ( dz ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    
    NuHm = NuHh * NDiffRatio
    NuVm = NuVh * NDiffRatio

    !-------------------------------------------------------------------
    ! $B=PNO(B
    !
    call MessageNotify( "M", module_name, "NuHh = %f", d=(/NuHh/) )
    call MessageNotify( "M", module_name, "NuVh = %f", d=(/NuVh/) )
    call MessageNotify( "M", module_name, "NuHm = %f", d=(/NuHm/) )
    call MessageNotify( "M", module_name, "NuVm = %f", d=(/NuVm/) )

  end subroutine advection_center4_mpdata_3d_init

!!!------------------------------------------------------------------------!!!
 
  subroutine advection_center4_mpdata_3d_dry(  &
    & pyz_VelXB,    pyz_VelXN,          & ! (in)
    & xqz_VelYB,    xqz_VelYN,          & ! (in)
    & xyr_VelZB,    xyr_VelZN,          & ! (in)
    & xyz_PTempB,   xyz_PTempN,         & ! (in)
    & xyz_ExnerB,   xyz_ExnerN,         & ! (in)
    & xyz_KmB,      xyz_KmN,            & ! (in)
    & pyz_DVelXDtAdv,  pyz_VelXnDiff,      & !(out)
    & xqz_DVelYDtAdv,  xqz_VelYnDiff,      & !(out)
    & xyr_DVelZDtAdv,  xyr_VelZnDiff,      & !(out)
    & xyz_DPTempDtAdv, xyz_PTempNDiff,     & !(out)
    & xyz_DExnerDtAdv, xyz_ExnerNDiff,     & !(out)
    & xyz_DKmDtAdv,    xyz_KmNDiff         & !(out)
    & )
    ! 
    ! $B0\N.7W;;(B ($B4%Ag(B)
    !
    !   $B0\N.(B: 4 $B<!Cf1{:9J,(B
    !   $B?tCM3H;6(B (4 $B3,(B): 2 $B<!Cf1{:9J,(B
    !
    ! $B%j!<%W%U%m%C%0$G(B, $B0\N.$rCf1{:9J,$G7W;;$9$k$?$a$K(B, 
    ! $B?tCM3H;69`$rDI2C$7$F$$$k(B. 
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x $BJ}8~$NG[Ns$N2<8B(B
      &                  imax,            &! x $BJ}8~$NG[Ns$N>e8B(B
      &                  jmin,            &! y $BJ}8~$NG[Ns$N2<8B(B
      &                  jmax,            &! y $BJ}8~$NG[Ns$N>e8B(B
      &                  kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                  kmax              ! z $BJ}8~$NG[Ns$N>e8B(B
    use axesset,  only : dx, dy, dz        ! $B3J;R4V3V(B
    use basicset, only : xyz_PTempBZ,    & ! $B290L$N4pK\>l(B
      &                  xyz_ExnerBZ       ! $B290L$N4pK\>l(B

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BG[Ns$NDj5A(B
    !
    real(DP), intent(in)    :: pyz_VelXB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmN(imin:imax,jmin:jmax,kmin:kmax)

    real(DP), intent(out)   :: pyz_DVelXDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_DVelYDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_DVelZDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DKmDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DExnerDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_DPTempDtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: pyz_VelXnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZnDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_KmNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerNDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_PTempNDiff(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)                :: xyz_VarB(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_VarN(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)                :: fct1, fct2
    integer                 :: i, j, k

      
    ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
    !
    fct1 = 9.0d0 / 8.0d0
    fct2 = 1.0d0 / 24.0d0

    !---------------------------------------------------------------------
    ! $BB.EY(B U
    ! 

    ! $B0\N.(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          pyz_DVelXDtAdv(i,j,k) =                                             &
            & - pyz_VelXN(i,j,k)                                           &
            &   * (                                                        &
            &       fct1 * (   pyz_VelXN(i+1,j,k) - pyz_VelXN(i-1,j,k) )   &
            &     - fct2 * (   pyz_VelXN(i+2,j,k) + pyz_VelXN(i+1,j,k)     &
            &                - pyz_VelXN(i-1,j,k) - pyz_VelXN(i-2,j,k) )   &
            &     ) * 5.0d-1 / dx                                          &
            & - (                                                          &
            &   + ( xqz_VelYN(i+1,j,k) + xqz_VelYN(i,j,k) )                &  
            &     * (                                                      &  
            &         fct1 * ( pyz_VelXN(i,j+1,k) - pyz_VelXN(i,j,k)   )   &
            &       - fct2 * ( pyz_VelXN(i,j+2,k) - pyz_VelXN(i,j-1,k) )   &
            &       )                                                      &
            &   + ( xqz_VelYN(i+1,j-1,k) + xqz_VelYN(i,j-1,k) )            & 
            &     * (                                                      &
            &         fct1 * ( pyz_VelXN(i,j,k)   - pyz_VelXN(i,j-1,k) )   &
            &       - fct2 * ( pyz_VelXN(i,j+1,k) - pyz_VelXN(i,j-2,k) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dy                                            &
            & - (                                                          &
            &   + ( xyr_VelZN(i+1,j,k) + xyr_VelZN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( pyz_VelXN(i,j,k+1) - pyz_VelXN(i,j,k)   )   &
            &       - fct2 * ( pyz_VelXN(i,j,k+2) - pyz_VelXN(i,j,k-1) )   &
            &       )                                                      &
            &   + ( xyr_VelZN(i+1,j,k-1) + xyr_VelZN(i,j,k-1) )            & 
            &     * (                                                      &
            &         fct1 * ( pyz_VelXN(i,j,k)   - pyz_VelXN(i,j,k-1) )   &
            &       - fct2 * ( pyz_VelXN(i,j,k+1) - pyz_VelXN(i,j,k-2) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dz
          
        end do
      end do
    end do

    ! $B?tCM3H;6(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          pyz_VelXnDiff(i,j,k) =                        &
            & - (                                       &
            &     + pyz_VelXB(i+2,j,k)                  &
            &     + pyz_VelXB(i-2,j,k)                  &
            &     - pyz_VelXB(i+1,j,k) * 4.0d0          &
            &     - pyz_VelXB(i-1,j,k) * 4.0d0          &
            &     + pyz_VelXB(i  ,j,k) * 6.0d0          &
            &   ) * NuHm / ( dx ** 4.0d0 )              &
            & - (                                       &
            &     + pyz_VelXB(i,j+2,k)                  &
            &     + pyz_VelXB(i,j-2,k)                  &
            &     - pyz_VelXB(i,j+1,k) * 4.0d0          &
            &     - pyz_VelXB(i,j-1,k) * 4.0d0          &
            &     + pyz_VelXB(i,j  ,k) * 6.0d0          &
            &   ) * NuHm / ( dy ** 4.0d0 )              &
            & - (                                       & 
            &     + pyz_VelXB(i,j,k+2)                  &
            &     + pyz_VelXB(i,j,k-2)                  &
            &     - pyz_VelXB(i,j,k+1) * 4.0d0          &
            &     - pyz_VelXB(i,j,k-1) * 4.0d0          &
            &     + pyz_VelXB(i,j,k  ) * 6.0d0          &
            &   ) * NuVm / ( dz ** 4.0d0 )
          
        end do
      end do
    end do


    ! $BCM$N3NDj(B
    !
    pyz_DVelXDtAdv(imin:imin+1,:,:) = 1.0d10
    pyz_DVelXDtAdv(imax-1:imax,:,:) = 1.0d10
    pyz_DVelXDtAdv(:,jmin:jmin+1,:) = 1.0d10
    pyz_DVelXDtAdv(:,jmax-1:jmax,:) = 1.0d10
    pyz_DVelXDtAdv(:,:,kmin:kmin+1) = 1.0d10
    pyz_DVelXDtAdv(:,:,kmax-1:kmax) = 1.0d10

    pyz_VelXnDiff(imin:imin+1,:,:) = 1.0d10
    pyz_VelXnDiff(imax-1:imax,:,:) = 1.0d10
    pyz_VelXnDiff(:,jmin:jmin+1,:) = 1.0d10
    pyz_VelXnDiff(:,jmax-1:jmax,:) = 1.0d10
    pyz_VelXnDiff(:,:,kmin:kmin+1) = 1.0d10
    pyz_VelXnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! $BB.EY(B V
    !       

    ! $B0\N.(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xqz_DVelYDtAdv(i,j,k) =                                             &
            & - (                                                          &
            &   + ( pyz_VelXN(i,j+1,k) + pyz_VelXN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( xqz_VelYN(i+1,j,k) - xqz_VelYN(i,j,k)   )   &
            &       - fct2 * ( xqz_VelYN(i+2,j,k) - xqz_VelYN(i-1,j,k) )   &
            &       )                                                      &
            &   + ( pyz_VelXN(i-1,j+1,k) + pyz_VelXN(i-1,j,k) )            & 
            &     * (                                                      &
            &         fct1 * ( xqz_VelYN(i,j,k)   - xqz_VelYN(i-1,j,k) )   &
            &       - fct2 * ( xqz_VelYN(i+1,j,k) - xqz_VelYN(i-2,j,k) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dx                                            &
            & - xqz_VelYN(i,j,k)                                           &
            &   * (                                                        &
            &       fct1 * (   xqz_VelYN(i,j+1,k) - xqz_VelYN(i,j-1,k) )   &
            &     - fct2 * (   xqz_VelYN(i,j+2,k) + xqz_VelYN(i,j+1,k)     &
            &                - xqz_VelYN(i,j-1,k) - xqz_VelYN(i,j-2,k) )   &
            &     ) * 5.0d-1 / dy                                          &
            & - (                                                          &
            &   + ( xyr_VelZN(i,j+1,k) + xyr_VelZN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( xqz_VelYN(i,j,k+1) - xqz_VelYN(i,j,k)   )   &
            &       - fct2 * ( xqz_VelYN(i,j,k+2) - xqz_VelYN(i,j,k-1) )   &
            &       )                                                      &
            &   + ( xyr_VelZN(i,j+1,k-1) + xyr_VelZN(i,j,k-1) )            & 
            &     * (                                                      &
            &         fct1 * ( xqz_VelYN(i,j,k)   - xqz_VelYN(i,j,k-1) )   &
            &       - fct2 * ( xqz_VelYN(i,j,k+1) - xqz_VelYN(i,j,k-2) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dz
          
        end do
      end do
    end do

    ! $B?tCM3H;6(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xqz_VelYnDiff(i,j,k) =                        &
            & - (                                       &
            &     + xqz_VelYB(i+2,j,k)                  &
            &     + xqz_VelYB(i-2,j,k)                  &
            &     - xqz_VelYB(i+1,j,k) * 4.0d0          &
            &     - xqz_VelYB(i-1,j,k) * 4.0d0          &
            &     + xqz_VelYB(i  ,j,k) * 6.0d0          &
            &   ) * NuHm / ( dx ** 4.0d0 )              &
            & - (                                       &
            &     + xqz_VelYB(i,j+2,k)                  &
            &     + xqz_VelYB(i,j-2,k)                  &
            &     - xqz_VelYB(i,j+1,k) * 4.0d0          &
            &     - xqz_VelYB(i,j-1,k) * 4.0d0          &
            &     + xqz_VelYB(i,j  ,k) * 6.0d0          &
            &   ) * NuHm / ( dy ** 4.0d0 )              &
            & - (                                       & 
            &     + xqz_VelYB(i,j,k+2)                  &
            &     + xqz_VelYB(i,j,k-2)                  &
            &     - xqz_VelYB(i,j,k+1) * 4.0d0          &
            &     - xqz_VelYB(i,j,k-1) * 4.0d0          &
            &     + xqz_VelYB(i,j,k  ) * 6.0d0          &
            &   ) * NuVm / ( dz ** 4.0d0 )

        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xqz_DVelYDtAdv(imin:imin+1,:,:) = 1.0d10
    xqz_DVelYDtAdv(imax-1:imax,:,:) = 1.0d10
    xqz_DVelYDtAdv(:,jmin:jmin+1,:) = 1.0d10
    xqz_DVelYDtAdv(:,jmax-1:jmax,:) = 1.0d10
    xqz_DVelYDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xqz_DVelYDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xqz_VelYnDiff(imin:imin+1,:,:) = 1.0d10
    xqz_VelYnDiff(imax-1:imax,:,:) = 1.0d10
    xqz_VelYnDiff(:,jmin:jmin+1,:) = 1.0d10
    xqz_VelYnDiff(:,jmax-1:jmax,:) = 1.0d10
    xqz_VelYnDiff(:,:,kmin:kmin+1) = 1.0d10
    xqz_VelYnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! $BB.EY(B W
    !     

    ! $B0\N.(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyr_DVelZDtAdv(i,j,k) =                                             &
            & - (                                                          &
            &   + ( pyz_VelXN(i,j,k+1) + pyz_VelXN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i+1,j,k) - xyr_VelZN(i,j,k)   )   &
            &       - fct2 * ( xyr_VelZN(i+2,j,k) - xyr_VelZN(i-1,j,k) )   &
            &       )                                                      &
            &   + ( pyz_VelXN(i-1,j,k+1) + pyz_VelXN(i-1,j,k) )            & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i,j,k)   - xyr_VelZN(i-1,j,k) )   &
            &       - fct2 * ( xyr_VelZN(i+1,j,k) - xyr_VelZN(i-2,j,k) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dx                                            &
            & - (                                                          &
            &   + ( xqz_VelYN(i,j,k+1) + xqz_VelYN(i,j,k) )                & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i,j+1,k) - xyr_VelZN(i,j,k)   )   &
            &       - fct2 * ( xyr_VelZN(i,j+2,k) - xyr_VelZN(i,j-1,k) )   &
            &       )                                                      &
            &   + ( xqz_VelYN(i,j-1,k+1) + xqz_VelYN(i,j-1,k) )            & 
            &     * (                                                      &
            &         fct1 * ( xyr_VelZN(i,j,k)   - xyr_VelZN(i,j-1,k) )   &
            &       - fct2 * ( xyr_VelZN(i,j+1,k) - xyr_VelZN(i,j-2,k) )   &
            &       )                                                      &
            &   ) * 2.5d-1 / dy                                            &
            & - xyr_VelZN(i,j,k)                                           &
            &   * (                                                        &
            &       fct1 * (   xyr_VelZN(i,j,k+1) - xyr_VelZN(i,j,k-1) )   &
            &     - fct2 * (   xyr_VelZN(i,j,k+2) + xyr_VelZN(i,j,k+1)     &
            &                - xyr_VelZN(i,j,k-1) - xyr_VelZN(i,j,k-2) )   &
            &     ) * 5.0d-1 / dz

        end do
      end do
    end do
      
    ! $B?tCM3H;6(B
    !
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyr_VelZnDiff(i,j,k) =                        &
            & - (                                       &
            &     + xyr_VelZB(i+2,j,k)                  &
            &     + xyr_VelZB(i-2,j,k)                  &
            &     - xyr_VelZB(i+1,j,k) * 4.0d0          &
            &     - xyr_VelZB(i-1,j,k) * 4.0d0          &
            &     + xyr_VelZB(i  ,j,k) * 6.0d0          &
            &   ) * NuHm / ( dx ** 4.0d0 )              &
            & - (                                       &
            &     + xyr_VelZB(i,j+2,k)                  &
            &     + xyr_VelZB(i,j-2,k)                  &
            &     - xyr_VelZB(i,j+1,k) * 4.0d0          &
            &     - xyr_VelZB(i,j-1,k) * 4.0d0          &
            &     + xyr_VelZB(i,j  ,k) * 6.0d0          &
            &   ) * NuHm / ( dy ** 4.0d0 )              &
            & - (                                       &
            &     + xyr_VelZB(i,j,k+2)                  &
            &     + xyr_VelZB(i,j,k-2)                  &
            &     - xyr_VelZB(i,j,k+1) * 4.0d0          &
            &     - xyr_VelZB(i,j,k-1) * 4.0d0          &
            &     + xyr_VelZB(i,j,k  ) * 6.0d0          &
            &   ) * NuVm / ( dz ** 4.0d0 )

        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xyr_DVelZDtAdv(imin:imin+1,:,:) = 1.0d10
    xyr_DVelZDtAdv(imax-1:imax,:,:) = 1.0d10
    xyr_DVelZDtAdv(:,jmin:jmin+1,:) = 1.0d10
    xyr_DVelZDtAdv(:,jmax-1:jmax,:) = 1.0d10
    xyr_DVelZDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyr_DVelZDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyr_VelZnDiff(imin:imin+1,:,:) = 1.0d10
    xyr_VelZnDiff(imax-1:imax,:,:) = 1.0d10
    xyr_VelZnDiff(:,jmin:jmin+1,:) = 1.0d10
    xyr_VelZnDiff(:,jmax-1:jmax,:) = 1.0d10
    xyr_VelZnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyr_VelZnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! $B290L(B
    !     

    ! $BA4NL$r7W;;(B
    !
    xyz_VarB = xyz_PTempB
    xyz_VarN = xyz_PTempN + xyz_PTempBZ


    fs_FluxX = 0.0d0 ; sf_FluxZ = 0.0d0 ; ss_AdvScalar = 0.0d0

    ! $B>eN.:9J,K!$rMQ$$$F%9%+%i!<$N%U%i%C%/%9$r7W;;(B
    !
    pyz_FluxX( imin:imax-1, :, : ) =                  &
      &  (                                            &
      &   + (                                         &
      &      +      pyz_VelX( imin:imax-1, :, : )     &
      &      + abs( pyz_VelX( imin:imax-1, :, : ) )   &
      &     )  *    xyz_VarN( imin:imax-1, :, : )     &
      &   + (                                         &
      &      +      pyz_VelX( imin:imax-1, :, : )     &
      &      - abs( pyz_VelX( imin:imax-1, :, : ) )   &
      &     )  *    xyz_VarN( imin+1:imax, :, : )     &
      &  ) * 0.5d0
    pyz_FluxX( imax, :, : ) = 0.0d0

    xqz_FluxY( :, jmin:jmax-1, : ) =                  &
      &  (                                            &
      &   + (                                         &
      &      +      xqz_VelY( :, jmin:jmax-1, : )     &
      &      + abs( xqz_VelY( :, jmin:jmax-1, : ) )   &
      &     )  *    xyz_VarN( :, jmin:jmax-1, : )     &
      &   + (                                         &
      &      +      xqz_VelY( :, jmin:jmax-1, : )     &
      &      - abs( xqz_VelY( :, jmin:jmax-1, : ) )   &
      &     )  *    xyz_VarN( :, jmin+1:jmax, : )     &
      &  ) * 0.5d0
    xqz_FluxY( :, jmax, : ) = 0.0d0

    xyr_FluxZ( :, :, kmin:kmax-1 ) =                  &
      &  (                                            &
      &   + (                                         &
      &      +      xyr_VelZ( :, :, kmin:kmax-1 )     &
      &      + abs( xyr_VelZ( :, :, kmin:kmax-1 ) )   &
      &     )  *    xyz_VarN( :, :, kmin:kmax-1 )     &
      &   + (                                         &
      &      +      xyr_VelZ( :, :, kmin:kmax-1 )     &
      &      - abs( xyr_VelZ( :, :, kmin:kmax-1 ) )   &
      &     )  *    xyz_VarN( :, :, kmin+1:kmax )     &
      &  ) * 0.5d0
    xyr_FluxZ( :, :, kmax ) = 0.0d0
    
    ! $BB.EY$N<}B+$r7W;;(B
    !
    do k = kmin + 1, kmax 
      do j = jmin + 1, jmax 
        do i = imin + 1, imax 

          xyz_VelDivN(i,j,k) =         &
            & + (                      &
            &     pyz_VelXN(i,   j, k) &
            &   - pyz_VelXN(i-1, j, k) &
            &   ) / dx                 &
            & + (                      &
            &     xqz_VelYN(i, j,   k) &
            &   - xqz_VelYN(i, j-1, k) &
            &   ) / dy                 &
            & + (                      &
            &     xyr_VelZN(i, j, k)   &
            &   - xyr_VelZN(i, j, k-1) &
            &   ) / dz

        end do
      end do
    end do

    ! $BCM$r3NDj$5$;$k(B
    !
    xyz_VelDivN(imin,:,:) = 1.0d10
    xyz_VelDivN(:,jmin,:) = 1.0d10  
    xyz_VelDivN(:,:,kmin) = 1.0d10

    
    ! $B0\N.$K$h$kA}J,$r7W;;(B
    !
    ss_AdvScalar = ss_dx_fs(fs_FluxX) + ss_dz_sf(sf_FluxZ)  !&
!      &                - ss_Scalar * ss_DivVel

    ! $B%9%+%i!<$N2>CM$r7W;;(B
    !
    ss_Scalar_tmp = ss_Scalar + DelTime * ( - ss_AdvScalar )
    







    ! $B0\N.(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_DPTempDtAdv(i,j,k) =                                           &
            & - (                                                         &
            &      pyz_VelXN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i+2,j,k) - xyz_VarN(i-1,j,k) ) &
            &          )                                                  &
            &    + pyz_VelXN(i-1,j,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i-1,j,k) ) &
            &          - fct2 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i-2,j,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dx                                           &
            & - (                                                         &
            &      xqz_VelYN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j+2,k) - xyz_VarN(i,j-1,k) ) &
            &          )                                                  &
            &    + xqz_VelYN(i,j-1,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j-1,k) ) &
            &          - fct2 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j-2,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dy                                           &
            & - (                                                         &
            &      xyr_VelZN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j,k+2) - xyz_VarN(i,j,k-1) ) &
            &          )                                                  &
            &    + xyr_VelZN(i,j,k-1)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j,k-1) ) &
            &          - fct2 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k-2) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dz

        end do
      end do
    end do
    
    ! $B?tCM3H;6(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_PTempNDiff(i,j,k) =                      &
            & - (                                      &
            &     + xyz_VarB(i+2,j,k)                  &
            &     + xyz_VarB(i-2,j,k)                  &
            &     - xyz_VarB(i+1,j,k) * 4.0d0          &
            &     - xyz_VarB(i-1,j,k) * 4.0d0          &
            &     + xyz_VarB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j+2,k)                  &
            &     + xyz_VarB(i,j-2,k)                  &
            &     - xyz_VarB(i,j+1,k) * 4.0d0          &
            &     - xyz_VarB(i,j-1,k) * 4.0d0          &
            &     + xyz_VarB(i,j  ,k) * 6.0d0          &
            &   ) * NuHh / ( dy ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j,k+2)                  &
            &     + xyz_VarB(i,j,k-2)                  &
            &     - xyz_VarB(i,j,k+1) * 4.0d0          &
            &     - xyz_VarB(i,j,k-1) * 4.0d0          &
            &     + xyz_VarB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )

        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xyz_DPTempDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DPTempDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DPTempDtAdv(:,jmin:jmin+1,:) = 1.0d10
    xyz_DPTempDtAdv(:,jmax-1:jmax,:) = 1.0d10
    xyz_DPTempDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DPTempDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_PTempnDiff(imin:imin+1,:,:) = 1.0d10
    xyz_PTempnDiff(imax-1:imax,:,:) = 1.0d10
    xyz_PTempnDiff(:,jmin:jmin+1,:) = 1.0d10
    xyz_PTempnDiff(:,jmax-1:jmax,:) = 1.0d10
    xyz_PTempnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_PTempnDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! $B%(%/%9%J!<4X?t(B
    !     

    ! $BA4NL$r7W;;(B
    !
    xyz_VarB = xyz_ExnerB
    xyz_VarN = xyz_ExnerN + xyz_ExnerBZ

    ! $B0\N.(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_DExnerDtAdv(i,j,k) =                                           &
            & - (                                                         &
            &      pyz_VelXN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i+2,j,k) - xyz_VarN(i-1,j,k) ) &
            &          )                                                  &
            &    + pyz_VelXN(i-1,j,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i-1,j,k) ) &
            &          - fct2 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i-2,j,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dx                                           &
            & - (                                                         &
            &      xqz_VelYN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j+2,k) - xyz_VarN(i,j-1,k) ) &
            &          )                                                  &
            &    + xqz_VelYN(i,j-1,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j-1,k) ) &
            &          - fct2 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j-2,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dy                                           &
            & - (                                                         &
            &      xyr_VelZN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j,k+2) - xyz_VarN(i,j,k-1) ) &
            &          )                                                  &
            &    + xyr_VelZN(i,j,k-1)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j,k-1) ) &
            &          - fct2 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k-2) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dz
          
        end do
      end do
    end do
    
    ! $B?tCM3H;6(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_ExnerNDiff(i,j,k) =                      &
            & - (                                      &
            &     + xyz_VarB(i+2,j,k)                  &
            &     + xyz_VarB(i-2,j,k)                  &
            &     - xyz_VarB(i+1,j,k) * 4.0d0          &
            &     - xyz_VarB(i-1,j,k) * 4.0d0          &
            &     + xyz_VarB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j+2,k)                  &
            &     + xyz_VarB(i,j-2,k)                  &
            &     - xyz_VarB(i,j+1,k) * 4.0d0          &
            &     - xyz_VarB(i,j-1,k) * 4.0d0          &
            &     + xyz_VarB(i,j  ,k) * 6.0d0          &
            &   ) * NuHh / ( dy ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j,k+2)                  &
            &     + xyz_VarB(i,j,k-2)                  &
            &     - xyz_VarB(i,j,k+1) * 4.0d0          &
            &     - xyz_VarB(i,j,k-1) * 4.0d0          &
            &     + xyz_VarB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )

        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xyz_DExnerDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DExnerDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DExnerDtAdv(:,jmin:jmin+1,:) = 1.0d10
    xyz_DExnerDtAdv(:,jmax-1:jmax,:) = 1.0d10
    xyz_DExnerDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DExnerDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_ExnernDiff(imin:imin+1,:,:) = 1.0d10
    xyz_ExnernDiff(imax-1:imax,:,:) = 1.0d10
    xyz_ExnernDiff(:,jmin:jmin+1,:) = 1.0d10
    xyz_ExnernDiff(:,jmax-1:jmax,:) = 1.0d10
    xyz_ExnernDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_ExnernDiff(:,:,kmax-1:kmax) = 1.0d10

    !---------------------------------------------------------------------
    ! $BMpN.3H;678?t(B
    !     

    ! $B4pK\>l$J$7(B
    !
    xyz_VarB = xyz_KmB 
    xyz_VarN = xyz_KmN 

    ! $B0\N.(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_DKmDtAdv(i,j,k) =                                              &
            & - (                                                         &
            &      pyz_VelXN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i+2,j,k) - xyz_VarN(i-1,j,k) ) &
            &          )                                                  &
            &    + pyz_VelXN(i-1,j,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i-1,j,k) ) &
            &          - fct2 * ( xyz_VarN(i+1,j,k) - xyz_VarN(i-2,j,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dx                                           &
            & - (                                                         &
            &      xqz_VelYN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j+2,k) - xyz_VarN(i,j-1,k) ) &
            &          )                                                  &
            &    + xqz_VelYN(i,j-1,k)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j-1,k) ) &
            &          - fct2 * ( xyz_VarN(i,j+1,k) - xyz_VarN(i,j-2,k) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dy                                           &
            & - (                                                         &
            &      xyr_VelZN(i,j,k)                                       &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k)   ) &
            &          - fct2 * ( xyz_VarN(i,j,k+2) - xyz_VarN(i,j,k-1) ) &
            &          )                                                  &
            &    + xyr_VelZN(i,j,k-1)                                     &
            &        * (                                                  &
            &            fct1 * ( xyz_VarN(i,j,k)   - xyz_VarN(i,j,k-1) ) &
            &          - fct2 * ( xyz_VarN(i,j,k+1) - xyz_VarN(i,j,k-2) ) &
            &          )                                                  &
            &   ) * 5.0d-1 / dz

        end do
      end do
    end do
    
    ! $B?tCM3H;6(B
    ! 
    do k = kmin + 2, kmax - 2
      do j = jmin + 2, jmax - 2
        do i = imin + 2, imax - 2
          
          xyz_KmNDiff(i,j,k) =                         &
            & - (                                      &
            &     + xyz_VarB(i+2,j,k)                  &
            &     + xyz_VarB(i-2,j,k)                  &
            &     - xyz_VarB(i+1,j,k) * 4.0d0          &
            &     - xyz_VarB(i-1,j,k) * 4.0d0          &
            &     + xyz_VarB(i  ,j,k) * 6.0d0          &
            &   ) * NuHh / ( dx ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j+2,k)                  &
            &     + xyz_VarB(i,j-2,k)                  &
            &     - xyz_VarB(i,j+1,k) * 4.0d0          &
            &     - xyz_VarB(i,j-1,k) * 4.0d0          &
            &     + xyz_VarB(i,j  ,k) * 6.0d0          &
            &   ) * NuHh / ( dy ** 4.0d0 )             &
            & - (                                      &
            &       xyz_VarB(i,j,k+2)                  &
            &     + xyz_VarB(i,j,k-2)                  &
            &     - xyz_VarB(i,j,k+1) * 4.0d0          &
            &     - xyz_VarB(i,j,k-1) * 4.0d0          &
            &     + xyz_VarB(i,j,k  ) * 6.0d0          &
            &   ) * NuVh / ( dz ** 4.0d0 )

        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xyz_DKmDtAdv(imin:imin+1,:,:) = 1.0d10
    xyz_DKmDtAdv(imax-1:imax,:,:) = 1.0d10
    xyz_DKmDtAdv(:,jmin:jmin+1,:) = 1.0d10
    xyz_DKmDtAdv(:,jmax-1:jmax,:) = 1.0d10
    xyz_DKmDtAdv(:,:,kmin:kmin+1) = 1.0d10
    xyz_DKmDtAdv(:,:,kmax-1:kmax) = 1.0d10

    xyz_KmnDiff(imin:imin+1,:,:) = 1.0d10
    xyz_KmnDiff(imax-1:imax,:,:) = 1.0d10
    xyz_KmnDiff(:,jmin:jmin+1,:) = 1.0d10
    xyz_KmnDiff(:,jmax-1:jmax,:) = 1.0d10
    xyz_KmnDiff(:,:,kmin:kmin+1) = 1.0d10
    xyz_KmnDiff(:,:,kmax-1:kmax) = 1.0d10

  end subroutine advection_center4_mpdata_3d_dry


!!!------------------------------------------------------------------------!!!
 
  subroutine advection_center4_mpdata_3d_tracer( &
    & pyz_VelXN, xqz_VelYN, xyr_VelZN,   & ! (in)
    & xyzf_QMixB,   xyzf_QMixN,          & ! (in)
    & xyzf_QMixAdv, xyzf_QMixNDiff       & !(out)
    & )
    ! 
    ! $B0\N.7W;;(B ($BJ*<A(B)
    !
    !   $B0\N.(B: 4 $B<!Cf1{:9J,(B
    !   $B?tCM3H;6(B (4 $B3,(B): 2 $B<!Cf1{:9J,(B
    !
    ! $B%j!<%W%U%m%C%0$G(B, $B0\N.$rCf1{:9J,$G7W;;$9$k$?$a$K(B, 
    ! $B?tCM3H;69`$rDI2C$7$F$$$k(B. 
    !

    !$B%b%8%e!<%kFI$_9~$_(B
    !
    use dc_types, only: DP
    use gridset, only :  imin,            &! x $BJ}8~$NG[Ns$N2<8B(B
      &                  imax,            &! x $BJ}8~$NG[Ns$N>e8B(B
      &                  jmin,            &! y $BJ}8~$NG[Ns$N2<8B(B
      &                  jmax,            &! y $BJ}8~$NG[Ns$N>e8B(B
      &                  kmin,            &! z $BJ}8~$NG[Ns$N2<8B(B
      &                  kmax,            &! z $BJ}8~$NG[Ns$N>e8B(B
      &                  ncmax
    use axesset,  only : dx, dy, dz        ! $B3J;R4V3V(B
    use basicset, only : xyzf_QMixBZ       ! $B:.9gHf$N4pK\>l(B

    ! $B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    ! $BG[Ns$NDj5A(B
    !
    real(DP), intent(in)    :: pyz_VelXN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYN(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZN(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyzf_QMixB(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixN(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)

    real(DP), intent(out)   :: xyzf_QMixAdv(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out)   :: xyzf_QMixNDiff(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)

    real(DP)                :: xyzf_VarB(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP)                :: xyzf_VarN(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)

    real(DP)                :: fct1, fct2
    integer                 :: i, j, k, s

      
    ! $BHyJ,$KMQ$$$k78?t$rM=$a7W;;(B
    !
    fct1 = 9.0d0 / 8.0d0
    fct2 = 1.0d0 / 24.0d0

    !---------------------------------------------------------------------
    ! $B:.9gHf(B
    !     

    ! $BA4NL$r7W;;(B
    !
    xyzf_VarB = xyzf_QMixB
    xyzf_VarN = xyzf_QMixN + xyzf_QMixBZ
   
    ! $B0\N.(B
    ! 
    do s = 1, ncmax
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyzf_QMixAdv(i,j,k,s) =                                               &
              & - (                                                               &
              &      pyz_VelXN(i,j,k)                                             &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i+1,j,k,s) - xyzf_VarN(i,j,k,s)   ) &
              &          - fct2 * ( xyzf_VarN(i+2,j,k,s) - xyzf_VarN(i-1,j,k,s) ) &
              &          )                                                        &
              &    + pyz_VelXN(i-1,j,k)                                           &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i,j,k,s)   - xyzf_VarN(i-1,j,k,s) ) &
              &          - fct2 * ( xyzf_VarN(i+1,j,k,s) - xyzf_VarN(i-2,j,k,s) ) &
              &          )                                                        &
              &   ) * 5.0d-1 / dx                                                 &
              & - (                                                               &
              &      xqz_VelYN(i,j,k)                                             &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i,j+1,k,s) - xyzf_VarN(i,j,k,s)   ) &
              &          - fct2 * ( xyzf_VarN(i,j+2,k,s) - xyzf_VarN(i,j-1,k,s) ) &
              &          )                                                        &
              &    + xqz_VelYN(i,j-1,k)                                           &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i,j,k,s)   - xyzf_VarN(i,j-1,k,s) ) &
              &          - fct2 * ( xyzf_VarN(i,j+1,k,s) - xyzf_VarN(i,j-2,k,s) ) &
              &          )                                                        &
              &   ) * 5.0d-1 / dy                                                 &
              & - (                                                               &
              &      xyr_VelZN(i,j,k)                                             &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i,j,k+1,s) - xyzf_VarN(i,j,k,s)   ) &
              &          - fct2 * ( xyzf_VarN(i,j,k+2,s) - xyzf_VarN(i,j,k-1,s) ) &
              &          )                                                        &
              &    + xyr_VelZN(i,j,k-1)                                           &
              &        * (                                                        &
              &            fct1 * ( xyzf_VarN(i,j,k,s)   - xyzf_VarN(i,j,k-1,s) ) &
              &          - fct2 * ( xyzf_VarN(i,j,k+1,s) - xyzf_VarN(i,j,k-2,s) ) &
              &          )                                                        &
              &   ) * 5.0d-1 / dz

          end do
        end do
      end do
    end do
      
    ! $B?tCM3H;6(B
    ! 
    do s = 1, ncmax
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyzf_QMixNDiff(i,j,k,s) =                       &
              & - (                                         &
              &       xyzf_VarB(i+2,j,k,s)                  &
              &     + xyzf_VarB(i-2,j,k,s)                  &
              &     - xyzf_VarB(i+1,j,k,s) * 4.0d0          &
              &     - xyzf_VarB(i-1,j,k,s) * 4.0d0          &
              &     + xyzf_VarB(i  ,j,k,s) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )                &
              & - (                                         &
              &       xyzf_VarB(i,j+2,k,s)                  &
              &     + xyzf_VarB(i,j-2,k,s)                  &
              &     - xyzf_VarB(i,j+1,k,s) * 4.0d0          &
              &     - xyzf_VarB(i,j-1,k,s) * 4.0d0          &
              &     + xyzf_VarB(i,j  ,k,s) * 6.0d0          &
              &   ) * NuHh / ( dy ** 4.0d0 )                &
              & - (                                         &
              &       xyzf_VarB(i,j,k+2,s)                  &
              &     + xyzf_VarB(i,j,k-2,s)                  &
              &     - xyzf_VarB(i,j,k+1,s) * 4.0d0          &
              &     - xyzf_VarB(i,j,k-1,s) * 4.0d0          &
              &     + xyzf_VarB(i,j,k  ,s) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
            
          end do
        end do
      end do
    end do

    ! $BCM$N3NDj(B
    !
    xyzf_QMixAdv(imin:imin+1,:,:,:) = 1.0d10
    xyzf_QMixAdv(imax-1:imax,:,:,:) = 1.0d10
    xyzf_QMixAdv(:,jmin:jmin+1,:,:) = 1.0d10
    xyzf_QMixAdv(:,jmax-1:jmax,:,:) = 1.0d10
    xyzf_QMixAdv(:,:,kmin:kmin+1,:) = 1.0d10
    xyzf_QMixAdv(:,:,kmax-1:kmax,:) = 1.0d10

    xyzf_QMixnDiff(imin:imin+1,:,:,:) = 1.0d10
    xyzf_QMixnDiff(imax-1:imax,:,:,:) = 1.0d10
    xyzf_QMixnDiff(:,jmin:jmin+1,:,:) = 1.0d10
    xyzf_QMixnDiff(:,jmax-1:jmax,:,:) = 1.0d10
    xyzf_QMixnDiff(:,:,kmin:kmin+1,:) = 1.0d10
    xyzf_QMixnDiff(:,:,kmax-1:kmax,:) = 1.0d10

  end subroutine advection_center4_mpdata_3d_tracer
  
end module advection_center4_mpdata_3d
