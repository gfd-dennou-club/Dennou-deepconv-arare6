!= Module FillNegative
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: fillnegative.f90,v 1.16 2014/07/08 01:06:28 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module FillNegative
  !
  !$BIi$N1@?eNL$J$I$r7jKd$a$9$k$?$a$N%k!<%A%s(B
  !  * $BIi$H$J$C$?>l9g$K$O(B, $B<~0O$N(B 8 $B3J;RE@$NCM$r:o$C$F7jKd$a$9$k(B.
  !  * F77$BHG(B deepconv $B$N(B QFILL.f $B$r2~JQ(B 

  !$B%b%8%e!<%kFI$_9~$_(B
  use dc_types,   only: DP

  !$B0EL[$N7?@k8@6X;_(B
  implicit none

  !$B8x3+MWAG(B
  public FillNegative_Init
  public xyza_FillNegative_xyza
  public xyz_FillNegativeDensity_xyz

  !$BJQ?tDj5A(B
  real(DP), save, private, allocatable  :: xyza_Basic(:,:,:,:)
  real(DP), save, private, allocatable  :: xyz_Dens(:,:,:)
  real(DP), save, private, allocatable  :: xza_Basic(:,:,:)
  real(DP), save, private, allocatable  :: xz_Dens(:,:)

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine FillNegative_Init( )    
    !
    ! $B4p=`CM$r<hF@(B. 
    ! $B4pK\>l$H>qMp@.J,$NOB$,%<%m$h$j$b>.$5$/$J$k$3$H$O5vMF$7$J$$(B
    !
    
    !$B%b%8%e!<%kFI$_9~$_(B
    use gtool_historyauto, &
      &              only: HistoryAutoAddVariable
    use gridset,     only: imin,  &! $BG[Ns$N(B X $BJ}8~$N2<8B(B
      &                    imax,  &! $BG[Ns$N(B X $BJ}8~$N>e8B(B
      &                    jmin,  &! $BG[Ns$N(B Y $BJ}8~$N2<8B(B
      &                    jmax,  &! $BG[Ns$N(B Y $BJ}8~$N>e8B(B
      &                    kmin,  &! $BG[Ns$N(B Z $BJ}8~$N2<8B(B
      &                    kmax,  &! $BG[Ns$N(B Z $BJ}8~$N>e8B(B
      &                    ncmax
    use basicset,    only: xyzf_QMixBZ, xyz_DensBZ
    use composition, only: SpcWetSymbol

    !$B0EL[$N7?@k8@6X;_(B
    implicit none

    !$BF~NOJQ?t(B
    integer  :: l

    !$B=i4|2=(B
    allocate( xyza_Basic(imin:imax,jmin:jmax,kmin:kmax, ncmax) )
    allocate( xyz_Dens  (imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xza_Basic(imin:imax,kmin:kmax, ncmax) )
    allocate( xz_Dens(imin:imax,kmin:kmax ) )
    
    !$BCM$NBeF~(B
    xyza_Basic = xyzf_QMixBZ
    xyz_Dens   = xyz_DensBZ

    xza_Basic = xyzf_QMixBZ(:,1,:,:)
    xz_Dens   = xyz_DensBZ(:,1,:)

    call HistoryAutoAddVariable(  &
      & varname="CDensFill", & 
      & dims=(/'x','y','z','t'/),     &
      & longname='Filling Negative term of cloud density ',  &
      & units='kg.m-3.s-1',    &
      & xtype='float')
    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtFill', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Filling Negative term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='float')
    end do

  end subroutine FillNegative_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyza_FillNegative_xyza( xyza_Var )
    
    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,    only: DP
    use gridset,     only: imin,  &! $BG[Ns$N(B X $BJ}8~$N2<8B(B
      &                    imax,  &! $BG[Ns$N(B X $BJ}8~$N>e8B(B
      &                    jmin,  &! $BG[Ns$N(B Y $BJ}8~$N2<8B(B
      &                    jmax,  &! $BG[Ns$N(B Y $BJ}8~$N>e8B(B
      &                    kmin,  &! $BG[Ns$N(B Z $BJ}8~$N2<8B(B
      &                    kmax,  &! $BG[Ns$N(B Z $BJ}8~$N>e8B(B
      &                    nx,ny,nz,ncmax

    implicit none

    real(DP), intent(in) :: xyza_Var(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_FillNegative_xyza(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_DQFILL(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_QSUMPN(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xza_Var(imin:imax, kmin:kmax, ncmax)
    real(DP) :: xza_DQFILL(imin:imax, kmin:kmax, ncmax)
    real(DP) :: xza_QSUMPN(imin:imax, kmin:kmax, ncmax)
    real(DP), parameter  :: EPS = 1.0d-60  !$BNm$G$N3d$j;;$rKI$0(B
    integer              :: i, j, k, s

    
    if (jmin==jmax) then 

      !$B=i4|2=(B
      xza_QSUMPN = 0.0d0
      xza_DQFILL = 0.0d0
      xza_Var(:,:,:) = xyza_Var(:,1,:,:)
      
      do s = 1, ncmax
        do k = kmin+2, kmax-2
          do i = imin+2, imax-2
            xza_QSUMPN(i,k,s) = 1.0d0 /                                     &
              & ( (    MAX( 0.0d0, xza_Basic(i-1,k,s) + xza_Var(i-1,k,s) )  &
              &         * xz_Dens(i-1,k)                                    &
              &      + MAX( 0.0d0, xza_Basic(i+1,k,s) + xza_Var(i+1,k,s) )  &
              &         * xz_Dens(i+1,k)                                    &
              &      + MAX( 0.0d0, xza_Basic(i,k-1,s) + xza_Var(i,k-1,s) )  &
              &         * xz_Dens(i,k-1)                                    &
              &      + MAX( 0.0d0, xza_Basic(i,k+1,s) + xza_Var(i,k+1,s) )  &
              &         * xz_Dens(i,k+1)                                    &
              &     ) * 0.75d0                                              &
              &  + (   MAX( 0.0d0, xza_Basic(i-2,k,s) + xza_Var(i-2,k,s) )  &
              &         * xz_Dens(i-2,k)                                    &
              &      + MAX( 0.0d0, xza_Basic(i+2,k,s) + xza_Var(i+2,k,s) )  &
              &         * xz_Dens(i+2,k)                                    &
              &      + MAX( 0.0d0, xza_Basic(i,k-2,s) + xza_Var(i,k-2,s) )  &
              &         * xz_Dens(i,k-2)                                    &
              &      + MAX( 0.0d0, xza_Basic(i,k+2,s) + xza_Var(i,k+2,s) )  &
              &         * xz_Dens(i,k+2)                                    &
              &     ) * 0.25d0                                              &
              &  + EPS ) 
          end do
        end do
      end do
      
      do s = 1, ncmax
        do k = 1, nz
          do i = 1, nx
            xza_DQFILL(i,k,s) =                                               &
              &  - MIN( 0.0d0, xza_Basic(i,k,s) + xza_Var(i,k,s) )            &
              &  + MAX( 0.0d0, xza_Basic(i,k,s) + xza_Var(i,k,s) )            &
              &    * ( ( MIN( 0.0d0, xza_Basic(i-1,k,s) + xza_Var(i-1,k,s) )  &
              &           * xza_QSUMPN(i-1,k,s)                               &
              &           * xz_Dens(i-1,k)                                    &
              &        + MIN( 0.0d0, xza_Basic(i+1,k,s) + xza_Var(i+1,k,s) )  &
              &           * xza_QSUMPN(i+1,k,s)                               &
              &           * xz_Dens(i+1,k)                                    &
              &        + MIN( 0.0d0, xza_Basic(i,k-1,s) + xza_Var(i,k-1,s) )  &
              &           * xza_QSUMPN(i,k-1,s)                               &
              &           * xz_Dens(i,k-1)                                    &
              &        + MIN( 0.0d0, xza_Basic(i,k+1,s) + xza_Var(i,k+1,s) )  &
              &           * xza_QSUMPN(i,k+1,s)                               &
              &           * xz_Dens(i,k+1)                                    &  
              &       ) * 0.75d0                                              &
              &      + ( MIN( 0.0d0, xza_Basic(i-2,k,s) + xza_Var(i-2,k,s) )  &
              &           * xza_QSUMPN(i-2,k,s)                               &
              &           * xz_Dens(i-2,k)                                    &
              &        + MIN( 0.0d0, xza_Basic(i+2,k,s) + xza_Var(i+2,k,s) )  &
              &           * xza_QSUMPN(i+2,k,s)                               &
              &           * xz_Dens(i+2,k)                                    &
              &        + MIN( 0.0d0, xza_Basic(i,k-2,s) + xza_Var(i,k-2,s) )  &
              &           * xza_QSUMPN(i,k-2,s)                               &
              &           * xz_Dens(i,k-2)                                    &
              &        + MIN( 0.0d0, xza_Basic(i,k+2,s) + xza_Var(i,k+2,s) )  &
              &           * xza_QSUMPN(i,k+2,s)                               &
              &           * xz_Dens(i,k+2)                                    &
              &       ) * 0.25d0                                              &
              &     )
          end do
        end do
      end do
      
      !$B=PNO(B
      xyza_FillNegative_xyza(:,1,:,:) = xza_Var + xza_DQFILL
      
    else
      
      !$B=i4|2=(B
      xyza_QSUMPN = 0.0d0
      xyza_DQFILL = 0.0d0
      
      do s = 1, ncmax
        do k = kmin+2, kmax-2
          do j = jmin+2, jmax-2
            do i = imin+2, imax-2
              xyza_QSUMPN(i,j,k,s) = 1.0d0 /                                       &
                & ( (    MAX( 0.0d0, xyza_Basic(i-1,j,k,s) + xyza_Var(i-1,j,k,s) ) &
                &         * xyz_Dens(i-1,j,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i+1,j,k,s) + xyza_Var(i+1,j,k,s) ) &
                &         * xyz_Dens(i+1,j,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j-1,k,s) + xyza_Var(i,j-1,k,s) ) &
                &         * xyz_Dens(i,j-1,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j+1,k,s) + xyza_Var(i,j+1,k,s) ) &
                &         * xyz_Dens(i,j+1,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j,k-1,s) + xyza_Var(i,j,k-1,s) ) &
                &         * xyz_Dens(i,j,k-1)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j,k+1,s) + xyza_Var(i,j,k+1,s) ) &
                &         * xyz_Dens(i,j,k+1)                                      &
                &     ) * 0.75d0                                                   &
                &  + (   MAX( 0.0d0, xyza_Basic(i-2,j,k,s) + xyza_Var(i-2,j,k,s) ) &
                &         * xyz_Dens(i-2,j,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i+2,j,k,s) + xyza_Var(i+2,j,k,s) ) &
                &         * xyz_Dens(i+2,j,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j-2,k,s) + xyza_Var(i,j-2,k,s) ) &
                &         * xyz_Dens(i,j-2,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j+2,k,s) + xyza_Var(i,j+2,k,s) ) &
                &         * xyz_Dens(i,j+2,k)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j,k-2,s) + xyza_Var(i,j,k-2,s) ) &
                &         * xyz_Dens(i,j,k-2)                                      &
                &      + MAX( 0.0d0, xyza_Basic(i,j,k+2,s) + xyza_Var(i,j,k+2,s) ) &
                &         * xyz_Dens(i,j,k+2)                                      &
                &     ) * 0.25d0                                              &
                &  + EPS ) 
            end do
          end do
        end do
      end do
      
      do s = 1, ncmax
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              xyza_DQFILL(i,j,k,s) =                                                 &
                &  - MIN( 0.0d0, xyza_Basic(i,j,k,s) + xyza_Var(i,j,k,s) )           &
                &  + MAX( 0.0d0, xyza_Basic(i,j,k,s) + xyza_Var(i,j,k,s) )           &
                &    * ( ( MIN( 0.0d0, xyza_Basic(i-1,j,k,s) + xyza_Var(i-1,j,k,s) ) &
                &           * xyza_QSUMPN(i-1,j,k,s) * xyz_Dens(i-1,j,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i+1,j,k,s) + xyza_Var(i+1,j,k,s) ) &
                &           * xyza_QSUMPN(i+1,j,k,s) * xyz_Dens(i+1,j,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j-1,k,s) + xyza_Var(i,j-1,k,s) ) &
                &           * xyza_QSUMPN(i,j-1,k,s) * xyz_Dens(i,j-1,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j+1,k,s) + xyza_Var(i,j+1,k,s) ) &
                &           * xyza_QSUMPN(i,j+1,k,s) * xyz_Dens(i,j+1,k)             &  
                &        + MIN( 0.0d0, xyza_Basic(i,j,k-1,s) + xyza_Var(i,j,k-1,s) ) &
                &           * xyza_QSUMPN(i,j,k-1,s) * xyz_Dens(i,j,k-1)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j,k+1,s) + xyza_Var(i,j,k+1,s) ) &
                &           * xyza_QSUMPN(i,j,k+1,s) * xyz_Dens(i,j,k+1)             &  
                &       ) * 0.75d0                                                   &
                &      + ( MIN( 0.0d0, xyza_Basic(i-2,j,k,s) + xyza_Var(i-2,j,k,s) ) &
                &           * xyza_QSUMPN(i-2,j,k,s) * xyz_Dens(i-2,j,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i+2,j,k,s) + xyza_Var(i+2,j,k,s) ) &
                &           * xyza_QSUMPN(i+2,j,k,s) * xyz_Dens(i+2,j,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j-2,k,s) + xyza_Var(i,j-2,k,s) ) &
                &           * xyza_QSUMPN(i,j-2,k,s) * xyz_Dens(i,j-2,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j+2,k,s) + xyza_Var(i,j+2,k,s) ) &
                &           * xyza_QSUMPN(i,j+2,k,s) * xyz_Dens(i,j+2,k)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j,k-2,s) + xyza_Var(i,j,k-2,s) ) &
                &           * xyza_QSUMPN(i,j,k-2,s) * xyz_Dens(i,j,k-2)             &
                &        + MIN( 0.0d0, xyza_Basic(i,j,k+2,s) + xyza_Var(i,j,k+2,s) ) &
                &           * xyza_QSUMPN(i,j,k+2,s) * xyz_Dens(i,j,k+2)             &
                &       ) * 0.25d0                                                   &
                &     )
            end do
          end do
        end do
      end do

      !$B=PNO(B
      xyza_FillNegative_xyza = xyza_Var + xyza_DQFILL
      
    end if

  end function xyza_FillNegative_xyza

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyz_FillNegativeDensity_xyz( xyz_Var )
    ! $B<g@.J,6E7kBPN.MQ(B. 
    ! $B1@L)EY$KBP$7$F$N$_E,MQ2DG=(B. 
    ! $B9T$J$C$F$$$k$3$H$OLZ@1HG$HF1$8(B. 

    !$B%b%8%e!<%kFI$_9~$_(B
    use dc_types,    only: DP
    use gridset,     only: imin,  &! $BG[Ns$N(B X $BJ}8~$N2<8B(B
      &                    imax,  &! $BG[Ns$N(B X $BJ}8~$N>e8B(B
      &                    jmin,  &! $BG[Ns$N(B Y $BJ}8~$N2<8B(B
      &                    jmax,  &! $BG[Ns$N(B Y $BJ}8~$N>e8B(B
      &                    kmin,  &! $BG[Ns$N(B Z $BJ}8~$N2<8B(B
      &                    kmax    ! $BG[Ns$N(B Z $BJ}8~$N>e8B(B

    implicit none

    real(DP), intent(inout) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_FillNegativeDensity_xyz(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_DQFILL(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_QSUMPN(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xz_Var(imin:imax, kmin:kmax)
    real(DP)                :: xz_DQFILL(imin:imax, kmin:kmax)
    real(DP)                :: xz_QSUMPN(imin:imax, kmin:kmax)
    real(DP), parameter     :: EPS = 1.0d-60  !$BNm$G$N3d$j;;$rKI$0(B
    integer                 :: i, j, k

    if (jmin==jmax) then     

      !$B=i4|2=(B
!      xz_QSUMPN = 0.0d0
!      xz_DQFILL = 0.0d0
      xz_Var = xyz_Var(:,1,:)

      do k = kmin+2, kmax-2
        do i = imin+2, imax-2
          xz_QSUMPN(i,k) = 1.0d0 /                &
            & ( (    MAX( 0.0d0, xz_Var(i-1,k) )  &
            &      + MAX( 0.0d0, xz_Var(i+1,k) )  &
            &      + MAX( 0.0d0, xz_Var(i,k-1) )  &
            &      + MAX( 0.0d0, xz_Var(i,k+1) )  &
            &     ) * 0.75d0                      &
            &  + (   MAX( 0.0d0, xz_Var(i-2,k) )  &
            &      + MAX( 0.0d0, xz_Var(i+2,k) )  &
            &      + MAX( 0.0d0, xz_Var(i,k-2) )  &
            &      + MAX( 0.0d0, xz_Var(i,k+2) )  &
            &     ) * 0.25d0                      &
            &  + EPS ) 
        end do
      end do

      do k = kmin+2, kmax-2
        do i = imin+2, imax-2
          xz_DQFILL(i,k) =                            &
            &  - MIN( 0.0d0, xz_Var(i,k) )            &
            &  + MAX( 0.0d0, xz_Var(i,k) )            &
            &    * ( ( MIN( 0.0d0, xz_Var(i-1,k) ) * xz_QSUMPN(i-1,k) &
            &        + MIN( 0.0d0, xz_Var(i+1,k) ) * xz_QSUMPN(i+1,k) &
            &        + MIN( 0.0d0, xz_Var(i,k-1) ) * xz_QSUMPN(i,k-1) &
            &        + MIN( 0.0d0, xz_Var(i,k+1) ) * xz_QSUMPN(i,k+1) &
            &       ) * 0.75d0                                        &
            &      + ( MIN( 0.0d0, xz_Var(i-2,k) ) * xz_QSUMPN(i-2,k) &
            &        + MIN( 0.0d0, xz_Var(i+2,k) ) * xz_QSUMPN(i+2,k) &
            &        + MIN( 0.0d0, xz_Var(i,k-2) ) * xz_QSUMPN(i,k-2) &
            &        + MIN( 0.0d0, xz_Var(i,k+2) ) * xz_QSUMPN(i,k+2) &
            &       ) * 0.25d0                                        &
            &     )
        end do
      end do
      !$B=PNO(B
      xyz_FillNegativeDensity_xyz(:,1,:) = xz_Var + xz_DQFILL

    else

      !$B=i4|2=(B
!      xyz_QSUMPN = 0.0d0
!      xyz_DQFILL = 0.0d0
      
      do k = kmin+2, kmax-2
        do j = jmin+2, jmax-2
          do i = imin+2, imax-2
            xyz_QSUMPN(i,j,k) = 1.0d0 /               &
              & ( (    MAX( 0.0d0, xyz_Var(i-1,j,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i+1,j,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j-1,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j+1,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j,k-1) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j,k+1) ) &
              &     ) * 0.75d0                        &
              &  + (   MAX( 0.0d0, xyz_Var(i-2,j,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i+2,j,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j-2,k ) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j+2,k) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j,k-2) ) &
              &      + MAX( 0.0d0, xyz_Var(i,j,k+2) ) &
              &     ) * 0.25d0                        &
              &  + EPS ) 
          end do
        end do
      end do

      do k = kmin+2, kmax-2
        do j = jmin+2, jmax-2
          do i = imin+2, imax-2
            xyz_DQFILL(i,j,k) =                                                &
              &  - MIN( 0.0d0, xyz_Var(i,j,k) )                                 &
              &  + MAX( 0.0d0, xyz_Var(i,j,k) )                                 &
              &    * ( ( MIN( 0.0d0, xyz_Var(i-1,j,k) ) * xyz_QSUMPN(i-1,j,k)  &
              &        + MIN( 0.0d0, xyz_Var(i+1,j,k) ) * xyz_QSUMPN(i+1,j,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j-1,k) ) * xyz_QSUMPN(i,j-1,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j+1,k) ) * xyz_QSUMPN(i,j+1,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j,k-1) ) * xyz_QSUMPN(i,j,k-1)  &
              &        + MIN( 0.0d0, xyz_Var(i,j,k+1) ) * xyz_QSUMPN(i,j,k+1)  &
              &       ) * 0.75d0                                                &
              &      + ( MIN( 0.0d0, xyz_Var(i-2,j,k) ) * xyz_QSUMPN(i-2,j,k)  &
              &        + MIN( 0.0d0, xyz_Var(i+2,j,k) ) * xyz_QSUMPN(i+2,j,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j-2,k) ) * xyz_QSUMPN(i,j-2,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j+2,k) ) * xyz_QSUMPN(i,j+2,k)  &
              &        + MIN( 0.0d0, xyz_Var(i,j,k-2) ) * xyz_QSUMPN(i,j,k-2)  &
              &        + MIN( 0.0d0, xyz_Var(i,j,k+2) ) * xyz_QSUMPN(i,j,k+2)  &
              &       ) * 0.25d0                                                &
              &     )
          end do
        end do
      end do
      
      !$B=PNO(B
      xyz_FillNegativeDensity_xyz = xyz_Var + xyz_DQFILL

    end if
    
  end function xyz_FillNegativeDensity_xyz
    
end module FillNegative
