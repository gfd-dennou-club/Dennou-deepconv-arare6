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
  !負の雲水量などを穴埋めするためのルーチン
  !  * 負となった場合には, 周囲の 8 格子点の値を削って穴埋めする.
  !  * F77版 deepconv の QFILL.f を改変 

  !モジュール読み込み
  use dc_types,   only: DP

  !暗黙の型宣言禁止
  implicit none

  !公開要素
  public FillNegative_Init
  public xyza_FillNegative_xyza
  public xyz_FillNegativeDensity_xyz

  !変数定義
  real(DP), save, private, allocatable  :: xyza_Basic(:,:,:,:)
  real(DP), save, private, allocatable  :: xyz_Dens(:,:,:)
  real(DP), save, private, allocatable  :: xza_Basic(:,:,:)
  real(DP), save, private, allocatable  :: xz_Dens(:,:)

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine FillNegative_Init( )    
    !
    ! 基準値を取得. 
    ! 基本場と擾乱成分の和がゼロよりも小さくなることは許容しない
    !
    
    !モジュール読み込み
    use gtool_historyauto, &
      &              only: HistoryAutoAddVariable
    use gridset,     only: imin,  &! 配列の X 方向の下限
      &                    imax,  &! 配列の X 方向の上限
      &                    jmin,  &! 配列の Y 方向の下限
      &                    jmax,  &! 配列の Y 方向の上限
      &                    kmin,  &! 配列の Z 方向の下限
      &                    kmax,  &! 配列の Z 方向の上限
      &                    ncmax
    use basicset,    only: xyzf_QMixBZ, xyz_DensBZ
    use composition, only: SpcWetSymbol

    !暗黙の型宣言禁止
    implicit none

    !入力変数
    integer  :: l

    !初期化
    allocate( xyza_Basic(imin:imax,jmin:jmax,kmin:kmax, ncmax) )
    allocate( xyz_Dens  (imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xza_Basic(imin:imax,kmin:kmax, ncmax) )
    allocate( xz_Dens(imin:imax,kmin:kmax ) )
    
    !値の代入
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
    
    !モジュール読み込み
    use dc_types,    only: DP
    use gridset,     only: imin,  &! 配列の X 方向の下限
      &                    imax,  &! 配列の X 方向の上限
      &                    jmin,  &! 配列の Y 方向の下限
      &                    jmax,  &! 配列の Y 方向の上限
      &                    kmin,  &! 配列の Z 方向の下限
      &                    kmax,  &! 配列の Z 方向の上限
      &                    nx,ny,nz,ncmax

    implicit none

    real(DP), intent(in) :: xyza_Var(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_FillNegative_xyza(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_DQFILL(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyza_QSUMPN(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xza_Var(imin:imax, kmin:kmax, ncmax)
    real(DP) :: xza_DQFILL(imin:imax, kmin:kmax, ncmax)
    real(DP) :: xza_QSUMPN(imin:imax, kmin:kmax, ncmax)
    real(DP), parameter  :: EPS = 1.0d-60  !零での割り算を防ぐ
    integer              :: i, j, k, s

    
    if (jmin==jmax) then 

      !初期化
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
      
      !出力
      xyza_FillNegative_xyza(:,1,:,:) = xza_Var + xza_DQFILL
      
    else
      
      !初期化
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

      !出力
      xyza_FillNegative_xyza = xyza_Var + xyza_DQFILL
      
    end if

  end function xyza_FillNegative_xyza

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyz_FillNegativeDensity_xyz( xyz_Var )
    ! 主成分凝結対流用. 
    ! 雲密度に対してのみ適用可能. 
    ! 行なっていることは木星版と同じ. 

    !モジュール読み込み
    use dc_types,    only: DP
    use gridset,     only: imin,  &! 配列の X 方向の下限
      &                    imax,  &! 配列の X 方向の上限
      &                    jmin,  &! 配列の Y 方向の下限
      &                    jmax,  &! 配列の Y 方向の上限
      &                    kmin,  &! 配列の Z 方向の下限
      &                    kmax    ! 配列の Z 方向の上限

    implicit none

    real(DP), intent(inout) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_FillNegativeDensity_xyz(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_DQFILL(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_QSUMPN(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xz_Var(imin:imax, kmin:kmax)
    real(DP)                :: xz_DQFILL(imin:imax, kmin:kmax)
    real(DP)                :: xz_QSUMPN(imin:imax, kmin:kmax)
    real(DP), parameter     :: EPS = 1.0d-60  !零での割り算を防ぐ
    integer                 :: i, j, k

    if (jmin==jmax) then     

      !初期化
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
      !出力
      xyz_FillNegativeDensity_xyz(:,1,:) = xz_Var + xz_DQFILL

    else

      !初期化
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
      
      !出力
      xyz_FillNegativeDensity_xyz = xyz_Var + xyz_DQFILL

    end if
    
  end function xyz_FillNegativeDensity_xyz
    
end module FillNegative
