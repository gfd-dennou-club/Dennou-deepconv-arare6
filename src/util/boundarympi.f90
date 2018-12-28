!= Module BoundaryMPI
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: boundarympi.f90,v 1.1 2008-09-22 07:33:59 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module BoundaryMPI

  !暗黙の型宣言禁止
  implicit none

  !private 属性にする
  private
	
  !関数を public にする
  public BoundaryXCyc_xza
  public BoundaryXCyc_xz
  public BoundaryXCyc_pz
  public BoundaryXCyc_xr

  public BoundaryZSym_xza
  public BoundaryZSym_xz
  public BoundaryZSym_pz
  public BoundaryZSym_xr

  public BoundaryZAntiSym_xza
  public BoundaryZAntiSym_xz
  public BoundaryZAntiSym_pz
  public BoundaryZAntiSym_xr

  !関数の定義. 2 度 3 度同じような関数を定義するのは避けるため.
  interface BoundaryXCyc_xza
    module procedure BoundaryXCyc_aaa
  end interface
  interface BoundaryXCyc_xz
    module procedure BoundaryXCyc_aa
  end interface
  interface BoundaryXCyc_pz
    module procedure BoundaryXCyc_aa
  end interface 
  interface BoundaryXCyc_xr
    module procedure BoundaryXCyc_aa
  end interface 

  interface BoundaryZSym_xza
    module procedure BoundaryZSym_aza
  end interface 
  interface BoundaryZSym_xz
    module procedure BoundaryZSym_az
  end interface 
  interface BoundaryZSym_pz
    module procedure BoundaryZSym_az
  end interface
  interface BoundaryZSym_xr
    module procedure BoundaryZSym_ar
  end interface 
  
  interface BoundaryZAntiSym_xza
    module procedure BoundaryZAntiSym_aza
  end interface 
  interface BoundaryZAntiSym_xz
    module procedure BoundaryZAntiSym_az
  end interface 
  interface BoundaryZAntiSym_pz
    module procedure BoundaryZAntiSym_az
  end interface
  interface BoundaryZAntiSym_xr
    module procedure BoundaryZAntiSym_ar
  end interface

contains

!!!---------------------------------------------------------------------!!!
  subroutine BoundaryXCyc_aa(aa_Var)
    
    implicit none

    !モジュール読み込み
    use gridset, only: MarginX,       &! x 方向の境界のグリッド数
      &                MarginZ,       &! z 方向の境界のグリッド数
      &                DimXMin,       &! x 方向の配列の下限
      &                DimXMax,       &! x 方向の配列の上限
      &                DimZMin,       &! z 方向の配列の下限
      &                DimZMax,       &! z 方向の配列の上限
      &                RegXMin,       &! x 方向の物理領域の下限
      &                RegXMax,       &! x 方向の物理領域の上限
      &                RegZMin,       &! z 方向の物理領域の下限
      &                RegZMax,       &! z 方向の物理領域の上限
      &                SpcNum          ! 化学種の数
    use mpiset,  only: myrank,        &!
      &                nprocs,        &!
      &                mpii_isend,    &!
      &                mpii_irecv,    &!
      &                mpii_wait
    
    !変数定義
    real(8), intent(inout)  :: aa_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    integer            :: idest_a, idep_a, idest_b, idep_b
    integer            :: im_a, km_a, im_b, km_b
    integer, parameter :: nvars = 1
    real(8)            :: &
      sbuf_a( MarginX+1, RegZMax-RegZMin, nvars ), &
      rbuf_a( MarginX+1, RegZMax-RegZMin, nvars ), &
      sbuf_b( MarginX  , RegZMax-RegZMin, nvars ), &
      rbuf_b( MarginX  , RegZMax-RegZMin, nvars )
    integer            :: ireqs_a, ireqr_a, ireqs_b, ireqr_b
    integer            :: i, k
    
    
    !----- 変数の初期化 -----
    im_a = MarginX+1
    km_a = RegZMax-RegZMin
    
    im_b = MarginX
    km_b = RegZMax-RegZMin
    
    !-------------------------------
    ! 配列の右側を, ノード間で通信する.

    !送信する配列(sbuf_a)を用意する
    do k = RegZMin+1, RegZMax
      do i = 0, MarginX
        sbuf_a( i + 1, k, 1 ) = aa_var(RegXMax - i , k ) 
      end do
    end do
    
    idest_a = mod(( myrank+1 )       , nprocs)   !送信先
    idep_a  = mod(( myrank-1 )+nprocs, nprocs)   !受信元
    call mpii_isend( idest_a, im_a, km_a, nvars, sbuf_a, ireqs_a )
    call mpii_irecv( idep_a , im_a, km_a, nvars, rbuf_a, ireqr_a )
    
    !-------------------------------
    do k = RegZMin+1, RegZMax
      do i = 1, MarginX
        sbuf_b( i, k, 1 ) = aa_var(RegXMin + i , k ) 
      end do
    end do
    
    idest_b = mod(( myrank-1 )+nprocs, nprocs)    !送信先
    idep_b  = mod(( myrank+1 )       , nprocs)    !受信元
    call mpii_isend( idest_b, im_b, km_b, nvars, sbuf_b, ireqs_b )
    call mpii_irecv( idep_b , im_b, km_b, nvars, rbuf_b, ireqr_b )

    !-------------------------------
    call mpii_wait( ireqs_a )
    call mpii_wait( ireqr_a )
  
    do k = RegZMin+1, RegZMax
      do i = 0, MarginX
        aa_var(RegXMin - i , k ) = rbuf_a( i + 1, k, 1 ) 
      end do
    end do

    !-------------------------------
    call mpii_wait( ireqs_b )
    call mpii_wait( ireqr_b )
    
    do k = RegZMin+1, RegZMax
      do i = 1, MarginX
        aa_var(RegXMax + i , k ) = rbuf_b( i, k, 1 ) 
      end do
    end do

  end subroutine BoundaryXCyc_aa


!!!---------------------------------------------------------------------!!!
  subroutine BoundaryXCyc_aaa(aaa_Var)
    
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: aaa_Var(DimXMin:DimXMax, DimZMin:DimZMax, 1:SpcNum)
    
    integer            :: idest_a, idep_a, idest_b, idep_b
    integer            :: im_a, km_a, im_b, km_b
!    integer, parameter :: nvars = SpcNum
    real(8)            :: &
      sbuf_a( MarginX+1, RegZMax-RegZMin, SpcNum ), &
      rbuf_a( MarginX+1, RegZMax-RegZMin, SpcNum ), &
      sbuf_b( MarginX  , RegZMax-RegZMin, SpcNum ), &
      rbuf_b( MarginX  , RegZMax-RegZMin, SpcNum )
    integer            :: ireqs_a, ireqr_a, ireqs_b, ireqr_b
    integer            :: i, k
    
    
    !----- 変数の初期化 -----
    im_a = MarginX+1
    km_a = RegZMax-RegZMin
    
    im_b = MarginX
    km_b = RegZMax-RegZMin
    
    !-------------------------------
    ! 配列の右側を, ノード間で通信する.

    !送信する配列(sbuf_a)を用意する
    do k = RegZMin+1, RegZMax
      do i = 0, MarginX
        sbuf_a( i + 1, k, 1:SpcNum ) = aaa_var(RegXMax - i , k, 1:SpcNum ) 
      end do
    end do
    
    idest_a = mod(( myrank+1 )       , nprocs)   !送信先
    idep_a  = mod(( myrank-1 )+nprocs, nprocs)   !受信元
    call mpii_isend( idest_a, im_a, km_a, SpcNum, sbuf_a, ireqs_a )
    call mpii_irecv( idep_a , im_a, km_a, SpcNum, rbuf_a, ireqr_a )
    
    !-------------------------------
    do k = RegZMin+1, RegZMax
      do i = 1, MarginX
        sbuf_b( i, k, 1:SpcNum ) = aaa_var(RegXMin + i , k, 1:SpcNum ) 
      end do
    end do
    
    idest_b = mod(( myrank-1 )+nprocs, nprocs)    !送信先
    idep_b  = mod(( myrank+1 )       , nprocs)    !受信元
    call mpii_isend( idest_b, im_b, km_b, SpcNum, sbuf_b, ireqs_b )
    call mpii_irecv( idep_b , im_b, km_b, SpcNum, rbuf_b, ireqr_b )

    !-------------------------------
    call mpii_wait( ireqs_a )
    call mpii_wait( ireqr_a )
  
    do k = RegZMin+1, RegZMax
      do i = 0, MarginX
        aaa_var(RegXMin - i , k, 1:SpcNum) = rbuf_a( i + 1, k, 1:SpcNum ) 
      end do
    end do

    !-------------------------------
    call mpii_wait( ireqs_b )
    call mpii_wait( ireqr_b )
    
    do k = RegZMin+1, RegZMax
      do i = 1, MarginX
        aaa_var(RegXMax + i , k, 1:SpcNum ) = rbuf_b( i, k, 1:SpcNum ) 
      end do
    end do

  end subroutine BoundaryXCyc_aaa

!  subroutine BoundaryXCyc_aa( aa_Var )
!    !
!    ! x 方向に「周期境界条件」を適用する. 
!    ! 格子点, 半格子点においても, 関数の形式は同じ
!    !
!
!    !暗黙の型宣言禁止
!    implicit none
!    
!    !変数定義
!    real(8), intent(inout)  :: aa_Var(DimXMin:DimXMax, DimZMin:DimZMax)
!    real(8)                 :: aa_Work(DimXMin:DimXMax, DimZMin:DimZMax)
!    integer                 :: i
!
!    aa_Work = aa_Var
!    
!    do i = 0, MarginX
!      aa_Var(RegXMin - i, :) = aa_Work(RegXMax - i, :)
!    end do
!    do i = 1, MarginX
!      aa_Var(RegXMax + i, :) = aa_Work(RegXMin + i, :) 
!    end do
!    
!    
!  end subroutine BoundaryXCyc_aa


!!!---------------------------------------------------------------------!!!
!  subroutine BoundaryXCyc_aaa( aaa_Var )
!    !
!    ! x 方向に「周期境界条件」を適用する. 
!    ! 格子点, 半格子点においても, 関数の形式は同じ
!    !
!
!    !暗黙の型宣言禁止
!    implicit none
!    
!    !変数定義
!    real(8), intent(inout)  :: aaa_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
!    real(8)                 :: aaa_Work(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
!    integer                 :: i, s
!    
!    
!    do s = 1, SpcNum
!      call BoundaryXCyc_aa(aaa_Var(:,:,s))
!    end do
!!
!!!    aaa_Work = aaa_Var
!!
!!    do i = 0, MarginX
!!      aaa_Var(RegXMin - i, :, :) = aaa_Work(RegXMax - i, :, :)
!!    end do
!!    do i = 1, MarginX
!!      aaa_Var(RegXMax + i, :, :) = aaa_Work(RegXMin + i, :, :) 
!!    end do
!    
!  end subroutine BoundaryXCyc_aaa

!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZSym_az( az_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: az_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8)                 :: az_Work(DimXMin:DimXMax, DimZMin:DimZMax)
    integer                 :: k
  
    az_Work = az_Var

    do k = 0, MarginZ
      az_Var( :, RegZMin - k ) = az_Work( :, RegZMin + 1 + k )
    end do
    do k = 1, MarginZ
      az_Var( :, RegZMax + k ) = az_Work( :, RegZMax + 1 - k )
    end do
    
  end subroutine BoundaryZSym_az


!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZSym_aza( aza_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: aza_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8)                 :: aza_Work(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    integer                 :: k
  
    aza_Work = aza_Var

    do k = 0, MarginZ
      aza_Var( :, RegZMin - k, : ) = aza_Work( :, RegZMin + 1 + k, : )
    end do
    do k = 1, MarginZ
      aza_Var( :, RegZMax + k, : ) = aza_Work( :, RegZMax + 1 - k, : )
    end do
    
  end subroutine BoundaryZSym_aza


!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZAntiSym_az( az_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: az_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8)                 :: az_Work(DimXMin:DimXMax, DimZMin:DimZMax)
    integer                 :: k
    
    az_Work = az_Var
    
    do k = 0, MarginZ
      az_Var( :, RegZMin - k ) = - az_Work( :, RegZMin + 1 + k )
    end do
    do k = 1, MarginZ
      az_Var( :, RegZMax + k ) = - az_Work( :, RegZMax + 1 - k )
    end do
    
  end subroutine BoundaryZAntiSym_az
  

!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZAntiSym_aza( aza_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: aza_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8)                 :: aza_Work(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    integer                 :: k
    
    aza_Work = aza_Var
    
    do k = 0, MarginZ
      aza_Var( :, RegZMin - k, : ) = - aza_Work( :, RegZMin + 1 + k, : )
    end do
    do k = 1, MarginZ
      aza_Var( :, RegZMax + k, : ) = - aza_Work( :, RegZMax + 1 - k, : )
    end do
    
  end subroutine BoundaryZAntiSym_aza
  

!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZSym_ar( ar_Var )
    !
    ! z 方向の格子点上に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: ar_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8)                 :: ar_Work(DimXMin:DimXMax, DimZMin:DimZMax)
    integer                 :: k
  
    ar_Work = ar_Var

    !境界での速度はゼロ
    ar_Var( :, RegZMin ) = 0.0d0
    ar_Var( :, RegZMax ) = 0.0d0
    
    do k = 1, MarginZ
      ar_Var( :, RegZMin - k ) = ar_Work( :, RegZMin + k )
      ar_Var( :, RegZMax + k ) = ar_Work( :, RegZMax - k )
    end do
    
  end subroutine BoundaryZSym_ar


!!!---------------------------------------------------------------------!!!
  subroutine BoundaryZAntiSym_ar( ar_Var )
    !
    ! z 方向の格子点上に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(8), intent(inout)  :: ar_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8)                 :: ar_Work(DimXMin:DimXMax, DimZMin:DimZMax)
    integer                 :: k
  
    ar_Work = ar_Var

    !境界での速度はゼロ
    ar_Var( :, RegZMin ) = 0.0d0
    ar_Var( :, RegZMax ) = 0.0d0

    do k = 1, MarginZ
      ar_Var( :, RegZMin - k ) = - ar_Work( :, RegZMin + k )
      ar_Var( :, RegZMax + k ) = - ar_Work( :, RegZMax - k )
    end do
    
  end subroutine BoundaryZAntiSym_ar


end module BoundaryMPI
