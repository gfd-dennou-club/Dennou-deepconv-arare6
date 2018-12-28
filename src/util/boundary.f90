!= 「そで領域」に値を代入するためのモジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: boundary.f90,v 1.7 2009-03-10 16:17:25 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module boundary 
  !
  ! 有限差分モデル用 境界条件モジュール
  ! 並列化には対応していない (テスト用として残す). 
  !

  !暗黙の型宣言禁止
  implicit none

  !private 属性にする
  private

  !関数を public にする
  public BoundaryXCyc_xza
  public BoundaryXCyc_xz
  public BoundaryXCyc_pz
  public BoundaryXCyc_xr

  public BoundaryZCyc_xz
  public BoundaryZCyc_pz
  public BoundaryZCyc_xr

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

  interface BoundaryZCyc_xz
    module procedure BoundaryZCyc_aa
  end interface
  interface BoundaryZCyc_pz
    module procedure BoundaryZCyc_aa
  end interface 
  interface BoundaryZCyc_xr
    module procedure BoundaryZCyc_aa
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

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryXCyc_aa( aa_Var )
    !
    ! x 方向に「周期境界条件」を適用する. 
    ! 格子点, 半格子点においても, 関数の形式は同じ
    !
    
    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginX,                & !x 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegXMin, RegXMax          !物理領域のサイズ

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: aa_Var(imin:imax, kmin:kmax)
    real(DP)                 :: aa_Work(imin:imax, kmin:kmax)
    integer                 :: i
    
    aa_Work = aa_Var

    !壁の値を決める
    aa_Var(RegXMin, :) = aa_Work(RegXMax, :)    

    !壁の外側の値を決める
    do i = 1, MarginX
      aa_Var(RegXMin - i, :) = aa_Work(RegXMax - i, :)    
      aa_Var(RegXMax + i, :) = aa_Work(RegXMin + i, :) 
    end do
    
  end subroutine BoundaryXCyc_aa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZCyc_aa( aa_Var )
    !
    ! z 方向に「周期境界条件」を適用する. 
    ! 格子点, 半格子点においても, 関数の形式は同じ
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax          !物理領域のサイズ

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: aa_Var(imin:imax, kmin:kmax)
    real(DP)                 :: aa_Work(imin:imax, kmin:kmax)
    integer                 :: i
    
    aa_Work = aa_Var

    !壁の値を決める
    aa_Var(:,RegZMin) = aa_Work(:,RegZMax)    

    !壁の外側の値を決める
    do i = 1, MarginZ
      aa_Var(:,RegZMin - i) = aa_Work(:,RegZMax - i)    
      aa_Var(:,RegZMax + i) = aa_Work(:,RegZMin + i) 
    end do
    
  end subroutine BoundaryZCyc_aa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryXCyc_aaa( aaa_Var )
    !
    ! x 方向に「周期境界条件」を適用する. 
    ! 格子点, 半格子点においても, 関数の形式は同じ
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginX,                & !x 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegXMin, RegXMax,       & !x 方向の物理領域の上限
      &                 SpcNum                    !化学種の数

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: aaa_Var(imin:imax, kmin:kmax, SpcNum)
    real(DP)                 :: aaa_Work(imin:imax, kmin:kmax, SpcNum)
    integer                 :: i
    
    aaa_Work = aaa_Var

    !壁の値を決める
    aaa_Var(RegXMin, :, :) = aaa_Work(RegXMax, :, :)

    !壁の外側の値を決める
    do i = 1, MarginX
      aaa_Var(RegXMin - i, :, :) = aaa_Work(RegXMax - i, :, :)
      aaa_Var(RegXMax + i, :, :) = aaa_Work(RegXMin + i, :, :) 
    end do
    
  end subroutine BoundaryXCyc_aaa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZSym_az( az_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax          !物理領域のサイズ

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: az_Var(imin:imax, kmin:kmax)
    real(DP)                 :: az_Work(imin:imax, kmin:kmax)
    integer                 :: k
  
    az_Work = az_Var

    do k = 0, MarginZ
      az_Var( :, RegZMin - k ) = az_Work( :, RegZMin + 1 + k )
    end do
    do k = 1, MarginZ
      az_Var( :, RegZMax + k ) = az_Work( :, RegZMax + 1 - k )
    end do
    
  end subroutine BoundaryZSym_az

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZSym_aza( aza_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset, only: MarginZ,                 & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax,       & !物理領域のサイズ
      &                 SpcNum                    !化学種の数

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: aza_Var(imin:imax, kmin:kmax, SpcNum)
    real(DP)                 :: aza_Work(imin:imax, kmin:kmax, SpcNum)
    integer                 :: k
  
    aza_Work = aza_Var

    do k = 0, MarginZ
      aza_Var( :, RegZMin - k, : ) = aza_Work( :, RegZMin + 1 + k, : )
    end do
    do k = 1, MarginZ
      aza_Var( :, RegZMax + k, : ) = aza_Work( :, RegZMax + 1 - k, : )
    end do
    
  end subroutine BoundaryZSym_aza

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZAntiSym_az( az_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax          !z 方向の物理領域の上限

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: az_Var(imin:imax, kmin:kmax)
    real(DP)                 :: az_Work(imin:imax, kmin:kmax)
    integer                 :: k
    
    az_Work = az_Var
    
    do k = 0, MarginZ
      az_Var( :, RegZMin - k ) = - az_Work( :, RegZMin + 1 + k )
    end do
    do k = 1, MarginZ
      az_Var( :, RegZMax + k ) = - az_Work( :, RegZMax + 1 - k )
    end do
    
  end subroutine BoundaryZAntiSym_az
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZAntiSym_aza( aza_Var )
    !
    ! z 方向に半格子ずれた点に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax,       & !物理領域のサイズ
      &                 SpcNum                    !化学種の数

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: aza_Var(imin:imax, kmin:kmax, SpcNum)
    real(DP)                 :: aza_Work(imin:imax, kmin:kmax, SpcNum)
    integer                 :: k
    
    aza_Work = aza_Var
    
    do k = 0, MarginZ
      aza_Var( :, RegZMin - k, : ) = - aza_Work( :, RegZMin + 1 + k, : )
    end do
    do k = 1, MarginZ
      aza_Var( :, RegZMax + k, : ) = - aza_Work( :, RegZMax + 1 - k, : )
    end do
    
  end subroutine BoundaryZAntiSym_aza
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZSym_ar( ar_Var )
    !
    ! z 方向の格子点上に存在する変数に対し, 
    ! z 方向に「対称境界条件」を適用する. 
    !
    
    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !x 方向の配列の下限
      &                 RegZMin, RegZMax          !物理領域のサイズ

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: ar_Var(imin:imax, kmin:kmax)
    real(DP)                 :: ar_Work(imin:imax, kmin:kmax)
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

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine BoundaryZAntiSym_ar( ar_Var )
    !
    ! z 方向の格子点上に存在する変数に対し, 
    ! z 方向に「反対称境界条件」を適用する. 
    !
    
    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: MarginZ,                & !z 方向の境界のグリッド数
      &                 imin, imax, kmin, kmax, & !配列サイズ
      &                 RegZMin, RegZMax          !物理領域のサイズ

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout)  :: ar_Var(imin:imax, kmin:kmax)
    real(DP)                 :: ar_Work(imin:imax, kmin:kmax)
    integer                  :: k
  
    ar_Work = ar_Var

    !境界での速度はゼロ
    ar_Var( :, RegZMin ) = 0.0d0
    ar_Var( :, RegZMax ) = 0.0d0

    do k = 1, MarginZ
      ar_Var( :, RegZMin - k ) = - ar_Work( :, RegZMin + k )
      ar_Var( :, RegZMax + k ) = - ar_Work( :, RegZMax - k )
    end do
    
  end subroutine BoundaryZAntiSym_ar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module boundary
