!= Module BasicEnvInit
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_basic.f90,v 1.9 2011/06/17 19:07:25 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
!デフォルトの基本場を設定するための変数参照型モジュール
!   * BasicEnvFile_init: 基本場の値を netCDF ファイルから取得
!   * BasicEnvCalc_Init: 基本場の情報を Namelist から取得して値を計算
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!== Future Plans
!
module initialdata_basic
  !
  ! デフォルトの基本場を作るためのルーチン

    
  !暗黙の型宣言禁止
  implicit none

  !デフォルトは private
  private

  real(8), save  :: Humidity = 0.0d0        !基本場の湿度 
  real(8), save  :: TempStrat = 100.0d0     !成層圏の温度 [k]
  real(8), save  :: HeightStrat = 18.0d3    !成層圏の高度 [k]
  real(8), save  :: Dhight = 5.0d3          !重み関数のパラメータ [m]

  ! public
  public initialdata_basic_nml
  public initialdata_basic_dry
  public initialdata_basic_moist
  public initialdata_basic_strat
  public initialdata_basic_isothermal
    
contains

  subroutine initialdata_basic_dry( z_TempBZ, z_PressBZ, za_MolFr )

    implicit none

    real(8), intent(out) :: z_TempBZ(kmin:kmax)
    real(8), intent(out) :: z_PressBZ(kmin:kmax)
    real(8), intent(out) :: za_MolFr(kmin:kmax, 1:ncmax)
    


  end subroutine initialdata_basic_dry


  subroutine initialdata_basic_moist( z_TempBZ, z_PressBZ, za_MolFr )

    implicit none

    real(8), intent(out) :: z_TempBZ(kmin:kmax)
    real(8), intent(out) :: z_PressBZ(kmin:kmax)
    real(8), intent(out) :: za_MolFr(kmin:kmax, 1:ncmax)
    
    
    
  end subroutine initialdata_basic_moist

  

 
end module initialdata_basic
