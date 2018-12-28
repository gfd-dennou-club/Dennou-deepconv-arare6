!= Module initialdata_Toon2002
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_toon2002.f90,v 1.3 2014/03/04 04:49:40 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module initialdata_Toon2002
  !
  ! Toon et al. (2002) を模した初期場を与える
  !

  !モジュール読み込み
  use dc_types,  only: DP

  !暗黙の型宣言禁止
  implicit none

  !内部変数
  real(DP), parameter, private :: AntA    = 27.4d0 
  real(DP), parameter, private :: AntB    = 3103.0d0
  real(DP), parameter, private :: TempLTP = 135.0d0
  real(DP), parameter, private :: Dhight  = 6.0d2

  !初期化だけ公開
  public  initialdata_Toon2002_basic

contains

!!!------------------------------------------------------------------------------!!!

  subroutine initialdata_toon2002_basic( z_Temp, z_Press )
    
    !モジュール読み込み
    use dc_types,  only: DP
    use gridset,   only: kmin, kmax,  &!配列サイズ (Z 方向)
      &                  nz            !物理領域のサイズ (Z 方向)
    use axesset,   only: z_Z,         &!スカラー格子点での高度
      &                  z_dz          !Z 方向の格子点間隔
    use constants, only: GasRDry,     &!乾燥成分の定圧比熱
      &                  CpDry,       &!乾燥成分の定圧比熱
      &                  Grav,        &!重力加速度
      &                  TempSfc,     &!地表面温度
      &                  PressSfc      !地表面圧力

    implicit none
    
    real(DP), intent(out):: z_Press(kmin:kmax)           !圧力
    real(DP), intent(out):: z_Temp(kmin:kmax)            !温度
    real(DP)             :: TempLCL
    real(DP)             :: Temp_0,  Temp_1
    real(DP)             :: Press_0, Press_1
    real(DP)             :: Weight1, Weight2
    real(DP)             :: LCL, LTP
    integer              :: k
    

    ! 乾燥断熱線, 湿潤断熱線, 等温線が交わる高度を計算し,
    ! 各領域で成り立つ式を用いて温度, 圧力を計算
    ! 乾燥断熱線と湿潤断熱線が交わる高度(LCL)を反復法で計算
    !
    Press_0 = PressSfc
    Temp_0 = TempSfc
    do
      ! 飽和温度 (press0 に対する): 最初は地表面での飽和温度.
      ! ln(p) = A - B/T
      Temp_1 = AntB / (AntA - dlog(Press_0))
      
      ! 乾燥断熱的に決めた圧力 
      !
      Press_1 = PressSfc * (Temp_1/TempSfc) **(CpDry / GasRDry)

      ! 絶対誤差が閾値より小さくなれば終了. 
      !
      if (abs(Temp_1 - Temp_0) < epsilon(0.0d0)) then
        LCL = (TempSfc * CpDry) / Grav &
          & * (1.0d0 - (Press_1 / PressSfc)**(GasRDry / CpDry))
        TempLCL = temp_1

        exit
      else
        Temp_0 = Temp_1
        Press_0 = Press_1
      end if
    end do


    ! 湿潤断熱線と等温線が交わる高度(LTP)を計算
    !
    LTP = LCL + GasRDry * AntB / Grav * dlog(TempLCL / TempLTP)
  
    ! 温度圧力を決める.
    !
    z_Temp(1)  = TempSfc  - Grav * z_Z(1) / CpDry 
    z_Press(1) = PressSfc - (Grav * PressSfc * z_dz(1) * 5.0d-1) / (GasRDry * TempSfc)
    do k = 2, nz

      !重みつけの関数を用意. tanh を用いる
      Weight1 = ( tanh( (z_Z(k) - LCL ) / Dhight ) + 1.0d0 ) * 5.0d-1

      !重みつけの関数を用意. tanh を用いる
      Weight2 = ( tanh( (z_Z(k) - LTP ) / Dhight ) + 1.0d0 ) * 5.0d-1

      !乾燥断熱
      if (z_z(k) < LCL) then 
        z_Temp(k) = TempSfc - Grav * z_Z(k) / CpDry 

      !湿潤断熱
      elseif (z_z(k) >= LCL .AND. z_z(k) < LTP) then 
        Temp_0 = TempSfc - Grav * z_Z(k) / CpDry 
        Temp_1 = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))
        z_Temp(k) = Temp_0 * ( 1.0d0 - Weight1 ) + Temp_1 * Weight1
!        z_Temp(k)  = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))

      !等温
      elseif (z_z(k) >= LTP) then 
        Temp_0 = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))
        z_Temp(k) = Temp_0 * ( 1.0d0 - Weight2 ) + TempLTP * Weight2
!        z_Temp(k) = TempLTP

      end if
    end do

    ! 静水圧平衡から圧力を決める
    !
    do k = 2, nz
      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * z_dz(k-1)) &
        & / ( GasRDry * z_Temp(k-1) )
    end do


  end subroutine initialdata_toon2002_basic
  
end module initialdata_Toon2002
