program main
  ! 
  ! Imamura and Hashimoto (1998) に基づいた計算モジュールの確認計算
  ! 物性値を確認するためのプログラム
  ! 
  
  use Imamura1998
  
  
  !暗黙の型宣言禁止
  implicit none
  
  !変数の宣言
  real(8) :: Temp, Press
  real(8) :: x1, x2
  real(8) :: GibbsRDivRT
  real(8) :: DelChemPotRDivRT1
  real(8) :: DelChemPotRDivRT2
  real(8) :: SatPress1
  real(8) :: SatPress2
  real(8) :: ChemData_SatPress2(16)
  real(8) :: ChemData_SatPress2_Temp(16)
  real(8) :: TempList(5), TempList2(7)
  integer :: i, j
  real(8) :: tmp1, tmp2
  real(8) :: NumDens1, NumDens2
  real(8) :: molfr1, molfr2
  
  real(8), parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  
!!!
!!! 実験データ
!!!
  !H2O 飽和蒸気圧 [Pa]
  ChemData_SatPress2 = (/ &
    & 610.749769736842d0, 1227.63236842105d0, 2337.94105263158d0, &
    & 4243.25101973684d0, 7376.72664473684d0, 12338.9851973684d0, &
    & 19923.6947368421d0, 31166.7700657895d0, 47364.1046052632d0, &
    & 70110.2338815789d0, 101325d0,           143268.217105263d0, &
    & 198530.338815789d0, 270111.118421053d0, 361356.947368421d0, &
    & 475974.1875d0 /)
  
  !H2O 飽和蒸気圧のデータに対応する温度
  ChemData_SatPress2_Temp = (/ &
    & 273.15d0, 283.15d0, 293.15d0, 303.15d0, 313.15d0, &
    & 323.15d0, 333.15d0, 343.15d0, 353.15d0, 363.15d0, &
    & 373.15d0, 383.15d0, 393.15d0, 403.15d0, 413.15d0, &
    & 423.15d0 /)
  
!!!
!!! H2SO4 の飽和蒸気圧の確認.
!!! 計算した結果を Excel でグラフ化し, Kulmala and Laaksonen (1990) Fig.1 と比較.
!!! 得られた結果は低温側で Fig.1 とファクター倍のずれがあるが気にしないことにする.
!!!
!!! Kulmala and Laaksonen, The Journal of Chemical Physics 93, 696 (1990);
!!! doi: 10.1063/1.459519
!!!
  write(*,*) "===================================="
  write(*,*) "  H2SO4 saturation vapor pressure   "
  write(*,*) "===================================="
  LoopH2SO4_1: do Temp = 150.0d0, 350.0d0, 1.0d0
    x1   = 1.0d0 - epsilon(1.0d0)
    call Imamura1998_SatPressRef( Temp, SatPress1, SatPress2 )
    call Imamura1998_SatPress( Temp, x1, tmp1, tmp2 )
!     write(*,*) "-----------------------------------"
!     write(*,*) "  Temp: ", Temp, "molfr : ", x1
!     write(*,*) "-----------------------------------"     
!     write(*,*) "SatPress = ", SatPress1, "( ", Press, " )"
!     write(*,*) "Ratio =         ", abs(SatPress1 / Press)
!     write(*,*) "AbsoluteError = ", abs(SatPress1 - Press)
!     write(*,*) "RelativeError = ", abs((SatPress1 - Press)/SatPress1)
!     write(*,*) ""
    write(*,*) int(Temp), ",", real(SatPress1), ",", real(tmp1), ","
  end do LoopH2SO4_1

  
!!!
!!! H2O の飽和蒸気圧の確認.
!!! 計算した値を飽和蒸気圧のデータ (化学工学便覧) と比較する.
!!! 得られた結果は化学工学便覧のデータと一致する. 
!!! 
  write(*,*) "===================================="
  write(*,*) "  H2O saturation vapor pressure     "
  write(*,*) "===================================="
  LoopH2O: do i = 1, 16
    Temp = ChemData_SatPress2_Temp(i)
    call Imamura1998_SatPressRef( Temp, SatPress1, SatPress2 )
    write(*,*) "-----------------------------------"
    write(*,*) "  Temp: ", temp
    write(*,*) "-----------------------------------"     
    write(*,*) "SatPress = ", SatPress2, "( ", ChemData_SatPress2(i), " )"
    write(*,*) "AbsoluteError = ", abs(SatPress2 - ChemData_SatPress2(i))
    write(*,*) "RelativeError = ", abs((SatPress2 - ChemData_SatPress2(i))/SatPress2)
    write(*,*) ""
  end do LoopH2O

!  write(*,*) "===================================="
!  write(*,*) "  H2O saturation vapor pressure (2) "
!  write(*,*) "===================================="
!  LoopH2O_1: do Temp = 100.0d0, 1000.0d0, 50.0d0
!     write(*,*) Temp
!     call Imamura1998_SatPressRef( Temp, SatPress1, SatPress2 )
!     write(*,*) "H2SO4 : ", real(Temp), real(SatPress1), real(SatPress1 / Temp / Boltz)
!     write(*,*) "H2O   : ", real(Temp), real(SatPress2), real(SatPress2 / Temp / Boltz)
!  end do LoopH2O_1
  
!!! 
!!! CHECK: relative gibbs free energy, relative potential energy
!!! Zeleznik (1999) Table 7 の数字との比較用
!!!
  write(*,*) "===================================="
  write(*,*) "  -G^r/RT, -\mu^r/RT     "
  write(*,*) "===================================="
  TempList = (/ 200.0d0, 250.0d0, 298.15d0, 300.0d0, 350.0d0/)
  LoopTemp1: do j = 1, 5
    Temp = TempList(j)
    write(*,*) "-----------------------------------"
    write(*,*) "  Temp: ", temp
    write(*,*) "-----------------------------------"
    write(*,*) "  x1                -G^r/RT,           -\mu^r_1/RT ,    -\mu^r_2/RT,"
    loop1: do i = 0, 50
      x1 = 0.02d0 * i
      x2 = 1.0d0 - x1
      
      call Imamura1998_GibbsRDivRT( temp, x1, GibbsRDivRT ) 
      call Imamura1998_DelChemPotRDivRT( temp, x1, DelChemPotRDivRT1, DelChemPotRDivRT2 )
      write(*,*) real(x1), ",",                     &
        & real(-1.0d0 * GibbsRDivRT), ",",       &
        & real(-1.0d0 * DelChemPotRDivRT1), ",", &
        & real(-1.0d0 * DelChemPotRDivRT2), ","
    end do loop1
  end do LoopTemp1
  
end program main
