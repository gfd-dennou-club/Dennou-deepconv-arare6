program main

  use ChemData
  use ChemCalc
  implicit none
  
  integer, parameter :: SpcNum = 3
  integer, parameter :: XMin = 0
  integer, parameter :: XMax = 40
  integer, parameter :: ZMin = 0
  integer, parameter :: ZMax = 1
  integer            :: SpcID(SpcNum)
  character(20)      :: SpcSymbol(SpcNum)
  real(8)            :: Temp(XMin:Xmax,ZMin:ZMax)
  real(8)            :: Svap(XMin:Xmax,ZMin:ZMax,SpcNum)
  real(8)            :: SvapAMP(XMin:Xmax,SpcNum)
  real(8)            :: tetens(2)
!  real(8)           :: LatentHeat, LatentHeat2
  integer            :: i, s


  call ChemData_Init

  SpcSymbol = (/ 'H2O-l  ',   &
    &            'H2O-s  ',   &
    &            'NH3-s  ' /)
  SpcID(:)  = ChemData_SpcID(SpcNum, SpcSymbol)

  write(*,*) SpcSymbol
  write(*,*) SpcID

  call ChemCalc_Init(XMin, XMax, ZMin, ZMax, SpcNum, SpcID)  

  !計算に利用する飽和蒸気圧と, 飽和蒸気圧との測定値をチェック. 
  call CheckSvapPress()

  do i = XMin, XMax    
    Temp(i,:) = i * 10.0d0 + 50.0d0
  end do

  do s = 1, SpcNum
    Svap(:,:,s) = xz_SvapPress( SpcID(s), Temp )
!    write(*,*) temp, Svap(1:2), SvapAMP(1:2), tetens(1:2)
  end do

  do s = 1, SpcNum
    do i = XMin, XMax    
!!      SvapAMP(i,s) = SvapPressAMP( SpcID(s), Temp(i,1) )
    end do
  end do

!  Temp2 = 273.15d0


  
!  write(*,*) xz_LatentHeat( SpcID(1), Temp )

  do i = XMin, XMax  
    tetens(1) = 6.11d0 * ( 10 ** ( 7.5d0 * (Temp(i,1) - 273.15d0) / ( 237.3d0 + Temp(i,1) - 273.15d0  ) ) ) * 1.0d2      
    tetens(2) = 6.11d0 * ( 10 ** ( 9.5d0 * (Temp(i,1) - 273.15d0) / ( 265.3d0 + Temp(i,1) - 273.15d0  ) ) ) * 1.0d2      
!    tetens(1) = 6.11d0 * dexp( 17.269 * (Temp(i,1) - 273.15d0) / (237.3d0 + Temp(i,1) - 273.15d0) ) * 1.0d2      
    write(*,*) Temp(i,1), Svap(i,1,:), SvapAMP(i,:), tetens
 end do

!  do s = 1, SpcNum
!    LatentHeat  = LatentHeatPerMol(SpcID(s), Temp2)
!    LatentHeat2 = LatentHeat(SpcID(s), Temp2)
!    SvapPress2  = SvapPress(SpcID(s), Temp2)
!    write(*,*) SvapPress2 
!    write(*,*) LatentHeat, LatentHeat2
!  end do
  

  contains
!!!
!!! 飽和蒸気圧チェックルーチン
!!!
!!!==========================================================================
  subroutine CheckSvapPress()
    !
    !飽和蒸気圧のチェックサブルーチン. 
    !蒸気圧の測定値と計算値の差を出力する.
    !

    !モジュール読み込み
    use ChemData, only : ChemData_DataNum,         & !データ点数 
      &                  ChemData_SVapPress,       & !飽和蒸気圧
      &                  ChemData_SVapPress_Temp     !温度
    
    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    integer             :: SNum
    integer             :: DNum
    real(8)             :: SVapPressData(SpcNum, ChemData_DataNum)
    real(8)             :: Temp(SpcNum, ChemData_DataNum)
    real(8)             :: SVapPressCalc
    real(8)             :: SVapPressCalc2

    !データベースから情報を取得する
    do SNum = 1, SpcNum
      SvapPressData(SNum,:) = ChemData_SvapPress(SpcID(SNum),:) 
      Temp(SNum,:)          = ChemData_SvapPress_Temp(SpcID(SNum),:) 
    end do

    do SNum = 1, SpcNum
      if (Phase(SNum) == 'Gas') cycle   !気相の場合
      
      write(*,*) "SpcName:      ", SpcID(SNum)
      
      do DNum = 1, ChemData_DataNum
        if ( Temp(SNum,DNum) == 0.0d0) cycle
        SvapPressCalc = SvapPress( SpcID(SNum), Temp(SNum,DNum) )
!!        SvapPressCalc2 = SvapPressAMP( SpcID(SNum), Temp(SNum,DNum) )
        
        write(*,*) "Temp:         ", Temp(SNum,DNum)
        write(*,*) "SVapPress:    ", SVapPressData(SNum,DNum), SVapPressCalc, SVapPressCalc2
        write(*,*) "DelSVapPress: ",                                                   &
          & abs(SVapPressData(SNum,DNum) - SVapPressCalc ) / SVapPressData(SNum,DNum), &
          & abs(SVapPressData(SNum,DNum) - SVapPressCalc2) / SVapPressData(SNum,DNum)
      end do
      SVapPressCalc = 0.0d0 !初期化
    end do
    
  end subroutine CheckSvapPress
  


end program main
