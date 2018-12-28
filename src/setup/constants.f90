!= 定数モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: constants.f90,v 1.12 2014/01/21 05:00:57 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module constants
  !
  != 定数モジュール
  !

  !モジュール読み込み
  !
  use dc_types,      only: DP, STRING

  !暗黙の型宣言禁止
  !
  implicit none

  ! private 属性
  !
  private

  !Public Interface
  !
  real(DP), save, public :: Grav = 9.8d0          !重力 [m/s^2]
  real(DP), save, public :: PressBasis = 965.0d0  !温位の基準圧力 [Pa]
  real(DP), save, public :: TempSfc = 0.0d0       !地表面温度 [K]
  real(DP), save, public :: PressSfc = 0.0d0      !地表面圧力 [Pa]
  real(DP), save, public :: TempTop = 0.0d0       !上部境界の温度 [K]
  real(DP), save, public :: PressTop = 0.0d0      !上部境界での圧力 [Pa]
  real(DP), save, public :: CpDry  = 0.0d0        !乾燥成分の定圧比熱 [J/K kg]
  real(DP), save, public :: CpDryMol = 0.0d0      !乾燥成分の定圧比熱 [J/K kg]
  real(DP), save, public :: CvDry = 0.0d0         !乾燥成分の定積比熱 [J/K kg]
  real(DP), save, public :: MolWtDry = 0.0d0      !乾燥成分の分子量   [kg/mol]
  real(DP), save, public :: GasRDry  = 0.0d0      !乾燥成分の気体定数 [J/K kg]
  real(DP), save, public :: DayTime = 86400.0d0   ! 1 日の長さ [s]
  real(DP), save, public :: FactorJ = 1.0d0       !雲物理過程のパラメータ
                                                  !木星では 3.0d0
                                                  !地球では 1.0d0 とする
  ! サブルーチンの公開
  !
  public constants_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine constants_Init
    !
    != 初期化ルーチン
    !
    ! namelist の設定に基づいて物性を決める. 順序は以下の通り. 
    !
    ! * namelist から定圧比熱 (CpDry) と平均分子量 (MolWtDry) が与えられた場合は, 
    !   それらを元に気体定数(GasRDry), 定積比熱 (CvDry), 定圧モル比熱 (CpDryMol) を決める. 
    ! * namelist から定圧比熱 (CpDry) と気体定数 (GasRDry) が与えられた場合は, 
    !   それらを元に平均分子量(MolWtDry), 定積比熱 (CvDry), 定圧モル比熱 (CpDryMol) を決める. 
    ! * namelist から乾燥成分のモル比 (SpcDryMolFr) が与えられた場合には, 
    !   それと熱力学テーブルより, 定圧比熱 (CpDry), 気体定数 (GasRDry), 
    !   平均分子量(MolWtDry), 定積比熱 (CvDry), 定圧モル比熱 (CpDryMol) を決める. 
    ! 
    
    !モジュール読み込み
    use dc_types,      only: DP, STRING
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify          !メッセージ出力
    use ChemData,      only: GasRUniv,             &!普遍気体定数
      &                      ChemData_OneSpcID,    &!化学種の ID
      &                      ChemData_CpPerMolRef, &!標準状態での比熱
      &                      ChemData_MolWt         !分子量
    use namelist_util, only: namelist_filename

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    integer                  :: SpcDryNum = 1   !乾燥成分の化学種の数
    character(20)            :: SpcDrySymbol(5) !乾燥成分の化学種名
    real(DP)                 :: SpcDryMolFr(5)  !乾燥成分の化学種の存在度
    integer, allocatable     :: SpcDryID(:)     !乾燥成分の化学種のID
    real(DP), allocatable    :: PropertyDry(:)  !作業配列
    integer                  :: s               !作業変数
    integer                  :: unit            !装置番号
    logical                  :: flag = .false.
    character(STRING)        :: Planet = ""
     
    !NAMELIST の定義
    NAMELIST /constants_nml/ &
      & Planet, Grav, PressBasis, TempSfc, PressSfc, TempTop, PressTop, & 
      & SpcDrySymbol, SpcDryMolFr, DayTime, CpDry, MolWtDry, GasRDry
 
    SpcDrySymbol = '' 
    SpcDryMolFr  = 0.0d0
    
    !ファイルオープン. 情報取得. 
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=constants_nml)
    close(unit)

    ! 重力加速度
    !
    if (trim(Planet) == "Earth") then 
      FactorJ = 1.0d0
      Grav    = 9.8d0
    elseif (trim(Planet) == "Jupiter") then 
      FactorJ = 3.0d0
      Grav    = 23.1d0
    end if

    ! 比熱, 分子量, 気体定数 
    !
    if (CpDry /= 0.0d0 .AND. MolWtDry /= 0.0d0) then 
    ! namelist から CpDry の値が入力された場合

      !定圧モル比熱
      CpDryMol = CpDry *  MolWtDry 
      
      !気体定数
      GasRDry = GasRUniv / MolWtDry
      
      !定積比熱
      CvDry    = CpDry - GasRDry

    elseif (CpDry /= 0.0d0 .AND. GasRDry /= 0.0d0) then 
    ! namelist から CpDry と GasRDry の値が入力された場合

      !乾燥成分の分子量
      MolWtDry = GasRUniv / GasRDry 

      !定圧モル比熱
      CpDryMol = CpDry *  MolWtDry 
      
      !定積比熱
      CvDry    = CpDry - GasRDry

    elseif (SpcDryMolFr(1) /= 0.0d0) then 
    ! namelist から モル比が入力された場合

      flag = .true.

      !----------------------------------------------------------
      ! 乾燥成分の物性値の初期化
      !
      !乾燥成分の個数を数える
      SpcDryNum = count(SpcDrySymbol /= "")
      
      !化学種の ID を取得    
      allocate(SpcDryID(SpcDryNum))    
      do s = 1, SpcDryNum
        SpcDryID(s) = ChemData_OneSpcID( SpcDrySymbol(s) )
      end do
      
      !作業配列の準備
      allocate(PropertyDry(SpcDryNum))
      
      !分子量
      do s = 1, SpcDryNum
        PropertyDry(s) = ChemData_MolWt(SpcDryID(s))
      end do
      MolWtDry = dot_product(PropertyDry, SpcDryMolFr(1:SpcDryNum)) 
      
      !定圧比熱(モル当量)
      do s = 1, SpcDryNum    
        PropertyDry(s) = ChemData_CpPerMolRef(SpcDryID(s))
      end do
      CpDryMol = dot_product(PropertyDry, SpcDryMolFr(1:SpcDryNum)) 
      
      !定圧比熱
      CpDry    = CpDryMol / MolWtDry
      
      !気体定数
      GasRDry = GasRUniv / MolWtDry
      
      !定積比熱
      CvDry    = CpDry - GasRDry

    else
      call MessageNotify( "E", "constants_init", "Enough Variables are not set" )
    end if
    
    !----------------------------------------------------------
    ! 確認
    !----------------------------------------------------------
    call MessageNotify( "M", &
      & "constants_init", "Grav = %f", d=(/Grav/) )
    call MessageNotify( "M", &
      &  "constants_init", "FactorJ = %f",  d=(/FactorJ/) )
    call MessageNotify( "M", &
      & "constants_init", "PressBasis = %f", d=(/PressBasis/))
    if (TempSfc /= 0) then 
      ! 地表面温度が与えられた場合 
      !
      call MessageNotify( "M", &
        & "constants_init", "TempSfc = %f",  d=(/TempSfc/) )
      call MessageNotify( "M", &
        & "constants_init", "PressSfc = %f", d=(/PressSfc/) )
    end if
    if (TempTop /= 0) then 
      ! 上部境界の温度が与えられた場合 
      !
      call MessageNotify( "M", &
        & "constants_init", "TempTop = %f",  d=(/TempTop/) )
      call MessageNotify( "M", &
        & "constants_init", "PressTop = %f", d=(/PressTop/) )
    end if
    if (flag) then 
      ! モル比が与えられた場合にはプロットする. 
      !
      do s = 1, SpcDryNum
        call MessageNotify( "M", &
          &  "constants_init", "SpcDryID = %d",      i=(/SpcDryID(s)/))
        call MessageNotify( "M", &
          &  "constants_init", "SpcDrySymbol = %c", c1=trim(SpcDrySymbol(s)))
        call MessageNotify( "M", &
          &  "constants_init", "SpcDryMolFr = %f",   d=(/SpcDryMolFr(s)/))
      end do
    end if
    
    call MessageNotify( "M", "constants_init", "CpDry    = %f",    d=(/CpDry/) )
    call MessageNotify( "M", "constants_init", "CpDryMol = %f", d=(/CpDryMol/) )
    call MessageNotify( "M", "constants_init", "CvDry    = %f",    d=(/CvDry/) )
    call MessageNotify( "M", "constants_init", "GasRDry  = %f",  d=(/GasRDry/) )
    call MessageNotify( "M", "constants_init", "MolWtDry = %f", d=(/MolWtDry/) )
    call MessageNotify( "M", "constants_init", "DayTime  = %f",  d=(/DayTime/)  )    

  end subroutine Constants_Init

end module Constants
