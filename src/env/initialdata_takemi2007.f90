!= Module initialdata_takemi2007
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_takemi2007.f90,v 1.8 2014/07/08 00:59:09 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_takemi2007
  !
  ! Takemi (2007) を模したの基本場・初期場を設定する

  !モジュール読み込み
  use dc_types, only: DP

  !暗黙の型宣言禁止
  implicit none

  !共通変数の定義
  real(DP), parameter, private :: PressSfcTakemi = 1.0d5   ! 地表の圧力
  real(DP), parameter, private :: PTempSfcTakemi = 300.0d0 ! 地表の温位
  real(DP), parameter, private :: AltTr    = 1.2d4   ! 対流圏界面高度
  real(DP), parameter, private :: HumMin   = 0.25d0  ! 混合比一定な高度
  integer,  parameter, private :: SpcID = 6          ! 水の番号
  real(DP), save, private      :: QMixSfc            ! 地表面でのモル比
  real(DP), save, private      :: VelXSfc            ! 地表の速度
  real(DP), save, private      :: PTempTr = 0.0d0    ! 
  real(DP), save, private      :: DryFact = 0.0d0    ! 
  real(DP), save, private      :: Alt1    = 0.0d0    ! 
  real(DP), save, private      :: Alt2    = 0.0d0    ! 
  integer,  save, private      :: amin, amax

  !初期化だけ公開
  public  initialdata_takemi2007_init
  public  initialdata_takemi2007_basic
  public  initialdata_takemi2007_wind

contains

!!!------------------------------------------------------------------------------!!!
  subroutine initialdata_takemi2007_init
    !
    !設定ファイルから出力ファイルに記載する情報を読み込む
    !

    !モジュール読み込み
    use dc_types,     only: STRING, DP
    use dc_iounit,    only: FileOpen 
    use dc_message,   only: MessageNotify
    use namelist_util,only: namelist_filename
    use gridset,      only: nz
    use axesset,      only: z_Z                !スカラー格子点での高度
    use constants,    only: PressBasis,       &!温位の基準圧力
      &                     TempSfc,          &!地表面温度
      &                     PressSfc           !地表面圧力
    
    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    real(DP)            :: Alt = 0.0d0  !高度
    integer             :: unit         !設定ファイル用装置番号  
    integer             :: k
    integer             :: ID_BasicZ = 0
    integer             :: ID_Wind   = 0
    integer, parameter  :: ID_MidLat_Q10     = 1
    integer, parameter  :: ID_MidLat_Q12     = 2
    integer, parameter  :: ID_MidLat_Q14     = 3
    integer, parameter  :: ID_MidLat_Q16     = 4
    integer, parameter  :: ID_MidLat_Q16DRY1 = 5
    integer, parameter  :: ID_MidLat_Q16DRY2 = 6
    integer, parameter  :: ID_MidLat_Q18     = 7
    integer, parameter  :: ID_Tropic_Q18     = 8
    integer, parameter  :: ID_Tropic_Q18DRY1 = 9
    integer, parameter  :: ID_Tropic_Q18DRY2 = 10
    integer, parameter  :: ID_Tropic_Q18DRY3 = 11
    integer, parameter  :: ID_Wind_LowLevel    = 1
    integer, parameter  :: ID_Wind_MiddleLevel = 2
    integer, parameter  :: ID_Wind_HighLevel   = 3

    character(STRING)   :: FlagEnv = ""
    character(STRING)   :: FlagWind = ""


    !設定ファイルから読み込む出力ファイル情報
    NAMELIST /initialdata_takemi2007_nml/ FlagEnv, FlagWind, VelXSfc
    
    !設定ファイルから出力ファイルに記載する情報を読み込む
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_takemi2007_nml)
    close(unit)

    ! 確認. 地表面温度・圧力を指定するのは別モジュールなので. 
    !
    if (     PressBasis /= PressSfcTakemi &
      & .OR. PressSfc /= PressSfcTakemi   &
      & .OR. TempSfc  /= PTempSfcTakemi  ) then 
      
      call MessageNotify( "E", "initaldata_takemi2007_init", &
        & "Constants are wrong. please PressSfc = 1.0d5, TempSfc = 300.0d0")
    end if
    
    call MessageNotify( "M", "initaldata_takemi2007_init", &
      & "VelXSfc= %f", d=(/VelXSfc/) ) 

    !基本場の選択
    !
    if (FlagEnv == "MidLat_Q10") then 
      ID_BasicZ = ID_MidLat_Q10
    elseif (FlagEnv == "MidLat_Q12") then 
      ID_BasicZ = ID_MidLat_Q12
    elseif (FlagEnv == "MidLat_Q14") then 
      ID_BasicZ = ID_MidLat_Q14
    elseif (FlagEnv == "MidLat_Q16") then 
      ID_BasicZ = ID_MidLat_Q16
    elseif (FlagEnv == "MidLat_Q16DRY1") then 
      ID_BasicZ = ID_MidLat_Q16DRY1
    elseif (FlagEnv == "MidLat_Q16DRY2") then 
      ID_BasicZ = ID_MidLat_Q16DRY2
    elseif (FlagEnv == "MidLat_Q18") then 
      ID_BasicZ = ID_MidLat_Q18
    elseif (FlagEnv == "Tropic_Q18") then 
      ID_BasicZ = ID_Tropic_Q18
    elseif (FlagEnv == "Tropic_Q18DRY1") then 
      ID_BasicZ = ID_Tropic_Q18DRY1
    elseif (FlagEnv == "Tropic_Q18DRY2") then 
      ID_BasicZ = ID_Tropic_Q18DRY2
    elseif (FlagEnv == "Tropic_Q18DRY3") then 
      ID_BasicZ = ID_Tropic_Q18DRY3
    end if

    ! 速度場の選択
    !
    if (FlagWind == "LowLevel") then 
      ID_Wind = ID_Wind_LowLevel
    elseif (FlagWind == "MiddleLevel") then 
      ID_Wind = ID_Wind_MiddleLevel
    elseif (FlagWind == "HighLevel") then 
      ID_Wind = ID_Wind_HighLevel
    end if

    ! 温位場
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q10:ID_MidLat_Q18)
      PTempTr = 343.0d0
    case (ID_Tropic_Q18:ID_Tropic_Q18DRY3)
      PTempTr = 358.0d0
    end select

    ! 混合比
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q10) 
      QMixSfc = 0.010d0
    case (ID_MidLat_Q12) 
      QMixSfc = 0.012d0
    case (ID_MidLat_Q14) 
      QMixSfc = 0.014d0
    case (ID_MidLat_Q16:ID_MidLat_Q16DRY2) 
      QMixSfc = 0.016d0
    case (ID_MidLat_Q18:ID_Tropic_Q18DRY3)
      QMixSfc = 0.018d0
    end select

    ! 湿度変化
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q16DRY1)
      DryFact = - 0.13d0
      Alt = 2.5d3
    case (ID_MidLat_Q16DRY2)
      DryFact = - 0.30d0
      Alt = 2.5d3
    case (ID_Tropic_Q18DRY1)
      DryFact = - 0.20d0
      Alt = 2.5d3
    case (ID_Tropic_Q18DRY2)
      DryFact = - 0.20d0
      Alt = 5.0d3
    case (ID_Tropic_Q18DRY3)
      DryFact = - 0.20d0
      Alt = 7.5d3
    end select

    do k = 1, nz
      if (z_Z(k) < AltTr .AND. AltTr <= z_Z(k+1)) then 
        amax = k
      end if
      if (z_Z(k) < Alt .AND. Alt <= z_Z(k+1)) then 
        amin = k
      end if
    end do

    ! 水平風速
    !
    select case (ID_Wind)
    case (ID_Wind_LowLevel) 
      Alt1 = 0.0d0
      Alt2 = 2.5d3
    case (ID_Wind_MiddleLevel)
      Alt1 = 2.5d3
      Alt2 = 5.0d3
    case (ID_Wind_HighLevel)
      Alt1 = 5.0d3
      Alt2 = 7.5d3
    end select

  end subroutine initialdata_takemi2007_init


!!!------------------------------------------------------------------------------!!!
  subroutine  initialdata_takemi2007_basic( z_Temp, z_Press, zf_MolFr )
    !
    !== 概要
    ! * deepconv の地球用のテスト計算としてTakemi(2007)の再現計算を
    !   するための湿度の基本場を作成する
    !   * 基本場の温度の式が温位で与えられているため, 温度に変換する必要がある
    !

    !モジュール読み込み
    use dc_types,     only: DP
    use gridset,      only: kmin, kmax,       &!配列サイズ (Z 方向)
      &                     ncmax,            &!凝縮成分の数
      &                     nz                 !物理領域の大きさ (Z方向)
    use axesset,      only: z_Z,              &!スカラー格子点での高度
      &                     dz                 !鉛直格子間隔
    use constants,    only: PressBasis,       &!温位の基準圧力
      &                     GasRDry,          &!乾燥成分の定圧比熱
      &                     CpDry,            &!乾燥成分の定圧比熱
      &                     Grav,             &!重力加速度
      &                     TempSfc,          &!地表面温度
      &                     PressSfc,         &!地表面圧力
      &                     MolWtDry  
    use composition,  only: MolWtWet
    use chemcalc,     only: SvapPress 

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out):: z_Press(kmin:kmax)           !圧力
    real(DP), intent(out):: z_Temp(kmin:kmax)            !温度
    real(DP), intent(out):: zf_MolFr(kmin:kmax, 1:ncmax) !モル比
    real(DP) :: z_Hum(kmin:kmax)                         ! 相対湿度
    real(DP) :: z_PTemp(kmin:kmax)                       ! 温位
    real(DP) :: QMix
    integer  :: k

    !-------------------------------------------
    ! 初期化
    !
    z_Temp   = 1.0d-60
    z_PTemp  = 1.0d-60
    z_Press  = 1.0d-60
    z_Hum    = 0.0d0
    zf_MolFr = 0.0d0

    !----------------------------------------------
    ! 湿度
    !   論文中の式 (2) より計算
    !    
    do k = 1, nz
      if (z_Z(k) <= AltTr) then 
        z_Hum(k) = 1.0d0 - 0.75d0 * (z_Z(k) / AltTr) ** 1.25d0
      elseif (z_Z(k) > AltTr) then 
        z_Hum(k) = HumMin
      end if
    end do

    ! Fig.2b は明らかに相対湿度が 95% 程度で打ち止めになっているので, 
    ! 上限値を設けてみた
    !
    where (z_Hum > 0.95d0) 
      z_Hum = 0.95d0
    end where
    
    ! DRY ケースの場合. 
    !   DRY 以外では, DryFact = 0.0 になっている. 
    !   相対湿度の下限値を下回らないよう調整している.
    !
    do k = 1, nz
      if (amin < k .AND. k <= amax) then
        z_Hum(k) = z_Hum(k) + DryFact
        if (z_Hum(k) <= HumMin) then 
          z_Hum(k) = HumMin
        end if
      end if
    end do

    !----------------------------------------------
    ! 温位, 圧力, 温度
    !   温位は, 論文中の式 (1) より計算
    !
    !    
    z_PTemp(1) = TempSfc + (PTempTr - TempSfc) * ((z_Z(1) / AltTr) ** 1.25d0)
    z_Press(1) = PressSfc - (Grav * PressSfc * dz * 5.0d-1) / (GasRDry * TempSfc)
    z_Temp(1)  = z_PTemp(1) * (z_Press(1) / PressBasis) ** (GasRDry / CpDry)

    do k = 2, nz
      if (k <= amax) then 
        z_PTemp(k) = TempSfc + (PTempTr - TempSfc) * ((z_Z(k) / AltTr) ** 1.25d0)
      elseif (k > amax) then 
        z_PTemp(k) = z_PTemp(amax) * exp( Grav * (z_Z(k) - AltTr) / (CpDry * z_Temp(amax)))
      end if

      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * dz) / (GasRDry * z_Temp(k-1))
      
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressBasis) ** (GasRDry / CpDry))
    end do
    
    !----------------------------------------------
    ! モル比
    !   論文中には, in the lowest 1.5km depth で混合比の値を設定するという
    !   記述があったが, Fig2b から判断するに, 式 (2) で与えられる相対湿度から
    !   混合比を計算し, それが地表での混合比 (nml で与える) を超えたら, 
    !   地表での混合比と同じとしているのだろう. 
    !
    do k = 1, nz
      zf_MolFr(k,1) = SvapPress(SpcID, z_Temp(k)) * z_Hum(k) / z_Press(k)
      QMix = zf_MolFr(k,1) / MolWtDry * MolWtWet(1)
      
      if (QMix > QMixSfc) then 
        zf_MolFr(k,1) = QMixSfc * MolWtDry / MolWtWet(1)
      end if
    end do
    
  end subroutine Initialdata_takemi2007_basic


  subroutine initialdata_takemi2007_wind(pyz_VelX)
    !
    !-------------------------------------------------------------!
    ! シアーの設定 (Takemi,2007)                                  !
    !-------------------------------------------------------------!
    !
    != 概要
    !* case "Takemi2007" での計算時に鉛直シアーのある風を与える時に使用する
    !* 風の与え方には, 以下のようなバリエーションがある
    !  (1) シアーを与える高度を変える
    !  (2) シアーのある風の最大風速 (U_s) を変える
    !
    !  (1) については, (a) 0 - 2.5 km, (b) 2.5 - 5.0 km, (c) 5.0 - 7.5 km の
    !  三パターンがある
    !  (2) については, Takemi (2007) では熱帯場と中緯度場の温度場毎に
    !  異なる値を設定している
    !
    !  その強度(Us)は, 以下の通り
    !  <熱帯場>   (1) 5 m/s, (2) 10 m/s, (3) 15 m/s
    !  <中緯度場> (1) 10 m/s, (2) 15 m/s, (3) 20 m/s
    !
    !* シアーの形の模式図 (Takemi, 2007)   |
    !                                     /| 7.5 km
    !                                    / |
    !                                   /  |
    !                                  / ←|
    !                                 ｜  /| 5.0 km
    !                                 ｜ / |
    !                                 ｜/  |
    !                                  / ←|
    !                                 ｜  /| 2.5 km
    !                                 ｜ / |
    !                                 ｜/  |
    !                                  / ←|
    !---------------------------------+------------- 0.0 km
    !                                Us (m/s)
    !

    !モジュール読み込み
    use dc_types,     only: DP
    use gridset,      only: imin, imax,       &!配列の X 方向の上限
      &                     jmin, jmax,       &!配列の Y 方向の上限
      &                     kmin, kmax,       &!配列の Z 方向の上限
      &                     nz
    use axesset,      only: z_Z                !スカラー格子点での高度

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    integer               :: k
    
    !初期化
    pyz_VelX = 0.0d0
    
    do k = 1, nz
      if (z_Z(k) <= Alt1) then 
        pyz_VelX(:,:,k) = - VelXSfc
      elseif (z_Z(k) > Alt1 .AND. z_Z(k) <= Alt2) then 
        pyz_VelX(:,:,k) = - VelXSfc + (VelXSfc / (Alt2 - Alt1)) * (z_Z(k) - Alt1)
      end if
    end do
    
  end subroutine initialdata_takemi2007_wind
  
end module initialdata_takemi2007
