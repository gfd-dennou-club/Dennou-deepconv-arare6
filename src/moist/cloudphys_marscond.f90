!---------------------------------------------------------------------
!     Copyright (C) GFD Dennou Club, 2004. All rights reserved.
!---------------------------------------------------------------------
!= Module CloudPhys_marscond
!
!   * Developer: KITAMORI Taichi, YAMASHITA Tatsuya, SUGIYAMA Ko-ichiro
!   * Version: $ 
!   * Tag Name: $Name:  $
!   * Change History: 


module Cloudphys_MarsCond
  !
  ! 拡散によって雲粒が成長する場合における成長方程式を解く. 
  ! 蒸発量が存在する雲の量よりも多い場合, 蒸発量を存在する雲の量にする. 
  ! 雲密度を時間発展(凝結部分)を計算. 
  ! 雲が非常に小さくなったときに蒸発が起きなくなってしまっていたので,
  ! 条件分岐を書き直した.

  !モジュール呼び出し
  use dc_types,   only: DP, STRING

  !暗黙の型宣言禁止
  implicit none

  !変数定義
  real(DP), save, private :: DensIce     = 1.565d3  ! 固相の密度 [kg/m^3]
  real(DP), save, private :: NumAerosol  = 0.0d0  ! エアロゾルの数密度 [1/kg]
  real(DP), save, private :: RadiAerosol = 0.0d0  ! エアロゾルの数密度 [1/kg]
  real(DP), save, private :: Kd          = 0.0d0  ! 大気の熱伝導係数 [W/K m]
  real(DP), save, private :: SatRatioCr  = 0.0d0  ! 臨界飽和比 []
  real(DP), save, private :: SatRtWetAdia = 0.0d0 ! 湿潤断熱線の飽和比 []
  real(DP), save, private :: CO2LatHeat  = 0.0d0  ! 単位質量あたりの凝結熱 [J/kg]
  real(DP), save, private :: AntA        = 27.4d0 
  real(DP), save, private :: AntB        = 3103.0d0
  real(DP), save, private :: Pi          = 3.1415926535897932385d0 ! 円周率
  real(DP), save, private :: CDensCr     = 5.0d-5

  !公開要素
  public cloudphys_marscond_init
  public cloudphys_marscond_forcing

contains
  
  subroutine cloudphys_marscond_init
    !
    ! NAMELIST から物理定数を読み込む
    !
    
    !モジュール呼び出し
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use constants,     only : GasRDry           ! 気体定数
    use namelist_util, only : namelist_filename

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer       :: unit

    !-------------------------------------------------------------
    ! NAMELIST
    !
    NAMELIST /cloudphys_marscond_nml/         &
      & DensIce, NumAerosol, RadiAerosol, Kd, & 
      & SatRatioCr, SatRtWetAdia, CDensCr
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=cloudphys_marscond_nml)
    close(unit)    

    !-------------------------------------------------------------
    ! 潜熱. 定数と見なせる. 
    !   クラウジウスクラペイロンの関係式 dp/dT = L p / (R T^2) と 
    !   Antoine の式 \ln p = A - B/T を利用. 
    !
    CO2LatHeat = AntB * GasRDry

    !-------------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", &
      & "cloudset_init", "DensIce = %f",  d=(/DensIce/) )
    call MessageNotify( "M", &
      & "cloudset_init", "NumAerosol = %f",  d=(/NumAerosol/) )
    call MessageNotify( "M", &
      & "cloudset_init", "RadiAerosol = %f",  d=(/RadiAerosol/) )
    call MessageNotify( "M", &
      & "cloudset_init", "Kd = %f",  d=(/Kd/) )
    call MessageNotify( "M", &
      & "cloudset_init", "SatRatioCr = %f",  d=(/SatRatioCr/) )
    call MessageNotify( "M", &
      & "cloudset_init", "SatRtWetAdia = %f",  d=(/SatRtWetAdia/) )

    !------------------------------------------------------------
    ! 変数定義
    !
    call cloudphys_marscond_historyauto

  end subroutine cloudphys_marscond_init
  

!!!----------------------------------------------------
  subroutine cloudphys_marscond_forcing(  &
    &  xyz_PTempNs,         &  !(in) 温位
    &  xyz_ExnerNs,         &  !(in) エクスナー関数
    &  xyz_CDensNs,         &  !(in) 
    &  xyz_DPTempDtNl,      &  !(in)    
    &  xyz_DExnerDtNl,      &  !(in)    
    &  xyz_DCDensDtNl,      &  !(in)    
    &  xyz_PTempAs,         &  !(out) 
    &  xyz_CDensAs,         &  !(out) 雲密度
    &  xyz_DExnerDtNs       &  !(out) 
    & )
    !
    ! tendency の計算 & 温位と雲密度の更新
    !
          
    !モジュール呼び出し
    use dc_types,      only: DP
    use gtool_historyauto,                    &
      &               only : HistoryAutoPut
    use gridset,      only : imin, imax,      &
      &                      jmin, jmax,      &
      &                      kmin, kmax,      &
      &                      nx, ny, nz
    use basicset,     only : xyz_PTempBZ,     &! 温位基本場
      &                      xyz_ExnerBZ,     &! 無次元圧力
      &                      xyz_VelSoundBZ,  &! 音速
      &                      xyz_VPTempBZ,    &! 仮温位
      &                      xyz_DensBZ        ! 密度
    use constants,    only : GasRDry,         & ! 気体定数
      &                      PressBasis,      & ! 温位の基準圧力
      &                      CpDry              ! 定圧比熱
    use SetMargin,    only : SetMargin_xyz
    use timeset,      only : TimeN, DelTimeShort


    ! 暗黙の型宣言を禁止
    implicit none
    
    ! Input
    real(DP), intent(in)   :: xyz_PTempNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 無次元圧力 [1]
    real(DP), intent(in)   :: xyz_ExnerNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 
    real(DP), intent(in)   :: xyz_CDensNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 
    real(DP), intent(in)   :: xyz_DPTempDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(in)   :: xyz_DExnerDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(in)   :: xyz_DCDensDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(out)  :: xyz_PTempAs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 温位 [K]
    real(DP), intent(out)  :: xyz_CDensAs(imin:imax, jmin:jmax, kmin:kmax)   
                                        ! 雲の密度   [kg/m^3]
    real(DP), intent(out)  :: xyz_DExnerDtNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 凝結熱による温度変化率 [K/s]
  
    ! Work
    real(DP)               :: xyz_RadiCloud(imin:imax, jmin:jmax, kmin:kmax)
                                        ! 雲粒の半径 [m]
    real(DP)               :: xyz_SatRatio(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Rh(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_TempAll(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Mcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Qcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_McondTmp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_PIcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Zero(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_PTempTmp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_CDensTmp(imin:imax, jmin:jmax, kmin:kmax)
    logical                :: xyz_Mask(imin:imax, jmin:jmax, kmin:kmax)
    integer                :: i,j,k        ! ループ変数
 

!    ! 時間発展
!    !
!    xyz_PTempTmp = xyz_PTempNs + DelTimeShort * xyz_DPTempDtNl
!    xyz_CDensTmp = xyz_CDensNs + DelTimeShort * xyz_DCDensDtNl

!    ! Set Margin
!    !
!    call SetMargin_xyz(xyz_PTempTmp)
!    call SetMargin_xyz(xyz_CDensTmp)

!    ! 移流で負になった部分を穴埋め
!    ! 
!    call FillNegativeDensity(xyz_CDensTmp)

!    ! Set Margin
!    !
!    call SetMargin_xyz(xyz_CDensTmp)

    ! tendency は Ns の値で計算
    !
    xyz_TempAll = (xyz_ExnerBZ + xyz_ExnerNs) * (xyz_PTempBZ + xyz_PTempNs)
    xyz_Mask = .false. 
    xyz_Zero = 0.0d0
   
    ! 熱輸送に関する係数
    ! 
    xyz_Rh = (CO2LatHeat**2.0d0) / (Kd * GasRDry * (xyz_TempAll**2.0d0))

    ! 飽和比 (1.36) 式. 
    !
    xyz_SatRatio =                                                   &
      &  PressBasis * (xyz_ExnerBZ + xyz_ExnerNs)**(CpDry / GasRDry) &
      &  / exp( AntA - AntB / xyz_TempAll )  

    ! 雲粒半径. 
    !
    xyz_RadiCloud =   &
      & (  &
      &     RadiAerosol**3.0d0  &
      &   + 3.0d0 * xyz_CDensNs / (4.0d0 * Pi * DensIce * xyz_DensBZ * NumAerosol) & 
      &  ) ** (1.0d0 / 3.0d0)

    ! 凝結量の計算
    !
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          if (xyz_SatRatio(i,j,k) - SatRatioCr > epsilon(0.0d0)) then
            ! 飽和比が定めた臨界飽和比よりも大きい場合
            ! CO2 の飽和比の式が地球大気で用いられる飽和比の式と等価なものなのか確認する必要がある. 
            !
            xyz_Mask(i,j,k) = .true. 
            
          else if (xyz_CDensNs(i,j,k) > CDensCr ) then
            ! 臨界飽和比を超えていないが, 雲密度が閾値以上である場合
            ! 飽和比が 1 以上のときに凝結, 1 以下のときに蒸発. 
            !
            xyz_Mask(i,j,k) = .true. 
            
          else if (xyz_CDensNs(i,j,k) /= 0.0d0 .and. xyz_SatRatio(i,j,k) < 1.0d0 ) then
            ! 雲密度が閾値未満で, 飽和比 1.0 未満である場合に蒸発. 
            !
            xyz_Mask(i,j,k) = .true. 

          end if
        end do
      end do
    end do

    ! 凝結量の仮値を計算. 
    !
    xyz_McondTmp =                                                &
      & max( - xyz_CDensTmp / DelTimeShort,                       &
      &      4.0d0 * Pi * xyz_RadiCloud * xyz_DensBZ * NumAerosol &
      &        / xyz_Rh * (xyz_SatRatio - 1.0d0)                  &
      &   )

    ! 凝結量を計算. Mask が .false. な要素にはゼロを入れる
    !
    xyz_Mcond = merge(xyz_McondTmp, xyz_Zero, xyz_Mask)

    ! 密度の計算 
    !
    xyz_CDensAs = xyz_CDensTmp + DelTimeShort * xyz_Mcond

    ! Set Margin
    !
    call SetMargin_xyz(xyz_CDensAs)

    ! 温位変化の計算
    ! 
    xyz_Qcond = CO2LatHeat * xyz_Mcond / (CpDry * xyz_DensBZ * xyz_ExnerBZ)
    xyz_PTempAs = xyz_PTempTmp + DelTimeShort * xyz_QCond

    ! Set Margin
    !
    call SetMargin_xyz(xyz_PTempAs)
    
    ! エクスナー関数の時間微分
    ! 
    xyz_PIcond =                                                   &
      &  ( xyz_VelSoundBZ ** 2.0d0 )                               &
      &     / (CpDry * xyz_DensBZ * (xyz_VPTempBZ ** 2.0d0))        &
      &     * xyz_Mcond                                            &
      &     * ( CO2LatHeat / (CpDry * xyz_ExnerBZ) - xyz_VPTempBZ ) 

    xyz_DExnerDtNs = xyz_DExnerDtNl + xyz_PIcond
    
    call HistoryAutoPut(TimeN, 'PTempCond', xyz_Qcond(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'ExnerCond', xyz_PIcond(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'CDensCond', xyz_Mcond(1:nx, 1:ny, 1:nz))

    call SetMargin_xyz(xyz_PTempAs)
    call SetMargin_xyz(xyz_CDensAs)
!    call SetMargin_xyz(xyz_DExnerDtNs)

  end subroutine Cloudphys_marscond_forcing


  subroutine Cloudphys_marscond_historyauto
    !
    ! tendency の出力設定
    !

    !モジュール呼び出し
    use gtool_historyauto, only: HistoryAutoAddVariable

    !暗黙の型宣言禁止
    implicit none

    !------------------------------------------------------
    ! tendency の定義
    !
    call HistoryAutoAddVariable(  &
      & varname='PTempCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of exner function', &
      & units='s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='CDensCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of cloud density', &
      & units='K.m-3.s-1',    &
      & xtype='float')

  end subroutine Cloudphys_marscond_historyauto
  
end module Cloudphys_MarsCond
