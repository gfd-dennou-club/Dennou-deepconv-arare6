!= Module HeatFlux
!
! Authors::   ODAKA Masatsugu
! Version::   $Id: surfaceflux_bulk_l1979.f90,v 1.1 2012/07/30 08:08:29 odakker Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module Surfaceflux_bulk_l1979
  !
  != 下部境界でのフラックスの計算モジュール
  !
  != Surface flux
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! Louis et al. (1979) の方法に基づいて地表面フラックスを計算. 
  !
  ! Surface fluxes are calculated by using the scheme by Louis et al. (1979). 
  !
  !== References
  !
  ! Louis, J-F., 1979:
  ! A parametric model of vertical eddy fluxes in the atmosphere. 
  ! Boundary Layer Meteorology, 17, 187-202.
  !== Variable List
  !
  !== Procedures List
  !

  ! モジュール引用; USE statement
  !  

  ! GTOOL 変数と手続き
  ! GTOOL variables and procedures
  !
  use dc_types, only: DP, STRING
  use dc_iounit, only: FileOpen
  use dc_message, only: MessageNotify
  use gtool_historyauto, only: HistoryAutoAddVariable, HistoryAutoPut

  ! 並列処理用変数
  ! Parallel processing variable 
  !       

  ! 格子点設定
  ! Grid points settings
  !
  use gridset, only: imin, & ! x 方向配列下限 
    &                      & ! Upper limit of array in x
    &                imax, & ! x 方向配列上限 
    &                      & ! Lower limit of array in x
    &                jmin, & ! y 方向配列下限 
    &                      & ! Upper limit of array in y
    &                jmax, & ! y 方向配列上限 
    &                      & ! Lower limit of array in y
    &                kmin, & ! z 方向配列下限 
    &                      & ! Upper limit of array in z
    &                kmax, & ! z 方向配列上限 
    &                      & ! Lower limit of array in z
    &                nx,   & ! x 方向格子点数 
    &                      & ! Number of grid point in x
    &                ny,   & ! y 方向格子点数 
    &                      & ! Number of grid point in y
    &                nz,   & ! z 方向格子点数 
    &                      & ! Number of grid point in z
    &                ncmax   ! 組成の数       
                             ! Number of spices

  ! 座標軸と演算子
  ! Axes and operator settings
  !
  use axesset, only: z_dz,        & ! z 方向格子点間隔 
    &                             & ! Grid size in z
    &                xyz_avr_pyz, & ! 平均操作
    &                             & ! Average operator
    &                xyz_avr_xqz, & ! 平均操作
    &                             & ! Average operator
    &                pyz_avr_xyz, & ! 平均操作
    &                             & ! Average operator
    &                xqz_avr_xyz    ! 平均操作
                                    ! Average operator

  ! 基本場変数
  ! Basic state variables
  !
  use basicset, only: xyz_ExnerBZ,  & ! 圧力関数
    &                               & ! Exner function
    &                 xyz_PressBZ,  & ! 圧力
    &                               & ! Pressure
    &                 xyz_PTempBZ,  & ! 温位
    &                               & ! Potential temperature
    &                 xyz_TempBZ,   & ! 温度
    &                               & ! Temperature
    &                 xyzf_QMixBZ,  & ! 混合比
    &                               & ! Mixing ration
    &                 xyz_DensBZ      ! 密度    
                                      ! Density

  ! 定数
  ! Constatns
  !
  use constants, only: Grav,       & ! 重力加速度
    &                              & ! Gravity
    &                  MolWtDry,   & ! 乾燥空気分子量 
    &                              & ! Molecular weight of dry air
    &                  PressBasis, & ! 温位の基準圧力
    &                              & ! Reference pressure
    &                  TempSfc,    & ! 地表面温度
    &                              & ! Surface temperature
    &                  PressSfc,   & ! 地表面気圧
    &                              & ! Surface pressure
    &                  CpDry,      & ! 乾燥空気定圧比熱 
    &                              & ! Specific heat of dry air
    &                  GasRDry       ! 乾燥空気気体定数
                                     ! Gas constant of dry air

  use constants0, only: FKarm        ! カルマン定数
                                     ! Karmann constant

  ! 組成に関する変数
  ! Setting for atmospheric composition
  ! 
  use composition, only: IdxCG,       & ! 凝結過程(気体)の配列添字 
    &                                 & ! Index of vapor
    &                    IdxCC,       & ! 凝結過程(雲)の配列添字 
    &                                 & ! Index of cloud
    &                    SpcWetID,    & !凝結成分の化学種ID
    &                                 & ! ID number of moist condensation spices
    &                    CondNum,     & ! 雲の数
    &                                 & ! Number of cloud component
    &                    MolWtWet,    & ! 凝結成分分子量
    &                                 & ! Molecular weight of moist air
    &                    SpcWetSymbol   ! 凝結成分の化学種名
                                        ! Chemical formula of condensation spices

  ! 化学量計算に関する設定
  ! Settings for chemical calculation
  !
  use chemcalc, only: SvapPress ! 飽和蒸気圧
                                ! Saturation vapor pressure

  ! NAMELIST に関する設定
  ! Settings for NAMELIST
  !
  use namelist_util, only: namelist_filename ! NAMELIST ファイル名
                                             ! NAMELIST file name
  ! 時刻に関する設定
  ! Setting for time
  !
  use timeset, only:  TimeN ! 時刻 t 
                            ! Time "t"

  ! 時間変化項
  ! Tendency of variable
  use DExnerDt, only: xy_DExnerDt_xy_xyf ! D$\pi$/Dt


  ! 暗黙の型宣言禁止
  ! Implicit none
  !
  implicit none

  ! 属性の指定
  !
  private

  ! 公開手続き
  ! Public procedure
  !
  public surfaceflux_bulk_init
  public surfaceflux_bulk_forcing

  ! 非公開変数
  ! Privete variables
  ! 
  logical, save:: FlagConstBulkCoef
                            ! Flag for using constant bulk coefficient
  logical, save:: FlagUseOfBulkCoefInNeutralCond
                            ! Flag for using bulk coefficient in neutral condition
  real(DP), save:: ConstBulkCoef
                            ! バルク係数一定値. 
                            ! Steady value of bulk coefficient
  real(DP), save :: VelMinForRi = 1.0d-8 
                            ! リチャード数計算用速度下限値
                            ! Lower limit of velocity for Ri
  real(DP), save :: SfcRoughLength = 1.0d-2
                            ! 祖度長さ
                            ! Roughness length
  real(DP), save :: VelBulkCoefMin = 0.0d0
                            ! $ u $ バルク係数最小値. 
                            ! Minimum value of $ u $ bulk coefficient
  real(DP), save :: TempBulkCoefMin = 0.0d0
                            ! $ T $ バルク係数最小値. 
                            ! Minimum value of $ T $ bulk coefficient
  real(DP), save :: QmixBulkCoefMin = 0.0d0
                            ! $ q $ バルク係数最小値. 
                            ! Minimum value of $ q $ bulk coefficient
  real(DP), save :: VelBulkCoefMax = 1.0d2
                            ! $ u $ バルク係数最大値. 
                            ! Maximum value of $ u $ bulk coefficient
  real(DP), save :: TempBulkCoefMax = 1.0d2
                            ! $ T $ バルク係数最大値. 
                            ! Maximum value of $ T $ bulk coefficient
  real(DP), save :: QmixBulkCoefMax = 1.0d2
                            ! $ q $ バルク係数最大値. 
                            ! Maximum value of $ q $ bulk coefficient

  real(DP), save :: DPTempDtFluxMin = 0.0d0
                            ! 温位フラックス最小値. 
                            ! Minimum value of potential temp. flux
  real(DP), save :: ExnerFluxMin = 0.0d0
                            ! 圧力関数フラックス最小値. 
                            ! Minimum value of exner function flux
  real(DP), save :: QmixFluxMin = 0.0d0
                            ! 混合比フラックス最小値. 
                            ! Minimum value of mixing ratio flux
  real(DP), save  :: Vel0 = 0.0d0  ! 下層での水平速度嵩上げ値
                                   ! 

  character(*), parameter:: module_name = 'surfaceflux_bulk_l1979'
                                   ! モジュールの名称.
                                   ! Module name

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Bulk_init
    !
    ! NAMELIST から必要な情報を読み取り, 時間関連の変数の設定を行う. 
    !

    ! 暗黙の型宣言禁止
    ! Implicit none
    !
    implicit none

    ! 作業変数
    ! Work variables
    !
    integer    :: l,   & ! 組成方向に回る DO ループ用作業変数
      &                & ! Work variables for DO loop in dimension of constituents
      &           unit   ! 出力装置番号
                         ! Device number 

    !---------------------------------------------------------------
    ! NAMELIST から情報を取得
    !
    NAMELIST /surfaceflux_bulk_nml/ &
      &  FlagConstBulkCoef,                                 &
      &  FlagUseOfBulkCoefInNeutralCond, ConstBulkCoef,     &
      ! 
      &  VelMinForRi, SfcRoughLength, Vel0,                 &
      !
      &  VelBulkCoefMin, TempBulkCoefMin, QmixBulkCoefMin,  &
      &  VelBulkCoefMax, TempBulkCoefMax, QmixBulkCoefMax             


    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=surfaceflux_bulk_nml)
    close(unit)  

    call MessageNotify( "M", module_name, "SfcRoughLength = %f",  &
      &  d=(/SfcRoughLength/))
    call MessageNotify( "M", module_name, "VelMinForRi = %f",     &
      &  d=(/VelMinForRi/) )
    call MessageNotify( "M", module_name, "Vel0 = %f", d=(/Vel0/))
    
    call MessageNotify( 'M', module_name, "FlagConstBulkCoef              = %b", l = (/ FlagConstBulkCoef /) )
    call MessageNotify( 'M', module_name, "FlagUseOfBulkCoefInNeutralCond = %b", l = (/ FlagUseOfBulkCoefInNeutralCond /) )
    call MessageNotify( 'M', module_name, "ConstBulkCoef   = %f", &
      &  d = (/ ConstBulkCoef   /) )
    call MessageNotify( 'M', module_name, "VelBulkCoefMin  = %f", &
      &  d = (/ VelBulkCoefMin  /) )
    call MessageNotify( 'M', module_name, "TempBulkCoefMin = %f", &
      &  d = (/ TempBulkCoefMin /) )
    call MessageNotify( 'M', module_name, "QmixBulkCoefMin = %f", &
      &  d = (/ QmixBulkCoefMin /) )
    call MessageNotify( "M", module_name, "VelBulkCoefMax = %f",  &
      &  d=(/VelBulkCoefMax/))
    call MessageNotify( "M", module_name, "TempBulkCoefMax = %f", &
      &  d=(/TempBulkCoefMax/))
    call MessageNotify( "M", module_name, "QmixBulkCoefMax = %f", &
      &  d=(/QmixBulkCoefMax/))

    call HistoryAutoAddVariable(      &
      & varname='PTempSfcFlux',       &
      & dims=(/'x','y','t'/),         &
      & longname='surface potential temperature flux (heat flux divided by density and specific heat)', &
      & units='K.m.s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(      &
      & varname='ExnerSfcFlux',       &
      & dims=(/'x','y','t'/),         &
      & longname='surface exner function flux (heat flux divided by density and specific heat)', &
      & units='s-1',             &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfcFlux',    &
      & dims=(/'x','y','t'/),     &
      & longname='surface flux of x-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfcFlux',    &
      & dims=(/'x','y','t'/),     &
      & longname='surface flux of y-component of velocity (momentum flux divided by density)', &
      & units='m2.s-2',           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfcFlux', & 
        & dims=(/'x','y','t'/),     &
        & longname='surface flux of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio (mass flux divided by density)',  &
        & units='m.s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(  &
      & varname='PTempSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='potential temperature tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='exner function tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelXSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='x-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='VelYSfc',         &
      & dims=(/'x','y','z','t'/), &
      & longname='y-component velocity tendency by surface flux', &
      & units='m.s-2',            &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfc', & 
        & dims=(/'x','y','z','t'/),     &
        & longname=trim(SpcWetSymbol(l))//' mixing ratio tendency by surface flux',  &
        & units='s-1',    &
        & xtype='float')
    end do


    call HistoryAutoAddVariable(       &
      & varname='SfcHeatFlux',         &
      & dims=(/'x','y','t'/),          &
      & longname='surface heat flux',  &
      & units='W.m-2',                 &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcXMomFlux',                        &
      & dims=(/'x','y','t'/),                         &
      & longname='surface x-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    call HistoryAutoAddVariable(                      &
      & varname='SfcYMomFlux',                        &
      & dims=(/'x','y','t'/),                         &
      & longname='surface y-component momentum flux', &
      & units='kg.m-2.s-1',                           &
      & xtype='float')

    do l = 1, ncmax
      call HistoryAutoAddVariable(                               &
        & varname='D'//trim(SpcWetSymbol(l))//'DtSfcMassFlux',         &
        & dims=(/'x','y','t'/),                              &
        & longname=trim(SpcWetSymbol(l))//' surface mass flux',  &
        & units='kg.m-2.s-1',                                    &
        & xtype='float')
    end do

  end subroutine Surfaceflux_Bulk_init
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_Bulk_forcing( &
    &   pyz_VelX, xqz_VelY, xyz_PTemp, xyz_Exner, xyzf_QMix, &
    &   pyz_DVelXDt, xqz_DVelYDt, xyz_DPTempDt, xyz_DExnerDt, xyzf_DQMixDt &
    & )
    ! 
    ! 下部境界からのフラックスによる温度の変化率を,
    ! バルク方法に基づいて計算する.
    !

    ! 暗黙の型宣言禁止
    ! Implicit none
    !
    implicit none

    ! 変数
    ! variables
    !
    real(DP), intent(in)   :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速
                              ! X-component velocity
    real(DP), intent(in)   :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速
                              ! Y-component velocity
    real(DP), intent(in)   :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
                              ! 温位
                              ! Potential temperature
    real(DP), intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                              ! 圧力関数
                              ! Exner function
    real(DP), intent(in)   :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! 混合比
                              ! Mixing ration
    real(DP), intent(inout):: pyz_DVelXDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! 風速時間変化率
                              ! X-component velocity tendency
    real(DP), intent(inout):: xqz_DVelYDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! 風速時間変化率
                              ! Y-component velocity tendency
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! 温位時間変化率
                              ! Potential tempreture tendency
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)
                              ! 圧力関数時間変化率
                              ! Exner function tendency
    real(DP), intent(inout):: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! 混合比時間変化率
                              ! Mixing ratio tendency

    ! 作業変数
    ! Work variables
    real(DP) :: xy_SurfBulkRiNum(imin:imax,jmin:jmax)
                              ! バルクリチャードソン数
                              ! Bulk Richardson number
    real(DP) :: xy_SurfRoughLength(imin:imax,jmin:jmax)
                              ! 祖度長さ
                              ! Roughness length
    real(DP) :: xy_SurfVelBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(運動量)
                              ! Bulk coefficient for momentum
    real(DP) :: xy_SurfTempBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(熱)
                              ! Bulk coefficient for heat
    real(DP) :: xy_SurfQmixBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(混合比)
                              ! Bulk coefficient for mixing ratio
    real(DP) :: py_VelXflux (imin:imax,jmin:jmax)
                              ! x 方向速度フラックス
                              ! velocity flux in x
    real(DP) :: xq_VelYflux (imin:imax,jmin:jmax)
                              ! y 方向速度フラックス
                              ! celocity flux in y
    real(DP) :: xy_DPTempDtFlux(imin:imax,jmin:jmax)
                              ! 温位フラックス
                              ! potential temperature flux
    real(DP) :: xy_ExnerFlux(imin:imax,jmin:jmax)
                              !
                              !
    real(DP) :: xyf_QMixFlux(imin:imax,jmin:jmax,ncmax)
                              ! 凝結成分混合比フラックス
                              ! Mixing ratio flux
    real(DP) :: xyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速 (xyz 格子)
                              ! X-component velocity (xyz grid)
    real(DP) :: xyz_VelY(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速 (xyz 格子)
                              ! Y-component velocity (xyz grid)
    real(DP) :: xyz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速絶対値 (xyz 格子)
                              ! Absolute value of horizontal velocity (xyz grid)
    real(DP) :: pyz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速絶対値 (pyz 格子)
                              ! Absolute value of horizontal velocity (pyz grid)
    real(DP) :: xqz_AbsVel(imin:imax,jmin:jmax,kmin:kmax)
                              ! 水平風速絶対値 (xqz 格子)
                              ! Absolute value of horizontal velocity (xqz grid)
    real(DP) :: xyz_PTempAll(imin:imax,jmin:jmax,kmin:kmax)
                              ! 温位(基本場 + 擾乱)
                              ! Total value of potential temperature
    real(DP) :: xyz_ExnerAll (imin:imax,jmin:jmax,kmin:kmax)
                              ! 圧力関数(基本場 + 擾乱)
                              ! Total value of exner function
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! 混合比(基本場 + 擾乱)
                              ! Total value of mixing ratios
    real(DP) :: xy_DPTempDtBulk(imin:imax,jmin:jmax)
                              ! 時間変化(温位)
                              ! potential temperature tendency by surface flux
    real(DP) :: xy_DExnerDtBulk(imin:imax,jmin:jmax)
                              ! 時間変化(圧力関数)
                              ! Exner function tendency by surface flux
    real(DP) :: xyf_DQMixDtBulk(imin:imax,jmin:jmax, ncmax)
                              ! 時間変化(混合比)
                              ! Mixing ratio tendency by surface flux
    real(DP) :: py_DVelXDtBulk (imin:imax,jmin:jmax)
                              ! 時間変化(U)
                              ! x-component velocity tendency by surface flux
    real(DP) :: xq_DVelYDtBulk (imin:imax,jmin:jmax)
                              ! 時間変化(V)
                              ! y-component velocity tendency by surface flux
    real(DP) :: xyz_DPTempDtBulk(imin:imax,jmin:jmax,kmin:kmax)
                              ! 時間変化(温位)
                              ! potential temperature tendency by surface flux
    real(DP) :: xyz_DExnerDtBulk(imin:imax,jmin:jmax,kmin:kmax)
                              ! 時間変化(圧力関数)
                              ! Exner function tendency by surface flux
    real(DP) :: xyzf_DQMixDtBulk(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                              ! 時間変化(混合比)
                              ! Mixing ratio tendency by surface flux
    real(DP) :: pyz_DVelXDtBulk (imin:imax,jmin:jmax,kmin:kmax)
                              ! 時間変化(U)
                              ! x-component velocity tendency by surface flux
    real(DP) :: xqz_DVelYDtBulk (imin:imax,jmin:jmax,kmin:kmax)
                              ! 時間変化(V)
                              ! y-component velocity tendency by surface flux
    real(DP) :: ExnerBZSfc    ! 地表面圧力関数 
                              ! Basic state Exner function at the surface
    real(DP) :: xy_PressSfc(imin:imax,jmin:jmax)
                              ! Total pressure at the surface
    integer  :: kz            ! 配列添字
                              ! Arrzy index
    integer  :: s             ! 組成方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in dimension of constituen                         

    ! 初期化
    ! Initialization
    ! 
    kz = 1

    ! 祖度長さの指定
    ! Specify surface length
    xy_SurfRoughLength = SfcRoughLength

    ! 全量の計算
    ! Calculate total value of thermodynamic variables
    ! 
    xyz_PTempAll  = xyz_PTemp + xyz_PTempBZ
    xyz_ExnerAll  = xyz_Exner + xyz_ExnerBZ
    xyzf_QMixAll  = xyzf_QMix + xyzf_QMixBZ

    ! Perturbation component of Exner function at the surface is assumed 
    ! to be same as that at the lowest layer. (YOT, 2011/09/03)
    ! 
    ExnerBZSfc    = (PressSfc / PressBasis) ** (GasRDry / CpDry)
    xy_PressSfc   = (PressBasis * xyz_ExnerAll(:,:,kz))**(CpDry / GasRDry)

    ! xyz 格子点の速度の計算
    ! Calculate velocities at xyz grid points
    !
    xyz_VelX = xyz_avr_pyz(pyz_VelX)
    xyz_VelY = xyz_avr_xqz(xqz_VelY)

    ! 水平風速の絶対値の計算
    ! Calculate of absoluto horizontal velocities
    ! 
    xyz_AbsVel = SQRT( xyz_VelX**2 + xyz_VelY**2 + Vel0**2 )
    pyz_AbsVel = pyz_avr_xyz(xyz_AbsVel)
    xqz_AbsVel = xqz_avr_xyz(xyz_AbsVel)

    ! バルク $ R_i $ 数算出
    ! Calculate bulk $ R_i $
    !
    xy_SurfBulkRiNum =                                 &
      &   Grav / ( xyz_PTempAll(:,:,1) )               & 
      &   * ( xyz_PTempAll(:,:,1)                      &
      &     -TempSfc / (ExnerBZSfc + xyz_Exner(:,:,1)))&
      &   / max( xyz_AbsVel(:,:,1), VelMinForRi )**2   &
      &   * z_dz(1) * 0.5d0

    ! バルク係数の計算
    ! Bulk coefficients are calculated.
    ! 
    call BulkCoef( &
      & xy_SurfBulkRiNum,    & ! (in)
      & xy_SurfRoughLength,  & ! (in)
      & xy_SurfVelBulkCoef,  & ! (out)
      & xy_SurfTempBulkCoef, & ! (out)
      & xy_SurfQmixBulkCoef  & ! (out)
      & )
    ! フラックスの計算
    ! Surface fluxes are calculated. 
    !
    xy_DPTempDtFlux = - xy_SurfTempBulkCoef * xyz_AbsVel(:,:,kz)  &
      & * ( xyz_PTempAll(:,:,kz) - TempSfc / ( ExnerBZSfc + xyz_Exner(:,:,kz) ) )

    xyf_QMixFlux = 0.0d0
    do s = 1, CondNum
      xyf_QMixFlux(:,:,IdxCG(s)) =                                 &
        &     - xy_SurfQmixBulkCoef * xyz_AbsVel(:,:,kz)           &
        &       * (                                                &
        &            xyzf_QMixAll(:,:,kz,s)                        &
        &          - SvapPress( SpcWetID(IdxCC(s)), TempSfc )      &
        &             / ( xy_PressSfc )                            &
        &             * (MolWtWet(IdxCG(s)) / MolWtDry)            &
        &         )
    end do
    !
    xy_ExnerFlux = xy_DExnerDt_xy_xyf(xy_DPTempDtFlux, xyf_QMixFlux, kz)
    !
    py_VelXFlux = - xy_SurfVelBulkCoef * pyz_AbsVel(:,:,kz) * pyz_VelX(:,:,kz)
    !
    xq_VelYFlux = - xy_SurfVelBulkCoef * xqz_AbsVel(:,:,kz) * xqz_VelY(:,:,kz)


    ! フラックスの下限値を設定
    ! Set lower limit of surface fluxes
    ! 
    xy_DPTempDtFlux = max( DPTempDtFluxMin, xy_DPTempDtFlux )
    xy_ExnerFlux = max( ExnerFluxMin, xy_ExnerFlux )
    xyf_QMixFlux = max( QMixFluxMin, xyf_QMixFlux )


    ! 地表フラックスによる時間変化を計算
    ! Tendencies by surface fluxes (convergences of fluxes) are calculated.
    !
    xy_DPTempDtBulk = - ( 0.0d0 - xy_DPTempDtFlux ) / z_dz(kz)
    !
    xy_DExnerDtBulk = - ( 0.0d0 - xy_ExnerFlux ) / z_dz(kz)
    !
    xyf_DQMixDtBulk = - ( 0.0d0 - xyf_QMixFlux ) / z_dz(kz)
    !
    py_DVelXDtBulk = - ( 0.0d0 - py_VelXFlux ) / z_dz(kz)
    !
    xq_DVelYDtBulk = - ( 0.0d0 - xq_VelYFlux ) / z_dz(kz)

    
    ! 仮引数配列へ格納
    ! Add tendency by surface flux convergence
    !
    xyz_DPTempDt(:,:,kz) = xyz_DPTempDt(:,:,kz) + xy_DPTempDtBulk
    xyz_DExnerDt(:,:,kz) = xyz_DExnerDt(:,:,kz) + xy_DExnerDtBulk
    do s = 1, ncmax
      xyzf_DQMixDt(:,:,kz,s) = xyzf_DQMixDt(:,:,kz,s) + xyf_DQMixDtBulk(:,:,s)
    end do
    pyz_DVelXDt (:,:,kz) = pyz_DVelXDt (:,:,kz) + py_DVelXDtBulk
    xqz_DVelYDt (:,:,kz) = xqz_DVelYDt (:,:,kz) + xq_DVelYDtBulk

    ! 出力
    ! Output
    !
    xyz_DPTempDtBulk = 0.0d0
    xyz_DExnerDtBulk = 0.0d0
    xyzf_DQMixDtBulk = 0.0d0
    pyz_DVelXDtBulk  = 0.0d0
    xqz_DVelYDtBulk  = 0.0d0

    xyz_DPTempDtBulk(:,:,kz) = xy_DPTempDtBulk
    xyz_DExnerDtBulk(:,:,kz) = xy_DExnerDtBulk
    do s = 1, ncmax
      xyzf_DQMixDtBulk(:,:,kz,s) = xyf_DQMixDtBulk(:,:,s)
    end do
    pyz_DVelXDtBulk (:,:,kz) = py_DVelXDtBulk
    xqz_DVelYDtBulk (:,:,kz) = xq_DVelYDtBulk
    !
    call HistoryAutoPut(TimeN, 'PTempSfc', xyz_DPTempDtBulk(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerSfc', xyz_DExnerDtBulk(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelXSfc',  pyz_DVelXDtBulk (1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelYSfc',  xqz_DVelYDtBulk (1:nx,1:ny,1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfc', &
        & xyzf_DQMixDtBulk(1:nx,1:ny,1:nz,s))
    end do


    call HistoryAutoPut(TimeN, 'PTempSfcFlux', xy_DPTempDtFlux(1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'ExnerSfcFlux', xy_ExnerFlux(1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'VelXSfcFlux',  py_VelXFlux (1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'VelYSfcFlux',  xq_VelYFlux (1:nx,1:ny))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfcFlux', &
        & xyf_QMixFlux(1:nx,1:ny,s))
    end do

    call HistoryAutoPut(TimeN, 'SfcHeatFlux', &
      & CpDry * xyz_DensBZ(1:nx,1:ny,1) * xy_DPTempDtFlux(1:nx,1:ny) &
      &   * ( ExnerBZSfc + xyz_Exner(1:nx,1:ny,1) ) )
    call HistoryAutoPut(TimeN, 'SfcXMomFlux', &
      & xyz_DensBZ(1:nx,1:ny,1) * py_VelXFlux (1:nx,1:ny))
    call HistoryAutoPut(TimeN, 'SfcYMomFlux', &
      & xyz_DensBZ(1:nx,1:ny,1) * xq_VelYFlux (1:nx,1:ny))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtSfcMassFlux', &
        & xyz_DensBZ(1:nx,1:ny,1) * xyf_QMixFlux(1:nx,1:ny,s))
    end do

  end subroutine Surfaceflux_Bulk_forcing
!!!------------------------------------------------------------------------!!!
  subroutine BulkCoef(     & 
    & xy_SurfBulkRiNum,    & ! (in)
    & xy_SurfRoughLength,  & ! (in)
    & xy_SurfVelBulkCoef,  & ! (out)
    & xy_SurfTempBulkCoef, & ! (out)
    & xy_SurfQmixBulkCoef  & ! (out)
    & )


    implicit none

    real(DP),intent(in) :: xy_SurfBulkRiNum(imin:imax,jmin:jmax)
                              ! バルクリチャードソン数
                              ! Bulk Richardson number
    real(DP),intent(in) :: xy_SurfRoughLength(imin:imax,jmin:jmax)
                              ! 祖度長さ
                              ! Roughness length
    real(DP),intent(out) :: xy_SurfVelBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(運動量)
                              ! Bulk coefficient for momentum
    real(DP),intent(out) :: xy_SurfTempBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(熱)
                              ! Bulk coefficient for heat
    real(DP),intent(out) :: xy_SurfQmixBulkCoef(imin:imax,jmin:jmax)
                              ! バルク係数(混合比)
                              ! Bulk coefficient for mixing ratio

    ! 作業変数
    ! Work variables
    !
    real(DP) :: xy_SurfBulkCoefInNeutCond(imin:imax, jmin:jmax)

    integer:: i               ! 経度方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in longitude
    integer:: j               ! 緯度方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in latitude

    ! 実行文 ; Executable statement
    !
    if ( FlagConstBulkCoef ) then

      ! Use of constant bulk coefficient
      !

      xy_SurfVelBulkCoef  = ConstBulkCoef
      xy_SurfTempBulkCoef = ConstBulkCoef
      xy_SurfQmixBulkCoef = ConstBulkCoef

    else

      ! Parameterization by Louis et al. (1981)
      !
      ! 中立バルク係数の計算
      ! Calculate bulk coefficient in neutral condition
      !
      xy_SurfBulkCoefInNeutCond  = &
        & ( FKarm &
        & / log ( z_dz(1) / xy_SurfRoughLength ) )**2

      if ( FlagUseOfBulkCoefInNeutralCond ) then

        ! 中立条件でのバルク係数の設定
        ! Set bulk coefficient in neutral condition
        !

        xy_SurfVelBulkCoef  = xy_SurfBulkCoefInNeutCond
        xy_SurfTempBulkCoef = xy_SurfBulkCoefInNeutCond

        xy_SurfQmixBulkCoef = xy_SurfTempBulkCoef

      else

        do j = jmin, jmax
          do i = imin, imax

            if ( xy_SurfBulkRiNum(i,j) > 0.0_DP ) then 

              xy_SurfVelBulkCoef(i,j) =                                        &
                &   xy_SurfBulkCoefInNeutCond(i,j)                             &
                &   / (   1.0_DP                                               &
                &       + 4.7_DP * xy_SurfBulkRiNum(i,j)                      &
                &     )**2

              xy_SurfTempBulkCoef(i,j) = xy_SurfVelBulkCoef(i,j) 

              xy_SurfQmixBulkCoef(i,j) = xy_SurfTempBulkCoef(i,j)

            else

              xy_SurfVelBulkCoef(i,j) =                                       &
                &   xy_SurfBulkCoefInNeutCond(i,j)                            &
                &   * (   1.0_DP                                              &
                &       - 9.4_DP * xy_SurfBulkRiNum(i,j)                     &
                &           / (   1.0_DP                                      &
                &             + 69.56_DP * xy_SurfBulkCoefInNeutCond(i,j)      &
                &               * sqrt( - z_dz(1)  &
                &                         / xy_SurfRoughLength(i,j)           &
                &                         * xy_SurfBulkRiNum(i,j)             &
                &                     )                                       &
                &           )                                                 &
                &     )

              xy_SurfTempBulkCoef(i,j) =                                      &
                &   xy_SurfBulkCoefInNeutCond(i,j)                            &
                &   * (   1.0_DP                                              &
                &       - 9.4_DP * xy_SurfBulkRiNum(i,j)                     &
                &         / (   1.0_DP                                        &
                &             + 49.82_DP * xy_SurfBulkCoefInNeutCond(i,j)      &
                &               * sqrt( - z_dz(1) &
                &                         / xy_SurfRoughLength(i,j)           &
                &                         * xy_SurfBulkRiNum(i,j)             &
                &                     )                                       &
                &           )                                                 &
                &     )

              xy_SurfQmixBulkCoef(i,j) = xy_SurfTempBulkCoef(i,j)
    
           end if
         end do
        end do
      end if
    end if

    ! 最大/最小 判定
    ! Measure maximum/minimum
    !
    do i = imin, imax
      do j = jmin, jmax

        xy_SurfVelBulkCoef(i,j)  = &
          & max( min( xy_SurfVelBulkCoef(i,j), VelBulkCoefMax ), &
          &      VelBulkCoefMin )

        xy_SurfTempBulkCoef(i,j) = &
          & max( min( xy_SurfTempBulkCoef(i,j), TempBulkCoefMax ), &
          &      TempBulkCoefMin )

        xy_SurfQmixBulkCoef(i,j) = &
          & max( min( xy_SurfQmixBulkCoef(i,j), QmixBulkCoefMax ), &
          &      QmixBulkCoefMin )

      end do
    end do

  end subroutine BulkCoef
!!!------------------------------------------------------------------------!!!  
end module Surfaceflux_bulk_l1979
