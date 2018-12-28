!= 簡単放射モジュール: とある高度領域を一様冷却・加熱する
!
! Authors::   YAMASHITA Tatsuya, SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_heatbalance.f90,v 1.13 2014/07/08 01:04:18 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_HeatBalance
  !
  ! 簡単放射モジュール: とある高度領域を一様冷却・加熱する
  !

  !モジュール読み込み
  use dc_types, only: DP, STRING
  
  !暗黙の型宣言禁止
  implicit none

  !private 属性
  private

  !変数定義
  real(DP),save   :: RadCoolRate = 0.0d0    !一様放射加熱率 [K/day]
  integer, save   :: IdxHeatUp   = 0        !加熱域上限の鉛直座標に対応する整数値
  integer, save   :: IdxHeatDown = 0        !冷却域上限の鉛直座標に対応する整数値
  integer, save   :: IdxCoolUp   = 0        !加熱域上限の鉛直座標に対応する整数値
  integer, save   :: IdxCoolDown = 0        !冷却域上限の鉛直座標に対応する整数値
  real(DP),save, allocatable :: xyz_RadHeightHeat(:,:,:)  !放射加熱が存在する領域
  real(DP),save, allocatable :: xyz_RadHeightCool(:,:,:)  !放射加熱が存在する領域

  real(DP), save  :: FactorDExnerDtRad = 1.0d0

  character(STRING), parameter:: module_name = 'Radiation_Heatbalance'

  !関数を public にする. 
  public Radiation_heatbalance_init
  public Radiation_heatbalance_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_heatbalance_init
    !
    !NAMELIST から冷却率, 冷却する高度領域, 加熱する高度領域, を設定する. 
    !加熱率は毎ステップ計算するので, 初期化ルーチンの中では決めない. 
    !

    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use namelist_util,     only : namelist_filename
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax        !z 方向の配列の上限
    use axesset,           only : z_Z         !Z 座標軸(スカラー格子点)
    
    !暗黙の型宣言禁止
    implicit none

    !入力変数
    real(8) :: HeightHeatUp   = 0.0d0  !加熱領域の上端の高度
    real(8) :: HeightHeatDown = 0.0d0  !加熱領域の下端の高度
    real(8) :: HeightCoolUp   = 0.0d0  !冷却領域の上端の高度
    real(8) :: HeightCoolDown = 0.0d0  !冷却領域の下端の高度
    integer :: k                       !ループ変数
    integer :: unit

    ! NAMELIST から情報を取得
    NAMELIST /radiation_heatbalance_nml/ &
      & RadCoolRate, HeightHeatUp, HeightHeatDown, HeightCoolUp, HeightCoolDown, &
      & FactorDExnerDtRad

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_heatbalance_nml)
    close(unit)

    do k = kmin, kmax-1
      if ( z_Z(k) <= HeightHeatUp .AND. HeightHeatUp < z_Z(k+1) ) then 
        IdxHeatUp = k
      elseif ( z_Z(k) <= HeightHeatDown .AND. HeightHeatDown < z_Z(k+1) ) then 
        IdxHeatDown = k
      elseif ( z_Z(k) <= HeightCoolUp .AND. HeightCoolUp < z_Z(k+1) ) then 
        IdxCoolUp = k
      elseif ( z_Z(k) <= HeightCoolDown .AND. HeightCoolDown < z_Z(k+1) ) then 
        IdxCoolDown = k
      end if
    end do

    allocate( xyz_RadHeightHeat(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_RadHeightCool(imin:imax, jmin:jmax, kmin:kmax) )

    xyz_RadHeightHeat = 1.0d0
    xyz_RadHeightHeat(:,:,IdxHeatDown:IdxHeatUp) = 1.0d0
    xyz_RadHeightCool = 0.0d0
    xyz_RadHeightCool(:,:,IdxCoolDown:IdxCoolUp) = 1.0d0

    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "RadCoolRate = %f", d=(/RadCoolRate/))
    call MessageNotify( "M", &
      & module_name, "HeightHeatUp = %f", d=(/HeightHeatUp/))
    call MessageNotify( "M", &
      & module_name, "HeightHeatDown= %f", d=(/HeightHeatDown/))
    call MessageNotify( "M", &
      & module_name, "HeightCoolUp = %f", d=(/HeightCoolUp/))
    call MessageNotify( "M", &
      & module_name, "HeightCoolDown= %f", d=(/HeightCoolDown/))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtRad= %f", d=(/ FactorDExnerDtRad /))

    ! ヒストリデータ定義
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1"',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of exner function', &
      & units='s-1"',    &
      & xtype='float')

  end subroutine Radiation_heatbalance_init

!!!------------------------------------------------------------------------!!!
  subroutine radiation_heatbalance_forcing(xyz_Exner, xyz_PTemp, xyz_DPTempDt, xyz_DExnerDt)
    !
    ! 温位の放射強制項. 
    ! 地表面から RadHeight で指定された高度までの間で空間的に一様な加熱, 
    ! RadHeight から RadHeight2 までの間で空間的な一様な冷却を与える. 
    ! 放射強制全体として加熱と冷却がバランスするように時々刻々一様加熱率
    ! および一様冷却率を変化させる. 
    ! 冷却の振幅を与え, 加熱の振幅を調節する. 

    ! モジュール呼び出し
    !
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax,     & !z 方向の配列の上限
      &                           nx, ny, nz
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use axesset,           only : z_dz  
    
    use constants,         only : DayTime,  & ! 1 日の長さ [s]
      &                           PressBasis, & !温位の基準圧力
      &                           CpDry,      & !定圧比熱
      &                           CvDry,      & !定積比熱
      &                           GasRDry       !気体定数
    use basicset,          only : xyz_ExnerBZ,  &!エクスナー関数の基本場
      &                           xyz_PTempBZ    !温位の基本場


    !暗黙の型宣言を禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: xyz_Exner(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)  :: xyz_PTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DPTempDt0(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DExnerDt0(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Rad(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_RadPI(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DensSum(imin:imax, jmin:jmax, kmin:kmax) 
                                        !密度(基本場成分と擾乱成分の和)
    real(DP)              :: HeatSum    !ある時刻での単位時間当たりの加熱量
    real(DP)              :: CoolSum    !ある時刻での単位時間当たりの冷却量
    real(DP)              :: RadHeatRate
    integer               :: i, j, k

    ! 初期化
    !
    HeatSum = 0.0d0
    CoolSum = 0.0d0
    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt

    ! 密度の算出
    xyz_DensSum = PressBasis                            &
      & * (xyz_ExnerBZ + xyz_Exner)**( CvDry /GasRDry ) &
      & / (GasRDry * (xyz_PTempBZ + xyz_PTemp ) )

    ! ある時刻での (加熱量/加熱率) を計算
    do k = IdxHeatDown, IdxHeatUp
      do j = 1, ny
        do i = 1, nx
          HeatSum = HeatSum + z_dz(k) * CpDry * xyz_DensSum(i,j,k)
        end do
      end do
    end do
    
    ! ある時刻での (冷却量/冷却率) を計算
    do k = IdxCoolDown, IdxCoolUp
      do j = 1, ny
        do i = 1, nx
          CoolSum = CoolSum + z_dz(k) * CpDry * xyz_DensSum(i,j,k)
        end do
      end do
    end do

    ! 加熱率の算出
    ! RadCoolRate が負値であることに注意. 
    RadHeatRate = - RadCoolRate * CoolSum / HeatSum
    
    xyz_Rad = &
      &   xyz_RadHeightHeat * RadHeatRate / ( xyz_ExnerBZ  * DayTime ) &
      & + xyz_RadHeightCool * RadCoolRate / ( xyz_ExnerBZ  * DayTime )

    xyz_DPTempDt = xyz_DPTempDt0 + xyz_Rad

    ! 圧力変化
    !
    xyz_RadPI = xyz_DExnerDt_xyz( xyz_Rad ) * FactorDExnerDtRad

    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_Rad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_RadPI(1:nx,1:ny,1:nz))
   
  end subroutine radiation_heatbalance_forcing
  
end module Radiation_HeatBalance
