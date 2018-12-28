!= Module initialdata_yamasaki1983
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_yamasaki1983.f90,v 1.8 2014/03/04 04:49:40 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_yamasaki1983
  !
  ! Yamasaki (1983) の基本場・擾乱場を設定する

   
  !暗黙の型宣言禁止
  implicit none

  !初期化だけ公開
  public initialdata_yamasaki1983_basic

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_yamasaki1983_basic( z_Temp, z_Press, zf_MolFr )
    !
    !== 概要
    !  * Yamasaki, 1983 の温度と相対湿度の観測値で基本場を作成する
    !    * 温度の基本場
    !      * 観測データをNetCDFファイル化したものから値を読み込む
    !        * 読み込んだ値を線形補間して温度の基本場を作成する
    !    * 湿度の基本場
    !      * 相対湿度は subroutine HUM で作成済み
    !        * ここでは相対湿度をモル比に変換して湿度の基本場を作成する 
    !    * 気圧の基本場
    !      * 湿潤成分と乾燥成分の分子量差を考慮した静水圧平衡から計算する
    !

    !モジュール読み込み
    use dc_types,   only: STRING, DP
    use dc_message, only: MessageNotify
    use gridset,    only: kmin, kmax,  &!配列サイズ (Z 方向)
      &                   nz,          &!物理領域のサイズ (Z 方向)
      &                   ncmax         !凝縮成分の数
    use axesset,    only: z_Z,         &!スカラー格子点での高度
      &                   dz            !格子間隔
    use constants,  only: GasRDry,       &!乾燥成分の定圧比熱
      &                   TempSfc,       &!地表面温度
      &                   PressSfc,      &!地表面圧力
      &                   PressBasis,    &!地表面圧力
      &                   Grav            !重力
    use chemcalc,   only: SvapPress       !飽和蒸気圧
    !  use composition

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer, parameter :: mmin = 1
    integer, parameter :: mmax1 = 37
    integer, parameter :: mmax2 = 31
    integer :: ID = 6    ! 水の物質番号

    real(DP), intent(out) :: z_Press(kmin:kmax)!圧力
    real(DP), intent(out) :: z_Temp(kmin:kmax) !温度
    real(DP), intent(out) :: zf_MolFr(kmin:kmax, 1:ncmax) !モル比
    real(DP)              :: Ob_Alt1(mmin:mmax1)    !観測データに対応する高度
    real(DP)              :: Ob_Alt2(mmin:mmax2)    !観測データに対応する高度
    real(DP)              :: Ob_Temp(mmin:mmax1)    !NetCDF から読み込んだ気温の観測値(1次元)
    real(DP)              :: Ob_Hum(mmin:mmax2)     !相対湿度の観測値(1次元)

    real(DP)              :: z_Hum(kmin:kmax)        !相対湿度の基本場の値(%)
    integer               :: k, m

    ! 初期化
    !
    z_Press  = 0.0d0
    z_Temp   = 0.0d0
    zf_MolFr = 0.0d0

    ! 確認. 地表面温度・圧力を指定するのは別モジュールなので. 
    !
    if (     PressBasis /= 1010.0d2  &
      & .OR. PressSfc   /= 1010.0d2  &
      & .OR. TempSfc    /= 302.0d0  ) then 
      
      call MessageNotify( "E", "initaldata_yamasaki1983_init", &
        & "Constants are wrong. please PressSfc = 1.01d5, TempSfc = 302.0d0")
    end if

    ! 観測データ. 
    !
    Ob_Alt1 = (/ &
      & 0.05d3, 0.16d3, 0.29d3, 0.44d3, 0.61d3,      &
      & 0.80d3, 1.02d3, 1.28d3, 1.58d3, 1.95d3,      &
      & 2.37d3, 2.82d3, 3.33d3, 3.90d3, 4.50d3,      &
      & 5.10d3, 5.70d3, 6.30d3, 6.90d3, 7.50d3,      &
      & 8.10d3, 8.70d3, 9.30d3, 9.90d3, 10.50d3,     &
      & 11.10d3, 11.70d3, 12.35d3, 13.05d3, 13.80d3, &
      & 14.65d3, 15.50d3, 16.55d3, 17.70d3, 18.95d3, &
      & 20.30d3, 21.80d3                             &
      &/)

    Ob_Temp = (/ &
      & 299.60d0, 298.72d0, 297.68d0, 296.48d0, 295.13d0, &
      & 293.90d0, 292.47d0, 290.90d0, 289.40d0, 287.55d0, &
      & 285.45d0, 283.13d0, 280.38d0, 277.25d0, 273.95d0, &
      & 270.65d0, 267.35d0, 263.90d0, 260.30d0, 256.55d0, &
      & 252.50d0, 248.15d0, 243.50d0, 238.70d0, 233.90d0, &
      & 229.10d0, 224.30d0, 219.25d0, 213.85d0, 208.50d0, &
      & 203.70d0, 200.55d0, 199.60d0, 201.40d0, 205.15d0, &
      & 209.20d0, 212.90d0                                &
      & /)

    Ob_Alt2 = (/ &
      & 0.00d3, 0.60d3, 1.20d3, 1.80d3, 2.40d3, &
      & 3.00d3, 3.60d3, 4.20d3, 4.80d3, 5.40d3, &
      & 6.00d3, 6.60d3, 7.20d3, 7.80d3, 8.40d3, &
      & 9.00d3, 9.60d3, 10.2d3, 10.8d3, 11.4d3, &
      & 12.0d3, 12.7d3, 13.4d3, 14.2d3, 15.1d3, &
      & 16.0d3, 17.1d3, 18.3d3, 19.6d3, 21.0d3, 22.6d3 &
      & /)

    Ob_Hum = (/ &
      & 83.0d-2, 91.0d-2, 95.0d-2, 95.0d-2, 91.0d-2, &
      & 85.0d-2, 80.0d-2, 75.0d-2, 71.0d-2, 68.0d-2, &
      & 66.0d-2, 65.0d-2, 64.0d-2, 63.0d-2, 62.0d-2, &
      & 61.0d-2, 61.0d-2, 62.0d-2, 62.0d-2, 63.0d-2, &
      & 63.0d-2, 64.0d-2, 67.0d-2, 74.0d-2, 75.0d-2, &
      & 46.0d-2, 26.0d-2,  7.0d-2,  1.0d-2,  0.5d-2,  0.2d-2 &
      & /)

    ! 線形補完
    !
    ! Ob_Hum を線形補間して z_Hum を作成する
    ! Ob_Alt(m) < z_Z(k) < Ob_Alt(m+1),
    ! Ob_Alt = z_Z(k),
    ! Ob_Alt < z_Z(k)
    ! の3つで場合分け
    !
    do k = 1, nz
      do m = mmin, mmax1 - 1
        ! z_Z が Ob_altitude のある区間に挟まれるとき,
        ! その区間の Obaltitude(m), Obaltitude(m+1) を結ぶ
        ! 直線で Ob_TempZ を線形補間する
        if (Ob_Alt1(m) /=  z_Z(k) .AND. z_Z(k) > Ob_Alt1(m) .AND. z_Z(k) < Ob_Alt1(m+1)) then

          z_Temp(k) = Ob_Temp(m) &
            & + ( ( Ob_Temp(m+1) - Ob_Temp(m) ) / ( Ob_Alt1(m+1) - Ob_Alt1(m) ) ) &
            &   * ( z_Z(k) - Ob_Alt1(m) )
          
        else if (Ob_Alt1(m) == z_Z(k)) then
          z_Temp(k)  = Ob_Temp(m)

        ! z_Z(k) > Ob_altitude では観測データが無いので等温大気にする          
        else if (Ob_Alt1(m) < z_Z(k)) then
          z_Temp(k)  = z_Temp(k-1)
          
        end if
      end do
    end do

    do k = 1, nz
      do m = mmin, mmax2 - 1 
        if (Ob_Alt2(m) /=  z_Z(k) .AND. z_Z(k) > Ob_Alt2(m) .AND. z_Z(k) < Ob_Alt2(m+1)) then

          z_Hum(k) = Ob_Hum(m) &
            & + ( ( Ob_Hum(m+1) - Ob_Hum(m) ) / ( Ob_Alt2(m+1) - Ob_Alt2(m) ) ) &
            &   * (z_Z(k) - Ob_Alt2(m))
          
        else if (Ob_Alt2(m) == z_Z(k)) then
          z_Hum(k) = Ob_Hum(m)

        ! z_Z(k) > Ob_altitude では観測データが無いので相対湿度一定にする 
        else if (Ob_Alt2(m) < z_Z(k)) then
          z_Hum(k) = z_Hum(k-1)

        end if
      end do
    end do

    ! 初期化
    z_Press = 1.0d-60

    ! 地表面気圧と温度から大気最下層の気圧を求める
    ! 静水圧の式 dP/dz = - \rho * g を使用
    z_Press(1) = PressSfc - (Grav * PressSfc * dz * 5.0d-1) / (GasRDry * TempSfc)

    ! 大気最下層のモル比を計算する
    zf_MolFr(1,1) = SvapPress( ID, z_Temp(1)) * z_Hum(1) / z_Press(1)
    
    ! 湿度の基本場の計算
    do k = 2, nz
      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * dz) &
        & / ( GasRDry * z_Temp(k-1) )
      
      zf_MolFr(k,1) = SvapPress( ID, z_Temp(k)) * z_Hum(k) / z_Press(k) 
    end do
    
  end subroutine INITIALDATA_YAMASAKI1983_basic
  
end module Initialdata_yamasaki1983
