!= Module EnergyMonit
!
! Authors::   ODAKA Masatsugu 
! Version::   $Id: energymonit.f90,v 1.4 2014/07/11 08:02:38 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module EnergyMonit
  !
  ! 保存量の計算と出力を行うためのモジュール
  !

  !モジュール読み込み
  use dc_types,   only: DP, STRING
  
  ! 変数定義
  !
  real(DP), private                     :: MassTotalBZ   !全質量(基本場)
  character(STRING), private, parameter :: module_name = 'energymonit'
                                                         ! モジュールの名称.
                                                         ! Module name
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine EnergyMonit_init
    !
    !出力変数の設定を行う. 
    !
    ! * DensDev2 という変数名はダサイが, そのままに. total を DensDev ではなく
    !   Mass とか何とかな名前に変更する方が良いのではないか
    !

    ! モジュール読み込み
    !
    use dc_types,   only : DP, STRING
    use gtool_historyauto, &
         &          only : HistoryAutoAddVariable
    use constants,  only : Grav                   !重力加速度
    use basicset,   only : xyz_PressBZ            !基本場の気圧
    use gridset,    only : nz                     !格子点数
    use axesset,    only : xmin, xmax, ymin, ymax !計算領域
    
    !暗黙の型宣言禁止
    !
    implicit none

    ! 共通変数の設定
    !   grav = 0 の時にゼロ割になってしまうので, 細工している. 
    !
    MassTotalBZ  = &
      &  ( xyz_PressBz(1,1,1) - xyz_PressBz(1,1,nz) ) &
      &  / max( Grav, 1.0d-40 )                       &
      &  * ( xmax - xmin )                            &
      &  * ( ymax - ymin )

    ! 出力設定
    !
    call HistoryAutoAddVariable(             &
      & varname='DensDev',                   &
      & dims=(/'t'/),                        &
      & longname='Total density deviation',  &
      & units='kg',                          &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='KinEnrgy',                  &
      & dims=(/'t'/),                        &
      & longname='Total kinetic energy',     &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='ElstEnrgy',                 &
      & dims=(/'t'/),                        &
      & longname='Total elastic energy',     &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='PotEnrgy',                  &
      & dims=(/'t'/),                        &
      & longname='Total potential energy',   &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='DensDev2',                  &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Density deviation',        &
      & units='kg m-3',                      &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='KinEnrgyDens',              &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Kinetic energy density',   &
      & units='J m-3',                       &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='ElstEnrgyDens',             &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Elastic energy density',   &
      & units='J m-3',                       &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='PotEnrgyDens',              &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Potential energy density', &
      & units='J m-3',                       &
      & xtype='float')

  end subroutine EnergyMonit_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine EnergyMonit_exec( &
    & pyz_VelX, xqz_VelY, xyr_VelZ, xyz_Exner, xyz_PTemp &
    )
    !
    !保存量の計算と出力を行う. 計算する量は
    !質量, 運動エネルギー, 弾性エネルギー, ポテンシャルエネルギー
    !
    ! 追加すべき量(要検討)
    ! * 全運動量
    ! * エンタルピー
    ! * トレーサーの総量
    ! * 潜熱

    use dc_types,   only : DP
    use gtool_historyauto, &
         &          only : HistoryAutoPut
    use gridset,    only : imin,         & !x 方向の配列の下限
         &                 imax,         & !x 方向の配列の上限
         &                 jmin,         & !y 方向の配列の下限
         &                 jmax,         & !y 方向の配列の上限
         &                 kmin,         & !z 方向の配列の下限
         &                 kmax,         & !z 方向の配列の上限
         &                 nx, ny, nz
    use axesset,    only : dx, dy, dz,   & !格子点間隔
         &                 xyz_Z           !座標
    use average,    only : xyz_pyz,      &
      &                    xyz_xqz,      &
      &                    xyz_xyr
    use constants,  only : CpDry,        & !定圧比熱
         &                 CvDry,        & !定積比熱
         &                 GasRDry,      & !期待定数
         &                 PressSfc,     & !地表気圧
         &                 Grav            !重力加速度
    use basicset,   only : xyz_ExnerBZ,  & !基本場の圧力関数
         &                 xyz_PTempBZ,  & !基本場の温位
         &                 xyz_DensBZ,   & !基本場の密度
         &                 xyz_VelSoundBZ  !音速
    use timeset,    only : TimeN

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)   !速度(x成分)
    real(DP),intent(in) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)   !速度(y成分)
    real(DP),intent(in) :: xyr_VelZ(imin:imax,jmin:jmax,kmin:kmax)   !速度(z成分)
    real(DP),intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)  !圧力関数
    real(DP),intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)  !温位
    
    real(DP) :: xyz_DensDev(imin:imax,jmin:jmax,kmin:kmax) 
                                     !密度(偏差)
    real(DP) :: xyz_KinEnrgyDens(imin:imax,jmin:jmax,kmin:kmax) 
                                     !運動エネルギー密度
    real(DP) :: xyz_ElstEnrgyDens(imin:imax,jmin:jmax,kmin:kmax)
                                     !弾性エネルギー密度
    real(DP) :: xyz_PotEnrgyDens(imin:imax,jmin:jmax,kmin:kmax)
                                     !ポテンシャルエネルギー密度

    real(DP) :: MassTotalDev         !全質量(偏差)
    real(DP) :: KinEnrgyTotal        !全運動エネルギー
    real(DP) :: ElstEnrgyTotal       !全弾性エネルギー
    real(DP) :: PotEnrgyTotal        !全ポテンシャルエネルギー


    ! 各格子点の大気質量密度偏差
    ! 
    xyz_DensDev =                                             &
      &  PressSfc                                             &
      &  * (                                                  &  
      &      (xyz_ExnerBZ + xyz_Exner) ** ( CvDry / GasRDry ) &
      &       - (xyz_ExnerBZ ) ** ( CvDry / GasRDry )         &
      &    ) / xyz_PTempBZ 
    
    MassTotalDev = sum( dx * dy * dz * xyz_DensDev(1:nx,1:ny,1:nz) )

   
    ! 各格子点の運動エネルギー密度
    ! * 密度は基本場の値で評価
    ! 
    xyz_KinEnrgyDens =                     &
      &   0.5d0 * xyz_DensBZ               &
      &   * (                              &
      &        xyz_pyz(pyz_VelX) ** 2.0d0  &
      &      + xyz_xqz(xqz_VelY) ** 2.0d0  &
      &      + xyz_xyr(xyr_VelZ) ** 2.0d0  &
      &     )

    !各格子点の弾性エネルギー密度
    ! * 密度は基本場の値で評価
    !  
    xyz_ElstEnrgyDens = 0.5d0 * xyz_DensBZ     &
      &  * ( CpDry * xyz_PTempBZ * xyz_Exner / xyz_VelSoundBZ ) ** 2.0d0 
      
    ! 領域全体の運動/弾性エネルギー
    !
    KinEnrgyTotal = sum( dx * dy * dz * xyz_KinEnrgyDens(1:nx,1:ny,1:nz)  )
    ElstEnrgyTotal= sum( dx * dy * dz * xyz_ElstEnrgyDens(1:nx,1:ny,1:nz) )


    !各格子点のポテンシャルエネルギー密度
    ! 
    xyz_PotEnrgyDens = &
      &  - Grav * xyz_DensBZ * xyz_PTemp * xyz_Z / xyz_PTempBZ

    ! 領域全体のポテンシャルエネルギー
    !
    PotEnrgyTotal = sum( dx * dy * dz * xyz_PotEnrgyDens(1:nx,1:ny,1:nz))


    !ファイルへの出力
    !
    call HistoryAutoPut(TimeN, 'DensDev',   MassTotalDev   )
    call HistoryAutoPut(TimeN, 'KinEnrgy',  KinEnrgyTotal  )
    call HistoryAutoPut(TimeN, 'ElstEnrgy', ElstEnrgyTotal )
    call HistoryAutoPut(TimeN, 'PotEnrgy',  PotEnrgyTotal  )

    call HistoryAutoPut(TimeN, 'DensDev2',      xyz_DensDev(1:nx,1:ny,1:nz)      )
    call HistoryAutoPut(TimeN, 'KinEnrgyDens',  xyz_KinEnrgyDens(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'ElstEnrgyDens', xyz_ElstEnrgyDens(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'PotEnrgyDens',  xyz_PotEnrgyDens(1:nx,1:ny,1:nz) )

  end subroutine EnergyMonit_exec

end module EnergyMonit
