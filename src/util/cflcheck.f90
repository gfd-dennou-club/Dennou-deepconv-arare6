!= CFL 条件のチェックをするためのパッケージ型モジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: cflcheck.f90,v 1.3 2014/03/04 05:55:06 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module CFLCheck
  !
  ! CFL 条件のチェックをするためのパッケージ型モジュール
  !  * 音波に対して CFL 条件をチェック
  !  * 入力された速度に対して CFL 条件をチェック
  !

  !暗黙の型宣言禁止
  implicit none

  !private 属性の指定
  private

  !型宣言
  character(*), parameter :: module_name = 'cflcheck'
                              ! モジュールの名称.
                              ! Module name
  
  !関数を public 属性に設定
  public CFLCheckTimeShort
  public CFLCheckTimeLongVelX
  public CFLCheckTimeLongVelY
  public CFLCheckTimeLongVelZ
  
contains  

!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeShort( xyz_VelSound )
    !
    !音波に対して CFL 条件をチェック
    ! 

    !モジュール読み込み
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: FlagCalc3D, &
      &                   imin,        &! x 方向の配列の下限
      &                   imax,        &! x 方向の配列の上限
      &                   jmin,        &! z 方向の配列の下限
      &                   jmax,        &! z 方向の配列の上限
      &                   kmin,        &! z 方向の配列の下限
      &                   kmax          ! z 方向の配列の上限
    use axesset,    only: dx,          &! x 方向の格子点間隔
      &                   dy            ! y 方向の格子点間隔
    use timeset,    only: DelTimeShort  !短い時間ステップ
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xyz_VelSound(imin:imax,jmin:jmax, kmin:kmax)
                                        !音速
    real(DP)             :: xyz_CFL(imin:imax,jmin:jmax, kmin:kmax)
                                        !クーラン数
    
    !音速と CFL 条件を求める
    !
    if ( FlagCalc3D ) then 
      xyz_CFL = DelTimeShort * xyz_VelSound       &
        &       * ((1.0d0 / (dx * dx) + 1.0d0 / (dy * dy)) ** 0.5d0)
    else
      xyz_CFL = DelTimeShort * xyz_VelSound / dx
    end if

    !メッセージ
    !
    call MessageNotify( "M", &
      & module_name, &
      & "Sound Wave Velocity = %f", d=(/maxval(xyz_VelSound)/) )
    call MessageNotify( "M", &
      & module_name, &
      & "DelTimeShort = %f", d=(/DelTimeShort/) )
    
    !警告メッセージ
    if ( maxval(xyz_CFL) >= 1.0) then 
      call MessageNotify( "E", &
        & module_name, &
        & "CFL Condition is broken, DelTimeShort * VelSound > min(DelX, DelZ)")
    else
      call MessageNotify( "M", &
        & module_name, &
        & "Courant number for DelTimeSort = %f", d=(/maxval(xyz_CFL)/) )
    end if

  end subroutine CFLCheckTimeShort


!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelX( pyz_VelX )
    !
    !水平速度に対して CFL 条件をチェック. 
    ! 

    !モジュール読み込み
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x 方向の配列の下限
      &                   imax,      &! x 方向の配列の上限
      &                   jmin,      &! z 方向の配列の下限
      &                   jmax,      &! z 方向の配列の上限
      &                   kmin,      &! z 方向の配列の下限
      &                   kmax        ! z 方向の配列の上限
    use axesset,    only: dx          ! x 方向の格子点間隔
    use timeset,    only: DelTimeLong !長い時間ステップ

    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    real(DP), intent(in) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL

    !CFL 条件を求める
    CFL = ( 2.0d0 * DelTimeLong ) * maxval( abs( pyz_VelX / dx ) )
  
    !メッセージ出力
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelX for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelX)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelX
    

!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelY( xqz_VelY )
    !
    !水平速度に対して CFL 条件をチェック. 
    ! 

    !モジュール読み込み
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x 方向の配列の下限
      &                   imax,      &! x 方向の配列の上限
      &                   jmin,      &! z 方向の配列の下限
      &                   jmax,      &! z 方向の配列の上限
      &                   kmin,      &! z 方向の配列の下限
      &                   kmax        ! z 方向の配列の上限
    use axesset,    only: dy          ! y 方向の格子点間隔
    use timeset,    only: DelTimeLong !長い時間ステップ

    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    real(DP), intent(in) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL

    !CFL 条件を求める
    CFL = (2.0d0 * DelTimeLong) * maxval( abs( xqz_VelY / dy ) )
  
    !メッセージ出力
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelY for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelY)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelY


!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelZ( xyr_VelZ )
    !
    !水平速度に対して CFL 条件をチェック. 
    ! 

    !モジュール読み込み
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x 方向の配列の下限
      &                   imax,      &! x 方向の配列の上限
      &                   jmin,      &! z 方向の配列の下限
      &                   jmax,      &! z 方向の配列の上限
      &                   kmin,      &! z 方向の配列の下限
      &                   kmax        ! z 方向の配列の上限
    use axesset,    only: dz          ! z 方向の格子点間隔
    use timeset,    only: DelTimeLong !長い時間ステップ

    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    real(DP), intent(in) :: xyr_VelZ(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL
    
    !CFL 条件を求める
    CFL = ( 2.0d0 * DelTimeLong ) * maxval( abs( xyr_VelZ / dz ) ) 
    
    !メッセージ出力
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelZ for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelZ)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelZ
  
end module CFLCheck
