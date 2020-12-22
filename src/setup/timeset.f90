!= 時刻設定用モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: timeset.f90,v 1.11 2014/01/21 05:00:57 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module timeset
  !
  != 時刻設定用モジュール
  !
  !引数に与えられた NAMELIST ファイルから, 時刻に関する情報を取得し, 
  !保管するための変数型モジュール
  !
  !== Procedures List
  ! timeset_init       :: 初期化ルーチン
  ! TimesetDelTimeHalf :: リープフロッグスキームの最初の 1 ステップを回すためのルーチン
  ! TimesetProgress    :: リープフロッグスキームで 1 ステップを回すためのルーチン
  
  !モジュール読み込み
  use dc_types,      only: DP
  use dc_iounit,     only: FileOpen
  use dc_message,    only: MessageNotify
  use namelist_util, only: namelist_filename
  
  !暗黙の型宣言禁止
  implicit none
  
  ! Namelist から取得する変数
  !
  real(DP), save, public  :: DelTimeLong  = 2.0d0    !長いタイムステップ
  real(DP), save, public  :: DelTimeShort = 2.0d-1   !短いタイムステップ
  real(DP), save, public  :: IntegPeriod  = 3600.0d0 !積分時間
  real(DP), save, public  :: RestartTime  = 0.0d0    !計算開始時刻
  real(DP), save, private :: DelTimeOutput= 2.0d0    !出力タイムステップ

  ! public な変数
  !
  real(DP), save, public  :: TimeA                   !時刻 t + \del t
  real(DP), save, public  :: TimeN                   !時刻 t
  real(DP), save, public  :: TimeB                   !時刻 t - \del t
  real(DP), save, public  :: tfil           = 1.0d-1 !アセリンの時間フィルタの係数
  integer,  save, public  :: NstepLong      = 20     !長いタイムステップのステップ数
  integer,  save, public  :: NstepShort     = 20     !短いタイムステップのステップ数
  integer,  save, public  :: NstepOutput    = 20     !リスタートファイルへの出力
  integer,  save, public  :: Nstep                   !ループ回数
  logical,  save, public  :: FlagInitialRun = .false.!t=0 か否かのフラグ. 

  ! private な変数
  !
  real(DP), parameter, private :: tfilSave = 1.0d-1  !アセリンの時間フィルタの係数 (保管)
  real(DP), save, private :: DelTimeLongSave         !長いタイムステップ(保管)
  real(DP), save, private :: DelTimeShortSave        !短いタイムステップのステップ数(保管)
  integer,  save, private :: NstepInit      = 0      !リスタートファイルへの出力

  ! arare4 との互換性のために残しておいた変数
  !
  real(DP), save, public  :: TimeInt   = 3600.0d0    !積分時間
  integer,  save, public  :: NstepDisp = 20          !リスタートファイルへの出力  

  ! 公開要素
  public timeset_init, TimesetDelTimeHalf, TimesetProgress

contains
   
  subroutine timeset_init
    !
    != 初期化ルーチン
    !
    ! NAMELIST から必要な情報を読み取り, 時間関連の変数の設定を行う. 
    !

    !暗黙の型宣言禁止
    !
    implicit none

    !内部変数
    !
    real(DP), parameter :: order = 1.0
    integer             :: unit

    !---------------------------------------------------------------    
    ! NAMELIST から情報を取得
    !
    NAMELIST /timeset_nml/ &
      & DelTimeLong, DelTimeShort, IntegPeriod, RestartTime, DelTimeOutput
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=timeset_nml)
    close(unit)
    
    ! 保管
    !
    DelTimeLongSave  = DelTimeLong
    DelTimeShortSave = DelTimeShort
   
    !---------------------------------------------------------------
    ! 時間刻みなどを実数で与えているので, 整数に変換した上で, 
    ! ループを回す回数を決定することにした. 
    ! 
    
    ! 長いタイムステップ(前進差分)のループ数
    !    
    NstepLong = nint( IntegPeriod * order ) / nint( DelTimeLong * order ) 
    
    ! 短いタイムステップ(前進差分)のループ数
    !    
    NstepShort = nint( 2.0d0 * DelTimeLong * order ) / nint( DelTimeShort * order ) 
    
    ! リスタートファイルを書き出すタイミング
    !
    NstepOutput = nint( DelTimeOutput * order ) / nint( DelTimeLong * order )
    
    ! 計算開始時のループ数
    !
    NstepInit = nint( RestartTime * order ) / nint( DelTimeLong * order )

    ! 現在のループ回数
    !
    Nstep = 0
    
    ! 割り切れない場合はエラーを出す
    ! 積分時間 / dt /= 0.0
    !
    if ( mod( nint( IntegPeriod * order ), nint( DelTimeLong * order ) ) /= 0 ) then
      call MessageNotify( "E", "timeset_init", "mod( IntegPeriod, DelTimeLong ) /= 0")
    end if

    ! 割り切れない場合はエラーを出す
    ! dt / d\tau /= 0.0
    !
    if ( mod( nint( DelTimeLong * order ), nint( DelTimeShort * order ) ) /= 0 ) then
      call MessageNotify( "E", "timeset_init", "mod( DelTimeLong, DelTimeShort ) /= 0")
    end if

    ! 割り切れない場合はエラーを出す
    ! リスタートファイル書き出しのタイミング (DelTimeOutput) / dt /= 0.0
    !
    if ( mod( nint( DelTimeOutput * order ), nint( DelTimeLong * order ) ) /= 0 ) then
      call MessageNotify( "E", "timeset_init", "mod( DelTimeOutput, DelTimeLong ) /= 0")
    end if

    ! 割り切れない場合はエラーを出す
    ! リスタート時刻 / dt /= 0.0
    !
    if ( mod( nint( RestartTime * order ), nint( DelTimeLong * order ) ) /= 0 ) then
      call MessageNotify( "E", "timeset_init", "mod( RestartTime, DelTimeLong ) /= 0")
    end if

    ! リスタートファイルの出力のタイミングが 1 以上でないとエラーを出す. 
    !
    if ( NstepOutput < 1 ) then 
      call MessageNotify( "E", "timeset_init", "NstepOutput is < 1, %d", i=(/NstepOutput/) )
    end if

    !---------------------------------------------------------------
    ! リスタートか否かを判断した上で時刻を設定する. 
    !
    if ( NstepInit == 0 ) then 
      FlagInitialRun = .true.
      Nstep = 0
      TimeB = 0.0d0
      TimeN = 0.0d0
      TimeA = DelTimeLong
    else
      FlagInitialRun = .false.
      Nstep = 1
      TimeB = ( NstepInit + Nstep - 1 ) * DelTimeLong
      TimeN = ( NstepInit + Nstep     ) * DelTimeLong
      TimeA = ( NstepInit + Nstep + 1 ) * DelTimeLong
    end if


    !---------------------------------------------------------------
    ! 確認
    !
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeShort = %f", d=(/DelTimeShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "IntegPeriod  = %f", d=(/IntegPeriod/) )
    call MessageNotify( "M", &
      & "timeset_init", "Restarttime  = %f", d=(/Restarttime/)  )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeOutput= %f", d=(/DelTimeOutput/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepLong    = %d", i=(/NstepLong/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepShort   = %d", i=(/NstepShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepOutput  = %d", i=(/NstepOutput/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepInit    = %d", i=(/NstepInit/) )
    call MessageNotify( "M", &
      & "timeset_init", "tfil  = %f", d=(/tfil/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeA = %f", d=(/TimeA/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeN = %f", d=(/TimeN/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeB = %f", d=(/TimeB/) )
    call MessageNotify( "M", &
      & "timeset_init", "FlagInitialRun  = %b", L=(/FlagInitialRun/) )
    
    !---------------------------------------------------------------
    ! arare4 用の変数の設定
    !
    NStepDisp = NStepOutput
    TimeInt   = IntegPeriod

  end subroutine timeset_init


  subroutine TimesetDelTimeHalf    
    !
    != リープフロッグスキームの最初の 1 ステップを回すためのルーチン
    !
    ! 最初の一回目はオイラースキームで回すための処置
    ! * 短い時間ステップのループ回数も半分にする. 
    ! * 長い時間ステップの DelTimeLong を半分にする.
    ! * Asselin の時間フィルタの係数をゼロとする. 
    
    ! 暗黙の型宣言禁止
    !
    implicit none

    ! FlagInitialRun = .false. の時はこのルーチンを呼ぶ必要は無い. 
    !
    if ( .NOT. FlagInitialRun ) then 
       call MessageNotify( "E", "timeset_init_TimesetDelTimeHalf", "FlagInitialRun = false" )
    end if

    ! 時刻刻み幅を直す 
    !
    DelTimeShort = DelTimeShortSave * 5.0d-1   
    DelTimeLong  = DelTimeLongSave  * 5.0d-1   

    ! Asselin の時間フィルタの係数
    !
    tfil = 0.0d0

    !---------------------------------------------------------------
    ! 確認
    !
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial DelTimeShort = %f", d=(/DelTimeShort/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Asselin Time Filter coefficient = %f", d=(/tfil/) )

  end subroutine TimesetDelTimeHalf
  

  subroutine TimesetProgress
    !
    != リープフロッグスキームで 1 ステップを回すためのルーチン
    !
    ! 時刻を進め, 時刻の刻み幅と Asselin の時間フィルタの係数を直す. 
    ! 

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! call された回数を 1 増加させる
    !
    Nstep = Nstep + 1

    ! 時刻刻み幅を直す 
    !
    DelTimeShort = DelTimeShortSave
    DelTimeLong  = DelTimeLongSave

    ! ループ回数から時刻を決める. 
    !
    TimeB = ( NstepInit + Nstep - 1 ) * DelTimeLong
    TimeN = ( NstepInit + Nstep     ) * DelTimeLong
    TimeA = ( NstepInit + Nstep + 1 ) * DelTimeLong
    
    ! Asselin の時間フィルタの係数
    !
    tfil = tfilSave
  
  end subroutine TimesetProgress
  
end module timeset
