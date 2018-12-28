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

  !属性
  private
  
  ! Public Interface
  real(DP), save, public  :: TimeA                   !時刻 t + \del t
  real(DP), save, public  :: TimeN                   !時刻 t
  real(DP), save, public  :: TimeB                   !時刻 t - \del t
  real(DP), save, public  :: DelTimeLong  = 2.0d0    !長いタイムステップ
  real(DP), save          :: DelTimeLongSave  = 2.0d0 !長いタイムステップ
  real(DP), save, public  :: DelTimeShort = 2.0d-1   !短いタイムステップ
  real(DP), save, public  :: RestartTime  = 0.0d0    !計算開始時刻
  real(DP), save, public  :: IntegPeriod  = 3600.0d0 !積分時間
  real(DP), save, public  :: EndTime      = 3600.0d0 !計算終了時刻
  real(DP), save, public  :: DelTimeOutput= 2.0d0    !出力タイムステップ
  real(DP), save, public  :: tfil         = 1.0d-1   !アセリンの時間フィルタの係数
  real(DP), parameter     :: tfilSave     = 1.0d-1   !アセリンの時間フィルタの係数
  integer,  save, public  :: NstepShort = 20         !短いタイムステップのステップ数
  integer,  save          :: NstepShortSave = 20     !短いタイムステップのステップ数
  integer,  save, public  :: NstepOutput    = 20     !リスタートファイルへの出力
  logical,  save, public  :: FlagInitialRun = .false.!t=0 か否かのフラグ. 

  integer,  save, public  :: NstepDisp = 20          !リスタートファイルへの出力

  real(DP), save, public  :: TimeInt  = 3600.0d0     !積分時間
  
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
    implicit none

    !内部変数
    integer    :: unit

    !---------------------------------------------------------------    
    ! NAMELIST から情報を取得
    !
    NAMELIST /timeset_nml/ &
      & DelTimeLong, DelTimeShort, IntegPeriod, RestartTime, DelTimeOutput
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=timeset_nml)
    close(unit)

    ! 計算終了時刻
    !
    EndTime = RestartTime + IntegPeriod
    
    ! 時刻・タイムステップの設定
    !   実数の割り算なので, 念の為に. 
    !    
    NstepShort = ( nint( 2.0d0 * DelTimeLong * 1.0d3 ) / nint( DelTimeShort * 1.0d3 ) )

    ! 時間刻みの設定
    !   元々の値を保管しておく
    DelTimeLongSave = DelTimeLong
    NstepShortSave  = NstepShort

    ! リスタートファイルを書き出すタイミング
    !
    NstepOutput = nint( DelTimeOutput * 1.0d2 ) / nint( DelTimeLong * 1.0d2 )
    if ( NstepOutput < 0 ) then 
       NstepOutput = nint( DelTimeOutput * 1.0d1 ) / nint( DelTimeLong * 1.0d1 )
    end if
    if ( NstepOutput < 0 ) then 
       NstepOutput = nint( DelTimeOutput ) / nint( DelTimeLong )
    end if
    if ( NstepOutput < 0 ) then 
      call MessageNotify( "E", &
        & "timeset_init", "NstepOutput is negative, %d", i=(/NstepOutput/) )
    end if


    ! 時刻の設定
    !
    TimeB = RestartTime - DelTimeLong
    TimeN = RestartTime
    TimeA = RestartTime + DelTimeLong      
    
    ! リスタートか否か. リスタートなら  .false.
    !
    if ( nint( RestartTime * 1.0d3 ) == 0 ) then 
      FlagInitialRun = .true.
    else
      FlagInitialRun = .false.
    end if

    !---------------------------------------------------------------
    ! 確認
    !
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeShort = %f", d=(/DelTimeShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "Restarttime  = %f", d=(/Restarttime/)  )
    call MessageNotify( "M", &
      & "timeset_init", "IntegPeriod  = %f", d=(/IntegPeriod/) )
    call MessageNotify( "M", &
      & "timeset_init", "EndTime      = %f", d=(/EndTime/) )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeOutput= %f", d=(/DelTimeOutput/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepShort   = %d", i=(/NstepShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepOutput  = %d", i=(/NstepOutput/) )
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
    ! 最初の一回目はオイラースキームで回すので, 
    ! 長い時間刻みを半分にし, 短い時間で回すループ回数も半分にする.
    ! Asselin の時間フィルタの係数をゼロとする. 

    implicit none

    TimeB       = TimeN
    DelTimeLong = DelTimeLongSave * 0.5d0
    NstepShort  = NstepShortSave  / 2
    tfil        = 0.0d0

    !---------------------------------------------------------------
    ! 確認
    !
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial NstepShort   = %d", i=(/NstepShort/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Asselin Time Filter coefficient = %f", d=(/tfil/) )

  end subroutine TimesetDelTimeHalf
  

  subroutine TimesetProgress
    !
    ! = リープフロッグスキームで 1 ステップを回すためのルーチン
    !
    ! 時刻を進め, 時刻の刻み幅と Asselin の時間フィルタの係数を直す. 
    ! 

    implicit none

    ! 時刻刻み幅を直す 
    !
    DelTimeLong = DelTimeLongSave
    NstepShort  = NstepShortSave

    ! 時刻を進める
    !
    TimeB = TimeN
    TimeN = TimeA
    TimeA = TimeA + DelTimeLong

    ! Asselin の時間フィルタの係数
    tfil = tfilSave

  end subroutine TimesetProgress
  
end module timeset
