! CPU 時間の計測を行うためのモジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: clockset.f90,v 1.3 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module clockset
  !
  != CPU 時間の計測を行うためのモジュール
  !
  ! CPU 時間の計測を行うためのサブルーチンを束ねたモジュール
  !
  !== Procedures List
  ! ClocksetInit      :: 初期化ルーチン
  ! ClocksetClose     :: CPU 時間計測の終了処理
  ! ClocksetPredict   :: プログラムが終了までの予測 CPU 時間・日時を表示
  ! ClockSetPreStart  :: CPU 時間計測開始 (初期化ルーチン)
  ! ClockSetPreStop   :: CPU 時間計測終了 (初期化ルーチン)
  ! ClockSetLoopStart :: CPU 時間計測開始 (時間ループ)
  ! ClocksetLoopStop  :: CPU 時間計測終了 (時間ループ)


  !モジュール読み込み
  !
  use dc_clock,   only : CLOCK

  !暗黙の型宣言禁止
  !
  implicit none

  !属性の指定
  !
  private
  
  ! Public Interface
  !
  type(CLOCK), save :: clock_init, clock_loop  ! Variables for CPU time counting 
                                               ! CPU 時間計測用変数
  public ClocksetInit, ClocksetClose, ClocksetPredict
  public ClockSetPreStart, ClockSetPreStop, ClockSetLoopStart, ClocksetLoopStop

contains
   
  subroutine ClocksetInit
    ! 
    ! CPU 時間計測 初期化ルーチン
    !

    !モジュール読み込み
    use dc_clock,    only : DCClockCreate

    implicit none
    
    ! 初期化ルーチン用の時刻計測初期化
    !
    call DCClockCreate( &           ! Initialize (初期化)
      & clk = clock_init, &         ! (out)
      & name = 'initialization' )   ! (in)

    ! 時間発展ループ用の時刻計測初期化
    !
    call DCClockCreate( &           ! Initialize (初期化)
      & clk = clock_loop, &         ! (out)
      & name = 'time-integration' ) ! (in)

  end subroutine ClocksetInit


  subroutine ClocksetPreStart
    !
    ! CPU 時間計測開始 (初期化ルーチン)
    !

    !モジュール読み込み
    use dc_clock,   only : DCClockStart
    
    !暗黙の型宣言禁止
    implicit none

    ! Start CPU time counting     
    call DCClockStart(clk = clock_init) ! (inout)
    
  end subroutine ClocksetPreStart
  

  subroutine ClocksetPreStop
    !
    ! CPU 時間計測修了 (初期化ルーチン)
    !
    
    !モジュール読み込み
    use dc_clock,   only : DCClockStop

    !暗黙の型宣言禁止
    implicit none
    
    ! Stop CPU time counting 
    ! (CPU 時間計測終了)    
    call DCClockStop(clk = clock_init) ! (inout)

  end subroutine ClocksetPreStop


  subroutine ClocksetLoopStart
    !
    ! CPU 時間計測開始 (時間ループ)
    !
    
    !モジュール読み込み
    use dc_clock,   only : DCClockStart

    !暗黙の型宣言禁止
    implicit none

    ! Start CPU time counting 
    call DCClockStart(clk = clock_loop) ! (inout) 

  end subroutine ClocksetLoopStart


  subroutine ClocksetLoopStop
    !
    ! CPU 時間計測終了 (時間ループ)
    !

    !モジュール読み込み
    use dc_clock,   only : DCClockStop

    !暗黙の型宣言禁止
    implicit none

    !Stop CPU time counting 
    !(CPU 時間計測終了)    
    call DCClockStop(clk = clock_loop) ! (inout)

  end subroutine ClocksetLoopStop


  subroutine ClocksetPredict
    !
    ! プログラムが終了するまでの予測 CPU 時間, および日時を表示
    ! 一旦計測を止めてから, 予想値を表示し, 計測を再開する. 
    !

    !モジュール読み込み
    use dc_message,  only : MessageNotify
    use dc_clock,    only : DCClockStop, DCClockPredict, operator(+)
    use timeset,     only : TimeA, RestartTime, IntegPeriod
    
    !暗黙の型宣言禁止
    implicit none

    !作業変数
    real(4) :: progress 

    ! 計算の進行の割合を計算
    progress = real((TimeA - RestartTime) / IntegPeriod, 4)

    ! CPU time measurement stops, temporarily.
    call ClocksetLoopStop

    call MessageNotify( "M", "ClockSet", "Time = %f", d=(/TimeA/) )

    call DCClockPredict( &                       ! 残り時間の予測
      &   clk = clock_init + clock_loop,       & ! (in)
      &   progress = progress                  & ! (in) 
      &  )

    ! CPU time measurement starts again.
    call ClocksetLoopStart

  end subroutine ClocksetPredict


  subroutine ClocksetClose
    !
    ! CPU 時間計測の終了処理
    ! 計測を終了し, CPU 時間を表示する.
    !

    !モジュール読み込み
    use dc_clock,    only : DCClockClose, DCClockResult

    !暗黙の型宣言禁止
    implicit none

    call DCClockResult( &                    ! 全 CPU 時間の表示
      & clks = (/clock_init, clock_loop/), & ! (in)
      & total_auto = .true. )                ! (in)

    call DCClockClose( clk = clock_init )      ! (inout) ! Finalize (後処理)
    call DCClockClose( clk = clock_loop )      ! (inout) ! Finalize (後処理)
    
  end subroutine ClocksetClose
  
end module clockset
