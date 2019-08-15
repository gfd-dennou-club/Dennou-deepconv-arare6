!= Module HistoryFileIO
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: historyfileio.f90,v 1.9 2014/03/04 05:55:04 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module HistoryFileIO
  !
  !ファイル出力. 長い時間ステップの値を出力.
  !

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types, only: STRING     ! 文字列. Strings. 

  !暗黙の型宣言禁止
  implicit none
  
  !属性の指定
  private

  !公開手続き
  public HistoryFileIO_init
  public HistoryFileIO_finalize

  ! 非公開変数
  ! Private variables
  !
  character(STRING), parameter, private :: module_name = 'historyfileio'
                              ! モジュールの名称. 
                              ! Module name
  character(STRING), parameter, private :: version = &
    & '$Name:  $' // &
    & '$Id: historyfileio.f90,v 1.9 2014/03/04 05:55:04 sugiyama Exp $'
                              ! モジュールのバージョン
                              ! Module version

contains 

!!!------------------------------------------------------------------------
  subroutine HistoryFileIO_init
    !
    ! history_file_io モジュールの初期化を行います. 
    !--
    ! NAMELIST#history_file_io_nml の読み込みはこの手続きで行われます. 
    !++
    !
    ! "history_file_io" module is initialized. 
    !--
    ! "NAMELIST#history_file_io_nml" is loaded in this procedure. 
    !++
    !

    ! モジュール引用 ; USE statements
    !

    ! gtool5 netCDF データの入出力インターフェース (大規模モデル用)
    ! Interface of Input/Output of gtool5 netCDF data (For large models)
    !
    use gtool_historyauto, only: HistoryAutoCreate,  &
      &                          HistoryAutoAddAttr, &
      &                          HistoryAutoPutAxis, &
      &                          HistoryAutoAddVariable

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen
    
    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types,      only: DP, &              ! 倍精度実数型. Double precision. 
      &                      STRING             ! 文字列.       Strings. 
    use mpi_wrapper,   only: FLAG_LIB_MPI
    use namelist_util, only: namelist_filename    
    use axesset,       only: x_X, y_Y, z_Z
    use gridset,       only: nx, ny, nz, ncmax
    use fileset,       only: filetitle,        &!データの表題
      &                      filesource,       &!データを作成する手順
      &                      FileInstitution    !最終変更者・組織
    use timeset,       only: Restarttime, IntegPeriod
    use composition,   only: SpcWetSymbol, GasNum
    
    ! 宣言文 ; Declaration statements
    !
    implicit none
    
    !変数定義
    real(DP), parameter :: TimeDisp = 1.0e5 !出力間隔のデフォルト値
    real(DP) :: EndTime
    integer  :: l, s

    EndTime = RestartTime + IntegPeriod

    !-----------------------------------------------------------
    ! ヒストリー作成
    !-----------------------------------------------------------
    call HistoryAutoCreate(                             &
      & title = FileTitle,                              &
      & source = FileSource,                            &
      & institution = FileInstitution,                  &
      & dims=(/'x','y','z','t'/),                       &
      & dimsizes=(/nx, ny, nz, 0/),                     &
      & longnames=(/'x-coordinate',                     &
      &             'y-coordinate',                     &
      &             'z-coordinate',                     &
      &             'time        '/),                   &
      & units=(/'m  ','m  ','m  ','sec'/),              &
      & xtypes=(/'double','double','double','double'/), &
      & origin   = Restarttime,                         & 
      & terminus = EndTime,                             &
      & interval = TimeDisp,                            &       
      & flag_mpi_split = FLAG_LIB_MPI,                  &
      & namelist_filename = namelist_filename)  

    call HistoryAutoAddAttr( &
      & varname = 'x', attrname = 'standard_name', &   ! (in)
      & value = 'x-coordinate' )                       ! (in)
    call HistoryAutoAddAttr( &
      & varname = 'y', attrname = 'standard_name', &   ! (in)
      & value = 'y-coordinate' )                       ! (in)
    call HistoryAutoAddAttr( &
      & varname = 'z', attrname = 'standard_name', &   ! (in)
      & value = 'z-coordinate' )                       ! (in)
    
    call HistoryAutoPutAxis('x', x_X(1:nx))
    call HistoryAutoPutAxis('y', y_Y(1:ny))
    call HistoryAutoPutAxis('z', z_Z(1:nz))

    ! 印字 ; Print
    !
!    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
!    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    !無次元圧力の擾乱
    call HistoryAutoAddVariable(                           &
      & varname='Exner',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='disturbunce of nondimensional pressure', &
      & units=' ',                                         &
      & xtype='float' )

    !無次元圧力
    call HistoryAutoAddVariable(                           &
      & varname='ExnerAll',                                &
      & dims=(/'x','y','z','t'/),                          &
      & longname='nondimensional pressure',                &
      & units=' ',                                         &
      & xtype='float' )
    
    call HistoryAutoAddVariable(                           &
      & varname='PTemp',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='disturbunce of potential temperature',   &
      & units='K',                                         &
      & xtype='float' )

    !温位の擾乱
    call HistoryAutoAddVariable(                           &
      & varname='PTempAll',                                &
      & dims=(/'x','y','z','t'/),                          &
      & longname='potential temperature',                  &
      & units='K',                                         &
      & xtype='float' )

    !水平速度
    call HistoryAutoAddVariable(                           &
      & varname='VelX',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='zonal velocity',                         &
      & units='m.s-1',                                     &
      & xtype='float' )

    !水平速度
    call HistoryAutoAddVariable(                           &
      & varname='VelY',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='meridional velocity',                    &
      & units='m.s-1',                                     &
      & xtype='float' )

    !鉛直速度
    call HistoryAutoAddVariable(                           &
      & varname='VelZ',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='vertical velocity',                      &
      & units='m.s-1',                                     &
      & xtype='float' )

    !渦粘性係数(運動量)
    call HistoryAutoAddVariable(                           &
      & varname='Km',                                      &
      & dims=(/'x','y','z','t'/),                          &
      & longname='turbulet diffusion coefficient',         &
      & units='m2.s-1',                                    &
      & xtype='float' )
  
    !渦粘性係数(熱)
    call HistoryAutoAddVariable(                           &
      & varname='Kh',                                      &
      & dims=(/'x','y','z','t'/),                          &
      & longname='turbulet diffusion coefficient for heat',&
      & units='m2.s-1',                                    &
      & xtype='float')

    !雲密度
    call HistoryAutoAddVariable(                           &
      & varname='CDens',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='Cloud density',                          &
      & units='kg.m-3',                                    &
      & xtype='float')
    
    ! 混合比
    do l = 1, ncmax
      call HistoryAutoAddVariable(                         &
        & varname=trim(SpcWetSymbol(l)),                   &
        & dims=(/'x','y','z','t'/),                        &
        & longname=trim(SpcWetSymbol(l))//' Mixing Ratio', &
        & units='kg.kg-1',                                 &
        & xtype='float')
    end do

    do l = 1, GasNum      
      call HistoryAutoAddVariable(                         &
        & varname=trim(SpcWetSymbol(l))//'All',            &
        & dims=(/'x','y','z','t'/),                        &
        & longname=trim(SpcWetSymbol(l))//' Mixing Ratio', &
        & units='kg.kg-1',                                 &
        & xtype='float')
    end do


    !----------------------------------------------------------------
    ! Mixing Ratio time change
    !----------------------------------------------------------------
    do s = 1, ncmax
    
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(s))//'DtFill1',       &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Filling Negative term 1 of '               &
        &           //trim(SpcWetSymbol(s))//' mixing ratio',  &
        & units='kg.kg-1.s-1',                                 &
        & xtype='float' )
    
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(s))//'DtFill2',             &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Filling Negative term 2 of '               &
        &           //trim(SpcWetSymbol(s))//' mixing ratio',  &
        & units='kg.kg-1.s-1',                                 &
        & xtype='float' )
    
    end do
   
  end subroutine HistoryFileIO_init


  subroutine HistoryFileIO_finalize
    !
    ! ヒストリデータファイル出力の終了処理を行います. 
    !
    ! Terminate history data files output. 

    ! モジュール引用 ; USE statements
    !
    use gtool_historyauto, only: HistoryAutoClose

    ! 宣言文 ; Declaration statements
    !
    implicit none

    ! 作業変数
    ! Work variables
    !

    ! 実行文 ; Executable statement
    !

    call HistoryAutoClose

  end subroutine HistoryFileIO_finalize
  
end module HistoryFileIO


