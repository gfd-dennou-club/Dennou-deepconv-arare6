!= Module ReStartFileIO
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: restartfileio.f90,v 1.18 2015/02/19 02:17:23 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module ReStartFileIO
  !
  !リスタート用の場の情報を netCDF ファイルに出力するためのルーチン
  !

  !モジュール読み込み
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING

  !暗黙の型宣言禁止
  implicit none
  
  !関数を public に指定
  public ReStartFileIO_init
  public ReStartFileio_Finalize
  public ReStartFileio_BasicZ_Get
  public ReStartFileio_Var_Get

  type(GT_HISTORY), save, public   :: rstat
  character(STRING), save, private :: InitialFile = ""
  character(STRING), save, private :: InputFile   = ""
  character(STRING), save, private :: OutputFile  = "output.nc"
  
contains

  subroutine ReStartFileio_Init ( FlagInitData )
    !
    !リスタートファイルの書き出し
    !

    !モジュール読み込み
    use gtool_history, only : HistoryCreate,  &
      &                       HistoryPut,     &
      &                       HistoryAddVariable
    use dc_message,    only : MessageNotify
    use dc_iounit,     only : FileOpen
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use axesset,       only : x_X,            &!X 座標軸(スカラー格子点)
      &                       y_Y,            &!Y 座標軸(スカラー格子点)
      &                       z_Z              !Z 座標軸(スカラー格子点)
    use gridset,       only : ncmax            !凝縮成分の数
    use fileset,       only : filetitle,      &!データの表題
      &                       filesource,     &!データを作成する手順
      &                       FileInstitution  !最終変更者・組織
    use namelist_util, only: namelist_filename
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    logical, intent(in), optional :: FlagInitData  ! 初期値生成プログラムの場合は .true.

    real(4)            :: SpcID(ncmax)
    integer            :: N, L, M
    integer            :: s    
    integer            :: unit     !装置番号
     
    !NAMELIST から情報を取得
    NAMELIST /restartfileio_nml/ InitialFile, InputFile, OutputFile
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=restartfileio_nml)
    close(unit)

    !確認
    call MessageNotify( "M", &
      & "restartfileioIO_init", "InitialFile  = %c", c1=trim(InitialFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "InputFile  = %c", c1=trim(InputFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "OutputFile = %c", c1=trim(OutputFile) )

    ! InputFile が空の場合は, 初期値ファイルを使う. 
    if ( InputFile == "" ) then 
      InputFile = InitialFile
    end if

    !初期値生成の場合には, outputfile 名を変更
    if ( present( FlagInitData ) ) then 
      if ( FlagInitData ) then 
        if ( .NOT. InitialFile == "" ) then  ! ファイル名が空ならば
          InputFile  = ""
          OutputFile = InitialFile
        end if
      end if
    end if

    !確認
    call MessageNotify( "M", &
      & "restartfileioIO_init", "INPUT  = %c", c1=trim(InputFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "OUTPUT = %c", c1=trim(OutputFile) )

    SpcID = 0.0d0
    do s = 1, ncmax
      SpcID(s) = real( s, 4 )
    end do
    
    N = size(x_X, 1)
    L = size(y_Y, 1)
    M = size(z_Z, 1)
    
    !-------------------------------------------------------------    
    ! ヒストリー作成
    !-------------------------------------------------------------  
    call HistoryCreate(                              &
      & file = Outputfile,                           &
      & title = filetitle,                           &
      & source = filesource,                         &
      & institution = FileInstitution,               &
      & dims=(/'x','y','z','s','t'/),                &
      & dimsizes=(/N, L, M, ncmax, 0/),              &
      & longnames=(/'X-coordinate',                  &
      &             'Y-coordinate',                  &
      &             'Z-coordinate',                  &
      &             'Species Num ',                  &
      &             'Time        '/),                &
      & units=(/'m  ','m  ','m  ','1  ','sec'/),     &
      & xtypes=(/'double', 'double', 'double', 'double', 'double'/), &
      & flag_mpi_split = FLAG_LIB_MPI,               &
      & origin=0.0, interval=1.0,                    &
      & history=rstat, quiet=.true. )
    
    !-------------------------------------------------------------  
    ! 変数出力
    !-------------------------------------------------------------
    call HistoryPut('x', x_X, rstat )
    call HistoryPut('y', y_Y, rstat )
    call HistoryPut('z', z_Z, rstat )
    call HistoryPut('s', real(SpcID, 4), rstat )

    !無次元圧力の基本場
    call HistoryAddVariable(                           &
      & varname='ExnerBZ', dims=(/'x','y','z'/),       &
      & longname='nondimensional pressure', units='1', &
      & xtype='double', history=rstat )
    
    !温位の基本場
    call HistoryAddVariable(                           &
      & varname='PTempBZ', dims=(/'x','y','z'/),       &
      & longname='potential temperature',              &
      & units='K', xtype='double', history=rstat ) 

    !温位の基本場
    call HistoryAddVariable(                           &
      & varname='VPTempBZ', dims=(/'x','y','z'/),      &
      & longname='virtual potential temperature',      &
      & units='K', xtype='double', history=rstat ) 
    
    !密度の基本場
    call HistoryAddVariable(                           &
      & varname='DensBZ', dims=(/'x','y','z'/),        &
      & longname='density',                            &
      & units='Kg.m-3', xtype='double', history=rstat )
    
    !音波速度の基本場
    call HistoryAddVariable(                           &
      & varname='VelSoundBZ', dims=(/'x','y','z'/),    &
      & longname='sound velocity',                     &
      & units='m.s-2', xtype='double', history=rstat )

    !温度の基本場
    call HistoryAddVariable(                           &
      & varname='TempBZ', dims=(/'x','y','z'/),        &
      & longname='Temperature of basic state',         &
      & units='K', xtype='double', history=rstat ) 
    
    !圧力の基本場
    call HistoryAddVariable(                           &
      & varname='PressBZ', dims=(/'x','y','z'/),       &
      & longname='Pressure of basic state',            &
      & units='Pa', xtype='double', history=rstat ) 
    
    !水蒸気混合比の基本場
    call HistoryAddVariable(                              &
      & varname='QMixBZ', dims=(/'x','y','z','s'/),       &
      & longname='Mixing ratio of Condensible volatiles', &
      & units='kg.kg-1', xtype='double', history=rstat ) 
    
    !分子量効果
    call HistoryAddVariable(                         &
      & varname='EffMolWtBZ', dims=(/'x','y','z'/),  &
      & longname='Effect of Mole Weight',            &
      & units='1', xtype='double', history=rstat ) 

    !湿度
    call HistoryAddVariable(                        &
      & varname='HumBZ', dims=(/'x','y','z','s'/),  &
      & longname='Humidity',                        &
      & units='1', xtype='double', history=rstat ) 
    
    !速度
    call HistoryAddVariable(                         &
      & varname='VelX', dims=(/'x','y','z','t'/),    &
      & longname='zonal velocity',                   &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !速度
    call HistoryAddVariable(                         &
      & varname='VelY', dims=(/'x','y','z','t'/),    &
      & longname='meridional velocity',              &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !速度
    call HistoryAddVariable(                         &
      & varname='VelZ', dims=(/'x','y','z','t'/),    &
      & longname='vertical velocity',                &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !無次元圧力
    call HistoryAddVariable(                         &
      & varname='Exner', dims=(/'x','y','z','t'/),   &
      & longname='nondimensional pressure',          &
      & units='1',                                   &
      & xtype='double', history=rstat )
    
    !温位の擾乱
    call HistoryAddVariable(                         &
      & varname='PTemp', dims=(/'x','y','z','t'/),   &
      & longname='virtual potential temperature',    &
      & units='K',                                   &
      & xtype='double', history=rstat )
    
    !渦粘性係数
    call HistoryAddVariable(                         &
      & varname='Km', dims=(/'x','y','z','t'/),      &
      & longname='Km',                               &
      & units='m2.s-1',                              &
      & xtype='double', history=rstat )
    
    !渦粘性係数
    call HistoryAddVariable(                         &
      & varname='Kh', dims=(/'x','y','z','t'/),      &
      & longname='Kh',                               &
      & units='m2.s-1',                              &
      & xtype='double', history=rstat )

    !雲密度
    call HistoryAddVariable(                         &
      & varname='CDens', dims=(/'x','y','z','t'/),   &
      & longname='CDens',                            &
      & units='kg.m-3',                              &
      & xtype='double', history=rstat )
    
    !混合比
    call HistoryAddVariable(                         &
      & varname='QMix', dims=(/'x','y','z','s','t'/),&
      & longname='Mixing Ratio',                     &
      & units='kg.kg-1"',                            & 
      & xtype='double', history=rstat )

!    !IH1998
!    call HistoryAddVariable(                         &
!      & varname='H2O', dims=(/'x','y','z','t'/),     &
!      & longname='Number Density of H2O',            &
!      & units='m-3',                                 &
!      & xtype='double', history=rstat )

!    !IH1998
!    call HistoryAddVariable(                         &
!      & varname='H2SO4', dims=(/'x','y','z','t'/),   &
!      & longname='Number Density of H2SO4',          &
!      & units='m-3',                                 &
!      & xtype='double', history=rstat )
    
  end subroutine ReStartFileio_Init
  
  
  subroutine ReStartFileIO_Finalize
    !
    !リスタートファイルのクローズ
    !

    !モジュール読み込み
    use gtool_history, only: HistoryClose
    
    !暗黙の型宣言禁止
    implicit none
    
    !ファイルを閉じる
    call HistoryClose(rstat, quiet=.true.)
    
  end subroutine ReStartFileIO_Finalize
  

  subroutine ReStartFileio_Var_Get(   &
    & pyz_VelXB,  pyz_VelXN,    & ! (out)
    & xqz_VelYB,  xqz_VelYN,    & ! (out)
    & xyr_VelZB,  xyr_VelZN,    & ! (out)
    & xyz_PTempB, xyz_PTempN,   & ! (out)
    & xyz_ExnerB, xyz_ExnerN,   & ! (out)
    & xyzf_QMixB, xyzf_QMixN,   & ! (out)
    & xyz_KmB,    xyz_KmN,      & ! (out)
    & xyz_KhB,    xyz_KhN,      & ! (out)
    & xyz_CDensB, xyz_CDensN  )   ! (out)
!    & xyz_H2SO4b, xyz_H2SO4n,   & ! (out) !IH1998
!    & xyz_H2Ob,   xyz_H2On    )   ! (out) !IH1998
    !
    !リスタートファイルから情報取得
    !

    !モジュール読み込み
    use dc_message,    only : MessageNotify
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use dc_string,     only : toChar
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!配列サイズ (X 方向)
      &                       jmin, jmax,    &!配列サイズ (Y 方向)
      &                       kmin, kmax,    &!配列サイズ (Z 方向)
      &                       ncmax           !凝縮成分の数
    use timeset,       only : DelTimeLong,   &
      &                       RestartTime 

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out) :: pyz_VelXN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelYN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_ExnerN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTempN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KmN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KhN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDensN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMixN &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(out) :: pyz_VelXB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelYB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_ExnerB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTempB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KmB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KhB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDensB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMixB &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)

!    !IH1998
!    real(DP), intent(out) :: xyz_H2SO4b(imin:imax,jmin:jmax,kmin:kmax)
!    real(DP), intent(out) :: xyz_H2SO4n(imin:imax,jmin:jmax,kmin:kmax)
!    real(DP), intent(out) :: xyz_H2Ob(imin:imax,jmin:jmax,kmin:kmax)
!    real(DP), intent(out) :: xyz_H2On(imin:imax,jmin:jmax,kmin:kmax)

    character(STRING)    :: name               !変数名
    character(STRING)    :: TimeN = ""
    character(STRING)    :: TimeB = ""

    real(DP)  :: RTimeN, RTimeB, RTimeRtnB, RTimeRtnN
    logical   :: flag_time_exist
    logical   :: flag_time_err

    !-------------------------------------------------------------
    ! Get a Value from netCDF File
    !-------------------------------------------------------------    
    RTimeN = RestartTime
    if (RestartTime /= 0.0d0) then
     RTimeB = RestartTime - DelTimeLong 
    else 
     RTimeB = RestartTime 
    end if

    timeB = 't=' // toChar( RTimeB )
    timeN = 't=' // toChar( RTimeN )
    
!!!
!!! 確認
!!!
    name = "PTemp"
    call HistoryGet( InputFile, name, xyz_PTempB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI, &
      & returned_time = RTimeRtnB, flag_time_exist = flag_time_exist, err = flag_time_err )
    call HistoryGet( InputFile, name, xyz_PTempN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI, &
      & returned_time = RTimeRtnN, flag_time_exist = flag_time_exist, err = flag_time_err )

    if ( .NOT. flag_time_exist .OR. flag_time_err ) then 
     call MessageNotify( "M", "restartfileio", "flag_time_exist = %b", L=(/flag_time_exist/) )
     call MessageNotify( "M", "restartfileio", "flag_time_err   = %b", L=(/flag_time_err/)   )
!     call MessageNotify( "E", "restartfileio", "restart error occure" )
    end if

    if ( RTimeRtnB - RTimeB /= 0.0d0 .OR. RTimeRtnN - RTimeN /= 0.0d0 ) then 
      call MessageNotify( "M", "restartfileio", "RTimeRtnB  = %f", d=(/RTimeRtnB/)  )
      call MessageNotify( "M", "restartfileio", "RTimeRtnN  = %f", d=(/RTimeRtnN/)  )
      call MessageNotify( "M", "restartfileio", "RTimeB     = %f", d=(/RTimeB/) )
      call MessageNotify( "M", "restartfileio", "RTimeN     = %f", d=(/RTimeN/) )
!      call MessageNotify( "E", "restartfileio", "restart time is incorrect " )
    end if


!!!
!!! ファイルオープン & 値の取り出し
!!!
    name = "PTemp"
    call HistoryGet( InputFile, name, xyz_PTempB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_PTempN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelX"
    call HistoryGet( InputFile, name, pyz_VelXB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, pyz_VelXN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelY"
    call HistoryGet( InputFile, name, xqz_VelYB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xqz_VelYN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelZ"
    call HistoryGet( InputFile, name, xyr_VelZB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyr_VelZN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Exner"
    call HistoryGet( InputFile, name, xyz_ExnerB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_ExnerN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Km"
    call HistoryGet( InputFile, name, xyz_KmB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_KmN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
    
    name = "Kh"
    call HistoryGet( InputFile, name, xyz_KhB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_KhN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "CDens"
    call HistoryGet( InputFile, name, xyz_CDensB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_CDensN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
      
    name = "QMix"
    call HistoryGet( InputFile, name, xyzf_QMixB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyzf_QMixN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

!    !IH1998
!    name = "H2SO4"
!    call HistoryGet( InputFile, name, xyz_H2SO4b, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
!    call HistoryGet( InputFile, name, xyz_H2SO4n, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

!    !IH1998
!    name = "H2O"
!    call HistoryGet( InputFile, name, xyzf_H2Ob, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
!    call HistoryGet( InputFile, name, xyzf_H2On, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

  end subroutine ReStartFileio_Var_Get
       


  subroutine ReStartFileio_BasicZ_Get()
    !
    !リスタートファイルから情報取得
    !

    !モジュール読み込み
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!配列サイズ (X 方向)
      &                       jmin, jmax,    &!配列サイズ (Y 方向)
      &                       kmin, kmax,    &!配列サイズ (Z 方向)
      &                       ncmax           !凝縮成分の数
    use basicset,      only : basicset_init 

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var3D &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: Var4D &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)              :: xyz_ExnerBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場のエクスナー関数
    real(DP)              :: xyz_DensBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の密度
    real(DP)              :: xyz_PTempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温位
    real(DP)              :: xyz_VelSoundBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の音速
    real(DP)              :: xyz_PressBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の圧力
    real(DP)              :: xyz_TempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温度
    real(DP)              :: xyzf_QMixBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !基本場の混合比
    real(DP)              :: xyz_EffMolWtBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の分子量効果
    character(STRING)    :: name               !変数名

    !-------------------------------------------------------------
    ! 基本場の取得
    !
    name = "DensBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_DensBZ = Var3D
    
    name = "ExnerBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_ExnerBZ = Var3D
    
    name = "PTempBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTempBZ = Var3D
    
    name = "VelSoundBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_VelSoundBZ = Var3D
    
    name = "TempBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_TempBZ = Var3D
    
    name = "PressBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PressBZ = Var3D
    
    name = "QMixBZ"
    call HistoryGet( InputFile, name, Var4D, flag_mpi_split = FLAG_LIB_MPI )
    xyzf_QMixBZ = Var4D
    
    name = "EffMolWtBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_EffMolWtBZ = Var3D
    
    !----------------------------------------------------------
    ! BasicSet モジュールに値を設定
    !----------------------------------------------------------
    call basicset_init( &
      & xyz_PressBZ,    &
      & xyz_ExnerBZ,    &
      & xyz_TempBZ,     &
      & xyz_PTempBZ,    &
      & xyz_DensBZ,     &
      & xyz_VelSoundBZ, &
      & xyzf_QMixBZ,    &
      & xyz_EffMolWtBZ  &
      & )
    
  end subroutine ReStartFileio_BasicZ_Get
       
end module ReStartFileIO
