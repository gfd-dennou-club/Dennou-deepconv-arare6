
module ReStartFileIO_IH1998
  !
  !リスタート用の場の情報を netCDF ファイルに出力するためのルーチン
  ! * Imamura and Hashimoto (1998) に基づく計算ルーチンについては, 
  !   数密度で計算しているので, リスタートファイルを別に分ける.
  !

  !モジュール読み込み
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING

  !暗黙の型宣言禁止
  implicit none
  
  !関数を public に指定
  public ReStartFileIO_IH1998_init
  public ReStartFileio_IH1998_Finalize
  public ReStartFileio_IH1998_Var_Get

  integer, private, save           :: ncmax = 3
  type(GT_HISTORY), save, public   :: rstat2
  character(STRING), save, private :: InputFile2
  character(STRING), save, private :: OutputFile2
  
contains

  subroutine ReStartFileio_IH1998_Init ( FlagInitData )
    !
    !リスタートファイルの初期化
    !

    !モジュール読み込み
    use gtool_history, only : HistoryCreate,  &
      &                       HistoryPut,     &
      &                       HistoryAddVariable
    use dc_message,    only : MessageNotify
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use axesset,       only : x_X,            &!X 座標軸(スカラー格子点)
      &                       y_Y,            &!Y 座標軸(スカラー格子点)
      &                       z_Z              !Z 座標軸(スカラー格子点)
    use fileset,       only : filetitle,      &!データの表題
      &                       filesource,     &!データを作成する手順
      &                       FileInstitution  !最終変更者・組織
    use restartfileio, only : Inputfile, Outputfile  !リスタートファイル名(デフォルト)
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    logical, intent(in), optional :: FlagInitData  ! 初期値生成プログラムの場合は .true.

    real(4)            :: SpcID(ncmax)
    integer            :: N, L, M
    integer            :: s    
    integer            :: num

    !入力ファイル名の設定
    if (Inputfile /= "") then 
      num = len_trim( Inputfile ) - 3
      Inputfile2 = Inputfile(1:num) // '_IH1998.nc'
      call MessageNotify( "M", &
        & "restartfileioIO_IH1998_init", "INPUT = %c", c1=trim(InputFile2) )
    end if

    !出力ファイル名の設定
    if (Outputfile /= "") then 
      num = len_trim( Outputfile ) - 3
      Outputfile2 = Outputfile(1:num) // '_IH1998.nc'
      call MessageNotify( "M", &
        & "restartfileioIO_IH1998_init", "OUTPUT = %c", c1=trim(OutputFile2) )
    end if

    !パラメタ
    N = size(x_X, 1)
    L = size(y_Y, 1)
    M = size(z_Z, 1)
    SpcID = 0.0d0
    do s = 1, ncmax
      SpcID(s) = real( s, 4 )
    end do    
    
    !-------------------------------------------------------------    
    ! ヒストリー作成
    !-------------------------------------------------------------  
    call HistoryCreate(                              &
      & file = Outputfile2,                           &
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
      & history=rstat2, quiet=.true. )
    
    !-------------------------------------------------------------  
    ! 変数出力
    !-------------------------------------------------------------
    call HistoryPut('x', x_X, rstat2 )
    call HistoryPut('y', y_Y, rstat2 )
    call HistoryPut('z', z_Z, rstat2 )
    call HistoryPut('s', real(SpcID, 4), rstat2 )

    !H2SO4 数密度
    call HistoryAddVariable(                            &
      & varname='NDens1', dims=(/'x','y','z','s','t'/), &
      & longname='Number Density of H2SO4 (all)',       &
      & units='m-3',                                    &
      & xtype='float', history=rstat2 )

    !H2O 数密度
    call HistoryAddVariable(                            &
      & varname='NDens2', dims=(/'x','y','z','s','t'/), &
      & longname='Number Density of H2O (all)',         &
      & units='m-3',                                    &
      & xtype='float', history=rstat2 )

  end subroutine ReStartFileio_IH1998_Init
  
  
  subroutine ReStartFileIO_IH1998_Finalize
    !
    !リスタートファイルのクローズ
    !

    !モジュール読み込み
    use gtool_history, only: HistoryClose
    
    !暗黙の型宣言禁止
    implicit none
    
    !ファイルを閉じる
    call HistoryClose(rstat2, quiet=.true.)
    
  end subroutine ReStartFileIO_IH1998_Finalize
  

  subroutine ReStartFileio_IH1998_Var_Get(  &
    &   xyzf_NDens1B, xyzf_NDens1N,         & ! (out)
    &   xyzf_NDens2B, xyzf_NDens2N          & ! (out)
    & )
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
      &                       kmin, kmax      !配列サイズ (Z 方向)
    use timeset,       only : DelTimeLong,   &
      &                       RestartTime 
    use restartfileio, only : TimeN, TimeB
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out) :: xyzf_NDens1B(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out) :: xyzf_NDens1N(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out) :: xyzf_NDens2B(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    real(DP), intent(out) :: xyzf_NDens2N(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
    character(STRING)     :: name

    ! ファイルオープン & 値の取り出し
    name = "NDens1"
    call HistoryGet( InputFile2, name, xyzf_NDens1B, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile2, name, xyzf_NDens1N, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "NDens2"
    call HistoryGet( InputFile2, name, xyzf_NDens2B, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile2, name, xyzf_NDens2N, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

  end subroutine ReStartFileio_IH1998_Var_Get
       
end module ReStartFileIO_IH1998
