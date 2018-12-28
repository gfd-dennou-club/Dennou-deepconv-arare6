!= Module Arare4InitFileio
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: arare4fileio.f90,v 1.12 2014/03/04 05:55:04 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module Arare4InitFileio
  !
  ! deepconv/arare4 のリスタートファイルから t = 0 での初期値を
  ! 作成するためのモジュール
  !

  !モジュール読み込み
  use dc_types,  only : STRING

  !暗黙の型宣言禁止
  implicit none
  
  !関数を public に指定
  public Arare4InitFileio_init
  public Arare4InitFileio_BZ_Get
  public Arare4InitFileio_Var_Get
  public Arare4InitFileio_MMC_BZ_Get
  public Arare4InitFileio_MMC_Var_Get

  character(STRING), save, private :: Arare4Prefix  = "arare4_"
  
contains

  subroutine Arare4InitFileio_Init
    !
    !リスタートファイル名の取得
    !

    !モジュール読み込み
    use dc_types,      only : DP
    use dc_message,    only : MessageNotify
    use dc_iounit,     only : FileOpen
    use namelist_util, only : namelist_filename

    !暗黙の型宣言禁止
    implicit none

    !変数
    integer            :: unit     !装置番号
     
    !NAMELIST から情報を取得
    NAMELIST /arare4fileio_nml/ Arare4Prefix
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=arare4fileio_nml)
    close(unit)
    
    !確認
    call MessageNotify( "M", &
      & "arare4FileIO_init", "Arare4Prefix = %c", c1=trim(Arare4Prefix) )

  end subroutine Arare4InitFileio_Init
  
  
  subroutine Arare4InitFileio_Var_Get(   &
    & pyz_VelX,     & ! (out)
    & xqz_VelY,     & ! (out)
    & xyr_VelZ,     & ! (out)
    & xyz_PTemp,    & ! (out)
    & xyz_Exner,    & ! (out)
    & xyzf_QMix,    & ! (out)
    & xyz_Km,       & ! (out)
    & xyz_Kh,       & ! (out)
    & xyz_CDens   )   ! (out)
    !
    !リスタートファイルから情報取得
    !

    !モジュール読み込み
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use dc_string,     only : toChar
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!配列サイズ (X 方向)
      &                       jmin, jmax,    &!配列サイズ (Y 方向)
      &                       kmin, kmax,    &!配列サイズ (Z 方向)
      &                       ncmax,         &!凝縮成分の数
      &                       nx, nz          !物理領域の大きさ
    use composition,   only : SpcWetSymbol
    use setmargin,     only : SetMargin_xyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_xyr, &
      &                       SetMargin_xyzf

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var3D(1:nx, 1, 1:nz)
    real(DP), intent(out) :: pyz_VelX &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelY &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZ &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Exner &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTemp &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Km &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Kh &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDens &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMix &
      &                      (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)             :: RestartTime
    integer              :: f

    character(STRING)    :: name               !変数名
    character(STRING)    :: Time
    character(STRING)    :: InputFile

    !-------------------------------------------------------------
    ! Get a Value from netCDF File
    !-------------------------------------------------------------    
    RestartTime = 0.0d0
    time = 't=' // trim(toChar( RestartTime ))

    name = "PotTemp"
!    name1 = "PotTempDist"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTemp(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_PTemp )

    name = "VelX"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    pyz_VelX(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_pyz( pyz_VelX )

    name = "VelY"
    xqz_VelY = 0.0d0

    name = "VelZ"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyr_VelZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyr( xyr_VelZ )

    name = "Exner"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Exner(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_Exner )

    name = "Km"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Km(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_Km )

    name = "Kh"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Kh(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_Kh )
      
!    name = "MixRt"
!    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
!    call HistoryGet( InputFile, name, Var4D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
!    xyzf_QMix(1:nx,1,1:nz,1:ncmax) = Var4D(1:nx,1,1:nz,1:ncmax)
    do f = 1, ncmax
      name = trim(SpcWetSymbol(f))
      InputFile = trim(Arare4Prefix) // trim(SpcWetSymbol(f)) // ".nc"
      call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
      xyzf_QMix(1:nx,1,1:nz,f) = Var3D(1:nx,1,1:nz)
    end do

    call SetMargin_xyzf( xyzf_QMix )

    name = "CloudDensity"    
    xyz_CDens = 0.0d0

  end subroutine Arare4InitFileio_Var_Get
       

  subroutine Arare4InitFileio_BZ_Get(   &
    & xyz_PressBZ,    & ! (out)
    & xyz_ExnerBZ,    & ! (out)
    & xyz_TempBZ,     & ! (out)
    & xyz_PTempBZ,    & ! (out)
    & xyz_DensBZ,     & ! (out)
    & xyz_VelSoundBZ, & ! (out)
    & xyzf_QMixBZ,    & ! (out)
    & xyz_EffMolWtBZ  & ! (out)
    & )
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
      &                       ncmax,         &!凝縮成分の数
      &                       nx, nz          !物理領域のサイズ
    use composition,   only : SpcWetSymbol
    use setmargin,     only : SetMargin_xyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_xyr, &
      &                       SetMargin_xyzf

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var3D(1:nx,1,1:nz)
    real(DP), intent(out) :: xyz_ExnerBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場のエクスナー関数
    real(DP), intent(out) :: xyz_DensBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の密度
    real(DP), intent(out) :: xyz_PTempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温位
    real(DP), intent(out) :: xyz_VelSoundBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の音速
    real(DP), intent(out) :: xyz_PressBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の圧力
    real(DP), intent(out) :: xyz_TempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温度
    real(DP), intent(out) :: xyzf_QMixBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !基本場の混合比
    real(DP), intent(out) :: xyz_EffMolWtBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の分子量効果
    character(STRING)    :: name               !変数名
    character(STRING)    :: InputFile
    integer              :: f
 
    !-------------------------------------------------------------
    ! 基本場の取得
    !
    InputFile = trim(Arare4Prefix) // "BasicZ.nc"

    name = "DensBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_DensBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_DensBZ )

    name = "ExnerBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_ExnerBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_ExnerBZ )
    
    name = "PotTempBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTempBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_PTempBZ )
    
    name = "VelSoundBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_VelSoundBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_VelSoundBZ )
    
    name = "TempBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_TempBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_TempBZ )
    
    name = "PressBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_PressBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_PressBZ)

!    name = "MixRtBasicZ"
!    call HistoryGet( InputFile, name, Var4D(:,1,:,:), flag_mpi_split = FLAG_LIB_MPI )
!    xyzf_QMixBZ(1:nx,1,1:nz,1:ncmax) = Var4D(1:nx,1,1:nz,1:ncmax)

    do f = 1, ncmax
      name = trim(SpcWetSymbol(f))//"BasicZ"
      call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
      xyzf_QMixBZ(1:nx,1,1:nz,f) = Var3D(1:nx,1,1:nz)
    end do
    call SetMargin_xyzf( xyzf_QMixBZ )
   
    name = "EffMolWtBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_EffMolWtBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_EffMolWtBZ )

  end subroutine Arare4InitFileio_BZ_Get


  subroutine Arare4InitFileio_MMC_Var_Get(   &
    & pyz_VelX,     & ! (out)
    & xqz_VelY,     & ! (out)
    & xyr_VelZ,     & ! (out)
    & xyz_PTemp,    & ! (out)
    & xyz_Exner,    & ! (out)
    & xyzf_QMix,    & ! (out)
    & xyz_Km,       & ! (out)
    & xyz_Kh,       & ! (out)
    & xyz_CDens   )   ! (out)
    !
    !リスタートファイルから情報取得
    !

    !モジュール読み込み
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use dc_string,     only : toChar
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!配列サイズ (X 方向)
      &                       jmin, jmax,    &!配列サイズ (Y 方向)
      &                       kmin, kmax,    &!配列サイズ (Z 方向)
      &                       ncmax,         &!凝縮成分の数
      &                       nx, nz          !物理領域のサイズ
    use setmargin,     only : SetMargin_xyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_xyr, &
      &                       SetMargin_xyzf

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var3D(1:nx, 1, 1:nz)
!    real(DP)              :: Var4D(1:nx, 1, 1:nz, 1:ncmax)
    real(DP), intent(out) :: pyz_VelX &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelY &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZ &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Exner &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTemp &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Km &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_Kh &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDens &
      &                      (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMix &
      &                      (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)             :: RestartTime

    character(STRING)    :: name               !変数名
    character(STRING)    :: name1
    character(STRING)    :: Time
    character(STRING)    :: InputFile

    !-------------------------------------------------------------
    ! Get a Value from netCDF File
    !-------------------------------------------------------------    
    RestartTime = 0.0d0
    time = 't=' // trim(toChar( RestartTime ))

    name = "PotTemp"
    name1 = "PotTempDist"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name1, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTemp(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_PTemp )

    name = "VelX"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    pyz_VelX(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_pyz( pyz_VelX )

    name = "VelY"
    xqz_VelY = 0.0d0

    name = "VelZ"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyr_VelZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyr( xyr_VelZ )

    name = "Exner"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Exner(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_Exner )

    name = "Km"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Km(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_Km )

    name = "Kh"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_Kh(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_Kh )

    name = "DensCloud"
    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
    call HistoryGet( InputFile, name, Var3D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    xyz_CDens(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_CDens )
      
!    name = "MixRt"
!    InputFile = trim(Arare4Prefix) // trim(name) // ".nc"
!    call HistoryGet( InputFile, name, Var4D, range=Time, flag_mpi_split = FLAG_LIB_MPI )
!    xyzf_QMix(1:nx,1,1:nz,1:ncmax) = Var4D(1:nx,1,1:nz,1:ncmax)
    xyzf_QMix(1:nx,1,1:nz,1:ncmax) = 0.0d0
    
  end subroutine Arare4InitFileio_MMC_Var_Get
       

  subroutine Arare4InitFileio_MMC_BZ_Get(   &
    & xyz_PressBZ,   &
    & xyz_ExnerBZ,   &
    & xyz_TempBZ,    &
    & xyz_PTempBZ,   &
    & xyz_DensBZ,    &
    & xyz_VelSoundBZ,&
    & xyzf_QMixBZ,   &
    & xyz_EffMolWtBZ &
    & )
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
      &                       ncmax,         &!凝縮成分の数
      &                       nx, nz          !物理領域のサイズ
    use setmargin,     only : SetMargin_xyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_pyz, &
      &                       SetMargin_xyr, &
      &                       SetMargin_xyzf

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var3D(1:nx,1,1:nz)
!    real(DP)              :: Var4D(1:nx,1,1:nz,1:ncmax)
    real(DP), intent(out) :: xyz_ExnerBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場のエクスナー関数
    real(DP), intent(out) :: xyz_DensBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の密度
    real(DP), intent(out) :: xyz_PTempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温位
    real(DP), intent(out) :: xyz_VelSoundBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の音速
    real(DP), intent(out) :: xyz_PressBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の圧力
    real(DP), intent(out) :: xyz_TempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の温度
    real(DP), intent(out) :: xyzf_QMixBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !基本場の混合比
    real(DP), intent(out) :: xyz_EffMolWtBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !基本場の分子量効果
    character(STRING)    :: name               !変数名
    character(STRING)    :: InputFile
 
    !-------------------------------------------------------------
    ! 基本場の取得
    !
    InputFile = trim(Arare4Prefix) // "BasicZ.nc"

    name = "DensBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_DensBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_DensBZ )

    name = "ExnerBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_ExnerBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_ExnerBZ )
    
    name = "PotTempBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTempBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_PTempBZ )
    
    name = "VelSoundBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_VelSoundBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_VelSoundBZ )
    
    name = "TempBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_TempBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_TempBZ )
    
    name = "PressBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_PressBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)
    
    call SetMargin_xyz( xyz_PressBZ)

!    name = "MixRtBasicZ"
!    call HistoryGet( InputFile, name, Var4D(:,1,:,:), flag_mpi_split = FLAG_LIB_MPI )
!    xyzf_QMixBZ(1:nx,1,1:nz) = Var4D(1:nx,1,1:nz)
     xyzf_QMixBZ = 0.0d0
    
    name = "EffMolWtBasicZ"
    call HistoryGet( InputFile, name, Var3D(:,1,:), flag_mpi_split = FLAG_LIB_MPI )
    xyz_EffMolWtBZ(1:nx,1,1:nz) = Var3D(1:nx,1,1:nz)

    call SetMargin_xyz( xyz_EffMolWtBZ )

  end subroutine Arare4InitFileio_MMC_BZ_Get
       
end module Arare4InitFileio

