!= Module Arare4ReStartFileIO
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: restartfileio.f90,v 1.17 2014/07/08 00:59:55 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module Arare4ReStartFileIO
  !
  !deepconv/arare4 のリスタートファイルを読み込むためのモジュール.
  !

  !モジュール読み込み
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING
 
  !暗黙の型宣言禁止
  implicit none
  
  !関数を public に指定
  public Arare4ReStartFileio_BasicZ_Get
  public Arare4ReStartFileio_Var_Get
  
  type(GT_HISTORY),  save, public  :: rstat
  character(STRING), save, private :: InitialFile = ""
  character(STRING), save, private :: InputFile   = ""
  character(STRING), save, private :: OutputFile  = "output.nc"
  
contains
  
  subroutine Arare4ReStartFileio_Var_Get(   &
    & pyz_VelXB,  pyz_VelXN,    & ! (out)
    & xqz_VelYB,  xqz_VelYN,    & ! (out)
    & xyr_VelZB,  xyr_VelZN,    & ! (out)
    & xyz_PTempB, xyz_PTempN,   & ! (out)
    & xyz_ExnerB, xyz_ExnerN,   & ! (out)
    & xyzf_QMixB, xyzf_QMixN,   & ! (out)
    & xyz_KmB,    xyz_KmN,      & ! (out)
    & xyz_KhB,    xyz_KhN,      & ! (out)
    & xyz_CDensB, xyz_CDensN  )   ! (out)

    !
    !リスタートファイルから情報取得
    !

    use gtool_history, only : HistoryGet
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use dc_types,      only : DP, STRING
    use dc_string,     only : toChar
    use gridset,       only : imin, imax, &!X 方向の配列サイズ
      &                       jmin, jmax, &!Y 方向の配列サイズ
      &                       kmin, kmax, &!Z 方向の配列サイズ
      &                       ncmax
    use timeset,       only : DelTimeLong, &
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
      &              (imin:imax,jmin:jmax,kmin:kmax,ncmax)
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
      &              (imin:imax,jmin:jmax,kmin:kmax,ncmax)

    ! deepconv では配列のサイズが異なるので調整する
    !
    real(DP)  :: pz_VelXN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xr_VelZN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_ExnerN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_PTempN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_KmN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_KhN(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xzf_QMixN(imin-1:imax,kmin-1:kmax,ncmax)
    real(DP)  :: pz_VelXB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xr_VelZB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_ExnerB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_PTempB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_KmB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xz_KhB(imin-1:imax,kmin-1:kmax)
    real(DP)  :: xzf_QMixB(imin-1:imax,kmin-1:kmax,ncmax)

    integer    :: j

    character(STRING)    :: name               !変数名
    character(STRING)    :: TimeN = ""
    character(STRING)    :: TimeB = ""

    !-------------------------------------------------------------
    ! Get a Value from netCDF File
    !-------------------------------------------------------------    
    if (RestartTime /= 0.0d0) then
      timeB = 't=' // toChar( RestartTime - DelTimeLong)
      timeN = 't=' // toChar( RestartTime  )
    else
      timeB = 't=' // toChar( RestartTime )
      timeN = 't=' // toChar( RestartTime )
    end if

    name = "PotTemp"
    call HistoryGet( InputFile, name, xz_PTempB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xz_PTempN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelX"
    call HistoryGet( InputFile, name, pz_VelXB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, pz_VelXN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

!    name = "VelY"
!    call HistoryGet( InputFile, name, xqz_VelYB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
!    call HistoryGet( InputFile, name, xqz_VelYN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
    xqz_VelYB = 0.0d0
    xqz_VelYN = 0.0d0

    name = "VelZ"
    call HistoryGet( InputFile, name, xr_VelZB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xr_VelZN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Exner"
    call HistoryGet( InputFile, name, xz_ExnerB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xz_ExnerN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Km"
    call HistoryGet( InputFile, name, xz_KmB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xz_KmN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
    
    name = "Kh"
    call HistoryGet( InputFile, name, xz_KhB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xz_KhN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

!    name = "CDens"
!    call HistoryGet( InputFile, name, xyz_CDensB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
!    call HistoryGet( InputFile, name, xyz_CDensN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
    xyz_CDensB = 0.0d0
    xyz_CDensN = 0.0d0
      
    name = "MixRt"
    call HistoryGet( InputFile, name, xzf_QMixB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xzf_QMixN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    do j = jmin, jmax
      xyz_PTempB(imin:imax,j,kmin:kmax) = xz_PTempB(imin:imax,kmin:kmax)
      xyz_PTempN(imin:imax,j,kmin:kmax) = xz_PTempN(imin:imax,kmin:kmax)
      pyz_VelXB(imin:imax,j,kmin:kmax) = pz_VelXB(imin:imax,kmin:kmax)
      pyz_VelXN(imin:imax,j,kmin:kmax) = pz_VelXN(imin:imax,kmin:kmax)
      xyr_VelZB(imin:imax,j,kmin:kmax) = xr_VelZB(imin:imax,kmin:kmax)
      xyr_VelZN(imin:imax,j,kmin:kmax) = xr_VelZN(imin:imax,kmin:kmax)
      xyz_ExnerB(imin:imax,j,kmin:kmax) = xz_ExnerB(imin:imax,kmin:kmax)
      xyz_ExnerN(imin:imax,j,kmin:kmax) = xz_ExnerN(imin:imax,kmin:kmax)
      xyz_KmB(imin:imax,j,kmin:kmax) = xz_KmB(imin:imax,kmin:kmax)
      xyz_KmN(imin:imax,j,kmin:kmax) = xz_KmN(imin:imax,kmin:kmax)
      xyz_KhB(imin:imax,j,kmin:kmax) = xz_KhB(imin:imax,kmin:kmax)
      xyz_KhN(imin:imax,j,kmin:kmax) = xz_KhN(imin:imax,kmin:kmax)
      xyzf_QMixB(imin:imax,j,kmin:kmax,:) = xzf_QMixB(imin:imax,kmin:kmax,:)
      xyzf_QMixN(imin:imax,j,kmin:kmax,:) = xzf_QMixN(imin:imax,kmin:kmax,:)
    end do

  end subroutine Arare4ReStartFileio_Var_Get
       


  subroutine Arare4ReStartFileio_BasicZ_Get()
    !
    !リスタートファイルから情報取得
    !

    !モジュール読み込み
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax, &!X 方向の配列サイズ
      &                       jmin, jmax, &!Y 方向の配列サイズ
      &                       kmin, kmax, &!Z 方向の配列サイズ
      &                       ncmax
    use basicset,      only : basicset_init 

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)              :: Var2D &
      &                     (imin-1:imax,kmin-1:kmax)
    real(DP)              :: Var3D &
      &                     (imin-1:imax,kmin-1:kmax,ncmax)
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
    integer              :: j
    character(STRING)    :: name               !変数名

    !-------------------------------------------------------------
    ! 基本場の取得
    !
    name = "DensBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_DensBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do

    name = "ExnerBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_ExnerBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do
    
    name = "PotTempBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_PTempBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do
    
    name = "VelSoundBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_VelSoundBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do

    name = "TempBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_TempBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do
    
    name = "PressBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_PressBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do
    
    name = "MixRtBasicZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyzf_QMixBZ(imin:imax,j,kmin:kmax,:) = Var3D(imin:imax,kmin:kmax,:)
    end do
    
    name = "EffMolWtBasicZ"
    call HistoryGet( InputFile, name, Var2D, flag_mpi_split = FLAG_LIB_MPI )
    do j = jmin, jmax
      xyz_EffMolWtBZ(imin:imax,j,kmin:kmax) = Var2D(imin:imax,kmin:kmax)
    end do

!    write(*,*) "PressBZ", minval(xyz_PressBZ), maxval(xyz_PressBZ)
!    write(*,*) "ExnerBZ", minval(xyz_ExnerBZ), maxval(xyz_ExnerBZ)
!    write(*,*) "TempBZ",  minval(xyz_TempBZ), maxval(xyz_TempBZ)
!    write(*,*) "PTempBZ", minval(xyz_PTempBZ), maxval(xyz_PTempBZ)
!    write(*,*) "DensBZ",  minval(xyz_DensBZ), maxval(xyz_DensBZ)
!    write(*,*) "VelSnd",  minval(xyz_VelSoundBZ), maxval(xyz_VelSoundBZ)
!    write(*,*) "QMixBZ", minval(xyzf_QMixBZ), maxval(xyzf_QMixBZ)
!    write(*,*) "EffMolWtBZ", minval(xyz_EffMolWtBZ), maxval(xyz_EffMolWtBZ)

    
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
    
  end subroutine Arare4ReStartFileio_BasicZ_Get
       
end module Arare4ReStartFileIO
