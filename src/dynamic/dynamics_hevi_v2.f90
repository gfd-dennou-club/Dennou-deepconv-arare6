!= Module DynamicsHEVI
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro), ODAKA Masatsugu 
! Version::   $Id: dynamics_hevi_v2.f90,v 1.8 2014/07/08 00:58:06 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module DynamicsHEVI_v2
  !
  ! 力学コア. 
  !   タイムスプリット法を利用. 音波モードとそれ以外を別々の時間刻みで解く. 
  !   短い時間ステップの計算には HE-VI 法を利用.
  !
  !   微分平均演算モジュールを使わずに書き下した版. 
  !   計算速度は v1 よりも 2 倍以上早いが, デバッグはしにくい. 
  !
  ! Note: 
  !  * エクスナー関数の空間方向の離散化は 2 次精度であるため, 気圧傾度
  !    力項の計算プログラムにおいて differentiate_center4 モジュールを
  !    指定することはできないので注意.
  !

  !モジュール読み込み
  use dc_types,   only : DP

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private
  
  real(DP), save,private :: beta  = 5.0d-1     !クランクニコルソン法なら 0.5
                                               !完全陰解法なら 1 
  integer, save,private  :: N = 10             !係数行列/改行列の次数, 整合寸法
  integer, save,private  :: M = 10             !方程式の組数
  integer, save,private  :: NUD = 1            !係数行列の上三角部分の帯幅
  integer, save,private  :: NLD = 1            !係数行列の下三角部分の帯幅
  integer, save,private  :: NAL = 1            !LU 分解の結果 L の整合寸法
  integer, save,private  :: NA = 3             !NUD + NLD + 1

  real(DP), allocatable, save,private :: xyz_F1BZ(:,:,:)
                                               !係数行列の計算に用いる配列
  real(DP), allocatable, save,private :: xyr_CpVPTempBZ(:,:,:)
                                               !係数行列の計算に用いる配列
  real(DP), allocatable, save,private :: xyr_CpDensVPTemp2BZ(:,:,:)
                                               !係数行列の計算に用いる配列
  real(DP), allocatable, save,private :: xyr_DensVPTempBZ(:,:,:) 
                                               !係数行列の計算に用いる配列

  real(DP), allocatable, save,private :: A(:)  !係数行列の対角成分
  real(DP), allocatable, save,private :: B(:)  !係数行列の上三角部分
  real(DP), allocatable, save,private :: C(:)  !係数行列の下三角部分
  real(DP), allocatable, save,private :: AL1(:)!LU 分解の結果 L (1 次元配列)
  integer,  allocatable, save,private :: IP(:) !部分ピボット交換の情報を格納

  real(DP), save,private :: AlphaH = 0.0d0     !音波減衰項の減衰係数 (水平方向)
  real(DP), save,private :: AlphaV = 0.0d0     !音波減衰項の減衰係数 (鉛直方向)
  real(DP), save,private :: NuHh   = 0.0d0     !熱に対する数値粘性の係数 (水平方向)
  real(DP), save,private :: NuVh   = 0.0d0     !熱に対する数値粘性の係数 (鉛直方向)
  real(DP), save,private :: NuHm   = 0.0d0     !運動量に対する数値粘性の係数 (水平方向)
  real(DP), save,private :: NuVm   = 0.0d0     !運動量に対する数値粘性の係数 (鉛直方向)

  character(*), parameter:: module_name = 'DynamicHEVI'
                                               ! モジュールの名称.
                                               ! Module name
  real(DP), save,private :: FactorBuoyTemp    = 1.0d0
                                               !浮力 (温度の寄与) の有無
                                               !考慮しない場合は値をゼロにする.
  real(DP), save,private :: FactorBuoyMolWt   = 1.0d0
                                               !浮力 (分子量効果) の有無
                                               !考慮しない場合は値をゼロにする.
  real(DP), save,private :: FactorBuoyLoading = 1.0d0
                                               !浮力 (荷重効果) の有無
                                               !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtAdv    = 1.0d0     !エクスナー関数の移流の有無
                                               !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtExpnd  = 1.0d0     !エクスナー関数の膨張項の有無
                                               !考慮しない場合は値をゼロにする.

  real(DP), allocatable, save,private :: xyr_Dummy(:,:,:) !ダミー

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:)   !基本場の鉛直微分
  real(DP), allocatable, save,private :: xyr_DExnerBZDz(:,:,:)   !基本場の鉛直微分
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:) !基本場の鉛直微分
  
  ! public 
  !
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Init
    !
    ! 力学コア 初期化ルーチン
    !

    !モジュール読み込み
    use dc_types,      only : DP
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use gridset,       only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
      &                       FlagCalc3D
    use timeset,       only : DelTimeShort, DelTimeLong
    use axesset,       only : dx, dy, dz        ! 格子間隔
    use namelist_util, only : namelist_filename
    use basicset,      only : xyz_PTempBZ,     &!基本場の温位
      &                       xyz_ExnerBZ,     &
      &                       xyzf_QMixBZ
    use differentiate_center4, &
      &                only : xyr_dz_xyz

    !暗黙の型宣言禁止
    implicit none
    
    real(DP)  :: AlphaSound = 5.0d-2  !音波減衰項の係数 (気象庁数値予報課報告・別冊49 より)
    real(DP)  :: AlphaNDiff = 1.0d-3  !4次の数値拡散の係数. CReSS マニュアルより
    real(DP)  :: NDiffRatio = 1.0d0   !速度に対する粘性を上げる場合は数字を 1 以上にする. 
    integer   :: unit                 !装置番号
    integer   :: f

    !-------------------------------------------------------------------
    ! Namelist から情報を取得する
    !
    NAMELIST /Dynamics_nml/                                    &
         & AlphaSound, AlphaNDiff, NDiffRatio, beta,           &
         & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
         & FactorDExnerDtAdv, FactorDExnerDtExpnd

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=dynamics_nml)
    close(unit)
   
    !-------------------------------------------------------------------
    ! 音波減衰項の減衰係数を決める
    ! 
    ! 気象庁予報課報告別冊 49 p53 に従い, 水平と鉛直とを分けて考える. 
    ! なお, 2 次元計算の場合には DelY に依存しないようにする. 
    !
    if ( FlagCalc3D ) then 
      AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
    else
      AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
      AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
    end if

    !-------------------------------------------------------------------
    ! 数値拡散係数を決める
    !
    ! CReSS マニュアルの記述に従って NuH, NuV を決める.
    ! 運動量と熱に対する数値拡散の大きさを変えられるように NDiffRatio を乗じている.
    ! 
    if ( FlagCalc3D ) then 
      NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    else
      NuHh = AlphaNDiff * ( dx ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    end if
    NuVh = AlphaNDiff * ( dz ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    NuHm = NuHh * NDiffRatio
    NuVm = NuVh * NDiffRatio

    !-------------------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", module_name, "AlphaH = %f", d=(/AlphaH/) )
    call MessageNotify( "M", module_name, "AlphaV = %f", d=(/AlphaV/) )
    call MessageNotify( "M", module_name, "NuHh = %f",   d=(/NuHh/)   )
    call MessageNotify( "M", module_name, "NuVh = %f",   d=(/NuVh/)   )
    call MessageNotify( "M", module_name, "NuHm = %f",   d=(/NuHm/)   )
    call MessageNotify( "M", module_name, "NuVm = %f",   d=(/NuVm/)   )
    call MessageNotify( "M", module_name, "FactorBuoyTemp = %f",    d=(/FactorBuoyTemp/)    )
    call MessageNotify( "M", module_name, "FactorBuoyMolWt = %f",   d=(/FactorBuoyMolWt/)   )
    call MessageNotify( "M", module_name, "FactorBuoyLoading = %f", d=(/FactorBuoyLoading/) )
    call MessageNotify( "M", module_name, "FactorDExnerDtAdv = %f",    d=(/FactorDExnerDtAdv/)    )
    call MessageNotify( "M", module_name, "FactorDExnerDtExpnd = %f",  d=(/FactorDExnerDtExpnd/)  )

    ! 鉛直陰解法を用いるために, 行列の準備を行う. 
    !
    call Dynamics_VI_init

    ! tendency の出力
    !
    call Dynamics_Tendency_output

    ! 配列の用意
    !
    allocate( xyr_Dummy(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DExnerBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_Dummy      = 0.0d0
    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
!    xyr_DExnerBZDz = xyr_dz_xyz( xyz_ExnerBZ )
    xyr_DExnerBZDz = 0.0d0
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Long_forcing(   &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,        & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
    & xyz_PTempBl, xyz_PTempNl,       & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
    & xyz_KmBl,    xyz_KmNl,          & ! (in)
    & pyz_DVelXDtNl,                  & ! (inout)
    & xqz_DVelYDtNl,                  & ! (inout)
    & xyr_DVelZDtNl,                  & ! (inout)
    & xyz_DPTempDtNl,                 & ! (inout)
    & xyz_DExnerDtNl,                 & ! (inout)
    & xyzf_DQMixDtNl,                 & ! (inout)
    & xyz_DKmDtNl                     & ! (inout)
    & )
    !
    ! 力学コア: 長い時間ステップで評価する項の計算.
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, ncmax, &
                                          !配列サイズ
      &                   FlagCalc3D      !フラグ

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(inout) :: xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax)


    !------------------------------------------------------
    ! 条件分岐
    !
    if ( FlagCalc3D ) then 
    
      call Dynamics3D_Long_forcing(       &
        & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
        & xqz_VelYBl,  xqz_VelYNl,        & ! (in)
        & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
        & xyz_PTempBl, xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
        & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
        & xyz_KmBl,    xyz_KmNl,          & ! (in)
        & pyz_DVelXDtNl,                  & ! (inout)
        & xqz_DVelYDtNl,                  & ! (inout)
        & xyr_DVelZDtNl,                  & ! (inout)
        & xyz_DPTempDtNl,                 & ! (inout)
        & xyz_DExnerDtNl,                 & ! (inout)
        & xyzf_DQMixDtNl,                 & ! (inout)
        & xyz_DKmDtNl                     & ! (inout)
        & )
     
    else

      call Dynamics2D_Long_forcing(       &
        & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
        & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
        & xyz_PTempBl, xyz_PTempNl,       & ! (in)
        & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
        & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
        & xyz_KmBl,    xyz_KmNl,          & ! (in)
        & pyz_DVelXDtNl,                  & ! (inout)
        & xqz_DVelYDtNl,                  & ! (inout)
        & xyr_DVelZDtNl,                  & ! (inout)
        & xyz_DPTempDtNl,                 & ! (inout)
        & xyz_DExnerDtNl,                 & ! (inout)
        & xyzf_DQMixDtNl,                 & ! (inout)
        & xyz_DKmDtNl                     & ! (inout)
        & )

    end if

  end subroutine Dynamics_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics3D_Long_forcing(   &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,        & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
    & xyz_PTempBl, xyz_PTempNl,       & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
    & xyz_KmBl,    xyz_KmNl,          & ! (in)
    & pyz_DVelXDtNl,                  & ! (inout)
    & xqz_DVelYDtNl,                  & ! (inout)
    & xyr_DVelZDtNl,                  & ! (inout)
    & xyz_DPTempDtNl,                 & ! (inout)
    & xyz_DExnerDtNl,                 & ! (inout)
    & xyzf_DQMixDtNl,                 & ! (inout)
    & xyz_DKmDtNl                     & ! (inout)
    & )
    ! 
    ! 移流計算 (3D 版)
    !
    !   移流: 4 次中央差分
    !   数値拡散 (4 階): 2 次中央差分
    !
    ! リープフロッグで, 移流を中央差分で計算するために, 
    ! 数値拡散項を追加している. 
    !

    !モジュール読み込み
    use dc_types,    only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x 方向の配列サイズ
      &                     jmin, jmax,       &! y 方向の配列サイズ
      &                     kmin, kmax,       &! z 方向の配列サイズ
      &                     nx, ny, nz,       &! 物理領域のサイズ
      &                     ncmax              ! 物質数
    use timeset,     only : TimeN
    use axesset,     only : dx, dy, dz         ! 格子間隔
    use composition, only : SpcWetSymbol

    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(inout) :: xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)             :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: pyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xqz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xqz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_nDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    integer              :: f


    !--------------------------------------------------------------------
    ! 乱流拡散係数

    ! 移流および数値拡散
    !    
    call AdvC4_nDiff_xyz( xyz_KmBl, xyz_KmNl, xyr_Dummy ) !(IN)
    
    ! tendency の更新
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_nDiff + xyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! 温位

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyz( xyz_PTempBl, xyz_PTempNl, xyr_DPTempBZDz ) !(IN)
    
    ! tendency の更新
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_nDiff + xyz_Adv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',   xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff',  xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! エクスナー関数

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyz( xyz_ExnerBl, xyz_ExnerNl, xyr_DExnerBZDz ) !(IN)
    
    ! tendency の更新
    !
    xyz_nDiff = 0.0d0
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_nDiff + xyz_Adv ) * FactorDExnerDtAdv
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))
    
    !--------------------------------------------------------------------
    ! 混合比
    ! 

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyzf( xyzf_QMixBl, xyzf_QMixNl, xyrf_DQMixBZDz ) !(IN)

    ! 落下項
    !
    call QMixFall

    ! tendency の更新    
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_nDiff + xyzf_Adv + xyzf_Fall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_nDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', &
        & xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------------------------------------------
    ! VelX, VelY, VelZ
    ! 

    ! 移流項・数値拡散項をまとめて計算
    !
    call AdvC4_nDiff_pyz_xqz_xyr

    ! tendency of VelX
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )

    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_nDiff(1:nx,1:ny,1:nz))

    ! tendency of VelY
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_NDiff + xqz_Adv )

    call HistoryAutoPut(TimeN, 'DVelYDtAdv',  xqz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_nDiff(1:nx,1:ny,1:nz))

    ! Buoyancy 
    ! 
    call BuoyancyLong_xyr

    ! tendency of VelZ
    !
    xyr_DVelZDtNl = xyr_DVelZDtNl                        &
      &             + (                                  &
      &                 xyr_NDiff                        &
      &               + xyr_Adv                          &
      &               + (                                &
      &                 + xyr_BuoyM * FactorBuoyMolWt    &
      &                 + xyr_BuoyD * FactorBuoyLoading  &
      &                 + xyr_BuoyT * FactorBuoyTemp     &
      &                 )                                &
      &               )

    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_nDiff(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyT', xyr_BuoyT(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyM', xyr_BuoyM(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyD', xyr_BuoyD(1:nx,1:ny,1:nz))    

  contains
    
    subroutine QmixFall
      !
      ! 雨粒の落下による移流を求める. 
      ! 

      ! モジュール呼び出し
      !
      use average,     only : xyr_xyz
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyzf_QMixBZ, &!基本場の混合比
        &                     xyz_DensBZ
      use differentiate_center4, &
        &              only : xyz_dz_xyr
      
      !暗黙の型宣言禁止
      !
      implicit none
      
      !変数定義
      !
      real(DP)  :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !蒸気混合比(擾乱 + 平均場)
      real(DP)  :: xyrf_QMixFallFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !雨粒の落下効果
      real(DP)  :: xyz_VelZFall(imin:imax,jmin:jmax,kmin:kmax)
                                                 !雨粒落下速度
      integer   :: s, iR

      ! 初期化
      !
      xyzf_QMixAll = max( 0.0d0, xyzf_QMixBl + xyzf_QMixBZ )
      xyrf_QMixFallFlux = 0.0d0
      xyzf_Fall = 0.0d0
      xyz_VelZFall = 0.0d0

      ! 落下項の計算
      !
      do s = 1, RainNum

        iR = IdxR(s)

        !雨粒終端速度
        xyz_VelZFall = - 12.2d0 * FactorJ * ( xyzf_QMixAll(:,:,:,iR) ** 0.125d0 )
        
        ! フラックスの計算
        !
        xyrf_QMixFallFlux(:,:,:,iR) =                                  &
          &  xyr_xyz (                                             &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall      &
          &  )

        ! 上端のフラックスはゼロ
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! 雨粒落下による時間変化
        !        
        xyzf_Fall(:,:,:,iR) =                                          &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall


    subroutine AdvC4_nDiff_pyz_xqz_xyr

      implicit none

      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            pyz_Adv(i,j,k) =                                                   &
              & - pyz_VelXNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   pyz_VelXNl(i+1,j,k) - pyz_VelXNl(i-1,j,k) )   &
              &     - fct2 * (   pyz_VelXNl(i+2,j,k) + pyz_VelXNl(i+1,j,k)     &
              &                - pyz_VelXNl(i-1,j,k) - pyz_VelXNl(i-2,j,k) )   &
              &     ) * 5.0d-1 / dx                                            &
              & - (                                                            &
              &   + ( xqz_VelYNl(i+1,j,k) + xqz_VelYNl(i,j,k) )                &  
              &     * (                                                        &  
              &         fct1 * ( pyz_VelXNl(i,j+1,k) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j+2,k) - pyz_VelXNl(i,j-1,k) )   &
              &       )                                                        &
              &   + ( xqz_VelYNl(i+1,j-1,k) + xqz_VelYNl(i,j-1,k) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j-1,k) )   &
              &       - fct2 * ( pyz_VelXNl(i,j+1,k) - pyz_VelXNl(i,j-2,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dy                                              &
              & - (                                                            &
              &   + ( xyr_VelZNl(i+1,j,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+2) - pyz_VelXNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i+1,j,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j,k-1) )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
            
          end do
        end do
      end do
      
      pyz_Adv(imin:imin+1,:,:) = 0.0d0
      pyz_Adv(imax-1:imax,:,:) = 0.0d0
      pyz_Adv(:,jmin:jmin+1,:) = 0.0d0
      pyz_Adv(:,jmax-1:jmax,:) = 0.0d0
      pyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      pyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xqz_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j+1,k) + pyz_VelXNl(i,j,k) )                &
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i+1,j,k) - xqz_VelYNl(i,j,k)   )   &
              &       - fct2 * ( xqz_VelYNl(i+2,j,k) - xqz_VelYNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j+1,k) + pyz_VelXNl(i-1,j,k) )            &
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k)   - xqz_VelYNl(i-1,j,k) )   &
              &       - fct2 * ( xqz_VelYNl(i+1,j,k) - xqz_VelYNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - xqz_VelYNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xqz_VelYNl(i,j+1,k) - xqz_VelYNl(i,j-1,k) )   &
              &     - fct2 * (   xqz_VelYNl(i,j+2,k) + xqz_VelYNl(i,j+1,k)     &
              &                - xqz_VelYNl(i,j-1,k) - xqz_VelYNl(i,j-2,k) )   &
              &     ) * 5.0d-1 / dy                                            &
              & - (                                                            &
              &   + ( xyr_VelZNl(i,j+1,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k+1) - xqz_VelYNl(i,j,k)   )   &
              &       - fct2 * ( xqz_VelYNl(i,j,k+2) - xqz_VelYNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i,j+1,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( xqz_VelYNl(i,j,k)   - xqz_VelYNl(i,j,k-1) )   &
              &       - fct2 * ( xqz_VelYNl(i,j,k+1) - xqz_VelYNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
          end do
        end do
      end do
      
      xqz_Adv(imin:imin+1,:,:) = 0.0d0
      xqz_Adv(imax-1:imax,:,:) = 0.0d0
      xqz_Adv(:,jmin:jmin+1,:) = 0.0d0
      xqz_Adv(:,jmax-1:jmax,:) = 0.0d0
      xqz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xqz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyr_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j,k+1) + pyz_VelXNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i+2,j,k) - xyr_VelZNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j,k+1) + pyz_VelXNl(i-1,j,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i-1,j,k) )   &
              &       - fct2 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - (                                                            &
              &   + ( xqz_VelYNl(i,j,k+1) + xqz_VelYNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j+1,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i,j+2,k) - xyr_VelZNl(i,j-1,k) )   &
              &       )                                                        &
              &   + ( xqz_VelYNl(i,j-1,k+1) + xqz_VelYNl(i,j-1,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i,j-1,k) )   &
              &       - fct2 * ( xyr_VelZNl(i,j+1,k) - xyr_VelZNl(i,j-2,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dy                                              &
              & - xyr_VelZNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xyr_VelZNl(i,j,k+1) - xyr_VelZNl(i,j,k-1) )   &
              &     - fct2 * (   xyr_VelZNl(i,j,k+2) + xyr_VelZNl(i,j,k+1)     &
              &                - xyr_VelZNl(i,j,k-1) - xyr_VelZNl(i,j,k-2) )   &
              &     ) * 5.0d-1 / dz
          end do
        end do
      end do
      
      xyr_Adv(imin:imin+1,:,:) = 0.0d0
      xyr_Adv(imax-1:imax,:,:) = 0.0d0
      xyr_Adv(:,jmin:jmin+1,:) = 0.0d0
      xyr_Adv(:,jmax-1:jmax,:) = 0.0d0
      xyr_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyr_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            pyz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + pyz_VelXBl(i+2,j,k)                  &
              &     + pyz_VelXBl(i-2,j,k)                  &
              &     - pyz_VelXBl(i+1,j,k) * 4.0d0          &
              &     - pyz_VelXBl(i-1,j,k) * 4.0d0          &
              &     + pyz_VelXBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + pyz_VelXBl(i,j+2,k)                  &
              &     + pyz_VelXBl(i,j-2,k)                  &
              &     - pyz_VelXBl(i,j+1,k) * 4.0d0          &
              &     - pyz_VelXBl(i,j-1,k) * 4.0d0          &
              &     + pyz_VelXBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        & 
              &     + pyz_VelXBl(i,j,k+2)                  &
              &     + pyz_VelXBl(i,j,k-2)                  &
              &     - pyz_VelXBl(i,j,k+1) * 4.0d0          &
              &     - pyz_VelXBl(i,j,k-1) * 4.0d0          &
              &     + pyz_VelXBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
            
          end do
        end do
      end do
      
      pyz_nDiff(imin:imin+1,:,:) = 0.0d0
      pyz_nDiff(imax-1:imax,:,:) = 0.0d0
      pyz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      pyz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      pyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      pyz_nDiff(:,:,kmax-1:kmax) = 0.0d0

      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xqz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xqz_VelYBl(i+2,j,k)                  &
              &     + xqz_VelYBl(i-2,j,k)                  &
              &     - xqz_VelYBl(i+1,j,k) * 4.0d0          &
              &     - xqz_VelYBl(i-1,j,k) * 4.0d0          &
              &     + xqz_VelYBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xqz_VelYBl(i,j+2,k)                  &
              &     + xqz_VelYBl(i,j-2,k)                  &
              &     - xqz_VelYBl(i,j+1,k) * 4.0d0          &
              &     - xqz_VelYBl(i,j-1,k) * 4.0d0          &
              &     + xqz_VelYBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        & 
              &     + xqz_VelYBl(i,j,k+2)                  &
              &     + xqz_VelYBl(i,j,k-2)                  &
              &     - xqz_VelYBl(i,j,k+1) * 4.0d0          &
              &     - xqz_VelYBl(i,j,k-1) * 4.0d0          &
              &     + xqz_VelYBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xqz_nDiff(imin:imin+1,:,:) = 0.0d0
      xqz_nDiff(imax-1:imax,:,:) = 0.0d0
      xqz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xqz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xqz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xqz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyr_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xyr_VelZBl(i+2,j,k)                  &
              &     + xyr_VelZBl(i-2,j,k)                  &
              &     - xyr_VelZBl(i+1,j,k) * 4.0d0          &
              &     - xyr_VelZBl(i-1,j,k) * 4.0d0          &
              &     + xyr_VelZBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j+2,k)                  &
              &     + xyr_VelZBl(i,j-2,k)                  &
              &     - xyr_VelZBl(i,j+1,k) * 4.0d0          &
              &     - xyr_VelZBl(i,j-1,k) * 4.0d0          &
              &     + xyr_VelZBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHm / ( dy ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j,k+2)                  &
              &     + xyr_VelZBl(i,j,k-2)                  &
              &     - xyr_VelZBl(i,j,k+1) * 4.0d0          &
              &     - xyr_VelZBl(i,j,k-1) * 4.0d0          &
              &     + xyr_VelZBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyr_nDiff(imin:imin+1,:,:) = 0.0d0
      xyr_nDiff(imax-1:imax,:,:) = 0.0d0
      xyr_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xyr_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xyr_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyr_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_pyz_xqz_xyr
    
    
    subroutine AdvC4_nDiff_xyz( xyz_VarBl, xyz_VarNl, xyr_DVarBZDz )
      
      implicit none
      
      real(DP), intent(in)  :: xyz_VarBl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyz_VarNl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyr_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax)
      
      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0

      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyz_Adv(i,j,k) =                                                  &
              & - (                                                           &
              &      pyz_VelXNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i+2,j,k) - xyz_VarNl(i-1,j,k) ) &
              &          )                                                    &
              &    + pyz_VelXNl(i-1,j,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i-1,j,k) ) &
              &          - fct2 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i-2,j,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dx                                             &
              & - (                                                           &
              &      xqz_VelYNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j+1,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j+2,k) - xyz_VarNl(i,j-1,k) ) &
              &          )                                                    &
              &    + xqz_VelYNl(i,j-1,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j-1,k) ) &
              &          - fct2 * ( xyz_VarNl(i,j+1,k) - xyz_VarNl(i,j-2,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dy                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+2) - xyz_VarNl(i,j,k-1) ) &
              &          )                                                    &
              &    + xyr_VelZNl(i,j,k-1)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j,k-1) ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k-2) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dz                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)   * xyr_DVarBZDz(i,j,k)                & 
              &    + xyr_VelZNl(i,j,k-1) * xyr_DVarBZDz(i,j,k-1)              &
              &   ) * 5.0d-1 
          end do
        end do
      end do
      
      xyz_Adv(imin:imin+1,:,:) = 0.0d0
      xyz_Adv(imax-1:imax,:,:) = 0.0d0
      xyz_Adv(:,jmin:jmin+1,:) = 0.0d0
      xyz_Adv(:,jmax-1:jmax,:) = 0.0d0
      xyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      ! 4 次の数値拡散: 2 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = jmin + 2, jmax - 2
          do i = imin + 2, imax - 2
            
            xyz_nDiff(i,j,k) =                            &
              & - (                                       &
              &     + xyz_VarBl(i+2,j,k)                  &
              &     + xyz_VarBl(i-2,j,k)                  &
              &     - xyz_VarBl(i+1,j,k) * 4.0d0          &
              &     - xyz_VarBl(i-1,j,k) * 4.0d0          &
              &     + xyz_VarBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j+2,k)                  &
              &     + xyz_VarBl(i,j-2,k)                  &
              &     - xyz_VarBl(i,j+1,k) * 4.0d0          &
              &     - xyz_VarBl(i,j-1,k) * 4.0d0          &
              &     + xyz_VarBl(i,j  ,k) * 6.0d0          &
              &   ) * NuHh / ( dy ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j,k+2)                  &
              &     + xyz_VarBl(i,j,k-2)                  &
              &     - xyz_VarBl(i,j,k+1) * 4.0d0          &
              &     - xyz_VarBl(i,j,k-1) * 4.0d0          &
              &     + xyz_VarBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyz_nDiff(imin:imin+1,:,:) = 0.0d0
      xyz_nDiff(imax-1:imax,:,:) = 0.0d0
      xyz_nDiff(:,jmin:jmin+1,:) = 0.0d0
      xyz_nDiff(:,jmax-1:jmax,:) = 0.0d0
      xyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyz
    
    
    subroutine AdvC4_nDiff_xyzf( xyzf_VarBl, xyzf_VarNl, xyrf_DVarBZDz ) 

      implicit none
      
      real(DP), intent(in)  :: xyzf_VarBl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyzf_VarNl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyrf_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      
      real(DP)              :: fct1, fct2
      integer               :: s, i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = jmin + 2, jmax - 2
            do i = imin + 2, imax - 2
              
              xyzf_Adv(i,j,k,s) =                                                     &
                & - (                                                                 &
                &      pyz_VelXNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i+2,j,k,s) - xyzf_VarNl(i-1,j,k,s) ) &
                &          )                                                          &
                &    + pyz_VelXNl(i-1,j,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i-1,j,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i-2,j,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dx                                                   &
                & - (                                                                 &
                &      xqz_VelYNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j+1,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j+2,k,s) - xyzf_VarNl(i,j-1,k,s) ) &
                &          )                                                          &
                &    + xqz_VelYNl(i,j-1,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j-1,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j+1,k,s) - xyzf_VarNl(i,j-2,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dy                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+2,s) - xyzf_VarNl(i,j,k-1,s) ) &
                &          )                                                          &
                &    + xyr_VelZNl(i,j,k-1)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j,k-1,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k-2,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dz                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)   * xyrf_DVarBZDz(i,j,k,s)                   &
                &    + xyr_VelZNl(i,j,k-1) * xyrf_DVarBZDz(i,j,k-1,s)                 & 
                &   ) * 5.0d-1 
            end do
          end do
        end do
      end do
      
      xyzf_Adv(imin:imin+1,:,:,:) = 0.0d0
      xyzf_Adv(imax-1:imax,:,:,:) = 0.0d0
      xyzf_Adv(:,jmin:jmin+1,:,:) = 0.0d0
      xyzf_Adv(:,jmax-1:jmax,:,:) = 0.0d0
      xyzf_Adv(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_Adv(:,:,kmax-1:kmax,:) = 0.0d0
      
      ! 数値拡散: 2 次精度中心差分
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = jmin + 2, jmax - 2
            do i = imin + 2, imax - 2
              
              xyzf_nDiff(i,j,k,s) =                            &
                & - (                                          &
                &       xyzf_VarBl(i+2,j,k,s)                  &
                &     + xyzf_VarBl(i-2,j,k,s)                  &
                &     - xyzf_VarBl(i+1,j,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i-1,j,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i  ,j,k,s) * 6.0d0          &
                &   ) * NuHh / ( dx ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j+2,k,s)                  &
                &     + xyzf_VarBl(i,j-2,k,s)                  &
                &     - xyzf_VarBl(i,j+1,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j-1,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j  ,k,s) * 6.0d0          &
                &   ) * NuHh / ( dy ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j,k+2,s)                  &
                &     + xyzf_VarBl(i,j,k-2,s)                  &
                &     - xyzf_VarBl(i,j,k+1,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j,k-1,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j,k  ,s) * 6.0d0          &
                &   ) * NuVh / ( dz ** 4.0d0 )
            end do
          end do
        end do
      end do
      
      xyzf_nDiff(imin:imin+1,:,:,:) = 0.0d0
      xyzf_nDiff(imax-1:imax,:,:,:) = 0.0d0
      xyzf_nDiff(:,jmin:jmin+1,:,:) = 0.0d0
      xyzf_nDiff(:,jmax-1:jmax,:,:) = 0.0d0
      xyzf_nDiff(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_nDiff(:,:,kmax-1:kmax,:) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyzf
  

    subroutine BuoyancyLong_xyr
      
      use composition, only: GasNum,       &! 
        &                    IdxG,         &!
        &                    MolWtWet       ! 湿潤成分の分子量
      use constants,   only: MolWtDry,     &! 乾燥成分の分子量
        &                    Grav           ! 重力加速度
      use basicset,    only: xyr_QMixBZ,         &
        &                    xyr_QMixBZPerMolWt, &
        &                    xyz_PTempBZ
      
      implicit none
      
      real(DP)              :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)              :: tmp1(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp2(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp3(imin:imax,jmin:jmax,kmin:kmax)
      integer               :: i, j, k, f, n
      
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do
      
      ! Buoyancy due to temperature disturbunce
      !
      do k = kmin, kmax - 1
        do j = jmin, jmax
          do i = imin, imax
            
            xyr_BuoyT(i,j,k) =                                  &
              & Grav                                            &
              & * (                                             &
              &     xyz_PTempNl(i,j,k+1) / xyz_PTempBZ(i,j,k+1) &
              &   + xyz_PTempNl(i,j,k)   / xyz_PTempBZ(i,j,k)   &
              &   ) * 5.0d-1
            
          end do
        end do
      end do
      
      xyr_BuoyT(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to molecular weight
      !
      tmp1 = sum(xyzf_QMixPerMolWt, 4) 
      tmp2 = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax - 1
        do j = jmin, jmax 
          do i = imin, imax 
            
            xyr_BuoyM(i,j,k) =                                       &
              & + Grav                                               &
              &   * ( tmp1(i,j,k+1) + tmp1(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) ) &
              & - Grav                                               &
              &   * ( tmp2(i,j,k+1) + tmp2(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 + xyr_QmixBZ(i,j,k) ) 
            
          end do
        end do
      end do
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      tmp3 = sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) 
      
      do k = kmin, kmax - 1
        do j = jmin, jmax 
          do i = imin, imax 
            
            xyr_BuoyD(i,j,k) =                                &
              & - Grav                                        &
              &   * ( tmp3(i,j,k+1) + tmp3(i,j,k) ) * 5.0d-1  &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do
      
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr
    
  end subroutine Dynamics3D_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics2D_Long_forcing( &
    & pyz_VelXBl,  pyz_VelXNl,        & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,        & ! (in)
    & xyz_PTempBl, xyz_PTempNl,       & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,       & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,       & ! (in)
    & xyz_KmBl,    xyz_KmNl,          & ! (in)
    & pyz_DVelXDtNl,                  & ! (inout)
    & xqz_DVelYDtNl,                  & ! (inout)
    & xyr_DVelZDtNl,                  & ! (inout)
    & xyz_DPTempDtNl,                 & ! (inout)
    & xyz_DExnerDtNl,                 & ! (inout)
    & xyzf_DQMixDtNl,                 & ! (inout)
    & xyz_DKmDtNl                     & ! (inout)
    & )
    ! 
    ! 移流計算 (2D 版)
    !
    !   移流: 4 次中央差分
    !   数値拡散 (4 階): 2 次中央差分
    !
    ! リープフロッグで, 移流を中央差分で計算するために, 
    ! 数値拡散項を追加している. 
    !

    !モジュール読み込み
    use dc_types,    only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x 方向の配列サイズ
      &                     jmin, jmax,       &! y 方向の配列サイズ
      &                     kmin, kmax,       &! z 方向の配列サイズ
      &                     nx, ny, nz,       &! 物理領域のサイズ
      &                     ncmax              ! 物質数
    use timeset,     only : TimeN
    use axesset,     only : dx, dz             ! 格子間隔
    use composition, only : SpcWetSymbol

    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(in)    :: xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP), intent(inout) :: xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP), intent(inout) :: xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)             :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: pyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_nDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)             :: xyzf_nDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    integer              :: f


    !--------------------------------------------------------------------
    ! 乱流拡散係数

    ! 移流および数値拡散
    !
    call AdvC4_nDiff_xyz( xyz_KmBl, xyz_KmNl, xyr_Dummy ) !(IN)
    
    ! tendency の更新
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_nDiff + xyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DKmDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))
    
    !--------------------------------------------------------------------
    ! 温位

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyz( xyz_PTempBl, xyz_PTempNl, xyr_DPTempBZDz ) !(IN)
 
    ! tendency の更新
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_nDiff + xyz_Adv )
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',   xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff',  xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! エクスナー関数

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyz( xyz_ExnerBl, xyz_ExnerNl, xyr_DExnerBZDz ) !(IN)

    ! tendency の更新
    !
    xyz_nDiff = 0.0d0
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_nDiff + xyz_Adv ) * FactorDExnerDtAdv
    
    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_nDiff(1:nx,1:ny,1:nz))

    !--------------------------------------------------------------------
    ! 混合比
    ! 

    ! 移流については, 基本場も考慮する.
    !
    call AdvC4_nDiff_xyzf( xyzf_QMixBl, xyzf_QMixNl, xyrf_DQMixBZDz ) !(IN)

    ! 落下項
    !
    call QMixFall

    ! tendency の更新    
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_nDiff + xyzf_Adv + xyzf_Fall )
    
    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',   &
        & xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', &
        & xyzf_nDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', &
        & xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------------------------------------------
    ! VelX, VelY, VelZ
    ! 

    ! 移流項・数値拡散項をまとめて計算
    !
    call AdvC4_nDiff_pyz_xyr

    ! tendency of VelX
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )
    
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',  pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_nDiff(1:nx,1:ny,1:nz))

    ! tendency of VelY
    !
    xqz_DVelYDtNl = 0.0d0

    ! Buoyancy 
    ! 
    call BuoyancyLong_xyr

    ! tendency of VelZ
    !
    xyr_DVelZDtNl = xyr_DVelZDtNl                        &
      &             + (                                  &
      &                 xyr_NDiff                        &
      &               + xyr_Adv                          &
      &               + (                                &
      &                 + xyr_BuoyM * FactorBuoyMolWt    &
      &                 + xyr_BuoyD * FactorBuoyLoading  &
      &                 + xyr_BuoyT * FactorBuoyTemp     &
      &                 )                                &
      &               )

    call HistoryAutoPut(TimeN, 'DVelZDtAdv',   xyr_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtDiff',  xyr_NDiff(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyT', xyr_BuoyT(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyM', xyr_BuoyM(1:nx,1:ny,1:nz))    
    call HistoryAutoPut(TimeN, 'DVelZDtBuoyD', xyr_BuoyD(1:nx,1:ny,1:nz))    

  contains

    subroutine QmixFall
      !
      ! 雨粒の落下による移流を求める. 
      ! 

      ! モジュール呼び出し
      !
      use average,     only : xyr_xyz
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyzf_QMixBZ, &!基本場の混合比
        &                     xyz_DensBZ
      use differentiate_center4, &
        &              only : xyz_dz_xyr
      
      !暗黙の型宣言禁止
      !
      implicit none
      
      !変数定義
      !
      real(DP)  :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !蒸気混合比(擾乱 + 平均場)
      real(DP)  :: xyrf_QMixFallFlux(imin:imax,jmin:jmax,kmin:kmax,ncmax)
                                                 !雨粒の落下効果
      real(DP)  :: xyz_VelZFall(imin:imax,jmin:jmax,kmin:kmax)
                                                 !雨粒落下速度
      integer   :: s, iR

      ! 初期化
      !
      xyzf_QMixAll = max( 0.0d0, xyzf_QMixBl + xyzf_QMixBZ )
      xyrf_QMixFallFlux = 0.0d0
      xyzf_Fall = 0.0d0
      xyz_VelZFall = 0.0d0

      ! 落下項の計算
      !
      do s = 1, RainNum

        iR = IdxR(s)

        !雨粒終端速度
        xyz_VelZFall = - 12.2d0 * FactorJ * ( xyzf_QMixAll(:,:,:,iR) ** 0.125d0 )
        
        ! フラックスの計算
        !
        xyrf_QMixFallFlux(:,:,:,iR) =                                  &
          &  xyr_xyz (                                                 &
          &    xyz_DensBZ * xyzf_QMixAll(:,:,:,iR) * xyz_VelZFall      &
          &  )

        ! 上端のフラックスはゼロ
        !        
        xyrf_QMixFallFlux(:,:,nz,iR) = 0.0d0

        ! 雨粒落下による時間変化 (DelTime をかけてある)
        !        
        xyzf_Fall(:,:,:,iR) =                                          &
          &  - xyz_dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall

    
    subroutine AdvC4_nDiff_pyz_xyr

      implicit none

      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            pyz_Adv(i,j,k) =                                                   &
              & - pyz_VelXNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   pyz_VelXNl(i+1,j,k) - pyz_VelXNl(i-1,j,k) )   &
              &     - fct2 * (   pyz_VelXNl(i+2,j,k) + pyz_VelXNl(i+1,j,k)     &
              &                - pyz_VelXNl(i-1,j,k) - pyz_VelXNl(i-2,j,k) )   &
              &     ) * 5.0d-1 / dx                                            &
              & - (                                                            &
              &   + ( xyr_VelZNl(i+1,j,k) + xyr_VelZNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k)   )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+2) - pyz_VelXNl(i,j,k-1) )   &
              &       )                                                        &
              &   + ( xyr_VelZNl(i+1,j,k-1) + xyr_VelZNl(i,j,k-1) )            & 
              &     * (                                                        &
              &         fct1 * ( pyz_VelXNl(i,j,k)   - pyz_VelXNl(i,j,k-1) )   &
              &       - fct2 * ( pyz_VelXNl(i,j,k+1) - pyz_VelXNl(i,j,k-2) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dz
            
          end do
        end do
      end do
      
      pyz_Adv(imin:imin+1,:,:) = 0.0d0
      pyz_Adv(imax-1:imax,:,:) = 0.0d0
      pyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      pyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyr_Adv(i,j,k) =                                                   &
              & - (                                                            &
              &   + ( pyz_VelXNl(i,j,k+1) + pyz_VelXNl(i,j,k) )                & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i,j,k)   )   &
              &       - fct2 * ( xyr_VelZNl(i+2,j,k) - xyr_VelZNl(i-1,j,k) )   &
              &       )                                                        &
              &   + ( pyz_VelXNl(i-1,j,k+1) + pyz_VelXNl(i-1,j,k) )            & 
              &     * (                                                        &
              &         fct1 * ( xyr_VelZNl(i,j,k)   - xyr_VelZNl(i-1,j,k) )   &
              &       - fct2 * ( xyr_VelZNl(i+1,j,k) - xyr_VelZNl(i-2,j,k) )   &
              &       )                                                        &
              &   ) * 2.5d-1 / dx                                              &
              & - xyr_VelZNl(i,j,k)                                            &
              &   * (                                                          &
              &       fct1 * (   xyr_VelZNl(i,j,k+1) - xyr_VelZNl(i,j,k-1) )   &
              &     - fct2 * (   xyr_VelZNl(i,j,k+2) + xyr_VelZNl(i,j,k+1)     &
              &                - xyr_VelZNl(i,j,k-1) - xyr_VelZNl(i,j,k-2) )   &
              &     ) * 5.0d-1 / dz
          end do
        end do
      end do
      
      xyr_Adv(imin:imin+1,:,:) = 0.0d0
      xyr_Adv(imax-1:imax,:,:) = 0.0d0
      xyr_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyr_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            pyz_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + pyz_VelXBl(i+2,j,k)                  &
              &     + pyz_VelXBl(i-2,j,k)                  &
              &     - pyz_VelXBl(i+1,j,k) * 4.0d0          &
              &     - pyz_VelXBl(i-1,j,k) * 4.0d0          &
              &     + pyz_VelXBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        & 
              &     + pyz_VelXBl(i,j,k+2)                  &
              &     + pyz_VelXBl(i,j,k-2)                  &
              &     - pyz_VelXBl(i,j,k+1) * 4.0d0          &
              &     - pyz_VelXBl(i,j,k-1) * 4.0d0          &
              &     + pyz_VelXBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
            
          end do
        end do
      end do
      
      pyz_nDiff(imin:imin+1,:,:) = 0.0d0
      pyz_nDiff(imax-1:imax,:,:) = 0.0d0
      pyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      pyz_nDiff(:,:,kmax-1:kmax) = 0.0d0

      do k = kmin + 2, kmax - 2
        do j = 1, ny          
          do i = imin + 2, imax - 2
            
            xyr_nDiff(i,j,k) =                             &
              & - (                                        &
              &     + xyr_VelZBl(i+2,j,k)                  &
              &     + xyr_VelZBl(i-2,j,k)                  &
              &     - xyr_VelZBl(i+1,j,k) * 4.0d0          &
              &     - xyr_VelZBl(i-1,j,k) * 4.0d0          &
              &     + xyr_VelZBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHm / ( dx ** 4.0d0 )               &
              & - (                                        &
              &     + xyr_VelZBl(i,j,k+2)                  &
              &     + xyr_VelZBl(i,j,k-2)                  &
              &     - xyr_VelZBl(i,j,k+1) * 4.0d0          &
              &     - xyr_VelZBl(i,j,k-1) * 4.0d0          &
              &     + xyr_VelZBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVm / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyr_nDiff(imin:imin+1,:,:) = 0.0d0
      xyr_nDiff(imax-1:imax,:,:) = 0.0d0
      xyr_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyr_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_pyz_xyr
    
    
    subroutine AdvC4_nDiff_xyz( xyz_VarBl, xyz_VarNl, xyr_DVarBZDz )
      
      implicit none
      
      real(DP), intent(in)  :: xyz_VarBl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyz_VarNl(imin:imax,jmin:jmax,kmin:kmax)
      real(DP), intent(in)  :: xyr_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax)
      
      real(DP)              :: fct1, fct2
      integer               :: i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0

      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyz_Adv(i,j,k) =                                                  &
              & - (                                                           &
              &      pyz_VelXNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i+2,j,k) - xyz_VarNl(i-1,j,k) ) &
              &          )                                                    &
              &    + pyz_VelXNl(i-1,j,k)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i-1,j,k) ) &
              &          - fct2 * ( xyz_VarNl(i+1,j,k) - xyz_VarNl(i-2,j,k) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dx                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)                                        &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k)   ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+2) - xyz_VarNl(i,j,k-1) ) &
              &          )                                                    &
              &    + xyr_VelZNl(i,j,k-1)                                      &
              &        * (                                                    &
              &            fct1 * ( xyz_VarNl(i,j,k)   - xyz_VarNl(i,j,k-1) ) &
              &          - fct2 * ( xyz_VarNl(i,j,k+1) - xyz_VarNl(i,j,k-2) ) &
              &          )                                                    &
              &   ) * 5.0d-1 / dz                                             &
              & - (                                                           &
              &      xyr_VelZNl(i,j,k)   * xyr_DVarBZDz(i,j,k)                &
              &    + xyr_VelZNl(i,j,k-1) * xyr_DVarBZDz(i,j,k-1)              &
              &   ) * 5.0d-1
          end do
        end do
      end do
      
      xyz_Adv(imin:imin+1,:,:) = 0.0d0
      xyz_Adv(imax-1:imax,:,:) = 0.0d0
      xyz_Adv(:,:,kmin:kmin+1) = 0.0d0
      xyz_Adv(:,:,kmax-1:kmax) = 0.0d0
      
      ! 4 次の数値拡散: 2 次精度中心差分
      ! 
      do k = kmin + 2, kmax - 2
        do j = 1, ny
          do i = imin + 2, imax - 2
            
            xyz_nDiff(i,j,k) =                            &
              & - (                                       &
              &     + xyz_VarBl(i+2,j,k)                  &
              &     + xyz_VarBl(i-2,j,k)                  &
              &     - xyz_VarBl(i+1,j,k) * 4.0d0          &
              &     - xyz_VarBl(i-1,j,k) * 4.0d0          &
              &     + xyz_VarBl(i  ,j,k) * 6.0d0          &
              &   ) * NuHh / ( dx ** 4.0d0 )              &
              & - (                                       &
              &       xyz_VarBl(i,j,k+2)                  &
              &     + xyz_VarBl(i,j,k-2)                  &
              &     - xyz_VarBl(i,j,k+1) * 4.0d0          &
              &     - xyz_VarBl(i,j,k-1) * 4.0d0          &
              &     + xyz_VarBl(i,j,k  ) * 6.0d0          &
              &   ) * NuVh / ( dz ** 4.0d0 )
          end do
        end do
      end do
      
      xyz_nDiff(imin:imin+1,:,:) = 0.0d0
      xyz_nDiff(imax-1:imax,:,:) = 0.0d0
      xyz_nDiff(:,:,kmin:kmin+1) = 0.0d0
      xyz_nDiff(:,:,kmax-1:kmax) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyz
    
    
    subroutine AdvC4_nDiff_xyzf( xyzf_VarBl, xyzf_VarNl, xyrf_DVarBZDz ) 

      implicit none
      
      real(DP), intent(in)  :: xyzf_VarBl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyzf_VarNl(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      real(DP), intent(in)  :: xyrf_DVarBZDz(imin:imax,jmin:jmax,kmin:kmax,1:ncmax)
      
      real(DP)              :: fct1, fct2
      integer               :: s, i, j, k
      
      ! 微分に用いる係数を予め計算
      !
      fct1 = 9.0d0 / 8.0d0
      fct2 = 1.0d0 / 24.0d0
      
      ! 移流項の計算. 移流項: 4 次精度中心差分
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = 1, ny
            do i = imin + 2, imax - 2
              
              xyzf_Adv(i,j,k,s) =                                                     &
                & - (                                                                 &
                &      pyz_VelXNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i+2,j,k,s) - xyzf_VarNl(i-1,j,k,s) ) &
                &          )                                                          &
                &    + pyz_VelXNl(i-1,j,k)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i-1,j,k,s) ) &
                &          - fct2 * ( xyzf_VarNl(i+1,j,k,s) - xyzf_VarNl(i-2,j,k,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dx                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)                                              &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k,s)   ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+2,s) - xyzf_VarNl(i,j,k-1,s) ) &
                &          )                                                          &
                &    + xyr_VelZNl(i,j,k-1)                                            &
                &        * (                                                          &
                &            fct1 * ( xyzf_VarNl(i,j,k,s)   - xyzf_VarNl(i,j,k-1,s) ) &
                &          - fct2 * ( xyzf_VarNl(i,j,k+1,s) - xyzf_VarNl(i,j,k-2,s) ) &
                &          )                                                          &
                &   ) * 5.0d-1 / dz                                                   &
                & - (                                                                 &
                &      xyr_VelZNl(i,j,k)   * xyrf_DVarBZDz(i,j,k,s)                   &
                &    + xyr_VelZNl(i,j,k-1) * xyrf_DVarBZDz(i,j,k-1,s)                 &
                &   ) * 5.0d-1 
            end do
          end do
        end do
      end do
      
      xyzf_Adv(imin:imin+1,:,:,:) = 0.0d0
      xyzf_Adv(imax-1:imax,:,:,:) = 0.0d0
      xyzf_Adv(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_Adv(:,:,kmax-1:kmax,:) = 0.0d0
      
      ! 数値拡散: 2 次精度中心差分
      ! 
      do s = 1, ncmax
        do k = kmin + 2, kmax - 2
          do j = 1, ny
            do i = imin + 2, imax - 2
              
              xyzf_nDiff(i,j,k,s) =                            &
                & - (                                          &
                &       xyzf_VarBl(i+2,j,k,s)                  &
                &     + xyzf_VarBl(i-2,j,k,s)                  &
                &     - xyzf_VarBl(i+1,j,k,s) * 4.0d0          &
                &     - xyzf_VarBl(i-1,j,k,s) * 4.0d0          &
                &     + xyzf_VarBl(i  ,j,k,s) * 6.0d0          &
                &   ) * NuHh / ( dx ** 4.0d0 )                 &
                & - (                                          &
                &       xyzf_VarBl(i,j,k+2,s)                  &
                &     + xyzf_VarBl(i,j,k-2,s)                  &
                &     - xyzf_VarBl(i,j,k+1,s) * 4.0d0          &
                &     - xyzf_VarBl(i,j,k-1,s) * 4.0d0          &
                &     + xyzf_VarBl(i,j,k  ,s) * 6.0d0          &
                &   ) * NuVh / ( dz ** 4.0d0 )
            end do
          end do
        end do
      end do
      
      xyzf_nDiff(imin:imin+1,:,:,:) = 0.0d0
      xyzf_nDiff(imax-1:imax,:,:,:) = 0.0d0
      xyzf_nDiff(:,:,kmin:kmin+1,:) = 0.0d0
      xyzf_nDiff(:,:,kmax-1:kmax,:) = 0.0d0
      
    end subroutine AdvC4_nDiff_xyzf
  

    subroutine BuoyancyLong_xyr
      
      use composition, only: GasNum,       &! 
        &                    IdxG,         &!
        &                    MolWtWet       ! 湿潤成分の分子量
      use constants,   only: MolWtDry,     &! 乾燥成分の分子量
        &                    Grav           ! 重力加速度
      use basicset,    only: xyr_QMixBZ,         &
        &                    xyr_QMixBZPerMolWt, &
        &                    xyz_PTempBZ
      
      implicit none
      
      real(DP)              :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
      real(DP)              :: tmp1(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp2(imin:imax,jmin:jmax,kmin:kmax)
      real(DP)              :: tmp3(imin:imax,jmin:jmax,kmin:kmax)
      integer               :: i, j, k, f, n
      
      do f = 1, GasNum
        n = IdxG(f)
        xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,n) / MolWtWet(n)
      end do
      
      ! Buoyancy due to temperature disturbunce
      !
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax
            
            xyr_BuoyT(i,j,k) =                                  &
              & Grav                                            &
              & * (                                             &
              &     xyz_PTempNl(i,j,k+1) / xyz_PTempBZ(i,j,k+1) &
              &   + xyz_PTempNl(i,j,k)   / xyz_PTempBZ(i,j,k)   &
              &   ) * 5.0d-1
            
          end do
        end do
      end do
      
      xyr_BuoyT(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to molecular weight
      !
      tmp1 = sum(xyzf_QMixPerMolWt, 4) 
      tmp2 = sum(xyzf_QMixNl(:,:,:,1:GasNum), 4)
      
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax 
            
            xyr_BuoyM(i,j,k) =                                       &
              & + Grav                                               &
              &   * ( tmp1(i,j,k+1) + tmp1(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt(i,j,k) ) &
              & - Grav                                               &
              &   * ( tmp2(i,j,k+1) + tmp2(i,j,k) ) * 5.0d-1         &
              &   / ( 1.0d0 + xyr_QmixBZ(i,j,k) ) 
            
          end do
        end do
      end do
      xyr_BuoyM(:,:,kmax) = 0.0d0
      
      ! Buoyancy due to loading
      !
      tmp3 = sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) 
      
      do k = kmin, kmax - 1
        do j = 1, ny
          do i = imin, imax 
            
            xyr_BuoyD(i,j,k) =                                &
              & - Grav                                        &
              &   * ( tmp3(i,j,k+1) + tmp3(i,j,k) ) * 5.0d-1  &
              &   / ( 1.0d0 + xyr_QMixBZ(i,j,k) )
            
          end do
        end do
      end do
      
      xyr_BuoyD(:,:,kmax) = 0.0d0
      
    end subroutine BuoyancyLong_xyr
    
  end subroutine Dynamics2D_Long_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xqz_VelYNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xqz_DVelYDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! 力学コア: 短い時間ステップで評価する項の計算.
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                           !配列サイズ
      &                    FlagCalc3D      !フラグ

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)    :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)

    if ( FlagCalc3D ) then 
      call Dynamics3D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xqz_VelYNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xqz_DVelYDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    else
      call Dynamics2D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    end if

  end subroutine Dynamics_Short_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics3D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xqz_VelYNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xqz_DVelYDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! 力学コア (短い時間ステップ) 3D 版
    !

    !モジュール読み込み
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gtool_historyauto, only: HistoryAutoPut 
    
    use gridset, only: &
      &                 imin,            &! x 方向の配列の下限
      &                 imax,            &! x 方向の配列の上限
      &                 jmin,            &! y 方向の配列の下限
      &                 jmax,            &! y 方向の配列の上限
      &                 kmin,            &! z 方向の配列の下限
      &                 kmax,            &! z 方向の配列の上限
      &                 nx,              &! x 方向の物理領域の上限
      &                 ny,              &! x 方向の物理領域の上限
      &                 nz                ! y 方向の物理領域の上限
    use constants,only: CpDry, CvDry, GasRDry ! 乾燥成分の比熱
    use timeset, only:  DelTimeShort,  TimeN
    use axesset, only:  dx, dy, dz        ! 格子間隔
    use basicset, only: xyz_VelSW,      &!基本場の音速 
      &                 xyz_VPTempBZ,      &!基本場の温位
      &                 pyz_VPTempBZ,      &!基本場の温位
      &                 xqz_VPTempBZ,      &!基本場の温位
      &                 xyr_VPTempBZ        !基本場の温位
    use setmargin,only: SetMargin_xyzf, SetMargin_xyz, &
      &                 SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

    implicit none

    real(DP), intent(in)    :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)    :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)   :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_PGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_SWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------     
    ! initialize: Divergence of velocity
    !
    call VelDivC2
    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))

    !------------------------------------------------------------
    ! VelX, VelY
    !  水平方向は陽解法で解く. 
    !
    call PGrad_HE

    ! tendency 
    !
    pyz_DVelXDtNs = pyz_PGrad + pyz_SWF
    xqz_DVelYDtNs = xqz_PGrad + xqz_SWF

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * (pyz_DVelXDtNl + pyz_DVelXDtNs)
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * (xqz_DVelYDtNl + xqz_DVelYDtNs)

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)
    
    !------------------------------------------------------------
    ! Exner function
    !  積分値を返すことに注意. 
    !

    ! 短い時間ステップで評価する圧力の式の tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency の合計
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    !エクスナー関数の計算
    ! 
    call Exner_HEVI

    !------------------------------------------------------------
    ! VelZ
    !
    call PGrad_VI

    ! tendency
    !
    xyr_DVelZDtNs = xyr_PGrad + xyr_SWF

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNl + xyr_DVelZDtNs)

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  contains

    subroutine VelDivC2
      
      implicit none
      integer              :: i, j, k
      
      do  k = kmin + 1, kmax 
        do j = jmin + 1, jmax 
          do i = imin + 1, imax 
            xyz_VelDivNs(i,j,k) =        &
              & + (                      &
              &     pyz_VelXNs(i,j,k)      &
              &   - pyz_VelXNs(i-1,j,k)    &
              &   ) / dx                 &
              & + (                      &
              &     xqz_VelYNs(i,j,k)      &
              &   - xqz_VelYNs(i,j-1,k)    &
              &   ) / dy                 &
              & + (                      &
              &     xyr_VelZNs(i,j,k)      &
              &   - xyr_VelZNs(i,j,k-1)    &
              &   ) / dz
          end do
        end do
      end do
      
      xyz_VelDivNs(imin,:,:) = 0.0d0 
      xyz_VelDivNs(:,jmin,:) = 0.0d0  
      xyz_VelDivNs(:,:,kmin) = 0.0d0 
      
    end subroutine VelDivC2


    subroutine PGrad_HE
    
      implicit none
      integer              :: i, j, k
      
      !------------------------------------------------------------------
      ! X 方向
      
      do k = kmin, kmax
        do j = jmin, jmax
          do i = imin, imax - 1

            ! 音波減衰項
            !            
            pyz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i+1,j,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dx            
            
            ! 圧力傾度力
            !
            pyz_PGrad(i,j,k) =                &
              & - CpDry * pyz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i+1,j,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) / dx                                
          end do
        end do
      end do
      
      ! 穴埋め
      !
      pyz_SWF(imax,:,:)   = 0.0d0 
      pyz_PGrad(imax,:,:) = 0.0d0 
      

      !------------------------------------------------------------------
      ! Y 方向
      
      do k = kmin, kmax
        do j = jmin, jmax - 1
          do i = imin, imax

            ! 音波減衰項
            !            
            xqz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j+1,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dy
               
            ! 圧力傾度力
            !             
            xqz_PGrad(i,j,k) =                &
              & - CpDry * xqz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i,j+1,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) /dy                     
            
          end do
        end do
      end do
      
      ! 穴埋め
      !
      xqz_SWF(:,jmax,:) = 0.0d0 
      xqz_PGrad(:,jmax,:) = 0.0d0 
      
    end subroutine PGrad_HE


    subroutine PGrad_VI
    
      implicit none
      integer               :: i, j, k

      do k = kmin, kmax - 1
        do j = jmin, jmax
          do i = imin, imax

            ! 音波減衰項
            !            
            xyr_SWF(i,j,k) =                  &
              & + AlphaV                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j,k+1)   & 
              &     - xyz_VelDivNs(i,j,k)     & 
              &     ) / dz
            
            ! 圧力傾度力
            !
            xyr_PGrad(i,j,k) =                 &
              & - CpDry * xyr_VPTempBZ(i,j,k)  &
              &   * (                          &
              &       beta                     &
              &       * (                      &
              &           xyz_ExnerAs(i,j,k+1) &
              &         - xyz_ExnerAs(i,j,k)   &
              &         )                      &
              &     + (1.0d0 - beta)           &
              &       * (                      &
              &           xyz_ExnerNs(i,j,k+1) &
              &         - xyz_ExnerNs(i,j,k)   &
              &         )                      &
              &     ) / dz
            
          end do
        end do
      end do
      
      xyr_PGrad(:,:,kmax) = 0.0d0
      xyr_SWF(:,:,kmax)   = 0.0d0
      
    end subroutine PGrad_VI


    subroutine Exner_HEVI
      !
      !陰解法を用いたエクスナー関数の計算. 
      !

      !暗黙の型宣言禁止
      implicit none

      !作業変数定義
      real(DP)               :: D1(1:nx,1:ny,1:nz)  
      real(DP)               :: D(nx*ny,nz)
      real(DP)               :: E(1:nx,1:ny,0:nz)
      real(DP)               :: F(1:nx,1:ny,1:nz)
      real(DP)               :: F0(1:nx,1:ny,kmin:kmax-1)  
      real(DP)               :: dt ! 短い時間格子間隔
      integer                :: i, j, k
      
      real(DP)               :: X(M, N)     !定数/解行列
      real(DP)               :: TX(N, M)    !解行列を転置したもの
      integer                :: NRHS        
      integer                :: INFO
      integer                :: LDB
      character(1),parameter :: TRANS = 'N'
            
      ! Initialize
      !
      dt = DelTimeShort

      !---------------------------------------------------------------
      !行列計算のための係数

      !  添字の範囲は, 1:nx, 1:ny, 0:nz
      !  D(:,:,1) を求める時に D(:,:,0) の値が必要になるため. 
      !
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            
            E(i,j,k) =                            &
              & - ( 1.0d0 - beta )                &
              &   * (                             &
              &     + xyz_ExnerNs(i,j,k+1)        & !
              &     - xyz_ExnerNs(i,j,k)          & ! xyz => xyr
              &     ) / dz                        &
              & + (                               &
              &   + AlphaV                        &
              &     * (                           &
              &       + xyz_VelDivNs(i,j,k+1)     & !
              &       - xyz_VelDivNs(i,j,k)       & ! xyz => xyr
              &       ) / dz                      &
              &   + xyr_DVelZDtNl(i,j,k)          &
              &   ) / xyr_CpVPTempBZ(i,j,k)  
            
          end do
        end do
      end do
      
      ! 被微分関数
      !   配列 F0 の添字の範囲は, 1:nx, 1:ny, kmin:kmax
      !   配列 F を求める際に F0 を z 方向に微分するため. 
      !
      do k = kmin, kmax-1
        do j = 1, ny
          do i = 1, nx
            
            F0(i,j,k)  =                            &
              & + xyr_DensVPTempBZ(i,j,k)           &
              &   * (                               &
              &     + xyr_VelZNs(i,j,k)             & 
              &     - xyr_CpVPTempBZ(i,j,k)         &
              &       * (1.0d0 - beta)              &
              &       * (                           &
              &           xyz_ExnerNs(i,j,k+1)      &
              &         - xyz_ExnerNs(i,j,k)        &
              &         ) / dz * dt                 &
              &     + AlphaV                        &
              &       * (                           &
              &           xyz_VelDivNs(i,j,k+1)     &
              &         - xyz_VelDivNs(i,j,k)       &
              &         ) / dz * dt                 &
              &     + xyr_DVelZDtNl(i,j,k) * dt     &
              &     )
          end do
        end do
      end do
      
      !行列計算のための係数
      !  添字の範囲は, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            F(i,j,k) = &
              & - beta * xyz_F1BZ(i,j,k) * dt   &
              &   * (                           &
              &       F0(i,j,k)                 &
              &     - F0(i,j,k-1)               &
              &     ) / dz                      &
              & + xyz_DExnerDtNl(i,j,k) * dt 
            
          end do
        end do
      end do
    
      !行列計算のための係数
      !  添字の範囲は, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz 
        do j = 1, ny
          do i = 1, nx
            
            D1(i,j,k) =                                 &
              & + xyz_ExnerNs(i,j,k)                    &
              & - (1.0d0 - beta)                        &
              &   * xyz_F1BZ(i,j,k) * dt                &
              &   * (                                   &
              &       xyr_DensVPTempBZ(i,j,k)           &
              &       * xyr_VelZNs(i,j,k)               &
              &     - xyr_DensVPTempBZ(i,j,k-1)         &
              &       * xyr_VelZNs(i,j,k-1)             &
              &     ) / dz                              &
              & - (xyz_VelSW(i,j,k) ** 2.0d0) * dt      &
              &   / (CpDry * xyz_VPTempBZ(i,j,k))       &
              &   * (                                   &
              &     + (                                 &
              &         pyz_VelXAs(i,j,k)               &
              &       - pyz_VelXAs(i-1,j,k)             &
              &       ) / dx                            &
              &     + (                                 &
              &         xqz_VelYAs(i,j,k)               &
              &       - xqz_VelYAs(i,j-1,k)             &
              &       ) / dy                            &
              &     )                                   &
              & + F(i,j,k)
          
          end do
        end do
      end do
      
      ! 行列計算のための係数
      !
      do j = 1, ny
        do i = 1, nx

          D1(i,j,1) =                                    &
            & + D1(i,j,1)                                &
            & - beta * xyz_F1BZ(i,j,1) * (dt ** 2.0d0)   &
            &   * xyr_CpDensVPTemp2BZ(i,j,0)             &
            &   * E(i,j,0)                               &
            &   / dz
          
          D1(i,j,nz) =                                   &
            & + D1(i,j,nz)                               &
            & + beta * xyz_F1BZ(i,j,nz) * (dt ** 2.0d0)  &
            &   * xyr_CpDensVPTemp2BZ(i,j,nz)            &
            &   * E(i,j,nz)                              &
            &   / dz
        end do
      end do
      
      !-----------------------------------------------------------
      !連立一次方程式の解を求める
      
      !変数の初期化
      !
      NRHS = M
      INFO = 0
      LDB  = N
      
      ! LAPACK の仕様に合わせて変形 
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            D(i + nx * (j - 1), k) =  D1(i,j,k)
          end do
        end do
      end do
      
      TX = transpose( D )
      
      !解行列の計算. LAPACK を使用. 
      !
      call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)
      
      !解のコンディションをチェック. 
      !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

      !戻り値を出力
      !
      X = transpose( TX )
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyz_ExnerAs(i,j,k) = X(i + nx * (j - 1 ), k)
          end do
        end do
      end do

      xyz_ExnerAs(imin:0,:,:) = 0.0d0
      xyz_ExnerAs(:,jmin:0,:) = 0.0d0
      xyz_ExnerAs(:,:,kmin:0) = 0.0d0
      xyz_ExnerAs(nx+1:imax,:,:) = 0.0d0
      xyz_ExnerAs(:,ny+1:jmax,:) = 0.0d0
      xyz_ExnerAs(:,:,nz+1:kmax) = 0.0d0
      
      ! のり代に値を入れる
      !
      call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
      
    end subroutine Exner_HEVI
    
  end subroutine Dynamics3D_Short_forcing
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  subroutine Dynamics2D_Short_forcing(  &
        &  pyz_VelXNs,          & ! (in)
        &  xyr_VelZNs,          & ! (in)
        &  xyz_ExnerNs,         & ! (in)
        &  pyz_DVelXDtNl,       & ! (in)
        &  xyr_DVelZDtNl,       & ! (in)
        &  xyz_DExnerDtNl,      & ! (in)
        &  xyz_DExnerDtNs,      & ! (in)
        &  pyz_VelXAs,          & ! (out)
        &  xqz_VelYAs,          & ! (out)
        &  xyr_VelZAs,          & ! (out)
        &  xyz_ExnerAs          & ! (out)
        & )
    !
    ! 力学コア (短い時間ステップ) 2D 版
    !

    !モジュール読み込み
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gtool_historyauto, only: HistoryAutoPut 
    use gridset, only: &
      &                 imin,            &! x 方向の配列の下限
      &                 imax,            &! x 方向の配列の上限
      &                 jmin,            &! y 方向の配列の下限
      &                 jmax,            &! y 方向の配列の上限
      &                 kmin,            &! z 方向の配列の下限
      &                 kmax,            &! z 方向の配列の上限
      &                 nx,              &! x 方向の物理領域の上限
      &                 ny,              &! x 方向の物理領域の上限
      &                 nz              ! y 方向の物理領域の上限
    use constants,only: CpDry, CvDry, GasRDry ! 乾燥成分の比熱
    use timeset, only:  DelTimeShort, TimeN
    use axesset, only:  dx,  dz        ! 格子間隔
    use basicset, only: xyz_VelSW,      &!基本場の音速 
      &                 xyz_VPTempBZ,      &!基本場の温位
      &                 pyz_VPTempBZ,      &!基本場の温位
      &                 xyr_VPTempBZ      !基本場の温位
    use setmargin,only: SetMargin_xyzf, SetMargin_xyz, &
      &                 SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

    implicit none

    real(DP), intent(in)     :: pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)     :: xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout)  :: xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out)    :: xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_PGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_PGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_SWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_SWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !------------------------------------------------------------
    ! initialize: Divergence of velocity
    !
    call VelDivC2

    !------------------------------------------------------------
    ! VelX, VelY
    !  水平方向は陽解法で解く. 
    !
    call PGrad_HE

    ! tendency 
    !
    pyz_DVelXDtNs = pyz_PGrad + pyz_SWF

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * (pyz_DVelXDtNl + pyz_DVelXDtNs)

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)

    ! Y 方向には値がゼロ
    xqz_VelYAs = 0.0d0

    
    !------------------------------------------------------------
    ! Exner function
    !  積分値を返すことに注意. 
    !

    ! 短い時間ステップで評価する圧力の式の tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency の合計
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    ! エクスナー関数の計算
    !
    call Exner_HEVI

    !------------------------------------------------------------
    ! VelZ
    !
    call PGrad_VI

    ! tendency
    !
    xyr_DVelZDtNs = xyr_PGrad + xyr_SWF

    ! 値の保管
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_PGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_SWF(1:nx,1:ny,1:nz))

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNl + xyr_DVelZDtNs)

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  contains

    subroutine VelDivC2
      
      implicit none
      integer              :: i, j, k
      
      do  k = kmin + 1, kmax 
        do j = 1, ny
          do i = imin + 1, imax 
            xyz_VelDivNs(i,j,k) =          &
              & + (                        &
              &     pyz_VelXNs(i,j,k)      &
              &   - pyz_VelXNs(i-1,j,k)    &
              &   ) / dx                   &
              & + (                        &
              &     xyr_VelZNs(i,j,k)      &
              &   - xyr_VelZNs(i,j,k-1)    &
              &   ) / dz
          end do
        end do
      end do
      
      xyz_VelDivNs(imin,:,:) = 0.0d0 
      xyz_VelDivNs(:,:,kmin) = 0.0d0 
      
    end subroutine VelDivC2


    subroutine PGrad_HE
      
      implicit none
      integer              :: i, j, k
      
      !------------------------------------------------------------------
      ! X 方向
      
      do k = kmin, kmax
        do j = 1, ny          
          do i = imin, imax - 1

            ! 音波減衰項
            !            
            pyz_SWF(i,j,k) =                  &
              &   AlphaH                      &
              &   * (                         &
              &       xyz_VelDivNs(i+1,j,k)   &
              &     - xyz_VelDivNs(i,j,k)     &
              &     ) / dx            
            
            ! 圧力傾度力
            !
            pyz_PGrad(i,j,k) =                &
              & - CpDry * pyz_VPTempBZ(i,j,k) &
              &   * (                         &
              &       xyz_ExnerNs(i+1,j,k)    &
              &     - xyz_ExnerNs(i,j,k)      &
              &     ) / dx                                
            
          end do
        end do
      end do
      
      ! 穴埋め
      !
      pyz_SWF(imax,:,:)   = 0.0d0 
      pyz_PGrad(imax,:,:) = 0.0d0 
      
    end subroutine PGrad_HE


    subroutine PGrad_VI
    
      implicit none
      integer               :: i, j, k

      do k = kmin, kmax - 1
        do j = 1, ny          
          do i = imin, imax

            ! 音波減衰項
            !            
            xyr_SWF(i,j,k) =                  &
              & + AlphaV                      &
              &   * (                         &
              &       xyz_VelDivNs(i,j,k+1)   & 
              &     - xyz_VelDivNs(i,j,k)     & 
              &     ) / dz
            
            ! 圧力傾度力
            !
            xyr_PGrad(i,j,k) =                 &
              & - CpDry * xyr_VPTempBZ(i,j,k)  &
              &   * (                          &
              &       beta                     &
              &       * (                      &
              &           xyz_ExnerAs(i,j,k+1) &
              &         - xyz_ExnerAs(i,j,k)   &
              &         )                      &
              &     + (1.0d0 - beta)           &
              &       * (                      &
              &           xyz_ExnerNs(i,j,k+1) &
              &         - xyz_ExnerNs(i,j,k)   &
              &         )                      &
              &     ) / dz
            
          end do
        end do
      end do
      
      xyr_PGrad(:,:,kmax) = 0.0d0
      xyr_SWF(:,:,kmax)   = 0.0d0
      
    end subroutine PGrad_VI


    subroutine Exner_HEVI
      !
      !陰解法を用いたエクスナー関数の計算. 
      !

      !暗黙の型宣言禁止
      implicit none

      !作業変数定義
      real(DP)               :: D1(1:nx,1:ny,1:nz)  
      real(DP)               :: D(nx*ny,nz)
      real(DP)               :: E(1:nx,1:ny,0:nz)
      real(DP)               :: F(1:nx,1:ny,1:nz)
      real(DP)               :: F0(1:nx,1:ny,kmin:kmax-1)  
      real(DP)               :: dt ! 短い時間格子間隔
      integer                :: i, j, k
      
      real(DP)               :: X(M, N)     !定数/解行列
      real(DP)               :: TX(N, M)    !解行列を転置したもの
      integer                :: NRHS        
      integer                :: INFO
      integer                :: LDB
      character(1),parameter :: TRANS = 'N'
      
      
      ! Initialize
      !
      dt = DelTimeShort

      !---------------------------------------------------------------
      !行列計算のための係数

      !  添字の範囲は, 1:nx, 1:ny, 0:nz
      !  D(:,:,1) を求める時に D(:,:,0) の値が必要になるため. 
      !
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            
            E(i,j,k) =                            &
              & - ( 1.0d0 - beta )                &
              &   * (                             &
              &     + xyz_ExnerNs(i,j,k+1)        & !
              &     - xyz_ExnerNs(i,j,k)          & ! xyz => xyr
              &     ) / dz                        &
              & + (                               &
              &   + AlphaV                        &
              &     * (                           &
              &       + xyz_VelDivNs(i,j,k+1)     & !
              &       - xyz_VelDivNs(i,j,k)       & ! xyz => xyr
              &       ) / dz                      &
              &   + xyr_DVelZDtNl(i,j,k)          &
              &   ) / xyr_CpVPTempBZ(i,j,k)  
            
          end do
        end do
      end do
      
      ! 被微分関数
      !   配列 F0 の添字の範囲は, 1:nx, 1:ny, kmin:kmax
      !   配列 F を求める際に F0 を z 方向に微分するため. 
      !
      do k = kmin, kmax-1
        do j = 1, ny
          do i = 1, nx
            
            F0(i,j,k)  =                            &
              & + xyr_DensVPTempBZ(i,j,k)           &
              &   * (                               &
              &     + xyr_VelZNs(i,j,k)             & 
              &     - xyr_CpVPTempBZ(i,j,k)         &
              &       * (1.0d0 - beta)              &
              &       * (                           &
              &           xyz_ExnerNs(i,j,k+1)      & !
              &         - xyz_ExnerNs(i,j,k)        & ! xyz => xyr
              &         ) / dz * dt                 &
              &     + AlphaV                        &
              &       * (                           &
              &           xyz_VelDivNs(i,j,k+1)     &
              &         - xyz_VelDivNs(i,j,k)       &
              &         ) / dz * dt                 &
              &     + xyr_DVelZDtNl(i,j,k) * dt     &
              &     )
          end do
        end do
      end do
      
      !行列計算のための係数
      !  添字の範囲は, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            
            F(i,j,k) = &
              & - beta * xyz_F1BZ(i,j,k) * dt   &
              &   * (                           &
              &       F0(i,j,k)                 &
              &     - F0(i,j,k-1)               &
              &     ) / dz                      &
              & + xyz_DExnerDtNl(i,j,k) * dt 
            
          end do
        end do
      end do
    
      !行列計算のための係数
      !  添字の範囲は, 1:nx, 1:ny, 1:nz
      !
      do k = 1, nz 
        do j = 1, ny
          do i = 1, nx
            
            D1(i,j,k) =                                 &
              & + xyz_ExnerNs(i,j,k)                    &
              & - (1.0d0 - beta)                        &
              &   * xyz_F1BZ(i,j,k) * dt                &
              &   * (                                   &
              &       xyr_DensVPTempBZ(i,j,k)           &
              &       * xyr_VelZNs(i,j,k)               &
              &     - xyr_DensVPTempBZ(i,j,k-1)         &
              &       * xyr_VelZNs(i,j,k-1)             &
              &     ) / dz                              &
              & - (xyz_VelSW(i,j,k) ** 2.0d0) * dt      &
              &   / (CpDry * xyz_VPTempBZ(i,j,k))       &
              &   * (                                   &
              &     + (                                 &
              &         pyz_VelXAs(i,j,k)               &
              &       - pyz_VelXAs(i-1,j,k)             &
              &       ) / dx                            &
              &     )                                   &
              & + F(i,j,k)
          
          end do
        end do
      end do
      
      ! 行列計算のための係数
      !
      do j = 1, ny
        do i = 1, nx

          D1(i,j,1) =                                    &
            & + D1(i,j,1)                                &
            & - beta * xyz_F1BZ(i,j,1) * (dt ** 2.0d0)   &
            &   * xyr_CpDensVPTemp2BZ(i,j,0)             &
            &   * E(i,j,0)                               &
            &   / dz
          
          D1(i,j,nz) =                                   &
            & + D1(i,j,nz)                               &
            & + beta * xyz_F1BZ(i,j,nz) * (dt ** 2.0d0)  &
            &   * xyr_CpDensVPTemp2BZ(i,j,nz)            &
            &   * E(i,j,nz)                              &
            &   / dz
        end do
      end do
      
      !-----------------------------------------------------------
      !連立一次方程式の解を求める
      
      !変数の初期化
      !
      NRHS = M
      INFO = 0
      LDB  = N
      
      ! LAPACK の仕様に合わせて変形 
      !
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            D(i + nx * (j - 1), k) =  D1(i,j,k)
          end do
        end do
      end do
      
      TX = transpose( D )
      
      !解行列の計算. LAPACK を使用. 
      !
      call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)
      
      !解のコンディションをチェック. 
      !
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if

      !戻り値を出力
      !
      X = transpose( TX )
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyz_ExnerAs(i,j,k) = X(i + nx * (j - 1 ), k)
          end do
        end do
      end do
      xyz_ExnerAs(imin:0,:,:) = 0.0d0
      xyz_ExnerAs(:,:,kmin:0) = 0.0d0
      xyz_ExnerAs(nx+1:imax,:,:) = 0.0d0
      xyz_ExnerAs(:,:,nz+1:kmax) = 0.0d0
      
      ! のり代に値を入れる
      !
      call SetMargin_xyz( xyz_ExnerAs ) ! (inout)
      
    end subroutine Exner_HEVI
    
  end subroutine Dynamics2D_Short_forcing

  
!!!--------------------------------------------------------------------!!!
  subroutine Dynamics_VI_init()
    !
    ! 力学コア 陰解法部分の初期化ルーチン
    ! * エクスナー関数を陰解法で解く際に必要となる, 係数行列の要素を決め, 
    !   LU 分解を行う. 
    !

    !モジュール呼び出し
    use dc_message, only : MessageNotify
    use gridset, only:  kmin,            &! z 方向の配列の下限
      &                 kmax,            &! z 方向の配列の上限
      &                 nx,              &! x 方向の物理領域の上限
      &                 ny,              &! x 方向の物理領域の上限
      &                 nz                ! y 方向の物理領域の上限
    use axesset, only:  dz
    use constants,only: CpDry             ! 乾燥成分の比熱
    use timeset, only : DelTimeShort
    use basicset, only: xyz_VelSW,       &!基本場の音速 
      &                 xyz_DensBZ,      &!基本場の密度
      &                 xyz_VPTempBZ      !基本場の温位

    !暗黙の型宣言禁止
    implicit none

    real(DP)  :: dt      ! 短い時間格子
    integer   :: INFO    !解のコンディションチェック
    integer   :: i, j, k

    !----------------------------------------------------------------
    ! 初期化

    ! 変数名が長すぎたので, 名前を置き換える
    !
    dt = DelTimeShort

    ! 配列の割り付け
    !
    allocate( A(1:nz) )
    allocate( B(2:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_F1BZ(1:nx,1:ny,1:nz) )
    allocate( xyr_CpDensVPTemp2BZ(1:nx,1:ny,kmin:kmax) )
    allocate( xyr_DensVPTempBZ(1:nx,1:ny,kmin:kmax) )
    allocate( xyr_CpVPTempBZ(1:nx,1:ny,kmin:kmax) )

    !----------------------------------------------------------------
    ! 係数行列および共通して利用される配列の値を決める
    !   A, B, C を求める際, 基本場 (BZ) の量は X 方向に一様なので. 
    !   nx, ny の値で代表させることとした. 
    !
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyz_F1BZ(i,j,k) =                                                    &
            &  ( xyz_VelSW(i,j,k) ** 2.0d0 )                                   &
            &  / ( CpDry * xyz_DensBZ(i,j,k) * (xyz_VPTempBZ(i,j,k) ** 2.0d0) )
        end do
      end do
    end do

    do k = kmin, kmax - 1
      do j = 1, ny        
        do i = 1, nx
          xyr_CpDensVPTemp2BZ(i,j,k)=                    &
            &  CpDry                                     &
            &  * (                                       &
            &    + xyz_DensBZ(i,j,k+1)                   &
            &      * ( xyz_VPTempBZ(i,j,k+1) ** 2.0d0 )  &
            &    + xyz_DensBZ(i,j,k)                     &
            &      * ( xyz_VPTempBZ(i,j,k) ** 2.0d0 )    &
            &    ) * 5.0d-1
        end do
      end do
    end do
    xyr_CpDensVPTemp2BZ(:,:,kmax) = 0.0d0  !穴埋め
    
    do k = kmin, kmax-1
      do j = 1, ny
        do i = 1, nx
          xyr_DensVPTempBZ(i,j,k) =                             &
            & + (                                               &
            &   + xyz_DensBZ(i,j,k+1) * xyz_VPTempBZ(i,j,k+1)   &
            &   + xyz_DensBZ(i,j,k)   * xyz_VPTempBZ(i,j,k)     &
            &   ) * 5.0d-1    
        end do
      end do
    end do
    xyr_DensVPTempBZ(:,:,kmax) = 0.0d0  !穴埋め

    do k = kmin, kmax-1
      do j = 1, ny
        do i = 1, nx
          xyr_CpVPTempBZ(i,j,k) =          &
            &   CpDry                      &
            &   * (                        &
            &     + xyz_VPTempBZ(i,j,k+1)  &
            &     + xyz_VPTempBZ(i,j,k)    &
            &     ) * 5.0d-1
        end do
      end do
    end do
    xyr_CpVPTempBZ(:,:,kmax) = 0.0d0
          
    do k = 2, nz-1
      A(k) =                                        &
        & + 1.0d0                                   &
        & + ( beta ** 2.0d0 )                       &
        &    * xyz_F1BZ(nx,ny,k) * ( dt * dt )      &
        &    * (                                    &
        &         xyr_CpDensVPTemp2BZ(nx,ny,k)      &
        &       + xyr_CpDensVPTemp2BZ(nx,ny,k-1)    &
        &       )                                   &
        &    / ( dz * dz )
    end do

    A(1) =                                   &
      & + 1.0d0                              &
      & + ( beta ** 2.0d0 )                  &
      &   * xyz_F1BZ(nx,ny,1) * ( dt * dt )  &
      &   * xyr_CpDensVPTemp2BZ(nx,ny,1)     &
      &   / ( dz * dz ) 

    A(nz) =                                  &
      & + 1.0d0                              &
      & + ( beta ** 2.0d0 )                  &
      &   * xyz_F1BZ(nx,ny,nz) * ( dt * dt ) &
      &   * xyr_CpDensVPTemp2BZ(nx,ny,nz-1)  &
      &   / ( dz * dz )  

    do k = 2, nz
      B(k) =                                     &
        & - ( beta ** 2.0d0 )                    &
        &   * xyz_F1BZ(nx,ny,k-1) * ( dt * dt )  &
        &   * xyr_CpDensVPTemp2BZ(nx,ny,k-1)     &
        &   / ( dz * dz )
    end do
    
    do k = 1, nz-1
      C(k) =                                     &
        & - ( beta ** 2.0d0 )                    &
        &   * xyz_F1BZ(nx,ny,k+1) * (dt * dt )   &
        &   * xyr_CpDensVPTemp2BZ(nx,ny,k)       &
        &   / ( dz * dz )
    end do

    !----------------------------------------------------------------
    ! 係数行列を LU 分解
    !
    !配列の大きさを保管
    N    = nz                   !係数行列/改行列の次数, 整合寸法
    M    = nx * ny              !方程式の組数
    NUD  = 1                    !係数行列の上三角部分の帯幅
    NLD  = 1                    !係数行列の下三角部分の帯幅
    NAL  = NLD                  !LU 分解の結果 L の整合寸法
    NA   = NUD + NLD + 1
    INFO = 0

    ! 配列の割り当て
    !
    allocate( AL1(N), IP(N) )

    ! 解行列の計算. LAPACK を使用. 
    !
    call DGTTRF(N, C, A, B, AL1, IP, INFO)
    
    ! 解のコンディションをチェック. 
    !
    if (INFO /= 0) then
      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
      stop
    end if
    
  end subroutine Dynamics_VI_init
  


  subroutine Dynamics_Tendency_Output
    !
    ! ファイル出力の定義
    !

    !モジュール呼び出し
    use gtool_historyauto, only: HistoryAutoAddVariable
    use composition,       only: SpcWetSymbol
    use gridset,           only: ncmax          ! 物質数

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer :: f

    call HistoryAutoAddVariable(  &
      & varname='DPTempDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of potential temperature',  &
      & units='K.s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of potential temperature',&
      & units='K.s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of exner function',  &
      & units='s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DExnerDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of exner function',&
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='CDensAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of cloud density',  &
      & units='kg.m-3.s-1',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='CDensDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of cloud density',&
      & units='kg.m-3.s-1',    &
      & xtype='double')

    do f = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtAdv', &
        & dims=(/'x','y','z','t'/),     &
        & longname='Advection term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
      
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtDiff', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Diffusion term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(f))//'DtFall', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Fall term of '          &
        &           //trim(SpcWetSymbol(f))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

    end do

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (x)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelXDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of velocity (x)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (x)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelXDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (x)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (y)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelYDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of velocity (y)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (y)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelYDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (y)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection term of velocity (z)',  &
      & units='m.s-2',    &
      & xtype='double')
    
    call HistoryAutoAddVariable(  &
      & varname='DVelZDtDiff',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Numerical diffusion term of Velocity (z)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyT',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (Temperature)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyM',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (MolWt)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtBuoyD',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Buoyancy (Drag)',&
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtPGrad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Pressure gradient term of velocity (z)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DVelZDtSWF', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Filter for acoustic mode (z)',  &
      & units='m.s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='VelDiv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='velocity divergence',  &
      & units='s-2',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtAdv', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Advection of Km',  &
      & units='s-1',    &
      & xtype='double')

    call HistoryAutoAddVariable(  &
      & varname='DKmDtDiff', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Diffusion term of Km',  &
      & units='s-1',    &
      & xtype='double')

  end subroutine Dynamics_Tendency_Output

  
end module DynamicsHEVI_v2

