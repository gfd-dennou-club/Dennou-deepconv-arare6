!= Module DynamicsHEVI
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro), ODAKA Masatsugu 
! Version::   $Id: dynamicshevi.f90,v 1.21 2014/05/28 15:27:42 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module DynamicsHEVI
  !
  ! 力学コア. 
  !   タイムスプリット法を利用. 音波モードとそれ以外を別々の時間刻みで解く. 
  !   短い時間ステップの計算には HE-VI 法を利用.
  !
  ! Note: 
  !  * エクスナー関数の空間方向の離散化は 2 次精度であるため, 気圧傾度
  !    力項の計算プログラムにおいて differentiate_center4 モジュールを
  !    指定することはできないので注意.

  !モジュール読み込み
  use dc_types,   only : DP

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private
  
  real(DP), save, private  :: beta  = 5.0d-1     !クランクニコルソン法なら 0.5
                                                 !完全陰解法なら 1
  integer, save, private   :: N = 10             !係数行列/改行列の次数, 整合寸法
  integer, save, private   :: M = 10             !方程式の組数
  integer, save, private   :: NUD = 1            !係数行列の上三角部分の帯幅
  integer, save, private   :: NLD = 1            !係数行列の下三角部分の帯幅
  integer, save, private   :: NAL = 1            !LU 分解の結果 L の整合寸法
  integer, save, private   :: NA = 3             !NUD + NLD + 1

  real(DP), allocatable, save, private :: xyz_F1BZ(:,:,:)
                                                 !係数行列の計算に用いる配列
  real(DP), allocatable, save, private :: xyr_F2BZ(:,:,:)
                                                 !係数行列の計算に用いる配列

  real(DP), allocatable, save, private :: A(:)   !係数行列の対角成分
  real(DP), allocatable, save, private :: B(:)   !係数行列の上三角部分
  real(DP), allocatable, save, private :: C(:)   !係数行列の下三角部分
  real(DP), allocatable, save, private :: AL1(:) !LU 分解の結果 L (1 次元配列)
  integer, allocatable, save, private  :: IP(:)  !部分ピボット交換の情報を格納

  real(DP), save, private :: AlphaH = 0.0d0      !音波減衰項の減衰係数
  real(DP), save, private :: AlphaV = 0.0d0      !音波減衰項の減衰係数
  real(DP), save, private :: NuHh   = 0.0d0      !数値粘性の係数 (水平方向)
  real(DP), save, private :: NuVh   = 0.0d0      !数値粘性の係数 (鉛直方向)
  real(DP), save, private :: NuHm   = 0.0d0      !数値粘性の係数 (水平方向)
  real(DP), save, private :: NuVm   = 0.0d0      !数値粘性の係数 (鉛直方向)

  character(*), parameter :: module_name = 'DynamicHEVI'
                                                 ! モジュールの名称.
                                                 ! Module name
  real(DP), save, private :: FactorBuoyTemp    = 1.0d0
                                                 !浮力 (温度の寄与) の有無
                                                 !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorBuoyMolWt   = 1.0d0
                                                 !浮力 (分子量効果) の有無
                                                 !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorBuoyLoading = 1.0d0
                                                 !浮力 (荷重効果) の有無
                                                 !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtAdv    = 1.0d0 
                                                 !エクスナー関数の移流の有無
                                                 !考慮しない場合は値をゼロにする.
  real(DP), save, private :: FactorDExnerDtExpnd  = 1.0d0
                                                 !エクスナー関数の膨張項の有無
                                                 !考慮しない場合は値をゼロにする. 

  real(DP), allocatable, save,private :: xyr_DPTempBZDz(:,:,:)   !基本場の鉛直微分
!  real(DP), allocatable, save,private :: xyr_DExnerBZDz(:,:,:)   !基本場の鉛直微分
  real(DP), allocatable, save,private :: xyrf_DQMixBZDz(:,:,:,:) !基本場の鉛直微分

  !public 
  public Dynamics_Init
  public Dynamics_Long_forcing
  public Dynamics_Short_forcing

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Init

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
    real(DP)  :: AlphaNDiff  = 1.0d-3 !4次の数値拡散の係数. CReSS マニュアルより
    real(DP)  :: NDiffRatio = 1.0d0   !速度に対する粘性を上げる場合は数字を 1 以上にする. 
    integer   :: unit            !装置番号
    integer   :: f

    !-------------------------------------------------------------------
    ! Namelist から情報を取得する
    !
    NAMELIST /Dynamics_nml/                                    &
         & AlphaSound, AlphaNDiff, NDiffRatio, beta,           &
         & FactorBuoyTemp, FactorBuoyMolWt, FactorBuoyLoading, &
         & FactorDExnerDtAdv, FactorDExnerDtExpnd

    !ファイルオープン. 情報取得. 
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

    ! 陰解法の計算設定の初期化
    !
    call DynamicsVI_init()

    ! tendency の出力
    !
    call Dynamics_Tendency_output

    ! 配列の用意
    !
    allocate( xyr_DPTempBZDz(imin:imax, jmin:jmax, kmin:kmax) )
!    allocate( xyr_DExnerBZDz(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyrf_DQMixBZDz(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    xyr_DPTempBZDz = xyr_dz_xyz( xyz_PTempBZ )
!    xyr_DExnerBZDz = xyr_dz_xyz( xyz_ExnerBZ )
    do f = 1, ncmax
      xyrf_DQMixBZDz(:,:,:,f) = xyr_dz_xyz( xyzf_QMixBZ(:,:,:,f) )
    end do

  end subroutine Dynamics_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Dynamics_Long_forcing(        &
    & pyz_VelXBl,  pyz_VelXNl,    & ! (in)
    & xqz_VelYBl,  xqz_VelYNl,    & ! (in)
    & xyr_VelZBl,  xyr_VelZNl,    & ! (in)
    & xyz_PTempBl, xyz_PTempNl,   & ! (in)
    & xyz_ExnerBl, xyz_ExnerNl,   & ! (in)
    & xyzf_QMixBl, xyzf_QMixNl,   & ! (in)
    & xyz_KmBl,    xyz_KmNl,      & ! (in)
    & pyz_DVelXDtNl,              & ! (inout)
    & xqz_DVelYDtNl,              & ! (inout)
    & xyr_DVelZDtNl,              & ! (inout)
    & xyz_DPTempDtNl,             & ! (inout)
    & xyz_DExnerDtNl,             & ! (inout)
    & xyzf_DQMixDtNl,             & ! (inout)
    & xyz_DKmDtNl                 & ! (inout)
    & )
    !
    ! 力学コア: 長い時間ステップで評価する項の計算.
    !

    !モジュール読み込み
    use dc_types,  only : DP
    use gtool_historyauto, &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x 方向の配列サイズ
      &                     jmin, jmax,       &! y 方向の配列サイズ
      &                     kmin, kmax,       &! z 方向の配列サイズ
      &                     nx, ny, nz,       &! 物理領域のサイズ
      &                     ncmax              ! 物質数
    use timeset,     only : TimeN
    use composition, only : SpcWetSymbol
    use basicset,    only : xyz_PTempBZ,       & ! 温位の基本場
      &                     xyzf_QMixBZ,       & ! 混合比の基本場
      &                     xyr_QMixBZ,        &
      &                     xyr_QMixBZPerMolWt
    use average,     only : pyz_xyz, pyz_pqz, pqz_xqz, pyz_pyr, pyr_xyr, &
      &                     xqz_pqz, pqz_pyz, xqz_xyz, xqz_xqr, xqr_xyr, &
      &                     xyr_pyr, pyr_pyz, xyr_xqr, xqr_xqz, xyr_xyz, &
      &                     xyz_pyz, xyz_xqz, xyz_xyr
    use differentiate_center4,                                                                              &
      &              only : pyz_c4dx_xyz => pyz_dx_xyz, xqz_c4dy_xyz => xqz_dy_xyz, xyr_c4dz_xyz => xyr_dz_xyz, &
      &                     xyz_c4dx_pyz => xyz_dx_pyz, pqz_c4dy_pyz => pqz_dy_pyz, pyr_c4dz_pyz => pyr_dz_pyz, &
      &                     pqz_c4dx_xqz => pqz_dx_xqz, xyz_c4dy_xqz => xyz_dy_xqz, xqr_c4dz_xqz => xqr_dz_xqz, &
      &                     pyr_c4dx_xyr => pyr_dx_xyr, xqr_c4dy_xyr => xqr_dy_xyr, xyz_c4dz_xyr => xyz_dz_xyr
    use differentiate_center2,                              &
      &              only : pyz_dx_xyz, xyz_dx_pyz, pyz_dy_pqz, &
      &                     pqz_dy_pyz, pyz_dz_pyr, pyr_dz_pyz, &
      &                     xqz_dx_pqz, pqz_dx_xqz, xqz_dy_xyz, &
      &                     xyz_dy_xqz, xqz_dz_xqr, xqr_dz_xqz, &
      &                     xyr_dx_pyr, pyr_dx_xyr, xyr_dy_xqr, &
      &                     xqr_dy_xyr, xyr_dz_xyz, xyz_dz_xyr, &
      &                     xyz_dx_pyz, pyz_dx_xyz, xyz_dy_xqz, &
      &                     xqz_dy_xyz, xyz_dz_xyr, xyr_dz_xyz
    use composition,  only: GasNum,       &! 
      &                     IdxG,         &!
      &                     MolWtWet       ! 湿潤成分の分子量
    use constants,    only: MolWtDry,     &! 乾燥成分の分子量
      &                     Grav           ! 重力加速度



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

    real(DP)                :: pyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: pyz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xqz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyT(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyM(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyr_BuoyD(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_Adv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyz_NDiff(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)                :: xyzf_Adv(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_Fall(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_NDiff(imin:imax,jmin:jmax,kmin:kmax, 1:ncmax)
    real(DP)                :: xyzf_QMixPerMolWt(imin:imax,jmin:jmax,kmin:kmax, 1:GasNum)
    integer                 :: f

    !------------------------------
    ! tendency of Km
    ! 

    ! Advection term
    !
    xyz_Adv  =                                             &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_KmNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_KmNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_KmNl ))    

    ! Numerical diffusion term 
    !
    xyz_NDiff =                                                                &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_KmBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_KmBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_KmBl ))))) 

    ! tendency
    !
    xyz_DKmDtNl = xyz_DKmDtNl + ( xyz_NDiff + xyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DKmDtAdv',    xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DKmDtDiff',  xyz_NDiff(1:nx,1:ny,1:nz))


    !------------------------------
    ! tendency of potential temperature
    ! 

    ! Advection term
    !
    xyz_Adv = &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_PTempNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_PTempNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_PTempNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_DPTempBZDz )

    ! numerical diffusion term
    !
    xyz_NDiff = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_PTempBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_PTempBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_PTempBl ))))) 

    ! sum
    !
    xyz_DPTempDtNl = xyz_DPTempDtNl + ( xyz_NDiff + xyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtDiff', xyz_NDiff(1:nx,1:ny,1:nz))

    
    !--------------------------------------
    ! Exner function
    !

    ! フラックス項の計算. 4 次精度中心差分を用いて計算
    !
    xyz_Adv = &
      & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyz_ExnerNl ))  &
      & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyz_ExnerNl ))  &
      & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyz_ExnerNl ))
!      & - xyz_xyr( xyr_VelZNl * xyr_DExnerBZDz ) 
    
    ! numerical diffusion term
    !
    xyz_NDiff = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyz_ExnerBl ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyz_ExnerBl ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyz_ExnerBl ))))) 
    xyz_NDiff = 0.0d0

    ! sum
    !
    xyz_DExnerDtNl = xyz_DExnerDtNl + ( xyz_NDiff + xyz_Adv ) * FactorDExnerDtAdv

    ! output
    !
    call HistoryAutoPut(TimeN, 'DExnerDtAdv',  xyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtDiff', xyz_NDiff(1:nx,1:ny,1:nz))


    !------------------------------
    ! tendency of mixing ratio
    ! 

    do f = 1, ncmax
      ! Advection term
      !
      xyzf_Adv(:,:,:,f) = &
        & - xyz_pyz( pyz_VelXNl * pyz_c4dx_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xqz( xqz_VelYNl * xqz_c4dy_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xyr( xyr_VelZNl * xyr_c4dz_xyz( xyzf_QMixNl(:,:,:,f) )) &
        & - xyz_xyr( xyr_VelZNl * xyrf_DQMixBZDz(:,:,:,f) )

      ! numerical diffusion term
      !
      xyzf_NDiff(:,:,:,f) = &
        &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyzf_QMixBl(:,:,:,f) ))))) &
        &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyzf_QMixBl(:,:,:,f) ))))) &
        &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyzf_QMixBl(:,:,:,f) )))))

    end do
    
    ! 落下項
    !
    call QMixFall

    ! sum
    !
    xyzf_DQMixDtNl = xyzf_DQMixDtNl + ( xyzf_NDiff + xyzf_Adv + xyzf_Fall )

    ! output
    !
    do f = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtAdv',  xyzf_Adv(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtDiff', xyzf_NDiff(1:nx,1:ny,1:nz,f))
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(f))//'DtFall', xyzf_Fall(1:nx,1:ny,1:nz,f))
    end do

    !------------------------------
    ! tendency of VelX
    ! 

    ! Advection term
    !
    pyz_Adv  = &
      & - pyz_VelXNl * pyz_xyz( xyz_c4dx_pyz( pyz_VelXNl ) )            &
      & - pyz_pqz( pqz_xqz( xqz_VelYNl ) * pqz_c4dy_pyz( pyz_VelXNl ) ) &
      & - pyz_pyr( pyr_xyr( xyr_VelZNl ) * pyr_c4dz_pyz( pyz_VelXNl ) )

    ! Numerical diffusion term 
    !
    pyz_NDiff = &
      & - NuHm * ( pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz( pyz_VelXBl ))))) &
      & - NuHm * ( pyz_dy_pqz(pqz_dy_pyz(pyz_dy_pqz(pqz_dy_pyz( pyz_VelXBl ))))) &
      & - NuVm * ( pyz_dz_pyr(pyr_dz_pyz(pyz_dz_pyr(pyr_dz_pyz( pyz_VelXBl )))))

    ! sum
    !
    pyz_DVelXDtNl = pyz_DVelXDtNl + ( pyz_NDiff + pyz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelXDtAdv',   pyz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtDiff', pyz_NDiff(1:nx,1:ny,1:nz))

    !------------------------------
    ! tendency of VelY
    !     

    ! Advection term
    !
    xqz_Adv  = &
      & - xqz_pqz( pqz_pyz( pyz_VelXNl ) * pqz_c4dx_xqz( xqz_VelYNl ) ) &
      & - xqz_VelYNl * xqz_xyz( xyz_c4dy_xqz( xqz_VelYNl ) ) &
      & - xqz_xqr( xqr_xyr( xyr_VelZNl ) * xqr_c4dz_xqz( xqz_VelYNl ) )

    ! Numerical diffusion term
    !
    xqz_NDiff = &
      & - NuHm * ( xqz_dx_pqz(pqz_dx_xqz(xqz_dx_pqz(pqz_dx_xqz( xqz_VelYBl ))))) &
      & - NuHm * ( xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz( xqz_VelYBl ))))) &
      & - NuVm * ( xqz_dz_xqr(xqr_dz_xqz(xqz_dz_xqr(xqr_dz_xqz( xqz_VelYBl )))))

    ! sum
    !
    xqz_DVelYDtNl = xqz_DVelYDtNl + ( xqz_NDiff + xqz_Adv )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelYDtAdv',   xqz_Adv(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtDiff', xqz_NDiff(1:nx,1:ny,1:nz))

    !------------------------------
    ! tendency of VelZ
    ! 

    ! 作業変数
    !
    do f = 1, GasNum
      xyzf_QMixPerMolWt(:,:,:,f) = xyzf_QMixNl(:,:,:,IdxG(f)) / MolWtWet(IdxG(f))
    end do

    ! Advection term
    !
    xyr_Adv  = &
      & - xyr_pyr( pyr_pyz( pyz_VelXNl ) * pyr_c4dx_xyr( xyr_VelZNl ) ) &
      & - xyr_xqr( xqr_xqz( xqz_VelYNl ) * xqr_c4dy_xyr( xyr_VelZNl ) ) &
      & - xyr_VelZNl * xyr_xyz( xyz_c4dz_xyr( xyr_VelZNl ) )

    ! Numerical diffusion term
    !
    xyr_NDiff = &
      & - NuHm * ( xyr_dx_pyr(pyr_dx_xyr(xyr_dx_pyr(pyr_dx_xyr( xyr_VelZBl ))))) &
      & - NuHm * ( xyr_dy_xqr(xqr_dy_xyr(xyr_dy_xqr(xqr_dy_xyr( xyr_VelZBl ))))) & 
      & - NuVm * ( xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr( xyr_VelZBl )))))

    ! Buoyancy due to temperature disturbunce
    !
    xyr_BuoyT = Grav * xyr_xyz( xyz_PTempNl / xyz_PTempBZ)

    ! Buoyancy due to molecular weight
    !
    xyr_BuoyM =                                                 &
      & + Grav * xyr_xyz( sum(xyzf_QMixPerMolWt, 4) )           &
      &    / ( 1.0d0 / MolWtDry + xyr_QMixBZPerMolWt )          &
      & - Grav * xyr_xyz( sum(xyzf_QMixNl(:,:,:,1:GasNum), 4) ) &
      &    / ( 1.0d0 + xyr_QmixBZ ) 

    ! Buoyancy due to loading
    !
    xyr_BuoyD =                                                       &
      & - Grav * xyr_xyz( sum(xyzf_QMixNl(:,:,:,GasNum+1:ncmax), 4) ) &
      &    / ( 1.0d0 + xyr_QMixBZ )

    ! sum
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

    ! output
    !
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
      use constants,   only : FactorJ
      use composition, only : IdxR, RainNum
      use basicset,    only : xyz_DensBZ
      
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
        xyzf_Fall(:,:,:,iR) =                                      &
          &  - xyz_c4dz_xyr( xyrf_QMixFallFlux(:,:,:,iR) ) / xyz_DensBZ  
        
      end do
            
    end subroutine QmixFall

  end subroutine Dynamics_Long_forcing
 

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
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 
    use gridset,     only : imin, imax,       &! x 方向の配列サイズ
      &                     jmin, jmax,       &! y 方向の配列サイズ
      &                     kmin, kmax,       &! z 方向の配列サイズ
      &                     nx, ny, nz         ! 物理領域のサイズ
    use timeset,     only : TimeN, DelTimeShort
    use constants,   only : CpDry, CvDry, GasRDry ! 乾燥成分の比熱
    use basicset,    only : pyz_VPTempBZ,     &! 基本場の温位
      &                     xqz_VPTempBZ,     &! 基本場の温位
      &                     xyr_VPTempBZ       ! 基本場の温位
!    use basicset,    only : xyz_VPTempBZ       !基本場の仮温位
    use average,     only : xyr_xyz, xqz_xyz, pyz_xyz
    use differentiate_center2,                &
      &              only : xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                     xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                     pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz
    use setmargin,  only : SetMargin_xyzf, SetMargin_xyz, &
      &                    SetMargin_pyz, SetMargin_xqz, SetMargin_xyr

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

    real(DP) :: pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtPGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtPGrad(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtPGrad(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: pyz_DVelXDtSWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xqz_DVelYDtSWF(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyr_DVelZDtSWF(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_DExnerDtExpnd(imin:imax,jmin:jmax,kmin:kmax)

    !--------------------------------------
    ! initialize: Divergence of velocity
    !
    xyz_VelDivNs =                 &
      &   xyz_dx_pyz( pyz_VelXNs ) &
      & + xyz_dy_xqz( xqz_VelYNs ) &
      & + xyz_dz_xyr( xyr_VelZNs )
!    call HistoryAutoPut(TimeN, 'VelDiv', xyz_VelDivNs(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'VelDiv', xyr_VelZNs(1:nx,1:ny,1:nz))
    
    !--------------------------------------
    ! VelX
    !
    pyz_DVelXDtSWF   =   pyz_dx_xyz( AlphaH * xyz_VelDivNs ) 
!    pyz_DVelXDtPGrad = - pyz_xyz( CpDry * xyz_VPTempBZ ) * pyz_dx_xyz( xyz_ExnerNs )      
    pyz_DVelXDtPGrad = - CpDry * pyz_VPTempBZ * pyz_dx_xyz( xyz_ExnerNs )
    pyz_DVelXDtNs =   pyz_DVelXDtPGrad + pyz_DVelXDtSWF

    ! Time integration
    !
    pyz_VelXAs = pyz_VelXNs + DelTimeShort * ( pyz_DVelXDtNs + pyz_DVelXDtNl )

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelXDtPGrad', pyz_DVelXDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelXDtSWF',   pyz_DVelXDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_pyz( pyz_VelXAs ) ! (inout)

    !--------------------------------------
    ! VelY
    !    
    xqz_DVelYDtSWF   =   xqz_dy_xyz( AlphaH * xyz_VelDivNs ) 
!    xqz_DVelYDtPGrad = - xqz_xyz( CpDry * xyz_VPTempBZ ) * xqz_dy_xyz( xyz_ExnerNs )
    xqz_DVelYDtPGrad = - CpDry * xqz_VPTempBZ * xqz_dy_xyz( xyz_ExnerNs )
    xqz_DVelYDtNs =   xqz_DVelYDtPGrad + xqz_DVelYDtSWF
    
    ! Time integration
    !
    xqz_VelYAs = xqz_VelYNs + DelTimeShort * (xqz_DVelYDtNs + xqz_DVelYDtNl)

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelYDtPGrad', xqz_DVelYDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSWF',   xqz_DVelYDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_xqz( xqz_VelYAs ) ! (inout)

    !----------------------------------------
    ! 次の時刻の Exner 関数
    !

    ! 短い時間ステップで評価する圧力の式の tendency
    !
    xyz_DExnerDtExpnd = GasRDry * xyz_ExnerNs * xyz_VelDivNs / CvDry   &
      &              * FactorDExnerDtExpnd

    ! tendency の合計
    !
    xyz_DExnerDtNs = xyz_DExnerDtNs + xyz_DExnerDtExpnd

    xyz_ExnerAs = xyz_Exner( &
      & pyz_VelXAs,          &
      & xqz_VelYAs,          &
      & xyr_VelZNs,          &
      & xyz_VelDivNs,        &
      & xyz_ExnerNs,         &
      & xyr_DVelZDtNl,       &
      & xyz_DExnerDtNl,      &
      & xyz_DExnerDtNs       &
      & )

    ! Set Margin
    !
    call SetMargin_xyz( xyz_ExnerAs ) ! (inout)

    !--------------------------------------
    ! VelZ
    !
    xyr_DVelZDtSWF =  xyr_dz_xyz( AlphaV * xyz_VelDivNs ) 
    xyr_DVelZDtPGrad =                                          &
      & - CpDry * xyr_VPTempBZ                               &
      &   * (                                                &
      &         beta           * xyr_dz_xyz( xyz_ExnerAs )   &
      &       + (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerNs )   &
      &     )
    xyr_DVelZDtNs = xyr_DVelZDtPGrad + xyr_DVelZDtSWF

    ! Time integration
    !
    xyr_VelZAs = xyr_VelZNs + DelTimeShort * (xyr_DVelZDtNs + xyr_DVelZDtNl)

    ! output
    !
    call HistoryAutoPut(TimeN, 'DVelZDtPGrad', xyr_DVelZDtPGrad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSWF',   xyr_DVelZDtSWF(1:nx,1:ny,1:nz))

    ! Set Margin
    !
    call SetMargin_xyr( xyr_VelZAs ) ! (inout)

  end subroutine Dynamics_Short_forcing

  
!!!--------------------------------------------------------------------!!!
  subroutine DynamicsVI_init()
    !
    !エクスナー関数を陰解法で解く際に必要となる, 係数行列の要素を決め, 
    !LU 分解を行う. 
    !

    ! モジュール読み込み
    use dc_types,   only : DP
!    use dc_message, only : MessageNotify
    use gridset,    only : imin, imax,      &!
      &                    jmin, jmax,      &
      &                    kmin, kmax,      &
      &                    nx,              &! x 方向の物理領域の上限
      &                    ny,              &! x 方向の物理領域の上限
      &                    nz                ! y 方向の物理領域の上限
    use constants,  only : CpDry           ! 乾燥成分の比熱
    use timeset,    only : DelTimeShort
    use axesset,    only : dz        ! 格子間隔
    use basicset,   only : xyz_VelSoundBZ,  &!基本場の音速 
      &                    xyz_DensBZ,      &!基本場の密度
      &                    xyz_VPTempBZ      !基本場の温位
    use average,    only : xyr_xyz

    !暗黙の型宣言禁止
    implicit none

    real(DP)  :: DTS ! 短い時間格子

    DTS = DelTimeShort

    !配列の割り付け
    allocate( A(1:nz) )
    allocate( B(1+1:nz) )
    allocate( C(1:nz-1) )
    allocate( xyz_F1BZ(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_F2BZ(imin:imax,jmin:jmax,kmin:kmax) )

    !----------------------------------------------------------------
    ! 係数行列と共通して利用される配列の値を決める
    !----------------------------------------------------------------

    !係数行列の計算
    !  A, B, C を求める際, F1BZ と F2BZ は X 方向に一様なので. 
    !  nx, ny の値で代表させることとした. 
    xyz_F1BZ =                                                &
      &  ( xyz_VelSoundBZ ** 2.0d0 )                          &
      &   / (CpDry * xyz_DensBZ * (xyz_VPTempBZ ** 2.0d0))

    xyr_F2BZ =                                                &
      &  xyr_xyz(                                             &
      &    CpDry * xyz_DensBZ * ( xyz_VPTempBZ ** 2.0d0 )     &
      &   )
        
    A(1+1: nz-1) =                                &
      & (beta ** 2.0d0)                           &
      &    * xyz_F1BZ(nx,ny,1+1: nz-1)            &  
      &    * (DTS ** 2.0d0)                       &
      &    * (                                    &
      &          xyr_F2BZ(nx,ny,1+1: nz-1)        &
      &            / dz                           &
      &        + xyr_F2BZ(nx,ny,1  : nz-2)        &
      &            / dz                           &
      &       )                                   &
      &    / dz                                   &
      & + 1.0d0

    A(1) =                                        &
      & (beta ** 2.0d0)                           &
      &   * xyz_F1BZ(nx,ny,1)                     &
      &   * xyr_F2BZ(nx,ny,1)                     &
      &     / dz                                  &
      &   * (DTS ** 2.0d0)                        &
      &   / dz                                    &
      & + 1.0d0                                         

    A(nz) =                                       &
      & (beta ** 2.0d0)                           &
      &   * xyz_F1BZ(nx,ny,nz)                    &
      &   * xyr_F2BZ(nx,ny,nz-1)                  &
      &     / dz                                  &
      &   * (DTS ** 2.0d0)                        &
      &   / dz                                    &
      & + 1.0d0                                         
    
    B(1+1:nz) =                                   &
      & - (beta ** 2.0d0)                         &
      &   * xyz_F1BZ(nx,ny,1:nz-1)                & 
      &   * xyr_F2BZ(nx,ny,1:nz-1)                &
      &   * (DTS ** 2.0d0)                        &
      &   / ( dz * dz ) 
    
    C(1: nz-1) =                                  &
      & - ( beta ** 2.0d0 )                       &
      &   * xyz_F1BZ(nx,ny,1+1:nz)                &
      &   * xyr_F2BZ(nx,ny,1  :nz-1)              & 
      &   * (DTS ** 2.0d0) &
      &   / ( dz * dz )


    !----------------------------------------------------------------
    ! 係数行列を LU 分解
    !----------------------------------------------------------------
    !配列の大きさを保管
    N   = nz  !係数行列/改行列の次数, 整合寸法
    M   = nx * ny 
                               !方程式の組数
    NUD = 1                    !係数行列の上三角部分の帯幅
    NLD = 1                    !係数行列の下三角部分の帯幅
    NAL = NLD                  !LU 分解の結果 L の整合寸法
    NA  = NUD + NLD + 1

    !配列の割り当て
!    allocate( AL1(N), AL2(NAL, N), AU2(NA, N), IP(N) )
    allocate( AL1(N), IP(N) )

    !LU 分解の実行
    !  LAPACK の利用
    call ResolvLU_Lapack( )

   
  end subroutine DynamicsVI_init
  

!!!--------------------------------------------------------------------!!!
  function xyz_Exner(      &
    & pyz_VelXAs,          &
    & xqz_VelYAs,          &
    & xyr_VelZNs,          &
    & xyz_VelDivNs,        &
    & xyz_ExnerNs,         &
    & xyr_DVelZDtNl,       &
    & xyz_DExnerDtNl,      &
    & xyz_DExnerDtNs       &
    & )
    !
    !陰解法を用いたエクスナー関数の計算. 
    !
    
    ! モジュールの読み込み
    !
    use dc_types, only : DP
    use gridset,  only : imin,            &! x 方向の配列の下限
      &                  imax,            &! x 方向の配列の上限
      &                  jmin,            &! y 方向の配列の下限
      &                  jmax,            &! y 方向の配列の上限
      &                  kmin,            &! z 方向の配列の下限
      &                  kmax,            &! z 方向の配列の上限
      &                  nx, ny, nz        ! 物理領域の大きさ
    use constants,only : CpDry             ! 乾燥成分の比熱
    use timeset,  only : DelTimeShort, TimeN
    use basicset, only : xyz_VelSoundBZ,  &!基本場の音速
      &                  xyz_DensBZ,      &!基本場の密度
      &                  xyz_VPTempBZ,    &!基本場の仮温位
      &                  xyr_VPTempBZ      !基本場の仮温位
    use axesset,  only : dz
    use average,  only : xyr_xyz
    use differentiate_center2, &
      &           only : xyr_dz_xyz, xyz_dz_xyr, xyz_dx_pyz, xyz_dy_xqz
    use gtool_historyauto,                    &
      &              only : HistoryAutoPut 

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数
    real(DP), intent(in)   :: pyz_VelXAs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 u [τ+Δτ]
    real(DP), intent(in)   :: xqz_VelYAs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 v [τ+Δτ]
    real(DP), intent(in)   :: xyr_VelZNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !速度 w [τ]
    real(DP), intent(in)   :: xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)   :: xyr_DVelZDtNl &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z 方向の外力項[t]
    real(DP), intent(in)   :: xyz_DExnerDtNl &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z 方向の外力項[t]
    real(DP), intent(in)   :: xyz_DExnerDtNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !Z 方向の外力項[t]
    real(DP), intent(in)   :: xyz_ExnerNs &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !無次元圧力
    real(DP)               :: xyz_Exner &
      &                     (imin:imax,jmin:jmax,kmin:kmax) 
                                                           !無次元圧力[τ+Δτ]

    !変数定義
    real(DP)  :: D1(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: D2(1:nx,1:ny,1:nz)  
    real(DP)  :: D(nx*ny,nz)
    real(DP)  :: E(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: F(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)  :: DTS ! 短い時間格子間隔
    integer   :: ix, jy, kz

    ! Initialize
    !
    DTS = DelTimeShort
    xyz_Exner = 0.0d0

    !行列計算のための係数
    E =   &
      & - ( 1.0d0 - beta ) * xyr_dz_xyz( xyz_ExnerNs )              &
      & + ( AlphaV * xyr_dz_xyz( xyz_VelDivNs ) + xyr_DVelZDtNl )   &
      &    / ( CpDry * xyr_VPTempBZ ) 

     F = - beta * xyz_F1BZ * DTS                                     &
      & * xyz_dz_xyr(                                               &
      &    xyr_xyz( xyz_DensBZ * xyz_VPTempBZ)                      &
      &    * (                                                      &
      &         xyr_VelZNs                                          &
      &       - xyr_xyz(CpDry * xyz_VPTempBZ) * DTS                 &
      &         * (1.0d0 - beta) * xyr_dz_xyz( xyz_ExnerNs )        &
      &       + xyr_dz_xyz( AlphaV * xyz_VelDivNs ) * DTS           &
      &       + xyr_DVelZDtNl * DTS                                 &
      &      )                                                      &
      &   )                                                         &
      & + (xyz_DExnerDtNs + xyz_DExnerDtNl ) * DTS

    D1 = xyz_ExnerNs                                                &
      & - (1.0d0 - beta)                                            &
      &   * xyz_F1BZ * DTS                                          &
      &   * xyz_dz_xyr(                                             &
      &       xyr_xyz(xyz_DensBZ * xyz_VPTempBZ) * xyr_VelZNs       &
      &     )                                                       &
      & - (xyz_VelSoundBZ ** 2.0d0) * DTS                           &
      &   / ( CpDry * xyz_VPTempBZ )                                &
      &   * ( xyz_dx_pyz( pyz_VelXAs ) +  xyz_dy_xqz( xqz_VelYAs ) )&
      & + F

    D1(:,:,1) = D1(:,:,1)                           &
      & - beta * xyz_F1BZ(:,:,1) * (DTS ** 2.0d0)   &
      &   * xyr_F2BZ(:,:,1-1) * E(:,:,1-1)          &
      &   / dz
    
    D1(:,:,nz) = D1(:,:,nz)                         &
      & + beta * xyz_F1BZ(:,:,nz) * (DTS ** 2.0d0)  &
      &   * xyr_F2BZ(:,:,nz) * E(:,:,nz)            &
      &   / dz

!    call HistoryAutoPut(TimeN, 'D', D1(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'E', E(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'F', F(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENs', xyz_DExnerDtNs(1:nx,1:ny,1:nz))
!    call HistoryAutoPut(TimeN, 'ENl', xyz_DExnerDtNl(1:nx,1:ny,1:nz))

    D2 = D1(1:nx,1:ny,1:nz)

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          D(ix + nx * (jy - 1), kz) =  D2(ix,jy,kz)
        end do
      end do
    end do

    !-----------------------------------------------------------
    !連立一次方程式の解を求める
    !------------------------------------------------------------

    !解の計算
    !  LAPACK 利用
    call LinSolv_Lapack( D )

    !戻り値を出力
    do kz = 1, nz
      do jy = 1, ny 
        do ix = 1, nx 
          xyz_Exner(ix,jy,kz) = D(ix + nx * (jy - 1 ), kz)
        end do
      end do
    end do

  end function xyz_Exner

!!!--------------------------------------------------------------------!!!
  subroutine ResolvLU_Lapack(  )
    !
    !実 3 項行列の LU 分解(倍精度). LAPACK 利用
    !

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer    :: INFO  !解のコンディションチェック
    
    !変数の初期化
    INFO = 0
    
    !解行列の計算. LAPACK を使用. 
    call DGTTRF(N, C, A, B, AL1, IP, INFO)
    
    !解のコンディションをチェック. 
!    if (INFO /= 0) then
!      call MessgaeNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if
    
  end subroutine ResolvLU_Lapack
  

!!!--------------------------------------------------------------------!!!
  subroutine LinSolv_Lapack( X )
    !
    !LU 分解された実 3 項行列の連立 1 次方程式(倍精度). LAPACK 利用
    !

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(inout) :: X(M, N)     !定数/解行列
    real(DP)                :: TX(N, M)    !解行列を転置したもの
    integer                :: NRHS         !
    integer                :: INFO
    integer                :: LDB
    character(1),parameter :: TRANS = 'N'

    !変数の初期化
    NRHS = M
    INFO = 0
    LDB  = N
    TX = transpose( X )
    
    !解行列の計算. LAPACK を使用. 
    call DGTTRS(TRANS, N, NRHS, C, A, B, AL1, IP, TX, LDB, INFO)

    !解の出力
    X = transpose( TX )
     
    !解のコンディションをチェック. 
!    if (INFO /= 0) then
!      call MessageNotify("Error", "lapack_linear", "INFO is not 0")
!      stop
!    end if
     
  end subroutine LinSolv_Lapack



  subroutine  Dynamics_Tendency_output
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
    integer :: l

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

    do l = 1, ncmax
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtAdv', &
        & dims=(/'x','y','z','t'/),     &
        & longname='Advection term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')
      
      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtDiff', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Diffusion term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
        & units='kg.kg-1.s-1',    &
        & xtype='double')

      call HistoryAutoAddVariable(  &
        & varname='D'//trim(SpcWetSymbol(l))//'DtFall', & 
        & dims=(/'x','y','z','t'/),     &
        & longname='Fall term of '          &
        &           //trim(SpcWetSymbol(l))//' mixing ratio',  &
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

!    call HistoryAutoAddVariable(  &
!      & varname='VelXTndNs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='Velocity Tendency (x)',  &
!      & units='m.s-2',    &
!      & xtype='double')

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

!    call HistoryAutoAddVariable(  &
!      & varname='VelYTndNs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='Velocity Tendency (y)',  &
!      & units='m.s-2',    &
!      & xtype='double')

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

!    call HistoryAutoAddVariable(  &
!      & varname='D', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix D',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='E', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix E',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='F', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='matrix F',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENs', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNs',  &
!      & units='s-1',    &
!      & xtype='double')

!    call HistoryAutoAddVariable(  &
!      & varname='ENl', &
!      & dims=(/'x','y','z','t'/),     &
!      & longname='DExnerDtNl',  &
!      & units='s-1',    &
!      & xtype='double')

  end subroutine Dynamics_Tendency_output

end module DynamicsHEVI

