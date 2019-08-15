!= deepconv/arare 湿潤大気対流計算用主プログラム 
!
!= deepconv/arare main program for moist atmospheric convection 
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro), ODAKA Masatsugu
! Version::   $Id: arare.f90,v 1.49 2014/07/08 01:01:44 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

program deepconv_arare
  !
  ! 非静力学モデル deepconv/arare 湿潤大気対流計算用主プログラム (三次元版)
  !

  ! モジュール引用  use statement 
  !

  ! gtool5 関連 
  ! gtool5 modules
  !
  use dc_types,          only: STRING, DP
  use dc_message,        only: MessageNotify, MessageSuppressMPI
  use gtool_history,     only: HistoryPut
  use gtool_historyauto, only: HistoryAutoPut

  ! 初期設定モジュール
  ! Initialize module
  !
  use mpi_wrapper,     only : MPIWrapperInit,                          &
    &                         MPIWrapperFinalize
  use argset,          only : argset_init
  use clockset,        only : ClocksetInit,                            &
    &                         ClocksetClose,    ClocksetPredict,       &
    &                         ClockSetPreStart, ClockSetLoopStart,     &
    &                         ClockSetPreStop,  ClockSetLoopStop
  use gridset,         only : gridset_init,                            &
    &                         imin, imax, jmin, jmax, kmin, kmax,      &
    &                         nx, ny, nz, ncmax
  use timeset,         only : timeset_init,                            &
    &                         TimesetDelTimeHalf, TimesetProgress,     &
    &                         TimeA, TimeN, TimeB,                     &
    &                         Nstep, NstepShort, NstepLong,            &
    &                         NstepOutput, FlagInitialRun
!    &                         TimeA, TimeN, DelTimeLong,               &
!    &                         NstepShort, NstepOutput,                 &
!    &                         EndTime, FlagInitialRun
  use axesset,         only : axesset_init
  use average,         only : xyz_pyz, xyz_xyr, xyz_xqz
  use constants,       only : constants_init
  use composition,     only : composition_init,                        &
    &                         SpcWetSymbol, GasNum
  use fileset,         only : fileset_init
  use basicset,        only : xyzf_QMixBZ, xyz_PTempBZ, xyz_ExnerBZ, &
    &                         xyz_VelSoundBZ
  use basicset,        only : xyzf_QMixBZ, xyz_DensBZ, xyz_EffMolWtBZ, &
    &                         xyz_PTempBZ, xyz_TempBZ, xyz_PressBZ,    &
    &                         xyz_VelSoundBZ, xyz_ExnerBZ, xyz_VPTempBZ
!!use ChemCalc,      only : ChemCalc_init, ChemCalc_init2
  use ChemCalc,        only : ChemCalc_init
  use namelist_util,   only : NmlutilInit,                             &
    &                         namelist_filename

  ! 化学量計算モジュール 
  ! Chemical calculation modules
  !
  use chemdata,        only : chemdata_init

  ! 力学過程計算用関数モジュール
  ! Dynamical processes module
  !
!  use DynamicsHEVI, only : Dynamics_Init,                        &
!  use DynamicsHEVI_v2, only : Dynamics_Init,                        &
  use DynamicsHEVI_v3, only : Dynamics_Init,                        &
    &                         Dynamics_Long_forcing,                &
    &                         Dynamics_Short_forcing

  ! 乱流拡散計算用モジュール
  ! Turbulent diffusion module
  !
  use Turbulence_kw1978_v2, only : Turbulence_kw1978_Init,          &
!  use Turbulence_kw1978_v1, only : Turbulence_kw1978_Init,          &
    &                              Turbulence_KW1978_forcing,       &
    &                              KmMax
  use Turbulence_constKm, only : Turbulence_constKm_Init,           & 
    &                            Turbulence_constKm_forcing

  ! 境界からのフラックス計算用モジュール
  ! Surface flux module
  !
  use Surfaceflux_Diff,       only : Surfaceflux_Diff_init,         &
    &                                Surfaceflux_Diff_forcing
!!use Surfaceflux_Bulk_l1982, only : Surfaceflux_Bulk_init,         &
!!  &                                Surfaceflux_Bulk_forcing
  use Surfaceflux_Bulk,       only : Surfaceflux_Bulk_init,         &
    &                                Surfaceflux_Bulk_forcing
  use Surfaceflux_Const,      only : Surfaceflux_Const_init,        &
    &                                Surfaceflux_Const_forcing
  use Surfaceflux_baker1998,  only : Surfaceflux_baker1998_init,    &
    &                                Surfaceflux_baker1998_forcing

  ! 放射強制計算用モジュール
  ! Radiative forceing module
  !
  use Radiation_simple,      only : Radiation_simple_init,           &
    &                               Radiation_HeatConst_forcing,     &
    &                               Radiation_HeatVary_forcing
  use Radiation_sounding,    only : Radiation_sounding_init,         &
    &                               Radiation_sounding_forcing
  use Radiation_baker1998,   only : Radiation_baker1998_init,        &
    &                               Radiation_baker1998_forcing
  use Radiation_heatbalance, only : Radiation_heatbalance_init,      &
    &                               Radiation_heatbalance_forcing

  ! 湿潤過程計算用モジュール
  ! Moist processes modules
  !
  use Cloudphys_K1969,       only : Cloudphys_K1969_Init,            &
    &                               Cloudphys_K1969_forcing 
  use Cloudphys_marscond,    only : cloudphys_marscond_Init,         &
    &                               cloudphys_marscond_forcing

  ! 下請けモジュール
  ! Utility modules
  !
  use damping,               only : Damping_init,                    &
    &                               SpongeLayer_forcing,             &
    &                               SpongeLayer_MeanFlow
  use fillnegative,          only : FillNegative_init
  use energymonit,           only : EnergyMonit_init,                &
    &                               EnergyMonit_exec
  use cflcheck,              only : CFLCheckTimeShort,               &
    &                               CFLCheckTimeLongVelX,            &
    &                               CFLCheckTimeLongVelY,            &
    &                               CFLCheckTimeLongVelZ
  use setmargin,             only : SetMargin_init

  ! ファイル入出力モジュール
  ! File I/O module
  !
  use RestartFileIO,         only : ReStartFileio_init,              &
    &                               ReStartFileio_Finalize,          &
    &                               ReStartFileio_BasicZ_Get,        &
    &                               ReStartFileio_Var_Get,           &
    &                               rstat
  use Arare4RestartFileIO,   only : Arare4ReStartFileio_BasicZ_Get,  &
    &                               Arare4ReStartFileio_Var_Get 
  use HistoryFileIO,         only : HistoryFileio_init,              &
    &                               HistoryFileio_Finalize
  use BasicFileIO,           only : BasicFileio_Output


  ! 暗黙の型宣言禁止
  !
  implicit none

  ! 内部変数
  ! Internal variables
  !
  character(*), parameter:: prog_name = 'arare'
                            ! 主プログラム名.
                            ! Main program name

  real(DP), allocatable :: pyz_VelXBl(:,:,:)    
                             ! $ u (t-\Delta t) $ 東西風 ; zonal wind
  real(DP), allocatable :: pyz_VelXNl(:,:,:)    
                             ! $ u (t) $          東西風 ; zonal wind
  real(DP), allocatable :: xyz_VelXNl(:,:,:)    
                             ! $ u (t) $          東西風 ; zonal wind
  real(DP), allocatable :: pyz_VelXAl(:,:,:)    
                             ! $ u (t+\Delta t) $ 東西風 ; zonal wind
  real(DP), allocatable :: pyz_VelXNs(:,:,:)    
                             ! $ u (\tau) $ 東西風 ; zonal wind
  real(DP), allocatable :: pyz_VelXAs(:,:,:)    
                             ! $ u (\tau +\Delta \tau) $ 東西風 ; zonal wind
  real(DP), allocatable :: xqz_VelYBl(:,:,:)    
                             ! $ v (t-\Delta t) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xqz_VelYNl(:,:,:)    
                             ! $ v (t) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xyz_VelYNl(:,:,:)    
                             ! $ v (t) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xqz_VelYAl(:,:,:)    
                             ! $ v (t+\Delta t) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xqz_VelYNs(:,:,:)   
                             ! $ v (\tau -\tau) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xqz_VelYAs(:,:,:)
                             ! $ v (\tau) $ 南北風 ; meridonal wind
  real(DP), allocatable :: xyr_VelZBl(:,:,:)    
                             ! $ w (t-\Delta t) $ 鉛直風 ; vertical wind
  real(DP), allocatable :: xyr_VelZNl(:,:,:)    
                             ! $ w (t) $ 鉛直風 ; vertical wind
  real(DP), allocatable :: xyz_VelZNl(:,:,:)    
                             ! $ w (t) $ 鉛直風 ; vertical wind
  real(DP), allocatable :: xyr_VelZAl(:,:,:)    
                             ! $ w (t+\Delta t) $ 鉛直風 ; vertical wind
  real(DP), allocatable :: xyr_VelZNs(:,:,:)    
                             ! $ w (\tau) $ 鉛直風 ; vertical wind
  real(DP), allocatable :: xyr_VelZAs(:,:,:) 
                             ! $ w (\tau +\Delta \tau)  鉛直風 ; vertical wind
  real(DP), allocatable :: xyz_ExnerBl(:,:,:)   
                             ! $ \pi (t-\Delta t) $ 圧力関数 ; Exner function
  real(DP), allocatable :: xyz_ExnerNl(:,:,:)   
                             ! $ \pi (t) $ 圧力関数 ; Exner function
  real(DP), allocatable :: xyz_ExnerAl(:,:,:)
                             ! $ \pi (t+\Delta t) $ 圧力関数 ; Exner function
  real(DP), allocatable :: xyz_ExnerNs(:,:,:)   
                             ! $ \pi (\tau -\Delta \tau) $ 圧力関数 ; Exner function
  real(DP), allocatable :: xyz_ExnerAs(:,:,:)   
                             ! $ \pi (\tau) $ 圧力関数 ; Exner function
  real(DP), allocatable :: xyz_PTempBl(:,:,:) 
                             ! $ \theta (t-\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_PTempNl(:,:,:) 
                             ! $ \theta (t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_PTempAl(:,:,:) 
                             ! $ \theta (t+\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_PTempNs(:,:,:) 
                             ! $ \theta (t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_PTempAs(:,:,:) 
                             ! $ \theta (t+\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_CDensBl(:,:,:) 
                             ! $ \theta (t-\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_CDensNl(:,:,:) 
                             ! $ \theta (t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_CDensAl(:,:,:) 
                             ! $ \theta (t+\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_CDensNs(:,:,:) 
                             ! $ \theta (t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_CDensAs(:,:,:) 
                             ! $ \theta (t+\Delta t) $ 温位 ; Potential temp.
  real(DP), allocatable :: xyz_KmBl(:,:,:)
                             ! $ Km (t-\Delta t) $ 乱流拡散係数 
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyz_KmNl(:,:,:)
                             ! $ K_m (t) $ 乱流拡散係数 
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyz_KmAl(:,:,:)
                             ! $ K_m (t+\Delta t) $ 乱流拡散係数 
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyz_KhBl(:,:,:)      
                             ! $ K_h (t-\Delta t) $ 乱流拡散係数
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyz_KhNl(:,:,:)
                             ! $ K_h (t) $ 乱流拡散係数 
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyz_KhAl(:,:,:)
                             ! $ K_h (t+\Delta t) $ 乱流拡散係数
                             ! Turbulent diffusion coeff. 
  real(DP), allocatable :: xyzf_QMixBl(:,:,:,:) 
                             ! $ q (t-\Delta t) $ 湿潤量の混合比
                             ! Mixing ratio of moist variables.
  real(DP), allocatable :: xyzf_QMixNl(:,:,:,:) 
                             ! $ q (t) $ 湿潤量の混合比 
                             ! Mixing ratio of moist variables
  real(DP), allocatable :: xyzf_QMixAl(:,:,:,:) ! 
                             ! $ q (t+\Delta t) $ 湿潤量の混合比 
                             !Mixing ratio of moist variables

  real(DP), allocatable :: pyz_DVelXDtNl(:,:,:)
  real(DP), allocatable :: xqz_DVelYDtNl(:,:,:)
  real(DP), allocatable :: xyr_DVelZDtNl(:,:,:)
  real(DP), allocatable :: xyz_DPTempDtNl(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNl(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNs(:,:,:)
  real(DP), allocatable :: xyzf_DQMixDtNl(:,:,:,:)
  real(DP), allocatable :: xyz_DKmDtNl(:,:,:)
  real(DP), allocatable :: xyz_DCDensDtNl(:,:,:)

  integer :: s, t, tau  ! do ループ変数 ; do loop variable 

  integer            :: IDTurbMethod          = 0
  integer, parameter :: IDTurbKW1978          = 2
  integer, parameter :: IDTurbConstKm         = 3
  integer            :: IDRadMethod           = 0
  integer, parameter :: IDRadHeatConst        = 1
  integer, parameter :: IDRadHeatVary         = 2
  integer, parameter :: IDRadHeatBalance      = 3
  integer, parameter :: IDRadSounding         = 4
  integer, parameter :: IDRadBaker1998        = 5
  integer            :: IDSurfaceMethod       = 0
  integer, parameter :: IDSurfaceDiff         = 1
  integer, parameter :: IDSurfaceBulk         = 2
  integer, parameter :: IDSurfaceConst        = 3
  integer, parameter :: IDSurfaceBaker1998    = 4
  integer            :: IDCloudMethod         = 0
  integer, parameter :: IDCloudK1969          = 1
  integer, parameter :: IDCloudMarsCond       = 2
  integer            :: IDDebugMethod         = 0
  integer, parameter :: IDDebugNoTendencyLong = 1
  integer, parameter :: IDDebugWindConst      = 2
  integer            :: IDRestartMethod       = 0
  integer, parameter :: IDRestartArare6       = 1
  integer, parameter :: IDRestartArare4       = 2
  integer            :: IDDampingMethod       = 0
  integer, parameter :: IDDampingSpongeLayer         = 1
  integer, parameter :: IDDampingSpongeLayerMeanFlow = 2


  !------------------------------------------
  ! 初期化手続き ; Initialize procedure 
  !
  call MainInit

  !------------------------------------------
  ! 時間積分 time integration 
  !
  call MessageNotify( "M", "main", "Time Integration Start" )

  ! CPU 時間計測開始
  ! Start CPU time counting 
  !
  call ClocksetLoopStart
  
  ! 時間発展ループのスタート
  !
  do while ( Nstep <= NstepLong )

    ! clear tendency 
    !
    pyz_DVelXDtNl = 0.0d0 ; xqz_DVelYDtNl  = 0.0d0 ; xyr_DVelZDtNl  = 0.0d0
    xyz_DPTempDtNl= 0.0d0 ; xyz_DExnerDtNl = 0.0d0 ; xyzf_DQMixDtNl = 0.0d0
    xyz_DKmDtNl   = 0.0d0 ; xyz_DCDensDtNl = 0.0d0

    !------------------------------------------
    ! スポンジ層; sponge layer
    !
    select case ( IDDampingMethod )
       
    case ( IDDampingSpongeLayer )
       
       call SpongeLayer_forcing(                           &
         &   pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,         & !(in)
         &   xyz_PTempBl, xyz_ExnerBl,                     & !(in)
         &   pyz_DVelXDtNl,  xqz_DVelYDtNl, xyr_DVelZDtNl, & !(inout)
         &   xyz_DPTempDtNl, xyz_DExnerDtNl                & !(inout)
         & )

    case ( IDDampingSpongeLayerMeanFlow )

       call SpongeLayer_forcing(                           &
         &   pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,         & !(in)
         &   xyz_PTempBl, xyz_ExnerBl,                     & !(in)
         &   pyz_DVelXDtNl,  xqz_DVelYDtNl, xyr_DVelZDtNl, & !(inout)
         &   xyz_DPTempDtNl, xyz_DExnerDtNl                & !(inout)
         & )
    
       call SpongeLayer_MeanFlow(                       &
         &   pyz_VelXBl,  xqz_VelYBl,                   & !(in)
         &   pyz_DVelXDtNl,  xqz_DVelYDtNl              & !(inout)
         & )

    end select

    !-----------------------------------------
    ! 移流拡散.
    ! Advection and diffusion
    !
    call Dynamics_Long_forcing(       &
      &   pyz_VelXBl,   pyz_VelXNl,   & ! (in)
      &   xqz_VelYBl,   xqz_VelYNl,   & ! (in)
      &   xyr_VelZBl,   xyr_VelZNl,   & ! (in)
      &   xyz_PTempBl,  xyz_PTempNl,  & ! (in) 
      &   xyz_ExnerBl,  xyz_ExnerNl,  & ! (in)
      &   xyzf_QMixBl,  xyzf_QMixNl,  & ! (in)
      &   xyz_KmBl,     xyz_KmNl,     & ! (in) 
      &   pyz_DVelXDtNl,              & ! (inout)
      &   xqz_DVelYDtNl,              & ! (inout)
      &   xyr_DVelZDtNl,              & ! (inout)
      &   xyz_DPTempDtNl,             & ! (inout)
      &   xyz_DExnerDtNl,             & ! (inout)
      &   xyzf_DQMixDtNl,             & ! (inout)
      &   xyz_DKmDtNl                 & ! (inout)
      )

    !-------------------------------
    ! 物理過程: 乱流
    !
    select case ( IDTurbMethod )
       
    case ( IDTurbKW1978 )

      ! Klemp and Wilhelmson (1978)
      !
      call turbulence_KW1978_forcing(                      &
        &   pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,          &!(in)
        &   xyz_PTempBl, xyz_ExnerBl, xyzf_QMixBl,         &!(in)
        &   xyz_KmBl,    xyz_KhBl,    xyz_CDensBl,         &!(in)
        &   pyz_DVelXDtNl, xqz_DVelYDtNl,  xyr_DVelZDtNl,  &!(inout)
        &   xyz_DPTempDtNl,xyz_DExnerDtNl, xyzf_DQMixDtNl, &!(inout)
        &   xyz_DKmDtNl,   xyz_DCDensDtNl                  &!(inout)
        & )
      
      ! 拡散係数の積分
      !
      call turbulence_integrate
      
    case ( IDTurbConstKm )
    
      ! 拡散係数一定
      !
      call turbulence_ConstKm_forcing(                     &
        &   pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,          &!(in)
        &   xyz_PTempBl,                                   &!(in)
        &   pyz_DVelXDtNl,  xqz_DVelYDtNl, xyr_DVelZDtNl,  &!(inout)
        &   xyz_DPTempDtNl, xyz_DExnerDtNl,                &!(inout)
        &   xyz_KmAl, xyz_KhAl                             &!(out)
        & )
      
    end select
     
    !-------------------------------
    ! 物理過程: 放射
    !
    select case (IDRadMethod)
      
    case (IDRadHeatConst)

      call Radiation_HeatConst_forcing(    &
        & xyz_DPTempDtNl, xyz_DExnerDtNl   & !(inout)
        & )
            
    case (IDRadHeatVary)

      call Radiation_HeatVary_forcing(     &
        & xyz_DPTempDtNl, xyz_DExnerDtNl   & !(inout)
        & )
            
    case (IDRadSounding)

      call Radiation_Sounding_forcing(     &
        & xyz_DPTempDtNl, xyz_DExnerDtNl   & !(inout)
        & )

    case (IDRadBaker1998)

      call Radiation_baker1998_forcing(    &
        & xyz_DPTempDtNl, xyz_DExnerDtNl   & !(inout)
        & )

    case (IDRadHeatBalance)
      call radiation_heatbalance_forcing(  &
        &   xyz_ExnerNl, xyz_PTempNl,      & !(in)
        &   xyz_DPTempDtNl, xyz_DExnerDtNl & !(inout)
        & )

    end select

    !--------------------------------
    ! 境界からの熱・運動量輸送
    !
    select case (IDSurfaceMethod)

    case (IDSurfaceDiff)
      call Surfaceflux_Diff_forcing(             &
        &   xyz_PTempNl, xyzf_QMixNl,            & !(in)
        &   xyz_DPTempDtNl, xyz_DExnerDtNl,      & !(inout)
        &   xyzf_DQMixDtNl                       & !(inout)
        & )
     
    case (IDSurfaceBulk)
      call Surfaceflux_Bulk_forcing(             &
        &   pyz_VelXNl, xqz_VelYNl, xyz_PTempNl, &!(in)
        &   xyz_ExnerNl, xyzf_QMixNl,            &!(in)
        &   pyz_DVelXDtNl, xqz_DVelYDtNl,        &!(inout)
        &   xyz_DPTempDtNl, xyz_DExnerDtNl,      &!(inout)
        &   xyzf_DQMixDtNl                       &!(inout)
        & )

    case (IDSurfaceConst)
      call Surfaceflux_Const_forcing(            &
        &   pyz_DVelXDtNl, xqz_DVelYDtNl,        &!(inout)
        &   xyz_DPTempDtNl, xyz_DExnerDtNl,      &!(inout)
        &   xyzf_DQMixDtNl                       &!(inout)
        & )

    case (IDSurfaceBaker1998)
      call Surfaceflux_baker1998_forcing(        &
        &   xyz_DPTempDtNl, xyz_DExnerDtNl       &!(inout)
        & )
      
    end select
    
    !-----------------------------------------
    ! 温位と混合比の積分
    ! Integration
    !
    call PTemp_integrate
    call QMix_integrate
    
    !------------------------------------------
    ! 凝結過程. Al な値を入れ替え.
    ! 
    select case (IDCloudMethod)
    case (IDCloudK1969)
      call Cloudphys_K1969_forcing(   &
        &   xyz_ExnerNl,              &!(in)
        &   xyz_DExnerDtNl,           &!(inout)
        &   xyz_PTempAl,              &!(inout)
        &   xyzf_QMixAl               &!(inout)
        & )
    end select
    
    ! 短い時間ステップの初期値作成.
    ! Initial values set up for time integration with short time step.
    !
    pyz_VelXNs  = pyz_VelXBl
    xqz_VelYNs  = xqz_VelYBl
    xyr_VelZNs  = xyr_VelZBl
    xyz_ExnerNs = xyz_ExnerBl
    xyz_PTempNs = xyz_PTempBl
    xyz_CDensNs = xyz_CDensBl

    ! DEBUG: 長い時間ステップでの tendency を無視する場合
    !
    select case (IDDebugMethod)
    case (IDDebugNoTendencyLong)
      pyz_DVelXDtNl  = 0.0d0
      xqz_DVelYDtNl  = 0.0d0
      xyr_DVelZDtNl  = 0.0d0
      xyz_DPTempDtNl = 0.0d0
      xyz_DExnerDtNl = 0.0d0
      xyz_DExnerDtNs = 0.0d0
      xyzf_DQMixDtNl = 0.0d0
      xyz_DKmDtNl    = 0.0d0
      xyz_DCDensDtNl = 0.0d0
    end select
 
    ! 短い時間ステップの時間積分. オイラー法を利用.
    ! Time integration with short time step.
    !
    Euler: do tau = 1, NstepShort

      !
      ! clear tendency
      !
      xyz_DExnerDtNs = 0.0d0

      ! 火星計算の場合. 凝結量の評価はここで行う. 
      ! 
      select case (IDCloudMethod)
      case (IDCloudMarsCond)

        call cloudphys_marscond_forcing(  &
          &   xyz_PTempNs,         &  !(in) 温位
          &   xyz_ExnerNs,         &  !(in) エクスナー関数
          &   xyz_CDensNs,         &  !(in) 
          &   xyz_DPTempDtNl,      &  !(in)    
          &   xyz_DExnerDtNl,      &  !(in)    
          &   xyz_DCDensDtNl,      &  !(in)    
          &   xyz_PTempAs,         &  !(out) 温位
          &   xyz_CDensAs,         &  !(out) 雲密度
          &   xyz_DExnerDtNs       &  !(out) 
          & )

      end select

      ! HE-VI : 速度 u, v の計算.
      !
      call Dynamics_Short_forcing( &
        &   pyz_VelXNs,          & ! (in)
        &   xqz_VelYNs,          & ! (in)
        &   xyr_VelZNs,          & ! (in)
        &   xyz_ExnerNs,         & ! (in)
        &   pyz_DVelXDtNl,       & ! (in)
        &   xqz_DVelYDtNl,       & ! (in)
        &   xyr_DVelZDtNl,       & ! (in)
        &   xyz_DExnerDtNl,      & ! (in)
        &   xyz_DExnerDtNs,      & ! (in)
        &   pyz_VelXAs,          & ! (out)
        &   xqz_VelYAs,          & ! (out)
        &   xyr_VelZAs,          & ! (out)
        &   xyz_ExnerAs          & ! (out)
        & )

      ! 短い時間ステップのループを回すための処置
      ! Renew prognostic variables for next short time step integration.
      !
      xyz_ExnerNs  = xyz_ExnerAs
      pyz_VelXNs   = pyz_VelXAs
      xqz_VelYNs   = xqz_VelYAs
      xyr_VelZNs   = xyr_VelZAs
      xyz_PTempNs  = xyz_PTempAs
      xyz_CDensNs  = xyz_CDensAs

    end do Euler

    ! 最終的な短い時間ステップでの値を長い時間ステップでの値とみなす
    ! Renew prognostic variables for next long time step integration.
    !
    xyz_ExnerAl  = xyz_ExnerAs
    pyz_VelXAl   = pyz_VelXAs
    xqz_VelYAl   = xqz_VelYAs
    xyr_VelZAl   = xyr_VelZAs
    select case (IDCloudMethod)
    case (IDCloudMarsCond)
      xyz_PTempAl = xyz_PTempAs
      xyz_CDensAl = xyz_CDensAs
    end select

    ! 速度場を初期速度で固定する場合
    !
    select case (IDDebugMethod)
    case (IDDebugWindConst)
      pyz_VelXAl = pyz_VelXBl
      xqz_VelYAl = xqz_VelYBl
      xyr_VelZAl = xyr_VelZBl
    end select
   
    ! 時間フィルタ. 
    ! Time filter. 
    !
    call AsselinTimeFilter 
    
    ! ヒストリファイル出力. 時間フィルタをかけた後の数値を出力. 
    ! Out put to history file.
    !
    xyz_VelXNl = xyz_pyz(pyz_VelXNl)
    xyz_VelYNl = xyz_xqz(xqz_VelYNl)
    xyz_VelZNl = xyz_xyr(xyr_VelZNl)
    
    call HistoryAutoPut(TimeN, 'PTemp', xyz_PTempNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'PTempAll', xyz_PTempNl(1:nx, 1:ny, 1:nz) + xyz_PTempBZ(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'Exner', xyz_ExnerNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'ExnerAll', xyz_ExnerNl(1:nx, 1:ny, 1:nz) + xyz_ExnerBZ(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'VelX',  xyz_VelXNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'VelY',  xyz_VelYNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'VelZ',  xyz_VelZNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'Km',    xyz_KmNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'Kh',    xyz_KhNl(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'CDens', xyz_CDensNl(1:nx, 1:ny, 1:nz))
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s)), xyzf_QMixNl(1:nx, 1:ny, 1:nz, s))
    end do
    do s = 1, GasNum
      call HistoryAutoPut(TimeN, trim(SpcWetSymbol(s)), xyzf_QMixNl(1:nx, 1:ny, 1:nz, s))
    end do

    !------------------------------------------
    ! 保存量の出力
    ! Out put conservation variables
    !
    call EnergyMonit_exec(                    &
      &  pyz_VelXNl,  xqz_VelYNL, xyr_VelZNl, & !(in)
      &  xyz_ExnerNl, xyz_PTempNl             & !(in)
      & )
    
    !----------------------------------------------------
    ! リスタートファイルの作成
    ! Make restartfile.
    !    
    if (mod(Nstep, NstepOutput) == 0) then 

      ! 移流に対する CFL 条件のチェック 
      ! CFL condtion check for advection
      !
      call CFLCheckTimeLongVelX( pyz_VelXNl ) !(in)
      call CFLCheckTimeLongVelY( xqz_VelYNl ) !(in)
      call CFLCheckTimeLongVelZ( xyr_VelZNl ) !(in)
      
      ! リスタートファイルに変数を出力
      !
      call HistoryPut( 't',     TimeN,       rstat)
      call HistoryPut( 'VelX',  pyz_VelXNl,  rstat)
      call HistoryPut( 'VelY',  xqz_VelYNl,  rstat)
      call HistoryPut( 'VelZ',  xyr_VelZNl,  rstat)
      call HistoryPut( 'Exner', xyz_ExnerNl, rstat)
      call HistoryPut( 'PTemp', xyz_PTempNl, rstat)
      call HistoryPut( 'Km',    xyz_KmNl,    rstat)
      call HistoryPut( 'Kh',    xyz_KhNl,    rstat)
      call HistoryPut( 'CDens', xyz_CDensNl, rstat)
      call HistoryPut( 'QMix',  xyzf_QMixNl, rstat)    

      call HistoryPut( 't',     TimeA,       rstat)
      call HistoryPut( 'VelX',  pyz_VelXAl,  rstat)
      call HistoryPut( 'VelY',  xqz_VelYAl,  rstat)
      call HistoryPut( 'VelZ',  xyr_VelZAl,  rstat)
      call HistoryPut( 'Exner', xyz_ExnerAl, rstat)
      call HistoryPut( 'PTemp', xyz_PTempAl, rstat)
      call HistoryPut( 'Km',    xyz_KmAl,    rstat)
      call HistoryPut( 'Kh',    xyz_KhAl,    rstat)
      call HistoryPut( 'CDens', xyz_CDensAl, rstat)
      call HistoryPut( 'QMix',  xyzf_QMixAl, rstat) 
      
      ! 基本場のファイル出力
      !
      call HistoryPut( 'DensBZ',     xyz_DensBZ    , rstat)
      call HistoryPut( 'ExnerBZ',    xyz_ExnerBZ   , rstat)
      call HistoryPut( 'PTempBZ',    xyz_PTempBZ   , rstat)
      call HistoryPut( 'VPTempBZ',   xyz_VPTempBZ  , rstat)
      call HistoryPut( 'VelSoundBZ', xyz_VelSoundBZ, rstat)
      call HistoryPut( 'TempBZ',     xyz_TempBZ    , rstat)
      call HistoryPut( 'PressBZ',    xyz_PressBZ   , rstat)
      call HistoryPut( 'QMixBZ',     xyzf_QMixBZ   , rstat)
      call HistoryPut( 'EffMolWtBZ', xyz_EffMolWtBZ, rstat)
!      call HistoryPut( 'HumBZ',      xyzf_HumBZ    , rstat)
     
      ! Show CPU time 
      ! 
      call ClocksetPredict
      
    end if

    ! 長い時間ステップのループを回すための処置.
    ! Renew prognostic variables for next long time step integration.
    !
    pyz_VelXBl  = pyz_VelXNl;  pyz_VelXNl  = pyz_VelXAl
    xqz_VelYBl  = xqz_VelYNl;  xqz_VelYNl  = xqz_VelYAl
    xyr_VelZBl  = xyr_VelZNl;  xyr_VelZNl  = xyr_VelZAl
    xyz_PTempBl = xyz_PTempNl; xyz_PTempNl = xyz_PTempAl
    xyz_ExnerBl = xyz_ExnerNl; xyz_ExnerNl = xyz_ExnerAl
    xyz_KmBl    = xyz_KmNl;    xyz_KmNl    = xyz_KmAl
    xyz_KhBl    = xyz_KhNl;    xyz_KhNl    = xyz_KhAl
    xyz_CDensBl = xyz_CDensNl; xyz_CDensNl = xyz_CDensAl
    xyzf_QMixBl = xyzf_QMixNl; xyzf_QMixNl = xyzf_QMixAl

    ! 時刻の進行
    ! Progress time
    !
    call TimesetProgress

  end do

  !----------------------------------------------------
  ! 出力ファイルのクローズ
  ! Close out put files.
  !
  call HistoryFileio_finalize
  call ReStartFileio_finalize

  !----------------------------------------------------
  ! CPU 時間の計測終了
  ! 
  call ClocksetLoopStop
  call ClocksetClose

  !----------------------------------------------------
  ! MPI END
  !
  call MPIWrapperFinalize

  
contains
!-----------------------------------------------------------------------
  subroutine VariableAllocate
    !
    !初期化として, 配列を定義し, 値としてゼロを代入する.
    !

    !暗黙の型宣言禁止
    implicit none

    !配列割り当て
    allocate( pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax) )

    pyz_VelXBl  = 0.0d0;    pyz_VelXNl  = 0.0d0;    pyz_VelXAl  = 0.0d0
    pyz_VelXNs  = 0.0d0;    pyz_VelXAs = 0.0d0    
    xyz_VelXNl  = 0.0d0

    allocate( xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_VelYNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax) )

    xqz_VelYBl  = 0.0d0;    xqz_VelYNl  = 0.0d0;    xqz_VelYAl  = 0.0d0
    xqz_VelYNs  = 0.0d0;    xqz_VelYAs = 0.0d0    
    xyz_VelYNl  = 0.0d0

    allocate( xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_VelZNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax) )

    xyr_VelZBl  = 0.0d0;    xyr_VelZNl  = 0.0d0;    xyr_VelZAl  = 0.0d0
    xyr_VelZNs  = 0.0d0;    xyr_VelZAs = 0.0d0
    xyz_VelZNl  = 0.0d0

    allocate( xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax) )

    xyz_ExnerBl = 0.0d0;    xyz_ExnerNl = 0.0d0;    xyz_ExnerAl = 0.0d0
    xyz_ExnerNs = 0.0d0;    xyz_ExnerAs = 0.0d0

    allocate( xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAs(imin:imax,jmin:jmax,kmin:kmax) )

    xyz_PTempBl = 0.0d0;  xyz_PTempNl = 0.0d0;  xyz_PTempAl = 0.0d0
    xyz_PTempNs = 0.0d0;  xyz_PTempAs = 0.0d0

    allocate( xyz_CDensBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CDensNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CDensAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CDensNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_CDensAs(imin:imax,jmin:jmax,kmin:kmax) )

    xyz_CDensBl = 0.0d0;  xyz_CDensNl = 0.0d0;  xyz_CDensAl = 0.0d0
    xyz_CDensNs = 0.0d0;  xyz_CDensAs = 0.0d0

    allocate( xyz_KmBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_KmNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_KmAl(imin:imax,jmin:jmax,kmin:kmax) )

    xyz_KmBl    = 0.0d0;    xyz_KmNl    = 0.0d0;    xyz_KmAl    = 0.0d0

    allocate( xyz_KhBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_KhNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_KhAl(imin:imax,jmin:jmax,kmin:kmax) )

    xyz_KhBl    = 0.0d0;    xyz_KhNl    = 0.0d0;    xyz_KhAl    = 0.0d0

    allocate( xyzf_QMixBl(imin:imax,jmin:jmax,kmin:kmax,ncmax) )
    allocate( xyzf_QMixNl(imin:imax,jmin:jmax,kmin:kmax,ncmax) )
    allocate( xyzf_QMixAl(imin:imax,jmin:jmax,kmin:kmax,ncmax) )

    xyzf_QMixBl = 0.0d0;   xyzf_QMixNl = 0.0d0;   xyzf_QMixAl = 0.0d0

    allocate( pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DKmDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DCDensDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyzf_DQMixDtNl(imin:imax,jmin:jmax,kmin:kmax,ncmax) )

    pyz_DVelXDtNl = 0.0d0
    xqz_DVelYDtNl = 0.0d0
    xyr_DVelZDtNl = 0.0d0
    xyz_DKmDtNl   = 0.0d0
    xyz_DPTempDtNl = 0.0d0
    xyz_DExnerDtNl = 0.0d0
    xyz_DExnerDtNs = 0.0d0
    xyz_DCDensDtNl = 0.0d0
    xyzf_DQMixDtNl = 0.0d0

    ! 時間ループの初期化
    !
    t = 1

  end subroutine VariableAllocate


  !-----------------------------------------------------------------------
  subroutine AsselinTimeFilter
    !
    ! 時間フィルター; Asselin のタイムフィルターを利用
    !   t = 0.0 の場合には tfil = 0.0d0, それ以外は tfil = 1.0d-1
    !   (t = 0 の時はオイラー法で, それ以外はリープフロッグ法で積分するため)
    !
    use TimeSet, only: tfil

    ! 暗黙の型宣言禁止    
    !
    implicit none

    ! 作業変数
    !
    real(DP) :: tfil2

    tfil2 = 1.0d0 - 2.0d0 * tfil 

    pyz_VelXNl  = tfil * ( pyz_VelXBl  + pyz_VelXAl  ) + tfil2 * pyz_VelXNl
    xqz_VelYNl  = tfil * ( xqz_VelYBl  + xqz_VelYAl  ) + tfil2 * xqz_VelYNl
    xyr_VelZNl  = tfil * ( xyr_VelZBl  + xyr_VelZAl  ) + tfil2 * xyr_VelZNl
    xyz_PTempNl = tfil * ( xyz_PTempBl + xyz_PTempAl ) + tfil2 * xyz_PTempNl
    xyz_ExnerNl = tfil * ( xyz_ExnerBl + xyz_ExnerAl ) + tfil2 * xyz_ExnerNl 
    xyz_KmNl    = tfil * ( xyz_KmBl    + xyz_KmAl    ) + tfil2 * xyz_KmNl
    xyz_KhNl    = tfil * ( xyz_KhBl    + xyz_KhAl    ) + tfil2 * xyz_KhNl
    xyz_CDensNl = tfil * ( xyz_CDensBl + xyz_CDensAl ) + tfil2 * xyz_CDensNl
    xyzf_QMixNl = tfil * ( xyzf_QMixBl + xyzf_QMixAl ) + tfil2 * xyzf_QMixNl
    
  end subroutine AsselinTimeFilter
  

  !-----------------------------------------------------------------------
  subroutine MainInit
    ! 初期化手続き ; Initialize procedure 
    !    

    ! 暗黙の型宣言禁止    
    !
    implicit none

    character(STRING) :: cfgfile ! NAMELIST ファイル名 ; NAMELIST fine name
    integer, parameter:: OutputRank = 0

    ! gtool のメッセージの抑制
    ! rank0 のみ標準出力
    !
    call MessageSuppressMPI( OutputRank )

    ! MPI
    !
    call MPIWrapperInit

    ! 時間計測
    !
    call ClocksetInit      !初期化
    call ClocksetPreStart  !初期化ルーチン計測開始
    
    ! NAMELIST ファイル名の読み込み
    ! Loading NAMELIST file.
    !
    call argset_init( cfgfile ) !(out)

    ! NAMELIST ファイル名のモジュールへの登録
    ! Loading NAMELIST file.
    !
    call NmlutilInit( cfgfile ) !(in)
    
    ! 時間積分の初期化
    ! Initialization of time integration.
    !
    call timeset_init

    ! 格子点情報の初期化
    ! Initialization of grid arrangement.
    !
    call gridset_init
    
    ! 化学定数の初期化
    ! Initialization of chemical constatns.
    !
    call chemdata_init

    ! 定数の情報の初期化
    ! Initialization of constant variables.
    !
    call constants_init

    ! 軸の情報の初期化
    ! Initialization of axis variables.
    !
    call axesset_init
    
    ! I/O ファイル名の初期化
    ! Initialization of output file name. 
    !
    call fileset_init
    
    ! 湿潤過程共有変数の初期化
    ! Initialization of common variables for moist process.
    !
    call composition_init

    ! ヒストリーファイル・リスタートファイルの初期化
    ! Initialize restart & history files.
    !
    call HistoryFileio_init
    call ReStartFileio_init

    ! マージンの設定の初期化
    ! Initialization of margin
    !
    call SetMargin_Init
   
    ! 内部変数の初期化
    ! Initialization of internal variables.
    !
    call VariableAllocate

    ! フラグ処理
    !
    call CheckFlag

    ! 初期値の代入 
    ! * ReStartFile が設定されている場合にはファイルを読み込む. 
    !   設定されていない場合にはデフォルトの基本場と擾乱場を作る. 
    !
    ! Initial value set up.
    ! * Read restartfile if it is specified. If not, make default basic
    !   state and disturbance.
    !
    call MessageNotify( "M", "main", "Initial value setup." )

    ! 基本場, 擾乱場の初期値を netCDF ファイルから取得する.
    ! 
    select case (IDRestartMethod)

    case (IDRestartArare6)
      call ReStartFileio_BasicZ_Get
      call ReStartFileio_Var_Get(     &
        & pyz_VelXBl,  pyz_VelXNl,    & ! (out)
        & xqz_VelYBl,  xqz_VelYNl,    & ! (out)
        & xyr_VelZBl,  xyr_VelZNl,    & ! (out)
        & xyz_PTempBl, xyz_PTempNl,   & ! (out)
        & xyz_ExnerBl, xyz_ExnerNl,   & ! (out)
        & xyzf_QMixBl, xyzf_QMixNl,   & ! (out)
        & xyz_KmBl,    xyz_KmNl,      & ! (out)
        & xyz_KhBl,    xyz_KhNl,      & ! (out)
        & xyz_CDensBl, xyz_CDensNl    & ! (out)
        & )

    case (IDRestartArare4)
      call Arare4ReStartFileio_BasicZ_Get
      call Arare4ReStartFileio_Var_Get(  &
        & pyz_VelXBl,  pyz_VelXNl,    & ! (out)
        & xqz_VelYBl,  xqz_VelYNl,    & ! (out)
        & xyr_VelZBl,  xyr_VelZNl,    & ! (out)
        & xyz_PTempBl, xyz_PTempNl,   & ! (out)
        & xyz_ExnerBl, xyz_ExnerNl,   & ! (out)
        & xyzf_QMixBl, xyzf_QMixNl,   & ! (out)
        & xyz_KmBl,    xyz_KmNl,      & ! (out)
        & xyz_KhBl,    xyz_KhNl,      & ! (out)
        & xyz_CDensBl, xyz_CDensNl    & ! (out)
        & )

    end select

    ! 化学計算ルーチンの初期化
    ! Initialization of chemical routines.
    !
    call chemcalc_init
    
    ! 数値摩擦係数の初期化
    ! Initialization of numerical friction coefficient.
    !
    call Damping_Init

    ! 力学モジュールの初期化
    ! Initialization of dynamical modules
    !
    call Dynamics_init
    
    ! 負の湿潤量の補填計算の初期化
    ! Initialization of negative moist value correction.
    !
    call FillNegative_Init
   
    ! 物理過程の初期化. 
    ! フラグを元に必要なものだけ初期化する. 
    !
    call PhysicalProcess_init

    !-------------------------------------------------------------
    ! 基本場のファイル出力
    !

    ! モニタ用
    !
    call BasicFileio_Output
    
    ! 音波に対する CFL 条件のチェック
    ! CFL condtion check for sound wave.
    !
    call CFLCheckTimeShort( &
      & xyz_VelSoundBZ   & ! (in)
      & )
    
    ! 保存量のモニタリング
    !
    call EnergyMonit_init

    ! t=0 から数値積分を行う場合は, 最初の一回目はオイラー法を利用. 
    !
    if ( FlagInitialRun ) call TimesetDelTimeHalf

    ! 初期化ルーチンの CPU 時間の計測終了
    !
    call ClocksetPreStop
   
  end subroutine MainInit
  

  subroutine CheckFlag
    !
    ! 設定ファイルの指定に従って, 必要な物理過程の初期化ルーチンを呼ぶ
    !

    use dc_message,    only: MessageNotify
    use dc_iounit,     only : FileOpen    

    implicit none

    integer            :: unit                    !装置番号
    character(STRING)  :: FlagTurbMethod    = ""  !乱流拡散に関する設定
    character(STRING)  :: FlagRadMethod     = ""  !放射に関する設定
    character(STRING)  :: FlagCloudMethod   = ""  !雲微物理に関する設定
    character(STRING)  :: FlagSurfaceMethod = ""  !地表面過程に関する設定
    character(STRING)  :: FlagWindMethod    = ""  !速度場に関する設定
    character(STRING)  :: FlagDebugMethod   = ""  !デバッグ用のフラグ
    character(STRING)  :: FlagRestartMethod = ""  !読み込むリスタートファイルの形式
    character(STRING)  :: FlagDampingMethod = ""  !ダンピングに関する設定

    NAMELIST /deepconv_main_nml / &
      & FlagTurbMethod,           &
      & FlagRadMethod,            &
      & FlagCloudMethod,          &
      & FlagSurfaceMethod,        &
      & FlagWindMethod,           &
      & FlagDebugMethod,          &
      & FlagRestartMethod,        &
      & FlagDampingMethod

    ! デフォルト値の設定
    ! Default values settings
    !
    FlagTurbMethod       = "Nothing"
!!$    FlagTurbMethod    = "KW1978"
!!$    FlagTurbMethod    = "ConstKm"

    FlagRadMethod        = "Nothing"
!!$    FlagRadMethod     = "HeatConst"
!!$    FlagRadMethod     = "HeatVary"
!!$    FlagRadMethod     = "HeatBalance"

    FlagSurfaceMethod    = "Nothing"
!!$    FlagSurfaceMethod = "Diff"
!!$    FlagSurfaceMethod = "Bulk"
!!$    FlagSurfaceMethod = "Const"
!!$    FlagSurfaceMethod = "Baker1998"

    FlagCloudMethod      = "Nothing"
!!$    FlagCloudMethod   = "K1969"
!!$    FlagCloudMethod   = "MarsCond"

    FlagWindMethod       = "Nothing"
!!$    FlagWindMethod    = "Const"

    FlagDebugMethod      = "Nothing"
!!$    FlagDebugMethod   = "NoTendencyLong"

    FlagRestartMethod    = "Arare6"
!!$    FlagRestartMethod = "Arare4"

!!$    FlagDampingMethod    = "Nothing"
    !整合性のためにデフォルト値を "SpongeLayer" に設定. 
    FlagDampingMethod    = "SpongeLayer"  
!!$    FlagDampingMethod    = "SpongeLayerMeanFlow"

    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=deepconv_main_nml)
    close(unit)

    select case ( FlagTurbMethod )
    case ( "Nothing" )
    case ( "KW1978" )
      IDTurbMethod = IDTurbKW1978 
    case ( "ConstKm" )
      IDTurbMethod = IDTurbConstKm
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagTurbMethod=<%c> is not supported.', &
        & c1 = trim(FlagTurbMethod) )
    end select

    select case ( FlagRadMethod )
    case ( "Nothing" )
    case ( "HeatConst" )
      IDRadMethod = IDRadHeatConst 
    case ( "HeatVary" )
      IDRadMethod = IDRadHeatVary
    case ( "HeatBalance" )
      IDRadMethod = IDRadHeatBalance
    case ( "Sounding" )
      IDRadMethod = IDRadSounding
    case ( "Baker1998" )
      IDRadMethod = IDRadBaker1998
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagRadMethod=<%c> is not supported.', &
        & c1 = trim(FlagRadMethod) )
    end select

    select case ( FlagSurfaceMethod )
    case ( "Nothing" )
    case ( "Diff" )
      IDSurfaceMethod = IDSurfaceDiff
    case ( "Bulk" )
      IDSurfaceMethod = IDSurfaceBulk
    case ( "Const" )
      IDSurfaceMethod = IDSurfaceConst
    case ( "Baker1998" )
      IDSurfaceMethod = IDSurfaceBaker1998
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagSurfaceMethod=<%c> is not supported.', &
        & c1 = trim(FlagSurfaceMethod) )
    end select

    select case ( FlagCloudMethod )
    case ( "Nothing" )
    case ( "K1969" )
      IDCloudMethod = IDCloudK1969
    case ( "MarsCond" )
      IDCloudMethod = IDCloudMarsCond
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagCloudMethod=<%c> is not supported.', &
        & c1 = trim(FlagCloudMethod) )
    end select

    select case ( FlagWindMethod )
    case ( "Nothing" )
    case ( "Const" )
      IDDebugMethod = IDDebugWindConst
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagWindMethod=<%c> is not supported.', &
        & c1 = trim(FlagWindMethod) )
    end select

    select case ( FlagDebugMethod )
    case ( "Nothing" )
    case ( "WindConst" )
      IDDebugMethod = IDDebugWindConst
    case ( "NoTendencyLong" )
      IDDebugMethod = IDDebugNoTendencyLong
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagDebugdMethod=<%c> is not supported.', &
        & c1 = trim(FlagDebugMethod) )
    end select

    select case ( FlagRestartMethod )
    case ( "Arare4" )
      IDRestartMethod = IDRestartArare4
    case default
      IDRestartMethod = IDRestartArare6
    end select

    select case ( FlagDampingMethod )
    case ( "Nothing" )
    case ( "SpongeLayer" )
      IDDampingMethod = IDDampingSpongeLayer
    case ( "SpongeLayerMeanFlow" )
      IDDampingMethod = IDDampingSpongeLayerMeanFlow
    case default
      call MessageNotify( 'E', prog_name, &
        & 'FlagDampingMethod=<%c> is not supported.', &
        & c1 = trim(FlagDampingMethod) )
    end select
    
    call MessageNotify( "M", "main", "FlagTurbMethod    = %c", c1=trim(FlagTurbMethod))
    call MessageNotify( "M", "main", "IDTurbMethod      = %d", i=(/IDTurbMethod/))
    call MessageNotify( "M", "main", "FlagRadMethod     = %c", c1=trim(FlagRadMethod))
    call MessageNotify( "M", "main", "IDRadMethod       = %d", i=(/IDRadMethod/))
    call MessageNotify( "M", "main", "FlagSurfaceMethod = %c", c1=trim(FlagSurfaceMethod))
    call MessageNotify( "M", "main", "IDSurfaceMethod   = %d", i=(/IDSurfaceMethod/))
    call MessageNotify( "M", "main", "FlagCloudMethod   = %c", c1=trim(FlagCloudMethod))
    call MessageNotify( "M", "main", "IDCloudMethod     = %d", i=(/IDCloudMethod/))
    call MessageNotify( "M", "main", "FlagDebugMethod   = %c", c1=trim(FlagDebugMethod))
    call MessageNotify( "M", "main", "IDDebugMethod     = %d", i=(/IDDebugMethod/))
    call MessageNotify( "M", "main", "FlagRestartMethod = %c", c1=trim(FlagRestartMethod))
    call MessageNotify( "M", "main", "IDRestartMethod   = %d", i=(/IDRestartMethod/))
    call MessageNotify( "M", "main", "FlagDampingMethod = %c", c1=trim(FlagDampingMethod))
    call MessageNotify( "M", "main", "IDDampingMethod   = %d", i=(/IDDampingMethod/))
        
  end subroutine CheckFlag


  subroutine PhysicalProcess_init

    select case ( IDTurbMethod )
    case ( IDTurbKW1978 )
      call turbulence_kw1978_init
    case ( IDTurbConstKm )
      call turbulence_constKm_init
    end select

    select case ( IDRadMethod )
    case ( IDRadHeatConst )
      call Radiation_Simple_init
    case ( IDRadHeatVary )
      call Radiation_Simple_init
    case ( IDRadHeatBalance )
      call radiation_heatbalance_init
    case ( IDRadSounding )
      call radiation_sounding_init
    case ( IDRadBaker1998 )
      call Radiation_BAKER1998_init
    end select

    select case ( IDSurfaceMethod )
    case ( IDSurfaceDiff )
      call Surfaceflux_Diff_init
    case ( IDSurfaceBulk )
      call Surfaceflux_Bulk_init
    case ( IDSurfaceConst )
      call Surfaceflux_Const_init
    case ( IDSurfaceBaker1998 )
      call Surfaceflux_baker1998_init
    end select

    select case ( IDCloudMethod )
    case ( IDCloudK1969 )
      call cloudphys_K1969_init
    case ( IDCloudMarsCond )
      call cloudphys_marscond_init
    end select

  end subroutine PhysicalProcess_init


  subroutine turbulence_integrate

    ! モジュール呼び出し
    !
    use timeset,              only : DelTimeLong
    use setmargin,            only : SetMargin_xyz
    use Turbulence_kw1978_v2, only : KmMax

    ! 暗黙の型宣言禁止
    !
    implicit none

    !
    ! Km の積分 (leap-frog)
    !
    xyz_KmAl = xyz_KmBl + 2.0d0 * DelTimeLong * xyz_DKmDtNl

    ! 値の上限下限の設定
    !  * 値は正になることを保証する
    !
    ! Upper and lower bound value are specified.
    !
    xyz_KmAl = max( 0.0d0, min( xyz_KmAl, KmMax ) )
    
    ! 境界条件 ; Boundary condition
    !
    call SetMargin_xyz( xyz_KmAl )  

    ! スカラーに対する渦拡散係数の計算 
    ! Specify turbulent diffusion coefficient for scalar variables.
    !
    xyz_KhAl = 3.0d0 * xyz_KmAl
         
  end subroutine turbulence_integrate



  subroutine PTemp_integrate

    ! モジュール呼び出し
    !
    use timeset,   only: DelTimeLong
    use setmargin, only : SetMargin_xyz

    ! 暗黙の型宣言禁止
    !
    implicit none
    
    !
    ! 温位の積分 (leap-frog)
    !
    xyz_PTempAl = xyz_PTempBl + 2.0d0 * DelTimeLong * xyz_DPTempDtNl

    ! 境界条件 ; Boundary condition
    !
    call SetMargin_xyz( xyz_PTempAl )  

  end subroutine PTemp_integrate
  


  subroutine QMix_integrate

    ! モジュール呼び出し
    !
    use timeset,      only : DelTimeLong
    use setmargin,    only : SetMargin_xyzf
    use fillnegative, only : xyza_FillNegative_xyza

    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 作業変数
    !    
    real(DP) :: xyzf_QMixWork(imin:imax,jmin:jmax,kmin:kmax,ncmax) 
    real(DP) :: xyzf_DQMixDt1(imin:imax,jmin:jmax,kmin:kmax,ncmax) 

    !
    ! 混合比の積分 (leap-frog)
    !
    xyzf_QMixAl = xyzf_QMixBl + 2.0d0 * DelTimeLong * xyzf_DQMixDtNl
    
    ! 境界条件 ; Boundary condition
    !
    call SetMargin_xyzf( xyzf_QMixAl )  

    ! 移流によって負になった部分を埋める
    ! Negative values due to advection are corrected.
    !
    xyzf_QMixWork = xyzf_QMixAl
    xyzf_QMixAl   = xyza_FillNegative_xyza( xyzf_QMixWork ) 
    
    ! 埋めた/削った量を保管. 単位時間当たりの量に変換. 
    ! Correction value is stored.
    !
    xyzf_DQMixDt1 = (xyzf_QMixAl - xyzf_QMixWork) / ( 2.0d0 * DelTimeLong )
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtFill1', xyzf_DQMixDt1(1:nx,1:ny,1:nz,s))
    end do
    
    ! 埋め切れなかった部分をゼロにする. Fill2 に保管
    ! Negative values mixing ratios are corrected.
    !
    xyzf_QMixWork = xyzf_QMixAl
    xyzf_QMixAl = max( - xyzf_QMixBZ, xyzf_QMixWork )
    
    xyzf_DQMixDt1 = (xyzf_QMixAl - xyzf_QMixWork) / ( 2.0d0 * DelTimeLong )
    
    do s = 1, ncmax
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(s))//'DtFill2', xyzf_DQMixDt1(1:nx,1:ny,1:nz,s))
    end do
    
    ! 境界条件 Boundary condition
    !
    call SetMargin_xyzf( xyzf_QMixAl )  ! (inout)   
         
  end subroutine QMix_integrate

      
end program deepconv_arare
