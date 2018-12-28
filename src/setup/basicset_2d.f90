!= デフォルトの基本場を設定するための変数参照型モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: basicset.f90,v 1.18 2014/07/08 01:05:32 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module basicset
  !
  != デフォルトの基本場を設定するための変数参照型モジュール
  !
  ! hogeBZ となる配列は, 外部から値を入力することで初期化を行うことに注意. 
  !

  !モジュール読み込み
  use dc_types,   only:  DP

  !暗黙の型宣言禁止
  implicit none

  !save 属性
  private

  !Public Interface
  real(DP), allocatable, save, public :: xz_DensBZ(:,:)     !密度
  real(DP), allocatable, save, public :: pz_DensBZ(:,:)     !密度
  real(DP), allocatable, save, public :: xqz_DensBZ(:,:)    !密度
  real(DP), allocatable, save, public :: xr_DensBZ(:,:)     !密度
  real(DP), allocatable, save, public :: xz_PressBZ(:,:)    !無次元圧力
  real(DP), allocatable, save, public :: xz_ExnerBZ(:,:)    !無次元圧力
  real(DP), allocatable, save, public :: xz_TempBZ(:,:)     !温度
  real(DP), allocatable, save, public :: xz_PTempBZ(:,:)    !温位
  real(DP), allocatable, save, public :: xr_PTempBZ(:,:)    !温位
  real(DP), allocatable, save, public :: xz_VPTempBZ(:,:)   !仮温位
  real(DP), allocatable, save, public :: pz_VPTempBZ(:,:)   !仮温位
  real(DP), allocatable, save, public :: xqz_VPTempBZ(:,:)  !仮温位
  real(DP), allocatable, save, public :: xr_VPTempBZ(:,:)   !仮温位
  real(DP), allocatable, save, public :: xz_VelSoundBZ(:,:) !音速
  real(DP), allocatable, save, public :: xz_VelSW(:,:)      !音速
  real(DP), allocatable, save, public :: xzf_QMixBZ(:,:,:)  !凝縮成分混合比
  real(DP), allocatable, save, public :: xz_EffMolWtBZ(:,:) !分子量効果
  real(DP), allocatable, save, public :: xz_QMixBZPerMolWt(:,:)
                              !基本場の混合比 / 分子量 の総和
  real(DP), allocatable, save, public :: xr_QMixBZPerMolWt(:,:)
                              !基本場の混合比 / 分子量 の総和
  real(DP), allocatable, save, public :: xz_QMixBZ(:,:)
                              !基本場の混合比の総和
  real(DP), allocatable, save, public :: xr_QMixBZ(:,:)
                              !基本場の混合比の総和

  public basicset_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine basicset_init(                                &
    &   xz_Press, xz_Exner, xz_Temp, xz_PTemp, xz_Dens,    & 
    &   xz_VelSound, xzf_QMix, xz_EffMolWt                 &
    & )
    !
    != 基本場の値を外部から取得する. 
    !
    ! dry の場合は混合比や乾燥成分と凝結成分の存在比を陽に使わないので,
    ! これらの変数は optional にしている. 
    !

    !モジュール読み込み
    use dc_types,   only:  DP
    use gridset,    only:  imin,       &! 配列の X 方向の下限
      &                    imax,       &! 配列の X 方向の上限
      &                    jmin,       &! 配列の Z 方向の下限
      &                    jmax,       &! 配列の Z 方向の上限
      &                    kmin,       &! 配列の Z 方向の下限
      &                    kmax,       &! 配列の Z 方向の上限
      &                    ncmax        ! 化学種の数
    use composition, only: GasNum,     &!気体の数
      &                    IdxG,       &!気体の配列添え字
      &                    MolWtWet     !凝縮成分の分子量
    
    !暗黙の型宣言禁止
    implicit none
    
    !入力変数
    real(DP), intent(in) :: xz_Press(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Exner(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Temp(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_PTemp(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Dens(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_VelSound(imin:imax,kmin:kmax)
    real(DP), intent(in), optional :: xzf_QMix(imin:imax,kmin:kmax,ncmax)
    real(DP), intent(in), optional :: xz_EffMolWt(imin:imax,kmin:kmax)

    real(DP)             :: xzf_QMixBZPerMolWt(imin:imax,kmin:kmax,GasNum)
    integer              :: s

    ! 配列の初期化
    call basicset_array_init
    
    !値の入力
    xz_PressBZ    = xz_Press
    xz_ExnerBZ    = xz_Exner
    xz_TempBZ     = xz_Temp
    xz_PTempBZ    = xz_PTemp
    xz_DensBZ     = xz_Dens
    xz_VelSoundBZ = xz_VelSound
    xz_VelSW      = xz_VelSound  ! 配列名を短くするために. 

    if (present(xz_EffMolWt)) xz_EffMolWtBZ = xz_EffMolWt

    if (present(xzf_QMix)) then 
      xzf_QMixBZ = xzf_QMix
      xzf_QMixBZPerMolWt = 0.0d0

      do s = 1, GasNum
        xzf_QMixBZPerMolWt(:,:,s) = &
          & xzf_QMixBZ(:,:,IdxG(s)) / MolWtWet(IdxG(s))
      end do

      xz_QMixBZPerMolWt = sum(xzf_QMixBZPerMolWt, 3) 
      xz_QMixBZ         = sum(xzf_QMixBZ, 3) 
      
      xr_QMixBZPerMolWt = xr_xz( xz_QMixBZPerMolWt )
      xr_QMixBZ         = xr_xz( xz_QMixBZ )

    end if

    xz_VPTempBZ = xz_PTempBZ / xz_EffMolWtBZ

    xr_PTempBZ  = xr_xz( xz_PTempBZ )
    xr_DensBZ   = xr_xz( xz_DensBZ )
    xr_VPTempBZ = xr_xz( xz_VPTempBZ )

    ! 水平一様なので, 平均操作は必要ない
    !
    pz_DensBZ    = xz_Dens
    xqz_DensBZ   = xz_Dens
    pz_VPTempBZ  = xz_VPTempBZ 
    xqz_VPTempBZ = xz_VPTempBZ 

  contains

    subroutine basicset_array_init
      !
      ! *BasicZ な配列の初期化
      !

      !暗黙の型宣言禁止
      implicit none

      allocate( & 
        & xz_DensBZ(imin:imax,kmin:kmax),        &
        & pz_DensBZ(imin:imax,kmin:kmax),        &
        & xqz_DensBZ(imin:imax,kmin:kmax),       &
        & xr_DensBZ(imin:imax,kmin:kmax),        &
        & xz_PressBZ(imin:imax,kmin:kmax),       &
        & xz_ExnerBZ(imin:imax,kmin:kmax),       &
        & xz_TempBZ(imin:imax,kmin:kmax),        &
        & xz_PTempBZ (imin:imax,kmin:kmax),      &
        & xr_PTempBZ (imin:imax,kmin:kmax),      &
        & xz_VPTempBZ (imin:imax,kmin:kmax),     &
        & pz_VPTempBZ (imin:imax,kmin:kmax),     &
        & xqz_VPTempBZ (imin:imax,kmin:kmax),    &
        & xr_VPTempBZ (imin:imax,kmin:kmax),     &
        & xz_VelSoundBZ(imin:imax,kmin:kmax),    &
        & xz_VelSW(imin:imax,kmin:kmax),         &
        & xzf_QMixBZ(imin:imax,kmin:kmax,ncmax), &
        & xz_EffMolWtBZ(imin:imax,kmin:kmax)     &
        & )
      allocate( &
        & xr_QMixBZPerMolWt(imin:imax,kmin:kmax),&
        & xz_QMixBZPerMolWt(imin:imax,kmin:kmax),&
        & xr_QMixBZ(imin:imax,kmin:kmax),        &
        & xz_QMixBZ(imin:imax,kmin:kmax)         &
        & )
      
      ! 値の確定
      xz_DensBZ     = 0.0d0
      pz_DensBZ     = 0.0d0
      xqz_DensBZ    = 0.0d0
      xr_DensBZ     = 0.0d0
      xz_PressBZ    = 0.0d0
      xz_ExnerBZ    = 0.0d0
      xz_TempBZ     = 0.0d0
      xz_PTempBZ    = 0.0d0
      xr_PTempBZ    = 0.0d0
      xz_VPTempBZ   = 0.0d0
      pz_VPTempBZ   = 0.0d0
      xqz_VPTempBZ  = 0.0d0
      xr_VPTempBZ   = 0.0d0
      xz_VelSoundBZ = 0.0d0
      xz_VelSW      = 0.0d0
      xzf_QMixBZ    = 0.0d0
      xz_EffMolWtBZ = 1.0d0 ! dry の場合は必ず 1.0
      
      xr_QMixBZPerMolWt = 0.0d0
      xz_QMixBZPerMolWt = 0.0d0
      xr_QMixBZ = 0.0d0
      xz_QMixBZ = 0.0d0
      
    end subroutine basicset_array_init
    

    function xr_xz(xz_Var)
      !
      ! 平均操作を行い z 方向半整数格子点の配列値を整数格子点上へ返す
      
      implicit none
      
      real(DP),intent(in) :: xz_Var(imin:imax,kmin:kmax) 
      real(DP)            :: xr_xz(imin:imax,kmin:kmax)
      integer             :: kz
      
      do kz = kmin, kmax-1
        xr_xz(:,kz) = ( xz_Var(:,kz+1) + xz_Var(:,kz) ) * 5.0d-1
      end do
      
      xr_xz(:,kmax) = xr_xz(:,kmax-1)
      
    end function xr_xz
    
    
    function xz_xr(xr_Var)
      !
      ! 平均操作を行い z 方向半整数格子点の配列値を整数格子点上へ返す
      
      implicit none
      
      real(DP),intent(in) :: xr_Var(imin:imax,kmin:kmax) 
      real(DP)            :: xz_xr(imin:imax,kmin:kmax)
      integer             :: kz
      
      do kz = kmin+1, kmax
        xz_xr(:,kz) = ( xr_Var(:,kz) + xr_Var(:,kz-1) ) * 5.0d-1
      end do
      
      xz_xr(:,kmin) = xz_xr(:,kmin+1)
      
    end function xz_xr

  end subroutine basicset_init

end module basicset
