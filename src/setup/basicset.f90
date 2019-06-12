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
  use dc_types, only:  DP

  !暗黙の型宣言禁止
  implicit none

  !save 属性
  private

  !Public Interface
  real(DP), allocatable, save, public :: xyz_DensBZ(:,:,:)    !密度
  real(DP), allocatable, save, public :: pyz_DensBZ(:,:,:)    !密度
  real(DP), allocatable, save, public :: xqz_DensBZ(:,:,:)    !密度
  real(DP), allocatable, save, public :: xyr_DensBZ(:,:,:)    !密度
  real(DP), allocatable, save, public :: xyz_PressBZ(:,:,:)   !無次元圧力
  real(DP), allocatable, save, public :: xyz_ExnerBZ(:,:,:)   !無次元圧力
  real(DP), allocatable, save, public :: xyz_TempBZ(:,:,:)    !温度
  real(DP), allocatable, save, public :: xyz_PTempBZ(:,:,:)   !温位
  real(DP), allocatable, save, public :: xyr_PTempBZ(:,:,:)   !温位
  real(DP), allocatable, save, public :: xyz_VPTempBZ(:,:,:)  !仮温位
  real(DP), allocatable, save, public :: pyz_VPTempBZ(:,:,:)  !仮温位
  real(DP), allocatable, save, public :: xqz_VPTempBZ(:,:,:)  !仮温位
  real(DP), allocatable, save, public :: xyr_VPTempBZ(:,:,:)  !仮温位
  real(DP), allocatable, save, public :: xyz_VelSoundBZ(:,:,:)!音速
  real(DP), allocatable, save, public :: xyz_VelSW(:,:,:)     !音速
  real(DP), allocatable, save, public :: xyzf_QMixBZ(:,:,:,:) !凝縮成分混合比
  real(DP), allocatable, save, public :: xyz_EffMolWtBZ(:,:,:)!分子量効果
  real(DP), allocatable, save, public :: xyz_QMixBZPerMolWt(:,:,:)
                              !基本場の混合比 / 分子量 の総和
  real(DP), allocatable, save, public :: xyr_QMixBZPerMolWt(:,:,:)
                              !基本場の混合比 / 分子量 の総和
  real(DP), allocatable, save, public :: xyz_QMixBZ(:,:,:)
                              !基本場の混合比の総和
  real(DP), allocatable, save, public :: xyr_QMixBZ(:,:,:)
                              !基本場の混合比の総和

  public basicset_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine basicset_init(                                     &
    &   xyz_Press, xyz_Exner, xyz_Temp, xyz_PTemp, xyz_Dens,    & 
    &   xyz_VelSound, xyzf_QMix, xyz_EffMolWt                   &
    & )
    !
    != 基本場の値を外部から取得する. 
    !
    ! dry の場合は混合比や乾燥成分と凝結成分の存在比を陽に使わないので,
    ! これらの変数は optional にしている. 
    !

    !モジュール呼び出し
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
    real(DP), intent(in) :: xyz_Press(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Dens(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_VelSound(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in), optional :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(in), optional :: xyz_EffMolWt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)             :: xyzf_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax,GasNum)
    integer              :: s, k

    ! 配列の初期化
    call basicset_array_init

    !値の入力
    xyz_PressBZ    = xyz_Press
    xyz_ExnerBZ    = xyz_Exner
    xyz_TempBZ     = xyz_Temp
    xyz_PTempBZ    = xyz_PTemp
    xyz_DensBZ     = xyz_Dens
    xyz_VelSoundBZ = xyz_VelSound
    xyz_VelSW      = xyz_VelSound  ! 配列名を短くするために. 

    if (present(xyz_EffMolWt)) xyz_EffMolWtBZ = xyz_EffMolWt

    if (present(xyzf_QMix)) then 
      xyzf_QMixBZ = xyzf_QMix
      xyzf_QMixBZPerMolWt = 0.0d0
      
      do s = 1, GasNum
        xyzf_QMixBZPerMolWt(:,:,:,s) = &
          & xyzf_QMixBZ(:,:,:,IdxG(s)) / MolWtWet(IdxG(s))
      end do
      
      xyz_QMixBZPerMolWt = sum(xyzf_QMixBZPerMolWt, 4) 
      xyz_QMixBZ         = sum(xyzf_QMixBZ, 4) 

      do k = kmin, kmax-1
        xyr_QMixBZPerMolWt(:,:,k) &
          & = ( xyz_QMixBZPerMolWt(:,:,k+1) + xyz_QMixBZPerMolWt(:,:,k) ) * 0.5d0
        xyr_QMixBZ(:,:,k) = ( xyz_QMixBZ(:,:,k+1) + xyz_QMixBZ(:,:,k) )   * 0.5d0
      end do
    end if
    
    xyz_VPTempBZ = xyz_PTempBZ / xyz_EffMolWtBZ

    do k = kmin, kmax-1    
      xyr_PTempBZ(:,:,k)  = (  xyz_PTempBZ(:,:,k+1) +  xyz_PTempBZ(:,:,k) ) * 0.5d0
      xyr_VPTempBZ(:,:,k) = ( xyz_VPTempBZ(:,:,k+1) + xyz_VPTempBZ(:,:,k) ) * 0.5d0
      xyr_DensBZ(:,:,k)   = (   xyz_DensBZ(:,:,k+1) +   xyz_DensBZ(:,:,k) ) * 0.5d0
    end do

    ! 水平一様なので, 平均操作は必要ない
    !
    pyz_DensBZ   = xyz_Dens
    xqz_DensBZ   = xyz_Dens
    pyz_VPTempBZ = xyz_VPTempBZ 
    xqz_VPTempBZ = xyz_VPTempBZ 
  
  contains

    subroutine basicset_array_init
      !
      ! *BZ な配列の初期化
      !
      
      allocate( & 
        & xyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & pyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xqz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyr_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_PressBZ(imin:imax,jmin:jmax,kmin:kmax),       &
        & xyz_ExnerBZ(imin:imax,jmin:jmax,kmin:kmax),       &
        & xyz_TempBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_PTempBZ (imin:imax,jmin:jmax,kmin:kmax),      &
        & xyr_PTempBZ (imin:imax,jmin:jmax,kmin:kmax),      &
        & xyz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & pyz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xqz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xyr_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xyz_VelSoundBZ(imin:imax,jmin:jmax,kmin:kmax),    &
        & xyz_VelSW(imin:imax,jmin:jmax,kmin:kmax),         &
        & xyzf_QMixBZ(imin:imax,jmin:jmax,kmin:kmax,ncmax), &
        & xyz_EffMolWtBZ(imin:imax,jmin:jmax,kmin:kmax)     &
        & )
      allocate( &
        & xyr_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax),&
        & xyz_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax),&
        & xyr_QMixBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_QMixBZ(imin:imax,jmin:jmax,kmin:kmax)         &
        & )
      
      ! 値の確定
      xyz_DensBZ     = 0.0d0
      pyz_DensBZ     = 0.0d0
      xqz_DensBZ     = 0.0d0
      xyr_DensBZ     = 0.0d0
      xyz_PressBZ    = 0.0d0   
      xyz_ExnerBZ    = 0.0d0   
      xyz_TempBZ     = 0.0d0   
      xyz_PTempBZ    = 0.0d0
      xyz_VPTempBZ   = 0.0d0
      pyz_VPTempBZ   = 0.0d0
      xqz_VPTempBZ   = 0.0d0
      xyr_VPTempBZ   = 0.0d0
      xyz_VelSoundBZ = 0.0d0
      xyz_VelSW      = 0.0d0
      xyzf_QMixBZ    = 0.0d0
      xyz_EffMolWtBZ = 1.0d0 ! dry の場合は必ず 1.0
      
      xyr_QMixBZPerMolWt = 0.0d0
      xyz_QMixBZPerMolWt = 0.0d0
      xyr_QMixBZ = 0.0d0
      xyz_QMixBZ = 0.0d0
      
    end subroutine basicset_array_init

  end subroutine basicset_init
    
end module basicset
