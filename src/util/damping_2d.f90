!= スポンジ層モジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: damping.f90,v 1.8 2010-08-13 07:18:17 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module Damping
  !
  ! スポンジ層 (境界付近で波の反射を抑え吸収するための層) における
  ! 減衰率とその計算を行うためのパッケージ型モジュール
  !

  !モジュール読み込み
  use dc_types, only : DP

  !暗黙の型宣言禁止
  implicit none

  !private 属性を指定
  private 
  
  !関数には public 属性を指定
  public Damping_Init
  public DampSound_Init
  public DampSponge_Init
  public DampSponge_xz
  public DampSponge_xr
  public DampSponge_pz
  public xz_Sponge
  public xr_Sponge
  public pz_Sponge

  !変数定義
  real(DP), save, public   :: DampSound  = 0.0d0   !音波減衰項の減衰係数
  real(DP), save, private  :: EFTime     = 100.0d0 !スポンジ層の e-folding time
  real(DP), save, private  :: DampDepthH = 0.0d0   !スポンジ層の厚さ(水平方向)
  real(DP), save, private  :: DampDepthV = 0.0d0   !スポンジ層の厚さ(鉛直方向)
  real(DP), allocatable, save, private :: xz_DampRateH(:,:) !ss 格子点減衰係数(水平方向)
  real(DP), allocatable, save, private :: xz_DampRateV(:,:) !ss 格子点減衰係数(鉛直方向)
  real(DP), allocatable, save, private :: pz_DampRateH(:,:) !fs 格子点減衰係数(水平方向)
  real(DP), allocatable, save, private :: pz_DampRateV(:,:) !fs 格子点減衰係数(鉛直方向)
  real(DP), allocatable, save, private :: xr_DampRateH(:,:) !sf 格子点減衰係数(水平方向)
  real(DP), allocatable, save, private :: xr_DampRateV(:,:) !sf 格子点減衰係数(鉛直方向)
  
contains 
  
!!!------------------------------------------------------------------------!!!
  subroutine Damping_Init( cfgfile ) 
    !
    ! 音波減衰項とスポンジ層の減衰係数の初期化
    ! 

    !モジュール呼び出し
    use dc_types,  only : DP
 
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    character(*), intent(in) :: cfgfile
    real(DP)                 :: Alpha    !音波減衰項の係数
    real(DP)                 :: Time     !
    real(DP)                 :: DepthH   !スポンジ層の厚さ(水平方向)
    real(DP)                 :: DepthV   !スポンジ層の厚さ(鉛直方向)

    !NAMELIST から取得
    NAMELIST /damping/ Alpha, Time, DepthH, DepthV
    open (10, FILE=cfgfile)
    read(10, NML=damping)
    close(10)

    !初期化
    call DampSound_Init( Alpha ) 
    call DampSponge_Init( Time, DepthH, DepthV )

  end subroutine Damping_Init


!!!------------------------------------------------------------------------!!!
  subroutine DampSound_Init( damp ) 
    !
    ! 音波減衰項の減衰係数の初期化
    ! 

    !モジュール呼び出し
    use dc_types,    only : DP
    use dc_message,  only : MessageNotify
    use mpi_wrapper, only : myrank
    use axesset,     only : DelX,        &! x 方向の格子点間隔
      &                     DelZ          ! z 方向の格子点間隔
    use timeset,     only : DelTimeShort  !短い時間ステップ

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in)  :: damp

    !-------------------------------------------------------------------
    ! 音波減衰項の減衰率   Min(DelX, DelZ) ** 2.0 に比例
    !-------------------------------------------------------------------
    DampSound = Damp * ( Min(DelX, DelZ) ** 2.0d0 ) / DelTimeShort
    
    !-----------------------------------------------------------------    
    ! 値の確認
    !-----------------------------------------------------------------
    if (myrank == 0) then 
      call MessageNotify( "M", &
        & "DampSound_init", "DampSound = %f", d=(/DampSound/) )
    end if
  end subroutine DampSound_Init


!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_Init( Time, DepthH, DepthV )
    !
    ! スポンジ層の減衰係数を決める
    !

    !モジュール呼び出し
    use dc_types,    only: DP
    use dc_message,  only: MessageNotify
    use mpi_wrapper, only: myrank
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限
    use axesset,     only: DelX,          &! x 方向の格子点間隔
      &                    DelZ,          &! z 方向の格子点間隔
      &                    s_X,           &!X 座標軸(スカラー格子点)
      &                    s_Z,           &!Z 座標軸(スカラー格子点)
      &                    f_X,           &!X 座標軸(フラックス格子点)
      &                    f_Z,           &!Z 座標軸(フラックス格子点)
      &                    XMax,          &!X 座標の最大値
      &                    ZMax            !Z 座標の最大値

    !暗黙の型宣言禁止
    implicit none

    !入力変数
    real(DP), intent(in)   :: Time     !スポンジ層の e-folding time
    real(DP), intent(in)   :: DepthH   !スポンジ層の厚さ(水平方向)
    real(DP), intent(in)   :: DepthV   !スポンジ層の厚さ(鉛直方向)
    real(DP), parameter    :: Pi =3.1415926535897932385d0   !円周率
    integer               :: i, k

    !初期化
    allocate( &
      & xz_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xz_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax), &
      & pz_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & pz_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xr_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xr_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax)    )
    xz_DampRateH = 0.0d0
    xz_DampRateV = 0.0d0
    pz_DampRateH = 0.0d0
    pz_DampRateV = 0.0d0
    xr_DampRateH = 0.0d0
    xr_DampRateV = 0.0d0

    !値の入力
    EFTime     = Time
    DampDepthH = DepthH
    DampDepthV = DepthV
    
    !-----------------------------------------------------------------    
    ! スポンジ層の減衰率
    !-----------------------------------------------------------------
    !水平方向の東側・西側境界
    if ( DampDepthH < DelX ) then 
      call MessageNotify( "W", &
        & "DampSponge_Init", "DampDepthH is too thin. DelX is %f", d=(/DelX/))

    else
      do i = DimXMin, DimXMax
        !スカラー格子点の西側境界
        if ( s_X(i) < DampDepthH) then 
          xz_DampRateH(i,:) = ((1.0d0 - s_X(i) / DampDepthH) ** 3.0d0) / EFTime
        end if
        
        !フラックス格子点の西側境界
        if ( f_X(i) < DampDepthH) then 
          pz_DampRateH(i,:) = ((1.0d0 - f_X(i) / DampDepthH) ** 3.0d0) / EFTime
         end if
        
        !スカラー格子点の東側境界    
        if ( s_X(i) > ( XMax - DampDepthH ) ) then 
          xz_DampRateH(i,:) = &
            & ((1.0d0 - (XMax - s_X(i)) / DampDepthH) ** 3.0d0) / EFTime 
        end if
        
        !フラックス格子点の東側境界    
        if ( f_X(i) > ( XMax - DampDepthH ) ) then 
          pz_DampRateH(i,:) = &
            & ((1.0d0 - (XMax - f_X(i)) / DampDepthH) ** 3.0d0) / EFTime 
        end if
      end do
    end if
    !sf と ss は X 方向に関しては同じ
    xr_DampRateH  = xz_DampRateH
    
    !鉛直方向の上部境界    
    if ( DampDepthV < DelZ ) then 
      call MessageNotify( "W", &
        & "DampSponge_Init", "DampDepthV is too thin. DelZ is %f", d=(/DelZ/) )
      
    else
      do k = DimZMin, DimZMax
        !スカラー格子点
        if ( s_Z(k) >= ( ZMax - DampDepthV ) ) then 
          xz_DampRateV(:,k) =  &
            & (1.0d0 - dcos(Pi * (s_Z(k) - ZMax + DampDepthV) / DampDepthV)) &
            &  / EFTime 
        end if
        
        !フラックス格子点
        if ( f_Z(k) >= ( ZMax - DampDepthV ) ) then 
          xr_DampRateV(:,k) =  &
            & (1.0d0 - dcos(Pi * (f_Z(k) - ZMax + DampDepthV)/ DampDepthV)) &
            &  / EFTime 
        end if
      end do
    end if
    !fs と ss は Z 方向に関しては同じ
    pz_DampRateV  = xz_DampRateV
    
    !-----------------------------------------------------------------    
    ! 値の確認
    !-----------------------------------------------------------------
    if (myrank == 0) then 
      call MessageNotify( "M", "DampSponge_Init", "EFTime = %f", d=(/EFTime/) )
      call MessageNotify( "M", "DampSponge_Init", "DampDepthH = %f", d=(/DampDepthH/) )
      call MessageNotify( "M", "DampSponge_Init", "DampDepthV = %f", d=(/DampDepthV/) )  
    end if
  end subroutine DampSponge_Init
  

!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_xz(xz_VarA, xz_VarB, DelTime)
    !
    ! ss 格子点に対するスポンジ層
    !

    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(inout):: xz_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: xz_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !スポンジ層によるダンピングを計算
    xz_Var  =  xz_VarA - ( xz_DampRateH + xz_DampRateV ) * xz_VarB * DelTime
    xz_VarA = xz_Var
    
  end subroutine DampSponge_xz


!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_xr(xr_VarA, xr_VarB, DelTime)
    !
    ! sf 格子点に対するスポンジ層
    !
    
    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(inout):: xr_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: xr_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: xr_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !スポンジ層によるダンピングを計算  
    xr_Var  = xr_VarA - ( xr_DampRateH + xr_DampRateV )* xr_VarB * DelTime
    xr_VarA = xr_Var
    
  end subroutine DampSponge_xr
  

!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_pz(pz_VarA, pz_VarB, DelTime)
    !
    ! fs 格子点に対するスポンジ層
    !

    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(inout):: pz_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: pz_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: pz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !スポンジ層によるダンピングを計算  
    pz_Var  = pz_VarA - ( pz_DampRateH + pz_DampRateV ) * pz_VarB * DelTime
    pz_VarA = pz_Var
    
  end subroutine DampSponge_pz


!!!------------------------------------------------------------------------!!!
  function xz_Sponge( xz_Var )
    !
    ! ss 格子点に対するスポンジ層
    !

    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)             :: xz_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !スポンジ層によるダンピングを計算
    xz_Sponge = - ( xz_DampRateH + xz_DampRateV ) * xz_Var
    
  end function xz_Sponge
  

!!!------------------------------------------------------------------------!!!
  function xr_Sponge( xr_Var )
    !
    ! sf 格子点に対するスポンジ層
    !
    
    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)             :: xr_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: xr_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !スポンジ層によるダンピングを計算  
    xr_Sponge = - ( xr_DampRateH + xr_DampRateV ) * xr_Var
    
  end function xr_Sponge
  
!!!------------------------------------------------------------------------!!!
  function pz_Sponge( pz_Var )
    !
    ! fs 格子点に対するスポンジ層
    !

    !モジュール呼び出し
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x 方向の配列の下限
      &                    DimXMax,       &! x 方向の配列の上限
      &                    DimZMin,       &! z 方向の配列の下限
      &                    DimZMax         ! z 方向の配列の上限

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP)             :: pz_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: pz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !スポンジ層によるダンピングを計算  
    pz_Sponge = - ( pz_DampRateH + pz_DampRateV ) * pz_Var 
    
  end function pz_Sponge
  
end module Damping
