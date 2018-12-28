!= 微分演算モジュール (4 次精度中央差分)
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: differentiate_center4.f90,v 1.3 2006-09-21 03:01:02 odakker Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module differentiate_center4
  !
  !4 次精度の微分演算を行うための関数群をまとめたパッケージ型モジュール
  !水平 Arakawa-C, 鉛直 Lorentz グリッドでの微分計算
  !
  ! 差分式は以下の通り
  !
  !    du/dx(i) = 9/8*[u(i)-u(i-1)]/dx - 1/24*[u(i+1)-u(i-2)]/dx
  ! 

  !暗黙の型宣言禁止
  implicit none

  !属性指定
  private

  public xyz_dx_pyz, xqz_dx_pqz, xyr_dx_pyr
  public pyz_dx_xyz, pqz_dx_xqz, pyr_dx_xyr
  public xyz_dy_xqz, pyz_dy_pqz, xyr_dy_xqr
  public xqz_dy_xyz, pqz_dy_pyz, xqr_dy_xyr
  public xyz_dz_xyr, pyz_dz_pyr, xqz_dz_xqr
  public xyr_dz_xyz, pyr_dz_pyz, xqr_dz_xqz

  !----------------------------------------------
  interface xyz_dx_pyz
    module procedure xaa_dx_paa
  end interface xyz_dx_pyz
  
  interface xqz_dx_pqz
    module procedure xaa_dx_paa
  end interface xqz_dx_pqz
  
  interface xyr_dx_pyr
    module procedure xaa_dx_paa
  end interface xyr_dx_pyr
  !----------------------------------------------
  interface pyz_dx_xyz
    module procedure paa_dx_xaa
  end interface pyz_dx_xyz
  
  interface pqz_dx_xqz
    module procedure paa_dx_xaa
  end interface pqz_dx_xqz
  
  interface pyr_dx_xyr
    module procedure paa_dx_xaa
  end interface pyr_dx_xyr
  !----------------------------------------------
  interface xyz_dy_xqz
    module procedure aya_dy_aqa
  end interface xyz_dy_xqz
  
  interface pyz_dy_pqz
    module procedure aya_dy_aqa
  end interface pyz_dy_pqz
  
  interface xyr_dy_xqr
    module procedure aya_dy_aqa
  end interface xyr_dy_xqr
  !----------------------------------------------
  interface xqz_dy_xyz
    module procedure aqa_dy_aya
  end interface xqz_dy_xyz

  interface pqz_dy_pyz
    module procedure aqa_dy_aya
  end interface pqz_dy_pyz

  interface xqr_dy_xyr
    module procedure aqa_dy_aya
  end interface xqr_dy_xyr
  !----------------------------------------------
  interface xyz_dz_xyr
    module procedure aaz_dz_aar
  end interface xyz_dz_xyr
  
  interface pyz_dz_pyr
    module procedure aaz_dz_aar
  end interface pyz_dz_pyr

  interface xqz_dz_xqr
    module procedure aaz_dz_aar
  end interface xqz_dz_xqr
  !----------------------------------------------
  interface xyr_dz_xyz
    module procedure aar_dz_aaz
  end interface xyr_dz_xyz
  
  interface pyr_dz_pyz
    module procedure aar_dz_aaz
  end interface pyr_dz_pyz

  interface xqr_dz_xqz
    module procedure aar_dz_aaz
  end interface xqr_dz_xqz
  !----------------------------------------------
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xaa_dx_paa(paa_Var)
    !
    ! x 微分: 半整数格子点上の微分
    !

    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dx    !格子間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義    
    real(DP),intent(in) :: paa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: xaa_dx_paa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix                          
    
    ! 1 階微分の計算
    !
    do ix = imin+2, imax-1
      xaa_dx_paa(ix,:,:) = &
        &   ( (paa_Var(ix,:,:)   - paa_Var(ix-1,:,:) ) * 9.0d0 / ( 8.0d0 * dx ) ) &
        & - ( (paa_Var(ix+1,:,:) - paa_Var(ix-2,:,:) ) / ( 24.0d0 * dx ) )
    end do

    xaa_dx_paa(imin,  :,:) = xaa_dx_paa(imin+2,:,:)
    xaa_dx_paa(imin+1,:,:) = xaa_dx_paa(imin+2,:,:)
    xaa_dx_paa(imax,  :,:) = xaa_dx_paa(imax-1,:,:)
    
  end function xaa_dx_paa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function paa_dx_xaa(xaa_Var)
    !
    ! x 微分: 整数格子点上の微分
    !
    
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dx    !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP),intent(in) :: xaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: paa_dx_xaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix                          

    ! 初期化
    !
    paa_dx_xaa = 0.0d0
    
    ! 1 階微分の計算
    !
    do ix = imin+1, imax-2
      paa_dx_xaa(ix,:,:) = &
        &   ( (xaa_Var(ix+1,:,:) - xaa_Var(ix,:,:)   ) * 9.0d0 / ( 8.0d0 * dx ) ) &
        & - ( (xaa_Var(ix+2,:,:) - xaa_Var(ix-1,:,:) ) / ( 24.0d0 * dx ) )
    end do

    paa_dx_xaa(imin,  :,:) = paa_dx_xaa(imin+1,:,:)
    paa_dx_xaa(imax-1,:,:) = paa_dx_xaa(imax-2,:,:)
    paa_dx_xaa(imax,  :,:) = paa_dx_xaa(imax-2,:,:)
    
  end function paa_dx_xaa
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aya_dy_aqa(aqz_Var)
    !
    ! y 微分: 半整数格子上の微分
    !
    
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax, &
                                            !領域サイズ
      &                         FlagCalc3D  !3 次元計算フラグ
    use axesset,         only : dy          !格子間隔
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP),intent(in) :: aqz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aya_dy_aqa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: jy

    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      aya_dy_aqa = 0.0d0
      return
    end if

    ! 1 階微分の計算
    !
    do jy = jmin+2, jmax-1
      aya_dy_aqa(:,jy,:) = &
        &   ( ( aqz_Var(:,jy,:)   - aqz_Var(:,jy-1,:) ) * 9.0d0 / ( 8.0d0 * dy ) ) &
        & - ( ( aqz_Var(:,jy+1,:) - aqz_Var(:,jy-2,:) ) / ( 24.0d0 * dy ) )
    end do

    aya_dy_aqa(:, jmin,  :) = aya_dy_aqa(:, jmin+2, :)
    aya_dy_aqa(:, jmin+1,:) = aya_dy_aqa(:, jmin+2, :)
    aya_dy_aqa(:, jmax,  :) = aya_dy_aqa(:, jmax-1, :)
    
  end function aya_dy_aqa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aqa_dy_aya(aya_Var)
    !
    ! y 微分: 整数格子上の微分
    !
    
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax, &
                                            !領域サイズ
      &                         FlagCalc3D  !3 次元計算フラグ
    use axesset,         only : dy          !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP),intent(in) :: aya_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aqa_dy_aya(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: jy

    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      aqa_dy_aya = 0.0d0
      return 
    end if

    ! 1 階微分の計算
    !
    do jy = jmin+1, jmax-2
      aqa_dy_aya(:,jy,:) = &
        &   ( ( aya_Var(:,jy+1,:) - aya_Var(:,jy,:)   ) * 9.0d0 / ( 8.0d0 * dy ) ) &
        & - ( ( aya_Var(:,jy+2,:) - aya_Var(:,jy-1,:) ) / ( 24.0d0 * dy ) )
    end do

    aqa_dy_aya(:, jmin,  :) = aqa_dy_aya(:, jmin+1, :)
    aqa_dy_aya(:, jmax-1,:) = aqa_dy_aya(:, jmax-2, :)
    aqa_dy_aya(:, jmax,  :) = aqa_dy_aya(:, jmax-2, :)
    
  end function aqa_dy_aya

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aaz_dz_aar(aar_Var)
    !
    ! z 微分: 半整数格子上の微分
    !
  
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dz
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義    
    real(DP),intent(in) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aaz_dz_aar(imin:imax,jmin:jmax,kmin:kmax)
                                 !領域サイズ
    integer             :: kz    !格子間隔
    
    ! 1 階微分の計算
    !
    do kz = kmin+2, kmax-1
      aaz_dz_aar(:,:,kz) = &
        &   ( ( aar_Var(:,:,kz)   - aar_Var(:,:,kz-1) ) * 9.0d0 / ( 8.0d0 * dz ) ) &
        & - ( ( aar_Var(:,:,kz+1) - aar_Var(:,:,kz-2) ) / ( 24.0d0 * dz ) )
    end do
    
    aaz_dz_aar(:,:,kmin  ) = aaz_dz_aar(:,:,kmin+2)
    aaz_dz_aar(:,:,kmin+1) = aaz_dz_aar(:,:,kmin+2)
    aaz_dz_aar(:,:,kmax  ) = aaz_dz_aar(:,:,kmax-1)

  end function aaz_dz_aar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aar_dz_aaz(aaz_Var)
    !
    ! z 微分: 整数格子上の微分
    !
    
    !モジュール呼び出し
    use dc_types,        only : DP
    use gridset,         only : imin, imax, jmin, jmax, kmin, kmax
                                      !領域サイズ
    use axesset,         only : dz    !格子間隔

    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP),intent(in) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aar_dz_aaz(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: kz
    
    ! 1 階微分の計算
    !
    do kz = kmin+1, kmax-2
      aar_dz_aaz(:,:,kz) = &
        &   ( ( aaz_Var(:,:,kz+1) - aaz_Var(:,:,kz)   ) * 9.0d0 / ( 8.0d0 * dz ) ) &
        & - ( ( aaz_Var(:,:,kz+2) - aaz_Var(:,:,kz-1) ) / ( 24.0d0 * dz ) )
    end do

    aar_dz_aaz(:,:,kmin  ) = aar_dz_aaz(:,:,kmin+1)
    aar_dz_aaz(:,:,kmax-1) = aar_dz_aaz(:,:,kmax-2)
    aar_dz_aaz(:,:,kmax  ) = aar_dz_aaz(:,:,kmax-2)
    
  end function aar_dz_aaz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module differentiate_center4
