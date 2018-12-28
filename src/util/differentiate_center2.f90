!= 微分演算モジュール (2 次精度中央差分)
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: differentiate_center2.f90,v 1.4 2007-04-11 11:59:58 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module differentiate_center2
  !
  ! 2 次精度の微分演算を行うための関数群をまとめたパッケージ型モジュール
  ! 水平 Arakawa-C, 鉛直 Lorentz グリッドでの微分計算
  !

  !暗黙の型宣言禁止
  implicit none

  !private 属性
  private
 
  !公開要素
  public xyz_dx_pyz, xqz_dx_pqz, xyr_dx_pyr
  public pyz_dx_xyz, pqz_dx_xqz, pyr_dx_xyr
  public xyz_dy_xqz, pyz_dy_pqz, xyr_dy_xqr
  public xqz_dy_xyz, pqz_dy_pyz, xqr_dy_xyr
  public xyz_dz_xyr, pyz_dz_pyr, xqz_dz_xqr
  public xyr_dz_xyz, pyr_dz_pyz, xqr_dz_xqz

  public xyz_dx4_xyz, pyz_dx4_pyz, xqz_dx4_xqz, xyr_dx4_xyr
  public xyz_dy4_xyz, pyz_dy4_pyz, xqz_dy4_xqz, xyr_dy4_xyr
  public xyz_dz4_xyz, pyz_dz4_pyz, xqz_dz4_xqz, xyr_dz4_xyr

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
  interface xyz_dx4_xyz
    module procedure aaa_dx4_aaa
  end interface xyz_dx4_xyz
  
  interface pyz_dx4_pyz
    module procedure aaa_dx4_aaa
  end interface pyz_dx4_pyz
  
  interface xqz_dx4_xqz
    module procedure aaa_dx4_aaa
  end interface xqz_dx4_xqz
  
  interface xyr_dx4_xyr
    module procedure aaa_dx4_aaa
  end interface xyr_dx4_xyr
  !----------------------------------------------
  interface xyz_dy4_xyz
    module procedure aaa_dy4_aaa
  end interface xyz_dy4_xyz
  
  interface pyz_dy4_pyz
    module procedure aaa_dy4_aaa
  end interface pyz_dy4_pyz
  
  interface xqz_dy4_xqz
    module procedure aaa_dy4_aaa
  end interface xqz_dy4_xqz

  interface xyr_dy4_xyr
    module procedure aaa_dy4_aaa
  end interface xyr_dy4_xyr
  !----------------------------------------------
  interface xyz_dz4_xyz
    module procedure aaa_dz4_aaa
  end interface xyz_dz4_xyz
  
  interface pyz_dz4_pyz
    module procedure aaa_dz4_aaa
  end interface pyz_dz4_pyz
  
  interface xqz_dz4_xqz
    module procedure aaa_dz4_aaa
  end interface xqz_dz4_xqz
  
  interface xyr_dz4_xyr
    module procedure aaa_dz4_aaa
  end interface xyr_dz4_xyr
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xaa_dx_paa(paa_Var)
    !
    ! x 微分: 半整数格子点上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dx    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: paa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: xaa_dx_paa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix
    
    ! 1 階微分の計算
    !
    do ix = imin+1, imax
      xaa_dx_paa(ix,:,:) = ( paa_Var(ix,:,:) - paa_Var(ix-1,:,:) ) / dx
    end do
    
    xaa_dx_paa(imin,:,:) = xaa_dx_paa(imin+1,:,:)
    
  end function xaa_dx_paa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function paa_dx_xaa(xaa_Var)
    !
    ! x 微分: 整数格子点上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dx    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: xaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: paa_dx_xaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix 
    
    ! 1 階微分の計算
    !
    do ix = imin, imax-1
      paa_dx_xaa(ix,:,:) = ( xaa_Var(ix+1,:,:) - xaa_Var(ix,:,:) ) / dx
    end do
    
    paa_dx_xaa(imax,:,:) = paa_dx_xaa(imax-1,:,:)
    
  end function paa_dx_xaa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aya_dy_aqa(aqa_Var)
    !
    ! y 微分: 半整数格子上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                      !配列サイズ
      &                   FlagCalc3D  !3次元計算か否かのフラグ
    use axesset,   only : dy          !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aqa_Var(imin:imax,jmin:jmax,kmin:kmax) 
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
    do jy = jmin+1, jmax
      aya_dy_aqa(:,jy,:) = ( aqa_Var(:,jy,:) - aqa_Var(:,jy-1,:) ) / dy
    end do
    
    aya_dy_aqa(:,jmin,:) = aya_dy_aqa(:,jmin+1,:)
    
  end function aya_dy_aqa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aqa_dy_aya(aya_Var)
    !
    ! y 微分: 整数格子上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                      !配列サイズ
      &                   FlagCalc3D  !3次元計算か否かのフラグ
    use axesset,   only : dy          !格子点間隔
    
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
    do jy = jmin, jmax-1
      aqa_dy_aya(:,jy,:) = ( aya_Var(:,jy+1,:) - aya_Var(:,jy,:) ) / dy
    end do
    
    aqa_dy_aya(:,jmax,:) = aqa_dy_aya(:,jmax-1,:)
    
  end function aqa_dy_aya

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aaz_dz_aar(aar_Var)
    !
    ! z 微分: 半整数格子上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dz    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aaz_dz_aar(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: kz
    
    ! 1 階微分の計算
    !
    do kz = kmin+1, kmax
      aaz_dz_aar(:,:,kz) = ( aar_Var(:,:,kz) - aar_Var(:,:,kz-1) ) / dz
    end do
    
    aaz_dz_aar(:,:,kmin) = aaz_dz_aar(:,:,kmin+1)
    
  end function aaz_dz_aar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aar_dz_aaz(aaz_Var)
    !
    ! z 微分: 整数格子上の微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dz    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aar_dz_aaz(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: kz
    
    ! 1 階微分の計算
    !
    do kz = kmin, kmax-1
      aar_dz_aaz(:,:,kz) = ( aaz_Var(:,:,kz+1) - aaz_Var(:,:,kz) ) / dz
    end do
    
    aar_dz_aaz(:,:,kmax) = aar_dz_aaz(:,:,kmax-1)
    
  end function aar_dz_aaz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  function aaa_dx4_aaa( aaa_Var )
    !
    ! x 方向の 4 階微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dx    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aaa_dx4_aaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin + 2, imax - 2            

          aaa_dx4_aaa(i,j,k) =                   &
            &   (                                &
            &       aaa_Var(i+2,j,k)             &
            &     + aaa_Var(i-2,j,k)             &
            &     - aaa_Var(i+1,j,k) * 4.0d0     &
            &     - aaa_Var(i-1,j,k) * 4.0d0     &
            &     + aaa_Var(i  ,j,k) * 6.0d0     &
            &   ) / ( dx ** 4.0d0 )

        end do
      end do
    end do

    aaa_dx4_aaa(imin:imin+1,:,:) = 1.0d10
    aaa_dx4_aaa(imax-1:imax,:,:) = 1.0d10
    
  end function aaa_dx4_aaa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aaa_dy4_aaa( aaa_Var )
    !
    ! y 方向の 4 階微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax, &
                                      !配列サイズ
      &                   FlagCalc3D  !3次元計算か否かのフラグ
    use axesset,   only : dy    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aaa_dy4_aaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: i, j, k

    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      aaa_dy4_aaa = 0.0d0
      return 
    end if


    do k = kmin, kmax
      do j = jmin + 2, jmax - 2
        do i = imin, imax

          aaa_dy4_aaa(i,j,k) =                 &
            &   (                              &
            &       aaa_Var(i,j+2,k)           &
            &     + aaa_Var(i,j-2,k)           &
            &     - aaa_Var(i,j+1,k) * 4.0d0   &
            &     - aaa_Var(i,j-1,k) * 4.0d0   &
            &     + aaa_Var(i,j  ,k) * 6.0d0   &
            &   ) / ( dy ** 4.0d0 ) 

        end do
      end do
    end do
   
    aaa_dy4_aaa(:,jmin:jmin+1,:) = 1.0d10
    aaa_dy4_aaa(:,jmax-1:jmax,:) = 1.0d10

  end function aaa_dy4_aaa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aaa_dz4_aaa( aaa_Var )
    !
    ! 4 階微分
    !
    
    !モジュール読み込み
    use dc_types,  only : DP
    use gridset,   only : imin, imax, jmin, jmax, kmin, kmax
                                !配列サイズ
    use axesset,   only : dz    !格子点間隔
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: aaa_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aaa_dz4_aaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: i, j, k

    do k = kmin + 2, kmax - 2
      do j = jmin, jmax
        do i = imin, imax

          aaa_dz4_aaa(i,j,k) =                 &
            &   (                              &
            &       aaa_Var(i,j,k+2)           &
            &     + aaa_Var(i,j,k-2)           &
            &     - aaa_Var(i,j,k+1) * 4.0d0   &
            &     - aaa_Var(i,j,k-1) * 4.0d0   &
            &     + aaa_Var(i,j,k  ) * 6.0d0   &
            &   ) / ( dz ** 4.0d0 )
          
        end do
      end do
    end do
    
    aaa_dz4_aaa(:,:,kmin:kmin+1) = 1.0d10
    aaa_dz4_aaa(:,:,kmax-1:kmax) = 1.0d10
    
  end function aaa_dz4_aaa

end module differentiate_center2
