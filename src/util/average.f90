!= 平均演算モジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: average_3d.f90,v 1.5 2007-04-11 11:59:58 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2015. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module average
  !
  !平均操作を行うための関数を束ねたパッケージ型モジュール. 
  !水平 Arakawa-C, 鉛直 Lorentz グリッドとする. 
  !

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !関数を public にする
  public :: xyz_pyz, xyr_pyr, xqz_pqz, xqr_pqr
  public :: pyz_xyz, pyr_xyr, pqz_xqz, pqr_xqr
  public :: xyz_xqz, pyz_pqz, xyr_xqr, pyr_pqr
  public :: xqz_xyz, pqz_pyz, xqr_xyr, pqr_pyr
  public :: xyz_xyr, pyz_pyr, xqz_xqr, pqz_pqr
  public :: xyr_xyz, pyr_pyz, xqr_xqz, pqr_pqz
  public :: pqz_xyz, pyr_xyz, xqr_xyz
  public :: xyz_pqz, xyz_pyr, xyz_xqr

  !--------------------------------------------
  interface xyz_pyz
    module procedure xaa_paa
  end interface xyz_pyz

  interface xyr_pyr
    module procedure xaa_paa
  end interface xyr_pyr
  
  interface xqz_pqz
    module procedure xaa_paa
  end interface xqz_pqz

  interface xqr_pqr
    module procedure xaa_paa
  end interface xqr_pqr
  !--------------------------------------------
  interface pyz_xyz
    module procedure paa_xaa
  end interface pyz_xyz

  interface pqz_xqz
    module procedure paa_xaa
  end interface pqz_xqz

  interface pyr_xyr
    module procedure paa_xaa
  end interface pyr_xyr

  interface pqr_xqr
    module procedure paa_xaa
  end interface pqr_xqr
  !--------------------------------------------
  interface xyz_xqz
    module procedure aya_aqa
  end interface xyz_xqz

  interface pyz_pqz
    module procedure aya_aqa
  end interface pyz_pqz

  interface xyr_xqr
    module procedure aya_aqa
  end interface xyr_xqr

  interface pyr_pqr
    module procedure aya_aqa
  end interface pyr_pqr
  !--------------------------------------------
  interface xqz_xyz
    module procedure aqa_aya
  end interface xqz_xyz
  
  interface pqz_pyz
    module procedure aqa_aya
  end interface pqz_pyz
  
  interface xqr_xyr
    module procedure aqa_aya
  end interface xqr_xyr

  interface pqr_pyr
    module procedure aqa_aya
  end interface pqr_pyr
  !--------------------------------------------
  interface xyz_xyr
    module procedure aaz_aar
  end interface xyz_xyr

  interface pyz_pyr
    module procedure aaz_aar
  end interface pyz_pyr

  interface xqz_xqr
    module procedure aaz_aar
  end interface xqz_xqr

  interface pqz_pqr
    module procedure aaz_aar
  end interface pqz_pqr
  !--------------------------------------------
  interface xyr_xyz
    module procedure aar_aaz
  end interface xyr_xyz

  interface pyr_pyz
    module procedure aar_aaz
  end interface pyr_pyz

  interface xqr_xqz
    module procedure aar_aaz
  end interface xqr_xqz
  
  interface pqr_pqz
    module procedure aar_aaz
  end interface pqr_pqz
  !--------------------------------------------
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xaa_paa(paa_Var)
    !
    ! x 方向: 半整数格子から整数格子への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in) :: paa_Var(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: xaa_paa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix
    
    
    do ix = imin+1, imax
      xaa_paa(ix,:,:) = ( paa_Var(ix,:,:) + paa_Var(ix-1,:,:) ) * 0.5d0 
    end do
    
    ! imin 格子上の値
    xaa_paa(imin,:,:) = xaa_paa(imin+1,:,:)
    
  end function xaa_paa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function paa_xaa(xaa_Var)
    !
    ! x 方向: 整数格子から半整数格子への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: xaa_Var(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: paa_xaa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: ix
    

    do ix = imin, imax-1
      paa_xaa(ix,:,:) = ( xaa_Var(ix+1,:,:) + xaa_Var(ix,:,:) ) * 0.5d0
    end do
    
    paa_xaa(imax,:,:) = paa_xaa(imax-1,:,:) 
    
  end function paa_xaa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aya_aqa(aqa_Var)
    !
    ! y 方向: 半整数格子から整数格子への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: aqa_Var(imin:imax,jmin:jmax,kmin:kmax)   
    real(DP)            :: aya_aqa(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: jy
    

    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      aya_aqa = aqa_Var
      return 
    end if
    
    do jy = jmin+1, jmax
      aya_aqa(:,jy,:) = ( aqa_Var(:,jy,:) + aqa_Var(:,jy-1,:) ) * 0.5d0 
    end do
    
    aya_aqa(:,jmin,:) = aya_aqa(:,jmin+1,:)
    
  end function aya_aqa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aqa_aya(aya_Var)
    !
    ! y 方向: 整数格子から半整数格子への変換
    !

    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP),intent(in) :: aya_Var(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)            :: aqa_aya(imin:imax,jmin:jmax,kmin:kmax) 
    integer             :: jy
    
    
    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      aqa_aya = aya_Var
      return 
    end if
    
    do jy = jmin, jmax-1
      aqa_aya(:,jy,:) = ( aya_Var(:,jy+1,:) + aya_Var(:,jy,:) ) * 0.5d0
    end do
    
    aqa_aya(:,jmax,:) =  aqa_aya(:,jmax-1,:) 
    
  end function aqa_aya

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function aaz_aar(aar_Var)
    !
    ! z 方向: 半整数格子から整数格子への変換
    !
  
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: aar_Var(imin:imax,jmin:jmax,kmin:kmax)   
    real(DP)            :: aaz_aar(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: kz
    
    
    do kz = kmin+1, kmax
      aaz_aar(:,:,kz) = ( aar_Var(:,:,kz) + aar_Var(:,:,kz-1) ) * 0.5d0 
    end do
    
    aaz_aar(:,:,kmin) = aaz_aar(:,:,kmin+1)
    
  end function aaz_aar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  function aar_aaz(aaz_Var)
    !
    ! z 方向: 整数格子から半整数格子への変換
    !

    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: aaz_Var(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP)            :: aar_aaz(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: kz


    do kz = kmin, kmax-1
      aar_aaz(:,:,kz) = ( aaz_Var(:,:,kz+1) + aaz_Var(:,:,kz) ) * 0.5d0
    end do
    
    aar_aaz(:,:,kmax) = aar_aaz(:,:,kmax-1)
    
  end function aar_aaz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pqz_xyz(xyz_Var)
    !
    ! xy 方向: 半整数格子 (pz) から整数格子 (xy) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: xyz_Var(imin:imax,jmin:jmax,kmin:kmax)! 入力
    real(DP)             :: pqz_xyz(imin:imax,jmin:jmax,kmin:kmax)! 出力
    integer              :: ix, jy, kz                            ! ループ添字
    
    ! 初期化
    pqz_xyz = 0.0d0
    
    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      do ix = imin, imax-1
        pqz_xyz(ix,:,:) =            &
          &  (                       &
          &    + xyz_Var(ix+1,:,:)   &
          &    + xyz_Var(ix,  :,:)   &
          &  ) * 0.5d0 
      end do
      return 
    end if
    
    do kz = kmin, kmax
      do jy = jmin, jmax-1
        do ix = imin, imax-1
          pqz_xyz(ix,jy,kz) =              &
            &  (                           &
            &    + xyz_Var(ix+1, jy+1, kz) &
            &    + xyz_Var(ix+1, jy,   kz) &
            &    + xyz_Var(ix,   jy+1, kz) &
            &    + xyz_Var(ix,   jy,   kz) &
            &  ) * 0.25d0 
        end do
      end do
    end do
    
  end function pqz_xyz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pyr_xyz(xyz_Var)
    !
    ! xz 方向: 半整数格子 (pr) から整数格子 (xz) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP), intent(in) :: xyz_Var(imin:imax,jmin:jmax,kmin:kmax)! 入力
    real(DP)             :: pyr_xyz(imin:imax,jmin:jmax,kmin:kmax)! 出力
    integer              :: ix, jy, kz                            ! ループ添字
    
    ! 初期化
    pyr_xyz = 0.0d0
    
    do kz = kmin, kmax-1
      do jy = jmin, jmax
        do ix = imin, imax-1
          pyr_xyz(ix,jy,kz) =               &
            &  (                            &
            &    + xyz_Var(ix+1, jy, kz+1)  &
            &    + xyz_Var(ix+1, jy, kz)    &
            &    + xyz_Var(ix,   jy, kz+1)  &
            &    + xyz_Var(ix,   jy, kz)    &
            &  ) * 0.25d0 
        end do
      end do
    end do
    
  end function pyr_xyz
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  function xqr_xyz(xyz_Var)
    !
    ! yz 方向: 半整数格子 (qr) から整数格子 (yz) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義    
    real(DP), intent(in) :: xyz_Var(imin:imax,jmin:jmax,kmin:kmax)! 入力
    real(DP)             :: xqr_xyz(imin:imax,jmin:jmax,kmin:kmax)! 出力
    integer              :: ix, jy, kz                            ! ループ添字
    
    ! 初期化
    xqr_xyz = 0.0d0
    
    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      do kz = kmin, kmax-1
        xqr_xyz(:,:,kz) =              &
          &  (                         &
          &    + xyz_Var(:,:,kz+1)     &
          &    + xyz_Var(:,:,kz)       &
          &  ) * 0.5d0 
      end do
      return 
    end if
    
    do kz = kmin, kmax-1
      do jy = jmin, jmax-1
        do ix = imin, imax
          xqr_xyz(ix,jy,kz) =              &
            &  (                           &
            &    + xyz_Var(ix, jy+1, kz+1) &
            &    + xyz_Var(ix, jy+1, kz)   &
            &    + xyz_Var(ix, jy,   kz+1) &
            &    + xyz_Var(ix, jy,   kz)   &
            &  ) * 0.25d0 
        end do
      end do
    end do

  end function xqr_xyz
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  function xyz_pqz(pqz_Var)
    !
    ! xy 方向: 整数格子 (xy) から整数格子 (pq) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(in) :: pqz_Var(imin:imax,jmin:jmax,kmin:kmax) ! 入力
    real(DP)             :: xyz_pqz(imin:imax,jmin:jmax,kmin:kmax) ! 出力
    integer              :: ix, jy, kz                             ! ループ添字
    
    ! 初期化
    xyz_pqz = 0.0d0
    
    ! y 方向に格子点が 1 の場合の処理
    !
    if (.NOT. FlagCalc3D) then 
      do ix = imin+1, imax
        xyz_pqz(ix,:,:) =             &
          &  (                        &
          &    + pqz_Var(ix,  :, :)   &
          &    + pqz_Var(ix-1,:, :)   &
          &  ) * 0.5d0 
      end do
      return 
    end if
    
    do kz = kmin, kmax
      do jy = jmin+1, jmax
        do ix = imin+1, imax
          xyz_pqz(ix,jy,kz) =              &
            &  (                           &
            &    + pqz_Var(ix,   jy,   kz) &
            &    + pqz_Var(ix,   jy-1, kz) &
            &    + pqz_Var(ix-1, jy,   kz) &
            &    + pqz_Var(ix-1, jy-1, kz) &
            &  ) * 0.25d0 
        end do
      end do
    end do

  end function xyz_pqz

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xyz_pyr(pyr_Var)
    !
    ! xz 方向: 整数格子 (xz) から整数格子 (pr) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax
       
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: pyr_Var(imin:imax,jmin:jmax,kmin:kmax) ! 入力
    real(DP)             :: xyz_pyr(imin:imax,jmin:jmax,kmin:kmax) ! 出力
    integer              :: ix, jy, kz                             ! ループ添字
    
    ! 初期化
    xyz_pyr = 0.0d0
    
    do kz = kmin+1, kmax
      do jy = jmin, jmax
        do ix = imin+1, imax
          xyz_pyr(ix,jy,kz) =                &
            &  (                             &
            &    + pyr_Var(ix,   jy, kz)     &
            &    + pyr_Var(ix,   jy, kz-1)   &
            &    + pyr_Var(ix-1, jy, kz)     &
            &    + pyr_Var(ix-1, jy, kz-1)   &
            &  ) * 0.25d0 
        end do
      end do
    end do

  end function xyz_pyr

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xyz_xqr(xqr_Var)
    !
    ! yz 方向: 整数格子 (yz) から整数格子 (qr) への変換
    !
    
    !モジュール呼び出し
    use dc_types,  only : DP
    use gridset,   only : imin, imax, &
      &                   jmin, jmax, &
      &                   kmin, kmax, FlagCalc3D
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xqr_Var(imin:imax,jmin:jmax,kmin:kmax) ! 入力
    real(DP)             :: xyz_xqr(imin:imax,jmin:jmax,kmin:kmax) ! 出力
    integer              :: ix, jy, kz                             ! ループ添字
    
    ! 初期化
    xyz_xqr = 0.0d0
    
    if (.NOT. FlagCalc3D) then 
      do kz = kmin+1, kmax
        xyz_xqr(:,:,kz) =             &
          &  (                        &
          &    + xqr_Var(:, :, kz)    &
          &    + xqr_Var(:, :, kz-1)  &
          &  ) * 0.5d0 
      end do
      return 
    end if
    
    do kz = kmin+1, kmax
      do jy = jmin+1, jmax
        do ix = imin, imax
          xyz_xqr(ix,jy,kz) =              &
            &  (                           &
            &    + xqr_Var(ix, jy,   kz)   &
            &    + xqr_Var(ix, jy,   kz-1) &
            &    + xqr_Var(ix, jy-1, kz)   &
            &    + xqr_Var(ix, jy-1, kz-1) &
            &  ) * 0.25d0 
        end do
      end do
    end do

  end function xyz_xqr


end module average
