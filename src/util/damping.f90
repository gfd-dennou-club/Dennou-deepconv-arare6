!= スポンジ層モジュール
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: damping.f90,v 1.14 2014/11/07 06:46:45 sugiyama Exp $
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
  public SpongeLayer_forcing
  public SpongeLayer_MeanFlow
  
  !変数定義
  real(DP), save, private              :: EFTime = 5000.0d0!スポンジ層の e-folding time
  real(DP), save, private              :: DepthH = 0.0d0   !スポンジ層の厚さ(水平方向)
  real(DP), save, private              :: DepthV = 0.0d0   !スポンジ層の厚さ(鉛直方向) [上部境界]
  real(DP), save, private              :: DepthVb= 0.0d0   !スポンジ層の厚さ(鉛直方向) [下部境界]
  real(DP), allocatable, save, private :: xyz_Gamma(:,:,:) !xyz 格子減衰係数(水平方向)
  real(DP), allocatable, save, private :: pyz_Gamma(:,:,:) !pyz 格子減衰係数(鉛直方向)
  real(DP), allocatable, save, private :: xqz_Gamma(:,:,:) !xqz 格子減衰係数(鉛直方向)
  real(DP), allocatable, save, private :: xyr_Gamma(:,:,:) !xyr 格子減衰係数(鉛直方向)
  real(DP), save, private              :: FactorSpngVelX  = 1.0d0
  real(DP), save, private              :: FactorSpngVelY  = 1.0d0
  real(DP), save, private              :: FactorSpngVelZ  = 1.0d0
  real(DP), save, private              :: FactorSpngPTemp = 1.0d0
  real(DP), save, private              :: FactorSpngExner = 0.0d0
  real(DP), save, public               :: DampSound  = 0.0d0   !音波減衰項の減衰係数

  real(DP), allocatable, save, private :: pyz_SpngMF(:,:,:)
  real(DP), allocatable, save, private :: xqz_SpngMF(:,:,:)
  real(DP), allocatable, save, private :: z_VelX0(:)
  real(DP), allocatable, save, private :: z_VelY0(:)

  real(DP), save, private              :: EFTimeMF = 50000.0d0
  real(DP), save, private              :: DelTimeMF = 1.0d0
  real(DP), save, private              :: ZMinMF  = 0.0d0
  real(DP), save, private              :: ZMaxMF  = 1.0d10

contains 
  
!!!------------------------------------------------------------------------!!!
  subroutine Damping_Init
    !
    ! 音波減衰項とスポンジ層の減衰係数の初期化
    ! 
    use dc_types,   only: DP, STRING
    use dc_iounit,  only: FileOpen
    use dc_message, only: MessageNotify
    use gtool_historyauto, only: HistoryAutoAddVariable
    use gridset, only: imin,       &! x 方向の配列の下限
      &                imax,       &! x 方向の配列の上限
      &                jmin,       &! y 方向の配列の下限
      &                jmax,       &! y 方向の配列の上限
      &                kmin,       &! z 方向の配列の下限
      &                kmax         ! z 方向の配列の上限
    use axesset, only: x_X,        &!X 座標軸(スカラー格子点)
      &                y_Y,        &!Y 座標軸(スカラー格子点)
      &                z_Z,        &!Z 座標軸(スカラー格子点)
      &                p_X,        &!X 座標軸(フラックス格子点)
      &                q_Y,        &!Y 座標軸(フラックス格子点)
      &                r_Z,        &!Z 座標軸(フラックス格子点)
      &                dx, dy, dz, &! 格子間隔
      &                XMax,       &!X 座標の最大値
      &                YMax,       &!Y 座標の最大値
      &                ZMin,       &!Z 座標の最小値 
      &                ZMax         !Z 座標の最大値 
    use namelist_util,  &
      &          only: namelist_filename
    use timeset, only: DelTimeShort    !短い時間ステップ
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), parameter  :: Pi =3.1415926535897932385d0   !円周率
    real(DP)             :: Alpha = 5.0d-2
    integer              :: unit ! 装置番号
    integer              :: i, j, k
    character(STRING)    :: InputFileMF = ""

    !NAMELIST から取得
    NAMELIST /damping_nml/ Alpha,                        &
      & EFTime, DepthH, DepthV, DepthVb,                 &
      & FactorSpngVelX, FactorSpngVelY, FactorSpngVelZ,  &
      & FactorSpngPTemp, FactorSpngExner,                &
      & DelTimeMF, ZMinMF, ZMaxMF, EFTimeMF, InputFileMF

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=damping_nml)
    close(unit)

    ! 音波減衰項
    DampSound = Alpha * ( Min(dx, dz) ** 2.0d0 ) / DelTimeShort

    !初期化
    allocate( &
      & xyz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & pyz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & xqz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & xyr_Gamma(imin:imax,jmin:jmax,kmin:kmax)    )

    !デフォルト値はゼロ. 
    xyz_Gamma = 0.0d0
    pyz_Gamma = 0.0d0
    xqz_Gamma = 0.0d0
    xyr_Gamma = 0.0d0

    !-----------------------------------------------------------------    
    ! 厚さのチェック
    ! DepthX > 0 の場合は, 厚さを dx, dy, dz より大きくする必要がある. 
    !
    if ( DepthH > 0.0d0 ) then 
      if ( DepthH < dx ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelX is %f", d=(/dx/))        
      else if ( DepthH < dx ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelX is %f", d=(/dx/))
      end if
      
      if ( DepthH < dy ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelY is %f", d=(/dy/))
      else if ( DepthH < dy ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelY is %f", d=(/dy/))
      end if
    end if

    if ( DepthV > 0.0d0 ) then 
      if ( DepthV < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthV is too thin. DelZ is %f", d=(/dz/) )      
      elseif ( DepthV < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthV is too thin. DelZ is %f", d=(/dz/) )      
      end if
    end if

    if ( DepthVb > 0.0d0 ) then 
      if ( DepthVb < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthVb is too thin. DelZ is %f", d=(/dz/) )      
      elseif ( DepthVb < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthVb is too thin. DelZ is %f", d=(/dz/) )      
      end if
    end if

    !-----------------------------------------------------------------    
    ! スポンジ層の減衰率
    ! スポンジ層の厚さ (DepthX) がゼロの場合はダンピングの係数をいじらない.
    !

    !水平方向の東側・西側境界
    !
    if ( DepthH > 0.0d0 ) then 
      do i = imin, imax
        !スカラー格子点の西側境界
        if ( x_X(i) < DepthH) then 
          xyz_Gamma(i,:,:) = xyz_Gamma(i,:,:) &
            & + ((1.0d0 - x_X(i) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !フラックス格子点の西側境界
        if ( p_X(i) < DepthH) then 
          pyz_Gamma(i,:,:) = pyz_Gamma(i,:,:) &
            & + ((1.0d0 - p_X(i) / DepthH) ** 3.0d0) / EFTime
        end if
        
      !スカラー格子点の東側境界    
        if ( x_X(i) > ( XMax - DepthH ) ) then 
          xyz_Gamma(i,:,:) = xyz_Gamma(i,:,:) &
            & + ((1.0d0 - (XMax - x_X(i)) / DepthH) ** 3.0d0) / EFTime 
        end if
        
        !フラックス格子点の東側境界    
        if ( p_X(i) > ( XMax - DepthH ) ) then 
          pyz_Gamma(i,:,:) = pyz_Gamma(i,:,:) &
            & + ((1.0d0 - (XMax - p_X(i)) / DepthH) ** 3.0d0) / EFTime 
        end if
      end do

      ! x 方向には同じ
      !
      xyr_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if
    
    !水平方向の南側・北側境界
    !
    if ( DepthH > 0.0d0 ) then 
      do j = jmin, jmax
        !スカラー格子点の南側境界
        if ( y_Y(j) < DepthH) then 
          xyz_Gamma(:,j,:) = xyz_Gamma(:,j,:) &
            & + ((1.0d0 - y_Y(j) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !フラックス格子点の南側境界
        if ( q_Y(j) < DepthH) then 
          xqz_Gamma(:,j,:) = xqz_Gamma(:,j,:) &
            & + ((1.0d0 - q_Y(j) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !スカラー格子点の北側境界    
        if ( y_Y(j) > ( YMax - DepthH ) ) then 
          xyz_Gamma(:,j,:) = xyz_Gamma(:,j,:) &
            & + ((1.0d0 - (YMax - y_Y(j)) / DepthH) ** 3.0d0) / EFTime 
        end if
        
        !フラックス格子点の北側境界    
        if ( q_Y(j) > ( YMax - DepthH ) ) then 
          xqz_Gamma(:,j,:) = xqz_Gamma(:,j,:)  &
            & + ((1.0d0 - (YMax - q_Y(j)) / DepthH) ** 3.0d0) / EFTime 
        end if
      end do

      ! y 方向には同じ
      !
      pyz_Gamma  = xyz_Gamma
      xyr_Gamma  = xyz_Gamma
    end if
    
    !鉛直方向の上部境界    
    !
    if ( DepthV > 0.0d0 ) then 
      do k = kmin, kmax
        !スカラー格子点
        if ( z_Z(k) >= ( ZMax - DepthV ) ) then 
          xyz_Gamma(:,:,k) = xyz_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * (z_Z(k) - ZMax + DepthV) / DepthV)) &
            &    / EFTime 
          
          ! 最大   z = ZMax
          !    gamma = 1.0 - cos( \pi ) = 1.0
          ! 最小   z = ZMax - DepthV
          !    gamma = 1.0 - cos( 0 ) = 0
          
        end if

        !フラックス格子点
        if ( r_Z(k) >= ( ZMax - DepthV ) ) then 
          xyr_Gamma(:,:,k) = xyr_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * (r_Z(k) - ZMax + DepthV)/ DepthV)) &
            &    / EFTime 
        end if
      end do
      
      ! z 方向には同じ
      !
      pyz_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if

    !鉛直方向の下部境界    
    !
    if ( DepthVb > 0.0d0 ) then 
      do k = kmin, kmax
        !スカラー格子点
        if ( z_Z(k) <= (ZMin + DepthVb) ) then 
          xyz_Gamma(:,:,k) = xyz_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * ( z_Z(k) - ( ZMin - DepthVb ) ) / DepthVb ) ) &
            &    / EFTime 
        end if

          ! 最大   z = ZMin + DepthVb
          !    gamma = 1.0 - cos( 0 ) = 0.0     ZMin
          ! 最小   z = ZMin
          !    gamma = 1.0 - cos( \pi ) = 1.0

        !フラックス格子点
        if ( r_Z(k) <= (ZMin + DepthVb) ) then 
          xyr_Gamma(:,:,k) = xyr_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * ( r_Z(k) - ( ZMin - DepthVb ) ) / DepthVb ) ) &
            &    / EFTime 
        end if
      end do
      
      ! z 方向には同じ
      !
      pyz_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if
      
    !-----------------------------------------------------------------    
    ! 値の確認
    !
    call MessageNotify( "M", "Damping_init", "EFTime = %f", d=(/EFTime/) )
    call MessageNotify( "M", "Damping_init", "DepthH = %f", d=(/DepthH/) )
    call MessageNotify( "M", "Damping_init", "DepthV = %f", d=(/DepthV/) )  
    call MessageNotify( "M", "Damping_init", "DepthVb= %f", d=(/DepthVb/) )  

!    write(*,*) minval( xyz_Gamma ), maxval( xyz_Gamma )
    
    !-----------------------------------------------------------------    
    ! 平均風のダンピング
    !
    allocate( &
      & pyz_SpngMF(imin:imax,jmin:jmax,kmin:kmax), &
      & xqz_SpngMF(imin:imax,jmin:jmax,kmin:kmax), &
      & z_VelX0(kmin:kmax),   &
      & z_VelY0(kmin:kmax)    )
    
    pyz_SpngMF = 0.0d0
    xqz_SpngMF = 0.0d0

    ! ダンピングさせる先の風速を得る
    !
    if ( InputFileMF == "" ) then 
      z_VelX0 = 0.0d0
      z_VelY0 = 0.0d0
    else
      call SpongeLayer_MeanFlow_Get( InputFileMF, z_VelX0, z_VelY0 )
    end if

    !-----------------------------------------------------------------    
    ! 値の確認
    !
    call MessageNotify( "M", "Damping_init", "DelTimeMF  = %f", d=(/DelTimeMF/) )
    call MessageNotify( "M", "Damping_init", "EFTimeMF   = %f", d=(/EFTimeMF/) )
    call MessageNotify( "M", "Damping_init", "ZMinMF     = %f", d=(/ZMinMF/) )
    call MessageNotify( "M", "Damping_init", "ZMaxMF     = %f", d=(/ZMaxMF/) )  
    call MessageNotify( "M", "Damping_init", "InputFileMF= %c", c1=trim(InputFileMF) )

    !-----------------------------------------------------------------    
    ! 出力
    !
    call HistoryAutoAddVariable(                          &
      & varname='DPTempDtSpng',                           &
      & dims=(/'x','y','z','t'/),                         &
      & longname='Damping term of potential temperature', &
      & units='K.s-1',                                    &
      & xtype='double')

    call HistoryAutoAddVariable(                    &
      & varname='DExnerDtSpng',                     &
      & dims=(/'x','y','z','t'/),                   &
      & longname='Damping term of exner function',  &
      & units='s-1',                                &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelXDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelX',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelYDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelY',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelZDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelZ',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(                      &
      & varname='DVelXDtSpngMF',                      &
      & dims=(/'x','y','z','t'/),                     &
      & longname='Damping term of VelX (Mean Flow)',  &
      & units='m.s-1',                                &
      & xtype='double')

    call HistoryAutoAddVariable(                      &
      & varname='DVelYDtSpngMF',                      &
      & dims=(/'x','y','z','t'/),                     &
      & longname='Damping term of VelY (Mean Flow)',  &
      & units='m.s-1',                                &
      & xtype='double')

  end subroutine Damping_Init


  subroutine SpongeLayer_forcing(                                         &
    & pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,  xyz_PTempBl,  xyz_ExnerBl,   & !(in) 
    & pyz_DVelXDt, xqz_DVelYDt, xyr_DVelZDt, xyz_DPTempDt, xyz_DExnerDt )   !(inout)

    use dc_types, only : DP
    use gtool_historyauto, only: HistoryAutoPut
    use timeset, only: TimeN        ! 現在の時刻
    use gridset, only: imin,       &! x 方向の配列の下限
      &                imax,       &! x 方向の配列の上限
      &                jmin,       &! y 方向の配列の下限
      &                jmax,       &! y 方向の配列の上限
      &                kmin,       &! z 方向の配列の下限
      &                kmax,       &! z 方向の配列の上限
      &                nx,         &! x 方向の物理領域の上限
      &                ny,         &! y 方向の物理領域の上限
      &                nz           ! z 方向の物理領域の上限

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax, jmin:jmax, kmin:kmax)   
    real(DP), intent(inout) :: pyz_DVelXDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyr_DVelZDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: pyz_SpngVelX(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xqz_SpngVelY(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyr_SpngVelZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_SpngPTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_SpngExner(imin:imax, jmin:jmax, kmin:kmax)

    !--------------------------------------------------------
    ! 擾乱場の値をゼロに戻すような強制を与える 
    !
    pyz_SpngVelX  =  - pyz_Gamma * pyz_VelXBl  * FactorSpngVelX
    xqz_SpngVelY  =  - xqz_Gamma * xqz_VelYBl  * FactorSpngVelY
    xyr_SpngVelZ  =  - xyr_Gamma * xyr_VelZBl  * FactorSpngVelZ
    xyz_SpngPTemp =  - xyz_Gamma * xyz_PTempBl * FactorSpngPTemp
    xyz_SpngExner =  - xyz_Gamma * xyz_ExnerBl * FactorSpngExner
    
    pyz_DVelXDt  = pyz_DVelXDt  + pyz_SpngVelX
    xqz_DVelYDt  = xqz_DVelYDt  + xqz_SpngVelY
    xyr_DVelZDt  = xyr_DVelZDt  + xyr_SpngVelZ
    xyz_DPTempDt = xyz_DPTempDt + xyz_SpngPTemp
    xyz_DExnerDt = xyz_DExnerDt + xyz_SpngExner
    
    call HistoryAutoPut(TimeN, 'DVelXDtSpng',  pyz_SpngVelX(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSpng',  xqz_SpngVelY(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSpng',  xyr_SpngVelZ(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtSpng', xyz_SpngPTemp(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtSpng', xyz_SpngExner(1:nx,1:ny,1:nz))

  end subroutine SpongeLayer_forcing

!!!
!!!================================================================================
!!!

  subroutine SpongeLayer_MeanFlow( &
    & pyz_VelX,    xqz_VelY,       & !(in) 
    & pyz_DVelXDt, xqz_DVelYDt     ) !(inout)

    use mpi_wrapper,       only: MPIWrapperAllreduce, nprocs
    use gtool_historyauto, only: HistoryAutoPut
    use dc_types, only: DP
    use timeset, only: TimeN, TimeA ! 現在の時刻
    use axesset, only: z_Z          !Z 座標軸(スカラー格子点)
    use gridset, only: imin,       &! x 方向の配列の下限
      &                imax,       &! x 方向の配列の上限
      &                jmin,       &! y 方向の配列の下限
      &                jmax,       &! y 方向の配列の上限
      &                kmin,       &! z 方向の配列の下限
      &                kmax,       &! z 方向の配列の上限
      &                nx,         &! x 方向の物理領域の上限
      &                ny,         &! y 方向の物理領域の上限
      &                nz           ! z 方向の物理領域の上限

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)    :: pyz_VelX(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xqz_VelY(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: z_VelX(nz)
    real(DP)                :: z_VelXSum(nz)
    real(DP)                :: z_VelXMean(nz)
    real(DP)                :: z_VelY(nz)
    real(DP)                :: z_VelYSum(nz)
    real(DP)                :: z_VelYMean(nz)
    integer                 :: i, j, k


    !------------------------------------------
    ! 平均流のダンプ係数を更新
    !
!!    write(*,*) TimeN, DelTimeMF, mod(TimeN, DelTimeMF)

    if ( mod(TimeA, DelTimeMF) == 0.0d0 ) then 

      !------------------------------------------
      ! ノード内で平均値の作成
      !
      z_VelX = 0.0d0
      z_VelY = 0.0d0
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            z_VelX(k) = z_VelX(k) + pyz_VelX(i,j,k)
            z_VelY(k) = z_VelY(k) + xqz_VelY(i,j,k)
          end do
        end do
      end do
      
      z_VelX = z_VelX / real( nx * ny, kind=DP)
      z_VelY = z_VelY / real( nx * ny, kind=DP)
      
      !--------------------------------------------------------
      ! MPI : 全ノードでの合計を得る
      !
      ! MPI を使わない場合は何も行わない 
      ! (z_VelXSum = z_VelX, z_VelYSum = z_VelY).
      !
      call MPIWrapperAllreduce(   &
        &  nz, z_VelX,        & !(in)
        &  z_VelXSum          & !(out)
        & )
      
      call MPIWrapperAllreduce(   &
        &  nz, z_VelY,        & !(in)
        &  z_VelYSum          & !(out)
        & )
      
      !--------------------------------------------------------
      ! ダンピング係数   
      !   - 1/T * ( \bar{u} / \bar{u}_max ) * ( \bar{u} - u0 )
      !   \bar{u} を \bar{u} の最大値で規格化する.
      !
      z_VelXMean = z_VelXSum / real(nprocs, kind=DP) 
      z_VelYMean = z_VelYSum / real(nprocs, kind=DP) 

      do k = 1, nz
        if ( z_Z(k) >= ZMinMF .OR.  z_Z(k) <= ZMaxMF ) then 
          pyz_SpngMF(:,:,k) = - ( z_VelXMean(k) - z_VelX0(k) ) / EFTimeMF 
          xqz_SpngMF(:,:,k) = - ( z_VelYMean(k) - z_VelY0(k) ) / EFTimeMF
        else
          pyz_SpngMF(:,:,k) = 0.0d0
          xqz_SpngMF(:,:,k) = 0.0d0
        end if
      end do
      
    end if

    !--------------------------------------------------------
    ! 擾乱場の値をゼロに戻すような強制を与える 
    !
    pyz_DVelXDt  = pyz_DVelXDt  + pyz_SpngMF
    xqz_DVelYDt  = xqz_DVelYDt  + xqz_SpngMF
    
    call HistoryAutoPut(TimeN, 'DVelXDtSpngMF',  pyz_SpngMF(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSpngMF',  xqz_SpngMF(1:nx,1:ny,1:nz))
    
  end subroutine SpongeLayer_MeanFlow


  subroutine SpongeLayer_MeanFlow_Get(   &
    & InputFileMF,                       & !(in)
    & z_VelX0, z_VelY0 )                   !(out)
    !
    !ファイルから情報取得
    !

    !モジュール読み込み
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!配列サイズ (X 方向)
      &                       jmin, jmax,    &!配列サイズ (Y 方向)
      &                       kmin, kmax      !配列サイズ (Z 方向)
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    character(STRING), intent(in) :: InputFileMF
    real(DP), intent(out) :: z_VelX0(kmin:kmax)
    real(DP), intent(out) :: z_VelY0(kmin:kmax)
    real(DP)              :: aaa_tmp1(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: aaa_tmp2(imin:imax,jmin:jmax,kmin:kmax)
    character(STRING)     :: name               !変数名
    character(STRING)     :: time = "t=^1"
    integer               :: k
    

    ! ファイルオープン & 値の取り出し
    ! 指定したファイルの先頭時刻の値を取り出す.
    !
    name = "VelX"
    call HistoryGet( InputFileMF, name, aaa_tmp1, range=Time, flag_mpi_split = FLAG_LIB_MPI )
    
    name = "VelY"
    call HistoryGet( InputFileMF, name, aaa_tmp2, range=Time, flag_mpi_split = FLAG_LIB_MPI )

    do k = kmin, kmax
      z_VelX0(k) = aaa_tmp1(1,1,k)
      z_VelY0(k) = aaa_tmp2(1,1,k)
    end do

  end subroutine SpongeLayer_MeanFlow_Get
  

end module Damping
