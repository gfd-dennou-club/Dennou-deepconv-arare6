!= Module cloudphys_k1969
!
! Authors::   杉山耕一朗(SUGIYAMA Ko-ichiro), 小高正嗣 (ODAKA Masatsugu), 高橋芳幸 (YOSHIYUKI Takahashi)
! Version::   $Id: cloudphys_k1969.f90,v 1.28 2014/03/04 05:55:05 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module cloudphys_k1969
  !
  !暖かい雨のバルク法を用いた, 水蒸気と雨, 雲と雨の混合比の変換係数を求める.
  !   * 中島健介 (1994) で利用した定式をそのまま利用. 
  ! 
  
  !モジュール読み込み
  use dc_types,   only : DP
  
  !暗黙の型宣言禁止
  implicit none
  
  !属性の指定
  private

  !関数を public にする
  public Cloudphys_K1969_Init
  public Cloudphys_K1969_forcing

  real(DP), save :: FactorCloud2Rain = 1.0d0 !雲から雨への変換の有無 
                                             !雨へ変換させない場合は値をゼロにする. 
  real(DP), save :: FactorRain2Gas = 1.0d0   !雨から蒸気への変換の有無 
                                             !蒸気へ変換させない場合は値をゼロにする. 
  real(DP), save, public :: FactorCloud2Gas = 1.0d0 
                                             !雲から蒸気への変換の有無 
                                             !蒸気へ変換させない場合は値をゼロにする. 

!  real(DP), save :: FactorJ      = 1.0d0 !雲物理過程のパラメータ
!                                         !木星では 3.0d0
!                                         !地球では 1.0d0 とする
  real(DP), save :: AutoConvTime = 1.0d3 !併合成長の時定数 [sec]
  real(DP), save :: QMixCr       = 1.0d-3 
                                         !併合成長を生じる臨界混合比 [kg/kg]
  real(DP), save :: FactorDExnerDtCloud = 1.0d0

contains  

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_Init
    !
    ! 初期化ルーチン
    !

    !モジュール呼び出し
    use dc_types,      only : DP, STRING
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use namelist_util, only : namelist_filename

    !暗黙の型宣言禁止
    implicit none

    !内部変数
    integer  :: unit    !装置番号
    character(*), parameter:: module_name = 'Cloudphys_K1969_Init'

    !-----------------------------------------------------------
    ! NAMELIST から情報を取得
    !
    NAMELIST /cloudphys_k1969_nml/                &
      & AutoConvTime, QMixCr,                     &
      & FactorDExnerDtCloud,                      &
      & FactorCloud2Rain, FactorRain2Gas, FactorCloud2Gas

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=cloudphys_k1969_nml)
    close(unit)

    !-----------------------------------------------------------
    ! 出力
    !
    call MessageNotify( "M", &
      &  module_name, "AutoConvTime = %f",  d=(/AutoConvTime/) )
    call MessageNotify( "M", &
      &  module_name, "QMixCr = %f",  d=(/QMixCr/) )
    call MessageNotify( "M", &
      &  module_name, "FactorCloud2Rain = %f",  d=(/FactorCloud2Rain/) )
    call MessageNotify( "M", &
      &  module_name, "FactorRain2Gas = %f",  d=(/FactorRain2Gas/) )
    call MessageNotify( "M", &
      &  module_name, "FactorCloud2Gas = %f",  d=(/FactorCloud2Gas/) )
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtCloud= %f", d=(/ FactorDExnerDtCloud /))

     !
     ! HistoryAuto
     !
     call Cloudphys_K1969_HistoryAuto

    
   end subroutine Cloudphys_K1969_Init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_forcing(         &
    & xyz_ExnerNl,                            &!(in)
    & xyz_DExnerDt, xyz_PTempAl, xyzf_QMixAl  &!(inout)
    & )

    !モジュール呼び出し
    use dc_types,only : DP, STRING
    use gtool_historyauto,                 &
      &          only : HistoryAutoPut
    use timeset, only : DelTimeLong, TimeN
    use gridset, only : imin,              &!x 方向の配列の下限
      &                 imax,              &!x 方向の配列の上限
      &                 jmin,              &!y 方向の配列の上限
      &                 jmax,              &!y 方向の配列の上限
      &                 kmin,              &!z 方向の配列の下限
      &                 kmax,              &!z 方向の配列の上限
      &                 nx, ny, nz, ncmax   !物理領域の大きさ
    use constants,only: FactorJ,           &!
      &                 PressBasis,        &!温位の基準圧力 
      &                 CpDry,             &!乾燥成分の比熱
      &                 MolWtDry,          &!
      &                 GasRDry             !乾燥成分の気体定数 
    use basicset, only: xyz_DensBZ,        &!基本場の密度
      &                 xyz_PTempBZ,       &!基本場の温位
      &                 xyz_ExnerBZ,       &!基本場の無次元圧力
      &                 xyzf_QMixBZ         !基本場の混合比
    use composition,                       &
      &           only:  MolWtWet,         &!
      &                 SpcWetID,          &!
      &                 SpcWetSymbol,      &!
      &                 CondNum,           &!凝結過程の数
      &                 IdxCG,             &!凝結過程(蒸気)の配列添え字
      &                 IdxCC,             &!凝結過程(雲)の配列添え字
      &                 IdxCR,             &!凝結過程(雲)の配列添え字
      &                 GasNum,            &!蒸気の数
      &                 CloudNum,          &!雲の数
      &                 RainNum,           &!雨の数
      &                 IdxG,              &!蒸気の配列添え字
      &                 IdxC,              &!雲の配列添え字
      &                 IdxR,              &!雨の配列添え字
      &                 IdxNH3,            &!NH3(蒸気)の配列添え字
      &                 IdxH2S,            &!H2S(蒸気)の配列添え字
      &                 IdxNH4SHr           !NH4SH(雨)の配列添え字
    use average, only : xyr_xyz
    use differentiate_center2,             &
      &          only : xyz_dz_xyr
    use ChemCalc,only : xyz_SvapPress, xyz_LatentHeat, ReactHeatNH4SH, xyz_DelQMixNH4SH    
    use MoistAdjust,                       &
      &          only : MoistAdjustSvapPress, MoistAdjustNH4SH
    use DExnerDt,only : xyz_DExnerDt_xyzf, xyz_DExnerDt_xyz
    use SetMargin,only: SetMargin_xyz, SetMargin_xyzf

    !暗黙の型宣言禁止
    implicit none

    real(DP), intent(in)    :: xyz_ExnerNl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_PTempAl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyzf_QMixAl(imin:imax, jmin:jmax, kmin:kmax, ncmax)

    real(DP) :: xyz_PTempOrig(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_PTempWork(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_DelPTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyz_PTempCond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP) :: xyzf_QMixOrig(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_QMixWork(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_DelQMix(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: xyzf_QMixCond(imin:imax, jmin:jmax, kmin:kmax, ncmax)
    real(DP) :: DelTime
    integer  :: s
    integer  :: iG, iC, iR

    real(DP) :: xyzf_Cloud2Rain(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                          !雲から雨への変換量
    real(DP) :: xyz_AutoConv(imin:imax,jmin:jmax,kmin:kmax)
                                          !飽和混合比
    real(DP) :: xyz_Collect(imin:imax,jmin:jmax,kmin:kmax)
                                          !規格化された潜熱
    real(DP) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                          !混合比の擾乱成分 + 平均成分
    real(DP) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                          !温度の擾乱成分 + 平均成分
    real(DP) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
                                          !全圧
    real(DP) :: xyz_ExnerAll(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_NonSaturate(imin:imax,jmin:jmax,kmin:kmax)
                                          !未飽和度(飽和混合比と蒸気の混合比の差)
    real(DP) :: xyzf_Rain2Gas(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyzf_Rain2GasNH4SH(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyzf_DelPTemp(imin:imax,jmin:jmax,kmin:kmax, ncmax)
    real(DP) :: xyz_DelPTempNH4SH(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtCondTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_DExnerDtCondQMix(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_QMixSat(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP) :: xyz_QMixHum(imin:imax,jmin:jmax,kmin:kmax) 


    !-------------------------------------------------------------
    ! 初期値を保管 Store Initial Value
    !
    xyz_PTempOrig = xyz_PTempAl
    xyzf_QMixOrig = xyzf_QMixAl

    !-------------------------------------------------------------
    ! 時間刻み幅. Leap-frog なので, 2 \del t
    !
    DelTime = 2.0d0 * DelTimeLong

    !------------------------------------------
    ! 全エクスナー関数・全圧を計算. サブルーチン内では変化しない.
    !
    xyz_ExnerAll = xyz_ExnerNl + xyz_ExnerBZ
    xyz_PressAll = PressBasis * (xyz_ExnerAll ** (CpDry / GasRDry))

    !------------------------------------------    
    ! 暖かい雨のパラメタリゼーション.
    ! * 雲<-->雨 の変換を行う.
    !
    ! Warm rain parameterization.
    ! * Conversion from cloud to rain.
    
    !これまでの値を作業配列に保管
    ! Previous values are stored to work area.
    !
    xyzf_QMixWork = xyzf_QMixAl
    
    !雨への変化量を計算
    ! Conversion values are calculated.
    !    
    xyzf_QMixAll = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyzf_Cloud2Rain = 0.0d0

    do s = 1, CloudNum

      ! 値を保管
      !
      iC = IdxC(s)
      iR = IdxR(s)

      !併合成長
      !
      xyz_AutoConv =                                           &
        & DelTime / AutoConvTime                               &
        & * max( 0.0d0, ( xyzf_QMixAll(:,:,:,iC) - QMixCr) )

      !衝突合体成長
      !
      xyz_Collect =                                            &
        &  DelTime                                             &
        &  * 2.2d0 * FactorJ * xyzf_QMixAll(:,:,:,iC)          &
        &  * (xyzf_QMixAll(:,:,:,iR) * xyz_DensBZ) ** 0.875d0  

      !雲の変換量: 併合成長と合体衝突の和
      !  元々の変化量を上限値として設定する. 負の値となる.
      !
      xyzf_Cloud2Rain(:,:,:,iC) =                                        &
        & - min( xyzf_QMixAll(:,:,:,iC), ( xyz_AutoConv + xyz_Collect ) )
      
      !雨の変換量. 符号は雲の変換量とは反対. 
      xyzf_Cloud2Rain(:,:,:,iR) = - xyzf_Cloud2Rain(:,:,:,iC) 

      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iC))//'DtC2R',   &
        & xyzf_Cloud2Rain(1:nx,1:ny,1:nz,iC) / DelTime)

      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtC2R',   &
        & xyzf_Cloud2Rain(1:nx,1:ny,1:nz,iR) / DelTime)

    end do

    ! 変化量を足し込む
    ! 雲から雨へ変換させない場合は FactorCloud2Rain = 0.0 とする. 
    !
    xyzf_QMixAl = xyzf_QMixWork + xyzf_Cloud2Rain * FactorCloud2Rain


    !-------------------------------------------    
    ! 暖かい雨のパラメタリゼーション.
    ! * 蒸気<-->雨 の変換を行う
    !
    ! Warm rain parameterization.
    ! * Conversion from rain to vapor.
    
    !これまでの値を作業配列に保管
    ! Previous values are stored to work area.
    !
    xyz_PTempWork = xyz_PTempAl
    xyzf_QMixWork = xyzf_QMixAl
    
    ! 雨から蒸気への混合比変化を求める
    ! * 温位の計算において, 混合比変化が必要となるため, 
    !   混合比変化を 1 つの配列として用意する.
    !
    ! Conversion values are calculated.
    !

    !温度, 圧力, 混合比の全量を求める
    !擾乱成分と平均成分の足し算
    !
    xyz_TempAll   = ( xyz_PTempAl + xyz_PTempBZ ) * ( xyz_ExnerNl + xyz_ExnerBZ )
    xyzf_QMixAll  = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry  = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )

    xyzf_Rain2Gas = 0.0d0
    xyzf_DelPTemp = 0.0d0
    xyzf_Rain2GasNH4SH = 0.0d0
    xyz_DelPTempNH4SH  = 0.0d0

    do s = 1, CondNum

       ! 値を保管
       !
       iG = IdxCG(s)
       iC = IdxCC(s)
       iR = IdxCR(s)

      !飽和蒸気圧と混合比の差(飽和度)を計算. 
      !  雨から蒸気への変換量は飽和度に比例する.
      !
      xyz_NonSaturate =                                   &
        & max(                                            &
        &   0.0d0,                                        &
        &   xyz_SvapPress(SpcWetID(iC), xyz_TempAll)      &
        &     * MolWtWet(iG) / ( MolWtDry * xyz_PressDry) &
        &     - xyzf_QMixAll(:,:,:,iG)                    &
        &    )

      !雨の変換量
      !  元々の雨粒の混合比以上に蒸発が生じないように上限値を設定
      !
      xyzf_Rain2Gas(:,:,:,iR) =                                    &
        & - min(                                                   &
        &    DelTime * 4.85d-2 * FactorJ * xyz_NonSaturate         &
        &     * ( xyzf_QMixAll(:,:,:,iR) * xyz_DensBZ )** 0.65d0,  &
        &    xyzf_QMixAll(:,:,:,iR)                                &
        &   ) 

      !蒸気の変換量
      !  雨粒の変換量とは符号が逆となる
      !
      xyzf_Rain2Gas(:,:,:,iG) = - xyzf_Rain2Gas(:,:,:,iR) 
    
      ! xyzf_DelQMix を元に潜熱を計算
      !
      xyzf_DelPTemp(:,:,:,s) =                          &
        & xyz_LatentHeat( SpcWetID(iR), xyz_TempAll )   &
        &  * xyzf_Rain2Gas(:,:,:,iR)                    &
        &  / (xyz_ExnerAll * CpDry) 

    end do

    !飽和蒸気圧と混合比の差(飽和度)を計算. 
    !  雨から蒸気への変換量は飽和度に比例する.
    !  未飽和度を求めたいので, マイナスをかけ算している
    !  (DelQMixNH4SH は, NH4SH が増加する方向, すなわち飽和度を正としている)
    !
    if (IdxNH4SHr /= 0) then 
      xyz_NonSaturate =                                                 &
        & max(                                                          &
        &  0.0d0,                                                       &
        &   - xyz_DelQMixNH4SH(                                         &  
        &       xyz_TempAll, xyz_PressAll, xyz_PressDry,                & 
!        &       xyz_TempAll, xyz_PressAll,                              & 
        &       xyzf_QMixAll(:,:,:,IdxNH3), xyzf_QMixAll(:,:,:,IdxH2S), &
        &       MolWtWet(IdxNH3), MolWtWet(IdxH2S)                      &
        &     )                                                         &
        &  )

      !雨の変換量
      !  元々の雨粒の混合比以上に蒸発が生じないように上限値を設定
      !
      xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) =                              &
        & - min(                                                         &
        &     DelTime * 4.85d-2 * FactorJ * xyz_NonSaturate              &
        &      * (xyzf_QMixAll(:,:,:,IdxNH4SHr) * xyz_DensBZ) ** 0.65d0, &
        &     xyzf_QMixAll(:,:,:,IdxNH4SHr)                              &
        &    ) 
     
      !蒸気の変換量
      !  雨粒の変換量とは符号が逆となる
      !
      xyzf_Rain2GasNH4SH(:,:,:,IdxNH3) =                           &
        & - xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) * MolWtWet(IdxNH3) &
        &   / MolWtWet(IdxNH4SHr)
      xyzf_Rain2GasNH4SH(:,:,:,IdxH2S) =                           &
        & - xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) * MolWtWet(IdxH2S) &
        &   / MolWtWet(IdxNH4SHr)

      xyz_DelPTempNH4SH                                          &
        & = ReactHeatNH4SH * xyzf_Rain2GasNH4SH(:,:,:,IdxNH4SHr) &
        &    / (xyz_ExnerAll * CpDry)

    end if

    !変化量を足し算
    !雨から蒸気への変換を切る場合は FactorRain2Gas = 0.0 とする. 
    !
    xyzf_DelQMix = xyzf_Rain2Gas + xyzf_Rain2GasNH4SH 
    xyz_DelPTemp = sum(xyzf_DelPTemp, 4) + xyz_DelPTempNH4SH 

    ! 温位と混合比の計算. 雨から蒸気への変換分を追加
    !
    xyz_PTempAl = xyz_PTempWork + xyz_DelPTemp * FactorRain2Gas
    xyzf_QMixAl = xyzf_QMixWork + xyzf_DelQMix * FactorRain2Gas

    call HistoryAutoPut(TimeN, 'PTempR2G',  &
      &  xyz_DelPTemp(1:nx,1:ny,1:nz) / DelTime)
    do s = 1, GasNum
      iG = IdxG(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtR2G',  &
        & xyzf_DelQMix(1:nx,1:ny,1:nz, iG) / DelTime)
    end do
    do s = 1, RainNum
      iR = IdxR(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iR))//'DtR2G',  &
        & xyzf_DelQMix(1:nx,1:ny,1:nz, iR) / DelTime)
    end do

    !-------------------------------------------
    ! 湿潤飽和調節
    ! * 蒸気<-->雲の変換を行う.
    !
    ! Moist adjustment.
    ! * Conversion from vapor to cloud.
    !

    !これまでの値を作業配列に保管
    ! Previous values are stored to work area.
    !
    xyz_PTempWork = xyz_PTempAl
    xyzf_QMixWork = xyzf_QMixAl
    
    ! 乾燥成分の圧力と混合比の全量を求める
    !
    xyzf_QMixAll = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )
    
    ! 湿潤飽和調節
    !
    call MoistAdjustSvapPress(   &
      & xyz_PressDry,            & ! (in) 
      & xyz_ExnerNl,             & ! (in)
      & xyz_PTempAl,             & ! (inout)
      & xyzf_QMixAl,             & ! (inout)
      & FactorCloud2Gas          & ! (in)
      & )
    if (IdxNH4SHr /= 0) then 
      call MoistAdjustNH4SH(     &
        & xyz_PressDry,          & !(in)
        & xyz_ExnerNl,           & !(in)
        & xyz_PTempAl,           & !(inout)
        & xyzf_QMixAl,           & !(inout)
        & FactorCloud2Gas        & !(in)
        & )
    end if

    ! Output
    !
    xyz_PTempCond = (xyz_PTempAl - xyz_PTempWork) / DelTime
    xyzf_QMixCond = (xyzf_QMixAl - xyzf_QMixWork) / DelTime

    call HistoryAutoPut(TimeN, 'PTempG2C', xyz_PTempCond(1:nx,1:ny,1:nz))
    do s = 1, GasNum
      iG = idxG(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtG2C',  &
        & xyzf_QMixCond(1:nx,1:ny,1:nz,iG) )
    end do
    do s = 1, CloudNum
      iC = idxC(s)
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iC))//'DtG2C',  &
        & xyzf_QMixCond(1:nx,1:ny,1:nz,iC) )
    end do


    !-----------------------------------
    ! エクスナー関数の tendency
    !

    ! 初期値からの差を取る
    !
    xyz_PTempCond = (xyz_PTempAl - xyz_PTempOrig) / DelTime
    xyzf_QMixCond = (xyzf_QMixAl - xyzf_QMixOrig) / DelTime

    xyz_DExnerDtCondTemp = xyz_DExnerDt_xyz( xyz_PTempCond )  * FactorDExnerDtCloud
    xyz_DExnerDtCondQMix = xyz_DExnerDt_xyzf( xyzf_QMixCond ) * FactorDExnerDtCloud

    xyz_DExnerDt  = xyz_DExnerDt + xyz_DExnerDtCondTemp + xyz_DExnerDtCondQMix

    call HistoryAutoPut(TimeN, 'DExnerDtCondTemp', xyz_DExnerDtCondTemp(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtCondQMix', xyz_DExnerDtCondQMix(1:nx,1:ny,1:nz))


    !----------------------------------------------------------------
    ! 飽和蒸気圧と平衡定数
    !
    xyz_TempAll   = ( xyz_PTempAl + xyz_PTempBZ ) * ( xyz_ExnerNl + xyz_ExnerBZ )
    xyzf_QMixAll  = max( 0.0d0, xyzf_QMixAl + xyzf_QMixBZ )
    xyz_PressDry  = xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )

    do s = 1, CondNum
      iG = IdxCG(s)
      iC = IdxCC(s)

      !飽和蒸気圧
      xyz_QMixSat =                                  &
        & xyz_SvapPress(SpcWetID(iC), xyz_TempAll)    &
        &  * MolWtWet(iC) / MolWtDry / xyz_PressDry

      !相対湿度
      xyz_QMixHum = &
        & xyzf_QMixAll(:,:,:,iG) / xyz_QMixSat * 100.0d0

      !出力
      call HistoryAutoPut(TimeN, 'D'//trim(SpcWetSymbol(iG))//'DtHum', xyz_QMixHum(1:nx, 1:ny, 1:nz))
    end do

    ! Set Margin
    !
    call SetMargin_xyz(xyz_PTempAl)
    call SetMargin_xyzf(xyzf_QMixAl)

  end subroutine Cloudphys_K1969_forcing

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xyz_PressDry_xyzf_xyz( xyzf_QMixAll, xyz_PressAll )
    !
    ! 乾燥成分の分圧の計算
    ! 

    !モジュール呼び出し
    use dc_types,    only : DP
    use gridset,     only : imin, imax, &
      &                     jmin, jmax, &
      &                     kmin, kmax, &
      &                     ncmax
    use composition, only : MolWtWet,   &!
      &                     IdxG,       &!蒸気の配列添え字
      &                     GasNum
    use constants,   only : MolWtDry

    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP), intent(in) :: xyzf_QMixAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(in) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_PressDry_xyzf_xyz(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_QMixAllPerMolWt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    integer              :: f, iG

    ! 混合比/分子量を計算
    ! 
    do f = 1, GasNum
      iG = IdxG(f)
      xyzf_QMixAllPerMolWt(:,:,:,f) = xyzf_QMixAll(:,:,:,iG) / MolWtWet(iG)
    end do

    ! 乾燥成分の分圧. 
    ! 
    xyz_PressDry_xyzf_xyz = xyz_PressAll / (1.0d0 + MolWtDry * sum( xyzf_QMixAllPerMolWt(:,:,:,1:GasNum), 4))

!!DBG
!!    xyz_PressDry_xyzf_xyz = xyz_PressAll

  end function xyz_PressDry_xyzf_xyz


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine Cloudphys_K1969_HistoryAuto
    !
    ! tendency 出力のための設定
    !

    !モジュール呼び出し
    use gtool_historyauto, only : HistoryAutoAddVariable
    use composition,       only : SpcWetSymbol,  & 
      &                           GasNum,        &!蒸気の数
      &                           CloudNum,      &!雲の数
      &                           RainNum,       &!雨の数  
      &                           IdxG,          &!蒸気の配列添え字
      &                           IdxC,          &!雲の配列添え字
      &                           IdxR            !雨の配列添え字

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer  :: l, iG, iC, iR

    call HistoryAutoAddVariable(                                     &
      & varname='PTempR2G',                                          &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Condensation term of potential temperature (R2G)', &
      & units='K.s-1',                                               &
      & xtype='double' )

    call HistoryAutoAddVariable(                                     &
      & varname='PTempG2C',                                          &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Condensation term of potential temperature (G2C)', &
      & units='K.s-1',                                               &
      & xtype='double' )

    call HistoryAutoAddVariable(                                     &
      & varname='DExnerDtCondTemp',                                     &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Latent heat term of exner function (Temp)',        &
      & units='K.s-1',                                               &
      & xtype='double')

    call HistoryAutoAddVariable(                                     &
      & varname='DExnerDtCondQMix',                                     &
      & dims=(/'x','y','z','t'/),                                    &
      & longname='Latent heat term of exner function (QMix)',        &
      & units='K.s-1',                                               &
      & xtype='double')

    do l = 1, GasNum
      iG = IdxG(l)

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtR2G',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iG))//' mixing ratio (R2G)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtG2C',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iG))//' mixing ratio (G2C)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iG))//'DtHum',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Humidity of ' //trim(SpcWetSymbol(iG)),          &
        & units='1',                                                 &
        & xtype='double')
      
    end do
    
    do l = 1, CloudNum
      iC = IdxC(l)
      
      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iC))//'DtC2R',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iC))//' mixing ratio (C2R)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )
      
      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iC))//'DtG2C',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iC))//' mixing ratio (G2C)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

    end do


    do l = 1, RainNum
      iR = IdxR(l)

      call HistoryAutoAddVariable(                                   &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtR2G',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio (R2G)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

     call HistoryAutoAddVariable(                                    &
        & varname='D'//trim(SpcWetSymbol(iR))//'DtC2R',                    &
        & dims=(/'x','y','z','t'/),                                  &
        & longname='Condensation term of '                           &
        &           //trim(SpcWetSymbol(iR))//' mixing ratio (C2R)', &
        & units='kg.kg-1.s-1',                                       &
        & xtype='double' )

    end do

  end subroutine Cloudphys_K1969_HistoryAuto
  
  
end module Cloudphys_k1969
