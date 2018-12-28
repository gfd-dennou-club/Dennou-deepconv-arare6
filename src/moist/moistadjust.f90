!= Module MoistAdjust
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: moistadjust.f90,v 1.7 2013/01/30 04:25:50 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

Module MoistAdjust
  !
  !湿潤飽和調節法を行うためのパッケージ型モジュール
  !  * 飽和蒸気圧を使う化学種は, 一元的に実行する
  !  * 化学反応については, それぞれの化学反応毎に実行する. 

  !暗黙の型宣言禁止
  implicit none
  
  !属性の指定
  private

  !関数の属性を public に変更
  public MoistAdjustSvapPress    !飽和蒸気圧を用いた飽和湿潤調節(簡易版)
  public MoistAdjustNH4SH        !化学反応の圧平衡定数を用いた飽和湿潤調節

contains

!!!------------------------------------------------------------------!!!
  subroutine MoistAdjustSvapPress(xyz_PressDry, xyz_Exner, xyz_PTemp, xyzf_QMix, FactorCloud2Gas)
    !
    ! 飽和蒸気圧を用いた湿潤飽和調整法の実行
    ! この副プログラムでは, 予め決めた回数だけ反復改良を行う. 
    !

    !モジュール読み込み
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : ncmax,             &!化学種の数
      &                    imin, imax,        &!x 方向の配列の上限と下限
      &                    jmin, jmax,        &!y 方向の配列の上限
      &                    kmin, kmax          !z 方向の配列の下限
    use basicset,   only : xyz_ExnerBZ,       &!無次元圧力(基本場)
      &                    xyz_PTempBZ,       &!温位(基本場)
      &                    xyzf_QMixBZ         !凝縮成分混合比(基本場)
    use constants,  only : CpDry,             &!乾燥成分の平均定圧比熱 [J/K kg]
      &                    MolWtDry            !乾燥成分の平均分子量   [kg/mol]
    use composition,only : MolWtWet,          &!湿潤成分の分子量  
      &                    SpcWetID,          &!湿潤成分の化学種のID 
      &                    CondNum,           &!凝結過程の数
      &                    IdxCG,             &!凝結過程(蒸気)の配列添え字
      &                    IdxCC               !凝結過程(雲)の配列添え字
      use ChemCalc,   only : xyz_SvapPress, xyz_LatentHeat, xyz_DQMixSatDPTemp
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in)   :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !乾燥成分の圧力
    real(DP),intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !エクスナー関数
    real(DP),intent(inout):: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !温位
    real(DP),intent(inout):: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                          !混合比
    real(DP),intent(in)   :: FactorCloud2Gas
    integer, parameter    :: ItrNum = 4                   !反復改良の回数

    real(DP):: xyz_QMixV_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixV_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixC_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixC_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_PTemp_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_PTemp_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_ExnerAll(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_DelQMix(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixSat(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_Cond(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP):: xyz_Evap(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP):: xyz_Gamma(imin:imax,jmin:jmax,kmin:kmax) 
    integer :: i, s                                       ! 添え字  
    integer :: IG, iC

    !---------------------------------------------------------------------
    ! 湿潤飽和調節法の実行
    !   ループを回すのは, 雲についてだけ.  
    !---------------------------------------------------------------------
    LoopSvapPress: do s = 1, CondNum

       ! 添字を保管
       iG = IdxCG(s)
       iC = IdxCC(s)
      
      !湿潤飽和法では圧力は変化させない. 
      xyz_ExnerAll = xyz_Exner + xyz_ExnerBZ
          
      !今までに得られた擾乱成分の値を暫定量とみなす. 添え字を追加
      xyz_QMixV_pre  = xyzf_QMix(:,:,:,iG) + xyzf_QMixBZ(:,:,:,iG)
      xyz_QMixC_pre  = xyzf_QMix(:,:,:,iC) + xyzf_QMixBZ(:,:,:,iC)
      xyz_PTemp_pre = xyz_PTemp

      Adjusting: do i = 1, ItrNum
        !---------------------------------------------------------------
        ! 飽和蒸気圧から飽和混合比を求める
        !---------------------------------------------------------------
        !温度
        xyz_TempAll = ( xyz_PTemp_pre + xyz_PTempBZ ) * xyz_ExnerAll

        !飽和蒸気圧から飽和混合比を計算(基本場からの差). 
        xyz_QMixSat =                                    &
          & xyz_SvapPress(SpcWetID(iC), xyz_TempAll)     &
          &  * MolWtWet(iC) / MolWtDry / xyz_PressDry

        !規格化された潜熱
        xyz_Gamma = xyz_LatentHeat(SpcWetID(iC), xyz_TempAll) &
          &        / (xyz_ExnerAll * CpDry)

        !凝結量を求める. 
        !  凝結が生じる場合には, xz_QMixV_pre - xz_QMixSat は必ず正の値となる.
        !  蒸発が生じる場合には, 蒸発量は - QMixC を超えることはない. 
        xyz_DelQMix =                                                    &
          & ( xyz_QMixV_pre - xyz_QMixSat )                              &
          &   / (1.0d0 + xyz_Gamma * xyz_DQMixSatDPTemp(                 &
          &        SpcWetID(iC), MolWtWet(iC), xyz_TempAll, xyz_ExnerAll & 
          &        ) )

        xyz_Cond = max( 0.0d0, min( xyz_QMixV_pre,   xyz_DelQMix ) )
        xyz_Evap = max( 0.0d0, min( xyz_QMixC_pre, - xyz_DelQMix ) ) * FactorCloud2Gas
        
        !より真に近い値を計算
        xyz_PTemp_nxt  = xyz_PTemp_pre + xyz_Gamma * ( xyz_Cond - xyz_Evap )
        xyz_QMixV_nxt  = xyz_QMixV_pre  - xyz_Cond + xyz_Evap
        xyz_QMixC_nxt  = xyz_QMixC_pre  + xyz_Cond - xyz_Evap

        !繰り返しのための変数定義
        xyz_PTemp_pre  = xyz_PTemp_nxt
        xyz_QMixV_pre  = xyz_QMixV_nxt
        xyz_QMixC_pre  = xyz_QMixC_nxt

      end do Adjusting
      
      xyz_PTemp           = xyz_PTemp_nxt                 
      xyzf_QMix(:,:,:,iG) = xyz_QMixV_nxt - xyzf_QMixBZ(:,:,:,iG)
      xyzf_QMix(:,:,:,iC) = xyz_QMixC_nxt - xyzf_QMixBZ(:,:,:,iC)
      
    end do LoopSvapPress
    
  end subroutine MoistAdjustSvapPress


!!!--------------------------------------------------------------------------!!!
  subroutine MoistAdjustNH4SH(xyz_PressDry, xyz_Exner, xyz_PTemp, xyzf_QMix, FactorCloud2Gas ) 
    !
    ! NH3 + H2S --> NH4SH の生成反応の圧平衡定数 Kp を用いた飽和湿潤調節法
    !

    !モジュール読み込み
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : ncmax,             &!化学種の数
      &                    imin, imax,        &!x 方向の配列の上限と下限
      &                    jmin, jmax,        &!y 方向の配列の上限
      &                    kmin, kmax          !z 方向の配列の下限
    use basicset,   only : xyz_ExnerBZ,       &!無次元圧力(基本場)
      &                    xyz_PTempBZ,       &!温位(基本場)
      &                    xyzf_QMixBZ         !凝縮成分混合比(基本場)
    use constants,  only : PressBasis,        &!温位の基準圧力         [Pa]
      &                    CpDry,             &!乾燥成分の平均定圧比熱 [J/K kg]
      &                    GasRDry             !乾燥成分の気体定数     [J/K kg]
    use composition,only : MolWtWet,          &!湿潤成分の分子量  
      &                    IdxNH3,            &!NH3(蒸気)の配列添え字
      &                    IdxH2S,            &!H2S(蒸気)の配列添え字
      &                    IdxNH4SHc           !NH4SH(雲)の配列添え字
    use ChemCalc,   only : ReactHeatNH4SH, xyz_DelQMixNH4SH
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数定義
    real(DP),intent(in)   :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !乾燥成分の圧力
    real(DP),intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !エクスナー関数
    real(DP),intent(inout):: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !温位
    real(DP),intent(inout):: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                          !凝縮成分の混合比
    real(DP),intent(in)   :: FactorCloud2Gas

    real(DP):: xyz_PTemp_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_PTemp_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixNH3_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixNH3_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixH2S_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixH2S_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixNH4SH_pre(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_QMixNH4SH_nxt(imin:imax,jmin:jmax,kmin:kmax) 
    
    real(DP):: xyz_ExnerAll(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_Gamma(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_DelQMix(imin:imax,jmin:jmax,kmin:kmax) 
    real(DP):: xyz_Cond(imin:imax,jmin:jmax,kmin:kmax)  
    real(DP)::xyz_Evap(imin:imax,jmin:jmax,kmin:kmax)  
    
    integer            :: i
    integer, parameter :: ItrNum = 2
            
    !---------------------------------------------------------------------
    ! 初期化
    !---------------------------------------------------------------------

    !湿潤飽和法では圧力は変化させない. 
    xyz_ExnerAll = xyz_Exner + xyz_ExnerBZ
    xyz_PressAll = PressBasis * (xyz_ExnerAll ** (CpDry / GasRDry))
        
    !今までに得られた擾乱成分の値を暫定量とみなす. 添え字を追加
    xyz_QMixNH3_pre(:,:,:)   = xyzf_QMix(:,:,:,IdxNH3) &
      &                           + xyzf_QMixBZ(:,:,:,IdxNH3)
    xyz_QMixH2S_pre(:,:,:)   = xyzf_QMix(:,:,:,IdxH2S) &
      &                           + xyzf_QMixBZ(:,:,:,IdxH2S)
    xyz_QMixNH4SH_pre(:,:,:) = xyzf_QMix(:,:,:,IdxNH4SHc) &
      &                           + xyzf_QMixBZ(:,:,:,IdxNH4SHc)
    xyz_PTemp_pre            = xyz_PTemp

    !規格化された反応熱 (NH4SH １kg に対する熱量)
    xyz_Gamma = ReactHeatNH4SH / ( xyz_ExnerAll * CpDry )

    AdjustNH4SH: do i = 1, ItrNum
      !---------------------------------------------------------------
      ! 変数の初期化
      !---------------------------------------------------------------
      !温度
      xyz_TempAll = ( xyz_PTemp_pre + xyz_PTempBZ ) * xyz_ExnerAll
      
      !NH4SH の生成量
      xyz_DelQMix =                                        &
        &    xyz_DelQMixNH4SH(                             &
        &        xyz_TempAll, xyz_PressAll, xyz_PressDry,  &
!        &        xyz_TempAll, xyz_PressAll,                &
        &        xyz_QMixNH3_pre, xyz_QMixH2S_pre,         &
        &        MolWtWet(IdxNH3), MolWtWet(IdxH2S)        &
        &      )

      xyz_Cond = max( 0.0d0, xyz_DelQMix )
      xyz_Evap = max( 0.0d0, min( - xyz_DelQMix, xyz_QMixNH4SH_pre ) ) * FactorCloud2Gas

      !---------------------------------------------------------------
      ! より真に近い値を求める飽和蒸気圧から飽和混合比を求める
      !---------------------------------------------------------------
      ! NH4SH の混合比を修正
      xyz_QMixNH4SH_nxt  = xyz_QMixNH4SH_pre + xyz_Cond - xyz_Evap
      
      ! DelPress を元に, NH3 と H2S の混合比を修正
      xyz_QMixNH3_nxt = xyz_QMixNH3_pre - ( xyz_Cond - xyz_Evap )      &
        &                 * MolWtWet(IdxNH3) / MolWtWet(IdxNH4SHc)
      xyz_QMixH2S_nxt = xyz_QMixH2S_pre - ( xyz_Cond - xyz_Evap )      &
        &                 * MolWtWet(IdxH2S) / MolWtWet(IdxNH4SHc)
          
      !温位を修正
      xyz_PTemp_nxt = xyz_PTemp_pre + xyz_Gamma * ( xyz_Cond - xyz_Evap )
      
      !ループを回すための変数変化
      xyz_PTemp_pre    = xyz_PTemp_nxt
      xyz_QMixNH3_pre   = xyz_QMixNH3_nxt 
      xyz_QMixH2S_pre   = xyz_QMixH2S_nxt 
      xyz_QMixNH4SH_pre = xyz_QMixNH4SH_nxt 

    end do AdjustNH4SH

    xyz_PTemp                  = xyz_PTemp_nxt                 
    xyzf_QMix(:,:,:,IdxNH3)    = xyz_QMixNH3_nxt   - xyzf_QMixBZ(:,:,:,IdxNH3)
    xyzf_QMix(:,:,:,IdxH2S)    = xyz_QMixH2S_nxt   - xyzf_QMixBZ(:,:,:,IdxH2S)
    xyzf_QMix(:,:,:,IdxNH4SHc) = xyz_QMixNH4SH_nxt - xyzf_QMixBZ(:,:,:,IdxNH4SHc)
    
  end subroutine MoistAdjustNH4SH

end Module MoistAdjust
