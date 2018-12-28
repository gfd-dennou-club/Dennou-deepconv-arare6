!= Module MoistAdjust
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: moistadjust.f90,v 1.7 2013/01/30 04:25:50 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

Module MoistAdjust
  !
  !����˰��Ĵ��ˡ��Ԥ�����Υѥå��������⥸�塼��
  !  * ˰�¾�������Ȥ����ؼ��, �층Ū�˼¹Ԥ���
  !  * ����ȿ���ˤĤ��Ƥ�, ���줾��β���ȿ����˼¹Ԥ���. 

  !���ۤη�����ػ�
  implicit none
  
  !°���λ���
  private

  !�ؿ���°���� public ���ѹ�
  public MoistAdjustSvapPress    !˰�¾��������Ѥ���˰�¼���Ĵ��(�ʰ���)
  public MoistAdjustNH4SH        !����ȿ���ΰ�ʿ��������Ѥ���˰�¼���Ĵ��

contains

!!!------------------------------------------------------------------!!!
  subroutine MoistAdjustSvapPress(xyz_PressDry, xyz_Exner, xyz_PTemp, xyzf_QMix, FactorCloud2Gas)
    !
    ! ˰�¾��������Ѥ�������˰��Ĵ��ˡ�μ¹�
    ! �������ץ����Ǥ�, ͽ���᤿�������ȿ�����ɤ�Ԥ�. 
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : ncmax,             &!���ؼ�ο�
      &                    imin, imax,        &!x ����������ξ�¤Ȳ���
      &                    jmin, jmax,        &!y ����������ξ��
      &                    kmin, kmax          !z ����������β���
    use basicset,   only : xyz_ExnerBZ,       &!̵��������(���ܾ�)
      &                    xyz_PTempBZ,       &!����(���ܾ�)
      &                    xyzf_QMixBZ         !�Ž���ʬ������(���ܾ�)
    use constants,  only : CpDry,             &!������ʬ��ʿ���갵��Ǯ [J/K kg]
      &                    MolWtDry            !������ʬ��ʿ��ʬ����   [kg/mol]
    use composition,only : MolWtWet,          &!������ʬ��ʬ����  
      &                    SpcWetID,          &!������ʬ�β��ؼ��ID 
      &                    CondNum,           &!�ŷ�����ο�
      &                    IdxCG,             &!�ŷ����(����)������ź����
      &                    IdxCC               !�ŷ����(��)������ź����
      use ChemCalc,   only : xyz_SvapPress, xyz_LatentHeat, xyz_DQMixSatDPTemp
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP),intent(in)   :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !������ʬ�ΰ���
    real(DP),intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !�������ʡ��ؿ�
    real(DP),intent(inout):: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !����
    real(DP),intent(inout):: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                          !������
    real(DP),intent(in)   :: FactorCloud2Gas
    integer, parameter    :: ItrNum = 4                   !ȿ�����ɤβ��

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
    integer :: i, s                                       ! ź����  
    integer :: IG, iC

    !---------------------------------------------------------------------
    ! ����˰��Ĵ��ˡ�μ¹�
    !   �롼�פ�󤹤Τ�, ���ˤĤ��Ƥ���.  
    !---------------------------------------------------------------------
    LoopSvapPress: do s = 1, CondNum

       ! ź�����ݴ�
       iG = IdxCG(s)
       iC = IdxCC(s)
      
      !����˰��ˡ�Ǥϰ��Ϥ��Ѳ������ʤ�. 
      xyz_ExnerAll = xyz_Exner + xyz_ExnerBZ
          
      !���ޤǤ�����줿������ʬ���ͤ�����̤Ȥߤʤ�. ź�������ɲ�
      xyz_QMixV_pre  = xyzf_QMix(:,:,:,iG) + xyzf_QMixBZ(:,:,:,iG)
      xyz_QMixC_pre  = xyzf_QMix(:,:,:,iC) + xyzf_QMixBZ(:,:,:,iC)
      xyz_PTemp_pre = xyz_PTemp

      Adjusting: do i = 1, ItrNum
        !---------------------------------------------------------------
        ! ˰�¾���������˰�º���������
        !---------------------------------------------------------------
        !����
        xyz_TempAll = ( xyz_PTemp_pre + xyz_PTempBZ ) * xyz_ExnerAll

        !˰�¾���������˰�º������׻�(���ܾ줫��κ�). 
        xyz_QMixSat =                                    &
          & xyz_SvapPress(SpcWetID(iC), xyz_TempAll)     &
          &  * MolWtWet(iC) / MolWtDry / xyz_PressDry

        !���ʲ����줿��Ǯ
        xyz_Gamma = xyz_LatentHeat(SpcWetID(iC), xyz_TempAll) &
          &        / (xyz_ExnerAll * CpDry)

        !�ŷ��̤����. 
        !  �ŷ뤬��������ˤ�, xz_QMixV_pre - xz_QMixSat ��ɬ�������ͤȤʤ�.
        !  ��ȯ����������ˤ�, ��ȯ�̤� - QMixC ��Ķ���뤳�ȤϤʤ�. 
        xyz_DelQMix =                                                    &
          & ( xyz_QMixV_pre - xyz_QMixSat )                              &
          &   / (1.0d0 + xyz_Gamma * xyz_DQMixSatDPTemp(                 &
          &        SpcWetID(iC), MolWtWet(iC), xyz_TempAll, xyz_ExnerAll & 
          &        ) )

        xyz_Cond = max( 0.0d0, min( xyz_QMixV_pre,   xyz_DelQMix ) )
        xyz_Evap = max( 0.0d0, min( xyz_QMixC_pre, - xyz_DelQMix ) ) * FactorCloud2Gas
        
        !��꿿�˶ᤤ�ͤ�׻�
        xyz_PTemp_nxt  = xyz_PTemp_pre + xyz_Gamma * ( xyz_Cond - xyz_Evap )
        xyz_QMixV_nxt  = xyz_QMixV_pre  - xyz_Cond + xyz_Evap
        xyz_QMixC_nxt  = xyz_QMixC_pre  + xyz_Cond - xyz_Evap

        !�����֤��Τ�����ѿ����
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
    ! NH3 + H2S --> NH4SH ������ȿ���ΰ�ʿ����� Kp ���Ѥ���˰�¼���Ĵ��ˡ
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only : DP
    use dc_message, only : MessageNotify
    use gridset,    only : ncmax,             &!���ؼ�ο�
      &                    imin, imax,        &!x ����������ξ�¤Ȳ���
      &                    jmin, jmax,        &!y ����������ξ��
      &                    kmin, kmax          !z ����������β���
    use basicset,   only : xyz_ExnerBZ,       &!̵��������(���ܾ�)
      &                    xyz_PTempBZ,       &!����(���ܾ�)
      &                    xyzf_QMixBZ         !�Ž���ʬ������(���ܾ�)
    use constants,  only : PressBasis,        &!���̤δ�వ��         [Pa]
      &                    CpDry,             &!������ʬ��ʿ���갵��Ǯ [J/K kg]
      &                    GasRDry             !������ʬ�ε������     [J/K kg]
    use composition,only : MolWtWet,          &!������ʬ��ʬ����  
      &                    IdxNH3,            &!NH3(����)������ź����
      &                    IdxH2S,            &!H2S(����)������ź����
      &                    IdxNH4SHc           !NH4SH(��)������ź����
    use ChemCalc,   only : ReactHeatNH4SH, xyz_DelQMixNH4SH
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP),intent(in)   :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !������ʬ�ΰ���
    real(DP),intent(in)   :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !�������ʡ��ؿ�
    real(DP),intent(inout):: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !����
    real(DP),intent(inout):: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                                          !�Ž���ʬ�κ�����
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
    ! �����
    !---------------------------------------------------------------------

    !����˰��ˡ�Ǥϰ��Ϥ��Ѳ������ʤ�. 
    xyz_ExnerAll = xyz_Exner + xyz_ExnerBZ
    xyz_PressAll = PressBasis * (xyz_ExnerAll ** (CpDry / GasRDry))
        
    !���ޤǤ�����줿������ʬ���ͤ�����̤Ȥߤʤ�. ź�������ɲ�
    xyz_QMixNH3_pre(:,:,:)   = xyzf_QMix(:,:,:,IdxNH3) &
      &                           + xyzf_QMixBZ(:,:,:,IdxNH3)
    xyz_QMixH2S_pre(:,:,:)   = xyzf_QMix(:,:,:,IdxH2S) &
      &                           + xyzf_QMixBZ(:,:,:,IdxH2S)
    xyz_QMixNH4SH_pre(:,:,:) = xyzf_QMix(:,:,:,IdxNH4SHc) &
      &                           + xyzf_QMixBZ(:,:,:,IdxNH4SHc)
    xyz_PTemp_pre            = xyz_PTemp

    !���ʲ����줿ȿ��Ǯ (NH4SH ��kg ���Ф���Ǯ��)
    xyz_Gamma = ReactHeatNH4SH / ( xyz_ExnerAll * CpDry )

    AdjustNH4SH: do i = 1, ItrNum
      !---------------------------------------------------------------
      ! �ѿ��ν����
      !---------------------------------------------------------------
      !����
      xyz_TempAll = ( xyz_PTemp_pre + xyz_PTempBZ ) * xyz_ExnerAll
      
      !NH4SH ��������
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
      ! ��꿿�˶ᤤ�ͤ����˰�¾���������˰�º���������
      !---------------------------------------------------------------
      ! NH4SH �κ��������
      xyz_QMixNH4SH_nxt  = xyz_QMixNH4SH_pre + xyz_Cond - xyz_Evap
      
      ! DelPress �򸵤�, NH3 �� H2S �κ��������
      xyz_QMixNH3_nxt = xyz_QMixNH3_pre - ( xyz_Cond - xyz_Evap )      &
        &                 * MolWtWet(IdxNH3) / MolWtWet(IdxNH4SHc)
      xyz_QMixH2S_nxt = xyz_QMixH2S_pre - ( xyz_Cond - xyz_Evap )      &
        &                 * MolWtWet(IdxH2S) / MolWtWet(IdxNH4SHc)
          
      !���̤���
      xyz_PTemp_nxt = xyz_PTemp_pre + xyz_Gamma * ( xyz_Cond - xyz_Evap )
      
      !�롼�פ�󤹤�����ѿ��Ѳ�
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
