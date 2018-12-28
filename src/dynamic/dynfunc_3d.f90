!= Module DynFunc_3D
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu 
! Version::   $Id: dynfunc_3d.f90,v 1.5 2008-06-26 09:24:27 odakker2 Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
!��ǥ���ϳز�����׻����뤿���ɬ�פȤʤ�ؿ�����«�ͤ��⥸�塼��
!����Ū�ˤϰʲ��ι��׻����뤿��δؿ����Ǽ����.  
!  * ��ή��
!  * ���Ϲ�
!  * ���������Ϲ�
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!  * �������ʡ��ؿ��ζ���������Υ�����ˤ�����, 2 �����٤�Υ�������ۤ����Ѥ��Ƥ��뤿��,
!    ���������Ϲ�η׻��ץ����ˤ����� differentiate_center4 �⥸�塼�����ꤹ��
!    ���ȤϤǤ��ʤ��Τ����.
!
!== Future Plans
!

module DynFunc_3d
  !
  !�۳������Ѥ����ϳز����γƹ�η׻��⥸�塼��. 
  !����Ū�ˤϰʲ��ι��׻����뤿��δؿ����Ǽ����.  
  !  * ��ή��
  !  * ���Ϲ�
  !  * ���������Ϲ�
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP
  
  use gridset, only:  DimXMin,           &! x ����������β���
    &                 DimXMax,           &! x ����������ξ��
    &                 DimYMin,           &! y ����������β���
    &                 DimYMax,           &! y ����������ξ��
    &                 DimZMin,           &! z ����������β���
    &                 DimZMax,           &! z ����������ξ��
    &                 SpcNum              !
  use damping_3d_v2,  only: DampSound           !���Ȥθ��그��
  use basicset_3d, only: xyz_EffMolWtBasicZ, &
    &                    xyz_PotTempBasicZ   !���ܾ�β���
  use constants, only: CpDry,             &!������ʬ����Ǯ
    &                  Grav                !���ϲ�®��
  use xyz_base_module, only: &
    &                xyz_avr_pyz, xyr_avr_pyr, xqz_avr_pqz, &
    &                pyz_avr_xyz, pyr_avr_xyr, pqz_avr_xqz, &
    &                xyz_avr_xqz, pyz_avr_pqz, xyr_avr_xqr, &
    &                xqz_avr_xyz, pqz_avr_pyz, xqr_avr_xyr, &
    &                xyz_avr_xyr, pyz_avr_pyr, xqz_avr_xqr, &
    &                xyr_avr_xyz, pyr_avr_pyz, xqr_avr_xqz 
!  use StorePotTemp_3d,  only: StorePotTempAdv
!  use StoreMixRt_3d,    only: StoreMixRtAdv

  !���ۤη�����ػ�
  implicit none

  !°���λ���
  private

  !��ή��׻��Τ���δؿ��� public �ˤ���
  public xyz_AdvScalar
  public xyz_AdvKm
  public xyza_AdvScalar
  public pyz_AdvVelX
  public xqz_AdvVelY
  public xyr_AdvVelZ

  !���Ϲ�׻��Τ���δؿ��� public �ˤ���
  public xyr_Buoy

  !���������Ϥη׻��Τ���δؿ��� public �ˤ���
  public pyz_GradPi
  public xqz_GradPi


contains


!!!------------------------------------------------------------------------!!!
  function xyz_AdvScalar(xyz_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP), intent(in) :: xyz_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !�����顼��
    real(DP)             :: xyz_AdvScalar &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
    
    xyz_AdvScalar =                                    &
      & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyz_Var))  &
      & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyz_Var))  &
      & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyz_Var))    

!    call StorePotTempAdv( xyz_AdvScalar )   

  end function xyz_AdvScalar


!!!------------------------------------------------------------------------!!!
  function xyz_AdvKm(xyz_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP), intent(in) :: xyz_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !�����顼��
    real(DP)             :: xyz_AdvKm &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
    
    xyz_AdvKm =                               &
      & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyz_Var))  &
      & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyz_Var))  &
      & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyz_Var))    

  end function xyz_AdvKm

!!!------------------------------------------------------------------------!!!
  function xyza_AdvScalar(xyza_Var, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz
!    use xyz_deriv_module, only: pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP), intent(in) :: xyza_Var &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax,SpcNum)
                                                        !�����顼��
    real(DP)             :: xyza_AdvScalar &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax,SpcNum)
                                                        !�����顼�̤ο�ʿ��ή
    integer :: s ! �롼���ѿ�

    do s = 1, SpcNum        
      xyza_AdvScalar(:,:,:,s) =                                  &
        & - xyz_avr_pyz(pyz_VelX * pyz_dx_xyz(xyza_Var(:,:,:,s))) &
        & - xyz_avr_xqz(xqz_VelY * xqz_dy_xyz(xyza_Var(:,:,:,s))) &
        & - xyz_avr_xyr(xyr_VelZ * xyr_dz_xyz(xyza_Var(:,:,:,s)))    
    end do

!    call StoreMixRtAdv( xyza_AdvScalar )   

  end function xyza_AdvScalar
!!!------------------------------------------------------------------------!!!
  function pyz_AdvVelX(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: xyz_dx_pyz, pqz_dy_pyz, pyr_dz_pyz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP)             :: pyz_AdvVelX(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ���ΰ�ή
    
!    pz_AdvVelX = 0.0d0  !�����

    pyz_AdvVelX =                                                 &
      & - pyz_VelX * pyz_avr_xyz( xyz_dx_pyz( pyz_VelX ) )        &
      & - pyz_avr_pqz( pqz_avr_xqz( xqz_VelY ) * pqz_dy_pyz( pyz_VelX ) ) &
      & - pyz_avr_pyr( pyr_avr_xyr( xyr_VelZ ) * pyr_dz_pyz( pyz_VelX ) )

    
  end function pyz_AdvVelX


!!!------------------------------------------------------------------------!!!
  function xyr_AdvVelZ(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: pyr_dx_xyr, xqr_dy_xyr, xyz_dz_xyr
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP)             :: xyr_AdvVelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ���ΰ�ή
  
!    xr_AdvVelZ = 0.0d0  !�����
    xyr_AdvVelZ =                                                 &
      & - xyr_avr_pyr( pyr_avr_pyz( pyz_VelX ) * pyr_dx_xyr( xyr_VelZ ) ) &
      & - xyr_avr_xqr( xqr_avr_xqz( xqz_VelY ) * xqr_dy_xyr( xyr_VelZ ) ) &
      & - xyr_VelZ * xyr_avr_xyz( xyz_dz_xyr( xyr_VelZ ) )
    
  end function xyr_AdvVelZ

!!!------------------------------------------------------------------------!!!
  function xqz_AdvVelY(pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! x ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use xyz_deriv_c4_module, only: pqz_dx_xqz, xyz_dy_xqz, xqr_dz_xqz
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ʿ��®
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ��®
    real(DP)             :: xqz_AdvVelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                                        !��ľ���ΰ�ή
  
!    xr_AdvVelZ = 0.0d0  !�����
    xqz_AdvVelY =                                                 &
      & - xqz_avr_pqz( pqz_avr_pyz( pyz_VelX ) * pqz_dx_xqz( xqz_VelY ) ) &
      & - xqz_VelY * xqz_avr_xyz( xyz_dy_xqz( xqz_VelY ) ) &
      & - xqz_avr_xqr( xqr_avr_xyr( xyr_VelZ ) * xqr_dz_xqz( xqz_VelY ) )
    
  end function xqz_AdvVelY

 
!!!------------------------------------------------------------------------!!!
  function xyr_Buoy(xyz_PotTemp)
    !
    ! ��ľ�����α�ư�������˸�������Ϲ��׻�

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in)  :: xyz_PotTemp(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !���̾���
    real(DP)              :: xyr_Buoy &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !���Ϲ�

!    !�����
!    xr_Buoy = 0.0d0

    !���Ϲ�η׻�
    xyr_Buoy = Grav * xyr_avr_xyz(xyz_PotTemp / xyz_PotTempBasicZ)

  end function xyr_Buoy


!!!------------------------------------------------------------------------!!!
  function pyz_GradPi(xyz_Exner, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! z ������Ⱦ�ʻҤ��줿���Ǥΰ��Ϸ����Ϲ�η׻�. 
    ! ���ȸ�����ޤ᤿�������꼰�����Ƥ��뤳�Ȥ����.
    !
    
    !�⥸�塼���ɤ߹���

    use xyz_deriv_module, only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, pyz_dx_xyz
        
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: xyz_Exner &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !�������ʡ��ؿ��ξ���
    real(DP), intent(in) :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ʿ®��
    real(DP), intent(in) :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ʿ®��
    real(DP), intent(in) :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ľ®��
    real(DP)             :: pyz_GradPi &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !���Ϸ�����
    real(DP)             :: xyz_DivVel &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !®�٤μ�«

    !®�٤μ�«
    xyz_DivVel = xyz_dx_pyz( pyz_VelX ) + xyz_dy_xqz( xqz_VelY ) &
      &          + xyz_dz_xyr( xyr_VelZ )
    
    !���Ϸ���
!    pyz_GradPi = 0.0d0
    pyz_GradPi =  &
      & + pyz_avr_xyz( CpDry * xyz_PotTempBasicZ / xyz_EffMolWtBasicZ )    &
      &   * pyz_dx_xyz( xyz_Exner ) &
      & - pyz_dx_xyz( DampSound * xyz_DivVel )

  end function pyz_GradPi

!!!------------------------------------------------------------------------!!!
  function xqz_GradPi(xyz_Exner, pyz_VelX, xqz_VelY, xyr_VelZ)
    !
    ! y ������Ⱦ�ʻҤ��줿���Ǥΰ��Ϸ����Ϲ�η׻�. 
    ! ���ȸ�����ޤ᤿�������꼰�����Ƥ��뤳�Ȥ����.
    !
    
    !�⥸�塼���ɤ߹���

    use xyz_deriv_module, only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, xqz_dy_xyz
        
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in)  :: xyz_Exner &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !�������ʡ��ؿ��ξ���
    real(DP), intent(in)  :: pyz_VelX &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ʿ®��
    real(DP), intent(in)  :: xqz_VelY &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ʿ®��
    real(DP), intent(in)  :: xyr_VelZ &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !��ľ®��
    real(DP)              :: xqz_GradPi &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !���Ϸ�����
    real(DP)              :: xyz_DivVel &
      &                    (DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
                                               !®�٤μ�«

    !®�٤μ�«
    xyz_DivVel = xyz_dx_pyz( pyz_VelX ) + xyz_dy_xqz( xqz_VelY ) &
      &          + xyz_dz_xyr( xyr_VelZ )
    
    !���Ϸ���
!    xqz_GradPi = 0.0d0
    xqz_GradPi =  &
      & + xqz_avr_xyz( CpDry * xyz_PotTempBasicZ / xyz_EffMolWtBasicZ )    &
      &   * xqz_dy_xyz( xyz_Exner )               &
      & - xqz_dy_xyz( DampSound * xyz_DivVel )
    
  end function xqz_GradPi
  
end module DynFunc_3d
