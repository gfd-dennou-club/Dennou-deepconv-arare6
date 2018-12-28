!= Module DynFunc
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: dynfunc.f90,v 1.12 2011-02-28 12:00:23 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
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
!    ���������Ϲ�η׻��ץ����ˤ����� differentiate_center4 �⥸�塼�����ꤹ�뤳�ȤϤǤ��ʤ��Τ����.
!
!== Future Plans
!

module DynFunc
  !
  !�۳������Ѥ����ϳز����γƹ�η׻��⥸�塼��. 
  !����Ū�ˤϰʲ��ι��׻����뤿��δؿ����Ǽ����.  
  !  * ��ή��
  !  * ���Ϲ�
  !  * ���������Ϲ�
  !

  !�⥸�塼���ɤ߹���
  use gridset, only:  DimXMin,           &! x ����������β���
    &                 DimXMax,           &! x ����������ξ��
    &                 DimZMin,           &! z ����������β���
    &                 DimZMax,           &! z ����������ξ��
    &                 SpcNum              !
  use damping,  only: DampSound           !���Ȥθ��그��
  use basicset, only: xz_PotTempBasicZ,  &!���ܾ�β���
    &                 xz_EffMolWtBasicZ   !���ܾ��ʬ���̸���
  use constants,only: CpDry,             &!������ʬ����Ǯ
    &                 Grav                !���ϲ�®��
  use average,  only: xz_avr_pz, xz_avr_xr, &
    &                 pz_avr_xz, pz_avr_pr, &
    &                 pr_avr_xr, pr_avr_pz, &
    &                 xr_avr_pr, xr_avr_xz
!  use StorePotTemp,  only: StorePotTempAdv
!  use StoreMixRt,    only: StoreMixRtAdv
!  use StoreMom,    only: StoreMomAdv
!  use StoreBuoy,    only: StoreBuoyTemp

  !���ۤη�����ػ�
  implicit none

  !°���λ���
  private

  !��ή��׻��Τ���δؿ��� public �ˤ���
  public xz_AdvScalar
  public xz_AdvKm
  public xza_AdvScalar
  public pz_AdvVelX
  public xr_AdvVelZ

  !���Ϲ�׻��Τ���δؿ��� public �ˤ���
  public xr_Buoy

  !���������Ϥη׻��Τ���δؿ��� public �ˤ���
  public pz_GradPi


contains


!!!------------------------------------------------------------------------!!!
  function xz_AdvScalar(xz_Var, pz_VelX, xr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ʿ��®
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ľ��®
    real(8), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼��
    real(8)             :: xz_AdvScalar(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
    
    xz_AdvScalar =                               &
      & - xz_avr_pz(pz_VelX * pz_dx_xz(xz_Var))  &
      & - xz_avr_xr(xr_VelZ * xr_dz_xz(xz_Var))    

!    call StorePotTempAdv( xz_AdvScalar )   

  end function xz_AdvScalar



!!!------------------------------------------------------------------------!!!
  function xz_AdvKm(xz_Var, pz_VelX, xr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ʿ��®
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ľ��®
    real(8), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼��
    real(8)             :: xz_AdvKm(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
    
    xz_AdvKm =                                   &
      & - xz_avr_pz(pz_VelX * pz_dx_xz(xz_Var))  &
      & - xz_avr_xr(xr_VelZ * xr_dz_xz(xz_Var))    
    
  end function xz_AdvKm



  function xza_AdvScalar(xza_Var, pz_VelX, xr_VelZ)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center4, only: pz_dx_xz, xr_dz_xz
!    use differentiate_center2, only: pz_dx_xz, xr_dz_xz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ʿ��®
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ľ��®
    real(8), intent(in) :: xza_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
                                                        !�����顼��
    real(8)             :: xza_AdvScalar(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
                                                        !�����顼�̤ο�ʿ��ή
    integer             :: s

    do s = 1, SpcNum
      xza_AdvScalar(:,:,s) =                               &
        & - xz_avr_pz(pz_VelX * pz_dx_xz(xza_Var(:,:,s)))  &
        & - xz_avr_xr(xr_VelZ * xr_dz_xz(xza_Var(:,:,s)))    
    end do

!    call StoreMixRtAdv( xza_AdvScalar )   
    
  end function xza_AdvScalar
  

!!!------------------------------------------------------------------------!!!
  function pz_AdvVelX(pz_VelX, xr_VelZ)
    !
    ! z ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center4, only: xz_dx_pz, pr_dz_pz
!    use differentiate_center2, only: xz_dx_pz, pr_dz_pz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ʿ��®
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ľ��®
    real(8)             :: pz_AdvVelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
    
!    pz_AdvVelX = 0.0d0  !�����
    pz_AdvVelX =                                                 &
      & - pz_VelX * pz_avr_xz( xz_dx_pz( pz_VelX ) )             &
      & - pz_avr_pr( pr_avr_xr( xr_VelZ ) * pr_dz_pz( pz_VelX ) )
    
!    call StoreMomAdv( pz_AdvVelX )   

  end function pz_AdvVelX


!!!------------------------------------------------------------------------!!!
  function xr_AdvVelZ(xr_VelZ, pz_VelX)
    !
    ! x ������Ⱦ�ʻҤ��줿���ˤ������ή��׻�
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center4, only: pr_dx_xr, xz_dz_xr
!    use differentiate_center2, only: pr_dx_xr, xz_dz_xr
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ʿ��®
    real(8), intent(in) :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !��ľ��®
    real(8)             :: xr_AdvVelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                                        !�����顼�̤ο�ʿ��ή
  
!    xr_AdvVelZ = 0.0d0  !�����
    xr_AdvVelZ =                                                 &
      & - xr_avr_pr( pr_avr_pz( pz_VelX ) * pr_dx_xr( xr_VelZ ) ) &
      & - xr_VelZ * xr_avr_xz( xz_dz_xr( xr_VelZ ) )
    
  end function xr_AdvVelZ
  

!!!------------------------------------------------------------------------!!!
  function xr_Buoy(xz_PotTemp)
    !
    ! ��ľ�����α�ư�������˸�������Ϲ��׻�
    !
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)  :: xz_PotTemp(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !���̾���
    real(8)              :: xr_Buoy(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !���Ϲ�

!    !�����
!    xr_Buoy = 0.0d0

    !���Ϲ�η׻�
    xr_Buoy = Grav * xr_avr_xz(xz_PotTemp / xz_PotTempBasicZ)

!    call StoreBuoyTemp(xz_avr_xr(xr_Buoy))

  end function xr_Buoy


!!!------------------------------------------------------------------------!!!
  function pz_GradPi(xz_Exner, pz_VelX, xr_VelZ)
    !
    ! z ������Ⱦ�ʻҤ��줿���Ǥΰ��Ϸ����Ϲ�η׻�. 
    ! ���ȸ�����ޤ᤿�������꼰�����Ƥ��뤳�Ȥ����.
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center2,  only: pz_dx_xz, xz_dx_pz, xz_dz_xr
        
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)  :: xz_Exner(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !�������ʡ��ؿ��ξ���
    real(8), intent(in)  :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !��ʿ®��
    real(8), intent(in)  :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !��ľ®��
    real(8)              :: pz_GradPi(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !���Ϸ�����
    real(8)              :: xz_DivVel(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !®�٤μ�«

    !®�٤μ�«
    xz_DivVel =  xz_dx_pz( pz_VelX ) + xz_dz_xr( xr_VelZ )
    
    !���Ϸ���
!    pz_GradPi = 0.0d0
!    pz_GradPi =  &
!      & pz_avr_xz( CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )    &
!      &   * ( pz_dx_xz( xz_Exner ) - pz_dx_xz( DampSound * xz_DivVel ) )  
    pz_GradPi =  &
      & pz_avr_xz( CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )    &
      &   * pz_dx_xz( xz_Exner ) &
      & - pz_dx_xz( DampSound * xz_DivVel ) 
    
  end function pz_GradPi

!!!--------------------------------------------------------------------!!!
  function xr_GradPi(xz_Exner, pz_VelX, xr_VelZ)
    !
    ! x ������Ⱦ�ʻҤ��줿���Ǥΰ��Ϸ����Ϲ�η׻�. 
    ! ���ȸ�����ޤ᤿�������꼰�����Ƥ��뤳�Ȥ����.
    !
    
    !�⥸�塼���ɤ߹���
    use differentiate_center2,  only: xr_dz_xz, xz_dx_pz, xz_dz_xr

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)  :: xz_Exner(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !�������ʡ��ؿ��ξ���
    real(8), intent(in)  :: pz_VelX(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !��ʿ®��
    real(8), intent(in)  :: xr_VelZ(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !��ľ®��
    real(8)              :: xr_GradPi(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !���Ϸ�����
    real(8)              :: xz_DivVel(DimXMin:DimXMax, DimZMin:DimZMax)
                                               !®�٤μ�«

    !®�٤μ�«
    xz_DivVel =  xz_dx_pz( pz_VelX ) + xz_dz_xr( xr_VelZ )
    
    !®�� w �ΰ��ϸ���
!    xr_GradPi = 0.0d0
!    xr_GradPi =   &
!      & xr_avr_xz(CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )   &
!      &   * ( xr_dz_xz( xz_Exner ) - xr_dz_xz( DampSound * xz_DivVel ) )
    xr_GradPi =   &
      & xr_avr_xz(CpDry * xz_PotTempBasicZ / xz_EffMolWtBasicZ )   &
      &   * xr_dz_xz( xz_Exner ) &
      & - xr_dz_xz( DampSound * xz_DivVel ) 
    
  end function xr_GradPi
  
end module DynFunc
