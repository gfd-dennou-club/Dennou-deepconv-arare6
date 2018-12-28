!= Module NumDiffusion
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: numdiffusion.f90,v 1.19 2011-02-28 12:28:38 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module NumDiffusion4th
  !
  !���ͳȻ���η׻��⥸�塼��
  ! 4 �������濴��ʬ������
  !

  !�⥸�塼��ƤӽФ�
  use dc_types,  only : DP

  !���ۤη�����ػ�
  implicit none

  !°���λ���
  private

  !�ؿ��� public �ˤ���
  public NumDiffusion_Init
  public xz_NumDiffScalar
  public xz_NumDiffKm
  public xza_NumDiffScalar
  public pz_NumDiffVelX
  public xr_NumDiffVelZ

  !�ѿ������
  real(DP), save  :: NuH   = 0.0d0   !����Ǵ���η��� (��ʿ����)
  real(DP), save  :: NuV   = 0.0d0   !����Ǵ���η��� (��ľ����)
! real(DP), save  :: Alpha = 1.0d-4
  real(DP), save  :: Alpha = 0.0d0
  real(DP), save  :: Alpha_Velocity = 1.0d0

contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine NumDiffusion_init(cfgfile)
    ! 
    ! NumDiffusion �⥸�塼��ν�����롼����
    !

    !�⥸�塼���ɤ߹���
    use dc_message,  only: MessageNotify
    use mpi_wrapper, only: myrank
    use axesset,     only: dx, dy, dz     !�ʻ����ֳ�
    use timeset,     only: DelTimeLong    !Ĺ�����֥��ƥå�

    !���ۤη�����ػ�
    implicit none

    ! �ѿ�������
    character(*), intent(in) :: cfgfile

    ! Namelist ���� Alpha ���ͤ���ľ��.
    NAMELIST /numdiffusion/ Alpha, Alpha_Velocity
    open (10, FILE=cfgfile)
    read(10, NML=numdiffusion)
    close(10)
    
    ! CReSS �ޥ˥奢��Ǥ�, 2 �������濴��ʬ�ξ��, Alpha < 1/8 ���餤��Ŭ���ȽҤ٤Ƥ���.
    NuH = Alpha * ( dx ** 4.0d0 ) / ( 2.0d0 * DelTimeLong )
    NuV = Alpha * ( dz ** 4.0d0 ) / ( 2.0d0 * DelTimeLong )
    
    !��ǧ
    if (myrank == 0) then
      call MessageNotify( "M", "NumDiffusion_init", "Alpha = %f", d=(/Alpha/) )
      call MessageNotify( "M", "NumDiffusion_init", "NuH = %f", d=(/NuH/) )
      call MessageNotify( "M", "NumDiffusion_init", "NuV = %f", d=(/NuV/) )
      call MessageNotify( "M", "NumDiffusion_init", "Alpha_Velocity = %f", d=(/Alpha_Velocity/) )
      if (Alpha > 0.125) then
        call MessageNotify( "E", "NumDiffusion_init", "Alpha = %f > 1/8", d=(/Alpha/) )
        stop
      end if
    end if

  end subroutine NumDiffusion_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xz_NumDiffScalar(xz_Scalar)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���Ǥο��ͳȻ����ɾ��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: xz_Scalar(imin:imax, kmin:kmax)
                                                    !�����顼��
    real(DP)             :: xz_NumDiffScalar(imin:imax, kmin:kmax)
                                                    !��ʿ�����ο��ͳȻ�
    
    xz_NumDiffScalar =   &
      &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xz_Scalar ))))) &
      &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xz_Scalar ))))) 
    
  end function xz_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xz_NumDiffKm(xz_Scalar)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���Ǥο��ͳȻ����ɾ��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: xz_Scalar(imin:imax, kmin:kmax)
                                                    !�����顼��
    real(DP)             :: xz_NumDiffKm(imin:imax, kmin:kmax)
                                                    !��ʿ�����ο��ͳȻ�
    
    xz_NumDiffKm =   &
      &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xz_Scalar ))))) &
      &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xz_Scalar ))))) 
    
  end function xz_NumDiffKm

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xza_NumDiffScalar(xza_Scalar)
    !
    ! x, z ������Ⱦ�ʻҤ��줿���Ǥο��ͳȻ����ɾ��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax, ncmax
    use differentiate_center2,                &
      &             only: xz_dx_pz, pz_dx_xz, &
      &                   xz_dz_xr, xr_dz_xz

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: xza_Scalar(imin:imax, kmin:kmax, ncmax)
                                                    !�����顼��
    real(DP)             :: xza_NumDiffScalar(imin:imax, kmin:kmax, ncmax)
                                                    !��ʿ�����ο��ͳȻ�
    integer             :: s
    

    do s = 1, ncmax
      xza_NumDiffScalar(:,:,s) =   &
        &  - NuH * (xz_dx_pz(pz_dx_xz(xz_dx_pz(pz_dx_xz( xza_Scalar(:,:,s) ))))) &
        &  - NuV * (xz_dz_xr(xr_dz_xz(xz_dz_xr(xr_dz_xz( xza_Scalar(:,:,s) ))))) 
    end do

  end function xza_NumDiffScalar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function pz_NumDiffVelX(pz_VarX)
    !
    ! z ������Ⱦ�ʻҤ��줿���Ǥο��ͳȻ����ɾ��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: pz_dx_xz, xz_dx_pz, &
      &                   pz_dz_pr, pr_dz_pz

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: pz_VarX(imin:imax, kmin:kmax)
                                                    !ʪ����
    real(DP)             :: pz_NumDiffVelX(imin:imax, kmin:kmax)
                                                    !���ͳȻ�
    
    !®�٤γȻ����������ѹ��������Alpha_Velocity���ͤ�1�ʳ����ͤˤ���
    pz_NumDiffVelX = &
      & -  Alpha_Velocity * NuH * ( pz_dx_xz( xz_dx_pz( pz_dx_xz( xz_dx_pz( pz_VarX ) ) ) ) ) &
      & - Alpha_Velocity * NuV * ( pz_dz_pr( pr_dz_pz( pz_dz_pr( pr_dz_pz( pz_VarX ) ) ) ) )
    
  end function pz_NumDiffVelX
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function xr_NumDiffVelZ(xr_VarZ)
    !
    ! x ������Ⱦ�ʻҤ��줿���Ǥο��ͳȻ����ɾ��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: imin, imax, kmin, kmax
    use differentiate_center2,                &
      &             only: xr_dx_pr, pr_dx_xr, &
      &                   xr_dz_xz, xz_dz_xr

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: xr_VarZ(imin:imax, kmin:kmax)
                                                    !ʪ����
    real(DP)             :: xr_NumDiffVelZ(imin:imax, kmin:kmax)
                                                    !���ͳȻ�

    !®�٤γȻ����������ѹ��������Alpha_Velocity���ͤ�1�ʳ����ͤˤ���
    xr_NumDiffVelZ = &
      & - Alpha_Velocity * NuH * ( xr_dx_pr( pr_dx_xr( xr_dx_pr( pr_dx_xr( xr_VarZ ) ) ) ) )&
      & - Alpha_Velocity * NuV * ( xr_dz_xz( xz_dz_xr( xr_dz_xz( xz_dz_xr( xr_VarZ ) ) ) ) )
    
  end function xr_NumDiffVelZ

end module NumDiffusion4th
