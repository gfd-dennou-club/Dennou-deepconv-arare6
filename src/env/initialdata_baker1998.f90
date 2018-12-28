!= Module initialdata_baker1998
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_baker1998.f90,v 1.2 2014/03/04 05:55:04 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_baker1998
  !
  ! Baker et al. (1998) ���Ϥ���������������뤿��Υ⥸�塼��

  !���ۤη�����ػ�
  implicit none

  !�ǥե���Ȥ� private
  private

  !�������������
  public initialdata_baker1998_basic

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_baker1998_basic( z_Temp, z_Press )
    !
    !== ����
    !  * Baker and Shubert (1989) �ǻȤ�줿����ͤ�Ƹ�����.
    !    * �絤�����٤��鲹�١����Ϥδ��ܾ����
    !      * ���̤Υǡ����ϰ����٤Ǽ�����Ƥ���Τ�, ����򲹰̤ǾƤ�ľ��
    !      * ���̤���, ���ܾ�β��١����Ϥ�׻�����. 
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,  only: STRING, DP
    use gridset,   only: kmin,          &!����� Z �����β���
      &                  kmax,          &!����� Z �����ξ��
      &                  nz              !����� Z �����ο�
    use axesset,   only: r_Z,           &!����
      &                  dz              !��ľ�ʻ����ֳ�
    use constants, only: GasRDry,       &!������ʬ���갵��Ǯ
      &                  TempTop,       &!��ɽ�̲���
!      &                  PressTop,      &!��ɽ�̰���
      &                  PressBasis,    &!��ɽ�̰���
      &                  CpDry,         &!������Ǯ��Ψ
      &                  Grav            !����

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(out) :: z_Press(kmin:kmax)!����
    real(DP), intent(out) :: z_Temp(kmin:kmax) !����
    real(DP)              :: r_DPTempDz(kmin:kmax)
    real(DP)              :: z_PTemp(kmin:kmax)
    real(DP), parameter   :: r_dPTempDz_45km = 3.9d-3
    real(DP), parameter   :: r_dPTempDz_48km = 0.0d0
    real(DP), parameter   :: r_dPTempDz_55km = 0.0d0
    real(DP), parameter   :: r_dPTempDz_60km = 8.35d-3
    integer               :: k

    ! �����
    !
    z_Press  = 0.0d0
    z_Temp   = 0.0d0

    ! �絤�����٤�Ϳ����.
    !   BS1998 ��ɽ 1 ���ܤ��ɤ߼��, �����ľ��������Ƥ���. 
    ! 
    do k = kmin, kmax
      if ( r_Z(k) < 45.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_45km + 0.30d-6 * ( r_Z(k) - 45.0d3 )
      elseif ( 45.0d3 <= r_Z(k) .AND. r_Z(k) < 48.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_45km - 1.30d-6 * ( r_Z(k) - 45.0d3 )
      elseif ( 48.0d3 <= r_Z(k) .AND. r_Z(k) < 55.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_48km
      elseif ( 55.0d3 <= r_Z(k) .AND. r_Z(k) < 60.0d3 ) then 
        r_DPTempDz(k) = r_DPTempDz_55km + 1.67d-6 * ( r_Z(k) - 55.0d3 )
      elseif ( 60.0d3 <= r_Z(k) ) then 
        r_DPTempDz(k) = r_DPTempDz_60km - 1.20d-6 * ( r_Z(k) - 60.0d3 )
      end if
    end do

    ! ���̸��� DPTempDz ����ʬ�����پ� z_PTemp ��׻�
    ! ������Ⱦ�����ʻ������ͤ��Ѥ���. 
    !
    z_PTemp(nz) = TempTop - r_DPTempDz(1) * dz * 0.5d0
    do k = nz-1, 1, -1
      z_PTemp(k) = z_PTemp(k+1) - r_DPTempDz(k) * dz
    end do

    !----------------------------------------------
    ! ���ϤȲ��٤η׻�
    !    
    z_Press(nz) = 2.20d4
    z_Temp(nz)  = 268.0d0
!    z_Press(nz) = PressTop + (Grav * PressTop * dz * 5.0d-1) / (GasRDry * TempTop)
!    z_Temp(nz)  = z_PTemp(nz) * (z_Press(nz) / PressBasis) ** (GasRDry / CpDry)

    do k = nz-1, 1, -1
      z_Press(k) = z_Press(k+1) + (Grav * z_Press(k+1) * dz) / (GasRDry * z_Temp(k+1))
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressBasis) ** (GasRDry / CpDry))
    end do
       
  end subroutine Initialdata_baker1998_basic
  
end module Initialdata_baker1998
