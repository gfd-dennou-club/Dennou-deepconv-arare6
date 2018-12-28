!= Module BasicEnvInit
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_basic.f90,v 1.9 2011/06/17 19:07:25 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
!�ǥե���Ȥδ��ܾ�����ꤹ�뤿����ѿ����ȷ��⥸�塼��
!   * BasicEnvFile_init: ���ܾ���ͤ� netCDF �ե����뤫�����
!   * BasicEnvCalc_Init: ���ܾ�ξ���� Namelist ������������ͤ�׻�
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!== Future Plans
!
module initialdata_basic
  !
  ! �ǥե���Ȥδ��ܾ���뤿��Υ롼����

    
  !���ۤη�����ػ�
  implicit none

  !�ǥե���Ȥ� private
  private

  real(8), save  :: Humidity = 0.0d0        !���ܾ�μ��� 
  real(8), save  :: TempStrat = 100.0d0     !���ط��β��� [k]
  real(8), save  :: HeightStrat = 18.0d3    !���ط��ι��� [k]
  real(8), save  :: Dhight = 5.0d3          !�Ťߴؿ��Υѥ�᡼�� [m]

  ! public
  public initialdata_basic_nml
  public initialdata_basic_dry
  public initialdata_basic_moist
  public initialdata_basic_strat
  public initialdata_basic_isothermal
    
contains

  subroutine initialdata_basic_dry( z_TempBZ, z_PressBZ, za_MolFr )

    implicit none

    real(8), intent(out) :: z_TempBZ(kmin:kmax)
    real(8), intent(out) :: z_PressBZ(kmin:kmax)
    real(8), intent(out) :: za_MolFr(kmin:kmax, 1:ncmax)
    


  end subroutine initialdata_basic_dry


  subroutine initialdata_basic_moist( z_TempBZ, z_PressBZ, za_MolFr )

    implicit none

    real(8), intent(out) :: z_TempBZ(kmin:kmax)
    real(8), intent(out) :: z_PressBZ(kmin:kmax)
    real(8), intent(out) :: za_MolFr(kmin:kmax, 1:ncmax)
    
    
    
  end subroutine initialdata_basic_moist

  

 
end module initialdata_basic
