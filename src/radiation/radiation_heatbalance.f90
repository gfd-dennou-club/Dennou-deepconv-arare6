!= ��ñ���ͥ⥸�塼��: �Ȥ�������ΰ�������ѡ���Ǯ����
!
! Authors::   YAMASHITA Tatsuya, SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_heatbalance.f90,v 1.13 2014/07/08 01:04:18 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_HeatBalance
  !
  ! ��ñ���ͥ⥸�塼��: �Ȥ�������ΰ�������ѡ���Ǯ����
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP, STRING
  
  !���ۤη�����ػ�
  implicit none

  !private °��
  private

  !�ѿ����
  real(DP),save   :: RadCoolRate = 0.0d0    !�������Ͳ�ǮΨ [K/day]
  integer, save   :: IdxHeatUp   = 0        !��Ǯ���¤α�ľ��ɸ���б�����������
  integer, save   :: IdxHeatDown = 0        !��Ѱ��¤α�ľ��ɸ���б�����������
  integer, save   :: IdxCoolUp   = 0        !��Ǯ���¤α�ľ��ɸ���б�����������
  integer, save   :: IdxCoolDown = 0        !��Ѱ��¤α�ľ��ɸ���б�����������
  real(DP),save, allocatable :: xyz_RadHeightHeat(:,:,:)  !���Ͳ�Ǯ��¸�ߤ����ΰ�
  real(DP),save, allocatable :: xyz_RadHeightCool(:,:,:)  !���Ͳ�Ǯ��¸�ߤ����ΰ�

  real(DP), save  :: FactorDExnerDtRad = 1.0d0

  character(STRING), parameter:: module_name = 'Radiation_Heatbalance'

  !�ؿ��� public �ˤ���. 
  public Radiation_heatbalance_init
  public Radiation_heatbalance_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_heatbalance_init
    !
    !NAMELIST �������Ψ, ��Ѥ�������ΰ�, ��Ǯ��������ΰ�, �����ꤹ��. 
    !��ǮΨ���襹�ƥå׷׻�����Τ�, ������롼�������ǤϷ��ʤ�. 
    !

    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use namelist_util,     only : namelist_filename
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax        !z ����������ξ��
    use axesset,           only : z_Z         !Z ��ɸ��(�����顼�ʻ���)
    
    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(8) :: HeightHeatUp   = 0.0d0  !��Ǯ�ΰ�ξ�ü�ι���
    real(8) :: HeightHeatDown = 0.0d0  !��Ǯ�ΰ�β�ü�ι���
    real(8) :: HeightCoolUp   = 0.0d0  !����ΰ�ξ�ü�ι���
    real(8) :: HeightCoolDown = 0.0d0  !����ΰ�β�ü�ι���
    integer :: k                       !�롼���ѿ�
    integer :: unit

    ! NAMELIST �����������
    NAMELIST /radiation_heatbalance_nml/ &
      & RadCoolRate, HeightHeatUp, HeightHeatDown, HeightCoolUp, HeightCoolDown, &
      & FactorDExnerDtRad

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_heatbalance_nml)
    close(unit)

    do k = kmin, kmax-1
      if ( z_Z(k) <= HeightHeatUp .AND. HeightHeatUp < z_Z(k+1) ) then 
        IdxHeatUp = k
      elseif ( z_Z(k) <= HeightHeatDown .AND. HeightHeatDown < z_Z(k+1) ) then 
        IdxHeatDown = k
      elseif ( z_Z(k) <= HeightCoolUp .AND. HeightCoolUp < z_Z(k+1) ) then 
        IdxCoolUp = k
      elseif ( z_Z(k) <= HeightCoolDown .AND. HeightCoolDown < z_Z(k+1) ) then 
        IdxCoolDown = k
      end if
    end do

    allocate( xyz_RadHeightHeat(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_RadHeightCool(imin:imax, jmin:jmax, kmin:kmax) )

    xyz_RadHeightHeat = 1.0d0
    xyz_RadHeightHeat(:,:,IdxHeatDown:IdxHeatUp) = 1.0d0
    xyz_RadHeightCool = 0.0d0
    xyz_RadHeightCool(:,:,IdxCoolDown:IdxCoolUp) = 1.0d0

    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "RadCoolRate = %f", d=(/RadCoolRate/))
    call MessageNotify( "M", &
      & module_name, "HeightHeatUp = %f", d=(/HeightHeatUp/))
    call MessageNotify( "M", &
      & module_name, "HeightHeatDown= %f", d=(/HeightHeatDown/))
    call MessageNotify( "M", &
      & module_name, "HeightCoolUp = %f", d=(/HeightCoolUp/))
    call MessageNotify( "M", &
      & module_name, "HeightCoolDown= %f", d=(/HeightCoolDown/))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtRad= %f", d=(/ FactorDExnerDtRad /))

    ! �ҥ��ȥ�ǡ������
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1"',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of exner function', &
      & units='s-1"',    &
      & xtype='float')

  end subroutine Radiation_heatbalance_init

!!!------------------------------------------------------------------------!!!
  subroutine radiation_heatbalance_forcing(xyz_Exner, xyz_PTemp, xyz_DPTempDt, xyz_DExnerDt)
    !
    ! ���̤����Ͷ�����. 
    ! ��ɽ�̤��� RadHeight �ǻ��ꤵ�줿���٤ޤǤδ֤Ƕ���Ū�˰��ͤʲ�Ǯ, 
    ! RadHeight ���� RadHeight2 �ޤǤδ֤Ƕ���Ū�ʰ��ͤ���Ѥ�Ϳ����. 
    ! ���Ͷ������ΤȤ��Ʋ�Ǯ����Ѥ��Х�󥹤���褦�˻�������Ͳ�ǮΨ
    ! ����Ӱ������Ψ���Ѳ�������. 
    ! ��Ѥο�����Ϳ��, ��Ǯ�ο�����Ĵ�᤹��. 

    ! �⥸�塼��ƤӽФ�
    !
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax,     & !z ����������ξ��
      &                           nx, ny, nz
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use axesset,           only : z_dz  
    
    use constants,         only : DayTime,  & ! 1 ����Ĺ�� [s]
      &                           PressBasis, & !���̤δ�వ��
      &                           CpDry,      & !�갵��Ǯ
      &                           CvDry,      & !������Ǯ
      &                           GasRDry       !�������
    use basicset,          only : xyz_ExnerBZ,  &!�������ʡ��ؿ��δ��ܾ�
      &                           xyz_PTempBZ    !���̤δ��ܾ�


    !���ۤη������ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: xyz_Exner(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)  :: xyz_PTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DPTempDt0(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DExnerDt0(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Rad(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_RadPI(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_DensSum(imin:imax, jmin:jmax, kmin:kmax) 
                                        !̩��(���ܾ���ʬ�Ⱦ�����ʬ����)
    real(DP)              :: HeatSum    !�������Ǥ�ñ�̻���������β�Ǯ��
    real(DP)              :: CoolSum    !�������Ǥ�ñ�̻���������������
    real(DP)              :: RadHeatRate
    integer               :: i, j, k

    ! �����
    !
    HeatSum = 0.0d0
    CoolSum = 0.0d0
    xyz_DPTempDt0 = xyz_DPTempDt
    xyz_DExnerDt0 = xyz_DExnerDt

    ! ̩�٤λ���
    xyz_DensSum = PressBasis                            &
      & * (xyz_ExnerBZ + xyz_Exner)**( CvDry /GasRDry ) &
      & / (GasRDry * (xyz_PTempBZ + xyz_PTemp ) )

    ! �������Ǥ� (��Ǯ��/��ǮΨ) ��׻�
    do k = IdxHeatDown, IdxHeatUp
      do j = 1, ny
        do i = 1, nx
          HeatSum = HeatSum + z_dz(k) * CpDry * xyz_DensSum(i,j,k)
        end do
      end do
    end do
    
    ! �������Ǥ� (�����/���Ψ) ��׻�
    do k = IdxCoolDown, IdxCoolUp
      do j = 1, ny
        do i = 1, nx
          CoolSum = CoolSum + z_dz(k) * CpDry * xyz_DensSum(i,j,k)
        end do
      end do
    end do

    ! ��ǮΨ�λ���
    ! RadCoolRate �����ͤǤ��뤳�Ȥ����. 
    RadHeatRate = - RadCoolRate * CoolSum / HeatSum
    
    xyz_Rad = &
      &   xyz_RadHeightHeat * RadHeatRate / ( xyz_ExnerBZ  * DayTime ) &
      & + xyz_RadHeightCool * RadCoolRate / ( xyz_ExnerBZ  * DayTime )

    xyz_DPTempDt = xyz_DPTempDt0 + xyz_Rad

    ! �����Ѳ�
    !
    xyz_RadPI = xyz_DExnerDt_xyz( xyz_Rad ) * FactorDExnerDtRad

    ! output
    !
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_Rad(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_RadPI(1:nx,1:ny,1:nz))
   
  end subroutine radiation_heatbalance_forcing
  
end module Radiation_HeatBalance
