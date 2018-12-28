!= ��ñ����: �Ȥ�����٤�������
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_simple.f90,v 1.10 2014/03/04 04:49:41 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_Simple
  !
  ! ��ñ����: �Ȥ�����٤�������

  !�⥸�塼���ɤ߹���
  !
  use dc_types, only: DP, STRING
  
  !���ۤη�����ػ�
  !
  implicit none

  !private °��
  !
  private

  !�ѿ����
  real(DP), save, allocatable, public :: xyz_DPTempDtRadVary(:,:,:)   
                                           !���Ͳ�Ǯ��
  real(DP), save, allocatable, public :: xyz_DPTempDtRadConst(:,:,:)  
                                           !���Ͳ�Ǯ��
  real(DP), save, allocatable, public :: xyz_ExnerRadVary(:,:,:)   
                                           !��������������Ǯ��Ǯ��
  real(DP), save, allocatable, public :: xyz_ExnerRadConst(:,:,:) 
                                           !��������������Ǯ��Ǯ�� 

  real(DP), save               :: FactorDExnerDtRad = 1.0d0  
  character(STRING), parameter :: module_name = 'radiation_simple'
                                           ! �⥸�塼���̾��.
                                           ! Module name

  public Radiation_Simple_init
  public Radiation_HeatConst_forcing
  public Radiation_HeatVary_forcing

contains

!!!------------------------------------------------------------------------!!!
  subroutine Radiation_Simple_init
    !
    !NAMELIST �������Ψ, �����ΰ������.
    !

    ! �⥸�塼����ɤ߹���
    !
    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax        !z ����������ξ��
    use axesset,           only : z_Z         !Z ��ɸ��(�����顼�ʻ���)
    use constants,         only : DayTime     ! 1 ����Ĺ�� [s]
    use basicset,          only : xyz_ExnerBZ !�������ʡ��ؿ��δ��ܾ�
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use namelist_util,     only : namelist_filename
    
    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP) :: HeightUp   = 0.0d0  !���Ͷ�����Ϳ�����ľ�ΰ�ξ��
    real(DP) :: HeightDown = 0.0d0  !���Ͷ�����Ϳ�����ľ�ΰ�β���
    real(DP) :: RadHeatRate = 0.0d0 !�������Ͳ�ǮΨ [K/day]
    integer  :: k                   !�롼���ѿ�
    integer  :: unit

    ! NAMELIST �����������
    NAMELIST /radiation_simple_nml/ &
      & RadHeatRate, HeightUp, HeightDown, &
      & FactorDExnerDtRad

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_simple_nml)
    close(unit)

    allocate( xyz_DPTempDtRadVary(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_DPTempDtRadConst(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadVary(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadConst(imin:imax, jmin:jmax, kmin:kmax) )

    
    ! ���̤μ������Ͷ�����. 
    ! ��ɽ�̤��� RadHeight �ǻ��ꤵ�줿���٤ޤǤδ֤ǰ���������Ѥ�Ϳ����. 
    !
    do k = kmin, kmax
      if ( z_Z(k) <= HeightDown  ) then
        xyz_DPTempDtRadConst(:,:,k) = 0.0d0 
      elseif( z_Z(k) >= HeightUp ) then
        xyz_DPTempDtRadConst(:,:,k) = 0.0d0 
      else
        xyz_DPTempDtRadConst(:,:,k) = RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k)
      end if
    end do    
    xyz_ExnerRadConst = xyz_DExnerDt_xyz(xyz_DPTempDtRadConst) * FactorDExnerDtRad
    

    ! ���̤����Ͷ�����.
    ! ��ɽ�̤��� HeightDown �ޤǤ� RadHeatRate �����. 
    ! HeightUp ������ϲ�ǮΨ����ˤʤ�褦�˲�ǮΨ�򸺾�������.
    ! Nakajima and Matsuno(1988),����(1944)�򻲹ͤˤ���
    !
    do k = kmin, kmax
      if ( z_Z(k) <= HeightDown ) then
        xyz_DPTempDtRadVary(:,:,k) = RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k)
        
      elseif ( z_Z(k) > HeightDown .AND. z_Z(k) <= HeightUP) then
        xyz_DPTempDtRadVary(:,:,k) =  &
          & RadHeatRate / DayTime / xyz_ExnerBZ(:,:,k) &
          & * (HeightUP - z_Z(k)) / (HeightUP - HeightDown) 
        
      else if (z_Z(k) > HeightUP) then
        xyz_DPTempDtRadVary(:,:,k) = 0.0d0
      end if
    end do
    xyz_ExnerRadVary = xyz_DExnerDt_xyz(xyz_DPTempDtRadVary) * FactorDExnerDtRad

    ! Output
    !
    call MessageNotify( "M", &
      & module_name, "RadHeatRate = %f", d=(/RadHeatRate/))
    call MessageNotify( "M", &
      & module_name, "HeightUp = %f", d=(/HeightUP/))
    call MessageNotify( "M", &
      & module_name, "HeightDown= %f", d=(/HeightDown/))
    call MessageNotify( "M", &
      & module_name, "FactorDExnerDtRad= %f", d=(/ FactorDExnerDtRad /))

    ! �ҥ��ȥ�ǡ������
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='float')

  end subroutine Radiation_Simple_init


  subroutine Radiation_HeatConst_forcing(  &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! �⥸�塼���ɤ߹���
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax,     & !z ����������ξ��
      &                           nx, ny, nz

    ! ���ۤη�����ػ�
    !
    implicit none

    ! �������ѿ�
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency �ι���
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadConst
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadConst

    ! �ե��������
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadConst(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadConst(1:nx,1:ny,1:nz))   

  end subroutine Radiation_HeatConst_forcing


  subroutine Radiation_HeatVary_forcing(   &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! �⥸�塼���ɤ߹���
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax,     & !z ����������ξ��
      &                           nx, ny, nz
    
    ! ���ۤη�����ػ�
    !
    implicit none

    ! �������ѿ�
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency �ι���
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadVary
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadVary
    
    ! �ե��������
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadVary(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadVary(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_HeatVary_forcing


end module Radiation_Simple
