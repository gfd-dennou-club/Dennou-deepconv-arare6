!= ����⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: constants.f90,v 1.12 2014/01/21 05:00:57 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module constants
  !
  != ����⥸�塼��
  !

  !�⥸�塼���ɤ߹���
  !
  use dc_types,      only: DP, STRING

  !���ۤη�����ػ�
  !
  implicit none

  ! private °��
  !
  private

  !Public Interface
  !
  real(DP), save, public :: Grav = 9.8d0          !���� [m/s^2]
  real(DP), save, public :: PressBasis = 965.0d0  !���̤δ�వ�� [Pa]
  real(DP), save, public :: TempSfc = 0.0d0       !��ɽ�̲��� [K]
  real(DP), save, public :: PressSfc = 0.0d0      !��ɽ�̰��� [Pa]
  real(DP), save, public :: TempTop = 0.0d0       !���������β��� [K]
  real(DP), save, public :: PressTop = 0.0d0      !���������Ǥΰ��� [Pa]
  real(DP), save, public :: CpDry  = 0.0d0        !������ʬ���갵��Ǯ [J/K kg]
  real(DP), save, public :: CpDryMol = 0.0d0      !������ʬ���갵��Ǯ [J/K kg]
  real(DP), save, public :: CvDry = 0.0d0         !������ʬ��������Ǯ [J/K kg]
  real(DP), save, public :: MolWtDry = 0.0d0      !������ʬ��ʬ����   [kg/mol]
  real(DP), save, public :: GasRDry  = 0.0d0      !������ʬ�ε������ [J/K kg]
  real(DP), save, public :: DayTime = 86400.0d0   ! 1 ����Ĺ�� [s]
  real(DP), save, public :: FactorJ = 1.0d0       !��ʪ�������Υѥ�᡼��
                                                  !�����Ǥ� 3.0d0
                                                  !�ϵ�Ǥ� 1.0d0 �Ȥ���
  ! ���֥롼����θ���
  !
  public constants_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine constants_Init
    !
    != ������롼����
    !
    ! namelist ������˴�Ť���ʪ�������. ����ϰʲ����̤�. 
    !
    ! * namelist �����갵��Ǯ (CpDry) ��ʿ��ʬ���� (MolWtDry) ��Ϳ����줿����, 
    !   �����򸵤˵������(GasRDry), ������Ǯ (CvDry), �갵�����Ǯ (CpDryMol) �����. 
    ! * namelist �����갵��Ǯ (CpDry) �ȵ������ (GasRDry) ��Ϳ����줿����, 
    !   �����򸵤�ʿ��ʬ����(MolWtDry), ������Ǯ (CvDry), �갵�����Ǯ (CpDryMol) �����. 
    ! * namelist ���鴥����ʬ�Υ���� (SpcDryMolFr) ��Ϳ����줿���ˤ�, 
    !   �����Ǯ�ϳإơ��֥���, �갵��Ǯ (CpDry), ������� (GasRDry), 
    !   ʿ��ʬ����(MolWtDry), ������Ǯ (CvDry), �갵�����Ǯ (CpDryMol) �����. 
    ! 
    
    !�⥸�塼���ɤ߹���
    use dc_types,      only: DP, STRING
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify          !��å���������
    use ChemData,      only: GasRUniv,             &!���׵������
      &                      ChemData_OneSpcID,    &!���ؼ�� ID
      &                      ChemData_CpPerMolRef, &!ɸ����֤Ǥ���Ǯ
      &                      ChemData_MolWt         !ʬ����
    use namelist_util, only: namelist_filename

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    integer                  :: SpcDryNum = 1   !������ʬ�β��ؼ�ο�
    character(20)            :: SpcDrySymbol(5) !������ʬ�β��ؼ�̾
    real(DP)                 :: SpcDryMolFr(5)  !������ʬ�β��ؼ��¸����
    integer, allocatable     :: SpcDryID(:)     !������ʬ�β��ؼ��ID
    real(DP), allocatable    :: PropertyDry(:)  !�������
    integer                  :: s               !����ѿ�
    integer                  :: unit            !�����ֹ�
    logical                  :: flag = .false.
    character(STRING)        :: Planet = ""
     
    !NAMELIST �����
    NAMELIST /constants_nml/ &
      & Planet, Grav, PressBasis, TempSfc, PressSfc, TempTop, PressTop, & 
      & SpcDrySymbol, SpcDryMolFr, DayTime, CpDry, MolWtDry, GasRDry
 
    SpcDrySymbol = '' 
    SpcDryMolFr  = 0.0d0
    
    !�ե����륪���ץ�. �������. 
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=constants_nml)
    close(unit)

    ! ���ϲ�®��
    !
    if (trim(Planet) == "Earth") then 
      FactorJ = 1.0d0
      Grav    = 9.8d0
    elseif (trim(Planet) == "Jupiter") then 
      FactorJ = 3.0d0
      Grav    = 23.1d0
    end if

    ! ��Ǯ, ʬ����, ������� 
    !
    if (CpDry /= 0.0d0 .AND. MolWtDry /= 0.0d0) then 
    ! namelist ���� CpDry ���ͤ����Ϥ��줿���

      !�갵�����Ǯ
      CpDryMol = CpDry *  MolWtDry 
      
      !�������
      GasRDry = GasRUniv / MolWtDry
      
      !������Ǯ
      CvDry    = CpDry - GasRDry

    elseif (CpDry /= 0.0d0 .AND. GasRDry /= 0.0d0) then 
    ! namelist ���� CpDry �� GasRDry ���ͤ����Ϥ��줿���

      !������ʬ��ʬ����
      MolWtDry = GasRUniv / GasRDry 

      !�갵�����Ǯ
      CpDryMol = CpDry *  MolWtDry 
      
      !������Ǯ
      CvDry    = CpDry - GasRDry

    elseif (SpcDryMolFr(1) /= 0.0d0) then 
    ! namelist ���� ����椬���Ϥ��줿���

      flag = .true.

      !----------------------------------------------------------
      ! ������ʬ��ʪ���ͤν����
      !
      !������ʬ�θĿ��������
      SpcDryNum = count(SpcDrySymbol /= "")
      
      !���ؼ�� ID �����    
      allocate(SpcDryID(SpcDryNum))    
      do s = 1, SpcDryNum
        SpcDryID(s) = ChemData_OneSpcID( SpcDrySymbol(s) )
      end do
      
      !�������ν���
      allocate(PropertyDry(SpcDryNum))
      
      !ʬ����
      do s = 1, SpcDryNum
        PropertyDry(s) = ChemData_MolWt(SpcDryID(s))
      end do
      MolWtDry = dot_product(PropertyDry, SpcDryMolFr(1:SpcDryNum)) 
      
      !�갵��Ǯ(�������)
      do s = 1, SpcDryNum    
        PropertyDry(s) = ChemData_CpPerMolRef(SpcDryID(s))
      end do
      CpDryMol = dot_product(PropertyDry, SpcDryMolFr(1:SpcDryNum)) 
      
      !�갵��Ǯ
      CpDry    = CpDryMol / MolWtDry
      
      !�������
      GasRDry = GasRUniv / MolWtDry
      
      !������Ǯ
      CvDry    = CpDry - GasRDry

    else
      call MessageNotify( "E", "constants_init", "Enough Variables are not set" )
    end if
    
    !----------------------------------------------------------
    ! ��ǧ
    !----------------------------------------------------------
    call MessageNotify( "M", &
      & "constants_init", "Grav = %f", d=(/Grav/) )
    call MessageNotify( "M", &
      &  "constants_init", "FactorJ = %f",  d=(/FactorJ/) )
    call MessageNotify( "M", &
      & "constants_init", "PressBasis = %f", d=(/PressBasis/))
    if (TempSfc /= 0) then 
      ! ��ɽ�̲��٤�Ϳ����줿��� 
      !
      call MessageNotify( "M", &
        & "constants_init", "TempSfc = %f",  d=(/TempSfc/) )
      call MessageNotify( "M", &
        & "constants_init", "PressSfc = %f", d=(/PressSfc/) )
    end if
    if (TempTop /= 0) then 
      ! ���������β��٤�Ϳ����줿��� 
      !
      call MessageNotify( "M", &
        & "constants_init", "TempTop = %f",  d=(/TempTop/) )
      call MessageNotify( "M", &
        & "constants_init", "PressTop = %f", d=(/PressTop/) )
    end if
    if (flag) then 
      ! ����椬Ϳ����줿���ˤϥץ�åȤ���. 
      !
      do s = 1, SpcDryNum
        call MessageNotify( "M", &
          &  "constants_init", "SpcDryID = %d",      i=(/SpcDryID(s)/))
        call MessageNotify( "M", &
          &  "constants_init", "SpcDrySymbol = %c", c1=trim(SpcDrySymbol(s)))
        call MessageNotify( "M", &
          &  "constants_init", "SpcDryMolFr = %f",   d=(/SpcDryMolFr(s)/))
      end do
    end if
    
    call MessageNotify( "M", "constants_init", "CpDry    = %f",    d=(/CpDry/) )
    call MessageNotify( "M", "constants_init", "CpDryMol = %f", d=(/CpDryMol/) )
    call MessageNotify( "M", "constants_init", "CvDry    = %f",    d=(/CvDry/) )
    call MessageNotify( "M", "constants_init", "GasRDry  = %f",  d=(/GasRDry/) )
    call MessageNotify( "M", "constants_init", "MolWtDry = %f", d=(/MolWtDry/) )
    call MessageNotify( "M", "constants_init", "DayTime  = %f",  d=(/DayTime/)  )    

  end subroutine Constants_Init

end module Constants
