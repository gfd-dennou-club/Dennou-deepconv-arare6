!= ʪ���������������
!
!= Physical and mathematical constants settings
!
! Authors::   Yasuhiro MORIKAWA, Yoshiyuki O. Takahashi
! Version::   $Id: constants0.f90,v 1.3 2012/07/18 07:30:40 odakker Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module constants0
  !
  != ʪ���������������
  !
  != Physical and mathematical constants settings
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! ʪ����������������ꤪ����ݴɤ�Ԥ��ޤ�. 
  ! �ǥե�����ͤ��ϵ��絤�����ꤷ���ͤ����ꤵ��Ƥ��ޤ�. 
  !
  ! Physical and mathematical constants are set and stored. 
  ! By default, values on atmosphere of earth are set. 
  !
  !== Procedures List
  !
  ! Constants0Init :: ʪ�����������
  ! ------------  :: ------------
  ! Constants0Init :: Settings of physical constants
  !
  !== NAMELIST
  !
  ! N/A
  !


  ! �⥸�塼����� ; USE statements
  !

  ! ���̷��ѥ�᥿
  ! Kind type parameter
  !
  use dc_types, only: DP     ! �����ټ¿���. Double precision. 

  ! ���ʸ ; Declaration statements
  !
  implicit none
  private

  ! ������³��
  ! Public procedure
  !
  public:: Constants0Init

  ! �����ѿ�
  ! Public variables
  !
  logical, save, public:: constants0_inited = .false.
                              ! �������ե饰. 
                              ! Initialization flag

  real(DP), parameter, public:: PI = 3.1415926535897932_DP
                              ! $ \pi $ .
                              ! �߼�Ψ.  Circular constant
  real(DP), parameter, public:: GasRUniv = 8.314_DP
                              ! $ R^{*} $ [J K-1 mol-1].
                              ! ���׵������.  Universal gas constant
  real(DP), parameter, public:: StB = 5.67e-8_DP
                              ! $ \sigma_{SB} $ . 
                              ! ���ƥե���ܥ�ĥޥ����. 
                              ! Stefan-Boltzmann constant
  real(DP), parameter, public :: FKarm = 0.4d0 
                              ! ����ޥ����
                              ! Karmann constant

  ! ������ѿ�
  ! Private variables
  !

  character(*), parameter:: module_name = 'constants0'
                              ! �⥸�塼���̾��. 
                              ! Module name
  character(*), parameter:: version = &
    & '$Name:  $' // &
    & '$Id: constants0.f90,v 1.3 2012/07/18 07:30:40 odakker Exp $'
                              ! �⥸�塼��ΥС������
                              ! Module version

contains

  subroutine Constants0Init
    !
    ! constants0 �⥸�塼��ν������Ԥ��ޤ�. 
    !
    ! "constants0" module is initialized. 
    !

    ! �⥸�塼����� ; USE statements
    !

    ! ��å���������
    ! Message output
    !
    use dc_message, only: MessageNotify

    ! ���ʸ ; Declaration statements
    !
    implicit none

    ! �¹�ʸ ; Executable statement
    !

    if ( constants0_inited ) return


    ! ���� ; Print
    !
!    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
!    call MessageNotify( 'M', module_name, '  PI               = %f', d = (/ PI               /) )
!    call MessageNotify( 'M', module_name, '  GasRUniv         = %f', d = (/ GasRUniv         /) )
!    call MessageNotify( 'M', module_name, '  StB              = %f', d = (/ StB              /) )
!    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    constants0_inited = .true.

  end subroutine Constants0Init

end module constants0
