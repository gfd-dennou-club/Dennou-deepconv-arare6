!= �ǥե���Ȥδ��ܾ�����ꤹ�뤿����ѿ����ȷ��⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: basicset.f90,v 1.18 2014/07/08 01:05:32 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module basicset
  !
  != �ǥե���Ȥδ��ܾ�����ꤹ�뤿����ѿ����ȷ��⥸�塼��
  !
  ! hogeBZ �Ȥʤ������, ���������ͤ����Ϥ��뤳�Ȥǽ������Ԥ����Ȥ����. 
  !

  !�⥸�塼���ɤ߹���
  use dc_types,   only:  DP

  !���ۤη�����ػ�
  implicit none

  !save °��
  private

  !Public Interface
  real(DP), allocatable, save, public :: xz_DensBZ(:,:)     !̩��
  real(DP), allocatable, save, public :: pz_DensBZ(:,:)     !̩��
  real(DP), allocatable, save, public :: xqz_DensBZ(:,:)    !̩��
  real(DP), allocatable, save, public :: xr_DensBZ(:,:)     !̩��
  real(DP), allocatable, save, public :: xz_PressBZ(:,:)    !̵��������
  real(DP), allocatable, save, public :: xz_ExnerBZ(:,:)    !̵��������
  real(DP), allocatable, save, public :: xz_TempBZ(:,:)     !����
  real(DP), allocatable, save, public :: xz_PTempBZ(:,:)    !����
  real(DP), allocatable, save, public :: xr_PTempBZ(:,:)    !����
  real(DP), allocatable, save, public :: xz_VPTempBZ(:,:)   !������
  real(DP), allocatable, save, public :: pz_VPTempBZ(:,:)   !������
  real(DP), allocatable, save, public :: xqz_VPTempBZ(:,:)  !������
  real(DP), allocatable, save, public :: xr_VPTempBZ(:,:)   !������
  real(DP), allocatable, save, public :: xz_VelSoundBZ(:,:) !��®
  real(DP), allocatable, save, public :: xz_VelSW(:,:)      !��®
  real(DP), allocatable, save, public :: xzf_QMixBZ(:,:,:)  !�Ž���ʬ������
  real(DP), allocatable, save, public :: xz_EffMolWtBZ(:,:) !ʬ���̸���
  real(DP), allocatable, save, public :: xz_QMixBZPerMolWt(:,:)
                              !���ܾ�κ����� / ʬ���� ������
  real(DP), allocatable, save, public :: xr_QMixBZPerMolWt(:,:)
                              !���ܾ�κ����� / ʬ���� ������
  real(DP), allocatable, save, public :: xz_QMixBZ(:,:)
                              !���ܾ�κ����������
  real(DP), allocatable, save, public :: xr_QMixBZ(:,:)
                              !���ܾ�κ����������

  public basicset_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine basicset_init(                                &
    &   xz_Press, xz_Exner, xz_Temp, xz_PTemp, xz_Dens,    & 
    &   xz_VelSound, xzf_QMix, xz_EffMolWt                 &
    & )
    !
    != ���ܾ���ͤ��������������. 
    !
    ! dry �ξ��Ϻ�����䴥����ʬ�ȶŷ���ʬ��¸������ۤ˻Ȥ�ʤ��Τ�,
    ! �������ѿ��� optional �ˤ��Ƥ���. 
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only:  DP
    use gridset,    only:  imin,       &! ����� X �����β���
      &                    imax,       &! ����� X �����ξ��
      &                    jmin,       &! ����� Z �����β���
      &                    jmax,       &! ����� Z �����ξ��
      &                    kmin,       &! ����� Z �����β���
      &                    kmax,       &! ����� Z �����ξ��
      &                    ncmax        ! ���ؼ�ο�
    use composition, only: GasNum,     &!���Το�
      &                    IdxG,       &!���Τ�����ź����
      &                    MolWtWet     !�Ž���ʬ��ʬ����
    
    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP), intent(in) :: xz_Press(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Exner(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Temp(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_PTemp(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_Dens(imin:imax,kmin:kmax)
    real(DP), intent(in) :: xz_VelSound(imin:imax,kmin:kmax)
    real(DP), intent(in), optional :: xzf_QMix(imin:imax,kmin:kmax,ncmax)
    real(DP), intent(in), optional :: xz_EffMolWt(imin:imax,kmin:kmax)

    real(DP)             :: xzf_QMixBZPerMolWt(imin:imax,kmin:kmax,GasNum)
    integer              :: s

    ! ����ν����
    call basicset_array_init
    
    !�ͤ�����
    xz_PressBZ    = xz_Press
    xz_ExnerBZ    = xz_Exner
    xz_TempBZ     = xz_Temp
    xz_PTempBZ    = xz_PTemp
    xz_DensBZ     = xz_Dens
    xz_VelSoundBZ = xz_VelSound
    xz_VelSW      = xz_VelSound  ! ����̾��û�����뤿���. 

    if (present(xz_EffMolWt)) xz_EffMolWtBZ = xz_EffMolWt

    if (present(xzf_QMix)) then 
      xzf_QMixBZ = xzf_QMix
      xzf_QMixBZPerMolWt = 0.0d0

      do s = 1, GasNum
        xzf_QMixBZPerMolWt(:,:,s) = &
          & xzf_QMixBZ(:,:,IdxG(s)) / MolWtWet(IdxG(s))
      end do

      xz_QMixBZPerMolWt = sum(xzf_QMixBZPerMolWt, 3) 
      xz_QMixBZ         = sum(xzf_QMixBZ, 3) 
      
      xr_QMixBZPerMolWt = xr_xz( xz_QMixBZPerMolWt )
      xr_QMixBZ         = xr_xz( xz_QMixBZ )

    end if

    xz_VPTempBZ = xz_PTempBZ / xz_EffMolWtBZ

    xr_PTempBZ  = xr_xz( xz_PTempBZ )
    xr_DensBZ   = xr_xz( xz_DensBZ )
    xr_VPTempBZ = xr_xz( xz_VPTempBZ )

    ! ��ʿ���ͤʤΤ�, ʿ������ɬ�פʤ�
    !
    pz_DensBZ    = xz_Dens
    xqz_DensBZ   = xz_Dens
    pz_VPTempBZ  = xz_VPTempBZ 
    xqz_VPTempBZ = xz_VPTempBZ 

  contains

    subroutine basicset_array_init
      !
      ! *BasicZ ������ν����
      !

      !���ۤη�����ػ�
      implicit none

      allocate( & 
        & xz_DensBZ(imin:imax,kmin:kmax),        &
        & pz_DensBZ(imin:imax,kmin:kmax),        &
        & xqz_DensBZ(imin:imax,kmin:kmax),       &
        & xr_DensBZ(imin:imax,kmin:kmax),        &
        & xz_PressBZ(imin:imax,kmin:kmax),       &
        & xz_ExnerBZ(imin:imax,kmin:kmax),       &
        & xz_TempBZ(imin:imax,kmin:kmax),        &
        & xz_PTempBZ (imin:imax,kmin:kmax),      &
        & xr_PTempBZ (imin:imax,kmin:kmax),      &
        & xz_VPTempBZ (imin:imax,kmin:kmax),     &
        & pz_VPTempBZ (imin:imax,kmin:kmax),     &
        & xqz_VPTempBZ (imin:imax,kmin:kmax),    &
        & xr_VPTempBZ (imin:imax,kmin:kmax),     &
        & xz_VelSoundBZ(imin:imax,kmin:kmax),    &
        & xz_VelSW(imin:imax,kmin:kmax),         &
        & xzf_QMixBZ(imin:imax,kmin:kmax,ncmax), &
        & xz_EffMolWtBZ(imin:imax,kmin:kmax)     &
        & )
      allocate( &
        & xr_QMixBZPerMolWt(imin:imax,kmin:kmax),&
        & xz_QMixBZPerMolWt(imin:imax,kmin:kmax),&
        & xr_QMixBZ(imin:imax,kmin:kmax),        &
        & xz_QMixBZ(imin:imax,kmin:kmax)         &
        & )
      
      ! �ͤγ���
      xz_DensBZ     = 0.0d0
      pz_DensBZ     = 0.0d0
      xqz_DensBZ    = 0.0d0
      xr_DensBZ     = 0.0d0
      xz_PressBZ    = 0.0d0
      xz_ExnerBZ    = 0.0d0
      xz_TempBZ     = 0.0d0
      xz_PTempBZ    = 0.0d0
      xr_PTempBZ    = 0.0d0
      xz_VPTempBZ   = 0.0d0
      pz_VPTempBZ   = 0.0d0
      xqz_VPTempBZ  = 0.0d0
      xr_VPTempBZ   = 0.0d0
      xz_VelSoundBZ = 0.0d0
      xz_VelSW      = 0.0d0
      xzf_QMixBZ    = 0.0d0
      xz_EffMolWtBZ = 1.0d0 ! dry �ξ���ɬ�� 1.0
      
      xr_QMixBZPerMolWt = 0.0d0
      xz_QMixBZPerMolWt = 0.0d0
      xr_QMixBZ = 0.0d0
      xz_QMixBZ = 0.0d0
      
    end subroutine basicset_array_init
    

    function xr_xz(xz_Var)
      !
      ! ʿ������Ԥ� z ����Ⱦ�����ʻ����������ͤ������ʻ�������֤�
      
      implicit none
      
      real(DP),intent(in) :: xz_Var(imin:imax,kmin:kmax) 
      real(DP)            :: xr_xz(imin:imax,kmin:kmax)
      integer             :: kz
      
      do kz = kmin, kmax-1
        xr_xz(:,kz) = ( xz_Var(:,kz+1) + xz_Var(:,kz) ) * 5.0d-1
      end do
      
      xr_xz(:,kmax) = xr_xz(:,kmax-1)
      
    end function xr_xz
    
    
    function xz_xr(xr_Var)
      !
      ! ʿ������Ԥ� z ����Ⱦ�����ʻ����������ͤ������ʻ�������֤�
      
      implicit none
      
      real(DP),intent(in) :: xr_Var(imin:imax,kmin:kmax) 
      real(DP)            :: xz_xr(imin:imax,kmin:kmax)
      integer             :: kz
      
      do kz = kmin+1, kmax
        xz_xr(:,kz) = ( xr_Var(:,kz) + xr_Var(:,kz-1) ) * 5.0d-1
      end do
      
      xz_xr(:,kmin) = xz_xr(:,kmin+1)
      
    end function xz_xr

  end subroutine basicset_init

end module basicset
