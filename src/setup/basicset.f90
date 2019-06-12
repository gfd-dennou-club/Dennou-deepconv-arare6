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
  use dc_types, only:  DP

  !���ۤη�����ػ�
  implicit none

  !save °��
  private

  !Public Interface
  real(DP), allocatable, save, public :: xyz_DensBZ(:,:,:)    !̩��
  real(DP), allocatable, save, public :: pyz_DensBZ(:,:,:)    !̩��
  real(DP), allocatable, save, public :: xqz_DensBZ(:,:,:)    !̩��
  real(DP), allocatable, save, public :: xyr_DensBZ(:,:,:)    !̩��
  real(DP), allocatable, save, public :: xyz_PressBZ(:,:,:)   !̵��������
  real(DP), allocatable, save, public :: xyz_ExnerBZ(:,:,:)   !̵��������
  real(DP), allocatable, save, public :: xyz_TempBZ(:,:,:)    !����
  real(DP), allocatable, save, public :: xyz_PTempBZ(:,:,:)   !����
  real(DP), allocatable, save, public :: xyr_PTempBZ(:,:,:)   !����
  real(DP), allocatable, save, public :: xyz_VPTempBZ(:,:,:)  !������
  real(DP), allocatable, save, public :: pyz_VPTempBZ(:,:,:)  !������
  real(DP), allocatable, save, public :: xqz_VPTempBZ(:,:,:)  !������
  real(DP), allocatable, save, public :: xyr_VPTempBZ(:,:,:)  !������
  real(DP), allocatable, save, public :: xyz_VelSoundBZ(:,:,:)!��®
  real(DP), allocatable, save, public :: xyz_VelSW(:,:,:)     !��®
  real(DP), allocatable, save, public :: xyzf_QMixBZ(:,:,:,:) !�Ž���ʬ������
  real(DP), allocatable, save, public :: xyz_EffMolWtBZ(:,:,:)!ʬ���̸���
  real(DP), allocatable, save, public :: xyz_QMixBZPerMolWt(:,:,:)
                              !���ܾ�κ����� / ʬ���� ������
  real(DP), allocatable, save, public :: xyr_QMixBZPerMolWt(:,:,:)
                              !���ܾ�κ����� / ʬ���� ������
  real(DP), allocatable, save, public :: xyz_QMixBZ(:,:,:)
                              !���ܾ�κ����������
  real(DP), allocatable, save, public :: xyr_QMixBZ(:,:,:)
                              !���ܾ�κ����������

  public basicset_init

contains

!!!-----------------------------------------------------------------!!!
  subroutine basicset_init(                                     &
    &   xyz_Press, xyz_Exner, xyz_Temp, xyz_PTemp, xyz_Dens,    & 
    &   xyz_VelSound, xyzf_QMix, xyz_EffMolWt                   &
    & )
    !
    != ���ܾ���ͤ��������������. 
    !
    ! dry �ξ��Ϻ�����䴥����ʬ�ȶŷ���ʬ��¸������ۤ˻Ȥ�ʤ��Τ�,
    ! �������ѿ��� optional �ˤ��Ƥ���. 
    !

    !�⥸�塼��ƤӽФ�
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
    real(DP), intent(in) :: xyz_Press(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_Dens(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyz_VelSound(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in), optional :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(in), optional :: xyz_EffMolWt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)             :: xyzf_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax,GasNum)
    integer              :: s, k

    ! ����ν����
    call basicset_array_init

    !�ͤ�����
    xyz_PressBZ    = xyz_Press
    xyz_ExnerBZ    = xyz_Exner
    xyz_TempBZ     = xyz_Temp
    xyz_PTempBZ    = xyz_PTemp
    xyz_DensBZ     = xyz_Dens
    xyz_VelSoundBZ = xyz_VelSound
    xyz_VelSW      = xyz_VelSound  ! ����̾��û�����뤿���. 

    if (present(xyz_EffMolWt)) xyz_EffMolWtBZ = xyz_EffMolWt

    if (present(xyzf_QMix)) then 
      xyzf_QMixBZ = xyzf_QMix
      xyzf_QMixBZPerMolWt = 0.0d0
      
      do s = 1, GasNum
        xyzf_QMixBZPerMolWt(:,:,:,s) = &
          & xyzf_QMixBZ(:,:,:,IdxG(s)) / MolWtWet(IdxG(s))
      end do
      
      xyz_QMixBZPerMolWt = sum(xyzf_QMixBZPerMolWt, 4) 
      xyz_QMixBZ         = sum(xyzf_QMixBZ, 4) 

      do k = kmin, kmax-1
        xyr_QMixBZPerMolWt(:,:,k) &
          & = ( xyz_QMixBZPerMolWt(:,:,k+1) + xyz_QMixBZPerMolWt(:,:,k) ) * 0.5d0
        xyr_QMixBZ(:,:,k) = ( xyz_QMixBZ(:,:,k+1) + xyz_QMixBZ(:,:,k) )   * 0.5d0
      end do
    end if
    
    xyz_VPTempBZ = xyz_PTempBZ / xyz_EffMolWtBZ

    do k = kmin, kmax-1    
      xyr_PTempBZ(:,:,k)  = (  xyz_PTempBZ(:,:,k+1) +  xyz_PTempBZ(:,:,k) ) * 0.5d0
      xyr_VPTempBZ(:,:,k) = ( xyz_VPTempBZ(:,:,k+1) + xyz_VPTempBZ(:,:,k) ) * 0.5d0
      xyr_DensBZ(:,:,k)   = (   xyz_DensBZ(:,:,k+1) +   xyz_DensBZ(:,:,k) ) * 0.5d0
    end do

    ! ��ʿ���ͤʤΤ�, ʿ������ɬ�פʤ�
    !
    pyz_DensBZ   = xyz_Dens
    xqz_DensBZ   = xyz_Dens
    pyz_VPTempBZ = xyz_VPTempBZ 
    xqz_VPTempBZ = xyz_VPTempBZ 
  
  contains

    subroutine basicset_array_init
      !
      ! *BZ ������ν����
      !
      
      allocate( & 
        & xyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & pyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xqz_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyr_DensBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_PressBZ(imin:imax,jmin:jmax,kmin:kmax),       &
        & xyz_ExnerBZ(imin:imax,jmin:jmax,kmin:kmax),       &
        & xyz_TempBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_PTempBZ (imin:imax,jmin:jmax,kmin:kmax),      &
        & xyr_PTempBZ (imin:imax,jmin:jmax,kmin:kmax),      &
        & xyz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & pyz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xqz_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xyr_VPTempBZ (imin:imax,jmin:jmax,kmin:kmax),     &
        & xyz_VelSoundBZ(imin:imax,jmin:jmax,kmin:kmax),    &
        & xyz_VelSW(imin:imax,jmin:jmax,kmin:kmax),         &
        & xyzf_QMixBZ(imin:imax,jmin:jmax,kmin:kmax,ncmax), &
        & xyz_EffMolWtBZ(imin:imax,jmin:jmax,kmin:kmax)     &
        & )
      allocate( &
        & xyr_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax),&
        & xyz_QMixBZPerMolWt(imin:imax,jmin:jmax,kmin:kmax),&
        & xyr_QMixBZ(imin:imax,jmin:jmax,kmin:kmax),        &
        & xyz_QMixBZ(imin:imax,jmin:jmax,kmin:kmax)         &
        & )
      
      ! �ͤγ���
      xyz_DensBZ     = 0.0d0
      pyz_DensBZ     = 0.0d0
      xqz_DensBZ     = 0.0d0
      xyr_DensBZ     = 0.0d0
      xyz_PressBZ    = 0.0d0   
      xyz_ExnerBZ    = 0.0d0   
      xyz_TempBZ     = 0.0d0   
      xyz_PTempBZ    = 0.0d0
      xyz_VPTempBZ   = 0.0d0
      pyz_VPTempBZ   = 0.0d0
      xqz_VPTempBZ   = 0.0d0
      xyr_VPTempBZ   = 0.0d0
      xyz_VelSoundBZ = 0.0d0
      xyz_VelSW      = 0.0d0
      xyzf_QMixBZ    = 0.0d0
      xyz_EffMolWtBZ = 1.0d0 ! dry �ξ���ɬ�� 1.0
      
      xyr_QMixBZPerMolWt = 0.0d0
      xyz_QMixBZPerMolWt = 0.0d0
      xyr_QMixBZ = 0.0d0
      xyz_QMixBZ = 0.0d0
      
    end subroutine basicset_array_init

  end subroutine basicset_init
    
end module basicset
