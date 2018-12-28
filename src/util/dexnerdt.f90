!= Module DExnerDt
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: dexnerdt.f90,v 1.2 2014/03/04 05:55:07 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module DExnerDt
  !
  ! ������������, �����Ѳ����ŷ���ʬ�κ�����λ����Ѳ��ˤ����ɾ��
  ! ���뤿��Υ⥸�塼��

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP, STRING
  use dc_iounit,  only: FileOpen
  use dc_message, only: MessageNotify

  use gridSet,  only: imin,     & !x ����������β���
    &                 imax,     & !x ����������ξ��
    &                 jmin,     & !y ����������β���
    &                 jmax,     & !y ����������ξ��
    &                 kmin,     & !z ����������β���
    &                 kmax,     & !z ����������ξ��
    &                 ncmax
  use basicset, only:  &
    &                 xyz_PTempBZ,        & !���̤δ��ܾ�
    &                 xyz_VPTempBZ,       & !���̤δ��ܾ�
    &                 xyz_VelSoundBZ,     & !
    &                 xyzf_QMixBZ,        & !
    &                 xyz_QMixBZ,         & !
    &                 xyz_QMixBZPerMolWt, & !
    &                 xyz_DensBZ            !���ܾ��̩��
  use constants,only: MolWtDry, CpDry, GasRDry
  use composition,only : MolWtWet, GasNum
  
  !���ۤη�����ػ�
  implicit none

  !private °��
  private

  public xyz_DExnerDt_xyz_xyzf, xyz_DExnerDt_xyz, xyz_DExnerDt_xyzf
  public xy_DExnerDt_xy_xyf

contains

!!!------------------------------------------------------------------------!!!
  function xyz_DExnerDt_xyz_xyzf( xyz_DPTempDt, xyzf_DQMixDt )

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), intent(in) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in) :: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)             :: xyz_DExnerDt_xyz_xyzf(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_DQMixDtPerMolWt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    integer              :: s
    
    do s = 1, ncmax
      xyzf_DQMixDtPerMolWt(:,:,:,s) = xyzf_DQMixDt(:,:,:,s) / MolWtWet(s)
    end do
    xyz_DExnerDt_xyz_xyzf =                                     &
      &  (xyz_VelSoundBZ ** 2.0d0) / (CpDry * xyz_VPTempBZ)     &
      &    * (                                                  &
      &          xyz_DPTempDt / xyz_PTempBZ                     &
      &        + sum(xyzf_DQMixDtPerMolWt(:,:,:,1:GasNum),4)    &
      &           / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt )   &
      &        - sum(xyzf_DQMixDt,4) / ( 1.0d0 + xyz_QMixBZ )   &
      &      )

  end function xyz_DExnerDt_xyz_xyzf

!!!------------------------------------------------------------------------!!!
  function xyz_DExnerDt_xyzf( xyzf_DQMixDt )

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), intent(in) :: xyzf_DQMixDt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)             :: xyz_DExnerDt_xyzf(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyzf_DQMixDtPerMolWt(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    integer              :: s
    
    do s = 1, ncmax
      xyzf_DQMixDtPerMolWt(:,:,:,s) = xyzf_DQMixDt(:,:,:,s) / MolWtWet(s)
    end do
    xyz_DExnerDt_xyzf =                                         &
      &  (xyz_VelSoundBZ ** 2.0d0) / (CpDry * xyz_VPTempBZ)     &
      &    * (                                                  &
      &        + sum(xyzf_DQMixDtPerMolWt(:,:,:,1:GasNum),4)    &
      &           / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt )   &
      &        - sum(xyzf_DQMixDt,4) / ( 1.0d0 + xyz_QMixBZ )   &
      &      )

  end function xyz_DExnerDt_xyzf

!!!------------------------------------------------------------------------!!!
  function xyz_DExnerDt_xyz( xyz_DPTempDt )

    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP), intent(in) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: xyz_DExnerDt_xyz(imin:imax,jmin:jmax,kmin:kmax)
    
    xyz_DExnerDt_xyz =                                          &
      &  (xyz_VelSoundBZ ** 2.0d0) / (CpDry * xyz_VPTempBZ)     &
      &    * (                                                  &
      &          xyz_DPTempDt / xyz_PTempBZ                     &
      &      )
    
  end function xyz_DExnerDt_xyz

!!!------------------------------------------------------------------------!!!
  function xy_DExnerDt_xy_xyf( xy_DPTempDt, xyf_DQMixDt, kz )

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), intent(in) :: xy_DPTempDt(imin:imax,jmin:jmax)
    real(DP), intent(in) :: xyf_DQMixDt(imin:imax,jmin:jmax,ncmax)
    integer,  intent(in) :: kz
    real(DP)             :: xy_DExnerDt_xy_xyf(imin:imax,jmin:jmax)
    real(DP)             :: xyf_DQMixDtPerMolWt(imin:imax,jmin:jmax,ncmax)
    integer              :: s

    do s = 1, ncmax
      xyf_DQMixDtPerMolWt(:,:,s) = xyf_DQMixDt(:,:,s) / MolWtWet(s)
    end do
    
    xy_DExnerDt_xy_xyf =                                                  &
      &  (xyz_VelSoundBZ(:,:,kz) ** 2.0d0) / (CpDry * xyz_VPTempBZ(:,:,kz))&
      &   * (                                                              &
      &       xy_DPTempDt / xyz_PTempBZ(:,:,kz)                            &
      &     + sum(xyf_DQMixDtPerMolWt(:,:,1:GasNum),3)                     &
      &        / ( 1.0d0 / MolWtDry + xyz_QMixBZPerMolWt(:,:,kz) )         &
      &     - sum(xyf_DQMixDt,3) / ( 1.0d0 + xyz_QMixBZ(:,:,kz) )          &
      &    )
    
  end function xy_DExnerDt_xy_xyf

end module DExnerDt

