!= ��ʬ�黻�⥸�塼�� (2 �����������ʬ)
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: differentiate_center2.f90,v 1.4 2007-04-11 11:59:58 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module differentiate_center2
  !
  ! 2 �����٤���ʬ�黻��Ԥ�����δؿ�����ޤȤ᤿�ѥå��������⥸�塼��
  ! ��ʿ Arakawa-C, ��ľ Lorentz ����åɤǤ���ʬ�׻�
  !

  !���ۤη�����ػ�
  implicit none

  !private °��
  private

  !��������
  public xz_dx_pz, xr_dx_pr
  public pz_dx_xz, pr_dx_xr
  public xz_dz_xr, pz_dz_pr
  public xr_dz_xz, pr_dz_pz

  !--------------------------------------------
  interface xz_dx_pz
    module procedure xa_dx_pa
  end interface xz_dx_pz

  interface xr_dx_pr
    module procedure xa_dx_pa
  end interface xr_dx_pr
  !--------------------------------------------
  interface pz_dx_xz
    module procedure pa_dx_xa
  end interface pz_dx_xz
  
  interface pr_dx_xr
    module procedure pa_dx_xa
  end interface pr_dx_xr
  !--------------------------------------------
  interface xz_dz_xr
    module procedure az_dz_ar
  end interface xz_dz_xr

  interface pz_dz_pr
    module procedure az_dz_ar
  end interface pz_dz_pr
  !--------------------------------------------
  interface xr_dz_xz
    module procedure ar_dz_az
  end interface xr_dz_xz

  interface pr_dz_pz
    module procedure ar_dz_az
  end interface pr_dz_pz
  !--------------------------------------------  

contains 

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function xa_dx_pa( pa_var ) 
    !
    ! x ��ʬ: Ⱦ�����ʻ��������ʬ
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !���󥵥���
    use axesset,   only : dx    !�ʻ����ֳ�
    
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: pa_var(imin:imax, kmin:kmax)
                                     !��ʬ�黻���оݤȤʤ��ѿ�
    real(DP)             :: xa_dx_pa(imin:imax, kmin:kmax)
                                     !��ʬ��
    integer             :: ix
    
    ! 1 ����ʬ�η׻�
    !
    do ix = imin+1, imax
      xa_dx_pa(ix,:) = ( pa_Var(ix,:) - pa_Var(ix-1,:) ) / dx
    end do
    
    xa_dx_pa(imin,:) = xa_dx_pa(imin+1,:)
    
  end function xa_dx_pa

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function pa_dx_xa( xa_var ) 
    !
    ! x ��ʬ: �����ʻ��������ʬ
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !���󥵥���
    use axesset,   only : dx    !�ʻ����ֳ�

    !���ۤη�����ػ�
    implicit none

    !�ѿ����    
    real(DP), intent(in) :: xa_var(imin:imax, kmin:kmax)
                                     !��ʬ�黻���оݤȤʤ��ѿ�
    real(DP)             :: pa_dx_xa(imin:imax, kmin:kmax)
                                     !��ʬ��    
    integer             :: ix

    ! 1 ����ʬ�η׻�
    !
    do ix = imin, imax-1
      pa_dx_xa(ix,:) = (xa_Var(ix+1,:) - xa_Var(ix,:))/dx
    end do
    
    pa_dx_xa(imax,:) = pa_dx_xa(imax-1,:)

  end function pa_dx_xa


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function az_dz_ar(ar_var) 
    !
    ! z ��ʬ: Ⱦ�����ʻ��������ʬ
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !���󥵥���
    use axesset,   only : dz    !�ʻ����ֳ�

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: ar_var(imin:imax, kmin:kmax)
                                     !��ʬ�黻���оݤȤʤ��ѿ�
    real(DP)             :: az_dz_ar(imin:imax, kmin:kmax)
                                     !��ʬ��
    integer             :: kz
    
    ! 1 ����ʬ�η׻�
    !
    do kz = kmin+1, kmax
      az_dz_ar(:,kz) = ( ar_Var(:,kz) - ar_Var(:,kz-1) ) / dz
    end do
    
    az_dz_ar(:,kmin) = az_dz_ar(:,kmin+1)
    
  end function az_dz_ar

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function ar_dz_az( az_var ) 
    !
    ! z ��ʬ: �����ʻ��������ʬ
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only : DP
    use gridset,   only : imin, imax, kmin, kmax
                                !���󥵥���
    use axesset,   only : dz    !�ʻ����ֳ�

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in) :: az_var(imin:imax, kmin:kmax)
                                     !��ʬ�黻���оݤȤʤ��ѿ�
    real(DP)             :: ar_dz_az(imin:imax, kmin:kmax)
                                     !��ʬ��
    integer             :: kz

    do kz = kmin, kmax-1
      ar_dz_az(:,kz) = ( az_Var(:,kz+1) - az_Var(:,kz) ) / dz
    end do
    
    ar_dz_az(:,kmax) = ar_dz_az(:,kmax-1)

  end function ar_dz_az

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
end module differentiate_center2
