!= Module TimeFilter
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: timefilter.f90,v 1.5 2007-04-19 14:29:33 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview
!
!���å����λ��֥ե��륿. leap-frog �θ��Ŭ��. 
! 
!== Error Handling
!
!== Bugs
!
!== Note
!
!== Future Plans
!
!

module TimeFilter

  !�⥸�塼���ɤ߹���
  use GridSet,  only: DimXMin,     & !x ����������β���
    &                 DimXMax,     & !x ����������ξ��
    &                 DimZMin,     & !z ����������β���
    &                 DimZMax,     & !z ����������ξ��
    &                 SpcNum         !���ؼ�ο�

  !���ۤη�����ػ�
  implicit none

  !°��������
  private

  !�ؿ��� public °��������
  public AsselinFilter_xz
  public AsselinFilter_pz
  public AsselinFilter_xr
  public AsselinFilter_xza

  interface AsselinFilter_xza
    module procedure AsselinFilter_aaa
  end interface 
  interface AsselinFilter_xz
    module procedure AsselinFilter_aa
  end interface 
  interface AsselinFilter_pz
    module procedure AsselinFilter_aa
  end interface
  interface AsselinFilter_xr
    module procedure AsselinFilter_aa
  end interface

  !�ѿ����
  real(8) :: tfil = 1.0d-1  !���å����λ��֥ե��륿�η���

  !�ͤ� save ����
  save tfil
  
contains
  
  subroutine AsselinFilter_aa(aa_VarA, aa_VarN, aa_VarB)
    !
    ! ���֥ե��륿��; Asselin �Υ�����ե��륿��������
    !

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)     :: aa_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8), intent(inout)  :: aa_VarN(DimXMin:DimXMax, DimZMin:DimZMax)
    real(8), intent(in)     :: aa_VarB(DimXMin:DimXMax, DimZMin:DimZMax)  
    real(8)                 :: aa_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !���֥ե��륿
    aa_Var  = aa_VarN + tfil * ( aa_VarB  - 2.0d0 * aa_VarN + aa_VarA ) 
    aa_VarN = aa_Var
    
  end subroutine AsselinFilter_aa
  

  subroutine AsselinFilter_aaa( aaa_VarA, aaa_VarN, aaa_VarB )
    !
    ! ���֥ե��륿��; Asselin �Υ�����ե��륿��������
    !

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)     :: aaa_VarA(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8), intent(inout)  :: aaa_VarN(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8), intent(in)     :: aaa_VarB(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(8)                 :: aaa_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)

    !���֥ե��륿
    aaa_Var  = aaa_VarN + tfil * ( aaa_VarB  - 2.0d0 * aaa_VarN + aaa_VarA ) 
    aaa_VarN = aaa_Var
    
  end subroutine AsselinFilter_aaa
  
end module TimeFilter
