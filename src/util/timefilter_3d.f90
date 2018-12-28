!= Module TimeFilter_3D
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: timefilter_3d.f90,v 1.4 2007-08-20 07:33:02 odakker Exp $
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

module TimeFilter_3d

  !�⥸�塼���ɤ߹���
  use dc_types, only : DP

  use GridSet_3d, only: DimXMin,     & !x ����������β���
    &                 DimXMax,     & !x ����������ξ��
    &                 DimYMin,     & !y ����������β���
    &                 DimYMax,     & !y ����������ξ��
    &                 DimZMin,     & !z ����������β���
    &                 DimZMax,     & !z ����������ξ��
    &                 SpcNum         !���ؼ�ο�

  !���ۤη�����ػ�
  implicit none

  !°��������
  private

  !�ؿ��� public °��������
  public AsselinFilter

  !�ѿ����
  real(DP) :: tfil = 1.0d-1  !���å����λ��֥ե��륿�η���

  !�ͤ� save ����
  save tfil
  
contains
  
  subroutine AsselinFilter(VarA, VarN, VarB)
    !
    ! ���֥ե��륿��; Asselin �Υ�����ե��륿��������
    !

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in)     :: VarA(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
    real(DP), intent(inout)  :: VarN(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)  
    real(DP), intent(in)     :: VarB(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)
    real(DP)                 :: Var(DimXMin:DimXMax,DimYMin:DimYMax,DimZMin:DimZMax)


    !���֥ե��륿
    Var  = VarN + tfil * ( VarB  - 2.0d0 * VarN + VarA ) 
    VarN = Var

    
  end subroutine AsselinFilter
  
end module TimeFilter_3d
