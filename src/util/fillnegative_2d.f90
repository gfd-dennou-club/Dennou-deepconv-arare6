!= ��α����̤ʤɤ����᤹�뤿��Υ롼����
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: fillnegative.f90,v 1.7 2006-11-10 04:48:41 odakker Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module FillNegative
  !
  !��α����̤ʤɤ����᤹�뤿��Υ롼����
  !  * ��Ȥʤä����ˤ�, ���Ϥ� 8 �ʻ������ͤ��äƷ���᤹��.
  !  * F77�� deepconv �� QFILL.f ����� 
  ! 

  !�⥸�塼���ɤ߹���
  use dc_types,   only: DP

  !���ۤη�����ػ�
  implicit none

  !��������
  public FillNegative_Init
  public xza_FillNegative_xza

  !�ѿ����
  real(DP), private, allocatable, save  :: xza_Basic(:,:,:)
  real(DP), private, allocatable, save  :: xz_Dens(:,:)

contains

  subroutine FillNegative_Init( xza_VarBasic, xz_DensBasic )
    !
    ! ����ͤ����. 
    ! ���ܾ�Ⱦ�����ʬ���¤�������⾮�����ʤ뤳�Ȥϵ��Ƥ��ʤ�
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! ����� X �����β���
      &                   DimXMax,  &! ����� X �����ξ��
      &                   DimZMin,  &! ����� Z �����β���
      &                   DimZMax,  &! ����� Z �����ξ��
      &                   SpcNum     ! ���ؼ�ο�

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), intent(in) :: xza_VarBasic(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(DP), intent(in) :: xz_DensBasic(DimXMin:DimXMax, DimZMin:DimZMax)

    !�����
    allocate( xza_Basic(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum), &
      &       xz_Dens(DimXMin:DimXMax, DimZMin:DimZMax )          )
    
    !�ͤ�����
    xza_Basic = xza_VarBasic
    xz_Dens    = xz_DensBasic

  end subroutine FillNegative_Init
  

  function xza_FillNegative_xza( xza_Var )

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use gridset,    only: DimXMin,  &! ����� X �����β���
      &                   DimXMax,  &! ����� X �����ξ��
      &                   DimZMin,  &! ����� Z �����β���
      &                   DimZMax,  &! ����� Z �����ξ��
      &                   SpcNum     ! ���ؼ�ο�

    !���ۤη�����ػ�    
    implicit none

    !�ѿ����
    real(DP), intent(inout) :: xza_Var(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(DP)                :: xza_FillNegative_xza(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(DP)                :: xza_DQFILL(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(DP)                :: xza_QSUMPN(DimXMin:DimXMax, DimZMin:DimZMax, SpcNum)
    real(DP), parameter     :: EPS = 1.0d-60  !��Ǥγ�껻���ɤ�
    integer                 :: i, k, s
    
    !�����
    xza_QSUMPN = 0.0d0
    xza_DQFILL = 0.0d0


    do s = 1, SpcNum
      do k = DimZMin+2, DimZMax-2
        do i = DimXMin+2, DimXMax-2
          xza_QSUMPN(i,k,s) = 1.0d0 /                                     &
            & ( (    MAX( 0.0d0, xza_Basic(i-1,k,s) + xza_Var(i-1,k,s) )  &
            &         * xz_Dens(i-1,k)                                    &
            &      + MAX( 0.0d0, xza_Basic(i+1,k,s) + xza_Var(i+1,k,s) )  &
            &         * xz_Dens(i+1,k)                                    &
            &      + MAX( 0.0d0, xza_Basic(i,k-1,s) + xza_Var(i,k-1,s) )  &
            &         * xz_Dens(i,k-1)                                    &
            &      + MAX( 0.0d0, xza_Basic(i,k+1,s) + xza_Var(i,k+1,s) )  &
            &         * xz_Dens(i,k+1)                                    &
            &     ) * 0.75d0                                              &
            &  + (   MAX( 0.0d0, xza_Basic(i-2,k,s) + xza_Var(i-2,k,s) )  &
            &         * xz_Dens(i-2,k)                                    &
            &      + MAX( 0.0d0, xza_Basic(i+2,k,s) + xza_Var(i+2,k,s) )  &
            &         * xz_Dens(i+2,k)                                    &
            &      + MAX( 0.0d0, xza_Basic(i,k-2,s) + xza_Var(i,k-2,s) )  &
            &         * xz_Dens(i,k-2)                                    &
            &      + MAX( 0.0d0, xza_Basic(i,k+2,s) + xza_Var(i,k+2,s) )  &
            &         * xz_Dens(i,k+2)                                    &
            &     ) * 0.25d0                                              &
            &  + EPS ) 
        end do
      end do
    end do

    do s = 1, SpcNum
      do k = DimZMin+2, DimZMax-2
        do i = DimXMin+2, DimXMax-2
          xza_DQFILL(i,k,s) =                                               &
            &  - MIN( 0.0d0, xza_Basic(i,k,s) + xza_Var(i,k,s) )            &
            &  + MAX( 0.0d0, xza_Basic(i,k,s) + xza_Var(i,k,s) )            &
            &    * ( ( MIN( 0.0d0, xza_Basic(i-1,k,s) + xza_Var(i-1,k,s) )  &
            &           * xza_QSUMPN(i-1,k,s)                               &
            &           * xz_Dens(i-1,k)                                    &
            &        + MIN( 0.0d0, xza_Basic(i+1,k,s) + xza_Var(i+1,k,s) )  &
            &           * xza_QSUMPN(i+1,k,s)                               &
            &           * xz_Dens(i+1,k)                                    &
            &        + MIN( 0.0d0, xza_Basic(i,k-1,s) + xza_Var(i,k-1,s) )  &
            &           * xza_QSUMPN(i,k-1,s)                               &
            &           * xz_Dens(i,k-1)                                    &
            &        + MIN( 0.0d0, xza_Basic(i,k+1,s) + xza_Var(i,k+1,s) )  &
            &           * xza_QSUMPN(i,k+1,s)                               &
            &           * xz_Dens(i,k+1)                                    &  
            &       ) * 0.75d0                                              &
            &      + ( MIN( 0.0d0, xza_Basic(i-2,k,s) + xza_Var(i-2,k,s) )  &
            &           * xza_QSUMPN(i-2,k,s)                               &
            &           * xz_Dens(i-2,k)                                    &
            &        + MIN( 0.0d0, xza_Basic(i+2,k,s) + xza_Var(i+2,k,s) )  &
            &           * xza_QSUMPN(i+2,k,s)                               &
            &           * xz_Dens(i+2,k)                                    &
            &        + MIN( 0.0d0, xza_Basic(i,k-2,s) + xza_Var(i,k-2,s) )  &
            &           * xza_QSUMPN(i,k-2,s)                               &
            &           * xz_Dens(i,k-2)                                    &
            &        + MIN( 0.0d0, xza_Basic(i,k+2,s) + xza_Var(i,k+2,s) )  &
            &           * xza_QSUMPN(i,k+2,s)                               &
            &           * xz_Dens(i,k+2)                                    &
            &       ) * 0.25d0                                              &
            &     )
        end do
      end do

!      write(*,*) sum( xza_DQFILL(DimXMin:DimXMax,DimZMin:DimZMax,s) )
    end do

    !����
    xza_FillNegative_xza = xza_Var + xza_DQFILL
    
  end function xza_FillNegative_xza

end module FillNegative
