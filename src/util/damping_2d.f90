!= ���ݥ��إ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: damping.f90,v 1.8 2010-08-13 07:18:17 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module Damping
  !
  ! ���ݥ��� (�����ն���Ȥ�ȿ�ͤ��ޤ��ۼ����뤿�����) �ˤ�����
  ! ����Ψ�Ȥ��η׻���Ԥ�����Υѥå��������⥸�塼��
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only : DP

  !���ۤη�����ػ�
  implicit none

  !private °�������
  private 
  
  !�ؿ��ˤ� public °�������
  public Damping_Init
  public DampSound_Init
  public DampSponge_Init
  public DampSponge_xz
  public DampSponge_xr
  public DampSponge_pz
  public xz_Sponge
  public xr_Sponge
  public pz_Sponge

  !�ѿ����
  real(DP), save, public   :: DampSound  = 0.0d0   !���ȸ����θ��그��
  real(DP), save, private  :: EFTime     = 100.0d0 !���ݥ��ؤ� e-folding time
  real(DP), save, private  :: DampDepthH = 0.0d0   !���ݥ��ؤθ���(��ʿ����)
  real(DP), save, private  :: DampDepthV = 0.0d0   !���ݥ��ؤθ���(��ľ����)
  real(DP), allocatable, save, private :: xz_DampRateH(:,:) !ss �ʻ������그��(��ʿ����)
  real(DP), allocatable, save, private :: xz_DampRateV(:,:) !ss �ʻ������그��(��ľ����)
  real(DP), allocatable, save, private :: pz_DampRateH(:,:) !fs �ʻ������그��(��ʿ����)
  real(DP), allocatable, save, private :: pz_DampRateV(:,:) !fs �ʻ������그��(��ľ����)
  real(DP), allocatable, save, private :: xr_DampRateH(:,:) !sf �ʻ������그��(��ʿ����)
  real(DP), allocatable, save, private :: xr_DampRateV(:,:) !sf �ʻ������그��(��ľ����)
  
contains 
  
!!!------------------------------------------------------------------------!!!
  subroutine Damping_Init( cfgfile ) 
    !
    ! ���ȸ����ȥ��ݥ��ؤθ��그���ν����
    ! 

    !�⥸�塼��ƤӽФ�
    use dc_types,  only : DP
 
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    character(*), intent(in) :: cfgfile
    real(DP)                 :: Alpha    !���ȸ����η���
    real(DP)                 :: Time     !
    real(DP)                 :: DepthH   !���ݥ��ؤθ���(��ʿ����)
    real(DP)                 :: DepthV   !���ݥ��ؤθ���(��ľ����)

    !NAMELIST �������
    NAMELIST /damping/ Alpha, Time, DepthH, DepthV
    open (10, FILE=cfgfile)
    read(10, NML=damping)
    close(10)

    !�����
    call DampSound_Init( Alpha ) 
    call DampSponge_Init( Time, DepthH, DepthV )

  end subroutine Damping_Init


!!!------------------------------------------------------------------------!!!
  subroutine DampSound_Init( damp ) 
    !
    ! ���ȸ����θ��그���ν����
    ! 

    !�⥸�塼��ƤӽФ�
    use dc_types,    only : DP
    use dc_message,  only : MessageNotify
    use mpi_wrapper, only : myrank
    use axesset,     only : DelX,        &! x �����γʻ����ֳ�
      &                     DelZ          ! z �����γʻ����ֳ�
    use timeset,     only : DelTimeShort  !û�����֥��ƥå�

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: damp

    !-------------------------------------------------------------------
    ! ���ȸ����θ���Ψ   Min(DelX, DelZ) ** 2.0 ������
    !-------------------------------------------------------------------
    DampSound = Damp * ( Min(DelX, DelZ) ** 2.0d0 ) / DelTimeShort
    
    !-----------------------------------------------------------------    
    ! �ͤγ�ǧ
    !-----------------------------------------------------------------
    if (myrank == 0) then 
      call MessageNotify( "M", &
        & "DampSound_init", "DampSound = %f", d=(/DampSound/) )
    end if
  end subroutine DampSound_Init


!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_Init( Time, DepthH, DepthV )
    !
    ! ���ݥ��ؤθ��그�������
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use dc_message,  only: MessageNotify
    use mpi_wrapper, only: myrank
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��
    use axesset,     only: DelX,          &! x �����γʻ����ֳ�
      &                    DelZ,          &! z �����γʻ����ֳ�
      &                    s_X,           &!X ��ɸ��(�����顼�ʻ���)
      &                    s_Z,           &!Z ��ɸ��(�����顼�ʻ���)
      &                    f_X,           &!X ��ɸ��(�ե�å����ʻ���)
      &                    f_Z,           &!Z ��ɸ��(�ե�å����ʻ���)
      &                    XMax,          &!X ��ɸ�κ�����
      &                    ZMax            !Z ��ɸ�κ�����

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), intent(in)   :: Time     !���ݥ��ؤ� e-folding time
    real(DP), intent(in)   :: DepthH   !���ݥ��ؤθ���(��ʿ����)
    real(DP), intent(in)   :: DepthV   !���ݥ��ؤθ���(��ľ����)
    real(DP), parameter    :: Pi =3.1415926535897932385d0   !�߼�Ψ
    integer               :: i, k

    !�����
    allocate( &
      & xz_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xz_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax), &
      & pz_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & pz_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xr_DampRateH(DimXMin:DimXMax, DimZMin:DimZMax), &
      & xr_DampRateV(DimXMin:DimXMax, DimZMin:DimZMax)    )
    xz_DampRateH = 0.0d0
    xz_DampRateV = 0.0d0
    pz_DampRateH = 0.0d0
    pz_DampRateV = 0.0d0
    xr_DampRateH = 0.0d0
    xr_DampRateV = 0.0d0

    !�ͤ�����
    EFTime     = Time
    DampDepthH = DepthH
    DampDepthV = DepthV
    
    !-----------------------------------------------------------------    
    ! ���ݥ��ؤθ���Ψ
    !-----------------------------------------------------------------
    !��ʿ��������¦����¦����
    if ( DampDepthH < DelX ) then 
      call MessageNotify( "W", &
        & "DampSponge_Init", "DampDepthH is too thin. DelX is %f", d=(/DelX/))

    else
      do i = DimXMin, DimXMax
        !�����顼�ʻ�������¦����
        if ( s_X(i) < DampDepthH) then 
          xz_DampRateH(i,:) = ((1.0d0 - s_X(i) / DampDepthH) ** 3.0d0) / EFTime
        end if
        
        !�ե�å����ʻ�������¦����
        if ( f_X(i) < DampDepthH) then 
          pz_DampRateH(i,:) = ((1.0d0 - f_X(i) / DampDepthH) ** 3.0d0) / EFTime
         end if
        
        !�����顼�ʻ�������¦����    
        if ( s_X(i) > ( XMax - DampDepthH ) ) then 
          xz_DampRateH(i,:) = &
            & ((1.0d0 - (XMax - s_X(i)) / DampDepthH) ** 3.0d0) / EFTime 
        end if
        
        !�ե�å����ʻ�������¦����    
        if ( f_X(i) > ( XMax - DampDepthH ) ) then 
          pz_DampRateH(i,:) = &
            & ((1.0d0 - (XMax - f_X(i)) / DampDepthH) ** 3.0d0) / EFTime 
        end if
      end do
    end if
    !sf �� ss �� X �����˴ؤ��Ƥ�Ʊ��
    xr_DampRateH  = xz_DampRateH
    
    !��ľ�����ξ�������    
    if ( DampDepthV < DelZ ) then 
      call MessageNotify( "W", &
        & "DampSponge_Init", "DampDepthV is too thin. DelZ is %f", d=(/DelZ/) )
      
    else
      do k = DimZMin, DimZMax
        !�����顼�ʻ���
        if ( s_Z(k) >= ( ZMax - DampDepthV ) ) then 
          xz_DampRateV(:,k) =  &
            & (1.0d0 - dcos(Pi * (s_Z(k) - ZMax + DampDepthV) / DampDepthV)) &
            &  / EFTime 
        end if
        
        !�ե�å����ʻ���
        if ( f_Z(k) >= ( ZMax - DampDepthV ) ) then 
          xr_DampRateV(:,k) =  &
            & (1.0d0 - dcos(Pi * (f_Z(k) - ZMax + DampDepthV)/ DampDepthV)) &
            &  / EFTime 
        end if
      end do
    end if
    !fs �� ss �� Z �����˴ؤ��Ƥ�Ʊ��
    pz_DampRateV  = xz_DampRateV
    
    !-----------------------------------------------------------------    
    ! �ͤγ�ǧ
    !-----------------------------------------------------------------
    if (myrank == 0) then 
      call MessageNotify( "M", "DampSponge_Init", "EFTime = %f", d=(/EFTime/) )
      call MessageNotify( "M", "DampSponge_Init", "DampDepthH = %f", d=(/DampDepthH/) )
      call MessageNotify( "M", "DampSponge_Init", "DampDepthV = %f", d=(/DampDepthV/) )  
    end if
  end subroutine DampSponge_Init
  

!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_xz(xz_VarA, xz_VarB, DelTime)
    !
    ! ss �ʻ������Ф��륹�ݥ���
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(inout):: xz_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: xz_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !���ݥ��ؤˤ�����ԥ󥰤�׻�
    xz_Var  =  xz_VarA - ( xz_DampRateH + xz_DampRateV ) * xz_VarB * DelTime
    xz_VarA = xz_Var
    
  end subroutine DampSponge_xz


!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_xr(xr_VarA, xr_VarB, DelTime)
    !
    ! sf �ʻ������Ф��륹�ݥ���
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(inout):: xr_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: xr_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: xr_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !���ݥ��ؤˤ�����ԥ󥰤�׻�  
    xr_Var  = xr_VarA - ( xr_DampRateH + xr_DampRateV )* xr_VarB * DelTime
    xr_VarA = xr_Var
    
  end subroutine DampSponge_xr
  

!!!------------------------------------------------------------------------!!!
  subroutine DampSponge_pz(pz_VarA, pz_VarB, DelTime)
    !
    ! fs �ʻ������Ф��륹�ݥ���
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(inout):: pz_VarA(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: pz_VarB(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in)   :: DelTime
    real(DP)               :: pz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !���ݥ��ؤˤ�����ԥ󥰤�׻�  
    pz_Var  = pz_VarA - ( pz_DampRateH + pz_DampRateV ) * pz_VarB * DelTime
    pz_VarA = pz_Var
    
  end subroutine DampSponge_pz


!!!------------------------------------------------------------------------!!!
  function xz_Sponge( xz_Var )
    !
    ! ss �ʻ������Ф��륹�ݥ���
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP)             :: xz_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: xz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !���ݥ��ؤˤ�����ԥ󥰤�׻�
    xz_Sponge = - ( xz_DampRateH + xz_DampRateV ) * xz_Var
    
  end function xz_Sponge
  

!!!------------------------------------------------------------------------!!!
  function xr_Sponge( xr_Var )
    !
    ! sf �ʻ������Ф��륹�ݥ���
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP)             :: xr_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: xr_Var(DimXMin:DimXMax, DimZMin:DimZMax)

    !���ݥ��ؤˤ�����ԥ󥰤�׻�  
    xr_Sponge = - ( xr_DampRateH + xr_DampRateV ) * xr_Var
    
  end function xr_Sponge
  
!!!------------------------------------------------------------------------!!!
  function pz_Sponge( pz_Var )
    !
    ! fs �ʻ������Ф��륹�ݥ���
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,    only: DP
    use gridset,     only: DimXMin,       &! x ����������β���
      &                    DimXMax,       &! x ����������ξ��
      &                    DimZMin,       &! z ����������β���
      &                    DimZMax         ! z ����������ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP)             :: pz_Sponge(DimXMin:DimXMax, DimZMin:DimZMax)
    real(DP), intent(in) :: pz_Var(DimXMin:DimXMax, DimZMin:DimZMax)
    
    !���ݥ��ؤˤ�����ԥ󥰤�׻�  
    pz_Sponge = - ( pz_DampRateH + pz_DampRateV ) * pz_Var 
    
  end function pz_Sponge
  
end module Damping
