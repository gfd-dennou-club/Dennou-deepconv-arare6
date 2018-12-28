!= Module HeatFluxBS1998
!
! Authors::   KAWABATA Takuya
! Version::   $Id: surfaceflux_baker1998.f90,v 1.4 2014/03/04 07:43:44 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2012. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Surfaceflux_Baker1998
  !
  ! Baker et al. 1989 �ξ岼�����Ǥ���ե�å����׻��⥸�塼��.
  !

  !�⥸�塼���ɤ߹���
  !
  use dc_types, only: DP, STRING

  !���ۤη�����ػ�
  !
  implicit none

  !°���λ���
  !
  private

  !�ѿ����
  !
  real(DP), save, allocatable :: xyz_DPTempDtFlux(:,:,:)
  real(DP), save, allocatable :: xyz_DExnerDtFlux(:,:,:)

  character(STRING), parameter:: module_name = 'surfaceflux_baker1998'
                              ! �⥸�塼���̾��.
                              ! Module name

  ! �ؿ��� public ������
  !
  public surfaceflux_baker1998_init
  public surfaceflux_baker1998_forcing

contains
!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_BAKER1998_init
    !
    ! Baker et al. 1989 �˽��ä�, �岼�����Ǥ�Ǯ�ե�å��������. 
    !

    ! �⥸�塼��ƤӽФ�
    use dc_types,          only : DP
    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable    
    use gridset,           only : imin,         & !x ����������β���
      &                           imax,         & !x ����������ξ��
      &                           jmin,         & !y ����������β���
      &                           jmax,         & !y ����������ξ��
      &                           kmin,         & !z ����������β���
      &                           kmax,         & !z ����������ξ��
      &                           nz
    use axesset,           only : dz              !z �����γʻ����ֳ�
    use basicset,          only : xyz_ExnerBZ,  & !�������ʡ��ؿ��δ��ܾ�   
      &                           xyz_DensBZ      !���ܾ��̩��
    use constants,         only : CpDry
    use axesset,           only : z_Z
    use DExnerDt,          only : xyz_DExnerDt_xyz
    
    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    real(DP), parameter :: Qsub_100km = 717.0d0  ! [W/m^2]
    real(DP), parameter :: Zref = 100.0d3
    real(DP)            :: QsubTop, QsubBtm
    real(DP)            :: z1, z2
    real(DP)            :: xyz_TempFlux(imin:imax,jmin:jmax,kmin:kmax)  !���٥ե�å��� [W/m^3]
    real(DP)            :: xyz_PTempFlux(imin:imax,jmin:jmax,kmin:kmax) !���̥ե�å��� [W/m^3]
    integer             :: k

    !---------------------------------------------------------------       
    ! �����
    !
    z1 = z_Z(nz) + dz       ! �ΰ�������ϴޤޤʤ�
    z2 = z_Z(1)             ! �ΰ�������ϴޤ�

    allocate( xyz_DPTempDtFlux(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtFlux(imin:imax,jmin:jmax,kmin:kmax) )
    
    xyz_DPTempDtFlux = 0.0d0
    xyz_DExnerDtFlux = 0.0d0
    xyz_TempFlux     = 0.0d0
    xyz_PTempFlux    = 0.0d0

    !---------------------------------------------------------------       
    ! �����Ǥβ�ǮΨ��׻�����. 
    ! 
    QsubTop = 0.0d0
    do while ( z1 <= Zref )
      QsubTop = QsubTop + cal_Qsub( z1 ) * dz
      z1 = z1 + dz
    end do
    xyz_TempFlux(:,:,nz)    = - Qsub_100km + QsubTop 
    xyz_TempFlux(:,:,nz+1)  = - Qsub_100km + QsubTop  ! ������
    xyz_PTempFlux(:,:,nz)   = - Qsub_100km + QsubTop 
    xyz_PTempFlux(:,:,nz+1) = - Qsub_100km + QsubTop  ! ������

    QsubBtm = 0.0d0
    do while ( z2 <= Zref )
      QsubBtm = QsubBtm + cal_Qsub( z2 ) * dz
      z2 = z2 + dz
    end do
    xyz_TempFlux(:,:,1) = Qsub_100km - QsubBtm
    xyz_TempFlux(:,:,0) = Qsub_100km - QsubBtm  ! ������
    
    QsubBtm = 0.0d0
    do k = 1, nz
      QsubBtm = QsubBtm + cal_Qsub( z_Z(k) ) / xyz_ExnerBZ(1,1,k) * dz
    end do
    xyz_PTempFlux(:,:,1) = Qsub_100km - QsubTop - QsubBtm
    xyz_PTempFlux(:,:,0) = Qsub_100km - QsubTop - QsubBtm

    !---------------------------------------------------------------
    ! �����Ѳ�Ψ�η׻�
    !   �Ǿ�������Ӻǲ��ؤγʻҤ� xyz_TempFlux ��Ϳ�������Ǯ [W/m^2] 
    !   ��Ŭ�Ѥ����褦�ˤ���. ���Τ���, �ʻҴֳ� dz �ǳ�껻����. 
    !
!    xyz_DPTempDtFlux = xyz_TempFlux / xyz_ExnerBZ / xyz_DensBZ / CpDry / dz 
    xyz_DPTempDtFlux = xyz_PTempFlux / xyz_DensBZ / CpDry / dz 
    xyz_DExnerDtFlux = xyz_DExnerDt_xyz( xyz_DPTempDtFlux )

    !---------------------------------------------------------------
    ! ��ǧ�Ѥν���
    !
    call MessageNotify( "M", module_name, "TempFluxTop = %f", &
      &                  d=(/xyz_TempFlux(1,1,nz)/) )
    call MessageNotify( "M", module_name, "PTempFluxTop = %f", &
      &                  d=(/xyz_PTempFlux(1,1,nz)/) )
    call MessageNotify( "M", module_name, "DPTempDtTop = %f", &
      &                  d=(/xyz_DPTempDtFlux(1,1,nz)/) )
    
    call MessageNotify( "M", module_name, "TempFluxBtm = %f", &
      &                  d=(/xyz_TempFlux(1,1,1)/) )
    call MessageNotify( "M", module_name, "PTempFluxBtm = %f", &
      &                  d=(/xyz_PTempFlux(1,1,1)/) )
    call MessageNotify( "M", module_name, "DPTempDtBtm = %f", &
      &                  d=(/xyz_DPTempDtFlux(1,1,1)/) )

    !---------------------------------------------------------------
    ! �ե���������
    !
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtSfc',    &
      & dims=(/'x','y','z','t'/), &
      & longname='potential temperature tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='DExnerDtSfc',    &
      & dims=(/'x','y','z','t'/), &
      & longname='exner function tendency by surface flux', &
      & units='K.s-1',            &
      & xtype='float')

  end subroutine Surfaceflux_BAKER1998_init

!!!------------------------------------------------------------------------!!!
  subroutine Surfaceflux_BAKER1998_forcing(  &
    &     xyz_DPTempDt, xyz_DExnerDt      &
    &  )

    ! �⥸�塼��ƤӽФ�
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,         & !x ����������β���
      &                           imax,         & !x ����������ξ��
      &                           jmin,         & !y ����������β���
      &                           jmax,         & !y ����������ξ��
      &                           kmin,         & !z ����������β���
      &                           kmax,         & !z ����������ξ��
      &                           nx, ny, nz

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(inout):: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout):: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! ����������ؤβû�
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtFlux
    xyz_DExnerDt = xyz_DExnerDt + xyz_DExnerDtFlux
    
    ! ����
    ! 
    call HistoryAutoPut(TimeN, 'DPTempDtSfc', xyz_DPTempDtFlux(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DExnerDtSfc', xyz_DExnerDtFlux(1:nx,1:ny,1:nz))

  end subroutine Surfaceflux_BAKER1998_forcing


!!!------------------------------------------------------------------------!!!
  function cal_Qsub( z0 )
    ! Ʊ���롼���� radiation_baker1998.f90 �ˤ�ޤޤ�Ƥ���Τ���դ��뤳��.
    
    use dc_types, only: DP

    implicit none

    real(DP), intent(in) :: z0
    real(DP)             :: cal_Qsub
    real(DP), parameter  :: z_Upper = 6.70d4      ! ����ĺ������
    real(DP), parameter  :: z_Lower = 2.70d4      ! ����ĺ������
    real(DP), parameter  :: c_Upper = 2.70d-2     ! �� 1 ����
    real(DP), parameter  :: c_Lower = 3.60d-3     ! �� 2 ����
    real(DP), parameter  :: s_Upper = 7.50d3      ! �� 1 ɸ���к�
    real(DP), parameter  :: s_Lower = 1.30d4      ! �� 2 ɸ���к�

    cal_Qsub =                                                                      &
      & (                                                                           &
      &     c_Lower * exp( - ( z0 - z_Lower )**2.0d0 / ( 2.0d0 * s_Lower **2.0d0) ) &
      &   + c_Upper * exp( - ( z0 - z_Upper )**2.0d0 / ( 2.0d0 * s_Upper **2.0d0) ) &
      &  ) 
   
    
  end function cal_Qsub
  
end module Surfaceflux_Baker1998
