!---------------------------------------------------------------------
!     Copyright (C) GFD Dennou Club, 2004. All rights reserved.
!---------------------------------------------------------------------
!= Module CloudPhys_marscond
!
!   * Developer: KITAMORI Taichi, YAMASHITA Tatsuya, SUGIYAMA Ko-ichiro
!   * Version: $ 
!   * Tag Name: $Name:  $
!   * Change History: 


module Cloudphys_MarsCond
  !
  ! �Ȼ��ˤ�äƱ�γ����Ĺ������ˤ�������Ĺ���������. 
  ! ��ȯ�̤�¸�ߤ�������̤���¿�����, ��ȯ�̤�¸�ߤ�������̤ˤ���. 
  ! ��̩�٤����ȯŸ(�ŷ���ʬ)��׻�. 
  ! �������˾������ʤä��Ȥ��˾�ȯ�������ʤ��ʤäƤ��ޤäƤ����Τ�,
  ! ���ʬ�����ľ����.

  !�⥸�塼��ƤӽФ�
  use dc_types,   only: DP, STRING

  !���ۤη�����ػ�
  implicit none

  !�ѿ����
  real(DP), save, private :: DensIce     = 1.565d3  ! �����̩�� [kg/m^3]
  real(DP), save, private :: NumAerosol  = 0.0d0  ! ��������ο�̩�� [1/kg]
  real(DP), save, private :: RadiAerosol = 0.0d0  ! ��������ο�̩�� [1/kg]
  real(DP), save, private :: Kd          = 0.0d0  ! �絤��Ǯ��Ƴ���� [W/K m]
  real(DP), save, private :: SatRatioCr  = 0.0d0  ! �׳�˰���� []
  real(DP), save, private :: SatRtWetAdia = 0.0d0 ! ������Ǯ����˰���� []
  real(DP), save, private :: CO2LatHeat  = 0.0d0  ! ñ�̼��̤�����ζŷ�Ǯ [J/kg]
  real(DP), save, private :: AntA        = 27.4d0 
  real(DP), save, private :: AntB        = 3103.0d0
  real(DP), save, private :: Pi          = 3.1415926535897932385d0 ! �߼�Ψ
  real(DP), save, private :: CDensCr     = 5.0d-5

  !��������
  public cloudphys_marscond_init
  public cloudphys_marscond_forcing

contains
  
  subroutine cloudphys_marscond_init
    !
    ! NAMELIST ����ʪ��������ɤ߹���
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_iounit,     only : FileOpen
    use dc_message,    only : MessageNotify
    use constants,     only : GasRDry           ! �������
    use namelist_util, only : namelist_filename

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    integer       :: unit

    !-------------------------------------------------------------
    ! NAMELIST
    !
    NAMELIST /cloudphys_marscond_nml/         &
      & DensIce, NumAerosol, RadiAerosol, Kd, & 
      & SatRatioCr, SatRtWetAdia, CDensCr
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=cloudphys_marscond_nml)
    close(unit)    

    !-------------------------------------------------------------
    ! ��Ǯ. ����ȸ��ʤ���. 
    !   ���饦����������ڥ����δط��� dp/dT = L p / (R T^2) �� 
    !   Antoine �μ� \ln p = A - B/T ������. 
    !
    CO2LatHeat = AntB * GasRDry

    !-------------------------------------------------------------
    ! ����
    !
    call MessageNotify( "M", &
      & "cloudset_init", "DensIce = %f",  d=(/DensIce/) )
    call MessageNotify( "M", &
      & "cloudset_init", "NumAerosol = %f",  d=(/NumAerosol/) )
    call MessageNotify( "M", &
      & "cloudset_init", "RadiAerosol = %f",  d=(/RadiAerosol/) )
    call MessageNotify( "M", &
      & "cloudset_init", "Kd = %f",  d=(/Kd/) )
    call MessageNotify( "M", &
      & "cloudset_init", "SatRatioCr = %f",  d=(/SatRatioCr/) )
    call MessageNotify( "M", &
      & "cloudset_init", "SatRtWetAdia = %f",  d=(/SatRtWetAdia/) )

    !------------------------------------------------------------
    ! �ѿ����
    !
    call cloudphys_marscond_historyauto

  end subroutine cloudphys_marscond_init
  

!!!----------------------------------------------------
  subroutine cloudphys_marscond_forcing(  &
    &  xyz_PTempNs,         &  !(in) ����
    &  xyz_ExnerNs,         &  !(in) �������ʡ��ؿ�
    &  xyz_CDensNs,         &  !(in) 
    &  xyz_DPTempDtNl,      &  !(in)    
    &  xyz_DExnerDtNl,      &  !(in)    
    &  xyz_DCDensDtNl,      &  !(in)    
    &  xyz_PTempAs,         &  !(out) 
    &  xyz_CDensAs,         &  !(out) ��̩��
    &  xyz_DExnerDtNs       &  !(out) 
    & )
    !
    ! tendency �η׻� & ���̤ȱ�̩�٤ι���
    !
          
    !�⥸�塼��ƤӽФ�
    use dc_types,      only: DP
    use gtool_historyauto,                    &
      &               only : HistoryAutoPut
    use gridset,      only : imin, imax,      &
      &                      jmin, jmax,      &
      &                      kmin, kmax,      &
      &                      nx, ny, nz
    use basicset,     only : xyz_PTempBZ,     &! ���̴��ܾ�
      &                      xyz_ExnerBZ,     &! ̵��������
      &                      xyz_VelSoundBZ,  &! ��®
      &                      xyz_VPTempBZ,    &! ������
      &                      xyz_DensBZ        ! ̩��
    use constants,    only : GasRDry,         & ! �������
      &                      PressBasis,      & ! ���̤δ�వ��
      &                      CpDry              ! �갵��Ǯ
    use SetMargin,    only : SetMargin_xyz
    use timeset,      only : TimeN, DelTimeShort


    ! ���ۤη������ػ�
    implicit none
    
    ! Input
    real(DP), intent(in)   :: xyz_PTempNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! ̵�������� [1]
    real(DP), intent(in)   :: xyz_ExnerNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 
    real(DP), intent(in)   :: xyz_CDensNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! 
    real(DP), intent(in)   :: xyz_DPTempDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(in)   :: xyz_DExnerDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(in)   :: xyz_DCDensDtNl(imin:imax, jmin:jmax, kmin:kmax) 

    real(DP), intent(out)  :: xyz_PTempAs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! ���� [K]
    real(DP), intent(out)  :: xyz_CDensAs(imin:imax, jmin:jmax, kmin:kmax)   
                                        ! ����̩��   [kg/m^3]
    real(DP), intent(out)  :: xyz_DExnerDtNs(imin:imax, jmin:jmax, kmin:kmax) 
                                        ! �ŷ�Ǯ�ˤ�벹���Ѳ�Ψ [K/s]
  
    ! Work
    real(DP)               :: xyz_RadiCloud(imin:imax, jmin:jmax, kmin:kmax)
                                        ! ��γ��Ⱦ�� [m]
    real(DP)               :: xyz_SatRatio(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Rh(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_TempAll(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Mcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Qcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_McondTmp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_PIcond(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_Zero(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_PTempTmp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)               :: xyz_CDensTmp(imin:imax, jmin:jmax, kmin:kmax)
    logical                :: xyz_Mask(imin:imax, jmin:jmax, kmin:kmax)
    integer                :: i,j,k        ! �롼���ѿ�
 

!    ! ����ȯŸ
!    !
!    xyz_PTempTmp = xyz_PTempNs + DelTimeShort * xyz_DPTempDtNl
!    xyz_CDensTmp = xyz_CDensNs + DelTimeShort * xyz_DCDensDtNl

!    ! Set Margin
!    !
!    call SetMargin_xyz(xyz_PTempTmp)
!    call SetMargin_xyz(xyz_CDensTmp)

!    ! ��ή����ˤʤä���ʬ������
!    ! 
!    call FillNegativeDensity(xyz_CDensTmp)

!    ! Set Margin
!    !
!    call SetMargin_xyz(xyz_CDensTmp)

    ! tendency �� Ns ���ͤǷ׻�
    !
    xyz_TempAll = (xyz_ExnerBZ + xyz_ExnerNs) * (xyz_PTempBZ + xyz_PTempNs)
    xyz_Mask = .false. 
    xyz_Zero = 0.0d0
   
    ! Ǯ͢���˴ؤ��뷸��
    ! 
    xyz_Rh = (CO2LatHeat**2.0d0) / (Kd * GasRDry * (xyz_TempAll**2.0d0))

    ! ˰���� (1.36) ��. 
    !
    xyz_SatRatio =                                                   &
      &  PressBasis * (xyz_ExnerBZ + xyz_ExnerNs)**(CpDry / GasRDry) &
      &  / exp( AntA - AntB / xyz_TempAll )  

    ! ��γȾ��. 
    !
    xyz_RadiCloud =   &
      & (  &
      &     RadiAerosol**3.0d0  &
      &   + 3.0d0 * xyz_CDensNs / (4.0d0 * Pi * DensIce * xyz_DensBZ * NumAerosol) & 
      &  ) ** (1.0d0 / 3.0d0)

    ! �ŷ��̤η׻�
    !
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          if (xyz_SatRatio(i,j,k) - SatRatioCr > epsilon(0.0d0)) then
            ! ˰���椬��᤿�׳�˰��������礭�����
            ! CO2 ��˰����μ����ϵ��絤���Ѥ�����˰����μ��������ʤ�ΤʤΤ���ǧ����ɬ�פ�����. 
            !
            xyz_Mask(i,j,k) = .true. 
            
          else if (xyz_CDensNs(i,j,k) > CDensCr ) then
            ! �׳�˰�����Ķ���Ƥ��ʤ���, ��̩�٤����Ͱʾ�Ǥ�����
            ! ˰���椬 1 �ʾ�ΤȤ��˶ŷ�, 1 �ʲ��ΤȤ��˾�ȯ. 
            !
            xyz_Mask(i,j,k) = .true. 
            
          else if (xyz_CDensNs(i,j,k) /= 0.0d0 .and. xyz_SatRatio(i,j,k) < 1.0d0 ) then
            ! ��̩�٤�����̤����, ˰���� 1.0 ̤���Ǥ�����˾�ȯ. 
            !
            xyz_Mask(i,j,k) = .true. 

          end if
        end do
      end do
    end do

    ! �ŷ��̤β��ͤ�׻�. 
    !
    xyz_McondTmp =                                                &
      & max( - xyz_CDensTmp / DelTimeShort,                       &
      &      4.0d0 * Pi * xyz_RadiCloud * xyz_DensBZ * NumAerosol &
      &        / xyz_Rh * (xyz_SatRatio - 1.0d0)                  &
      &   )

    ! �ŷ��̤�׻�. Mask �� .false. �����Ǥˤϥ���������
    !
    xyz_Mcond = merge(xyz_McondTmp, xyz_Zero, xyz_Mask)

    ! ̩�٤η׻� 
    !
    xyz_CDensAs = xyz_CDensTmp + DelTimeShort * xyz_Mcond

    ! Set Margin
    !
    call SetMargin_xyz(xyz_CDensAs)

    ! �����Ѳ��η׻�
    ! 
    xyz_Qcond = CO2LatHeat * xyz_Mcond / (CpDry * xyz_DensBZ * xyz_ExnerBZ)
    xyz_PTempAs = xyz_PTempTmp + DelTimeShort * xyz_QCond

    ! Set Margin
    !
    call SetMargin_xyz(xyz_PTempAs)
    
    ! �������ʡ��ؿ��λ�����ʬ
    ! 
    xyz_PIcond =                                                   &
      &  ( xyz_VelSoundBZ ** 2.0d0 )                               &
      &     / (CpDry * xyz_DensBZ * (xyz_VPTempBZ ** 2.0d0))        &
      &     * xyz_Mcond                                            &
      &     * ( CO2LatHeat / (CpDry * xyz_ExnerBZ) - xyz_VPTempBZ ) 

    xyz_DExnerDtNs = xyz_DExnerDtNl + xyz_PIcond
    
    call HistoryAutoPut(TimeN, 'PTempCond', xyz_Qcond(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'ExnerCond', xyz_PIcond(1:nx, 1:ny, 1:nz))
    call HistoryAutoPut(TimeN, 'CDensCond', xyz_Mcond(1:nx, 1:ny, 1:nz))

    call SetMargin_xyz(xyz_PTempAs)
    call SetMargin_xyz(xyz_CDensAs)
!    call SetMargin_xyz(xyz_DExnerDtNs)

  end subroutine Cloudphys_marscond_forcing


  subroutine Cloudphys_marscond_historyauto
    !
    ! tendency �ν�������
    !

    !�⥸�塼��ƤӽФ�
    use gtool_historyauto, only: HistoryAutoAddVariable

    !���ۤη�����ػ�
    implicit none

    !------------------------------------------------------
    ! tendency �����
    !
    call HistoryAutoAddVariable(  &
      & varname='PTempCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of exner function', &
      & units='s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='CDensCond',&
      & dims=(/'x','y','z','t'/),     &
      & longname='Latent heat term of cloud density', &
      & units='K.m-3.s-1',    &
      & xtype='float')

  end subroutine Cloudphys_marscond_historyauto
  
end module Cloudphys_MarsCond
