!= Module initialdata_takemi2007
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_takemi2007.f90,v 1.8 2014/07/08 00:59:09 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_takemi2007
  !
  ! Takemi (2007) ���Ϥ����δ��ܾ졦���������ꤹ��

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP

  !���ۤη�����ػ�
  implicit none

  !�����ѿ������
  real(DP), parameter, private :: PressSfcTakemi = 1.0d5   ! ��ɽ�ΰ���
  real(DP), parameter, private :: PTempSfcTakemi = 300.0d0 ! ��ɽ�β���
  real(DP), parameter, private :: AltTr    = 1.2d4   ! ��ή�����̹���
  real(DP), parameter, private :: HumMin   = 0.25d0  ! ���������ʹ���
  integer,  parameter, private :: SpcID = 6          ! ����ֹ�
  real(DP), save, private      :: QMixSfc            ! ��ɽ�̤ǤΥ����
  real(DP), save, private      :: VelXSfc            ! ��ɽ��®��
  real(DP), save, private      :: PTempTr = 0.0d0    ! 
  real(DP), save, private      :: DryFact = 0.0d0    ! 
  real(DP), save, private      :: Alt1    = 0.0d0    ! 
  real(DP), save, private      :: Alt2    = 0.0d0    ! 
  integer,  save, private      :: amin, amax

  !�������������
  public  initialdata_takemi2007_init
  public  initialdata_takemi2007_basic
  public  initialdata_takemi2007_wind

contains

!!!------------------------------------------------------------------------------!!!
  subroutine initialdata_takemi2007_init
    !
    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    !

    !�⥸�塼���ɤ߹���
    use dc_types,     only: STRING, DP
    use dc_iounit,    only: FileOpen 
    use dc_message,   only: MessageNotify
    use namelist_util,only: namelist_filename
    use gridset,      only: nz
    use axesset,      only: z_Z                !�����顼�ʻ����Ǥι���
    use constants,    only: PressBasis,       &!���̤δ�వ��
      &                     TempSfc,          &!��ɽ�̲���
      &                     PressSfc           !��ɽ�̰���
    
    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP)            :: Alt = 0.0d0  !����
    integer             :: unit         !����ե������������ֹ�  
    integer             :: k
    integer             :: ID_BasicZ = 0
    integer             :: ID_Wind   = 0
    integer, parameter  :: ID_MidLat_Q10     = 1
    integer, parameter  :: ID_MidLat_Q12     = 2
    integer, parameter  :: ID_MidLat_Q14     = 3
    integer, parameter  :: ID_MidLat_Q16     = 4
    integer, parameter  :: ID_MidLat_Q16DRY1 = 5
    integer, parameter  :: ID_MidLat_Q16DRY2 = 6
    integer, parameter  :: ID_MidLat_Q18     = 7
    integer, parameter  :: ID_Tropic_Q18     = 8
    integer, parameter  :: ID_Tropic_Q18DRY1 = 9
    integer, parameter  :: ID_Tropic_Q18DRY2 = 10
    integer, parameter  :: ID_Tropic_Q18DRY3 = 11
    integer, parameter  :: ID_Wind_LowLevel    = 1
    integer, parameter  :: ID_Wind_MiddleLevel = 2
    integer, parameter  :: ID_Wind_HighLevel   = 3

    character(STRING)   :: FlagEnv = ""
    character(STRING)   :: FlagWind = ""


    !����ե����뤫���ɤ߹�����ϥե��������
    NAMELIST /initialdata_takemi2007_nml/ FlagEnv, FlagWind, VelXSfc
    
    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_takemi2007_nml)
    close(unit)

    ! ��ǧ. ��ɽ�̲��١����Ϥ���ꤹ��Τ��̥⥸�塼��ʤΤ�. 
    !
    if (     PressBasis /= PressSfcTakemi &
      & .OR. PressSfc /= PressSfcTakemi   &
      & .OR. TempSfc  /= PTempSfcTakemi  ) then 
      
      call MessageNotify( "E", "initaldata_takemi2007_init", &
        & "Constants are wrong. please PressSfc = 1.0d5, TempSfc = 300.0d0")
    end if
    
    call MessageNotify( "M", "initaldata_takemi2007_init", &
      & "VelXSfc= %f", d=(/VelXSfc/) ) 

    !���ܾ������
    !
    if (FlagEnv == "MidLat_Q10") then 
      ID_BasicZ = ID_MidLat_Q10
    elseif (FlagEnv == "MidLat_Q12") then 
      ID_BasicZ = ID_MidLat_Q12
    elseif (FlagEnv == "MidLat_Q14") then 
      ID_BasicZ = ID_MidLat_Q14
    elseif (FlagEnv == "MidLat_Q16") then 
      ID_BasicZ = ID_MidLat_Q16
    elseif (FlagEnv == "MidLat_Q16DRY1") then 
      ID_BasicZ = ID_MidLat_Q16DRY1
    elseif (FlagEnv == "MidLat_Q16DRY2") then 
      ID_BasicZ = ID_MidLat_Q16DRY2
    elseif (FlagEnv == "MidLat_Q18") then 
      ID_BasicZ = ID_MidLat_Q18
    elseif (FlagEnv == "Tropic_Q18") then 
      ID_BasicZ = ID_Tropic_Q18
    elseif (FlagEnv == "Tropic_Q18DRY1") then 
      ID_BasicZ = ID_Tropic_Q18DRY1
    elseif (FlagEnv == "Tropic_Q18DRY2") then 
      ID_BasicZ = ID_Tropic_Q18DRY2
    elseif (FlagEnv == "Tropic_Q18DRY3") then 
      ID_BasicZ = ID_Tropic_Q18DRY3
    end if

    ! ®�پ������
    !
    if (FlagWind == "LowLevel") then 
      ID_Wind = ID_Wind_LowLevel
    elseif (FlagWind == "MiddleLevel") then 
      ID_Wind = ID_Wind_MiddleLevel
    elseif (FlagWind == "HighLevel") then 
      ID_Wind = ID_Wind_HighLevel
    end if

    ! ���̾�
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q10:ID_MidLat_Q18)
      PTempTr = 343.0d0
    case (ID_Tropic_Q18:ID_Tropic_Q18DRY3)
      PTempTr = 358.0d0
    end select

    ! ������
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q10) 
      QMixSfc = 0.010d0
    case (ID_MidLat_Q12) 
      QMixSfc = 0.012d0
    case (ID_MidLat_Q14) 
      QMixSfc = 0.014d0
    case (ID_MidLat_Q16:ID_MidLat_Q16DRY2) 
      QMixSfc = 0.016d0
    case (ID_MidLat_Q18:ID_Tropic_Q18DRY3)
      QMixSfc = 0.018d0
    end select

    ! �����Ѳ�
    !
    select case (ID_BasicZ)
    case (ID_MidLat_Q16DRY1)
      DryFact = - 0.13d0
      Alt = 2.5d3
    case (ID_MidLat_Q16DRY2)
      DryFact = - 0.30d0
      Alt = 2.5d3
    case (ID_Tropic_Q18DRY1)
      DryFact = - 0.20d0
      Alt = 2.5d3
    case (ID_Tropic_Q18DRY2)
      DryFact = - 0.20d0
      Alt = 5.0d3
    case (ID_Tropic_Q18DRY3)
      DryFact = - 0.20d0
      Alt = 7.5d3
    end select

    do k = 1, nz
      if (z_Z(k) < AltTr .AND. AltTr <= z_Z(k+1)) then 
        amax = k
      end if
      if (z_Z(k) < Alt .AND. Alt <= z_Z(k+1)) then 
        amin = k
      end if
    end do

    ! ��ʿ��®
    !
    select case (ID_Wind)
    case (ID_Wind_LowLevel) 
      Alt1 = 0.0d0
      Alt2 = 2.5d3
    case (ID_Wind_MiddleLevel)
      Alt1 = 2.5d3
      Alt2 = 5.0d3
    case (ID_Wind_HighLevel)
      Alt1 = 5.0d3
      Alt2 = 7.5d3
    end select

  end subroutine initialdata_takemi2007_init


!!!------------------------------------------------------------------------------!!!
  subroutine  initialdata_takemi2007_basic( z_Temp, z_Press, zf_MolFr )
    !
    !== ����
    ! * deepconv ���ϵ��ѤΥƥ��ȷ׻��Ȥ���Takemi(2007)�κƸ��׻���
    !   ���뤿��μ��٤δ��ܾ���������
    !   * ���ܾ�β��٤μ������̤�Ϳ�����Ƥ��뤿��, ���٤��Ѵ�����ɬ�פ�����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,     only: DP
    use gridset,      only: kmin, kmax,       &!���󥵥��� (Z ����)
      &                     ncmax,            &!�Ž���ʬ�ο�
      &                     nz                 !ʪ���ΰ���礭�� (Z����)
    use axesset,      only: z_Z,              &!�����顼�ʻ����Ǥι���
      &                     dz                 !��ľ�ʻҴֳ�
    use constants,    only: PressBasis,       &!���̤δ�వ��
      &                     GasRDry,          &!������ʬ���갵��Ǯ
      &                     CpDry,            &!������ʬ���갵��Ǯ
      &                     Grav,             &!���ϲ�®��
      &                     TempSfc,          &!��ɽ�̲���
      &                     PressSfc,         &!��ɽ�̰���
      &                     MolWtDry  
    use composition,  only: MolWtWet
    use chemcalc,     only: SvapPress 

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(out):: z_Press(kmin:kmax)           !����
    real(DP), intent(out):: z_Temp(kmin:kmax)            !����
    real(DP), intent(out):: zf_MolFr(kmin:kmax, 1:ncmax) !�����
    real(DP) :: z_Hum(kmin:kmax)                         ! ���м���
    real(DP) :: z_PTemp(kmin:kmax)                       ! ����
    real(DP) :: QMix
    integer  :: k

    !-------------------------------------------
    ! �����
    !
    z_Temp   = 1.0d-60
    z_PTemp  = 1.0d-60
    z_Press  = 1.0d-60
    z_Hum    = 0.0d0
    zf_MolFr = 0.0d0

    !----------------------------------------------
    ! ����
    !   ��ʸ��μ� (2) ���׻�
    !    
    do k = 1, nz
      if (z_Z(k) <= AltTr) then 
        z_Hum(k) = 1.0d0 - 0.75d0 * (z_Z(k) / AltTr) ** 1.25d0
      elseif (z_Z(k) > AltTr) then 
        z_Hum(k) = HumMin
      end if
    end do

    ! Fig.2b �����餫�����м��٤� 95% ���٤��Ǥ��ߤ�ˤʤäƤ���Τ�, 
    ! ����ͤ��ߤ��Ƥߤ�
    !
    where (z_Hum > 0.95d0) 
      z_Hum = 0.95d0
    end where
    
    ! DRY �������ξ��. 
    !   DRY �ʳ��Ǥ�, DryFact = 0.0 �ˤʤäƤ���. 
    !   ���м��٤β����ͤ򲼲��ʤ��褦Ĵ�����Ƥ���.
    !
    do k = 1, nz
      if (amin < k .AND. k <= amax) then
        z_Hum(k) = z_Hum(k) + DryFact
        if (z_Hum(k) <= HumMin) then 
          z_Hum(k) = HumMin
        end if
      end if
    end do

    !----------------------------------------------
    ! ����, ����, ����
    !   ���̤�, ��ʸ��μ� (1) ���׻�
    !
    !    
    z_PTemp(1) = TempSfc + (PTempTr - TempSfc) * ((z_Z(1) / AltTr) ** 1.25d0)
    z_Press(1) = PressSfc - (Grav * PressSfc * dz * 5.0d-1) / (GasRDry * TempSfc)
    z_Temp(1)  = z_PTemp(1) * (z_Press(1) / PressBasis) ** (GasRDry / CpDry)

    do k = 2, nz
      if (k <= amax) then 
        z_PTemp(k) = TempSfc + (PTempTr - TempSfc) * ((z_Z(k) / AltTr) ** 1.25d0)
      elseif (k > amax) then 
        z_PTemp(k) = z_PTemp(amax) * exp( Grav * (z_Z(k) - AltTr) / (CpDry * z_Temp(amax)))
      end if

      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * dz) / (GasRDry * z_Temp(k-1))
      
      z_Temp(k)  = z_PTemp(k) * ((z_Press(k) / PressBasis) ** (GasRDry / CpDry))
    end do
    
    !----------------------------------------------
    ! �����
    !   ��ʸ��ˤ�, in the lowest 1.5km depth �Ǻ�������ͤ����ꤹ��Ȥ���
    !   ���Ҥ����ä���, Fig2b ����Ƚ�Ǥ����, �� (2) ��Ϳ���������м��٤���
    !   �������׻���, ���줬��ɽ�Ǥκ����� (nml ��Ϳ����) ��Ķ������, 
    !   ��ɽ�Ǥκ������Ʊ���Ȥ��Ƥ���Τ���. 
    !
    do k = 1, nz
      zf_MolFr(k,1) = SvapPress(SpcID, z_Temp(k)) * z_Hum(k) / z_Press(k)
      QMix = zf_MolFr(k,1) / MolWtDry * MolWtWet(1)
      
      if (QMix > QMixSfc) then 
        zf_MolFr(k,1) = QMixSfc * MolWtDry / MolWtWet(1)
      end if
    end do
    
  end subroutine Initialdata_takemi2007_basic


  subroutine initialdata_takemi2007_wind(pyz_VelX)
    !
    !-------------------------------------------------------------!
    ! ������������ (Takemi,2007)                                  !
    !-------------------------------------------------------------!
    !
    != ����
    !* case "Takemi2007" �Ǥη׻����˱�ľ�������Τ�������Ϳ������˻��Ѥ���
    !* ����Ϳ�����ˤ�, �ʲ��Τ褦�ʥХꥨ������󤬤���
    !  (1) ��������Ϳ������٤��Ѥ���
    !  (2) �������Τ������κ�����® (U_s) ���Ѥ���
    !
    !  (1) �ˤĤ��Ƥ�, (a) 0 - 2.5 km, (b) 2.5 - 5.0 km, (c) 5.0 - 7.5 km ��
    !  ���ѥ����󤬤���
    !  (2) �ˤĤ��Ƥ�, Takemi (2007) �Ǥ�Ǯ�Ӿ������پ�β��پ����
    !  �ۤʤ��ͤ����ꤷ�Ƥ���
    !
    !  ���ζ���(Us)��, �ʲ����̤�
    !  <Ǯ�Ӿ�>   (1) 5 m/s, (2) 10 m/s, (3) 15 m/s
    !  <����پ�> (1) 10 m/s, (2) 15 m/s, (3) 20 m/s
    !
    !* �������η����ϼ��� (Takemi, 2007)   |
    !                                     /| 7.5 km
    !                                    / |
    !                                   /  |
    !                                  / ��|
    !                                 ��  /| 5.0 km
    !                                 �� / |
    !                                 ��/  |
    !                                  / ��|
    !                                 ��  /| 2.5 km
    !                                 �� / |
    !                                 ��/  |
    !                                  / ��|
    !---------------------------------+------------- 0.0 km
    !                                Us (m/s)
    !

    !�⥸�塼���ɤ߹���
    use dc_types,     only: DP
    use gridset,      only: imin, imax,       &!����� X �����ξ��
      &                     jmin, jmax,       &!����� Y �����ξ��
      &                     kmin, kmax,       &!����� Z �����ξ��
      &                     nz
    use axesset,      only: z_Z                !�����顼�ʻ����Ǥι���

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(out) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    integer               :: k
    
    !�����
    pyz_VelX = 0.0d0
    
    do k = 1, nz
      if (z_Z(k) <= Alt1) then 
        pyz_VelX(:,:,k) = - VelXSfc
      elseif (z_Z(k) > Alt1 .AND. z_Z(k) <= Alt2) then 
        pyz_VelX(:,:,k) = - VelXSfc + (VelXSfc / (Alt2 - Alt1)) * (z_Z(k) - Alt1)
      end if
    end do
    
  end subroutine initialdata_takemi2007_wind
  
end module initialdata_takemi2007
