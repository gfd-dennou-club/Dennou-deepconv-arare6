!= ���ش�Ϣ�ν��̤�׻����뤿��Υ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: chemcalc.f90,v 1.12 2014/07/08 01:05:32 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006-2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module ChemCalc
  !
  != ���ش�Ϣ�ν��̤�׻����뤿��Υ⥸�塼��. 
  !
  ! AMP �� Antoine ��˰�¾����������Ѥ��ưʲ������. 
  ! �ǥե���ȤǤ� AMP ����Ȥ��褦�ˤ��Ƥ���. 
  !  * ˰�¾�����
  !  * ˰�¾������β�����ʬ
  !  * ��Ǯ
  !

  !�⥸�塼��ƤӽФ�
  use dc_types,   only: DP              !���ٻ���
  use ChemData,   only: ChemData_SpcNum !�ǡ����١�����β��ؼ��

  ! ���ۤη�����ػ�
  !
  implicit none

  ! �ѿ������
  !
  real(DP), save, public  :: ReactHeatNH4SH       !NH4SH ����ȿ��Ǯ [J/K kg]
  real(DP), save, public  :: ReactHeatNH4SHPerMol !NH4SH ����ȿ��Ǯ [J/K mol]

  integer,  save, private :: a_kmin(ChemData_SpcNum)  !ʪ����˷�᤿����β���
  integer,  save, private :: a_kmax(ChemData_SpcNum)  !ʪ����˷�᤿����β���
  real(DP), save, private :: a_SwAmp(ChemData_SpcNum) !�����å�. AMP ��Ȥ����� 1.0, �����Ǥʤ���� 0.0  
  real(DP), save, private :: a_SwAnt(ChemData_SpcNum) !�����å�. Antoine ��Ȥ����� 1.0, �����Ǥʤ���� 0.0
  real(DP), save, private :: a_antA(ChemData_SpcNum)  !Antoine �ξ��������� A ����
  real(DP), save, private :: a_antB(ChemData_SpcNum)  !Antoine �ξ��������� B ����
  real(DP), save, private :: a_antC(ChemData_SpcNum)  !Antoine �ξ��������� C ����
  real(DP), save, private :: a_antU(ChemData_SpcNum)  !Antoine �ξ���������ñ�̴����Τ���η���
  real(DP), save, private :: a_ampA(ChemData_SpcNum)  !AMP ���ξ��������� A ����
  real(DP), save, private :: a_ampB(ChemData_SpcNum)  !AMP ���ξ��������� B ����
  real(DP), save, private :: a_ampC(ChemData_SpcNum)  !AMP ���ξ��������� C ����
  real(DP), save, private :: a_ampD(ChemData_SpcNum)  !AMP ���ξ��������� D ����
  real(DP), save, private :: a_ampE(ChemData_SpcNum)  !AMP ���ξ��������� E ����
  real(DP), save, private :: a_MolWt(ChemData_SpcNum) !ʬ����

  ! �������륵�֥롼����� public °�����դ���
  !
  public ChemCalc_Init                 !������롼����
  public MolWt                         !ʬ����
  public GasR                          !�������
  public CpRef, CpPerMolRef, CvRef     !�갵��Ǯ, ������Ǯ
  public SvapPress, xyz_SvapPress      !˰�¾����� [Pa]
  public xyz_LatentHeat                !��Ǯ [J/K kg]
  public LatentHeatPerMol              !��Ǯ [J/K mol]
  public xyz_DQMixSatDPTemp
  public xyz_DelQMixNH4SH
  public DelMolFrNH4SH
 
contains
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ChemCalc_Init
    !
    !������롼����
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP                               !���ٻ���
    use dc_message, only: MessageNotify                    !��å�����ɽ��
    use ChemData,   only: ChemData_init,                 & !�����
      &                   ChemData_SpcNum,               & !�ǡ����١�����β��ؼ��
      &                   ChemData_SvapPress_AntoineA,   & !Antoine ���� A ����
      &                   ChemData_SvapPress_AntoineB,   & !Antoine ���� B ����
      &                   ChemData_SvapPress_AntoineC,   & !Antoine ���� C ����
      &                   ChemData_SvapPress_AntoineUnit,& !ñ���Ѵ��ѷ���
      &                   ChemData_SvapPress_AMPA,       & !AMP ���� A ����
      &                   ChemData_SvapPress_AMPB,       & !AMP ���� B ����
      &                   ChemData_SvapPress_AMPC,       & !AMP ���� C ����
      &                   ChemData_SvapPress_AMPD,       & !AMP ���� D ����
      &                   ChemData_SvapPress_AMPE,       & !AMP ���� E ����
      &                   ChemData_MolWt,                & !ʬ����
      &                   GasRUniv,                      & !�������
      &                   ChemData_SpcSymbol,            & !ʬ��̾
      &                   ChemData_OneSpcID                !���ؼ�� ID ����     
    use gridset,    only: nz,                            & ! ʪ���ΰ���礭��
      &                   kmin, kmax                       ! ����� Z �����ξ�¡�����
    use constants,  only: PressSfc                         !���������Ǥΰ��� [Pa]
    use basicset,   only: xyz_TempBZ,                    & !���٤δ��ܾ�
      &                   xyz_PressBZ                      !���Ϥδ��ܾ�
    use axesset,    only: z_Z                              !z ��
    use namelist_util, &
      &             only: namelist_filename
    use dc_iounit,  only: FileOpen

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    character(20)      :: Name
    integer            :: id
    integer            :: k
    integer            :: unit, ierr
    real(DP)           :: Temp
    real(DP)           :: Press
    real(DP),parameter :: Temp0C = 273.15d0
    real(DP)           :: logsvap
    real(DP)           :: HeightUp = 0.0d0
    real(DP)           :: HeightDown = 0.0d0

    !NAMELIST �����
    NAMELIST /chemcalc_nml/ HeightUp, HeightDown
 
    !�ե����륪���ץ�. �������. 
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=chemcalc_nml, iostat=ierr, err=99)
    close(unit)
99  call MessageNotify( "M", "ChemCalc_Init", "No information of chemcalc_nml in config file; use default values")

    !-----------------------------------------------------------
    ! �����
    !

    ! �ǡ����١����ν����
    call chemdata_init

    !Antoine ��˰�¾��������η���
    a_antA = ChemData_SvapPress_AntoineA
    a_antB = ChemData_SvapPress_AntoineB
    a_antC = ChemData_SvapPress_AntoineC
    a_antU = ChemData_SvapPress_AntoineUnit

    !AMP ����˰�¾��������η���
    a_ampA = ChemData_SvapPress_AMPA
    a_ampB = ChemData_SvapPress_AMPB
    a_ampC = ChemData_SvapPress_AMPC
    a_ampD = ChemData_SvapPress_AMPD
    a_ampE = ChemData_SvapPress_AMPE

    !ʬ����
    a_MolWt = ChemData_MolWt
    
    !NH4SH ��ȿ��Ǯ�ν����
    !  NH4SH 1kg ���Ф���ȿ��Ǯ�ˤ���.
    Name = 'NH4SH-s'
    id   = ChemData_OneSpcID( Name )  
    
    ReactHeatNH4SHPerMol  = GasRUniv * 10834.0d0
    ReactHeatNH4SH = GasRUniv * 10834.0d0 / MolWt( id )

    !--------------------------------------------------------
    ! ʪ���ˤ�ä�, AMP ��Ȥ��� Antoine ��Ȥ�������. 
    ! AMP �η���������ʤ�� Antoine ��Ȥ����Ȥˤ���. 
    !     

    do ID = 1, ChemData_SpcNum
      if ( a_ampA(ID) /= 0.0d0 ) then 

        ! AMP �η�����Ϳ�����Ƥ�����
        !
        a_SwAmp(ID) = 1.0d0
        a_SwAnt(ID) = 0.0d0

      elseif ( a_antA(ID) /= 0.0d0 ) then 

        ! Antoine �η�������Ϳ�����Ƥ�����
        !
        a_SwAmp(ID) = 0.0d0
        a_SwAnt(ID) = 1.0d0

      else

        ! ���Τξ��
        !
        a_SwAmp(ID) = 0.0d0
        a_SwAnt(ID) = 0.0d0

      end if
    end do

    ! ����� (ʪ���ΰ�ξ�¤�Ϳ����) 
    a_kmax = nz

    ! ����� (ʪ���ΰ�β��¤�Ϳ����) 
    a_kmin = 1


    !!--------------------------------------------------------
    !! TempBZ, PressBZ �� allocated ������ͤ����ꤷ�Ƥ������, 
    !! �׻��򥵥ܤ뤿��ΰʲ��ν��֤�Ԥ�
    !!
    if ( allocated( xyz_TempBZ) ) then 

    !--------------------------------------------------------
    ! �׻���Ŭ���˥��ܤ뤿��ν��� (1) 
    !
    ! ˰�¾���������ʬ���������, �׻�����ɬ�פϤʤ��Ȥ����ɤ��Τ�����?
    ! �����׻��ξ��, ��ή�������ǤϿ�ξ������ϤۤȤ�ɥ������, ��ή
    ! �ˤ�äƱ����������ն�ޤǻ����夬��, ����ʤ����ȯ����.
    ! ��ŷŪ�˹��٤���ꤹ��Τ��񤷤����ʤΤ�, ���� (HeightUp) �����ꤵ�줿
    ! ���ˤ�, ���������η׻��ϹԤ�ʤ��Ȥ������Ȥ�. 

    ! HeightUp �����ꤵ��Ƥ�����˽�����Ԥ�. 
    ! ʪ����˰㤦�Ȥ������ȤϤ��ꤨ�ʤ��Τ�, ʪ�����Ф���롼�פϲ󤵤ʤ�. 
    !
    if ( HeightUp > 0.0d0 ) then 

      do k = kmin, kmax
        if ( z_Z(k) > HeightUp ) then 
          a_kmax = k
          exit
        end if
      end do
     
    end if

    !--------------------------------------------------------
    ! �׻���Ŭ���˥��ܤ뤿��ν��� (2) 
    ! 
    ! ˰�¾�������׻����뤿�������ź���β��¤����. 
    ! HeightDown �����ꤵ��Ƥ�����ˤϤ����ͥ�褷, 
    ! HeightDown ����ξ��ˤϰʲ��μ�³���ǲ������ꤹ��. 
    !
    ! * ˰�¾����������������Ǥΰ��Ϥ����Ϥ����ʤ��Ȥ������Ȥ���˷���.
    !   * �ɤ�ʪ�����Ф��Ƥ� Antoine �η�����Ϳ�����Ƥ��뤳�Ȥ�����. 
    ! 
    ! Fujitsu Fortran �Ǥ�, exp(logsvap) [logsvap > 700] �ǥ��顼���Ф�. 
    ! HeightUp �����Ǥ��äƤ�, logsvap > 700 �Ȥʤ���٤򲼸¤Ȥ���. 

    do ID = 1, ChemData_SpcNum
      ! �ŷ�ʪ�ξ��˷׻���Ԥ� (���Τξ��� a_antA = 0.0)
      !
      if ( a_antA(ID) /= 0.0d0 ) then 

        ! ��ľ�����Ͼ������롼�פ��.        
        do k = nz, 1, -1

          ! ���ܾ�β��٤��Ф���˰�¾������� log ��׻�
          !
          Temp  = xyz_TempBZ(1,1,k)
          Press = xyz_PressBZ(1,1,k)
          logsvap =                                    &
            &       (                                  &
            &           a_antA(ID)                     &
            &         - a_antB(ID)                     &
            &           / (a_antC(ID) + Temp - Temp0C) &
            &        ) * dlog(10.0d0)                  &
            &        + a_antU(ID)

!          write(*,*) ChemData_SpcSymbol(ID), k, dexp( logsvap ), PressSfc, Press

          ! logsvap > 700 �Ȥʤ��ź�����ݴ�
          ! 
          if ( logsvap > 700 ) then
            a_kmin(ID) = k
            exit

          ! ���¤����ꤵ�줿���
          !
          elseif ( z_Z(k) <= HeightDown ) then 
            a_kmin(ID) = k
            exit

          ! ˰�¾�������������
          !
!!          elseif( HeightDown < 0.0d0 .AND. logsvap >= dlog( PressSfc ) ) then 
          elseif( logsvap >= dlog( PressSfc ) ) then 
            a_kmin(ID) = k
            exit

          end if
        end do
        
      end if
    end do
    end if

    !--------------------------------------------------------
    ! �ͤγ�ǧ
    !
    call MessageNotify( "M", &
      & "ChemCalc_Init", "ReactHeatNH4SH = %f", d=(/ReactHeatNH4SH/) )
    id   = ChemData_OneSpcID( Name )  
    call MessageNotify( "M", &
      & "ChemCalc_Init", "NH4SH MolWt = %f", d=(/MolWt(id)/) )
    
    do k = 1, ChemData_SpcNum
      call MessageNotify( "M", "ChemCalc_Init", &
        & "%c : a_kmin= %d, a_kmax= %d, s_SwAmp= %f, s_SwAnt= %f", &
        & i=(/a_kmin(k), a_kmax(k)/), d=(/a_SwAMP(k), a_SwAnt(k)/), c1=trim(ChemData_SpcSymbol(k)))
    end do

  end subroutine ChemCalc_Init

!!!
!!! ˰�¾�����, ��Ǯ, etc. �δ��ܴؿ�. 
!!! ���ؼ�� ID �Ȳ��٤��Ф����ͤ��֤�
!!!  

!!!==========================================================================
  function CpRef(ID)
    !
    !������Ϳ����줿���ؼ���Ф���, ɸ����֤Ǥ�ñ�̼�����������갵��Ǯ��׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use chemdata,   only: ChemData_CpRef !ɸ����֤Ǥ�ñ�̼������������Ǯ

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: CpRef        !ɸ����֤Ǥ�ñ�̼������������Ǯ
    integer, intent(in) :: ID           !���ؼ�� ID

    
    !�ǡ����١�������������
    CpRef = ChemData_CpRef(ID)

  end function CpRef


!!!==========================================================================
  function CpPerMolRef(ID)
    !
    !������Ϳ����줿���ؼ���Ф���, ɸ����֤Ǥ�ñ�̥����������갵��Ǯ��׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP                   !���ٻ���
    use chemdata,   only: ChemData_CpPerMolRef !ɸ����֤Ǥ�ñ�̥�����������Ǯ

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: CpPerMolRef  !ɸ����֤Ǥ�ñ�̥�����������Ǯ
    integer, intent(in) :: ID           !���ؼ�� ID

    
    !�ǡ����١�������������
    CpPerMolRef = ChemData_CpPerMolRef(ID)

  end function CpPerMolRef


!!!==========================================================================
  function CvRef(ID)
    !
    !������Ϳ����줿���ؼ���Ф���, ɸ����֤Ǥ�ñ�̼�����������갵��Ǯ��׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use chemdata,   only: ChemData_CvRef !ɸ����֤Ǥ�ñ�̼������������Ǯ

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: CvRef       !ɸ����֤Ǥ�ñ�̼������������Ǯ
    integer, intent(in) :: ID          !���ؼ�� ID

    
    !�ǡ����١�������������
    CvRef = ChemData_CvRef(ID)

  end function CvRef


!!!==========================================================================
  function MolWt(ID)
    !
    !������Ϳ����줿���ؼ���Ф���, ʬ���̤�׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use chemdata,   only: ChemData_MolWt !ʬ����

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: MolWt         !ʬ����
    integer, intent(in) :: ID            !���ؼ�� ID

    
    !�ǡ����١�������������
    MolWt = ChemData_MolWt(ID)

  end function MolWt


!!!==========================================================================
  function GasR(ID)
    !
    !������Ϳ����줿���ؼ���Ф���, ���������׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use chemdata,   only: ChemData_GasR  !������� [J/K kg]

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: GasR          !ʬ����
    integer, intent(in) :: ID            !���ؼ�� ID
    
    
    !�ǡ����١�������������
    GasR = ChemData_GasR(ID)

  end function GasR


!!!
!!! ���� 3 �����δؿ���
!!!

!!!==========================================================================
  function xyz_SvapPress( ID, xyz_Temp )
    !
    != ������Ϳ����줿���ؼ�Ȳ��٤��Ф���, ˰�¾�������׻�. 
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use gridset,    only: nx, ny,      & !ʪ���ΰ���礭��
      &                   imin, imax,  & !����� X �����ξ�¡�����
      &                   jmin, jmax,  & !����� Y �����ξ�¡�����
      &                   kmin, kmax     !����� Z �����ξ�¡�����
    use constants,  only: PressSfc       !���������Ǥΰ���       [Pa]

    !���ۤη�����ػ�
    implicit none
    
    ! �������ѿ�  
    real(DP)            :: xyz_SvapPress(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !˰�¾�����
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)  
                                                          !����
    integer, intent(in) :: ID                             !���ؼ�� ID
  
    !�����ѿ�
    real(DP)            :: LogSvapPress
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! �����
    ! * ˰�¾������Ͻ�ʬ�礭���ͤˤ��Ƥ���.
    !
    xyz_SvapPress = PressSfc * 100.0d0

    ! ˰�¾������η׻�
    ! a_SwAmp, a_SwAnt ���Ѥ��뤳�Ȥ�, ���򤵤줿�׻���ˡ���Ѥ���.
    !
    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx
          
          ! ˰�¾������� log ��׻�
          !
          LogSvapPress =                                           &
            &      (                                               &
            &         a_ampA(ID) / xyz_Temp(i,j,k)                 &
            &       + a_ampB(ID)                                   &
            &       + a_ampC(ID) * dlog( xyz_Temp(i,j,k) )         &
            &       + a_ampD(ID) * xyz_Temp(i,j,k)                 &
            &       + a_ampE(ID) * ( xyz_temp(i,j,k) ** 2 )        &
            &       + dlog(1.0d-1)                                 &
            &      ) * a_SwAmp(ID)                                 &
            &    + (                                               &
            &       + (                                            &
            &          + a_antA(ID)                                &
            &          - a_antB(ID)                                &
            &            / (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) &
            &         ) * dlog(10.0d0)                             &
            &       + a_antU(ID)                                   &
            &      ) * a_SwAnt(ID)
          
          !˰�¾�������׻�
          !
          xyz_SvapPress(i,j,k) =  dexp( LogSvapPress )

        end do
      end do
    end do

  end function xyz_SvapPress  

!!!==========================================================================
  function xyz_LatentHeat(ID, xyz_Temp)
    !
    != ˰�¾����������Ǯ��׻�����. 
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use gridset,    only: nx, ny,      & !ʪ���ΰ���礭��
      &                   imin, imax,  & !����� X �����ξ�¡�����
      &                   jmin, jmax,  & !����� Y �����ξ�¡�����
      &                   kmin, kmax     !����� Z �����ξ�¡�����

    !���ۤη�����ػ�
    implicit none

    !�������ѿ�
    real(DP)            :: xyz_LatentHeat(imin:imax,jmin:jmax,kmin:kmax)
                                                            !��Ǯ[J/K kg]
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
                                                    !����[K]
    integer, intent(in) :: ID                       !���ؼ�� ID

    !�����ؿ�
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: GasRUniv = 8.314d0
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! �����
    !
    xyz_LatentHeat = 0.0d0

    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx

          ! ˰�¾������β�����ʬ
          ! a_SwAmp, a_SwAnt ���Ѥ��뤳�Ȥ�, ���򤵤줿�׻���ˡ���Ѥ���.
          !
          DLogSvapPressDTemp =                                             &
            &    (                                                         &
            &     - a_ampA(ID) / (xyz_Temp(i,j,k) ** 2.0d0)                &
            &     + a_ampC(ID) / xyz_Temp(i,j,k)                           &
            &     + a_ampD(ID)                                             &
            &     + a_ampE(ID) * 2.0d0 * xyz_Temp(i,j,k)                   &
            &    ) * a_SwAmp(ID)                                           &
            &  + (                                                         &
            &     + a_antB(ID) * dlog(10.0d0)                              &
            &       / ( (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) ** 2.0d0 ) &
            &    ) * a_SwAnt(ID)
          
          xyz_LatentHeat(i,j,k) =                                          &
            & DLogSvapPressDTemp * GasRUniv * (xyz_Temp(i,j,k) ** 2.0d0)   &
            &  / a_MolWt(ID)

        end do
      end do
    end do
    
  end function xyz_LatentHeat


!!!==========================================================================
  function SvapPress(ID, Temp)
    !
    != ������Ϳ����줿���ؼ�Ȳ��٤��Ф���, ˰�¾�������׻�. 
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���

    !���ۤη�����ػ�
    implicit none
    
    !�������ѿ�  
    real(DP)            :: SvapPress   !˰�¾�����
    real(DP),intent(in) :: Temp        !���� [K]
    integer, intent(in) :: ID          !���ؼ�� ID

    !�����ѿ�
    real(DP)            :: LogSvapPress
    real(DP), parameter :: Temp0C = 273.15d0

    ! ˰�¾������� log ��׻�
    ! �п����礭���ʤꤹ���ʤ��褦�ˤ���. 
    ! Fujitsu Fortran Compiler �Ǥ� 700 ����礭������ exp ����ȷٹ𤬽Ф�.
    !
    LogSvapPress =                               &
      & min(                                     &
      &      (                                   &
      &         a_ampA(ID) / Temp                &
      &       + a_ampB(ID)                       &
      &       + a_ampC(ID) * dlog( Temp )        &
      &       + a_ampD(ID) * Temp                &
      &       + a_ampE(ID) * ( Temp ** 2 )       &
      &       + dlog(1.0d-1)                     &
      &      ) * a_SwAmp(ID)                     &
      &    + (                                   &
      &        (                                 &
      &         + a_antA(ID)                     &
      &         - a_antB(ID)                     &
      &           / (a_antC(ID) + Temp - Temp0C) &
      &        ) * dlog(10.0d0)                  &
      &       + a_antU(ID)                       &
      &      ) * a_SwAnt(ID),                    &
      &   700.0d0                                &
      & )
          
    !˰�¾�������׻�
    !
    SvapPress =  dexp( LogSvapPress )

  end function SvapPress


!!!==========================================================================
  function LatentHeatPerMol(ID, Temp)
    !
    != ������Ϳ����줿���ؼ�Ȳ��٤��Ф���, ��Ǯ [J/K/mol] ��׻�
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���

    !���ۤη�����ػ�
    implicit none

    !�������ѿ�
    real(DP)            :: LatentHeatPerMol   !��Ǯ
    real(DP),intent(in) :: Temp               !����
    integer, intent(in) :: ID                 !���ؼ�̾
    
    !�����ѿ�
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: GasRUniv = 8.314d0  !���׵������
    real(DP),parameter  :: Temp0C = 273.15d0

    ! ˰�¾������β�����ʬ
    ! a_SwAmp, a_SwAnt ���Ѥ��뤳�Ȥ�, ���򤵤줿�׻���ˡ���Ѥ���.
    !
    DLogSvapPressDTemp =                                  &
      &    (                                              &
      &     - a_ampA(ID) / (Temp ** 2.0d0)                &
      &     + a_ampC(ID) / Temp                           &
      &     + a_ampD(ID)                                  &
      &     + a_ampE(ID) * 2.0d0 * Temp                   &
      &    ) * a_SwAmp(ID)                                &
      &  + (                                              &
      &     + a_antB(ID) * dlog(10.0d0)                   &
      &       / ( (a_antC(ID) + Temp - Temp0C) ** 2.0d0 ) &
      &    ) * a_SwAnt(ID)
          
    ! ��Ǯ�η׻�
    !
    LatentHeatPerMol =                                    &
      & DLogSvapPressDTemp * GasRUniv * (Temp ** 2.0d0)   

  end function LatentHeatPerMol

!!!-----------------------------------------------------------------------!!!
  function xyz_DQMixSatDPTemp(ID, MolWt, xyz_Temp, xyz_Exner)
    !
    !˰�¾������� �� ��ʬ��Ԥ�
    !�ºݤˤ�, dq/dp * dp/dT * dT/d�� ��¹�. (â�� p ��˰�¾�����)
    !
    ! * dq/dp =  Mv / (Md * p_all) 
    !   (q = p * Mv / (Md * p_all) )
    ! * dT/d��= \pi  (T = \pi \theta)
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use gridset,    only: nx, ny,      & !ʪ���ΰ���礭��
      &                   imin, imax,  & !����� X �����ξ�¡�����
      &                   jmin, jmax,  & !����� Y �����ξ�¡�����
      &                   kmin, kmax     !����� Z �����ξ�¡�����
    use constants,  only: PressBasis,  & !���̤�ɸ�వ��         [Pa]
      &                   CpDry,       & !������ʬ��ʿ���갵��Ǯ [J/K kg]
      &                   MolWtDry,    & !������ʬ��ʿ��ʬ����   [kg/mol]
      &                   GasRDry        !������ʬ�ε������     [J/K kg]

    !���ۤη�����ػ�
    implicit none 
    
    !�������ѿ�
    integer, intent(in) :: ID
    real(DP),intent(in) :: MolWt
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
                                            !����(���� + ���ܾ�)
    real(DP),intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                                            !�������ʡ��ؿ�(���� + ���ܾ�)
    real(DP)            :: xyz_DQMixSatDPTemp(imin:imax,jmin:jmax,kmin:kmax)
                           
    !�����ѿ�
    real(DP)            :: xyz_Press(imin:imax,jmin:jmax,kmin:kmax)
                                            !����(���� + ���ܾ�)
    real(DP)            :: DSvapPressDTemp
                                            !˰�¾������β�����ʬ [Pa/K]
    real(DP)            :: LogSvapPress
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! �����
    !
    xyz_DQMixSatDPTemp = 0.0d0
    xyz_Press = PressBasis * (xyz_Exner ** (CpDry / GasRDry))

    ! ˰�¾������β�����ʬ
    !
    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx
          ! ˰�¾������� log ��׻�
          !
          LogSvapPress =                                           &
            &      (                                               &
            &         a_ampA(ID) / xyz_Temp(i,j,k)                 &
            &       + a_ampB(ID)                                   &
            &       + a_ampC(ID) * dlog( xyz_Temp(i,j,k) )         &
            &       + a_ampD(ID) * xyz_Temp(i,j,k)                 &
            &       + a_ampE(ID) * ( xyz_temp(i,j,k) ** 2 )        &
            &       + dlog(1.0d-1)                                 &
            &      ) * a_SwAmp(ID)                                 &
            &    + (                                               &
            &       + (                                            &
            &          + a_antA(ID)                                &
            &          - a_antB(ID)                                &
            &            / (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) &
            &         ) * dlog(10.0d0)                             &
            &       + a_antU(ID)                                   &
            &      ) * a_SwAnt(ID)

          ! ˰�¾������β�����ʬ
          !
          DLogSvapPressDTemp =                                             &
            &    (                                                         &
            &     - a_ampA(ID) / (xyz_Temp(i,j,k) ** 2.0d0)                &
            &     + a_ampC(ID) / xyz_Temp(i,j,k)                           &
            &     + a_ampD(ID)                                             &
            &     + a_ampE(ID) * 2.0d0 * xyz_Temp(i,j,k)                   &
            &    ) * a_SwAmp(ID)                                           &
            &  + (                                                         &
            &     + a_antB(ID) * dlog(10.0d0)                              &
            &       / ( (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) ** 2.0d0 ) &
            &    ) * a_SwAnt(ID)
      
          DSvapPressDTemp = DLogSvapPressDTemp * dexp( LogSvapPress ) 

          xyz_DQMixSatDPTemp(i,j,k) =                               &
            &   MolWt / ( MolWtDry * xyz_Press(i,j,k) )             &
            &   * DSvapPressDTemp * xyz_Exner(i,j,k)   
          
        end do
      end do
    end do
    
  end function xyz_DQMixSatDPTemp


!!!-----------------------------------------------------------------------!!!
  function xyz_DelQMixNH4SH(xyz_TempAll, xyz_PressAll, xyz_PressDry, &
    &                       xyz_QMixNH3, xyz_QMixH2S, &
    &                       MolWtNH3, MolWtH2S)
    !
    ! NH4SH ����ȿ����ȼ��, NH4SH ��������(������)�����
    !

    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���
    use gridset,    only: imin, imax,  & ! ����� X �����ξ�¡�����
      &                   jmin, jmax,  & ! ����� Y �����ξ�¡�����
      &                   kmin, kmax     ! ����� Z �����ξ�¡�����
    use constants,  only: MolWtDry

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP),intent(in) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                         !����
    real(DP),intent(in) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
                                         !����
    real(DP),intent(in) :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax)
                                         !����
    real(DP),intent(in) :: xyz_QMixNH3(imin:imax,jmin:jmax,kmin:kmax)
                                         !NH3 �κ�����
    real(DP),intent(in) :: xyz_QMixH2S(imin:imax,jmin:jmax,kmin:kmax)
                                         !H2S �κ�����
    real(DP),intent(in) :: MolWtNH3      !NH3 ��ʬ����
    real(DP),intent(in) :: MolWtH2S      !H2S ��ʬ����

    real(DP) :: xyz_DelQMixNH4SH(imin:imax,jmin:jmax,kmin:kmax)
                                         !NH4SH �κ�����
    real(DP) :: xyz_EquivConst(imin:imax,jmin:jmax,kmin:kmax)
                                         !��ʿ�����
    real(DP) :: xyzf_PPress(imin:imax,jmin:jmax,kmin:kmax,2)
                                         !�������(ʬ��)
    real(DP) :: xyz_Solution(imin:imax,jmin:jmax,kmin:kmax)
                                         !�������(�������β�)

    !�����
!    xyz_DelQMixNH4SH = 0.0d0
    
    !�����˥���β�����Ǥ�ʬ��. 
    xyzf_PPress(:,:,:,1) = xyz_QMixNH3 * xyz_PressAll * MolWtDry / MolWtNH3 
    xyzf_PPress(:,:,:,2) = xyz_QMixH2S * xyz_PressAll * MolWtDry / MolWtH2S 

    !��ʿ�����
    xyz_EquivConst = 61.781d0 - 10834.0d0 / xyz_TempAll - dlog(1.0d2)

    !�����Ѳ������. 
    !  (P_NH3 - X) * (P_H2S - X) = exp(Kp)
    !  DelX^2 - (P_NH3 + P_H2S) * DelX + P_NH3 * P_H2S - exp( Kp ) = 0
    !  �Ȥ����������������ɬ�פ����뤬, (P_NH3 - X) > 0 ��
    !  (P_H2S - X) > 0 ������������ˤ�, ��θ����Τ���������������򤵤��.
    !
    xyz_Solution  =                                                       &
      & (                                                                 &
      &     sum(xyzf_PPress, 4)                                           &
      &   - dsqrt( (xyzf_PPress(:,:,:,1) - xyzf_PPress(:,:,:,2)) ** 2.0d0 &
      &            + 4.0d0 * dexp( min( 700.0d0, xyz_EquivConst ) ) )     &
      &  ) * 5.0d-1

    !�����̤����
    xyz_DelQMixNH4SH = xyz_Solution * ( MolWtNH3 + MolWtH2S ) &
      &                   / ( xyz_PressDry * MolWtDry )

  end function xyz_DelQMixNH4SH
  

!!!-----------------------------------------------------------------------!!!
  function DelMolFrNH4SH(TempAll, PressAll, MolFrNH3, MolFrH2S, Humidity)
    !
    ! NH4SH ����ȿ����ȼ�� H2S �� NH3 �Υ����θ���ʬ�����
    !
    
    !�⥸�塼��ƤӽФ�
    use dc_types,   only: DP             !���ٻ���

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP),intent(in) :: TempAll       !����
    real(DP),intent(in) :: PressAll      !����
    real(DP),intent(in) :: MolFrNH3      !NH3 �Υ����
    real(DP),intent(in) :: MolFrH2S      !H2S �Υ����
    real(DP),intent(in) :: Humidity      !˰����
    real(DP)            :: DelMolFrNH4SH !NH4SH ������ȼ��������Ѳ�
    real(DP)            :: EquivConst    !��ʿ�����
    real(DP)            :: PPress(2)     !�������(ʬ��)
    real(DP)            :: Solution      !�������(�������β�)

    !------------------------------------------------------------
    !NH4SH ��ʿ�վ��
    !------------------------------------------------------------
    !�����˥���β�����Ǥ�ʬ��
    PPress(1) = MolFrNH3 * PressAll
    PPress(2) = MolFrH2S * PressAll

    !��ʿ�����
    EquivConst = 61.781d0 - 10834.0d0 / TempAll - dlog(1.0d2) - 2.0d0 * dlog( Humidity )
    
    !�����Ѳ������������β�Ȥ��Ƶ���. 
    Solution = 5.0d-1 * (sum(PPress)                                        &
      &        - dsqrt( (PPress(1) - PPress(2))**2.0d0                      &
      &                    + 4.0d0 * dexp( min( 700.0d0, EquivConst ))) )
    
    !NH4SH ��������. 
    DelMolFrNH4SH = Solution / PressAll

  end function DelMolFrNH4SH

    
end module ChemCalc
