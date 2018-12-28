!= Program ArareInitData
!
! Authors::   SUGIYAMA Ko-ichiro (�����̰�ϯ), ODAKA Masatsugu (��������)
! Version::   $Id: arare_init-data.f90,v 1.30 2014/07/08 01:01:45 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

program ArareInitData
  !
  ! �����ϳإ�ǥ� deepconv/arare
  ! ����ͺ����ѥץ����
  !

  !----- �⥸�塼���ɤ߹��� ------

  !-----   �����, ʸ�������   ----
  use dc_types,       only : STRING, DP
  use dc_string,      only : StoA
  use gtool_history,  only : HistoryPut

  !-----   ��å���������   -----
  use dc_message,     only: MessageNotify, MessageSuppressMPI

  !  ���ޥ�ɥ饤��������
  use argset,        only : argset_init  

  !  MPI �ν����
  !  Initialize MPI wrapper
  !
  use mpi_wrapper,   only : MPIWrapperInit, MPIWrapperFinalize

  !-----    �����⥸�塼��   -----
  !  �����ϥե�����̾�����⥸�塼��
  use fileset,    only : fileset_init

  !  �ʻ��������⥸�塼�� 
  use gridset,    only : gridset_init, &
    &                    imin, imax, jmin, jmax, kmin, kmax, &
    &                    nx, ny, nz,  ncmax

  !  ���ܾ�����⥸�塼��
  use basicset,      only : basicset_init

  use axesset,       only : axesset_init

  use constants,     only : constants_init, molwtdry, &
    &                       Grav, TempSfc, PressSfc, &
    &                       pressbasis, gasrdry, cpdry, cvdry

  !  ����롼��������⥸�塼��
  use composition,      only: composition_init, molwtwet, SpcWetMolFr, &
    &                         CondNum, IdxCC, IdxCG, SpcWetID

  ! ���ط׻�
  use ChemCalc, only: ChemCalc_init, xyz_SvapPress

  !  �������Ŭ�ѥ⥸�塼��
  use setmargin

  use axesset, only:  z_Z,          &!�����顼�ʻ����Ǥι���
    &                 dz             !Z �����γʻ����ֳ�

  use namelist_util, only: NmlutilInit, namelist_filename

  use eccm,     only: ECCM_Dry,      &!
    &                 ECCM_Wet              

  use initialdata_disturb, only:                                  &
    &  initialdata_disturb_random,                                &
    &  initialdata_disturb_circleXZ,                              &
    &  initialdata_disturb_gaussXZ, initialdata_disturb_gaussXY,  &
    &  initialdata_disturb_gaussYZ, initialdata_disturb_gaussXYZ, &
    &  initialdata_disturb_cosXZ,   initialdata_disturb_cosXY,    &
    &  initialdata_disturb_cosYZ,   initialdata_disturb_cosXYZ,   &
    &  initialdata_disturb_coneXZ,  initialdata_disturb_coneXY,   &
    &  initialdata_disturb_coneYZ,                                &
    &  initialdata_disturb_tanh,    initialdata_disturb_tanh_sin, &
    &  initialdata_disturb_dryreg,  initialdata_disturb_moist,    &
    &  initialdata_disturb_square
  use initialdata_baker1998, only: initialdata_baker1998_basic
  use initialdata_Skamarock1994, only: initialdata_SK1994_basic, &
    &                                  initialdata_SK1994_disturbXZ, &
    &                                  initialdata_SK1994_disturbYZ
  use initialdata_yamasaki1983, only: initialdata_yamasaki1983_basic
  use initialdata_takemi2007, only: initialdata_takemi2007_init, &
    &  initialdata_takemi2007_basic, initialdata_takemi2007_wind
  use initialdata_toon2002, only: initialdata_toon2002_basic
  use initialdata_sounding, only: initialdata_sounding_init, &
    & initialdata_sounding_basic, initialdata_sounding_wind

  !-----    �����ϥ⥸�塼��   -----
  !  �ꥹ�����ȥե����������ϥ⥸�塼��
  use RestartFileIO,    only : ReStartFileIO_Init, ReStartFileIO_Finalize, rstat
  use Arare4InitFileIO, only : Arare4InitFileIO_Init, Arare4InitFileIO_Var_Get, Arare4InitFileIO_BZ_Get, &
    &                          Arare4InitFileIO_MMC_Var_Get, Arare4InitFileIO_MMC_BZ_Get

  !���ۤη�����ػ�
  implicit none

  !�����ѿ�
  real(DP), allocatable :: pyz_VelX(:,:,:)
  real(DP), allocatable :: xqz_VelY(:,:,:)
  real(DP), allocatable :: xyr_VelZ(:,:,:)
  real(DP), allocatable :: xyz_Exner(:,:,:)
  real(DP), allocatable :: xyz_PTemp(:,:,:)
  real(DP), allocatable :: xyz_Km(:,:,:)
  real(DP), allocatable :: xyz_Kh(:,:,:)
  real(DP), allocatable :: xyz_CDens(:,:,:)
  real(DP), allocatable :: xyzf_QMix(:,:,:,:)
  real(DP), allocatable :: xyz_DensBZ(:,:,:)
  real(DP), allocatable :: xyz_PressBZ(:,:,:)
  real(DP), allocatable :: xyz_ExnerBZ(:,:,:)
  real(DP), allocatable :: xyz_TempBZ(:,:,:)
  real(DP), allocatable :: xyz_PTempBZ(:,:,:)
  real(DP), allocatable :: xyz_VPTempBZ(:,:,:)
  real(DP), allocatable :: xyz_VelSoundBZ(:,:,:)
  real(DP), allocatable :: xyzf_QMixBZ(:,:,:,:)
  real(DP), allocatable :: xyz_EffMolWtBZ(:,:,:)
  real(DP), allocatable :: z_TempBZ(:)
  real(DP), allocatable :: z_PressBZ(:)
  real(DP), allocatable :: xyzf_MolFr(:,:,:,:)
  real(DP), allocatable :: zf_MolFr(:,:)
  real(DP), allocatable :: xyzf_QMixDivMolWt(:,:,:,:)
  real(DP), allocatable :: xyzf_HumBZ(:,:,:,:)
  real(DP)              :: Time
  integer               :: i, j, k, s
  logical, parameter    :: FlagInitData = .true.
  
  !�ѿ����
  ! Moist ��
  real(DP), save :: Humidity = 0.0d0 !���м���
  ! Gauss ��, cos ��, cone ��
  real(DP), save :: PTempMax = 0.0d0 !������
  real(DP), save :: ExnerMax = 0.0d0 !������
  real(DP), save :: QMixMax = 0.0d0  !������
  real(DP), save :: Xc = 0.0d0       !������濴����(X����)
  real(DP), save :: Yc = 0.0d0       !������濴����(Y����)
  real(DP), save :: Zc = 0.0d0       !������濴����(��ľ����)
  real(DP), save :: Xr = 0.0d0       !�����Ⱦ��(X����)
  real(DP), save :: Yr = 0.0d0       !�����Ⱦ��(Y����)
  real(DP), save :: Zr = 0.0d0       !�����Ⱦ��(��ľ����)
  ! tanh ��
  real(DP), save :: PTempMean= 0.0d0 !���������濴����
  real(DP), save :: VelMean  = 0.0d0 !���������濴����
  real(DP), save :: PTempDel = 0.0d0 !�������Ǥ��Ѳ���
  ! Therma-Random ��
  real(DP), save :: Zpos = 0.0d0     !����� Z ��ɸ [m] 
  ! Square ��
  real(DP), save :: XposMin = 0.0d0    !������ X ��ɸ [m] 
  real(DP), save :: YposMin = 0.0d0    !������ Y ��ɸ [m] 
  real(DP), save :: ZposMin = 0.0d0    !������ Z ��ɸ [m] 
  real(DP), save :: XposMax = 0.0d0    !������ X ��ɸ [m] 
  real(DP), save :: YposMax = 0.0d0    !������ Y ��ɸ [m] 
  real(DP), save :: ZposMax = 0.0d0    !������ Z ��ɸ [m] 
  ! WindConst ��
  real(DP), save :: VelX0 = 0.0d0   !����˰����®�٤�Ϳ������
  real(DP), save :: VelY0 = 0.0d0   !����˰����®�٤�Ϳ������
  real(DP), save :: VelZ0 = 0.0d0   !����˰����®�٤�Ϳ������

  real(DP), save :: TempTr  = 10.0d0    !��ή�����̤β��� [k]
  real(DP), save :: Dheight = 10.0d3    !�Ťߴؿ��Υѥ�᡼�� [m]
!  real(DP), save :: HeightTr= 10000.0d3 !��ή�����̤ι��� [m]

  integer               :: IDBasic                = 0
  integer, parameter    :: IDBasicArare4          = 1
  integer, parameter    :: IDBasicArare4mmc       = 2
  integer, parameter    :: IDBasicYamasaki1983    = 3
  integer, parameter    :: IDBasicTakemi2007      = 4
  integer, parameter    :: IDBasicIsoThermal      = 5
  integer, parameter    :: IDBasicDry             = 6
  integer, parameter    :: IDBasicMoist           = 7
  integer, parameter    :: IDBasicToon2002        = 8
  integer, parameter    :: IDBasicSounding        = 9
  integer, parameter    :: IDBasicBaker1998       = 10
  integer, parameter    :: IDBasicSK1994   = 11
  integer               :: IDDisturbPTemp         = 0
  integer, parameter    :: IDDisturbPTempGaussXY  = 1
  integer, parameter    :: IDDisturbPTempGaussXZ  = 2
  integer, parameter    :: IDDisturbPTempGaussYZ  = 3
  integer, parameter    :: IDDisturbPTempGaussXYZ = 4
  integer, parameter    :: IDDisturbPTempRandom   = 5
  integer, parameter    :: IDDisturbPTempSquare= 6
  integer, parameter    :: IDDisturbPTempCosXY    = 7
  integer, parameter    :: IDDisturbPTempCosXZ    = 8
  integer, parameter    :: IDDisturbPTempCosYZ    = 9
  integer, parameter    :: IDDisturbPTempCosXYZ   = 10
  integer, parameter    :: IDDisturbPTempConeXY   = 11
  integer, parameter    :: IDDisturbPTempConeXZ   = 12
  integer, parameter    :: IDDisturbPTempConeYZ   = 13
  integer, parameter    :: IDDisturbPTempTanh     = 14
  integer, parameter    :: IDDisturbPTempSK1994XZ = 15
  integer, parameter    :: IDDisturbPTempSK1994YZ = 16
  integer, parameter    :: IDDisturbPTempCircleXZ = 17
  integer               :: IDDisturbExner         = 0
  integer, parameter    :: IDDisturbExnerGaussXY  = 1
  integer, parameter    :: IDDisturbExnerGaussXZ  = 2
  integer, parameter    :: IDDisturbExnerGaussYZ  = 3
  integer, parameter    :: IDDisturbExnerGaussXYZ = 4
  integer, parameter    :: IDDisturbExnerCosXY    = 5
  integer, parameter    :: IDDisturbExnerCosXZ    = 6
  integer, parameter    :: IDDisturbExnerCosYZ    = 7
  integer, parameter    :: IDDisturbExnerCosXYZ   = 8
  integer, parameter    :: IDDisturbExnerConeXY   = 9
  integer, parameter    :: IDDisturbExnerConeXZ   = 10
  integer, parameter    :: IDDisturbExnerConeYZ   = 11
  integer               :: IDDisturbQMix          = 0
  integer, parameter    :: IDDisturbQMixGaussXY   = 1
  integer, parameter    :: IDDisturbQMixGaussXZ   = 2
  integer, parameter    :: IDDisturbQMixGaussYZ   = 3
  integer, parameter    :: IDDisturbQMixGaussXYZ  = 4
  integer, parameter    :: IDDisturbQMixDryreg    = 5
  integer, parameter    :: IDDisturbQMixMoist     = 6
  integer, parameter    :: IDDisturbQMixCosXY     = 7
  integer, parameter    :: IDDisturbQMixCosXZ     = 8
  integer, parameter    :: IDDisturbQMixCosYZ     = 9
  integer, parameter    :: IDDisturbQMixCosXYZ    = 10
  integer, parameter    :: IDDisturbQMixConeXY    = 11
  integer, parameter    :: IDDisturbQMixConeXZ    = 12
  integer, parameter    :: IDDisturbQMixConeYZ    = 13
  integer, parameter    :: IDDisturbQMixSquare = 14
  integer               :: IDDisturbWind          = 0
  integer, parameter    :: IDDisturbWindTakemi2007= 1
  integer, parameter    :: IDDisturbWindSounding  = 2
  integer, parameter    :: IDDisturbWindZonal     = 3
  integer, parameter    :: IDDisturbWindTanh      = 4
  integer               :: IDDisturb              = 0
  integer, parameter    :: IDDisturbArare4        = 1
  integer, parameter    :: IDDisturbArare4mmc     = 2
  

  !------------------------------------------
  ! �������³�� ; Initialize procedure 
  !
  call MainInit

  Time = 0.0d0
  call MessageNotify( "M", "main", "Making Initial data...." )

  !---------------------------------------------------------------  
  ! ���ܾ���������. 
  !
  select case(IDBasic)
    
  case (IDBasicArare4)
    ! deepconv/arare4 �Υҥ��ȥ꡼�ե����뤫�����ͤ��ɤ߹���
    !
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_BZ_Get( &
      & xyz_PressBZ,   & ! (out)
      & xyz_ExnerBZ,   & ! (out)
      & xyz_TempBZ,    & ! (out)
      & xyz_PTempBZ,   & ! (out)
      & xyz_DensBZ,    & ! (out)
      & xyz_VelSoundBZ,& ! (out)
      & xyzf_QMixBZ,   & ! (out)
      & xyz_EffMolWtBZ & ! (out)
      & )

  case (IDBasicArare4mmc)
    ! deepconv/arare4 �����ǤΥҥ��ȥ꡼�ե����뤫�����ͤ��ɤ߹���
    !    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )

    call Arare4Initfileio_MMC_BZ_Get( &
      & xyz_PressBZ,   & ! (out)
      & xyz_ExnerBZ,   & ! (out)
      & xyz_TempBZ,    & ! (out)
      & xyz_PTempBZ,   & ! (out)
      & xyz_DensBZ,    & ! (out)
      & xyz_VelSoundBZ,& ! (out)
      & xyzf_QMixBZ,   & ! (out)
      & xyz_EffMolWtBZ & ! (out)
      & )
    
  case(IDBasicIsoThermal) 
    ! ����/���Ϲ�θ���ʤ����
    !  
    call MessageNotify( "M", "main", "Making Initial data (basic) named IsoThermal..." )
    
    z_TempBZ = TempSfc
    z_PressBZ = PressSfc
    zf_MolFr = 0.0d0
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            

  case(IDBasicDry) 
    ! ������ǮŪ�ʽ����
    !
    call MessageNotify( "M", "main", "Making Initial data (basic) named DRY..." )
    call eccm_dry( SpcWetMolFr(1:ncmax), Humidity, z_TempBZ, z_PressBZ, zf_MolFr )
    if (minval(z_TempBZ(1:nz)) < TempTr) then 
!      call initialdata_basic_strat( z_TempBZ, z_PressBZ ) !(inout)
      call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)
    end if
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            
    
! �ǥХå��Ѥ�Ǥ��ʤ�
!   
!  case(IDBasicMoist) 
!    ! ������ǮŪ�ʽ����
!    !
!    call MessageNotify( "M", "main", "Making Initial data (basic) named MOIST..." )
!    call eccm_wet( SpcWetMolFr(1:ncmax), Humidity, z_TempBZ, z_PressBZ, zf_MolFr )
!    call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)
!
!    ! �Ĥ�δ��ܾ���ͤ����
!    call DetermineBZ            

  case (IDBasicSounding)
    ! ������ǥ��󥰥ե����뤫���ͤ��ɤ߹���
    !
    call MessageNotify( "M", "main", "Making Initial data using sounding files ...." )
    call initialdata_sounding_basic( z_TempBZ, z_PressBZ )  ! (out)
    zf_MolFr = 0.0d0
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            

  case(IDBasicToon2002) 
    ! ������
    !
    call MessageNotify( "M", "main", "Making Initial data (basic) named Toon et al. 2002..." )
    call initialdata_toon2002_basic( z_TempBZ, z_PressBZ ) !(out)
!    call initialdata_basic_strat_v2( z_TempBZ, z_PressBZ ) !(inout)

    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            
    
  case(IDBasicTakemi2007)
    ! Takemi(2007)�δ��ܾ����Ѥ�����
    !
    call MessageNotify( "M", "main", "Making Initial data named Takemi2007..." )
    call initialdata_takemi2007_basic( z_TempBZ, z_PressBZ, zf_MolFr )

    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            

  case(IDBasicYamasaki1983) 
    ! Yamasaki(1983)�β��٤����м��٤δ�¬�ͤ���Ѥ����� 
    !
    call MessageNotify( "M", "main", "Making Initial data named Yamasaki1983..." )
    call initialdata_yamasaki1983_basic( z_TempBZ, z_PressBZ, zf_MolFr )
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            

  case(IDBasicBaker1998) 
    ! Baker and Schubert (1998) �β��٤����м��٤δ�¬�ͤ���Ѥ����� 
    !
    call MessageNotify( "M", "main", "Making Initial data named Baker and Shubert 1998..." )
    call initialdata_baker1998_basic( z_TempBZ, z_PressBZ ) !(out)
    zf_MolFr = 0.0d0
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            

  case(IDBasicSK1994) 
    ! ����:  Skamarock & Klemp (1994) 
    !
    call MessageNotify( "M", "main", "Making Initial data named Skamarock & Klemp (1994) ..." )
    call initialdata_SK1994_basic( z_TempBZ, z_PressBZ ) !(out)
    zf_MolFr = 0.0d0
    
    ! �Ĥ�δ��ܾ���ͤ����
    call DetermineBZ            
    
  end select

  !---------------------------------------------------------
  ! �������ͤ����
  ! 

  ! ����: �������ͤ����
  !    
  select case(IDDisturbPTemp)

  case(IDDisturbPTempGaussXY)
    ! ����: �����������ʬ�ۤ�Ϳ���� (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussXZ)
    ! ����: �����������ʬ�ۤ�Ϳ���� (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussYZ)
    ! ����: �����������ʬ�ۤ�Ϳ���� (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempGaussXYZ)
    ! ����: �����������ʬ�ۤ�Ϳ���� (XYZ)
    !                   
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(PTempMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempRandom)
    ! ����: ���ꤵ�줿���٤˥������ʬ�ۤ�Ϳ����. 
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named random..." )
    call initialdata_disturb_random(PTempMax, Zpos, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSquare)
    ! ����: ����ʾ����Ϳ����
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Square..." )
    call initialdata_disturb_square(PTempMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXY)
    ! ����: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXZ)
    ! ����: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosYZ)
    ! ����: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCosXYZ)
    ! ����: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XYZ)
    !                   
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(PTempMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )
    
  case(IDDisturbPTempConeXY)
    ! ����: �߿��ʬ�ۤ�Ϳ���� (XY)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(PTempMax, Xc, Xr, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempConeXZ)
    ! ����: �߿��ʬ�ۤ�Ϳ���� (XZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(PTempMax, Xc, Xr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempConeYZ)
    ! ����: �߿��ʬ�ۤ�Ϳ���� (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(PTempMax, Yc, Yr, Zc, Zr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempTanh)
    ! ����: tanh ��ʬ�ۤ�Ϳ���� (YZ)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named tanh..." )
    call initialdata_disturb_tanh_sin(PTempMean, PTempDel, Zc, Zr, xyz_PTemp, xyz_PTempBZ)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSK1994XZ)
    ! ����:  Skamarock & Klemp (1994)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named SK1994..." )
    call initialdata_SK1994_disturbXZ(PTempMax, Xc, Xr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempSK1994YZ)
    ! ����:  Skamarock & Klemp (1994)
    !                
    call MessageNotify( "M", "main", "Making Initial data (disturb) named SK1994..." )
    call initialdata_SK1994_disturbYZ(PTempMax, Yc, Yr, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  case(IDDisturbPTempCircleXZ)
    ! ����: ���ΰ��Ʊ�����̾��� (XZ)
    !          
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CircleXZ..." )
    call initialdata_disturb_circleXZ(PTempMax, Xc, Xr, Zc, xyz_PTemp)
    call SetMargin_xyz( xyz_PTemp )

  end select

  
  ! �������ʡ��ؿ�: �������ͤ����    
  !                
  select case(IDDisturbExner)

  case(IDDisturbExnerGaussXY)
    ! �������ʡ��ؿ�: �����������ʬ�ۤ�Ϳ���� (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussXZ)
    ! �������ʡ��ؿ�: �����������ʬ�ۤ�Ϳ���� (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussYZ)
    ! �������ʡ��ؿ�: �����������ʬ�ۤ�Ϳ���� (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerGaussXYZ)      
    ! �������ʡ��ؿ�: �����������ʬ�ۤ�Ϳ���� (XYZ)
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(ExnerMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXY)
    ! �������ʡ��ؿ�: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXZ)
    ! �������ʡ��ؿ�: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosYZ)
    ! �������ʡ��ؿ�: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerCosXYZ)      
    ! �������ʡ��ؿ�: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XYZ)
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(ExnerMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeXY)
    ! �������ʡ��ؿ�: �߿��ʬ�ۤ�Ϳ���� (XY)
    !                         
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(ExnerMax, Xc, Xr, Yc, Yr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeXZ)
    ! �������ʡ��ؿ�: �߿��ʬ�ۤ�Ϳ���� (XZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(ExnerMax, Xc, Xr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  case(IDDisturbExnerConeYZ)
    ! �������ʡ��ؿ�: �߿��ʬ�ۤ�Ϳ���� (YZ)
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(ExnerMax, Yc, Yr, Zc, Zr, xyz_Exner)
    call SetMargin_xyz( xyz_Exner )

  end select


  ! ������: �������ͤ����
  !                    
  select case(IDDisturbQMix)

  case(IDDisturbQMixGaussXY)
    ! ������: �����������ʬ�ۤ�Ϳ���� (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXY..." )
    call initialdata_disturb_gaussXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixGaussXZ)
    ! ������: �����������ʬ�ۤ�Ϳ���� (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXZ..." )
    call initialdata_disturb_gaussXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixGaussYZ)
    ! ������: �����������ʬ�ۤ�Ϳ���� (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussYZ..." )
    call initialdata_disturb_gaussYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixGaussXYZ)
    ! ������: �����������ʬ�ۤ�Ϳ���� (XYZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named GaussXYZ..." )
    call initialdata_disturb_gaussXYZ(QMixMax, Xc, Xr, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixDryreg)
    ! ������: �����ΰ����
    !                            
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Dryreg..." )
    call initialdata_disturb_Dryreg(XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixMoist)
    ! ������: ������ǮŪ��ʬ�ۤ�Ϳ���� 
    !                                  
    call MessageNotify( "M", "main", "Making Initial data (disturb) named MOIST..." )
    call initialdata_disturb_moist(Humidity, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosXY)
    ! ������: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXY..." )
    call initialdata_disturb_cosXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixCosXZ)
    ! ������: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXZ..." )
    call initialdata_disturb_cosXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosYZ)
    ! ������: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosYZ..." )
    call initialdata_disturb_cosYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixCosXYZ)
    ! ������: A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0��ʬ�ۤ�Ϳ���� (XYZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named CosXYZ..." )
    call initialdata_disturb_cosXYZ(QMixMax, Xc, Xr, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixConeXY)
    ! ������: �߿��ʬ�ۤ�Ϳ���� (XY)
    !                               
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXY..." )
    call initialdata_disturb_coneXY(QMixMax, Xc, Xr, Yc, Yr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )
    
  case(IDDisturbQMixConeXZ)
    ! ������: �߿��ʬ�ۤ�Ϳ���� (XZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeXZ..." )
    call initialdata_disturb_coneXZ(QMixMax, Xc, Xr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixConeYZ)
    ! ������: �߿��ʬ�ۤ�Ϳ���� (YZ)
    !                                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named ConeYZ..." )
    call initialdata_disturb_coneYZ(QMixMax, Yc, Yr, Zc, Zr, xyzf_QMix(:,:,:,1))
    call SetMargin_xyzf( xyzf_QMix )

  case(IDDisturbQMixSquare)
    ! ������: ����ʾ����Ϳ����
    !                      
    call MessageNotify( "M", "main", "Making Initial data (disturb) named Square..." )
    call initialdata_disturb_square(QMixMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, xyzf_QMix)
    call SetMargin_xyzf( xyzf_QMix )

  end select

  ! ®�پ���ͤ����
  !                    
  select case(IDDisturbWind)  
  case (IDDisturbWindTakemi2007)  
    ! Takemi (2007) ����®ʬ��
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named takemi2007..." )
    call initialdata_takemi2007_wind( pyz_VelX ) 
    call SetMargin_pyz( pyz_VelX )   

  case (IDDisturbWindSounding)  
    ! 
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named Sounding..." )
    call initialdata_sounding_wind( pyz_VelX, xqz_VelY ) 
    call SetMargin_pyz( pyz_VelX )   
    call SetMargin_xqz( xqz_VelY )   

  case (IDDisturbWindZonal)  
    ! 
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named Zonal..." )

    ! �ͤ�����
    pyz_VelX = VelX0
    xqz_VelY = VelY0
    xyr_VelZ = VelZ0

    call SetMargin_pyz( pyz_VelX )
    call SetMargin_xqz( xqz_VelY )
    call SetMargin_xyr( xyr_VelZ )

  case (IDDisturbWindTanh)  
    ! tanh ������®
    !
    call MessageNotify( "M", "main", "Making Initial wind data (disturb) named tanh..." )
    call initialdata_disturb_tanh(VelMean, VelMean, Zc, Zr, pyz_VelX)
    call SetMargin_pyz( pyz_VelX )   

  end select

  ! ��礷�ƾ�����Ϳ������
  !
  select case(IDDisturb)  
  case (IDDisturbArare4)    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_Var_Get( &
      & pyz_VelX,     & ! (out)
      & xqz_VelY,     & ! (out)
      & xyr_VelZ,     & ! (out)
      & xyz_PTemp,    & ! (out)
      & xyz_Exner,    & ! (out)
      & xyzf_QMix,    & ! (out)
      & xyz_Km,       & ! (out)
      & xyz_Kh,       & ! (out)
      & xyz_CDens   )   ! (out)

  case (IDDisturbArare4mmc)    
    call MessageNotify( "M", "main", "Making Initial data using arare4 files ...." )
    call Arare4Initfileio_MMC_Var_Get( &
      & pyz_VelX,     & ! (out)
      & xqz_VelY,     & ! (out)
      & xyr_VelZ,     & ! (out)
      & xyz_PTemp,    & ! (out)
      & xyz_Exner,    & ! (out)
      & xyzf_QMix,    & ! (out)
      & xyz_Km,       & ! (out)
      & xyz_Kh,       & ! (out)
      & xyz_CDens   )   ! (out)
  end select

  !------------------------------------------------
  ! �ե��������
  !
  call MessageNotify( "M", "main", "Output variables into netCDF file..." )

  ! �ꥹ�����ȥե��������. ���ܾ�Ⱦ��������. 
  !
  call ReStartFileIO_Init( FlagInitData )

  call HistoryPut( 't',     Time,      rstat)
  call HistoryPut( 'VelX',  pyz_VelX,  rstat)
  call HistoryPut( 'VelY',  xqz_VelY,  rstat)
  call HistoryPut( 'VelZ',  xyr_VelZ,  rstat)
  call HistoryPut( 'Exner', xyz_Exner, rstat)
  call HistoryPut( 'PTemp', xyz_PTemp, rstat)
  call HistoryPut( 'Km',    xyz_Km,    rstat)
  call HistoryPut( 'Kh',    xyz_Kh,    rstat)
  call HistoryPut( 'QMix',  xyzf_QMix, rstat)    
  call HistoryPut( 'CDens', xyz_CDens, rstat)    

  ! ���ܾ�Υե��������
  !
  xyz_VPTempBZ = xyz_PTempBZ / xyz_EffMolWtBZ  ! �ե��������ϤΤ���

  call HistoryPut( 'DensBZ',     xyz_DensBZ    , rstat)
  call HistoryPut( 'ExnerBZ',    xyz_ExnerBZ   , rstat)
  call HistoryPut( 'PTempBZ',    xyz_PTempBZ   , rstat)
  call HistoryPut( 'VPTempBZ',   xyz_VPTempBZ  , rstat)
  call HistoryPut( 'VelSoundBZ', xyz_VelSoundBZ, rstat)
  call HistoryPut( 'TempBZ',     xyz_TempBZ    , rstat)
  call HistoryPut( 'PressBZ',    xyz_PressBZ   , rstat)
  call HistoryPut( 'QMixBZ',     xyzf_QMixBZ   , rstat)
  call HistoryPut( 'EffMolWtBZ', xyz_EffMolWtBZ, rstat)
  call HistoryPut( 'HumBZ',      xyzf_HumBZ    , rstat)
  
  call ReStartFileIO_Finalize

  !----------------------------------------------------
  ! MPI END
  !
  call MPIWrapperFinalize
  
contains

  subroutine ArareAlloc
    !
    !������Ȥ���, ����������, �ͤȤ��ƥ������������.
    !

    !���ۤη�����ػ�
    implicit none

    !����������
    allocate( pyz_VelX(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xqz_VelY(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_VelZ(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Exner(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_PTemp(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Km(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_Kh(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_CDens(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, ncmax) )

    allocate( xyz_DensBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_PressBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_ExnerBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_TempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_PTempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_VPTempBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyz_VelSoundBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( xyzf_QMixBZ(imin:imax,jmin:jmax,kmin:kmax,ncmax))
    allocate( xyz_EffMolWtBZ(imin:imax,jmin:jmax,kmin:kmax))
    allocate( z_TempBZ(kmin:kmax))
    allocate( z_PressBZ(kmin:kmax))
    allocate( xyzf_MolFr(imin:imax, jmin:jmax, kmin:kmax, ncmax))
    allocate( zf_MolFr(kmin:kmax, ncmax))
    allocate( xyzf_QMixDivMolWt(imin:imax,jmin:jmax,kmin:kmax, ncmax))
    allocate( xyzf_HumBZ(imin:imax,jmin:jmax,kmin:kmax, ncmax))

    pyz_VelX = 0.0d0
    xqz_VelY = 0.0d0
    xyr_VelZ = 0.0d0
    xyz_PTemp = 0.0d0
    xyz_Exner = 0.0d0
    xyz_Km = 0.0d0
    xyz_Kh = 0.0d0
    xyz_CDens = 0.0d0
    xyzf_QMix = 0.0d0

    xyz_DensBZ = 0.0d0
    xyz_PressBZ = 0.0d0
    xyz_ExnerBZ = 0.0d0
    xyz_TempBZ = 0.0d0
    xyz_PTempBZ = 0.0d0
    xyz_VPTempBZ = 0.0d0
    xyz_VelSoundBZ = 0.0d0
    xyzf_QMixBZ = 0.0d0
    xyz_EffMolWtBZ  = 0.0d0
    z_TempBZ  = 0.0d0
    z_PressBZ = 0.0d0
    xyzf_MolFr = 0.0d0
    zf_MolFr = 0.0d0 
    xyzf_QMixDivMolWt = 0.0d0
    xyzf_HumBZ = 0.0d0

  end subroutine ArareAlloc
    

  subroutine DetermineBZ

    !---------------------------------------------------------------
    ! ���١�����
    !

    ! 3 ��������˳�Ǽ
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyz_TempBZ(i,j,k)  = z_TempBZ(k)  
          xyz_PressBZ(i,j,k) = z_PressBZ(k)
        end do
      end do
    end do
    
    !�Τ�����ͤ����
    !
    call SetMargin_xyz( xyz_TempBZ )
    call SetMargin_xyz( xyz_PressBZ)

    !---------------------------------------------------------------
    ! ������
    !
    !��ʿ�����ˤϰ���
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyzf_MolFr(i,j,k,:) = zf_MolFr(k,:)
        end do
      end do
    end do

    !�Τ�����ͤ����
    ! 
    call SetMargin_xyzf( xyzf_MolFr )

    !����Υ����򺮹�����Ѵ�
    do s = 1, ncmax
      xyzf_QMixBZ(:,:,:,s) = xyzf_MolFr(:,:,:,s) * MolWtWet(s) / MolWtDry
    end do
    
    !  !�ͤ��������ʤꤹ���ʤ��褦�˺����ͤ�Ϳ����
    !  where (xyzf_QMixBZ <= 1.0d-20 )
    !    xyzf_QMixBZ = 1.0d-20
    !  end where
    
    !�Τ�����ͤ����
    !
    call SetMargin_xyzf( xyzf_QMixBZ )
    
    !---------------------------------------------------------------
    ! ʬ���̤θ���
    !
    do s = 1, ncmax
      xyzf_QMixDivMolWt(:,:,:,s) = xyzf_QMixBZ(:,:,:,s) / MolWtWet(s)
    end do
    
    xyz_EffMolWtBZ = &
      & (1.0d0 + sum(xyzf_QMixBZ,4) ) &
      & / ( MolWtDry * ((1.0d0 / MolWtDry) + sum(xyzf_QMixDivMolWt,4)) )
    
    !�Τ�����ͤ����  
    !
    call SetMargin_xyz( xyz_EffMolWtBZ )

    !---------------------------------------------------------------    
    ! ����
    !
    xyz_PTempBZ = &
      & xyz_TempBZ * (PressBasis / xyz_PressBZ) ** (GasRDry / CpDry) 

    !�Τ�����ͤ����  
    !
    call SetMargin_xyz( xyz_PTempBZ )

    
    !---------------------------------------------------------------    
    ! �������ʡ��ؿ�
    !
    xyz_ExnerBZ = xyz_TempBZ / xyz_PTempBZ    
    
    !�Τ�����ͤ����
    !
    call SetMargin_xyz( xyz_ExnerBZ )

    !---------------------------------------------------------------    
    ! ̩��
    !    
    xyz_DensBZ = &
      & PressBasis * (xyz_ExnerBZ ** (CvDry / GasRDry)) &
      &  / (GasRDry * xyz_PTempBZ / xyz_EffMolWtBZ)
    
    !�Τ�����ͤ����
    !
    call SetMargin_xyz( xyz_DensBZ )
    
    !---------------------------------------------------------------    
    ! ��®
    !
    xyz_VelSoundBZ = &
      & sqrt ( &
      &   CpDry * GasRDry * xyz_ExnerBZ * xyz_PTempBZ &
      &   / (CvDry * xyz_EffMolWtBZ) &
      & )
    
    !�Τ�����ͤ����
    !
    call SetMargin_xyz( xyz_VelSoundBZ )

    !---------------------------------------------------------------    
    ! ����
    !
    if (CondNum >= 1) then 
      do s = 1, CondNum
        xyzf_HumBZ(:,:,:,IdxCG(s)) = &
          & xyz_PressBZ * xyzf_MolFr(:,:,:,IdxCG(s)) &
          & / xyz_SvapPress(SpcWetID(IdxCC(s)), xyz_TempBZ) &
          & * 100.0d0
      end do
    end if

    !----------------------------------------------------------
    ! BasicSet �⥸�塼����ͤ�����
    !   
    call BasicSet_Init(                           &
      & xyz_PressBZ, xyz_ExnerBZ, xyz_TempBZ,     &
      & xyz_PTempBZ, xyz_DensBZ,  xyz_VelSoundBZ, &
      & xyzf_QMixBZ, xyz_EffMolWtBZ )
    
  end subroutine DetermineBZ


  subroutine Initialdata_init

    use dc_iounit,     only : FileOpen    

    implicit none

    integer                       :: unit     !�����ֹ�

    character(STRING) :: FlagBasic        = ""
    character(STRING) :: FlagDisturbPTemp = ""
    character(STRING) :: FlagDisturbExner = ""
    character(STRING) :: FlagDisturbQMix  = ""
    character(STRING) :: FlagDisturbWind  = ""
    character(STRING) :: FlagDisturb      = ""

    NAMELIST /initialdata_nml/ &
      & FlagBasic,        & 
      & FlagDisturb,      &
      & FlagDisturbPTemp, &    
      & FlagDisturbExner, &    
      & FlagDisturbQMix,  &    
      & FlagDisturbWind  

    NAMELIST /initialdata_basic_nml/ &
!      & Humidity, TempTr, DHeight, HeightTr
      & Humidity, TempTr, DHeight

    NAMELIST /initialdata_disturb_gauss_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_cos_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_cone_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Yr, Zc, Zr

    NAMELIST /initialdata_disturb_tanh_nml/ &
      & PTempMean, PTempDel, VelMean, Zc, Zr

    NAMELIST /initialdata_disturb_random_nml/ &
!      & PTempMax, ExnerMax, QMixMax, Zpos
      & PTempMax, Zpos

    NAMELIST /initialdata_disturb_square_nml/ &
      & QMixMax, PTempMax, XposMin, YposMin, ZposMin, XposMax, YposMax, ZposMax

!    NAMELIST /initialdata_disturb_dryreg_nml/ &
!      & XposMin, YposMin, ZposMin, XposMax, YposMax, ZposMax

    NAMELIST /initialdata_disturb_zonal_nml/ &
      & VelX0, VelY0, VelZ0

    NAMELIST /initialdata_disturb_moist_nml/ &
      & Humidity

    NAMELIST /initialdata_disturb_SK1994_nml/ &
      & PTempMax, Xc, Xr, Yc, Yr

    NAMELIST /initialdata_disturb_circle_nml/ &
      & PTempMax, ExnerMax, QMixMax, Xc, Xr, Yc, Zc

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_nml)
    close(unit)
    
    if (FlagBasic == "Arare4") then 
      IDBasic = IDBasicArare4
    elseif(FlagBasic == "Arare4mmc") then 
      IDBasic = IDBasicArare4mmc
    elseif(FlagBasic == "Yamasaki1983") then 
      IDBasic = IDBasicYamasaki1983
    elseif(FlagBasic == "Takemi2007") then 
      IDBasic = IDBasicTakemi2007
    elseif(FlagBasic == "IsoThermal") then 
      IDBasic = IDBasicIsoThermal
    elseif(FlagBasic == "Dry"       ) then 
      IDBasic = IDBasicDry
    elseif(FlagBasic == "Moist"     ) then 
      IDBasic = IDBasicMoist
    elseif(FlagBasic == "Toon2002"  ) then 
      IDBasic = IDBasicToon2002
    elseif(FlagBasic == "Sounding") then 
      IDBasic = IDBasicSounding
    elseif(FlagBasic == "Baker1998") then 
      IDBasic = IDBasicBaker1998
    elseif(FlagBasic == "SK1994") then 
      IDBasic = IDBasicSK1994
    end if

    if(FlagDisturbPTemp == "GaussXY" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussXY
    elseif(FlagDisturbPTemp == "GaussXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussXZ
    elseif(FlagDisturbPTemp == "GaussYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempGaussYZ
    elseif(FlagDisturbPTemp == "GaussXYZ") then     
      IDDisturbPTemp = IDDisturbPTempGaussXYZ
    elseif(FlagDisturbPTemp == "Random"  ) then     
      IDDisturbPTemp = IDDisturbPTempRandom
    elseif(FlagDisturbPTemp == "Square"  ) then     
      IDDisturbPTemp = IDDisturbPTempSquare
    elseif(FlagDisturbPTemp == "CosXY" ) then     
      IDDisturbPTemp = IDDisturbPTempCosXY
    elseif(FlagDisturbPTemp == "CosXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCosXZ
    elseif(FlagDisturbPTemp == "CosYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCosYZ
    elseif(FlagDisturbPTemp == "CosXYZ") then     
      IDDisturbPTemp = IDDisturbPTempCosXYZ
    elseif(FlagDisturbPTemp == "ConeXY" ) then     
      IDDisturbPTemp = IDDisturbPTempConeXY
    elseif(FlagDisturbPTemp == "ConeXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempConeXZ
    elseif(FlagDisturbPTemp == "ConeYZ" ) then     
      IDDisturbPTemp = IDDisturbPTempConeYZ
    elseif(FlagDisturbPTemp == "tanh" ) then     
      IDDisturbPTemp = IDDisturbPTempTanh
    elseif(FlagDisturbPTemp == "SK1994XZ" ) then     
      IDDisturbPTemp = IDDisturbPTempSK1994XZ
    elseif(FlagDisturbPTemp == "SK1994YZ" ) then     
      IDDisturbPTemp = IDDisturbPTempSK1994YZ
    elseif(FlagDisturbPTemp == "CircleXZ" ) then     
      IDDisturbPTemp = IDDisturbPTempCircleXZ
    end if

    if(FlagDisturbExner == "GaussXY" ) then     
      IDDisturbExner = IDDisturbExnerGaussXY 
    elseif(FlagDisturbExner == "GaussXZ" ) then     
      IDDisturbExner = IDDisturbExnerGaussXZ 
    elseif(FlagDisturbExner == "GaussYZ" ) then     
      IDDisturbExner = IDDisturbExnerGaussYZ 
    elseif(FlagDisturbExner == "GaussXYZ") then    
      IDDisturbExner = IDDisturbExnerGaussXYZ  
    elseif(FlagDisturbExner == "CosXY" ) then     
      IDDisturbExner = IDDisturbExnerCosXY 
    elseif(FlagDisturbExner == "CosXZ" ) then     
      IDDisturbExner = IDDisturbExnerCosXZ 
    elseif(FlagDisturbExner == "CosYZ" ) then     
      IDDisturbExner = IDDisturbExnerCosYZ 
    elseif(FlagDisturbExner == "CosXYZ") then    
      IDDisturbExner = IDDisturbExnerCosXYZ  
    elseif(FlagDisturbExner == "ConeXY" ) then     
      IDDisturbExner = IDDisturbExnerConeXY 
    elseif(FlagDisturbExner == "ConeXZ" ) then     
      IDDisturbExner = IDDisturbExnerConeXZ 
    elseif(FlagDisturbExner == "ConeYZ" ) then     
      IDDisturbExner = IDDisturbExnerConeYZ 
    end if

    if(FlagDisturbQMix == "GaussXY"  ) then     
      IDDisturbQMix = IDDisturbQMixGaussXY
    elseif(FlagDisturbQMix == "GaussXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixGaussXZ   
    elseif(FlagDisturbQMix == "GaussYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixGaussYZ   
    elseif(FlagDisturbQMix == "GaussXYZ" ) then 
      IDDisturbQMix = IDDisturbQMixGaussXYZ    
    elseif(FlagDisturbQMix == "Dryreg" ) then
      IDDisturbQMix = IDDisturbQMixDryreg
    elseif(FlagDisturbQMix == "Moist"    ) then    
      IDDisturbQMix = IDDisturbQMixMoist
    elseif(FlagDisturbQMix == "CosXY"  ) then     
      IDDisturbQMix = IDDisturbQMixCosXY
    elseif(FlagDisturbQMix == "CosXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixCosXZ   
    elseif(FlagDisturbQMix == "CosYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixCosYZ   
    elseif(FlagDisturbQMix == "CosXYZ" ) then 
      IDDisturbQMix = IDDisturbQMixCosXYZ  
    elseif(FlagDisturbQMix == "ConeXY"  ) then     
      IDDisturbQMix = IDDisturbQMixConeXY
    elseif(FlagDisturbQMix == "ConeXZ"  ) then 
      IDDisturbQMix = IDDisturbQMixConeXZ   
    elseif(FlagDisturbQMix == "ConeYZ"  ) then 
      IDDisturbQMix = IDDisturbQMixConeYZ     
    elseif(FlagDisturbQMix == "Square") then     
      IDDisturbQMix = IDDisturbQMixSquare
    end if

    if (FlagDisturbWind == "Takemi2007") then
      IDDisturbWind = IDDisturbWindTakemi2007
    elseif (FlagDisturbWind == "Sounding") then
      IDDisturbWind = IDDisturbWindSounding
    elseif (FlagDisturbWind == "Zonal") then
      IDDisturbWind = IDDisturbWindZonal
    elseif (FlagDisturbWind == "tanh") then
      IDDisturbWind = IDDisturbWindTanh
    end if

    if(FlagDisturb == "Arare4"        ) then    
      IDDisturb = IDDisturbArare4
    elseif(FlagDisturb == "Arare4mmc" ) then    
      IDDisturb = IDDisturbArare4mmc
    end if

    !----------------------------------------------------
    ! �⥸�塼��ν����
    !
    !   Yamasaki, baker, ����ץ������, ���������ɬ�פʤ�. 
    !
    select case (IDDisturbPTemp)
    case (IDDisturbPTempGaussXY:IDDisturbPTempGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case (IDDisturbPTempRandom)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_random_nml)
      close(unit)

    case (IDDisturbPTempSquare)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)

    case (IDDisturbPTempCosXY:IDDisturbPTempCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbPTempConeXY:IDDisturbPTempConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)

    case (IDDisturbPTempTanh)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_tanh_nml)
      close(unit)

    case (IDDisturbPTempSK1994XZ:IDDisturbPTempSK1994YZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_SK1994_nml)
      close(unit)

    case (IDDisturbPTempCircleXZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_circle_nml)
      close(unit)
    end select

    select case (IDDisturbExner)
    case (IDDisturbExnerGaussXY:IDDisturbExnerGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case (IDDisturbExnerCosXY:IDDisturbExnerCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbExnerConeXY:IDDisturbExnerConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)
    end select

    select case (IDDisturbQMix)
    case (IDDisturbQMixGaussXY:IDDisturbQMixGaussXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_gauss_nml)
      close(unit)

    case(IDDisturbQMixDryreg)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)

    case(IDDisturbQMixMoist)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_moist_nml)
      close(unit)

    case (IDDisturbQMixCosXY:IDDisturbQMixCosXYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cos_nml)
      close(unit)

    case (IDDisturbQMixConeXY:IDDisturbQMixConeYZ)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_cone_nml)
      close(unit)

    case (IDDisturbQMixSquare)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_square_nml)
      close(unit)
      
    end select

    select case (IDDisturbWind)
    case (IDDisturbWindZonal)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_zonal_nml)
      close(unit)

    case (IDDisturbWindTanh)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_disturb_tanh_nml)
      close(unit)
    end select

    select case (IDBasic)
    case (IDBasicDry, IDBasicMoist)
      call FileOpen(unit, file=namelist_filename, mode='r')
      read(unit, NML=initialdata_basic_nml)
      close(unit)
    end select

    if (FlagBasic == "Takemi2007" .OR. FlagDisturbWind == "Takemi2007") then 
      call initialdata_takemi2007_init
    end if

    if (FlagBasic == "Sounding" .OR. FlagDisturbWind == "Sounding") then 
      call initialdata_sounding_init
    end if

    if (     FlagBasic == "Arare4"        &
      & .OR. FlagDisturb == "Arare4"      &
      & .OR. FlagBasic   == "Arare4mmc"   &
      & .OR. FlagDisturb == "Arare4mmc") then 
      call Arare4InitFileio_init
    end if

    call MessageNotify( "M", "main", "FlagBasic        = %c", c1=trim(FlagBasic))
    call MessageNotify( "M", "main", "IDBasic          = %d", i=(/IDBasic/))
    call MessageNotify( "M", "main", "FlagDisturbPTemp = %c", c1=trim(FlagDisturbPTemp))
    call MessageNotify( "M", "main", "IDDisturbPTemp   = %d", i=(/IDDisturbPTemp/))
    call MessageNotify( "M", "main", "FlagDisturbExner = %c", c1=trim(FlagDisturbExner))
    call MessageNotify( "M", "main", "IDDisturbExner   = %d", i=(/IDDisturbExner/))
    call MessageNotify( "M", "main", "FlagDisturbWind  = %c", c1=trim(FlagDisturbWind))
    call MessageNotify( "M", "main", "IDDisturbWind    = %d", i=(/IDDisturbWind/))
    call MessageNotify( "M", "main", "FlagDisturbQMix  = %c", c1=trim(FlagDisturbQMix))
    call MessageNotify( "M", "main", "IDDisturbQMix    = %d", i=(/IDDisturbQMix/))
    call MessageNotify( "M", "main", "FlagDisturb      = %c", c1=trim(FlagDisturb))
    call MessageNotify( "M", "main", "IDDisturb        = %d", i=(/IDDisturb/))

  end subroutine Initialdata_init


!!!--------------------------------------------------------------
  subroutine initialdata_basic_strat(z_TempBZ, z_PressBZ)
    !
    implicit none
    
    real(DP), intent(inout) :: z_TempBZ(kmin:kmax)
    real(DP), intent(inout) :: z_PressBZ(kmin:kmax)
    real(DP)                :: HeightTr
    real(DP)                :: z_DTempDZ(kmin:kmax)
    real(DP)                :: DTempDZ
    real(DP)                :: Weight
    integer                 :: k    
    
    do k = 1, kmax
      z_DTempDZ(k) = (z_TempBZ(k) - z_TempBZ(k-1)) / dz
    end do

    ! ��ή������
    !
    HeightTr =  minval(z_z(1:nz), 1, z_TempBZ(1:nz) < TempTr)
    
    ! ��ή�����̤���ΰ���
    !   �����̤���������絤�Ȥ���. �������絤���������絤�ؤ����ܤ� 
    !   tanh ���Ѥ��Ƥʤ�餫�ˤĤʤ�
    do k = 2, kmax
      
      !�ŤߤĤ��δؿ����Ѱ�. tanh ���Ѥ���
      Weight = ( tanh( (z_Z(k) - HeightTr ) / DHeight ) + 1.0d0 ) * 5.0d-1
      
      !���ͤȤ��Ʋ��٤�׻�����. �����̤���Ǥ� TempTr �������絤�˶�Ť���
      z_TempBZ(k) = z_TempBZ(k) * ( 1.0d0 - Weight ) + TempTr * Weight
      
      !���ٸ�Ψ����Ǯ���ٸ�Ψ��꾮�����ʤ�ʤ��褦��
      DTempDZ = max( z_DTempDZ(k), (z_TempBZ(k) - z_TempBZ(k-1)) / dz )
      
      !���ܾ�β��٤����
      z_TempBZ(k) = z_TempBZ(k-1) + DTempDZ * dz
      
      !���Ϥ��ſ尵ʿ�դ���׻�
      z_PressBZ(k) =                                          &
        &  z_PressBZ(k-1) * ( ( z_TempBZ(k-1) / z_TempBZ(k) ) &
        &    ** (Grav / ( DTempDZ * GasRDry ) ) )
    end do
    
  end subroutine initialdata_basic_strat


!!!--------------------------------------------------------------
  subroutine initialdata_basic_strat_v2(z_TempBZ, z_PressBZ)
    !
    implicit none
    
    real(DP), intent(inout) :: z_TempBZ(kmin:kmax)
    real(DP), intent(inout) :: z_PressBZ(kmin:kmax)
    real(DP)                :: z_DTempDZ(kmin:kmax)
    real(DP)                :: DTempDZ
    integer                 :: k, k1, k2

    ! ���ߤβ��ٸ��ۤ�׻�. 
    !
    do k = 1, kmax
      z_DTempDZ(k) = (z_TempBZ(k) - z_TempBZ(k-1)) / dz
    end do

    ! ��ή�����̤���ΰ���. ���ꤵ�줿���٤Ȥʤ���٤�����̤ȸ��ʤ�
    !
!    k1 = minloc( z_Z, 1, z_Z > HeightTr ) 
    k1 = maxloc( z_TempBZ, 1, z_TempBZ <= TempTr ) 

    ! ���ٸ�Ψ
    !
    k2 = int( DHeight / dz )
    k1 = k1 - k2 / 2          ! ���ꤵ�줿���٤�ޤ���
    DTempDZ = z_DTempDZ(k1-1)

    do k = k1, kmax

      ! ���ٸ�Ψ�򥼥�˶�Ť���. 
      !
      DTempDZ = min( -1.0d-14, DTempDZ - DTempDZ / k2 * ( k - k1 ) )
      
      !���ܾ�β��٤����
      !
      z_TempBZ(k) = z_TempBZ(k-1) + DTempDZ * dz
      
      !���Ϥ��ſ尵ʿ�դ���׻�
      !
      z_PressBZ(k) =                                          &
        &  z_PressBZ(k-1) * ( ( z_TempBZ(k-1) / z_TempBZ(k) ) &
        &    ** (Grav / ( DTempDZ * GasRDry ) ) )

    end do
    
  end subroutine initialdata_basic_strat_v2


  subroutine MainInit

    implicit none

    character(STRING)  :: cfgfile = ""    
    integer, parameter :: OutputRank = 0

    ! gtool �Υ�å�����������
    ! rank0 �Τ�ɸ�����
    !
    call MessageSuppressMPI( OutputRank )

    ! MPI
    !
    call MPIWrapperInit
    
    !NAMELIST �ե�����̾���ɤ߹���
    !
    call argset_init( cfgfile ) !(out)
    
    ! NAMELIST �ե�����̾�Υ⥸�塼��ؤ���Ͽ
    ! Loading NAMELIST file.
    !
    call NmlutilInit( cfgfile ) !(in)
    
    !�ʻ�������ν����
    !  NAMELIST ������������, �ʻ�����׻�����
    !
    call gridset_init
    
    ! ���ط׻��롼����ν����
    ! Initialization of chemical routines.
    !
    call chemcalc_init
    
    ! ���ξ���ν����
    ! Initialization of axis variables.
    !
    call axesset_init
    
    ! ����ξ���ν����
    ! Initialization of constant variables.
    !
    call constants_init
    
    ! ���������ͭ�ѿ��ν����
    ! Initialization of common variables for moist process.
    !
    call composition_init
    
    ! I/O �ե�����̾�ν����
    ! Initialization of output file name. 
    !
    call fileset_init
    
    ! nml ����������Ф� (�������֥롼����)
    !
    call InitialData_init

    ! �ޡ����������ν����
    ! Initialization of margin
    !
    call SetMargin_Init

    !�����ѿ��ν����. �Ȥꤢ���������������ͤ���ꤵ���Ƥ���. 
    !
    call ArareAlloc  
    
  end subroutine MainInit

end program ArareInitData
