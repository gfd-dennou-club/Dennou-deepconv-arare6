!= Module ReStartFileIO
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: restartfileio.f90,v 1.18 2015/02/19 02:17:23 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module ReStartFileIO
  !
  !�ꥹ�������Ѥξ�ξ���� netCDF �ե�����˽��Ϥ��뤿��Υ롼����
  !

  !�⥸�塼���ɤ߹���
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING

  !���ۤη�����ػ�
  implicit none
  
  !�ؿ��� public �˻���
  public ReStartFileIO_init
  public ReStartFileio_Finalize
  public ReStartFileio_BasicZ_Get
  public ReStartFileio_Var_Get

  type(GT_HISTORY), save, public   :: rstat
  character(STRING), save, private :: InitialFile = ""
  character(STRING), save, private :: InputFile   = ""
  character(STRING), save, private :: OutputFile  = "output.nc"
  
contains

  subroutine ReStartFileio_Init ( FlagInitData )
    !
    !�ꥹ�����ȥե�����ν񤭽Ф�
    !

    !�⥸�塼���ɤ߹���
    use gtool_history, only : HistoryCreate,  &
      &                       HistoryPut,     &
      &                       HistoryAddVariable
    use dc_message,    only : MessageNotify
    use dc_iounit,     only : FileOpen
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use axesset,       only : x_X,            &!X ��ɸ��(�����顼�ʻ���)
      &                       y_Y,            &!Y ��ɸ��(�����顼�ʻ���)
      &                       z_Z              !Z ��ɸ��(�����顼�ʻ���)
    use gridset,       only : ncmax            !�Ž���ʬ�ο�
    use fileset,       only : filetitle,      &!�ǡ�����ɽ��
      &                       filesource,     &!�ǡ��������������
      &                       FileInstitution  !�ǽ��ѹ��ԡ��ȿ�
    use namelist_util, only: namelist_filename
    
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    logical, intent(in), optional :: FlagInitData  ! ����������ץ����ξ��� .true.

    real(4)            :: SpcID(ncmax)
    integer            :: N, L, M
    integer            :: s    
    integer            :: unit     !�����ֹ�
     
    !NAMELIST �����������
    NAMELIST /restartfileio_nml/ InitialFile, InputFile, OutputFile
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=restartfileio_nml)
    close(unit)

    !��ǧ
    call MessageNotify( "M", &
      & "restartfileioIO_init", "InitialFile  = %c", c1=trim(InitialFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "InputFile  = %c", c1=trim(InputFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "OutputFile = %c", c1=trim(OutputFile) )

    ! InputFile �����ξ���, ����ͥե������Ȥ�. 
    if ( InputFile == "" ) then 
      InputFile = InitialFile
    end if

    !����������ξ��ˤ�, outputfile ̾���ѹ�
    if ( present( FlagInitData ) ) then 
      if ( FlagInitData ) then 

        if ( .NOT. InitialFile == "" ) then  ! �ե�����̾�����Ǥʤ����
          OutputFile = InitialFile
          InputFile  = ""
        end if

        if ( InitialFile == "" .AND. (.NOT. InputFile == "") ) then  
          OutputFile = InputFile
          InputFile  = ""
        end if

      end if
    end if

    !��ǧ
    call MessageNotify( "M", &
      & "restartfileioIO_init", "INPUT  = %c", c1=trim(InputFile) )
    call MessageNotify( "M", &
      & "restartfileioIO_init", "OUTPUT = %c", c1=trim(OutputFile) )

    SpcID = 0.0d0
    do s = 1, ncmax
      SpcID(s) = real( s, 4 )
    end do
    
    N = size(x_X, 1)
    L = size(y_Y, 1)
    M = size(z_Z, 1)
    
    !-------------------------------------------------------------    
    ! �ҥ��ȥ꡼����
    !-------------------------------------------------------------  
    call HistoryCreate(                              &
      & file = Outputfile,                           &
      & title = filetitle,                           &
      & source = filesource,                         &
      & institution = FileInstitution,               &
      & dims=(/'x','y','z','s','t'/),                &
      & dimsizes=(/N, L, M, ncmax, 0/),              &
      & longnames=(/'X-coordinate',                  &
      &             'Y-coordinate',                  &
      &             'Z-coordinate',                  &
      &             'Species Num ',                  &
      &             'Time        '/),                &
      & units=(/'m  ','m  ','m  ','1  ','sec'/),     &
      & xtypes=(/'double', 'double', 'double', 'double', 'double'/), &
      & flag_mpi_split = FLAG_LIB_MPI,               &
      & origin=0.0, interval=1.0,                    &
      & history=rstat, quiet=.true. )
    
    !-------------------------------------------------------------  
    ! �ѿ�����
    !-------------------------------------------------------------
    call HistoryPut('x', x_X, rstat )
    call HistoryPut('y', y_Y, rstat )
    call HistoryPut('z', z_Z, rstat )
    call HistoryPut('s', real(SpcID, 4), rstat )

    !̵�������Ϥδ��ܾ�
    call HistoryAddVariable(                           &
      & varname='ExnerBZ', dims=(/'x','y','z'/),       &
      & longname='nondimensional pressure', units='1', &
      & xtype='double', history=rstat )
    
    !���̤δ��ܾ�
    call HistoryAddVariable(                           &
      & varname='PTempBZ', dims=(/'x','y','z'/),       &
      & longname='potential temperature',              &
      & units='K', xtype='double', history=rstat ) 

    !���̤δ��ܾ�
    call HistoryAddVariable(                           &
      & varname='VPTempBZ', dims=(/'x','y','z'/),      &
      & longname='virtual potential temperature',      &
      & units='K', xtype='double', history=rstat ) 
    
    !̩�٤δ��ܾ�
    call HistoryAddVariable(                           &
      & varname='DensBZ', dims=(/'x','y','z'/),        &
      & longname='density',                            &
      & units='Kg.m-3', xtype='double', history=rstat )
    
    !����®�٤δ��ܾ�
    call HistoryAddVariable(                           &
      & varname='VelSoundBZ', dims=(/'x','y','z'/),    &
      & longname='sound velocity',                     &
      & units='m.s-2', xtype='double', history=rstat )

    !���٤δ��ܾ�
    call HistoryAddVariable(                           &
      & varname='TempBZ', dims=(/'x','y','z'/),        &
      & longname='Temperature of basic state',         &
      & units='K', xtype='double', history=rstat ) 
    
    !���Ϥδ��ܾ�
    call HistoryAddVariable(                           &
      & varname='PressBZ', dims=(/'x','y','z'/),       &
      & longname='Pressure of basic state',            &
      & units='Pa', xtype='double', history=rstat ) 
    
    !�����������δ��ܾ�
    call HistoryAddVariable(                              &
      & varname='QMixBZ', dims=(/'x','y','z','s'/),       &
      & longname='Mixing ratio of Condensible volatiles', &
      & units='kg.kg-1', xtype='double', history=rstat ) 
    
    !ʬ���̸���
    call HistoryAddVariable(                         &
      & varname='EffMolWtBZ', dims=(/'x','y','z'/),  &
      & longname='Effect of Mole Weight',            &
      & units='1', xtype='double', history=rstat ) 

    !����
    call HistoryAddVariable(                        &
      & varname='HumBZ', dims=(/'x','y','z','s'/),  &
      & longname='Humidity',                        &
      & units='1', xtype='double', history=rstat ) 
    
    !®��
    call HistoryAddVariable(                         &
      & varname='VelX', dims=(/'x','y','z','t'/),    &
      & longname='zonal velocity',                   &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !®��
    call HistoryAddVariable(                         &
      & varname='VelY', dims=(/'x','y','z','t'/),    &
      & longname='meridional velocity',              &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !®��
    call HistoryAddVariable(                         &
      & varname='VelZ', dims=(/'x','y','z','t'/),    &
      & longname='vertical velocity',                &
      & units='m.s-1',                               &
      & xtype='double', history=rstat )
    
    !̵��������
    call HistoryAddVariable(                         &
      & varname='Exner', dims=(/'x','y','z','t'/),   &
      & longname='nondimensional pressure',          &
      & units='1',                                   &
      & xtype='double', history=rstat )
    
    !���̤ξ���
    call HistoryAddVariable(                         &
      & varname='PTemp', dims=(/'x','y','z','t'/),   &
      & longname='virtual potential temperature',    &
      & units='K',                                   &
      & xtype='double', history=rstat )
    
    !��Ǵ������
    call HistoryAddVariable(                         &
      & varname='Km', dims=(/'x','y','z','t'/),      &
      & longname='Km',                               &
      & units='m2.s-1',                              &
      & xtype='double', history=rstat )
    
    !��Ǵ������
    call HistoryAddVariable(                         &
      & varname='Kh', dims=(/'x','y','z','t'/),      &
      & longname='Kh',                               &
      & units='m2.s-1',                              &
      & xtype='double', history=rstat )

    !��̩��
    call HistoryAddVariable(                         &
      & varname='CDens', dims=(/'x','y','z','t'/),   &
      & longname='CDens',                            &
      & units='m2.s-1',                              &
      & xtype='double', history=rstat )
    
    !������
    call HistoryAddVariable(                         &
      & varname='QMix', dims=(/'x','y','z','s','t'/),&
      & longname='Mixing Ratio',                     &
      & units='kg.kg-1"',                            & 
      & xtype='double', history=rstat )
    
  end subroutine ReStartFileio_Init
  
  
  subroutine ReStartFileIO_Finalize
    !
    !�ꥹ�����ȥե�����Υ�����
    !

    !�⥸�塼���ɤ߹���
    use gtool_history, only: HistoryClose
    
    !���ۤη�����ػ�
    implicit none
    
    !�ե�������Ĥ���
    call HistoryClose(rstat, quiet=.true.)
    
  end subroutine ReStartFileIO_Finalize
  

  subroutine ReStartFileio_Var_Get(   &
    & pyz_VelXB,  pyz_VelXN,    & ! (out)
    & xqz_VelYB,  xqz_VelYN,    & ! (out)
    & xyr_VelZB,  xyr_VelZN,    & ! (out)
    & xyz_PTempB, xyz_PTempN,   & ! (out)
    & xyz_ExnerB, xyz_ExnerN,   & ! (out)
    & xyzf_QMixB, xyzf_QMixN,   & ! (out)
    & xyz_KmB,    xyz_KmN,      & ! (out)
    & xyz_KhB,    xyz_KhN,      & ! (out)
    & xyz_CDensB, xyz_CDensN  )   ! (out)
    !
    !�ꥹ�����ȥե����뤫��������
    !

    !�⥸�塼���ɤ߹���
    use dc_message,    only : MessageNotify
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use dc_string,     only : toChar
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!���󥵥��� (X ����)
      &                       jmin, jmax,    &!���󥵥��� (Y ����)
      &                       kmin, kmax,    &!���󥵥��� (Z ����)
      &                       ncmax           !�Ž���ʬ�ο�
    use timeset,       only : DelTimeLong,   &
      &                       RestartTime 

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(out) :: pyz_VelXN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelYN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_ExnerN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTempN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KmN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KhN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDensN &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMixN &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP), intent(out) :: pyz_VelXB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelYB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyr_VelZB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_ExnerB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_PTempB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KmB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_KhB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyz_CDensB &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xyzf_QMixB &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)

    character(STRING)    :: name               !�ѿ�̾
    character(STRING)    :: TimeN = ""
    character(STRING)    :: TimeB = ""

    integer, parameter  :: rtrn = 1000
    real(DP)            :: RTime(rtrn)
    real(DP)            :: xyz_var(imin:imax,jmin:jmax,kmin:kmax)
    integer             :: nt, nt2, rNum(1)
    character(STRING)   :: ct

    !-------------------------------------------------------------
    ! Get a Value from netCDF File
    !-------------------------------------------------------------    

    ! netCDF ���ݴɤ���Ƥ���������Ф�. �ѿ��Ȥ��Ʋ��̤�����. 
    !
    RTime = -999.0    ! ����� (��»��)
    name  = "PTemp"   ! �����
    !
    CHECK_RestartTime: do nt = 1, rtrn

      ! ����μ��Ф�
      !
      ct = 't=^' // toChar( nt )
      call HistoryGet( InputFile, name, xyz_var, range=ct, &
        & flag_mpi_split = FLAG_LIB_MPI, returned_time = RTime(nt) )

      ! ���Ф������郎���פ�����롼�׽�λ. nt2 �˥롼�ײ���ݴ�. 
      !
      if ( nt > 1 .AND. RTime(nt - 1) == RTime(nt) ) then 
        nt2 = nt - 1
        exit CHECK_RestartTime
      end if
    end do CHECK_RestartTime

    ! RestartTime �Ȥ��ƻ��ꤵ��Ƥ�����������ź���������. 
    ! gtool ���ʤ�����Τ�, ���и����Ǿ��ˤʤ�����ź���Ȥ���.
    !
    RTime(1:nt2) = abs( RTime(1:nt2) - RestartTime )
    rNum = minloc( RTime(1:nt2) )
    
    ! �ꥹ�����ȤǻȤ������ź������ꤹ��. 
    ! ���ߤλ��� (TimeN) �� 1 ���ƥå����λ��� (TimeB) ��ź��
    !
    if (rNum(1) == 1) then 
      timeB = 't=^' // toChar( rNum(1)     )
      timeN = 't=^' // toChar( rNum(1)     )
    else
      timeB = 't=^' // toChar( rNum(1)     )
      timeN = 't=^' // toChar( rNum(1) + 1 )
    end if

    ! ��ǧ
    !
    call MessageNotify( "M", "restartfileio_var_get", "rNum  = %d", i=(/rNum/) )
    call MessageNotify( "M", "restartfileio_var_get", "timeB = %c", c1=trim(timeB) )
    call MessageNotify( "M", "restartfileio_var_get", "timeN = %c", c1=trim(timeN) )

!!!
!!! �ե����륪���ץ� & �ͤμ��Ф�
!!!
    name = "PTemp"
    call HistoryGet( InputFile, name, xyz_PTempB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_PTempN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelX"
    call HistoryGet( InputFile, name, pyz_VelXB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, pyz_VelXN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelY"
    call HistoryGet( InputFile, name, xqz_VelYB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xqz_VelYN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "VelZ"
    call HistoryGet( InputFile, name, xyr_VelZB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyr_VelZN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Exner"
    call HistoryGet( InputFile, name, xyz_ExnerB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_ExnerN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "Km"
    call HistoryGet( InputFile, name, xyz_KmB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_KmN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
    
    name = "Kh"
    call HistoryGet( InputFile, name, xyz_KhB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_KhN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

    name = "CDens"
    call HistoryGet( InputFile, name, xyz_CDensB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyz_CDensN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )
      
    name = "QMix"
    call HistoryGet( InputFile, name, xyzf_QMixB, range=TimeB, flag_mpi_split = FLAG_LIB_MPI )
    call HistoryGet( InputFile, name, xyzf_QMixN, range=TimeN, flag_mpi_split = FLAG_LIB_MPI )

  end subroutine ReStartFileio_Var_Get


  subroutine ReStartFileio_BasicZ_Get()
    !
    !�ꥹ�����ȥե����뤫��������
    !

    !�⥸�塼���ɤ߹���
    use gtool_history, only : HistoryGet
    use dc_types,      only : DP, STRING
    use mpi_wrapper,   only : FLAG_LIB_MPI
    use gridset,       only : imin, imax,    &!���󥵥��� (X ����)
      &                       jmin, jmax,    &!���󥵥��� (Y ����)
      &                       kmin, kmax,    &!���󥵥��� (Z ����)
      &                       ncmax           !�Ž���ʬ�ο�
    use basicset,      only : basicset_init 

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP)              :: Var3D &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: Var4D &
      &                     (imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)              :: xyz_ExnerBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ�Υ������ʡ��ؿ�
    real(DP)              :: xyz_DensBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ��̩��
    real(DP)              :: xyz_PTempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ�β���
    real(DP)              :: xyz_VelSoundBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ�β�®
    real(DP)              :: xyz_PressBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ�ΰ���
    real(DP)              :: xyz_TempBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ�β���
    real(DP)              :: xyzf_QMixBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax, ncmax)
                                               !���ܾ�κ�����
    real(DP)              :: xyz_EffMolWtBZ &
      &                     (imin:imax,jmin:jmax,kmin:kmax)
                                               !���ܾ��ʬ���̸���
    character(STRING)    :: name               !�ѿ�̾

    !-------------------------------------------------------------
    ! ���ܾ�μ���
    !
    name = "DensBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_DensBZ = Var3D
    
    name = "ExnerBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_ExnerBZ = Var3D
    
    name = "PTempBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PTempBZ = Var3D
    
    name = "VelSoundBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_VelSoundBZ = Var3D
    
    name = "TempBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_TempBZ = Var3D
    
    name = "PressBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_PressBZ = Var3D
    
    name = "QMixBZ"
    call HistoryGet( InputFile, name, Var4D, flag_mpi_split = FLAG_LIB_MPI )
    xyzf_QMixBZ = Var4D
    
    name = "EffMolWtBZ"
    call HistoryGet( InputFile, name, Var3D, flag_mpi_split = FLAG_LIB_MPI )
    xyz_EffMolWtBZ = Var3D
    
    !----------------------------------------------------------
    ! BasicSet �⥸�塼����ͤ�����
    !----------------------------------------------------------
    call basicset_init( &
      & xyz_PressBZ,    &
      & xyz_ExnerBZ,    &
      & xyz_TempBZ,     &
      & xyz_PTempBZ,    &
      & xyz_DensBZ,     &
      & xyz_VelSoundBZ, &
      & xyzf_QMixBZ,    &
      & xyz_EffMolWtBZ  &
      & )
    
  end subroutine ReStartFileio_BasicZ_Get
       
end module ReStartFileIO
