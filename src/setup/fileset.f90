!= �ե�����̾�˴ؤ��������Ԥ�����Υ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: fileset.f90,v 1.10 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module fileset
  !
  != �ե�����̾�˴ؤ��������Ԥ�����Υ⥸�塼��
  !
  !������Ϳ����줿 NAMELIST �ե����뤫��, I/O �ե�����̾�������, 
  !�ݴɤ��뤿����ѿ����ȷ��⥸�塼��
  !

  !�⥸�塼���ɤ߹���
  !
  use dc_types,      only: STRING        

  !���ۤη�����ػ�
  !
  implicit none

  ! �ǥե���Ȥ�°��
  !
  private

  !�����ѿ�
  !
  character(STRING), save, public :: FileTitle = 'cloud moist convection experiment'
                              ! �¸�̾.
                              ! Title of experiment
  character(STRING), save, public :: FileSource = 'deepconv/arare5 (http://www.gfd-dennou.org/library/deepconv)'
                              ! �ǡ����ե���������ץ����̾. 
                              ! Source of data file
  character(STRING), save, public :: FileInstitution = 'GFD Dennou Club (http://www.gfd-dennou.org)'
                              ! �ǡ����ե����������/���롼��.
                              ! Institution or person that changes data files for the last time


  character(STRING), save, public :: ExpTitle          !�ǡ�����ɽ��
  character(STRING), save, public :: ExpSrc            !�ǡ��������������
  character(STRING), save, public :: ExpInst           !�ǽ��ѹ��ԡ��ȿ�

  ! ������롼����θ���
  !
  public fileset_init

contains

  subroutine fileset_init
    !
    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    !
    
    !�⥸�塼���ɤ߹���
    use dc_iounit,     only: FileOpen      
    use dc_message,    only: MessageNotify 
    use namelist_util, only: namelist_filename
    
    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    integer                       :: unit     !����ե������������ֹ�
     
    !����ե����뤫���ɤ߹�����ϥե��������
    !
    NAMELIST /fileset_nml/ FileTitle, FileSource, FileInstitution
    
    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=fileset_nml)
    close(unit)
    
    ! �ɤ߹������������
    !
    call MessageNotify( "M", &
      & "fileset_init", "FileTitle = %c",    c1=trim(FileTitle) )
    call MessageNotify( "M", &
      & "fileset_init", "FileSource = %c",   c1=trim(FileSource) )
    call MessageNotify( "M", &
      & "fileset_init", "FileInstitution = %c", c1=trim(FileInstitution) )

    !----------------------------------------------------
    ! arare4 �Ѥ��ѿ�
    !
    ExpTitle= FileTitle
    ExpSrc  = FileSource 
    ExpInst = FileInstitution

  end subroutine fileset_init

end module fileset
