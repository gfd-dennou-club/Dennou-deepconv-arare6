!= ���ޥ�ɥ饤������β���Ԥ��⥸�塼��
!
! Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
! Version::   $Id: argset.f90,v 1.4 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module argset
  !
  != ���ޥ�ɥ饤������β���Ԥ��⥸�塼��
  !
  != Procedures List
  !
  ! argset_init :: ������롼����
  !

  !�⥸�塼�����
  !
  use dc_types, only : STRING
  use dc_args,  only : ARGS
  
  !���ۤη�����ػ�
  !
  implicit none

  !�����ѿ�
  !-----   ���ޥ�ɥ饤������������ѿ�    -----
  type(ARGS), save, private        :: arg             ! ���ޥ�ɥ饤���������
  logical, save, private           :: OPT_namelist    ! ���ޥ�ɥ饤������������ѿ�
  character(STRING), save, private :: VAL_namelist    ! ���ޥ�ɥ饤���������

contains

  subroutine argset_init(namelist_filename)
    !
    !���ޥ�ɥ饤��������ᤷ��, ������Ϳ����줿 NAMELIST �ե�����̾
    !���֤�. 
    !

    !�⥸�塼�����
    !
    use dc_types,    only : STRING
    use dc_string,   only : StoA
    use dc_args,     only : ARGS, Open, Debug, Help, Strict, Close, Option
    use dc_message,  only : MessageNotify
    
    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    character(STRING), intent(out)     :: namelist_filename

    !NAMELIST �ե�����μ���
    !  gtool5 �饤�֥��� dc_args �⥸�塼�������.
    !  �����ǽ�ʥ��ץ����ϰʲ����̤�.
    !
    !  -N=VAL, --namelist=VAL
    !    specify Namelist file (default is 'arare.conf'). 
    !
    !  -D=VAL, --debug=VAL
    !    call dc_trace#SetDebug (display a lot of messages for debug).
    !    VAL is unit number (default is standard output)
    !
    !  -h=VAL, -H=VAL, --help=VAL
    !    display this help and exit. VAL is unit number (default is
    !    standard output)
    !
    call Open(arg)
    call Option(arg, StoA('-N', '--namelist'), &
      &         OPT_namelist, VAL_namelist, &
      &         help="specify Namelist file (default is 'arare.conf')." )
                      ! "-N/--namelist" ���ץ���������
    call Debug(arg)   ! �ǥХå����ץ����μ�ư����
    call Help(arg)    ! �إ�ץ��ץ����μ�ư����
    call Strict(arg)  ! ̵���ʥ��ץ���������˷ٹ��ɽ��

    !"-N/-namelist" ���ץ����β��
    !  Ϳ�����Ƥ��ʤ����ϥǥե������ (arare.conf) �� 
    !  NAMLIST �ե�����̾�Ȥ���.
    !
    if (OPT_namelist) then
      call MessageNotify( "M", "main", &
        &                 "Namelist file is '%c'", c1=trim(VAL_namelist) )
      namelist_filename=trim(VAL_namelist)
    else
      call MessageNotify( "W", "main", &
        &                 "Namelist file is not specified." )
      call MessageNotify( "M", "main", &
        &                 "Use default Namelist file (arare.conf)." )
      namelist_filename="arare.conf"
    end if

    call Close(arg)

    ! ��ǧ
    !
    call MessageNotify( "M", "argset_init",   &
      &                 "NAMELIST FILE = %c", c1=namelist_filename )

  end subroutine argset_init

end module argset
