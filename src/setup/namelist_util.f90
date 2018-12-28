!= NAMELIST �ե��������Ϥ˴ؤ���桼�ƥ���ƥ�
!
!= Utilities for NAMELIST file input
!
! Authors::   Yoshiyuki O. Takahashi, Yasuhiro MORIKAWA, SUGIYAMA Ko-ichiro
! Version::   $Id: namelist_util.f90,v 1.2 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module namelist_util
  !
  != NAMELIST �ե��������Ϥ˴ؤ���桼�ƥ���ƥ�
  !
  != Utilities for NAMELIST file input
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  !== Variables List
  !
  ! namelist_filename :: NAMELIST �ե������̾��. 
  ! MaxNmlArySize     :: NAMELIST �����ɤ߹�������κ��祵����. 
  ! ------------      :: ------------
  ! namelist_filename :: NAMELIST file name
  ! MaxNmlArySize     :: Maximum size of arrays loaded from NAMELIST
  !
  !== Procedures List
  !
  ! NmlutilInit      :: NAMELIST �ե�����̾������
  ! NmlutilMsg       :: NAMELIST �ե��������Ϥ˴ؤ����å�����ɽ��
  ! NmlutilAryValid  :: NAMELIST �ե����뤫���ɤ߹���������ͭ����������å�
  ! ------------  :: ------------
  ! NmlutilInit      :: Settings of NAMELIST file name
  ! NmlutilMsg       :: Print messages about NAMELIST file input
  ! NmlutilAryValid  :: Check validation of arrays loaded from NAMELIST file

  ! �⥸�塼����� ; USE statements
  !

  ! ���̷��ѥ�᥿
  ! Kind type parameter
  !
  use dc_types,    only: STRING   ! ʸ����. Strings. 

  ! ��å���������
  ! Message output
  !
  use dc_message, only: MessageNotify

  ! ���ʸ ; Declaration statements
  !
  implicit none
  private

  ! ������³��
  ! Public procedure
  !
  public:: NmlutilInit, NmlutilMsg, NmlutilAryValid

  ! �����ѿ�
  ! Public variables
  !
  logical, save, public:: namelist_util_inited = .false.
                              ! �������ե饰. 
                              ! Initialization flag
  character(STRING), save, public:: namelist_filename = ''
                              ! NAMELIST �ե������̾��. 
                              ! NAMELIST file name
  integer, parameter, public:: MaxNmlArySize = 256
                              ! NAMELIST �����ɤ߹�������κ��祵����. 
                              ! Maximum size of arrays loaded from NAMELIST

  ! ������ѿ�
  ! Private variables
  !

  !  NAMELIST �ѿ���
  !  NAMELIST group name
  !
!!$  namelist /namelist_util_nml/ MaxNmlArySize

  character(*), parameter:: module_name = 'namelist_util'
                              ! �⥸�塼���̾��. 
                              ! Module name
  character(*), parameter:: version = &
    & '$Name:  $' // &
    & '$Id: namelist_util.f90,v 1.2 2014/01/20 08:12:41 sugiyama Exp $'
                              ! �⥸�塼��ΥС������
                              ! Module version

contains

  subroutine NmlutilInit(  &
    & namelist_filename_in & ! (in)
    & )
    !
    ! namelist_util �⥸�塼��ν�������Ԥ��ޤ�. 
    !
    ! Initialize "namelist_util" module. 
    !

    ! �⥸�塼����� ; USE statements
    !

    ! �ե��������������
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! ���ʸ ; Declaration statements
    !
    implicit none

    character(*), intent(in ) :: namelist_filename_in
                              ! NAMELIST �ե������̾��. 
                              ! NAMELIST file name

!!$    integer:: unit_nml        ! NAMELIST �ե����륪���ץ��������ֹ�. 
!!$                              ! Unit number for NAMELIST file open
!!$    integer:: iostat_nml      ! NAMELIST �ɤ߹��߻��� IOSTAT. 
!!$                              ! IOSTAT of NAMELIST read

    ! �¹�ʸ ; Executable statement
    !

    if ( namelist_util_inited ) return
    namelist_filename = namelist_filename_in

    ! NAMELIST ���ɤ߹���
    ! NAMELIST is input
    !
!!$    call FileOpen( unit_nml, &          ! (out)
!!$      & namelist_filename, mode = 'r' ) ! (in)
!!$
!!$    rewind( unit_nml )
!!$    read( unit_nml, &               ! (in)
!!$      & nml = namelist_util_nml, &  ! (out)
!!$      & iostat = iostat_nml )       ! (out)
!!$    close( unit_nml )
!!$
!!$    call NmlutilMsg( iostat_nml, module_name ) ! (in)

    ! ���� ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  MaxNmlArySize = %d', i = (/ MaxNmlArySize /) )
    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    ! ������ե饰
    !
    namelist_util_inited = .true.

  end subroutine NmlutilInit


  subroutine NmlutilMsg( iostat, name & ! (in)
    & )
    !
    ! NAMELIST �ե��������ϥ��ơ���������
    ! Ŭ�ڤʥ�å�������ɽ�����ޤ�. 
    !
    ! Appropriate messages are output from 
    ! status of NAMELIST file loading. 
    !

    implicit none
    integer, intent(in):: iostat
                              ! NAMELIST �ɤ߹��߻��Υ��ơ�����. 
                              ! Status of NAMELIST loading
    character(*), intent(in):: name
                              ! �⥸�塼���̾��. 
                              ! Module name

    ! �¹�ʸ ; Executable statement
    !

    if ( iostat == 0 ) then
      call MessageNotify( 'M', name, &
        & 'NAMELIST group "%c" is loaded from "%c".', &
        & c1 = trim(name) // '_nml', &
        & c2 = trim(namelist_filename) )
    else
      call MessageNotify( 'W', name, &
        & 'NAMELIST group "%c" is not found in "%c" (iostat=%d).', &
        & c1 = trim(name) // '_nml', &
        & c2 = trim(namelist_filename), &
        & i = (/iostat/) )
    end if
  end subroutine NmlutilMsg

  subroutine NmlutilAryValid( name, &               ! (in)
    & array, array_name, need_num, need_num_name, & ! (in)
    & valid_limit &                                 ! (in) optional
    )
    !
    ! NAMELIST �����ɤ߹�������󷿥ǡ�������������
    ! �����å����ޤ�. 
    ! 
    ! �ǥե���ȤǤ�, �����ͤ�ͭ���Ȱ����ޤ�. 
    ! ̵���Ǥ���ȸ��ڤ��줿���ˤ�, ���顼��ȯ�������ޤ�. 
    ! 
    ! Check validation of array data loaded from NAMELIST. 
    !
    ! By defaut, positive values are treated as valid values.
    ! If invalidation is checked, an error is occurred. 
    ! 

    ! �⥸�塼����� ; USE statements
    !

    ! ���̷��ѥ�᥿
    ! Kind type parameter
    !
    use dc_types, only: DP      ! �����ټ¿���. Double precision. 

    ! ���ʸ ; Declaration statements
    !
    implicit none
    character(*), intent(in):: name
                              ! ���Υ��֥롼�����ƤӽФ��⥸�塼���̾��. 
                              ! Module name calling this subroutine
    real(DP), intent(in):: array(:)
                              ! ���ڤ��٤�����ǡ���. 
                              ! Checked array data 
    character(*), intent(in):: array_name
                              ! ���ڤ��٤�����ǡ�����̾��. 
                              ! Name of checked array data 
    integer, intent(in):: need_num
                              ! ɬ�פʥǡ�����. 
                              ! 0 ̤���ο���Ϳ����ȥ��顼�������ޤ�. 
                              ! 
                              ! Number of needed data. 
                              ! If number less than 0, an error is occurred. 
    character(*), intent(in):: need_num_name
                              ! ɬ�פʥǡ������򼨤��ѿ���̾��. 
                              ! Name of a variable that indicates number of needed data
    real(DP), intent(in), optional:: valid_limit
                              ! ͭ�������� (�ǥե���Ȥ� 0.0). 
                              ! Lower limit of validation (defalt is 0.0) 

    ! ����ѿ�
    ! Work variables
    !
    real(DP):: valid_limit_work
                              ! ͭ�������� (�ǥե���Ȥ� 0.0). 
                              ! Lower limit of validation (defalt is 0.0) 
    integer:: valid_count     ! ����ǡ�����ͭ���ͤο�. 
                              ! Number of valid values in an array
    integer:: size_array      ! ����ǡ����Υ�����
                              ! Size of array data

    ! �¹�ʸ ; Executable statement
    !

    ! need_num ����ǤϤʤ����Ȥ�����å�
    ! Check that "need_num" is not negative 
    !
    if ( need_num < 0 ) then
      call MessageNotify( 'E', name, '%c=<%d> must not be negative.', &
        &                 c1 = trim(need_num_name), i = (/ need_num /) )
    end if

    ! array �Υ���������ʬ�Ǥ��뤳�Ȥ�����å�
    ! Check that size of "array" is enough
    !
    size_array = size(array)
    if ( need_num > size_array ) then
      call MessageNotify( 'E', name, &
        &  'Maximum size=<%d> of "%c" is too smaller than %c=<%d>. ' // &
        &  'Please search for a statement "MaxNmlArySize = %d" in ' // &
        &  '"namelist_util.f90", and change it into "MaxNmlArySize = %d".', &
!!$        &  'Please add a following phrase to NAMELIST file. ' // &
!!$        &  ' "&namelist_util_nml  MaxNmlArySize=%d /"', &
        &  i = (/ size_array, need_num, MaxNmlArySize, need_num /), &
        &  c1 = trim(array_name), c2 = trim(need_num_name) )
    end if

    ! array ������å�
    ! Check "array"
    !
    if ( need_num > 0 ) then
      valid_limit_work = 0.0_DP
      if ( present( valid_limit ) ) valid_limit_work = valid_limit

      if ( any( array(1:need_num) < valid_limit_work ) ) then
        valid_count = count( .not. ( array(1:need_num) < valid_limit_work ) )
        if ( valid_count > 0 ) then
          call MessageNotify( 'E', name, &
            &   'Number of valid data of %c=<%*f> is %d. ' // &
            &   'Valid data is %c=<%d> necessary.', &
            &   c1 = trim( array_name ), c2 = trim( need_num_name ), &
            &   d = array(1:valid_count), n = (/ valid_count /), &
            &   i = (/ valid_count, need_num /) )
        else
          call MessageNotify( 'E', name, &
            &   'Valid data of %c is nothing. ' // &
            &   'Valid data is %c=<%d> necessary.', &
            &   c1 = trim( array_name ), c2 = trim( need_num_name ), &
            &   i = (/ need_num /) )
        end if
      end if
    end if

  end subroutine NmlutilAryValid

end module namelist_util
