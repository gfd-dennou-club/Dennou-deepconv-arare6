!= ���������ѥ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: timeset.f90,v 1.11 2014/01/21 05:00:57 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module timeset
  !
  != ���������ѥ⥸�塼��
  !
  !������Ϳ����줿 NAMELIST �ե����뤫��, ����˴ؤ������������, 
  !�ݴɤ��뤿����ѿ����⥸�塼��
  !
  !== Procedures List
  ! timeset_init       :: ������롼����
  ! TimesetDelTimeHalf :: �꡼�ץե�å���������κǽ�� 1 ���ƥåפ�󤹤���Υ롼����
  ! TimesetProgress    :: �꡼�ץե�å���������� 1 ���ƥåפ�󤹤���Υ롼����
  
  !�⥸�塼���ɤ߹���
  use dc_types,      only: DP
  use dc_iounit,     only: FileOpen
  use dc_message,    only: MessageNotify
  use namelist_util, only: namelist_filename
  
  !���ۤη�����ػ�
  implicit none

  !°��
  private
  
  ! Public Interface
  real(DP), save, public  :: TimeA                   !���� t + \del t
  real(DP), save, public  :: TimeN                   !���� t
  real(DP), save, public  :: TimeB                   !���� t - \del t
  real(DP), save, public  :: DelTimeLong  = 2.0d0    !Ĺ�������ॹ�ƥå�
  real(DP), save          :: DelTimeLongSave  = 2.0d0 !Ĺ�������ॹ�ƥå�
  real(DP), save, public  :: DelTimeShort = 2.0d-1   !û�������ॹ�ƥå�
  real(DP), save, public  :: RestartTime  = 0.0d0    !�׻����ϻ���
  real(DP), save, public  :: IntegPeriod  = 3600.0d0 !��ʬ����
  real(DP), save, public  :: EndTime      = 3600.0d0 !�׻���λ����
  real(DP), save, public  :: DelTimeOutput= 2.0d0    !���ϥ����ॹ�ƥå�
  real(DP), save, public  :: tfil         = 1.0d-1   !�������λ��֥ե��륿�η���
  real(DP), parameter     :: tfilSave     = 1.0d-1   !�������λ��֥ե��륿�η���
  integer,  save, public  :: NstepShort = 20         !û�������ॹ�ƥåפΥ��ƥå׿�
  integer,  save          :: NstepShortSave = 20     !û�������ॹ�ƥåפΥ��ƥå׿�
  integer,  save, public  :: NstepOutput    = 20     !�ꥹ�����ȥե�����ؤν���
  logical,  save, public  :: FlagInitialRun = .false.!t=0 ���ݤ��Υե饰. 

  integer,  save, public  :: NstepDisp = 20          !�ꥹ�����ȥե�����ؤν���

  real(DP), save, public  :: TimeInt  = 3600.0d0     !��ʬ����
  
  ! ��������
  public timeset_init, TimesetDelTimeHalf, TimesetProgress

contains
   
  subroutine timeset_init
    !
    != ������롼����
    !
    ! NAMELIST ����ɬ�פʾ�����ɤ߼��, ���ִ�Ϣ���ѿ��������Ԥ�. 
    !

    !���ۤη�����ػ�
    implicit none

    !�����ѿ�
    integer    :: unit

    !---------------------------------------------------------------    
    ! NAMELIST �����������
    !
    NAMELIST /timeset_nml/ &
      & DelTimeLong, DelTimeShort, IntegPeriod, RestartTime, DelTimeOutput
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=timeset_nml)
    close(unit)

    ! �׻���λ����
    !
    EndTime = RestartTime + IntegPeriod
    
    ! ��������ॹ�ƥåפ�����
    !   �¿��γ�껻�ʤΤ�, ǰ�ΰ٤�. 
    !    
    NstepShort = ( nint( 2.0d0 * DelTimeLong * 1.0d3 ) / nint( DelTimeShort * 1.0d3 ) )

    ! ���ֹ�ߤ�����
    !   �������ͤ��ݴɤ��Ƥ���
    DelTimeLongSave = DelTimeLong
    NstepShortSave  = NstepShort

    ! �ꥹ�����ȥե������񤭽Ф������ߥ�
    !
    NstepOutput = nint( DelTimeOutput * 1.0d2 ) / nint( DelTimeLong * 1.0d2 )
    if ( NstepOutput < 0 ) then 
       NstepOutput = nint( DelTimeOutput * 1.0d1 ) / nint( DelTimeLong * 1.0d1 )
    end if
    if ( NstepOutput < 0 ) then 
       NstepOutput = nint( DelTimeOutput ) / nint( DelTimeLong )
    end if
    if ( NstepOutput < 0 ) then 
      call MessageNotify( "E", &
        & "timeset_init", "NstepOutput is negative, %d", i=(/NstepOutput/) )
    end if


    ! ���������
    !
    TimeB = RestartTime - DelTimeLong
    TimeN = RestartTime
    TimeA = RestartTime + DelTimeLong      
    
    ! �ꥹ�����Ȥ��ݤ�. �ꥹ�����Ȥʤ�  .false.
    !
    if ( nint( RestartTime * 1.0d3 ) == 0 ) then 
      FlagInitialRun = .true.
    else
      FlagInitialRun = .false.
    end if

    !---------------------------------------------------------------
    ! ��ǧ
    !
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeShort = %f", d=(/DelTimeShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "Restarttime  = %f", d=(/Restarttime/)  )
    call MessageNotify( "M", &
      & "timeset_init", "IntegPeriod  = %f", d=(/IntegPeriod/) )
    call MessageNotify( "M", &
      & "timeset_init", "EndTime      = %f", d=(/EndTime/) )
    call MessageNotify( "M", &
      & "timeset_init", "DelTimeOutput= %f", d=(/DelTimeOutput/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepShort   = %d", i=(/NstepShort/) )
    call MessageNotify( "M", &
      & "timeset_init", "NstepOutput  = %d", i=(/NstepOutput/) )
    call MessageNotify( "M", &
      & "timeset_init", "tfil  = %f", d=(/tfil/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeA = %f", d=(/TimeA/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeN = %f", d=(/TimeN/) )
    call MessageNotify( "M", &
      & "timeset_init", "TimeB = %f", d=(/TimeB/) )
    call MessageNotify( "M", &
      & "timeset_init", "FlagInitialRun  = %b", L=(/FlagInitialRun/) )
    
    !---------------------------------------------------------------
    ! arare4 �Ѥ��ѿ�������
    !
    NStepDisp = NStepOutput
    TimeInt   = IntegPeriod

  end subroutine timeset_init


  subroutine TimesetDelTimeHalf    
    !
    != �꡼�ץե�å���������κǽ�� 1 ���ƥåפ�󤹤���Υ롼����
    !
    ! �ǽ�ΰ���ܤϥ����顼��������ǲ󤹤Τ�, 
    ! Ĺ�����ֹ�ߤ�Ⱦʬ�ˤ�, û�����֤ǲ󤹥롼�ײ����Ⱦʬ�ˤ���.
    ! Asselin �λ��֥ե��륿�η����򥼥�Ȥ���. 

    implicit none

    TimeB       = TimeN
    DelTimeLong = DelTimeLongSave * 0.5d0
    NstepShort  = NstepShortSave  / 2
    tfil        = 0.0d0

    !---------------------------------------------------------------
    ! ��ǧ
    !
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial DelTimeLong  = %f", d=(/DelTimeLong/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Initial NstepShort   = %d", i=(/NstepShort/) )
    call MessageNotify( "M",               &
      & "timeset_init_TimesetDelTimeHalf", &
      & "Asselin Time Filter coefficient = %f", d=(/tfil/) )

  end subroutine TimesetDelTimeHalf
  

  subroutine TimesetProgress
    !
    ! = �꡼�ץե�å���������� 1 ���ƥåפ�󤹤���Υ롼����
    !
    ! �����ʤ�, ����ι������ Asselin �λ��֥ե��륿�η�����ľ��. 
    ! 

    implicit none

    ! ����������ľ�� 
    !
    DelTimeLong = DelTimeLongSave
    NstepShort  = NstepShortSave

    ! �����ʤ��
    !
    TimeB = TimeN
    TimeN = TimeA
    TimeA = TimeA + DelTimeLong

    ! Asselin �λ��֥ե��륿�η���
    tfil = tfilSave

  end subroutine TimesetProgress
  
end module timeset
