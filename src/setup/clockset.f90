! CPU ���֤η�¬��Ԥ�����Υ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: clockset.f90,v 1.3 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module clockset
  !
  != CPU ���֤η�¬��Ԥ�����Υ⥸�塼��
  !
  ! CPU ���֤η�¬��Ԥ�����Υ��֥롼�����«�ͤ��⥸�塼��
  !
  !== Procedures List
  ! ClocksetInit      :: ������롼����
  ! ClocksetClose     :: CPU ���ַ�¬�ν�λ����
  ! ClocksetPredict   :: �ץ���ब��λ�ޤǤ�ͽ¬ CPU ���֡�������ɽ��
  ! ClockSetPreStart  :: CPU ���ַ�¬���� (������롼����)
  ! ClockSetPreStop   :: CPU ���ַ�¬��λ (������롼����)
  ! ClockSetLoopStart :: CPU ���ַ�¬���� (���֥롼��)
  ! ClocksetLoopStop  :: CPU ���ַ�¬��λ (���֥롼��)


  !�⥸�塼���ɤ߹���
  !
  use dc_clock,   only : CLOCK

  !���ۤη�����ػ�
  !
  implicit none

  !°���λ���
  !
  private
  
  ! Public Interface
  !
  type(CLOCK), save :: clock_init, clock_loop  ! Variables for CPU time counting 
                                               ! CPU ���ַ�¬���ѿ�
  public ClocksetInit, ClocksetClose, ClocksetPredict
  public ClockSetPreStart, ClockSetPreStop, ClockSetLoopStart, ClocksetLoopStop

contains
   
  subroutine ClocksetInit
    ! 
    ! CPU ���ַ�¬ ������롼����
    !

    !�⥸�塼���ɤ߹���
    use dc_clock,    only : DCClockCreate

    implicit none
    
    ! ������롼�����Ѥλ����¬�����
    !
    call DCClockCreate( &           ! Initialize (�����)
      & clk = clock_init, &         ! (out)
      & name = 'initialization' )   ! (in)

    ! ����ȯŸ�롼���Ѥλ����¬�����
    !
    call DCClockCreate( &           ! Initialize (�����)
      & clk = clock_loop, &         ! (out)
      & name = 'time-integration' ) ! (in)

  end subroutine ClocksetInit


  subroutine ClocksetPreStart
    !
    ! CPU ���ַ�¬���� (������롼����)
    !

    !�⥸�塼���ɤ߹���
    use dc_clock,   only : DCClockStart
    
    !���ۤη�����ػ�
    implicit none

    ! Start CPU time counting     
    call DCClockStart(clk = clock_init) ! (inout)
    
  end subroutine ClocksetPreStart
  

  subroutine ClocksetPreStop
    !
    ! CPU ���ַ�¬��λ (������롼����)
    !
    
    !�⥸�塼���ɤ߹���
    use dc_clock,   only : DCClockStop

    !���ۤη�����ػ�
    implicit none
    
    ! Stop CPU time counting 
    ! (CPU ���ַ�¬��λ)    
    call DCClockStop(clk = clock_init) ! (inout)

  end subroutine ClocksetPreStop


  subroutine ClocksetLoopStart
    !
    ! CPU ���ַ�¬���� (���֥롼��)
    !
    
    !�⥸�塼���ɤ߹���
    use dc_clock,   only : DCClockStart

    !���ۤη�����ػ�
    implicit none

    ! Start CPU time counting 
    call DCClockStart(clk = clock_loop) ! (inout) 

  end subroutine ClocksetLoopStart


  subroutine ClocksetLoopStop
    !
    ! CPU ���ַ�¬��λ (���֥롼��)
    !

    !�⥸�塼���ɤ߹���
    use dc_clock,   only : DCClockStop

    !���ۤη�����ػ�
    implicit none

    !Stop CPU time counting 
    !(CPU ���ַ�¬��λ)    
    call DCClockStop(clk = clock_loop) ! (inout)

  end subroutine ClocksetLoopStop


  subroutine ClocksetPredict
    !
    ! �ץ���ब��λ����ޤǤ�ͽ¬ CPU ����, �����������ɽ��
    ! ��ö��¬��ߤ�Ƥ���, ͽ���ͤ�ɽ����, ��¬��Ƴ�����. 
    !

    !�⥸�塼���ɤ߹���
    use dc_message,  only : MessageNotify
    use dc_clock,    only : DCClockStop, DCClockPredict, operator(+)
    use timeset,     only : TimeA, RestartTime, IntegPeriod
    
    !���ۤη�����ػ�
    implicit none

    !����ѿ�
    real(4) :: progress 

    ! �׻��οʹԤγ���׻�
    progress = real((TimeA - RestartTime) / IntegPeriod, 4)

    ! CPU time measurement stops, temporarily.
    call ClocksetLoopStop

    call MessageNotify( "M", "ClockSet", "Time = %f", d=(/TimeA/) )

    call DCClockPredict( &                       ! �Ĥ���֤�ͽ¬
      &   clk = clock_init + clock_loop,       & ! (in)
      &   progress = progress                  & ! (in) 
      &  )

    ! CPU time measurement starts again.
    call ClocksetLoopStart

  end subroutine ClocksetPredict


  subroutine ClocksetClose
    !
    ! CPU ���ַ�¬�ν�λ����
    ! ��¬��λ��, CPU ���֤�ɽ������.
    !

    !�⥸�塼���ɤ߹���
    use dc_clock,    only : DCClockClose, DCClockResult

    !���ۤη�����ػ�
    implicit none

    call DCClockResult( &                    ! �� CPU ���֤�ɽ��
      & clks = (/clock_init, clock_loop/), & ! (in)
      & total_auto = .true. )                ! (in)

    call DCClockClose( clk = clock_init )      ! (inout) ! Finalize (�����)
    call DCClockClose( clk = clock_loop )      ! (inout) ! Finalize (�����)
    
  end subroutine ClocksetClose
  
end module clockset
