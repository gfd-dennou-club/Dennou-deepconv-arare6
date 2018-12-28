!= ��ñ����: ������ǥ��󥰥ե����뤫���ǮΨ��Ϳ����
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_sounding.f90,v 1.2 2014/03/04 04:49:42 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_Sounding
  !
  ! ��ñ����: ������ǥ��󥰥ե����뤫���ǮΨ��Ϳ����

  ! �⥸�塼���ɤ߹���
  !
  use dc_types, only: DP
  
  ! ���ۤη�����ػ�
  !
  implicit none

  ! private °��
  !
  private

  !�ѿ����
  !
  real(DP), save, allocatable, public :: xyz_DPTempDtRadSndg(:,:,:)   
                                           !���Ͳ�Ǯ��
  real(DP), save, allocatable, public :: xyz_ExnerRadSndg(:,:,:)  
                                           !���Ͳ�Ǯ��
  real(DP), save :: FactorDExnerDtRad = 1.0d0

  public  radiation_sounding_init
  public  radiation_sounding_forcing

contains

!!!----------------------------------------------------------------------!!!
  subroutine radiation_sounding_init
    !
    ! ������ǥ��󥰥ե����뤫���ǮΨ���ɤ߹���
    !
    
    ! �⥸�塼���ɤ߹���
    !
    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
!    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax        !z ����������ξ��
    use axesset,           only : z_Z, r_Z    !Z ��ɸ��
    use constants,         only : DayTime     ! 1 ����Ĺ�� [s]
    use basicset,          only : xyz_ExnerBZ !�������ʡ��ؿ��δ��ܾ�
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use namelist_util,     only : namelist_filename
    
    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    integer             :: AltCol = 0       !�ֹ��١פ����ֹ� (������ǥ��󥰥ե�������)
    integer             :: SWaveCol = 0     !��û�����ͤˤ���ǮΨ�פ����ֹ� (������ǥ��󥰥ե�������)
    integer             :: LWaveCol = 0     !��Ĺ�����ͤˤ���ǮΨ�פ����ֹ� (������ǥ��󥰥ե�������)
    integer             :: unit             !����ե������������ֹ�  
    integer             :: k
    integer             :: io
    character(30)       :: SoundingFile     ! ������ǥ��󥰥ե�����
    
    integer, parameter  :: maxch=12
    character(len=200)  :: buf, eachcol(maxch)
    integer             :: i, num, MaxCol

    real(DP)            :: r_tmpQradSW(10000)
    real(DP)            :: r_tmpQradLW(10000)
    real(DP)            :: r_tmpAlt(10000)
    real(DP)            :: r_tmpQrad(10000)
    integer             :: NumRec = 0

    real(DP), allocatable :: xyr_DPTempDtRadSndg(:,:,:)   

    logical             :: flag

   
    !����ե����뤫���ɤ߹�����ϥե��������
    !
    NAMELIST /radiation_sounding_nml/ &
      & SoundingFile, AltCol, SWaveCol, LWaveCol, FactorDExnerDtRad
              
    !����ե����뤫����ϥե�����˵��ܤ��������ɤ߹���
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_sounding_nml)
    close(unit)
      
    ! �����
    !
    io = 0
    NumRec = 0
    MaxCol = max( AltCol, max( SWaveCol, LWaveCol ) )
!    write(*,*) "MaxCol", MaxCol
    r_tmpAlt  = 0.0d0; r_tmpQrad = 0.0d0; r_tmpQradSW = 0.0d0; r_tmpQradLW = 0.0d0

    ! �ե�����Υ����ץ�
    !
    open (17, file=SoundingFile, status='old')
    
    ! �ե�����ƤӽФ�
    !
    do while ( io == 0 ) 
      ! 1 ��ʬ�ɤ߽Ф�
      !
      read (17, '(a)', IOSTAT=io) buf
      
      ! �Ԥ򥫥�޶��ڤ��ʬ��
      !
      call devidecsv( buf, eachcol, maxch, num )
      
      ! ��ǧ
      !
!      write(*,*) num
!      do i=1, num
!        write(*,*) i, eachcol(i)(1:len_trim(eachcol(i)))
!      end do
      
      ! num ���ͤ���������Τϥإå��Ȥߤʤ�. 
      !
      if (num >= MaxCol) then 
        ! �Կ��η׻�
        !
        NumRec = NumRec + 1        
        
        ! �ͤ�����
        !
        if (AltCol > 0)   read( eachcol(AltCol)(1:len_trim(eachcol(AltCol))), *)   r_tmpAlt(NumRec) 
        if (SWaveCol > 0) read( eachcol(SWaveCol)(1:len_trim(eachcol(SWaveCol))), *) r_tmpQradSW(NumRec) 
        if (LWaveCol > 0) read( eachcol(LWaveCol)(1:len_trim(eachcol(LWaveCol))), *) r_tmpQradLW(NumRec) 

        r_tmpQrad(NumRec) = r_tmpQradSW(NumRec) + r_tmpQradLW(NumRec)

      end if

!      write(*,*) r_tmpAlt(NumRec), r_tmpQradSW(NumRec), r_tmpQradLW(NumRec), r_tmpQrad(NumRec) 

    end do

    ! �ե�����Υ�����
    !
    close (17)

    !�����
    !
    allocate( xyz_DPTempDtRadSndg(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DPTempDtRadSndg(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadSndg(imin:imax, jmin:jmax, kmin:kmax) )

    xyz_DPTempDtRadSndg = 0.0d0
    xyr_DPTempDtRadSndg = 0.0d0
    xyz_ExnerRadSndg = 0.0d0

    flag = .false.
    
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_Z(k) == r_tmpAlt(i) ) then           
          xyr_DPTempDtRadSndg(:,:,k) = r_tmpQrad(i)         
        end if
      end do
    end do

    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
!          write(*,*) r_tmpAlt(i), r_tmpQrad(i), '(', r_tmpQradSW(i), '+', r_tmpQradLW(i), ')'
          xyz_DPTempDtRadSndg(:,:,k) = r_tmpQrad(i)
        end if
      end do
    end do

    if (.NOT. flag) then 
      do k = kmin+1, kmax
        xyz_DPTempDtRadSndg(:,:,k) = ( xyr_DPTempDtRadSndg(:,:,k-1) + xyr_DPTempDtRadSndg(:,:,k) ) * 5.0d-1
      end do
    end if


    ! ñ�̴���
    !
    xyz_DPTempDtRadSndg = xyz_DPTempDtRadSndg / DayTime / xyz_ExnerBZ

    ! ���Ϥ� tendency
    !
    xyz_ExnerRadSndg = xyz_DExnerDt_xyz(xyz_DPTempDtRadSndg) * FactorDExnerDtRad

    ! �ҥ��ȥ�ǡ������
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='float')

  end subroutine radiation_sounding_init


  subroutine Radiation_Sounding_forcing(   &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! �⥸�塼���ɤ߹���
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x ����������β���
      &                           imax,     & !x ����������ξ��
      &                           jmin,     & !y ����������β���
      &                           jmax,     & !y ����������ξ��
      &                           kmin,     & !z ����������β���
      &                           kmax,     & !z ����������ξ��
      &                           nx, ny, nz
    
    ! ���ۤη�����ػ�
    !
    implicit none

    ! �������ѿ�
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency �ι���
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadSndg
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadSndg
    
    ! �ե��������
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadSndg(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadSndg(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_Sounding_forcing


  subroutine devidecsv( buf, eachcol, maxch, num )
    !
    ! ������ǥ��󥰥ե�������ɤ߹��ि��Υ롼����

    implicit none

    integer          :: maxch, num
    character(len=*) :: buf, eachcol(maxch)
    integer          :: i, j, prev, now, now1, ll
    logical          :: quote
    
    quote = .FALSE.
    
    ll = len_trim(buf)
    prev=0; now=1; i=1
    do j=1, ll
      if( buf(j:j) == '"' ) then
        quote = .NOT. quote
      end if
      if( quote .eqv. .FALSE. ) then
        if(buf(j:j) == ',') then
          now = j
          now1 = now -1
          if( now1 < prev+1 ) then
            eachcol(i) = ''
          else
            eachcol(i) = buf(prev+1:now1)
          end if
          i=i+1
          if( i > maxch ) then
            write(0,*) 'maxch is too small'
            stop
          end if
          prev=now
        end if
      end if
    end do
    
    if( prev < ll ) then
      eachcol(i) =  buf(prev+1:ll)
      num = i
    else
      num=i-1
    end if
    
  end subroutine devidecsv
  
end module Radiation_Sounding
