! ������ǥ��󥰥ե����뤫���������뤿��Υ��֥롼����
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_sounding.f90,v 1.8 2014/07/08 00:59:09 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module initialdata_sounding
  !
  ! ������ǥ��󥰥ե����뤫���������뤿��Υ��֥롼����
  ! ������ǥ��󥰥ե�����ϥƥ����ȥե�����ǽ񤫤�Ƥ���. 
  ! ����Ū�ˤ� netCDF ���ѹ�����ͽ��.
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only: DP
  
  !���ۤη�����ػ�
  implicit none

  !�����ѿ�
  real(DP), save, private :: r_tmpAlt(10000)
  real(DP), save, private :: r_tmpTemp(10000)
  real(DP), save, private :: r_tmpPress(10000)
  real(DP), save, private :: r_tmpPTemp(10000)
  real(DP), save, private :: r_tmpVelX(10000)
  real(DP), save, private :: r_tmpVelY(10000)
  real(DP), save, private :: r_tmpQrad(10000)
  real(DP), save, private :: TempTr = 0.0d0
  real(DP), save, private :: AltTr  = 0.0d0
  real(DP), save, private :: DelAlt = 4.0d3
  integer,  save, private :: NumRec = 0

  !�����ѿ�
  public  initialdata_sounding_init
  public  initialdata_sounding_basic
  public  initialdata_sounding_wind

contains

!!!------------------------------------------------------------------------------!!!
  subroutine initialdata_sounding_init
    !
    !�ե����뤫�饵����ǥ��󥰥ǡ������ɤ߹���
    !

    !�⥸�塼���ɤ߹���
    use dc_types,      only: DP
    use dc_iounit,     only: FileOpen      
    use namelist_util, only: namelist_filename
    
    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    integer             :: AltCol = 0
    integer             :: TempCol = 0
    integer             :: PressCol = 0
    integer             :: VelXCol = 0
    integer             :: VelYCol = 0
    integer             :: unit         !����ե������������ֹ�  
    integer             :: io
    character(30)       :: SoundingFile    
    
    integer, parameter  :: maxch=12
    character(len=100)  :: buf, eachcol(maxch)
    integer             :: num, MaxCol
    
    !����ե����뤫���ɤ߹�����ϥե��������
    !
    NAMELIST /initialdata_sounding_nml/ SoundingFile, AltCol, TempCol, PressCol, VelXCol, VelYCol, TempTr, AltTr, DelAlt    

    !����ե����뤫���ɤ߹���
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_sounding_nml)
    close(unit)

    ! �����
    !
    io = 0
    NumRec = 0
    MaxCol = max( AltCol, max( TempCol, PressCol ) )
    r_tmpAlt  = 0.0d0; r_tmpTemp = 0.0d0; r_tmpPress = 0.0d0; 
    r_tmpVelX = 0.0d0; r_tmpVelY = 0.0d0

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
      
      ! num ���ͤ���������Τϥإå��Ȥߤʤ�. 
      !
      if (num >= MaxCol) then 
        ! �Կ��η׻�
        !
        NumRec = NumRec + 1        

        ! �ͤ�����
        !
        if (AltCol > 0)   read( eachcol(AltCol)(1:len_trim(eachcol(AltCol))), *)    r_tmpAlt(NumRec) 
        if (TempCol > 0)  read( eachcol(TempCol)(1:len_trim(eachcol(TempCol))), *)   r_tmpTemp(NumRec)   
        if (PressCol > 0) read( eachcol(PressCol)(1:len_trim(eachcol(PressCol))), *) r_tmpPress(NumRec) 
        if (VelXCol > 0)  read( eachcol(VelXCol)(1:len_trim(eachcol(VelXCol))), *)   r_tmpVelX(NumRec)   
        if (VelYCol > 0)  read( eachcol(VelYCol)(1:len_trim(eachcol(VelYCol))), *)   r_tmpVelY(NumRec) 

      end if

    end do

    ! �ե�����Υ�����
    !
    close (17)    

  end subroutine initialdata_sounding_init


!!!------------------------------------------------------------------------------!!!
  subroutine  initialdata_sounding_basic( z_Temp, z_Press )
    !
    ! ������ǥ��󥰥ǡ���������ܾ���������
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use gridset,   only: kmin, kmax,       &!���󥵥��� (Z ����)
      &                  nz                 !ʪ���ΰ���礭�� (Z ����)
    use axesset,   only: r_Z, z_Z,         &!����
      &                  dz                 !�ʻҴֳ� (Z ����)
    use constants, only: Grav,             &
      &                  GasRDry
    
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(out):: z_Press(kmin:kmax)           !����
    real(DP), intent(out):: z_Temp(kmin:kmax)            !����
    real(DP)             :: r_Press(kmin:kmax)           !����
    real(DP)             :: r_Temp(kmin:kmax)            !����
    real(DP)             :: z_DTempDZ(kmin:kmax)
    real(DP)             :: DTempDZ
    real(DP)             :: ratio
    integer              :: i, k, k1, k2
    logical              :: flag

    ! �����
    !
    z_Temp  = 0.0d0
    r_Temp  = 0.0d0
    z_Press = 0.0d0
    r_Press = 0.0d0

    flag = .false. 

    ! �����䴰
    !
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_tmpAlt(i) <= r_Z(k) .AND. r_Z(k) < r_tmpAlt(i+1) ) then 
          ratio = ( r_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          r_Temp(k)  = r_tmpTemp(i)  + ( r_tmpTemp(i+1) - r_tmpTemp(i) )   * ratio
          r_Press(k) = r_tmpPress(i) + ( r_tmpPress(i+1) - r_tmpPress(i) ) * ratio
        end if
      end do
    end do

    ! �ǡ����ɤ߹���  
    ! ���ޤ�Ⱦ�ʻҤ˹礦���⤢��Τ�, ���ξ��ϥǡ�����ͥ��. 
    !    
    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
          z_Temp(k) = r_tmpTemp(i)
          z_Press(k) = r_tmpPress(i)
        end if
      end do
    end do
    
    ! r => z ���Ѵ�
    ! 
    if (.NOT. flag) then 
      do k = kmin+1, kmax
        z_Temp(k)  = ( r_Temp(k-1)  + r_Temp(k)  ) / 2.0d0
        z_Press(k) = ( r_Press(k-1) + r_Press(k) ) / 2.0d0
      end do
    end if
  
!!!
!!! �����̤�¸�ߤ�����ˤ�, �����Ǥβ��ٸ��ۤ�ˤ䤫�ˤ���.
!!!
    if ( AltTr > z_Z(1) .AND. AltTr < z_Z(nz) ) then

      ! ���ߤβ���ʬ�ۤθ��ۤ���
      !
      do k = 1, kmax
        z_DTempDZ(k) = (z_Temp(k) - z_Temp(k-1)) / dz
      end do
      
      ! ��ή�����̤���ΰ���. ���ꤵ�줿���٤β��ٸ�Ψ��Ȥ�³����. 
      !
      k1 = minloc( z_Z, 1, z_Z > AltTr ) 
      
      ! ���ٸ�Ψ
      !
      DTempDZ = z_DTempDZ(k1-1)
      k2 = int( DelAlt / dz )    !���������Ѵ�
      
      do k = k1, kmax
        
        ! ���ٸ�Ψ�򥼥�˶�Ť���. 
        !
        DTempDZ = min( -1.0d-14, DTempDZ - DTempDZ / k2 * ( k - k1 ) )
        
        !���ܾ�β��٤����
        z_Temp(k) = z_Temp(k-1) + DTempDZ * dz
        
        !���Ϥ��ſ尵ʿ�դ���׻�
        z_Press(k) =                                      &
          &  z_Press(k-1) * ( ( z_Temp(k-1) / z_Temp(k) ) &
          &    ** (Grav / ( DTempDZ * GasRDry ) ) )
        
      end do

    end if

  end subroutine Initialdata_sounding_basic


  subroutine initialdata_sounding_wind(pyz_VelX, xqz_VelY)

    !�⥸�塼���ɤ߹���
    use dc_types, only: DP
    use gridset,  only: imin, imax,       &!���󥵥��� (X ����)
      &                 jmin, jmax,       &!���󥵥��� (Y ����)
      &                 kmin, kmax         !���󥵥��� (Z ����)
    use axesset,  only: r_Z, z_Z           !����
      
    implicit none

    real(DP), intent(out) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: pyr_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: xqr_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: ratio
    integer               :: i, k
    logical               :: flag
    
    !�����
    pyz_VelX = 0.0d0
    pyr_VelX = 0.0d0
    xqz_VelY = 0.0d0
    xqr_VelY = 0.0d0

    flag = .false. 

    ! �ǡ����ɤ߹���  
    !        
!    do k = kmin, kmax
!      do i = 1, NumRec
!        if ( r_Z(k) == r_tmpAlt(i) ) then 
!          pyr_VelX(:,:,k) = r_tmpVelX(i)
!          xqr_VelY(:,:,k) = r_tmpVelY(i)
!        end if
!      end do
!    end do

    ! �����䴰
    !
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_tmpAlt(i) <= r_Z(k) .AND. r_Z(k) < r_tmpAlt(i+1) ) then 
          ratio = ( r_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          pyr_VelX(:,:,k) = r_tmpVelX(i) + ( r_tmpVelX(i+1) - r_tmpVelX(i) ) * ratio
          xqr_VelY(:,:,k) = r_tmpVelY(i) + ( r_tmpVelY(i+1) - r_tmpVelY(i) ) * ratio
        end if
      end do
    end do


    ! �ǡ����ɤ߹���  
    ! ���ޤ�Ⱦ�ʻҤ˹礦���⤢��Τ�, ���ξ��ϥǡ�����ͥ��. 
    !        
    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
          pyz_VelX(:,:,k) = r_tmpVelX(i)
          xqz_VelY(:,:,k) = r_tmpVelY(i)
        end if
      end do
    end do

    ! r => z ���Ѵ�
    ! 
    if (.NOT. flag) then 
      do k = kmin+1, kmax
        pyz_VelX(:,:,k) = ( pyr_VelX(:,:,k-1) + pyr_VelX(:,:,k) ) * 5.0d-1
        xqz_VelY(:,:,k) = ( xqr_VelY(:,:,k-1) + xqr_VelY(:,:,k) ) * 5.0d-1
      end do
    end if

  end subroutine initialdata_sounding_wind
  

  subroutine devidecsv( buf, eachcol, maxch, num )
    !
    ! csv �����Υե����뤫��ɬ�פʾ�����ɤ߽Ф�
    !

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
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

end module initialdata_sounding
