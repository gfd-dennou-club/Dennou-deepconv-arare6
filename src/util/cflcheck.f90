!= CFL ���Υ����å��򤹤뤿��Υѥå��������⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: cflcheck.f90,v 1.3 2014/03/04 05:55:06 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module CFLCheck
  !
  ! CFL ���Υ����å��򤹤뤿��Υѥå��������⥸�塼��
  !  * ���Ȥ��Ф��� CFL ��������å�
  !  * ���Ϥ��줿®�٤��Ф��� CFL ��������å�
  !

  !���ۤη�����ػ�
  implicit none

  !private °���λ���
  private

  !�����
  character(*), parameter :: module_name = 'cflcheck'
                              ! �⥸�塼���̾��.
                              ! Module name
  
  !�ؿ��� public °��������
  public CFLCheckTimeShort
  public CFLCheckTimeLongVelX
  public CFLCheckTimeLongVelY
  public CFLCheckTimeLongVelZ
  
contains  

!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeShort( xyz_VelSound )
    !
    !���Ȥ��Ф��� CFL ��������å�
    ! 

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: FlagCalc3D, &
      &                   imin,        &! x ����������β���
      &                   imax,        &! x ����������ξ��
      &                   jmin,        &! z ����������β���
      &                   jmax,        &! z ����������ξ��
      &                   kmin,        &! z ����������β���
      &                   kmax          ! z ����������ξ��
    use axesset,    only: dx,          &! x �����γʻ����ֳ�
      &                   dy            ! y �����γʻ����ֳ�
    use timeset,    only: DelTimeShort  !û�����֥��ƥå�
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: xyz_VelSound(imin:imax,jmin:jmax, kmin:kmax)
                                        !��®
    real(DP)             :: xyz_CFL(imin:imax,jmin:jmax, kmin:kmax)
                                        !��������
    
    !��®�� CFL �������
    !
    if ( FlagCalc3D ) then 
      xyz_CFL = DelTimeShort * xyz_VelSound       &
        &       * ((1.0d0 / (dx * dx) + 1.0d0 / (dy * dy)) ** 0.5d0)
    else
      xyz_CFL = DelTimeShort * xyz_VelSound / dx
    end if

    !��å�����
    !
    call MessageNotify( "M", &
      & module_name, &
      & "Sound Wave Velocity = %f", d=(/maxval(xyz_VelSound)/) )
    call MessageNotify( "M", &
      & module_name, &
      & "DelTimeShort = %f", d=(/DelTimeShort/) )
    
    !�ٹ��å�����
    if ( maxval(xyz_CFL) >= 1.0) then 
      call MessageNotify( "E", &
        & module_name, &
        & "CFL Condition is broken, DelTimeShort * VelSound > min(DelX, DelZ)")
    else
      call MessageNotify( "M", &
        & module_name, &
        & "Courant number for DelTimeSort = %f", d=(/maxval(xyz_CFL)/) )
    end if

  end subroutine CFLCheckTimeShort


!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelX( pyz_VelX )
    !
    !��ʿ®�٤��Ф��� CFL ��������å�. 
    ! 

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x ����������β���
      &                   imax,      &! x ����������ξ��
      &                   jmin,      &! z ����������β���
      &                   jmax,      &! z ����������ξ��
      &                   kmin,      &! z ����������β���
      &                   kmax        ! z ����������ξ��
    use axesset,    only: dx          ! x �����γʻ����ֳ�
    use timeset,    only: DelTimeLong !Ĺ�����֥��ƥå�

    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP), intent(in) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL

    !CFL �������
    CFL = ( 2.0d0 * DelTimeLong ) * maxval( abs( pyz_VelX / dx ) )
  
    !��å���������
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelX for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelX)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelX
    

!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelY( xqz_VelY )
    !
    !��ʿ®�٤��Ф��� CFL ��������å�. 
    ! 

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x ����������β���
      &                   imax,      &! x ����������ξ��
      &                   jmin,      &! z ����������β���
      &                   jmax,      &! z ����������ξ��
      &                   kmin,      &! z ����������β���
      &                   kmax        ! z ����������ξ��
    use axesset,    only: dy          ! y �����γʻ����ֳ�
    use timeset,    only: DelTimeLong !Ĺ�����֥��ƥå�

    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP), intent(in) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL

    !CFL �������
    CFL = (2.0d0 * DelTimeLong) * maxval( abs( xqz_VelY / dy ) )
  
    !��å���������
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelY for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelY)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelY


!!!-----------------------------------------------------------------!!!
  subroutine CFLCheckTimeLongVelZ( xyr_VelZ )
    !
    !��ʿ®�٤��Ф��� CFL ��������å�. 
    ! 

    !�⥸�塼���ɤ߹���
    use dc_types,   only: DP
    use dc_message, only: MessageNotify
    use gridset,    only: imin,      &! x ����������β���
      &                   imax,      &! x ����������ξ��
      &                   jmin,      &! z ����������β���
      &                   jmax,      &! z ����������ξ��
      &                   kmin,      &! z ����������β���
      &                   kmax        ! z ����������ξ��
    use axesset,    only: dz          ! z �����γʻ����ֳ�
    use timeset,    only: DelTimeLong !Ĺ�����֥��ƥå�

    !���ۤη�����ػ�
    implicit none
    
    !�����ѿ�
    real(DP), intent(in) :: xyr_VelZ(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)             :: CFL
    
    !CFL �������
    CFL = ( 2.0d0 * DelTimeLong ) * maxval( abs( xyr_VelZ / dz ) ) 
    
    !��å���������
    call MessageNotify( "M", &
      & module_name, &
      & "Courant number of VelZ for DelTimeLong = %f", d=(/CFL/) )

    if (CFL > 1.0d0) then 
      call MessageNotify( "E", "CFLCheck (VelZ)", "CFL > 1.0" )
    end if

  end subroutine CFLCheckTimeLongVelZ
  
end module CFLCheck
