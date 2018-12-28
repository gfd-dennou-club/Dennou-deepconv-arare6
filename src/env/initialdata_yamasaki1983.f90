!= Module initialdata_yamasaki1983
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_yamasaki1983.f90,v 1.8 2014/03/04 04:49:40 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module initialdata_yamasaki1983
  !
  ! Yamasaki (1983) �δ��ܾ졦���������ꤹ��

   
  !���ۤη�����ػ�
  implicit none

  !�������������
  public initialdata_yamasaki1983_basic

contains

!!------------------------------------------------------------------------------!!!
  subroutine initialdata_yamasaki1983_basic( z_Temp, z_Press, zf_MolFr )
    !
    !== ����
    !  * Yamasaki, 1983 �β��٤����м��٤δ�¬�ͤǴ��ܾ���������
    !    * ���٤δ��ܾ�
    !      * ��¬�ǡ�����NetCDF�ե����벽������Τ����ͤ��ɤ߹���
    !        * �ɤ߹�����ͤ�������֤��Ʋ��٤δ��ܾ���������
    !    * ���٤δ��ܾ�
    !      * ���м��٤� subroutine HUM �Ǻ����Ѥ�
    !        * �����Ǥ����м��٤�������Ѵ����Ƽ��٤δ��ܾ��������� 
    !    * �����δ��ܾ�
    !      * ������ʬ�ȴ�����ʬ��ʬ���̺����θ�����ſ尵ʿ�դ���׻�����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only: STRING, DP
    use dc_message, only: MessageNotify
    use gridset,    only: kmin, kmax,  &!���󥵥��� (Z ����)
      &                   nz,          &!ʪ���ΰ�Υ����� (Z ����)
      &                   ncmax         !�Ž���ʬ�ο�
    use axesset,    only: z_Z,         &!�����顼�ʻ����Ǥι���
      &                   dz            !�ʻҴֳ�
    use constants,  only: GasRDry,       &!������ʬ���갵��Ǯ
      &                   TempSfc,       &!��ɽ�̲���
      &                   PressSfc,      &!��ɽ�̰���
      &                   PressBasis,    &!��ɽ�̰���
      &                   Grav            !����
    use chemcalc,   only: SvapPress       !˰�¾�����
    !  use composition

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    integer, parameter :: mmin = 1
    integer, parameter :: mmax1 = 37
    integer, parameter :: mmax2 = 31
    integer :: ID = 6    ! ���ʪ���ֹ�

    real(DP), intent(out) :: z_Press(kmin:kmax)!����
    real(DP), intent(out) :: z_Temp(kmin:kmax) !����
    real(DP), intent(out) :: zf_MolFr(kmin:kmax, 1:ncmax) !�����
    real(DP)              :: Ob_Alt1(mmin:mmax1)    !��¬�ǡ������б��������
    real(DP)              :: Ob_Alt2(mmin:mmax2)    !��¬�ǡ������б��������
    real(DP)              :: Ob_Temp(mmin:mmax1)    !NetCDF �����ɤ߹���������δ�¬��(1����)
    real(DP)              :: Ob_Hum(mmin:mmax2)     !���м��٤δ�¬��(1����)

    real(DP)              :: z_Hum(kmin:kmax)        !���м��٤δ��ܾ����(%)
    integer               :: k, m

    ! �����
    !
    z_Press  = 0.0d0
    z_Temp   = 0.0d0
    zf_MolFr = 0.0d0

    ! ��ǧ. ��ɽ�̲��١����Ϥ���ꤹ��Τ��̥⥸�塼��ʤΤ�. 
    !
    if (     PressBasis /= 1010.0d2  &
      & .OR. PressSfc   /= 1010.0d2  &
      & .OR. TempSfc    /= 302.0d0  ) then 
      
      call MessageNotify( "E", "initaldata_yamasaki1983_init", &
        & "Constants are wrong. please PressSfc = 1.01d5, TempSfc = 302.0d0")
    end if

    ! ��¬�ǡ���. 
    !
    Ob_Alt1 = (/ &
      & 0.05d3, 0.16d3, 0.29d3, 0.44d3, 0.61d3,      &
      & 0.80d3, 1.02d3, 1.28d3, 1.58d3, 1.95d3,      &
      & 2.37d3, 2.82d3, 3.33d3, 3.90d3, 4.50d3,      &
      & 5.10d3, 5.70d3, 6.30d3, 6.90d3, 7.50d3,      &
      & 8.10d3, 8.70d3, 9.30d3, 9.90d3, 10.50d3,     &
      & 11.10d3, 11.70d3, 12.35d3, 13.05d3, 13.80d3, &
      & 14.65d3, 15.50d3, 16.55d3, 17.70d3, 18.95d3, &
      & 20.30d3, 21.80d3                             &
      &/)

    Ob_Temp = (/ &
      & 299.60d0, 298.72d0, 297.68d0, 296.48d0, 295.13d0, &
      & 293.90d0, 292.47d0, 290.90d0, 289.40d0, 287.55d0, &
      & 285.45d0, 283.13d0, 280.38d0, 277.25d0, 273.95d0, &
      & 270.65d0, 267.35d0, 263.90d0, 260.30d0, 256.55d0, &
      & 252.50d0, 248.15d0, 243.50d0, 238.70d0, 233.90d0, &
      & 229.10d0, 224.30d0, 219.25d0, 213.85d0, 208.50d0, &
      & 203.70d0, 200.55d0, 199.60d0, 201.40d0, 205.15d0, &
      & 209.20d0, 212.90d0                                &
      & /)

    Ob_Alt2 = (/ &
      & 0.00d3, 0.60d3, 1.20d3, 1.80d3, 2.40d3, &
      & 3.00d3, 3.60d3, 4.20d3, 4.80d3, 5.40d3, &
      & 6.00d3, 6.60d3, 7.20d3, 7.80d3, 8.40d3, &
      & 9.00d3, 9.60d3, 10.2d3, 10.8d3, 11.4d3, &
      & 12.0d3, 12.7d3, 13.4d3, 14.2d3, 15.1d3, &
      & 16.0d3, 17.1d3, 18.3d3, 19.6d3, 21.0d3, 22.6d3 &
      & /)

    Ob_Hum = (/ &
      & 83.0d-2, 91.0d-2, 95.0d-2, 95.0d-2, 91.0d-2, &
      & 85.0d-2, 80.0d-2, 75.0d-2, 71.0d-2, 68.0d-2, &
      & 66.0d-2, 65.0d-2, 64.0d-2, 63.0d-2, 62.0d-2, &
      & 61.0d-2, 61.0d-2, 62.0d-2, 62.0d-2, 63.0d-2, &
      & 63.0d-2, 64.0d-2, 67.0d-2, 74.0d-2, 75.0d-2, &
      & 46.0d-2, 26.0d-2,  7.0d-2,  1.0d-2,  0.5d-2,  0.2d-2 &
      & /)

    ! �����䴰
    !
    ! Ob_Hum ��������֤��� z_Hum ���������
    ! Ob_Alt(m) < z_Z(k) < Ob_Alt(m+1),
    ! Ob_Alt = z_Z(k),
    ! Ob_Alt < z_Z(k)
    ! ��3�ĤǾ��ʬ��
    !
    do k = 1, nz
      do m = mmin, mmax1 - 1
        ! z_Z �� Ob_altitude �Τ����֤˶��ޤ��Ȥ�,
        ! ���ζ�֤� Obaltitude(m), Obaltitude(m+1) ����
        ! ľ���� Ob_TempZ ��������֤���
        if (Ob_Alt1(m) /=  z_Z(k) .AND. z_Z(k) > Ob_Alt1(m) .AND. z_Z(k) < Ob_Alt1(m+1)) then

          z_Temp(k) = Ob_Temp(m) &
            & + ( ( Ob_Temp(m+1) - Ob_Temp(m) ) / ( Ob_Alt1(m+1) - Ob_Alt1(m) ) ) &
            &   * ( z_Z(k) - Ob_Alt1(m) )
          
        else if (Ob_Alt1(m) == z_Z(k)) then
          z_Temp(k)  = Ob_Temp(m)

        ! z_Z(k) > Ob_altitude �Ǥϴ�¬�ǡ�����̵���Τ������絤�ˤ���          
        else if (Ob_Alt1(m) < z_Z(k)) then
          z_Temp(k)  = z_Temp(k-1)
          
        end if
      end do
    end do

    do k = 1, nz
      do m = mmin, mmax2 - 1 
        if (Ob_Alt2(m) /=  z_Z(k) .AND. z_Z(k) > Ob_Alt2(m) .AND. z_Z(k) < Ob_Alt2(m+1)) then

          z_Hum(k) = Ob_Hum(m) &
            & + ( ( Ob_Hum(m+1) - Ob_Hum(m) ) / ( Ob_Alt2(m+1) - Ob_Alt2(m) ) ) &
            &   * (z_Z(k) - Ob_Alt2(m))
          
        else if (Ob_Alt2(m) == z_Z(k)) then
          z_Hum(k) = Ob_Hum(m)

        ! z_Z(k) > Ob_altitude �Ǥϴ�¬�ǡ�����̵���Τ����м��ٰ���ˤ��� 
        else if (Ob_Alt2(m) < z_Z(k)) then
          z_Hum(k) = z_Hum(k-1)

        end if
      end do
    end do

    ! �����
    z_Press = 1.0d-60

    ! ��ɽ�̵����Ȳ��٤����絤�ǲ��ؤε��������
    ! �ſ尵�μ� dP/dz = - \rho * g �����
    z_Press(1) = PressSfc - (Grav * PressSfc * dz * 5.0d-1) / (GasRDry * TempSfc)

    ! �絤�ǲ��ؤΥ�����׻�����
    zf_MolFr(1,1) = SvapPress( ID, z_Temp(1)) * z_Hum(1) / z_Press(1)
    
    ! ���٤δ��ܾ�η׻�
    do k = 2, nz
      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * dz) &
        & / ( GasRDry * z_Temp(k-1) )
      
      zf_MolFr(k,1) = SvapPress( ID, z_Temp(k)) * z_Hum(k) / z_Press(k) 
    end do
    
  end subroutine INITIALDATA_YAMASAKI1983_basic
  
end module Initialdata_yamasaki1983
