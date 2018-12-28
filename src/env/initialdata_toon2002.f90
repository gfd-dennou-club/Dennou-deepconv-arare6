!= Module initialdata_Toon2002
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_toon2002.f90,v 1.3 2014/03/04 04:49:40 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module initialdata_Toon2002
  !
  ! Toon et al. (2002) ���Ϥ���������Ϳ����
  !

  !�⥸�塼���ɤ߹���
  use dc_types,  only: DP

  !���ۤη�����ػ�
  implicit none

  !�����ѿ�
  real(DP), parameter, private :: AntA    = 27.4d0 
  real(DP), parameter, private :: AntB    = 3103.0d0
  real(DP), parameter, private :: TempLTP = 135.0d0
  real(DP), parameter, private :: Dhight  = 6.0d2

  !�������������
  public  initialdata_Toon2002_basic

contains

!!!------------------------------------------------------------------------------!!!

  subroutine initialdata_toon2002_basic( z_Temp, z_Press )
    
    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use gridset,   only: kmin, kmax,  &!���󥵥��� (Z ����)
      &                  nz            !ʪ���ΰ�Υ����� (Z ����)
    use axesset,   only: z_Z,         &!�����顼�ʻ����Ǥι���
      &                  z_dz          !Z �����γʻ����ֳ�
    use constants, only: GasRDry,     &!������ʬ���갵��Ǯ
      &                  CpDry,       &!������ʬ���갵��Ǯ
      &                  Grav,        &!���ϲ�®��
      &                  TempSfc,     &!��ɽ�̲���
      &                  PressSfc      !��ɽ�̰���

    implicit none
    
    real(DP), intent(out):: z_Press(kmin:kmax)           !����
    real(DP), intent(out):: z_Temp(kmin:kmax)            !����
    real(DP)             :: TempLCL
    real(DP)             :: Temp_0,  Temp_1
    real(DP)             :: Press_0, Press_1
    real(DP)             :: Weight1, Weight2
    real(DP)             :: LCL, LTP
    integer              :: k
    

    ! ������Ǯ��, ������Ǯ��, ��������������٤�׻���,
    ! ���ΰ������Ω�ļ����Ѥ��Ʋ���, ���Ϥ�׻�
    ! ������Ǯ���ȼ�����Ǯ�����������(LCL)��ȿ��ˡ�Ƿ׻�
    !
    Press_0 = PressSfc
    Temp_0 = TempSfc
    do
      ! ˰�²��� (press0 ���Ф���): �ǽ����ɽ�̤Ǥ�˰�²���.
      ! ln(p) = A - B/T
      Temp_1 = AntB / (AntA - dlog(Press_0))
      
      ! ������ǮŪ�˷�᤿���� 
      !
      Press_1 = PressSfc * (Temp_1/TempSfc) **(CpDry / GasRDry)

      ! ���и������ͤ�꾮�����ʤ�н�λ. 
      !
      if (abs(Temp_1 - Temp_0) < epsilon(0.0d0)) then
        LCL = (TempSfc * CpDry) / Grav &
          & * (1.0d0 - (Press_1 / PressSfc)**(GasRDry / CpDry))
        TempLCL = temp_1

        exit
      else
        Temp_0 = Temp_1
        Press_0 = Press_1
      end if
    end do


    ! ������Ǯ�������������������(LTP)��׻�
    !
    LTP = LCL + GasRDry * AntB / Grav * dlog(TempLCL / TempLTP)
  
    ! ���ٰ��Ϥ����.
    !
    z_Temp(1)  = TempSfc  - Grav * z_Z(1) / CpDry 
    z_Press(1) = PressSfc - (Grav * PressSfc * z_dz(1) * 5.0d-1) / (GasRDry * TempSfc)
    do k = 2, nz

      !�ŤߤĤ��δؿ����Ѱ�. tanh ���Ѥ���
      Weight1 = ( tanh( (z_Z(k) - LCL ) / Dhight ) + 1.0d0 ) * 5.0d-1

      !�ŤߤĤ��δؿ����Ѱ�. tanh ���Ѥ���
      Weight2 = ( tanh( (z_Z(k) - LTP ) / Dhight ) + 1.0d0 ) * 5.0d-1

      !������Ǯ
      if (z_z(k) < LCL) then 
        z_Temp(k) = TempSfc - Grav * z_Z(k) / CpDry 

      !������Ǯ
      elseif (z_z(k) >= LCL .AND. z_z(k) < LTP) then 
        Temp_0 = TempSfc - Grav * z_Z(k) / CpDry 
        Temp_1 = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))
        z_Temp(k) = Temp_0 * ( 1.0d0 - Weight1 ) + Temp_1 * Weight1
!        z_Temp(k)  = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))

      !����
      elseif (z_z(k) >= LTP) then 
        Temp_0 = TempLCL * exp(-Grav * (z_Z(k) - LCL) / (GasRDry * AntB))
        z_Temp(k) = Temp_0 * ( 1.0d0 - Weight2 ) + TempLTP * Weight2
!        z_Temp(k) = TempLTP

      end if
    end do

    ! �ſ尵ʿ�դ��鰵�Ϥ����
    !
    do k = 2, nz
      z_Press(k) = z_Press(k-1) - (Grav * z_Press(k-1) * z_dz(k-1)) &
        & / ( GasRDry * z_Temp(k-1) )
    end do


  end subroutine initialdata_toon2002_basic
  
end module initialdata_Toon2002
