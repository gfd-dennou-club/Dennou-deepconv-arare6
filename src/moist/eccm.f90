!= Module ECCM
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: eccm.f90,v 1.11 2014/03/04 05:55:05 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module ECCM
  !
  !��ǮŪ�˾徺���뵤���β��ٸ�Ψ��׻���, �ſ尵ʿ�դ��鰵�Ϥ����
  !

  !���ۤη�����ػ�
  implicit none

  !�ؿ��θ���
  public ECCM_MolFr
  public ECCM_Stab
  public ECCM_Dry
  public ECCM_Wet

contains

!!!----------------------------------------------------------------------!!!
  subroutine ECCM_Dry( a_MolFrIni, Humidity, z_Temp, z_Press, za_MolFr )
    !
    !== ����
    !  * ������Ǯ���ٸ�Ψ�˱�ä����١����Ϥ��ᡢ�������Ф���
    !    ���ꤵ�줿���м��٤Ȥʤ�褦�˶ŷ���ʬ�Υ��������
    !  * ��Ǯ�ϴ��絤���Τ�Τ���ɽ������
    !    * ή�Τ��������ˤ�������Ǯ�ϴ�����ʬ�Τ�Τ���ɽ�����Ƥ��뤿��
    !  * �絤��ʿ��ʬ���̤ˤϼ�����ʬ��ʬ���̤����
    !    * ή�Τ��������ˤ�����, ������ʬ��ʬ���̤Ϲ�θ���Ƥ��뤿��
    !
    
    !�⥸�塼���ɤ߹���
    use dc_types,    only: DP
    use gridset,     only: kmin, kmax,    &!���󥵥��� (Z ����)
      &                    ncmax           !ʪ����
    use axesset,     only: dz              !�ʻҥ�����
    use constants,   only: MolWtDry,      &!
      &                    CpDryMol,      &!
      &                    TempSfc,       &!
      &                    PressSfc,      &!
      &                    Grav            !
    use chemcalc,    only: SvapPress,     &!
      &                    DelMolFrNH4SH 
    use composition, only: CondNum,       &!�ŷ�����ο�
      &                    SpcWetID,      &!
      &                    IdxCG,         &!�ŷ����(����)������ź����
      &                    IdxCC,         &!�ŷ����(��)������ź����
      &                    IdxNH3,        &!NH3(����)������ź����
      &                    IdxH2S          !H2S(����)������ź����
    use ChemData,    only: GasRUniv        
    
    !���ۤη�����ػ�
    implicit none
    
    real(DP), intent(in) :: a_MolFrIni(1:ncmax)     !���������ǤΥ����
    real(DP), intent(in) :: Humidity                !���м��� ( Humidity <= 1.0 )
    real(DP), intent(out):: z_Temp(kmin:kmax)       !����
    real(DP), intent(out):: z_Press(kmin:kmax)      !����
    real(DP), intent(out):: za_MolFr(kmin:kmax, 1:ncmax) 
                                                    !���ʬΨ
    real(DP)             :: SatPress                !˰�¾�����
    real(DP)             :: VapPress                !������
    real(DP)             :: DelMolFr
    integer              :: k, s
    logical              :: a_FlagCond(1:CondNum)
    
    !-------------------------------------------------------------
    ! ����ν����
    !-------------------------------------------------------------
    !�����
    a_FlagCond = .false.

    !��ɽ�̤�ʬ���̤����
    za_MolFr  = 1.0d-60
    za_MolFr(1, 1:ncmax) = a_MolFrIni(1:ncmax) 

    !��ɽ�̤Ǥβ���(nx ��, ���� DelZ / 2 ������)
    z_Temp    = 1.0d-60
    z_Temp(1) = TempSfc - Grav * MolWtDry &
      &               / CpDryMol * ( dz * 5.0d-1 )

    !��ɽ�̤Ǥΰ���(1 ��, ���� DelZ / 2 ������)
    z_Press    = 1.0d-60
    z_Press(1) = PressSfc *((TempSfc / z_Temp(1)) ** (- CpDryMol /  GasRUniv))

    !-----------------------------------------------------------
    ! (1) ������Ǯ���˱�ä����٤����
    ! (2) �ſ尵ʿ�դ��鰵�Ϥ����
    ! (3) (1),(2) �β��ٰ��Ϥ��Ф���, �Ȥ������м��٤Ȥʤ���������
    !-----------------------------------------------------------    
    DtDz: do k = 1, kmax-1
      
      !(1)������Ǯ���˱�ä� k+1 �Ǥβ��٤�׻�
      z_Temp(k+1) = z_Temp(k) - Grav * MolWtDry / CpDryMol * dz
      
      !ǰ��
      if (z_Temp(k+1) <= 0.0d0 ) z_Temp(k+1) = z_Temp(k) 
      
      !(2)���Ϥ��ſ尵ʿ�դ���׻�
      z_Press(k+1) =                                                  &
        &  z_Press(k) * ((z_Temp(k) / z_Temp(k+1)) ** (- CpDryMol / GasRUniv)) 

      !(3)�����η׻�
      !  �ޤ��ϥ������Ѳ����ʤ���ΤȤ��ƥ�����Ϳ����
      !  ˰�¾�������ʿ������Ȥ�ʿ�վ�������Ŭ�Ѥ��Ƥ���
      za_MolFr(k+1,:) = za_MolFr(k,:)
      
      do s = 1, CondNum

        !˰�¾�����
        !
        SatPress = SvapPress( SpcWetID(IdxCC(s)), z_Temp(k+1) )        
        
        !�����Υ��ʬΨ���Ѥ��Ƹ��ߤξ�������׻�
        !
        VapPress = za_MolFr(k,IdxCG(s)) * z_Press(k+1)

        !�ŷ���٤�Ķ�������ݤ��Υ����å�
        !
        if (.NOT. a_FlagCond(s) ) then 
          if ( VapPress > SatPress ) then
            a_FlagCond(s) = .true.
          end if
        end if

        !�ŷ���٤��⤤����, ˰�¾������Ȱ��Ϥ��������׻�
        !
        if ( a_FlagCond(s) ) then 
          za_MolFr(k+1,IdxCG(s)) = max(SatPress * Humidity / z_Press(k+1), 1.0d-16)          
        end if

      end do
      
      !NH4SH ��ʿ�վ��
      if ( IdxNH3 /= 0 ) then 
        DelMolFr =                                              &
          & max (                                               &
          &    DelMolFrNH4SH(                                   &
          &         z_Temp(k+1), z_Press(k+1),                  &
          &         za_MolFr(k+1,IdxNH3), za_MolFr(k+1,IdxH2S), &
          &         Humidity                                    &
          &      ),                                             &
          &    0.0d0                                            &
          &  )
        za_MolFr(k+1,IdxNH3) = za_MolFr(k+1,IdxNH3) - DelMolFr
        za_MolFr(k+1,IdxH2S) = za_MolFr(k+1,IdxH2S) - DelMolFr
      end if
      
    end do DtDz
    
  end subroutine ECCM_Dry


!!-----------------------------------------------------------------------!!!
  subroutine ECCM_Wet( a_MolFrIni, Humidity, z_Temp, z_Press, za_MolFr )
    !
    !== ����
    !  * ������Ǯ���ٸ�Ψ�˱�ä����١����Ϥ��ᡢ�������Ф���
    !    ���ꤵ�줿���м��٤Ȥʤ�褦�˶ŷ���ʬ�Υ��������
    !  * ��Ǯ�ϴ��絤���Τ�Τ���ɽ������
    !    * ή�Τ��������ˤ�������Ǯ�ϴ�����ʬ�Τ�Τ���ɽ�����Ƥ��뤿��
    !  * �絤��ʿ��ʬ���̤ˤϼ�����ʬ��ʬ���̤����
    !    * ή�Τ��������ˤ�����, ������ʬ��ʬ���̤Ϲ�θ���Ƥ��뤿��
    !

    !�⥸�塼���ɤ߹���
    use dc_types,    only: DP
    use gridset,     only: kmin, kmax,    &!���󥵥��� (Z ����)
      &                    ncmax           !ʪ����
    use axesset,     only: dz              !�ʻҥ�����
    use constants,   only: MolWtDry,      &!
      &                    CpDryMol,      &!
      &                    TempSfc,       &!
      &                    PressSfc,      &!
      &                    Grav            !
    use chemcalc,    only: SvapPress,     &!
      &                    DelMolFrNH4SH 
    use composition, only: CondNum,       &!�ŷ�����ο�
      &                    SpcWetID,      &!
      &                    IdxCG,         &!�ŷ����(����)������ź����
      &                    IdxCC,         &!�ŷ����(��)������ź����
      &                    IdxNH3,        &!NH3(����)������ź����
      &                    IdxH2S          !H2S(����)������ź����
    use ChemData,    only: GasRUniv        
      
    !���ۤη�����ػ�
    implicit none
    
    real(DP), intent(in) :: a_MolFrIni(1:ncmax)    !���������ǤΥ����
    real(DP), intent(in) :: Humidity                !���м��� ( Humidity <= 1.0 )
    real(DP), intent(out):: z_Temp(kmin:kmax) !����
    real(DP), intent(out):: z_Press(kmin:kmax)!����
                                                   !ʿ��ʬ����
    real(DP), intent(out):: za_MolFr(kmin:kmax, 1:ncmax) 
                                                   !���ʬΨ
    real(DP)             :: SatPress                !˰�¾�����
    real(DP)             :: VapPress                !������
    real(DP)             :: DelMolFr
    real(DP)             :: a_MolFr(ncmax)         !�����κ������
    integer              :: k, s

    real(DP)             :: Temp1, Press1, DTempDZ1
    real(DP)             :: Temp2, Press2, DTempDZ2
    real(DP)             :: Temp3, Press3, DTempDZ3
    real(DP)             :: Temp4, Press4, DTempDZ4
    real(DP)             :: DTempDZ
    
    !-------------------------------------------------------------
    ! ����ν����
    !-------------------------------------------------------------
    !��ɽ�̤�ʬ���̤����
    za_MolFr  = 1.0d-60
    za_MolFr(1, 1:ncmax)   = a_MolFrIni(1:ncmax) 
    
    !��ɽ�̤Ǥβ���(1 ��, ���� DelZ / 2 ������)
    z_Temp          = 1.0d-60
    z_Temp(1) = TempSfc - Grav * MolWtDry  &
      &               / CpDryMol * ( dz * 5.0d-1 )
    
    !��ɽ�̤Ǥΰ���(1 ��, ���� DelZ / 2 ������)
    z_Press           = 1.0d-60
    z_Press(1)  = &
      & PressSfc *((TempSfc / z_Temp(1)) ** (- CpDryMol /  GasRUniv))
    
    !-----------------------------------------------------------
    ! ��Ǯ��Ψ dT/dz �η׻�. 
    !-----------------------------------------------------------    
    DtDz: do k = 1, kmax-1
      
      !�����
      za_MolFr(k+1,:) = za_MolFr(k,:)
      
      !----------------------------------------------------
      !������Ǯ��Ψ���󥲥��å�ˡ���Ѥ��Ʒ׻�
      !----------------------------------------------------
      ! (0) ���� k �Ǥ��ͤ���������ݴ�
      Temp1  = z_Temp(k)
      Press1 = z_Press(k)
      a_MolFr  = za_MolFr(k,:)

      ! (1) ���� k �Ǥ��ͤ��Ѥ��Ʋ����Ѳ���׻�
      call ECCM_DTempDZ( Temp1, Press1, dz, a_MolFr, DTempDZ1 )

      ! (2) (1) �ǵ�᤿�ͤ��Ѥ���, ���� k + ��k/2 �Ǥ��ͤ��Ѥ��Ʋ����Ѳ���׻�
      !     ���ΤȤ�, ʬ���̤��Ѳ����ʤ���ΤȤ���.
      Temp2  = Temp1 + DTempDZ1 * dz * 5.0d-1
      Press2 =                                     &
        & Press1 * ((Temp1 / Temp2) ** (Grav * MolWtDry / (GasRUniv * DTempDZ1))) 
      call ECCM_DTempDZ( Temp2, Press2, dz, a_MolFr, DTempDZ2 )

      ! (3) (2) �ǵ�᤿�ͤ��Ѥ���, ���� k + ��k/2 �Ǥ��ͤ��Ѥ��Ʋ����Ѳ���׻�
      !     ���ΤȤ�, ʬ���̤��Ѳ����ʤ���ΤȤ���.
      Temp3  = Temp1 + DTempDZ2 * dz * 5.0d-1
      Press3 =                                                                   &
        & Press1 * ((Temp1 / Temp3) ** (Grav * MolWtDry / (GasRUniv * DTempDZ2)))
      call ECCM_DTempDZ( Temp3, Press3, dz, a_MolFr, DTempDZ3 )
      
      ! (4) (3) �ǵ�᤿�ͤ��Ѥ���, ���� k + ��k �Ǥ��ͤ��Ѥ��Ʋ����Ѳ���׻�
      !     ���ΤȤ�, ʬ���̤��Ѳ����ʤ���ΤȤ���.
      Temp4  = Temp1 + DTempDZ3 * dz
      Press4 =                                               &
        & Press1 * ((Temp1 / Temp4) ** (Grav * MolWtDry / (GasRUniv * DTempDZ3)))
      call ECCM_DTempDZ( Temp4, Press4, dz, a_MolFr, DTempDZ4 )
      
      ! (5) �ǽ�Ū�ʷ��������
      DTempDZ = (DTempDZ1 + DTempDZ2 * 2.0d0 + DTempDZ3 * 2.0d0 + DTempDZ4) / 6.0d0

      !----------------------------------------------------
      !����줿���ٸ�Ψ��겹�٤Ȱ��Ϥ����
      !----------------------------------------------------
      !���٤�׻�
      z_Temp(k+1) = z_Temp(k) + DTempDz * dz

      !��ǰ
      if(z_Temp(k+1) < 0.0d0) z_Temp(k+1) = z_Temp(k) 
      
      !���Ϥ��ſ尵ʿ�դ���׻�
      z_Press(k+1) =                                                  &
        &  z_Press(k) * ( ( z_Temp(k) / z_Temp(k+1))                  &
        &    ** (Grav * MolWtDry / ( DTempDZ * GasRUniv ) ) )
      
      !----------------------------------------------------
      !�����η׻�
      !----------------------------------------------------
      do s = 1, CondNum      
        !˰�¾�����
        SatPress = SvapPress( SpcWetID(IdxCC(s)), z_Temp(k+1) )
        
        !�����Υ��ʬΨ���Ѥ��Ƹ��ߤξ�������׻�
        VapPress = za_MolFr(k,IdxCG(s)) * z_Press(k+1)
        
        !˰�¾������Ȱ��Ϥ��鸽�ߤΥ�����׻�
        if ( VapPress > SatPress ) then         
          za_MolFr(k+1,IdxCG(s)) = max(SatPress * Humidity / z_Press(k+1), 1.0d-16)
        end if
      end do
      
      !NH4SH ��ʿ�վ��
      if ( IdxNH3 /= 0 ) then 
        DelMolFr =                                              &
          & max (                                               &
          &    DelMolFrNH4SH(                                   &
          &         z_Temp(k+1), z_Press(k+1),                  &
          &         za_MolFr(k+1,IdxNH3), za_MolFr(k+1,IdxH2S),       &
          &         Humidity                                    &
          &      ),                                             &
          &    0.0d0                                            &
          &  )
        za_MolFr(k+1,IdxNH3) = za_MolFr(k+1,IdxNH3) - DelMolFr
        za_MolFr(k+1,IdxH2S) = za_MolFr(k+1,IdxH2S) - DelMolFr
      end if
      
    end do DtDz
    
  end subroutine ECCM_Wet

!!!------------------------------------------------------------------------------!!!
  subroutine ECCM_MolFr( a_MolFrIni, Humidity, z_Temp, z_Press, za_MolFr )
    !
    ! Ϳ����줿���٤��Ф�, ��������ǮŪ�˾徺�������˼¸������
    ! �����Υץ�ե���������
    !

    !�⥸�塼���ɤ߹���
    use dc_types,   only : DP
    use gridset,    only: kmin, kmax,    &!���󥵥��� (Z ����)
      &                   ncmax           !ʪ����
    use chemcalc,   only: SvapPress,     &!
      &                   DelMolFrNH4SH 
    use composition, only: CondNum,      &!�ŷ�����ο�
      &                    SpcWetID,     &!
      &                    IdxCG,        &!�ŷ����(����)������ź����
      &                    IdxCC,        &!�ŷ����(��)������ź����
      &                    IdxNH3,       &!NH3(����)������ź����
      &                    IdxH2S         !H2S(����)������ź����
        
    !���ۤη�����ػ�
    implicit none
    
    real(DP), intent(in) :: a_MolFrIni(1:ncmax)
    real(DP), intent(in) :: Humidity
    real(DP), intent(in) :: z_Temp(kmin:kmax)
    real(DP), intent(in) :: z_Press(kmin:kmax)
    real(DP), intent(out):: za_MolFr(kmin:kmax, 1:ncmax)
    
    real(DP)             :: DelMolFr
    integer              :: k, s
    
    !-----------------------------------------------------------
    ! ����ν����
    !-----------------------------------------------------------
    do s = 1, ncmax
      za_MolFr(:,s) = a_MolFrIni(s) 
    end do

    !-----------------------------------------------------------
    ! ��Ǯ��Ψ dT/dz �η׻�. 
    !-----------------------------------------------------------
    do k = 1, kmax

      za_MolFr(k,:) = za_MolFr(k-1,:)
      
      !------------------------------------------------------------
      !NH4SH �ʳ��β��ؼ��ʿ�վ��
      !------------------------------------------------------------
      do s = 1, CondNum

        !���������
        !���������Υ��ƥåפǤΥ�����Ķ���뤳�ȤϤʤ�
        za_MolFr(k,IdxCG(s)) =                                 &
          & min(                                                &
          &       za_MolFr(k-1,IdxCG(s)),                      &
          &       SvapPress( SpcWetID(IdxCC(s)), z_Temp(k) ) &
          &        * Humidity / z_Press(k)                      &
          &      )
        
      end do

      !------------------------------------------------------------
      !NH4SH ��ʿ�վ��
      !------------------------------------------------------------
      if ( IdxNH3 /= 0 ) then 
        
        !�������Ѳ�. 
        !�Ȥꤢ���� NH4SH ���Ф���˰����� 1.0 �Ȥ���(��ȴ��...).
        DelMolFr =                                            &
          & max (                                             &
          &    DelMolFrNH4SH(                                 &
          &      z_Temp(k), z_Press(k),                       &
          &      za_MolFr(k,IdxNH3), za_MolFr(k,IdxH2S), Humidity   &
          &     ),                                            &
          &    0.0d0                                          &
          &  )
        
        za_MolFr(k,IdxNH3) = za_MolFr(k,IdxNH3) - DelMolFr 
        za_MolFr(k,IdxH2S) = za_MolFr(k,IdxH2S) - DelMolFr
      end if
      
    end do
  end subroutine ECCM_MolFr

!!!----------------------------------------------------------------------!!!
  subroutine ECCM_DTempDZ( Temp, Press, DelZ, MolFr, DTempDZ )

    !�⥸�塼���ɤ߹���
    use dc_types,    only: DP
    use gridset,     only: ncmax           !ʪ����
    use constants,   only: MolWtDry,      &!
      &                    CpDryMol,      &!
      &                    Grav            !
    use chemcalc,    only: SvapPress,        &!
      &                    LatentHeatPerMol, &!
      &                    ReactHeatNH4SHPerMol, &
      &                    DelMolFrNH4SH 
    use composition, only: CondNum,          &!�ŷ�����ο�
      &                    SpcWetID,         &!
      &                    IdxCG,            &!�ŷ����(����)������ź����
      &                    IdxCC,            &!�ŷ����(��)������ź����
      &                    IdxNH3,           &!NH3(����)������ź����
      &                    IdxH2S             !H2S(����)������ź����
    use ChemData,    only: GasRUniv        
      
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in) :: Temp
    real(DP), intent(in) :: Press
    real(DP), intent(in) :: DelZ
    real(DP), intent(inout) :: MolFr(1:ncmax)    !���ʬΨ
    real(DP), intent(out):: DTempDZ
    real(DP)            :: ReactHeat
    real(DP)            :: Heat(ncmax)
    real(DP)            :: DelMolFr
    real(DP)            :: SatPress
    real(DP)            :: VapPress
    real(DP)            :: Humidity
    real(DP)            :: A, B
    integer             :: s

    !�����
    DTempDZ      = 0.0d0
    ReactHeat    = 0.0d0
    Heat         = 0.0d0
    DelMolFr     = 0.0d0
    SatPress     = 0.0d0
    VapPress     = 0.0d0

    !------------------------------------------------------------
    !NH4SH �ʳ��β��ؼ��ʿ�վ��
    !------------------------------------------------------------
    do s = 1, CondNum      
      
      !˰�¾�����
      SatPress = SvapPress( SpcWetID(IdxCC(s)), Temp )
      
      !��Ǯ. 
      Heat(IdxCG(s)) = LatentHeatPerMol( SpcWetID(IdxCC(s)), Temp )
      
      !�����Υ��ʬΨ���Ѥ��Ƹ��ߤξ�������׻�
      VapPress = MolFr(IdxCG(s)) * Press
      
      !˰�¾���������ŷ��̵ͭ�����
      if ( VapPress < SatPress ) then         
        !�ŷ뤷�Ƥ��ʤ��Τ���Ǯ�ʤ�.
        Heat(IdxCG(s)) = 0.0d0          

      else      

        !˰�¾������Ȱ��Ϥ��鸽�ߤΥ�����׻�
        MolFr(IdxCG(s)) = max(SatPress / Press, 1.0d-16)

      end if
    end do
    
    !------------------------------------------------------------
    !NH4SH ��ʿ�վ��
    !------------------------------------------------------------
    if ( IdxNH3 /= 0 ) then 
      
      Humidity = 1.0d0
      DelMolFr =                                            &
        & max (                                             &
        &    DelMolFrNH4SH(                                 &
        &         Temp, Press, MolFr(IdxNH3), MolFr(IdxH2S),&
        &         Humidity                                  &
        &      ),                                           &
        &    0.0d0                                          &
        &  )
      MolFr(IdxNH3) = MolFr(IdxNH3) - DelMolFr
      MolFr(IdxH2S) = MolFr(IdxH2S) - DelMolFr

      ReactHeat = ReactHeatNH4SHPerMol * DelMolFr
    end if
    
    !------------------------------------------------------------
    !���ٸ��ۤ�׻�
    !------------------------------------------------------------
    !����. ���� Temp(i) ��ɾ��
    A = dot_product( Heat(1:ncmax), MolFr(1:ncmax)) &
      &  / ( GasRUniv * Temp )
    B = dot_product(( Heat(1:ncmax) ** 2.0d0), MolFr(1:ncmax)) &
      &  / ( CpDryMol * GasRUniv * ( Temp ** 2.0d0 ) )
    
    !��Ǯ���ٸ�Ψ
    DTempDZ = - Grav * MolWtDry * (1.0d0 + A) / (CpDryMol * (1.0d0 + B))  &
      &       + ReactHeat / (CpDryMol * DelZ)
    
  end subroutine ECCM_DTempDZ

!!!----------------------------------------------------------------------!!!
  subroutine ECCM_Stab( xyz_PTemp, xyz_Exner, xyzf_QMix )

    !�⥸�塼���ɤ߹���
    use dc_types,    only: DP
    use gridset,     only: imin, imax,    &!���󥵥��� (X ����)
      &                    jmin, jmax,    &!���󥵥��� (Y ����)
      &                    kmin, kmax,    &!���󥵥��� (Z ����)
      &                    ncmax           !ʪ����
    use basicset,    only: xyz_ExnerBZ,   &!
      &                    xyz_PTempBZ,   &! 
      &                    xyzf_QMixBZ
    use constants,   only: MolWtDry,      &!
      &                    CpDry,         &!
      &                    Grav          
    use composition, only: GasNum,        &!���Το�
      &                    MolWtWet
    use average,     only: xyz_xyr
    use differentiate_center2, &
      &              only: xyr_dz_xyz
    
    implicit none

    real(DP), intent(in)  :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyzf_QMix(imin:imax,jmin:jmax,kmin:kmax, ncmax)

    real(DP) :: xyz_Stab(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_StabTemp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_StabMolWt(imin:imax,jmin:jmax,kmin:kmax)

    real(DP)   :: xyzf_MolFrAll(imin:imax,jmin:jmax,kmin:kmax,ncmax)
    real(DP)   :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)   :: xyz_MolWtWet(imin:imax,jmin:jmax,kmin:kmax)
    integer    :: i, j, k, s

    xyz_TempAll = (xyz_PTemp + xyz_PTempBZ) * (xyz_Exner + xyz_ExnerBZ)
    do s = 1, ncmax
      xyzf_MolFrAll(:,:,:,s) =                          &
        &   (xyzf_QMix(:,:,:,s) + xyzf_QMixBZ(:,:,:,s)) &
        &   * MolWtDry / MolWtWet(s) 
    end do
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_MolWtWet(i,j,k) = &
            &     dot_product( MolWtWet(1:GasNum), xyzf_MolFrAll(i,j,k,1:GasNum) )
        end do
      end do
    end do
    
    xyz_StabTemp =                                           &
      &         Grav / xyz_TempAll                           &
      &           * (   xyz_xyr( xyr_dz_xyz( xyz_TempAll ) ) &
      &               + Grav / CpDry ) 
    xyz_StabMolWt =                                          &
      &       - Grav * xyz_xyr( xyr_dz_xyz( xyz_MolWtWet ) ) &
      &         / MolWtDry 
    xyz_Stab = xyz_StabTemp + xyz_StabMolWt

    where (xyz_Stab < 1.0d-7) 
      xyz_Stab = 1.0d-7
    end where

  end subroutine ECCM_Stab

end module ECCM
