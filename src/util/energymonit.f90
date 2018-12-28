!= Module EnergyMonit
!
! Authors::   ODAKA Masatsugu 
! Version::   $Id: energymonit.f90,v 1.4 2014/07/11 08:02:38 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module EnergyMonit
  !
  ! ��¸�̤η׻��Ƚ��Ϥ�Ԥ�����Υ⥸�塼��
  !

  !�⥸�塼���ɤ߹���
  use dc_types,   only: DP, STRING
  
  ! �ѿ����
  !
  real(DP), private                     :: MassTotalBZ   !������(���ܾ�)
  character(STRING), private, parameter :: module_name = 'energymonit'
                                                         ! �⥸�塼���̾��.
                                                         ! Module name
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine EnergyMonit_init
    !
    !�����ѿ��������Ԥ�. 
    !
    ! * DensDev2 �Ȥ����ѿ�̾�ϥ�������, ���Τޤޤ�. total �� DensDev �ǤϤʤ�
    !   Mass �Ȥ����Ȥ���̾�����ѹ����������ɤ��ΤǤϤʤ���
    !

    ! �⥸�塼���ɤ߹���
    !
    use dc_types,   only : DP, STRING
    use gtool_historyauto, &
         &          only : HistoryAutoAddVariable
    use constants,  only : Grav                   !���ϲ�®��
    use basicset,   only : xyz_PressBZ            !���ܾ�ε���
    use gridset,    only : nz                     !�ʻ�����
    use axesset,    only : xmin, xmax, ymin, ymax !�׻��ΰ�
    
    !���ۤη�����ػ�
    !
    implicit none

    ! �����ѿ�������
    !   grav = 0 �λ��˥����ˤʤäƤ��ޤ��Τ�, �ٹ����Ƥ���. 
    !
    MassTotalBZ  = &
      &  ( xyz_PressBz(1,1,1) - xyz_PressBz(1,1,nz) ) &
      &  / max( Grav, 1.0d-40 )                       &
      &  * ( xmax - xmin )                            &
      &  * ( ymax - ymin )

    ! ��������
    !
    call HistoryAutoAddVariable(             &
      & varname='DensDev',                   &
      & dims=(/'t'/),                        &
      & longname='Total density deviation',  &
      & units='kg',                          &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='KinEnrgy',                  &
      & dims=(/'t'/),                        &
      & longname='Total kinetic energy',     &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='ElstEnrgy',                 &
      & dims=(/'t'/),                        &
      & longname='Total elastic energy',     &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='PotEnrgy',                  &
      & dims=(/'t'/),                        &
      & longname='Total potential energy',   &
      & units='J',                           &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='DensDev2',                  &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Density deviation',        &
      & units='kg m-3',                      &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='KinEnrgyDens',              &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Kinetic energy density',   &
      & units='J m-3',                       &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='ElstEnrgyDens',             &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Elastic energy density',   &
      & units='J m-3',                       &
      & xtype='float')

    call HistoryAutoAddVariable(             &
      & varname='PotEnrgyDens',              &
      & dims=(/'x', 'y', 'z', 't'/),         &
      & longname='Potential energy density', &
      & units='J m-3',                       &
      & xtype='float')

  end subroutine EnergyMonit_init

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine EnergyMonit_exec( &
    & pyz_VelX, xqz_VelY, xyr_VelZ, xyz_Exner, xyz_PTemp &
    )
    !
    !��¸�̤η׻��Ƚ��Ϥ�Ԥ�. �׻������̤�
    !����, ��ư���ͥ륮��, �������ͥ륮��, �ݥƥ󥷥�륨�ͥ륮��
    !
    ! �ɲä��٤���(�׸�Ƥ)
    ! * ����ư��
    ! * ���󥿥�ԡ�
    ! * �ȥ졼����������
    ! * ��Ǯ

    use dc_types,   only : DP
    use gtool_historyauto, &
         &          only : HistoryAutoPut
    use gridset,    only : imin,         & !x ����������β���
         &                 imax,         & !x ����������ξ��
         &                 jmin,         & !y ����������β���
         &                 jmax,         & !y ����������ξ��
         &                 kmin,         & !z ����������β���
         &                 kmax,         & !z ����������ξ��
         &                 nx, ny, nz
    use axesset,    only : dx, dy, dz,   & !�ʻ����ֳ�
         &                 xyz_Z           !��ɸ
    use average,    only : xyz_pyz,      &
      &                    xyz_xqz,      &
      &                    xyz_xyr
    use constants,  only : CpDry,        & !�갵��Ǯ
         &                 CvDry,        & !������Ǯ
         &                 GasRDry,      & !�������
         &                 PressSfc,     & !��ɽ����
         &                 Grav            !���ϲ�®��
    use basicset,   only : xyz_ExnerBZ,  & !���ܾ�ΰ��ϴؿ�
         &                 xyz_PTempBZ,  & !���ܾ�β���
         &                 xyz_DensBZ,   & !���ܾ��̩��
         &                 xyz_VelSoundBZ  !��®
    use timeset,    only : TimeN

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP),intent(in) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)   !®��(x��ʬ)
    real(DP),intent(in) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)   !®��(y��ʬ)
    real(DP),intent(in) :: xyr_VelZ(imin:imax,jmin:jmax,kmin:kmax)   !®��(z��ʬ)
    real(DP),intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)  !���ϴؿ�
    real(DP),intent(in) :: xyz_PTemp(imin:imax,jmin:jmax,kmin:kmax)  !����
    
    real(DP) :: xyz_DensDev(imin:imax,jmin:jmax,kmin:kmax) 
                                     !̩��(�к�)
    real(DP) :: xyz_KinEnrgyDens(imin:imax,jmin:jmax,kmin:kmax) 
                                     !��ư���ͥ륮��̩��
    real(DP) :: xyz_ElstEnrgyDens(imin:imax,jmin:jmax,kmin:kmax)
                                     !�������ͥ륮��̩��
    real(DP) :: xyz_PotEnrgyDens(imin:imax,jmin:jmax,kmin:kmax)
                                     !�ݥƥ󥷥�륨�ͥ륮��̩��

    real(DP) :: MassTotalDev         !������(�к�)
    real(DP) :: KinEnrgyTotal        !����ư���ͥ륮��
    real(DP) :: ElstEnrgyTotal       !���������ͥ륮��
    real(DP) :: PotEnrgyTotal        !���ݥƥ󥷥�륨�ͥ륮��


    ! �Ƴʻ������絤����̩���к�
    ! 
    xyz_DensDev =                                             &
      &  PressSfc                                             &
      &  * (                                                  &  
      &      (xyz_ExnerBZ + xyz_Exner) ** ( CvDry / GasRDry ) &
      &       - (xyz_ExnerBZ ) ** ( CvDry / GasRDry )         &
      &    ) / xyz_PTempBZ 
    
    MassTotalDev = sum( dx * dy * dz * xyz_DensDev(1:nx,1:ny,1:nz) )

   
    ! �Ƴʻ����α�ư���ͥ륮��̩��
    ! * ̩�٤ϴ��ܾ���ͤ�ɾ��
    ! 
    xyz_KinEnrgyDens =                     &
      &   0.5d0 * xyz_DensBZ               &
      &   * (                              &
      &        xyz_pyz(pyz_VelX) ** 2.0d0  &
      &      + xyz_xqz(xqz_VelY) ** 2.0d0  &
      &      + xyz_xyr(xyr_VelZ) ** 2.0d0  &
      &     )

    !�Ƴʻ������������ͥ륮��̩��
    ! * ̩�٤ϴ��ܾ���ͤ�ɾ��
    !  
    xyz_ElstEnrgyDens = 0.5d0 * xyz_DensBZ     &
      &  * ( CpDry * xyz_PTempBZ * xyz_Exner / xyz_VelSoundBZ ) ** 2.0d0 
      
    ! �ΰ����Τα�ư/�������ͥ륮��
    !
    KinEnrgyTotal = sum( dx * dy * dz * xyz_KinEnrgyDens(1:nx,1:ny,1:nz)  )
    ElstEnrgyTotal= sum( dx * dy * dz * xyz_ElstEnrgyDens(1:nx,1:ny,1:nz) )


    !�Ƴʻ����Υݥƥ󥷥�륨�ͥ륮��̩��
    ! 
    xyz_PotEnrgyDens = &
      &  - Grav * xyz_DensBZ * xyz_PTemp * xyz_Z / xyz_PTempBZ

    ! �ΰ����ΤΥݥƥ󥷥�륨�ͥ륮��
    !
    PotEnrgyTotal = sum( dx * dy * dz * xyz_PotEnrgyDens(1:nx,1:ny,1:nz))


    !�ե�����ؤν���
    !
    call HistoryAutoPut(TimeN, 'DensDev',   MassTotalDev   )
    call HistoryAutoPut(TimeN, 'KinEnrgy',  KinEnrgyTotal  )
    call HistoryAutoPut(TimeN, 'ElstEnrgy', ElstEnrgyTotal )
    call HistoryAutoPut(TimeN, 'PotEnrgy',  PotEnrgyTotal  )

    call HistoryAutoPut(TimeN, 'DensDev2',      xyz_DensDev(1:nx,1:ny,1:nz)      )
    call HistoryAutoPut(TimeN, 'KinEnrgyDens',  xyz_KinEnrgyDens(1:nx,1:ny,1:nz) )
    call HistoryAutoPut(TimeN, 'ElstEnrgyDens', xyz_ElstEnrgyDens(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'PotEnrgyDens',  xyz_PotEnrgyDens(1:nx,1:ny,1:nz) )

  end subroutine EnergyMonit_exec

end module EnergyMonit
