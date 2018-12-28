!= ���ݥ��إ⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: damping.f90,v 1.14 2014/11/07 06:46:45 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module Damping
  !
  ! ���ݥ��� (�����ն���Ȥ�ȿ�ͤ��ޤ��ۼ����뤿�����) �ˤ�����
  ! ����Ψ�Ȥ��η׻���Ԥ�����Υѥå��������⥸�塼��
  !

  !�⥸�塼���ɤ߹���
  use dc_types, only : DP

  !���ۤη�����ػ�
  implicit none

  !private °�������
  private 
  
  !�ؿ��ˤ� public °�������
  public Damping_Init
  public SpongeLayer_forcing
  public SpongeLayer_MeanFlow
  
  !�ѿ����
  real(DP), save, private              :: EFTime = 5000.0d0!���ݥ��ؤ� e-folding time
  real(DP), save, private              :: DepthH = 0.0d0   !���ݥ��ؤθ���(��ʿ����)
  real(DP), save, private              :: DepthV = 0.0d0   !���ݥ��ؤθ���(��ľ����) [��������]
  real(DP), save, private              :: DepthVb= 0.0d0   !���ݥ��ؤθ���(��ľ����) [��������]
  real(DP), allocatable, save, private :: xyz_Gamma(:,:,:) !xyz �ʻҸ��그��(��ʿ����)
  real(DP), allocatable, save, private :: pyz_Gamma(:,:,:) !pyz �ʻҸ��그��(��ľ����)
  real(DP), allocatable, save, private :: xqz_Gamma(:,:,:) !xqz �ʻҸ��그��(��ľ����)
  real(DP), allocatable, save, private :: xyr_Gamma(:,:,:) !xyr �ʻҸ��그��(��ľ����)
  real(DP), save, private              :: FactorSpngVelX  = 1.0d0
  real(DP), save, private              :: FactorSpngVelY  = 1.0d0
  real(DP), save, private              :: FactorSpngVelZ  = 1.0d0
  real(DP), save, private              :: FactorSpngPTemp = 1.0d0
  real(DP), save, private              :: FactorSpngExner = 0.0d0
  real(DP), save, public               :: DampSound  = 0.0d0   !���ȸ����θ��그��

  real(DP), save, private              :: EFTimeMF = 50000.0d0
  real(DP), save, private              :: DelTimeMF = 1.0d0
  real(DP), save, private              :: ZMinMF  = 0.0d0
  real(DP), save, private              :: ZMaxMF  = 1.0d10

contains 
  
!!!------------------------------------------------------------------------!!!
  subroutine Damping_Init
    !
    ! ���ȸ����ȥ��ݥ��ؤθ��그���ν����
    ! 
    use dc_types,   only: DP
    use dc_iounit,  only: FileOpen
    use dc_message, only: MessageNotify
    use gtool_historyauto, only: HistoryAutoAddVariable
    use gridset, only: imin,       &! x ����������β���
      &                imax,       &! x ����������ξ��
      &                jmin,       &! y ����������β���
      &                jmax,       &! y ����������ξ��
      &                kmin,       &! z ����������β���
      &                kmax         ! z ����������ξ��
    use axesset, only: x_X,        &!X ��ɸ��(�����顼�ʻ���)
      &                y_Y,        &!Y ��ɸ��(�����顼�ʻ���)
      &                z_Z,        &!Z ��ɸ��(�����顼�ʻ���)
      &                p_X,        &!X ��ɸ��(�ե�å����ʻ���)
      &                q_Y,        &!Y ��ɸ��(�ե�å����ʻ���)
      &                r_Z,        &!Z ��ɸ��(�ե�å����ʻ���)
      &                dx, dy, dz, &! �ʻҴֳ�
      &                XMax,       &!X ��ɸ�κ�����
      &                YMax,       &!Y ��ɸ�κ�����
      &                ZMin,       &!Z ��ɸ�κǾ��� 
      &                ZMax         !Z ��ɸ�κ����� 
    use namelist_util,  &
      &          only: namelist_filename
    use timeset, only: DelTimeShort    !û�����֥��ƥå�
    
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), parameter       :: Pi =3.1415926535897932385d0   !�߼�Ψ
    real(DP)                  :: Alpha = 5.0d-2
    integer                   :: unit ! �����ֹ�
    integer                   :: i, j, k

    !NAMELIST �������
    NAMELIST /damping_nml/ Alpha,                        &
      & EFTime, DepthH, DepthV, DepthVb,                 &
      & FactorSpngVelX, FactorSpngVelY, FactorSpngVelZ,  &
      & FactorSpngPTemp, FactorSpngExner,                &
      & DelTimeMF, ZMinMF, ZMaxMF, EFTimeMF

    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=damping_nml)
    close(unit)

    ! ���ȸ����
    DampSound = Alpha * ( Min(dx, dz) ** 2.0d0 ) / DelTimeShort

    !�����
    allocate( &
      & xyz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & pyz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & xqz_Gamma(imin:imax,jmin:jmax,kmin:kmax), &
      & xyr_Gamma(imin:imax,jmin:jmax,kmin:kmax)    )
    xyz_Gamma = 0.0d0
    pyz_Gamma = 0.0d0
    xqz_Gamma = 0.0d0
    xyr_Gamma = 0.0d0

    !-----------------------------------------------------------------    
    ! �����Υ����å�
    !
    if ( DepthH > 0.0d0 ) then 
      if ( DepthH < dx ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelX is %f", d=(/dx/))        
      else if ( DepthH < dx ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelX is %f", d=(/dx/))
      end if
      
      if ( DepthH < dy ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelY is %f", d=(/dy/))
      else if ( DepthH < dy ) then 
        call MessageNotify( "E", "Damping_init", "DepthH is too thin. DelY is %f", d=(/dy/))
      end if
    end if

    if ( DepthV > 0.0d0 ) then 
      if ( DepthV < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthV is too thin. DelZ is %f", d=(/dz/) )      
      elseif ( DepthV < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthV is too thin. DelZ is %f", d=(/dz/) )      
      end if
    end if

    if ( DepthVb > 0.0d0 ) then 
      if ( DepthVb < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthVb is too thin. DelZ is %f", d=(/dz/) )      
      elseif ( DepthVb < dz ) then 
        call MessageNotify( "E", "Damping_init", "DepthVb is too thin. DelZ is %f", d=(/dz/) )      
      end if
    end if

    !-----------------------------------------------------------------    
    ! ���ݥ��ؤθ���Ψ
    !

    !��ʿ��������¦����¦����
    !
    if ( DepthH > 0.0d0 ) then 
      do i = imin, imax
        !�����顼�ʻ�������¦����
        if ( x_X(i) < DepthH) then 
          xyz_Gamma(i,:,:) = xyz_Gamma(i,:,:) &
            & + ((1.0d0 - x_X(i) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !�ե�å����ʻ�������¦����
        if ( p_X(i) < DepthH) then 
          pyz_Gamma(i,:,:) = pyz_Gamma(i,:,:) &
            & + ((1.0d0 - p_X(i) / DepthH) ** 3.0d0) / EFTime
        end if
        
      !�����顼�ʻ�������¦����    
        if ( x_X(i) > ( XMax - DepthH ) ) then 
          xyz_Gamma(i,:,:) = xyz_Gamma(i,:,:) &
            & + ((1.0d0 - (XMax - x_X(i)) / DepthH) ** 3.0d0) / EFTime 
        end if
        
        !�ե�å����ʻ�������¦����    
        if ( p_X(i) > ( XMax - DepthH ) ) then 
          pyz_Gamma(i,:,:) = pyz_Gamma(i,:,:) &
            & + ((1.0d0 - (XMax - p_X(i)) / DepthH) ** 3.0d0) / EFTime 
        end if
      end do

      ! x �����ˤ�Ʊ��
      !
      xyr_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if
    
    !��ʿ��������¦����¦����
    !
    if ( DepthH > 0.0d0 ) then 
      do j = jmin, jmax
        !�����顼�ʻ�������¦����
        if ( y_Y(j) < DepthH) then 
          xyz_Gamma(:,j,:) = xyz_Gamma(:,j,:) &
            & + ((1.0d0 - y_Y(j) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !�ե�å����ʻ�������¦����
        if ( q_Y(j) < DepthH) then 
          xqz_Gamma(:,j,:) = xqz_Gamma(:,j,:) &
            & + ((1.0d0 - q_Y(j) / DepthH) ** 3.0d0) / EFTime
        end if
        
        !�����顼�ʻ�������¦����    
        if ( y_Y(j) > ( YMax - DepthH ) ) then 
          xyz_Gamma(:,j,:) = xyz_Gamma(:,j,:) &
            & + ((1.0d0 - (YMax - y_Y(j)) / DepthH) ** 3.0d0) / EFTime 
        end if
        
        !�ե�å����ʻ�������¦����    
        if ( q_Y(j) > ( YMax - DepthH ) ) then 
          xqz_Gamma(:,j,:) = xqz_Gamma(:,j,:)  &
            & + ((1.0d0 - (YMax - q_Y(j)) / DepthH) ** 3.0d0) / EFTime 
        end if
      end do

      ! y �����ˤ�Ʊ��
      !
      pyz_Gamma  = xyz_Gamma
      xyr_Gamma  = xyz_Gamma
    end if
    
    !��ľ�����ξ�������    
    !
    if ( DepthV > 0.0d0 ) then 
      do k = kmin, kmax
        !�����顼�ʻ���
        if ( z_Z(k) >= ( ZMax - DepthV ) ) then 
          xyz_Gamma(:,:,k) = xyz_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * (z_Z(k) - ZMax + DepthV) / DepthV)) &
            &    / EFTime 
        end if
      
        !�ե�å����ʻ���
        if ( r_Z(k) >= ( ZMax - DepthV ) ) then 
          xyr_Gamma(:,:,k) = xyr_Gamma(:,:,k) &
            & + (1.0d0 - dcos(Pi * (r_Z(k) - ZMax + DepthV)/ DepthV)) &
            &    / EFTime 
        end if
      end do
      
      ! z �����ˤ�Ʊ��
      !
      pyz_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if

    !��ľ�����β�������    
    !
    if ( DepthVb > 0.0d0 ) then 
      do k = kmin, kmax
        !�����顼�ʻ���
        if ( z_Z(k) <= (ZMin + DepthV) ) then 
          xyz_Gamma(:,:,k) = xyz_Gamma(:,:,k) &
            & + (1.0d0 + dcos(Pi * ( z_Z(k) - ZMin ) / DepthV) ) &
            &    / EFTime 
        end if
        
        !�ե�å����ʻ���
        if ( r_Z(k) <= (ZMin + DepthV) ) then 
          xyr_Gamma(:,:,k) = xyr_Gamma(:,:,k) &
            & + (1.0d0 + dcos(Pi * ( r_Z(k) - ZMin )/ DepthV ) ) &
            &    / EFTime 
        end if
      end do
      
      ! z �����ˤ�Ʊ��
      !
      pyz_Gamma  = xyz_Gamma
      xqz_Gamma  = xyz_Gamma
    end if
      
    ! ��ǧ
    ! 
!    i = 1; j = 1
!    write(*,*) 'z_Z,   xyz_Gamma'
!    do k = kmin, kmax      
!      write(*,*) z_Z(k), '=>', xyz_Gamma(i,j,k)
!    end do
!    write(*,*) 'r_Z,   xyr_Gamma'
!    do k = kmin, kmax      
!      write(*,*) r_Z(k), '=>', xyr_Gamma(i,j,k)
!    end do

    !-----------------------------------------------------------------    
    ! �ͤγ�ǧ
    !
    call MessageNotify( "M", "Damping_init", "EFTime = %f", d=(/EFTime/) )
    call MessageNotify( "M", "Damping_init", "DepthH = %f", d=(/DepthH/) )
    call MessageNotify( "M", "Damping_init", "DepthV = %f", d=(/DepthV/) )  
    call MessageNotify( "M", "Damping_init", "DepthVb= %f", d=(/DepthVb/) )  

    !-----------------------------------------------------------------    
    ! ����
    !
    call HistoryAutoAddVariable(                          &
      & varname='DPTempDtSpng',                           &
      & dims=(/'x','y','z','t'/),                         &
      & longname='Damping term of potential temperature', &
      & units='K.s-1',                                    &
      & xtype='double')

    call HistoryAutoAddVariable(                    &
      & varname='ExnerSpng',                        &
      & dims=(/'x','y','z','t'/),                   &
      & longname='Damping term of exner function',  &
      & units='s-1',                                &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelXDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelX',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelYDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelY',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(          &
      & varname='DVelZDtSpng',            &
      & dims=(/'x','y','z','t'/),         &
      & longname='Damping term of VelZ',  &
      & units='m.s-1',                    &
      & xtype='double')

    call HistoryAutoAddVariable(                      &
      & varname='DVelXDtSpngMF',                      &
      & dims=(/'x','y','z','t'/),                     &
      & longname='Damping term of VelX (Mean Flow)',  &
      & units='m.s-1',                                &
      & xtype='double')

    call HistoryAutoAddVariable(                      &
      & varname='DVelYDtSpngMF',                      &
      & dims=(/'x','y','z','t'/),                     &
      & longname='Damping term of VelY (Mean Flow)',  &
      & units='m.s-1',                                &
      & xtype='double')

  end subroutine Damping_Init


  subroutine SpongeLayer_forcing(                                         &
    & pyz_VelXBl,  xqz_VelYBl,  xyr_VelZBl,  xyz_PTempBl,  xyz_ExnerBl,   & !(in) 
    & pyz_DVelXDt, xqz_DVelYDt, xyr_DVelZDt, xyz_DPTempDt, xyz_DExnerDt )   !(inout)

    use dc_types, only : DP
    use gtool_historyauto, only: HistoryAutoPut
    use timeset, only: TimeN        ! ���ߤλ���
    use gridset, only: imin,       &! x ����������β���
      &                imax,       &! x ����������ξ��
      &                jmin,       &! y ����������β���
      &                jmax,       &! y ����������ξ��
      &                kmin,       &! z ����������β���
      &                kmax,       &! z ����������ξ��
      &                nx,         &! x ������ʪ���ΰ�ξ��
      &                ny,         &! y ������ʪ���ΰ�ξ��
      &                nz           ! z ������ʪ���ΰ�ξ��

    !���ۤη�����ػ�
    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyr_VelZBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyz_PTempBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xyz_ExnerBl(imin:imax, jmin:jmax, kmin:kmax)   
    real(DP), intent(inout) :: pyz_DVelXDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyr_DVelZDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: pyz_SpngVelX(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xqz_SpngVelY(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyr_SpngVelZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_SpngPTemp(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xyz_SpngExner(imin:imax, jmin:jmax, kmin:kmax)

    !--------------------------------------------------------
    ! �������ͤ򥼥���᤹�褦�ʶ�����Ϳ���� 
    !
    pyz_SpngVelX  =  - pyz_Gamma * pyz_VelXBl  * FactorSpngVelX
    xqz_SpngVelY  =  - xqz_Gamma * xqz_VelYBl  * FactorSpngVelY
    xyr_SpngVelZ  =  - xyr_Gamma * xyr_VelZBl  * FactorSpngVelZ
    xyz_SpngPTemp =  - xyz_Gamma * xyz_PTempBl * FactorSpngPTemp
    xyz_SpngExner =  - xyz_Gamma * xyz_ExnerBl * FactorSpngExner
    
    pyz_DVelXDt  = pyz_DVelXDt  + pyz_SpngVelX
    xqz_DVelYDt  = xqz_DVelYDt  + xqz_SpngVelY
    xyr_DVelZDt  = xyr_DVelZDt  + xyr_SpngVelZ
    xyz_DPTempDt = xyz_DPTempDt + xyz_SpngPTemp
    xyz_DExnerDt = xyz_DExnerDt + xyz_SpngExner
    
    call HistoryAutoPut(TimeN, 'DVelXDtSpng',  pyz_SpngVelX(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSpng',  xqz_SpngVelY(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelZDtSpng',  xyr_SpngVelZ(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DPTempDtSpng', xyz_SpngPTemp(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerSpng', xyz_SpngExner(1:nx,1:ny,1:nz))

  end subroutine SpongeLayer_forcing


  subroutine SpongeLayer_MeanFlow( &
    & pyz_VelXBl,  xqz_VelYBl,     & !(in) 
    & pyz_DVelXDt, xqz_DVelYDt     ) !(inout)

    use mpi_wrapper, only: nprocs, myrank, MPIWrapperGather, MPIWrapperAllreduce
    use gtool_historyauto, only: HistoryAutoPut
    use dc_types, only: DP
    use timeset, only: TimeN        ! ���ߤλ���
    use axesset, only: z_Z          !Z ��ɸ��(�����顼�ʻ���)
    use gridset, only: imin,       &! x ����������β���
      &                imax,       &! x ����������ξ��
      &                jmin,       &! y ����������β���
      &                jmax,       &! y ����������ξ��
      &                kmin,       &! z ����������β���
      &                kmax,       &! z ����������ξ��
      &                nx,         &! x ������ʪ���ΰ�ξ��
      &                ny,         &! y ������ʪ���ΰ�ξ��
      &                nz           ! z ������ʪ���ΰ�ξ��

    !���ۤη�����ػ�
    implicit none

    real(DP), intent(in)    :: pyz_VelXBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(in)    :: xqz_VelYBl(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: pyz_DVelXDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(inout) :: xqz_DVelYDt(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: pyz_SpngVelX(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: xqz_SpngVelY(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)                :: z_SpngVelX(kmin:kmax)
    real(DP)                :: z_SpngVelY(kmin:kmax)
    real(DP)                :: z_VelX(nz)
    real(DP)                :: z_VelXMean(nz)
    real(DP)                :: z_VelXSum(nz)
    real(DP)                :: a_VelX(nz * nprocs)
    real(DP)                :: za_VelX(nz, nprocs)
    real(DP)                :: z_VelY(nz)
    real(DP)                :: z_VelYMean(nz)
    real(DP)                :: z_VelYSum(nz)
    real(DP)                :: a_VelY(nz * nprocs)
    real(DP)                :: za_VelY(nz, nprocs)
    integer                 :: i, j, k, n, n0, n1


    !------------------------------------------
    ! ʿ��ή�����פ��뤫�Υ����å�. 
    !
    if (mod(TimeN, DelTimeMF) /= 0.0d0 ) return 
    

    !------------------------------------------
    ! �Ρ������ʿ���ͤκ���
    !
    z_VelX = 0.0d0
    z_VelY = 0.0d0

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          z_VelX(k) = z_VelX(k) + pyz_VelXBl(i,j,k)
          z_VelY(k) = z_VelY(k) + xqz_VelYBl(i,j,k)
        end do
      end do
    end do
    
    z_VelX = z_VelX / real( nx * ny, kind=8)
    z_VelY = z_VelY / real( nx * ny, kind=8)

    !--------------------------------------------------------
    ! MPI : ���Ρ��ɤǤι�פ�����
    !
    call MPIWrapperAllreduce(   &
      &  nz, z_VelX,        & !(in)
      &  z_VelXSum          & !(out)
      & )

    call MPIWrapperAllreduce(   &
      &  nz, z_VelY,        & !(in)
      &  z_VelYSum          & !(out)
      & )

    !--------------------------------------------------------
    ! ����ԥ󥰷���   
    !
    do k = kmin, kmax
      if ( z_Z(k) >= ZMinMF .OR.  z_Z(k) <= ZMaxMF ) then 
        pyz_SpngVelX(:,:,k) = - z_VelXSum(k) / real( nprocs, kind=8 ) / EFTimeMF
        xqz_SpngVelY(:,:,k) = - z_VelYSum(k) / real( nprocs, kind=8 ) / EFTimeMF
      else
        pyz_SpngVelX(:,:,k) = 0.0d0
        xqz_SpngVelY(:,:,k) = 0.0d0
      end if
    end do

!!    if ( myrank == 0 ) then 
!!      write(*,*) "z_VelXSum / nprocs : ", z_VelXSum(1:10) / real( nprocs, kind=8 )
!!      write(*,*) "z_VelYSum / nprocs : ", z_VelYSum / real( nprocs, kind=8 )
!!    end if

    !--------------------------------------------------------
    ! �������ͤ򥼥���᤹�褦�ʶ�����Ϳ���� 
    !
    pyz_DVelXDt  = pyz_DVelXDt  + pyz_SpngVelX
    xqz_DVelYDt  = xqz_DVelYDt  + xqz_SpngVelY
    
    call HistoryAutoPut(TimeN, 'DVelXDtSpngMF',  pyz_SpngVelX(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'DVelYDtSpngMF',  xqz_SpngVelY(1:nx,1:ny,1:nz))

  end subroutine SpongeLayer_MeanFlow
  
end module Damping
