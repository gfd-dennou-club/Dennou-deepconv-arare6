!---------------------------------------------------------------------
!     Copyright (C) GFD Dennou Club, 2004, 2005. All rights reserved.
!---------------------------------------------------------------------
!= Module DisturbEnv
!
!   * Developer: SUGIYAMA Ko-ichiro, ODAKA Masatsugu
!   * Version: $Id: initialdata_disturb.f90,v 1.17 2014/07/08 00:59:08 sugiyama Exp $ 
!   * Tag Name: $Name:  $
!   * Change History: 


module initialdata_disturb
  !
  !����Υǥե�����ͤ�Ϳ���뤿��Υ롼����. 
  !
  
  !���ۤη�����ػ�
  implicit none

  !��������
  public initialdata_disturb_random
  public initialdata_disturb_gaussXZ
  public initialdata_disturb_gaussXY
  public initialdata_disturb_gaussYZ
  public initialdata_disturb_gaussXYZ
  public initialdata_disturb_cosXZ
  public initialdata_disturb_cosXY
  public initialdata_disturb_cosYZ
  public initialdata_disturb_cosXYZ
  public initialdata_disturb_coneXZ
  public initialdata_disturb_coneXY
  public initialdata_disturb_coneYZ
  public initialdata_disturb_dryreg
  public initialdata_disturb_square
  public initialdata_disturb_moist
  public initialdata_disturb_tanh
  public initialdata_disturb_tanh_sin
  public initialdata_disturb_circleXZ

contains
    
  subroutine initialdata_disturb_random( DelMax, Zpos, xyz_Var )
    !
    ! �������������Ϳ����
    !

    !�⥸�塼���ɤ߹���
    use mpi_wrapper,only: myrank, nprocs
    use dc_types,   only: DP
    use axesset,    only: z_Z                   ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,    only: nx, imin, imax,      &! ���󥵥��� (X ����)
      &                   ny, jmin, jmax,      &! ���󥵥��� (Y ����)
      &                   kmin, kmax            ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Zpos
    real(DP), intent(out) :: xyz_Var(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: Random           !�ե����뤫������������
    real(DP)              :: Random1(imin:imax, jmin:jmax)
    integer               :: i, j, k, kpos, ix, jy

    ! �����
    xyz_Var = 0.0d0
    kpos    = 1

    ! 0.0--1.0 �ε������ȯ��
    !  mpi �ξ���, �� CPU �λ���������ۤʤ�褦Ĵ�����Ƥ���.  
    !
    do j = jmin, jmax + ( ny * nprocs )
      do i = imin, imax * ( nx * nprocs )
        call random_number(random)
        if (imin + nx * myrank <= i .AND. i <= imax + nx * myrank) then 
          if (jmin + ny * myrank <= j .AND. j <= jmax + ny * myrank) then 
            ix = i - nx * myrank
            jy = j - ny * myrank
            Random1(ix,jy) = random
          end if
        end if
      end do
    end do

    ! ���ꤵ�줿���٤�����ź�����Ѱ�   
    do k = kmin, kmax
      if ( z_Z(k) >= Zpos ) then 
        kpos = k
        exit
      end if
    end do

    ! �������ΤȤ��Ƥϥ���Ȥʤ�褦��Ĵ��. ʿ�Ѥ���κ��ˤ���. 
    do j = 1, ny
      do i = 1, nx
        xyz_Var(i, j, kpos) = &
          & DelMax * (Random1(i,j) - sum( Random1(1:nx,1:ny) ) / real((nx * ny),8))
      end do
    end do
    
  end subroutine initialdata_disturb_random

!!!
!!! Klemp and Wilhelmson (1978) �ν����
!!! 

!  subroutine initialdata_disturb_kw1978(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
!    !
!    ! Klemp and Wilhelmson (1978) �ν����
!
!    implicit none
!
!    real(DP), intent(in)   :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
!    real(DP), intent(out)  :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
!    real(DP)               :: xyz_beta(imin:imax, jmin:jmax, kmin:kmax)
!    integer                :: i, j, k
!
!    do k = kmin, kmax
!      do j = jmin, jmax
!        do i = imin, imax
!          beta(i,j,k) =                               &
!            & (                                       &
!            &      ( ( x_X(i) - Xc ) / Xr ) ** 2.0d0  &
!            &    + ( ( y_Y(j) - Yc ) / Yr ) ** 2.0d0  &
!            &    + ( ( z_Z(k) - Zc ) / Zr ) ** 2.0d0  &
!            &  ) ** 5.0d-1
!        end do
!      end do
!    end do
!    
!    where ( beta < 1.0d0 )
!      xyz_Var = DelMax * ( dcos( Pi * 5.0d-1 * beta ) ** 2.0d0 )
!    end where
!
!  end subroutine initialdata_disturb_kw1978



!!!
!!! tanh
!!! 

  subroutine initialdata_disturb_tanh(VarMean, VarDel, Zc, Zr, aaz_Var, aaz_VarBZ)
    !
    ! tanh ���Υ���
    !   A(z) = A0 + A1 \tanh( (z - Zc) / h )

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)            :: VarMean, VarDel, Zc, Zr
    real(DP), intent(in), optional  :: aaz_VarBZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(out)           :: aaz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer                         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          aaz_Var(i,j,k) = VarMean + VarDel * tanh( (z_Z(k) - Zc) / Zr ) 
        end do
      end do
    end do

    if ( present( aaz_VarBZ) ) then 
      aaz_Var = aaz_Var - aaz_VarBZ
    end if

  end subroutine initialdata_disturb_tanh


  subroutine initialdata_disturb_tanh_sin(VarMean, VarDel, Zc, Zr, aaz_Var, aaz_VarBZ)
    !
    ! tanh ���Υ���
    !   A(z) = A0 + A1 \tanh( (z - Zc) / h )

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,  XMax,    &! X ��ɸ��(�����顼�ʻ���)
      &                  z_Z,  ZMax      ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)            :: VarMean, VarDel, Zc, Zr
    real(DP), intent(in), optional  :: aaz_VarBZ(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), intent(out)           :: aaz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter             :: Pi = 3.1415926535897932385d0 
    integer                         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          aaz_Var(i,j,k) =                                  &
            &   VarMean                                     &
            & + VarDel * tanh( (z_Z(k) - Zc) / Zr )         &
            & + VarDel * sin( x_X(i) * 2.0d0 * Pi / XMax )  &
            &          * sin( z_Z(k) * 2.0d0 * Pi / ZMax )
        end do
      end do
    end do

    if ( present( aaz_VarBZ) ) then 
      aaz_Var = aaz_Var - aaz_VarBZ
    end if

  end subroutine initialdata_disturb_tanh_sin


!!!
!!! �߿��
!!! 

  subroutine initialdata_disturb_ConeXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! �߿��ν����

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  y_Y             ! Y ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((x_X(i) - Xc) / Xr )**2.0d0   &
            &                    +  ((y_Y(j) - Yc) / Yr )**2.0d0) )
          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeXY


  subroutine initialdata_disturb_ConeXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! �߿��ν����
    
    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((x_X(i) - Xc) / Xr )**2.0d0   &
            &                    +  ((z_Z(k) - Zc) / Zr )**2.0d0) )

          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeXZ 


  subroutine initialdata_disturb_ConeYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! �߿��ν����

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: y_Y,           &! Y ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����    
    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * (1.0 - sqrt( ((y_Y(i) - Yc) / Yr )**2.0d0   &
            &                    +  ((z_Z(k) - Zc) / Zr )**2.0d0) )

          xyz_Var(i,j,k) = MAX(xyz_Var(i,j,k), 0.0d0)
        end do
      end do
    end do

  end subroutine initialdata_disturb_ConeYZ

  
!!!
!!! ����
!!!  
  
  subroutine initialdata_disturb_circleXZ(DelMax, Xc, Xr, Zc, xyz_Var)
    !
    ! ���췿�ν����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Zc, Xr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k

    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          if ( ( (x_X(i) - Xc) ** 2.0d0 + (z_Z(k) - Zc) ** 2.0d0 ) < Xr*Xr ) then 
            xyz_Var(i,j,k) = DelMax
          else
            xyz_Var(i,j,k) = 0.0d0
          end if
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_circleXZ


!!!
!!! ����������
!!!  
  
  subroutine initialdata_disturb_gaussXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! ����������ʽ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k

    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1 ) 
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussXZ  

  subroutine initialdata_disturb_gaussXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! ����������ʽ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types, only: DP
    use axesset,  only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                 y_Y             ! Y ��ɸ��(�����顼�ʻ���)
    use gridset,  only: imin, imax,    &! ���󥵥��� (X ����)
      &                 jmin, jmax,    &! ���󥵥��� (Y ����)
      &                 kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1 )
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussXY


  subroutine initialdata_disturb_gaussYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ����������ʽ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: y_Y,           &! Y ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1 )
        end do
      end do
    end do

!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where
    
  end subroutine initialdata_disturb_gaussYZ


  subroutine initialdata_disturb_gaussXYZ(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ����������ʽ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  y_Y,           &! Y ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer               :: i, j, k
    
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          xyz_Var(i,j,k) = &
            & DelMax * dexp( - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1   &
            &                - ( (y_Y(j) - Yc) / Yr )**2.0d0 * 5.0d-1   &
            &                - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1 ) 
        end do
      end do
    end do
    
!    where ( xyz_Var < DelMax * 1.0d-2) 
!      xyz_Var = 0.0d0
!    end where

  end subroutine initialdata_disturb_gaussXYZ

!!!
!!! 
!!!
  subroutine initialdata_disturb_cosXZ(DelMax, Xc, Xr, Zc, Zr, xyz_Var)
    !
    ! ����ή�η׻������Ѥ������
    ! A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: xyz_X,         &! X ��ɸ��(�����顼�ʻ���)
      &                  xyz_Z           ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 

    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXZ
  
  
  subroutine initialdata_disturb_cosXY(DelMax, Xc, Xr, Yc, Yr, xyz_Var)
    !
    ! ����ή�η׻������Ѥ������
    ! A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: xyz_X,         &! X ��ɸ��(�����顼�ʻ���)
      &                  xyz_Y           ! Y ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ����� X �����ξ��
      &                  jmin, jmax,    &! ����� Y �����ξ��
      &                  kmin, kmax      ! ����� Z �����ξ��

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 
    
    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXY
  

  subroutine initialdata_disturb_cosYZ(DelMax, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ����ή�η׻������Ѥ������
    ! A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: xyz_Y,        &! X ��ɸ��(�����顼�ʻ���)
      &                  xyz_Z          ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,   &! ���󥵥��� (X ����)
      &                  jmin, jmax,   &! ���󥵥��� (Y ����)
      &                  kmin, kmax     ! ���󥵥��� (Z ����)
    
    !���ۤη�����ػ�
    implicit none

    real(DP), intent(in)  :: DelMax, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 

    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )

    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosYZ
  
  
  subroutine initialdata_disturb_cosXYZ(DelMax, Xc, Xr, Yc, Yr, Zc, Zr, xyz_Var)
    !
    ! ����ή�η׻������Ѥ������
    ! A [cos(��L) + 1]*0.5 ( L < 1.0 ) or 0.0

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: xyz_X,        &! X ��ɸ��(�����顼�ʻ���)
      &                  xyz_Y,        &! X ��ɸ��(�����顼�ʻ���)
      &                  xyz_Z          ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,   &! ���󥵥��� (X ����)
      &                  jmin, jmax,   &! ���󥵥��� (Y ����)
      &                  kmin, kmax     ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(DP), intent(in)  :: DelMax, Xc, Xr, Yc, Yr, Zc, Zr
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    real(DP)              :: xyz_Var2(imin:imax, jmin:jmax, kmin:kmax)
    real(DP), parameter   :: Pi = 3.1415926535897932385d0 
    
    xyz_Var = 0.0d0   
    xyz_Var2 = sqrt(                       &
         &   + ((xyz_X - Xc) / Xr )**2.0d0 &
         &   + ((xyz_Y - Yc) / Yr )**2.0d0 &
         &   + ((xyz_Z - Zc) / Zr )**2.0d0 &
         &  )
    
    where ( xyz_Var2 <= 1.0d0 )
       xyz_Var = DelMax * ( cos( Pi * xyz_Var2 ) + 1.0d0 ) * 0.5d0
    end where
    
  end subroutine initialdata_disturb_cosXYZ

!!!
!!!---------------------------------------------------------------------
!!!

  subroutine initialdata_disturb_dryreg( &
    & XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, &
    & xyzf_QMix)
    !
    ! ����ʴ����ΰ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,             &! X ��ɸ��(�����顼�ʻ���)
      &                  y_Y,             &! X ��ɸ��(�����顼�ʻ���)
      &                  z_Z               ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,      &! ���󥵥��� (X ����)
      &                  jmin, jmax,      &! ���󥵥��� (Y ����)
      &                  kmin, kmax,      &! ���󥵥��� (Z ����)
      &                  ncmax             ! �׻��ΰ�Υޡ�����
    use basicset, only: xyzf_QMixBZ
    
    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax
    real(DP), intent(out) :: xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, 1:ncmax)
    integer               :: i, j, k, s
    
    ! XposMin:XposMax,ZposMin:ZposMax �ǰϤޤ줿�ΰ�ν���μ��٤򥼥�ˤ��뤿���
    ! ���ܾ�ȵ����ο���������Ϳ����
    do s = 1, ncmax
      do k = kmin,kmax  
        do j = jmin, jmax
          do i = imin,imax
            if (z_Z(k) >= ZposMin .AND. z_Z(k) < ZposMax &
              & .AND. y_Y(j) >= YposMin .AND. y_Y(j) < YposMax &
              & .AND. x_X(i) >= XposMin .AND. x_X(i) < XposMax) then
              xyzf_QMix(i,j,k,s) = - xyzf_QMixBZ(i,j,k,s)
            end if
          end do
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_dryreg
  
  
  subroutine initialdata_disturb_square( &
    & DelMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax, &
    & xyz_Var)
    !
    ! Ω����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use axesset,   only: x_X,           &! X ��ɸ��(�����顼�ʻ���)
      &                  y_Y,           &! Y ��ɸ��(�����顼�ʻ���)
      &                  z_Z             ! Z ��ɸ��(�����顼�ʻ���)
    use gridset,   only: imin, imax,    &! ���󥵥��� (X ����)
      &                  jmin, jmax,    &! ���󥵥��� (Y ����)
      &                  kmin, kmax      ! ���󥵥��� (Z ����)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: DelMax, XposMin, XposMax, YposMin, YposMax, ZposMin, ZposMax
    real(DP), intent(out) :: xyz_Var(imin:imax, jmin:jmax, kmin:kmax)
    integer         :: i, j, k

  
    do k = kmin, kmax  
      do j = jmin, jmax
        do i = imin, imax
          if (z_Z(k) >= ZposMin .AND. z_Z(k) <= ZposMax &
            & .AND. y_Y(j) >= YposMin .AND. y_Y(j) <= YposMax &
            & .AND. x_X(i) >= XposMin .AND. x_X(i) <= XposMax) then
            xyz_Var(i,j,k) = DelMax
!            write(*,*) x_X(i), y_Y(j), z_Z(k), xyz_Var(i,j,k)
          end if
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_square
  
  
  subroutine initialdata_disturb_moist(Hum, xyzf_QMix)
    !
    ! ���� Hum �ʺ�����α�ľʬ�ۤ����
    !

    !�⥸�塼���ɤ߹���
    use dc_types,  only: DP
    use gridset,   only: nx, imin, imax,   &! ���󥵥��� (X ����)
      &                  ny, jmin, jmax,   &! ���󥵥��� (Y ����)
      &                  nz, kmin, kmax,   &! ���󥵥��� (Z ����)
      &                  ncmax              ! �׻��ΰ�Υޡ�����
    use basicset,  only: xyz_TempBZ,       &! ���ܾ�β���
      &                  xyz_PressBZ,      &! ���ܾ�ΰ���
      &                  xyzf_QMixBZ        ! ���ܾ�κ�����
    use composition, only:                 &
      &                  MolWtWet,         &!�Ž���ʬ��ʬ����
      &                  SpcWetMolFr        !�Ž���ʬ�ν�������
    use constants, only: MolWtDry           !������ʬ��ʬ����
    use eccm,      only: eccm_molfr

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(DP), intent(in)  :: Hum
    real(DP), intent(out) :: xyzf_QMix(imin:imax, jmin:jmax, kmin:kmax, 1:ncmax)
    real(DP)              :: zf_MolFr(kmin:kmax, 1:ncmax)
    integer               :: i, j, k, s
  
    ! ���٥���ʤ鲿�⤷�ʤ�
    if ( Hum == 0.0d0 ) return

    ! ��ʿ���ͤʤΤ�, i=0 �����׻�. 
    i = 1
    j = 1
    call eccm_molfr( SpcWetMolFr(1:ncmax), Hum, xyz_TempBZ(i,j,:), &
      &              xyz_PressBZ(i,j,:), zf_MolFr )
    
    !����Υ����򺮹�����Ѵ�
    do s = 1, ncmax
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xyzf_QMix(i,j,k,s) = zf_MolFr(k,s) * MolWtWet(s) / MolWtDry - xyzf_QMixBZ(i,j,k,s)
          end do
        end do
      end do
    end do
    
  end subroutine initialdata_disturb_moist
  
end module initialdata_disturb
