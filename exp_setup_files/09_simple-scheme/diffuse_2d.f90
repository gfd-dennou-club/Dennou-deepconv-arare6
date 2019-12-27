! Sample program 
! 
!
program diffuse_2d

  use gridset, only: DimXMin, DimXMax, DimZMin, DimZMax, &
    &                FileXMin, FileXMax, FileZMin, FileZMax, &
    &                x_X, z_Z, gridset_init2                              ! �ʻ���������
!  use average                                                            ! ʿ��
  use differentiate_center4, only: xz_dx_pz, pz_dx_xz, xz_dz_xr, xr_dz_xz ! ��ʬ
  use boundary, only: BoundaryXCyc_xz, BoundaryZCyc_xz                    ! �������
  use gt4_history

  implicit none

 !---- ���ֲ��������� ----
  integer, parameter :: im=40, km=40            ! �ʻ���������(X,Y)
  integer, parameter :: mgn = 5 
  
  !---- �ѿ� ----
  real(8) :: xz_ZetaN(-mgn:im+mgn,-mgn:km+mgn)
  real(8) :: xz_ZetaA(-mgn:im+mgn,-mgn:km+mgn)
  
  !---- ��ɸ�ѿ��ʤ� ----
  real(8), parameter :: xmin=0.0, xmax=1.0
  real(8), parameter :: zmin=0.0, zmax=1.0
  integer, parameter :: snum = 1                ! ���ߡ�

 !---- ������ʬ�ѥ�᥿�� ----
  real(8), parameter :: dt=1e-4                 ! ���֥��ƥå״ֳ�
  integer, parameter :: nt=1000, ndisp=20       ! ������ʬ��, ɽ�����ƥå�

 !---- ʪ���ѥ�᥿�� ----
  real(8), parameter :: nu=1.0                  ! Ǵ������
  real(8), parameter :: sigma=0.1               ! ���ʬ�ۤ��礭��
  real(8), parameter :: x1=5.0d-1      ! ���ʬ�� X ��ɸ
  real(8), parameter :: z1=5.0d-1      ! ���ʬ�� Y ��ɸ

  integer :: it                                 ! DO �ѿ�
  integer :: i, k


 !---------------- ��ɸ�ͤ����� ---------------------
  call gridset_init2(im,km,xmin,xmax,zmin,zmax,mgn,snum)    ! �ʻ����ν����

 !------------------- ��������� ----------------------
  do k = DimXMin, DimXMax
    do i = DimZMin, DimZMax
      xz_ZetaN(i,k) = dexp(-((x_X(i)-x1)**2 + (z_Z(k)-z1)**2)/ (2*sigma**2))
    end do
  end do

  it = 0
  call output_gtool4_init                            ! �ҥ��ȥ꡼�����
  call output_gtool4

 !------------------- ������ʬ ----------------------
  do it=1,nt    
    write(*,*) '*it = ',it
    xz_ZetaA = xz_ZetaN                 &
      &        + dt * nu * (            &
      &              xz_dx_pz(pz_dx_xz(xz_ZetaN)) &
      &            + xz_dz_xr(xr_dz_xz(xz_ZetaN)) &
      &          )                                   ! Euler ˡ�ˤ�������ʬ
    call BoundaryXcyc_xz( xz_ZetaA ) 
    call BoundaryZCyc_xz( xz_ZetaA ) 
    xz_ZetaN = xz_ZetaA

    if(mod(it,ndisp) .eq. 0)then                    ! ����
      call output_gtool4
    endif
  enddo

  call output_gtool4_close()                        ! �ҥ��ȥ꡼�����
  stop

contains

  subroutine output_gtool4_init
    call HistoryCreate( &                                  ! �ҥ��ȥ꡼����
      file='diffuse_2d.nc', title='2D diffusion model',   &
      source='Sample program of deepconv/arare4', &
      institution='GFD_Dennou Club davis project',     &
      dims=(/'x','z','t'/), dimsizes=(/im,km,0/),      &
      longnames=(/'X-coordinate','Z-coordinate','time        '/),&
      units=(/'1','1','1'/),                           &
      origin=0.0, interval=real(ndisp*dt) )
    
    call HistoryPut('x',x_X(FileXMin:FileXMax))               ! �ѿ�����
    call HistoryPut('z',z_Z(FileZMin:FileZMax))               ! �ѿ�����

    call HistoryAddVariable( &                                ! �ѿ����
      varname='zeta', dims=(/'x','z','t'/), & 
      longname='vorticity', units='1', xtype='double')
  end subroutine output_gtool4_init

  subroutine output_gtool4
    write(*,*) 'it = ',it
    call HistoryPut('zeta', xz_ZetaN(FileXMin:FileXMax, FileZMin:FIleZMax))
  end subroutine output_gtool4

  subroutine output_gtool4_close()
    call HistoryClose
  end subroutine output_gtool4_close
  
end program diffuse_2d
