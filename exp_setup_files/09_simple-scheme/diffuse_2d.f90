! Sample program 
! 
!
program diffuse_2d

  use gridset, only: DimXMin, DimXMax, DimZMin, DimZMax, &
    &                FileXMin, FileXMax, FileZMin, FileZMax, &
    &                x_X, z_Z, gridset_init2                              ! 格子点の設定
!  use average                                                            ! 平均
  use differentiate_center4, only: xz_dx_pz, pz_dx_xz, xz_dz_xr, xr_dz_xz ! 微分
  use boundary, only: BoundaryXCyc_xz, BoundaryZCyc_xz                    ! 境界条件
  use gt4_history

  implicit none

 !---- 空間解像度設定 ----
  integer, parameter :: im=40, km=40            ! 格子点の設定(X,Y)
  integer, parameter :: mgn = 5 
  
  !---- 変数 ----
  real(8) :: xz_ZetaN(-mgn:im+mgn,-mgn:km+mgn)
  real(8) :: xz_ZetaA(-mgn:im+mgn,-mgn:km+mgn)
  
  !---- 座標変数など ----
  real(8), parameter :: xmin=0.0, xmax=1.0
  real(8), parameter :: zmin=0.0, zmax=1.0
  integer, parameter :: snum = 1                ! ダミー

 !---- 時間積分パラメター ----
  real(8), parameter :: dt=1e-4                 ! 時間ステップ間隔
  integer, parameter :: nt=1000, ndisp=20       ! 時間積分数, 表示ステップ

 !---- 物理パラメター ----
  real(8), parameter :: nu=1.0                  ! 粘性係数
  real(8), parameter :: sigma=0.1               ! 初期分布の大きさ
  real(8), parameter :: x1=5.0d-1      ! 初期分布 X 座標
  real(8), parameter :: z1=5.0d-1      ! 初期分布 Y 座標

  integer :: it                                 ! DO 変数
  integer :: i, k


 !---------------- 座標値の設定 ---------------------
  call gridset_init2(im,km,xmin,xmax,zmin,zmax,mgn,snum)    ! 格子点の初期化

 !------------------- 初期値設定 ----------------------
  do k = DimXMin, DimXMax
    do i = DimZMin, DimZMax
      xz_ZetaN(i,k) = dexp(-((x_X(i)-x1)**2 + (z_Z(k)-z1)**2)/ (2*sigma**2))
    end do
  end do

  it = 0
  call output_gtool4_init                            ! ヒストリー初期化
  call output_gtool4

 !------------------- 時間積分 ----------------------
  do it=1,nt    
    write(*,*) '*it = ',it
    xz_ZetaA = xz_ZetaN                 &
      &        + dt * nu * (            &
      &              xz_dx_pz(pz_dx_xz(xz_ZetaN)) &
      &            + xz_dz_xr(xr_dz_xz(xz_ZetaN)) &
      &          )                                   ! Euler 法による時間積分
    call BoundaryXcyc_xz( xz_ZetaA ) 
    call BoundaryZCyc_xz( xz_ZetaA ) 
    xz_ZetaN = xz_ZetaA

    if(mod(it,ndisp) .eq. 0)then                    ! 出力
      call output_gtool4
    endif
  enddo

  call output_gtool4_close()                        ! ヒストリー後処理
  stop

contains

  subroutine output_gtool4_init
    call HistoryCreate( &                                  ! ヒストリー作成
      file='diffuse_2d.nc', title='2D diffusion model',   &
      source='Sample program of deepconv/arare4', &
      institution='GFD_Dennou Club davis project',     &
      dims=(/'x','z','t'/), dimsizes=(/im,km,0/),      &
      longnames=(/'X-coordinate','Z-coordinate','time        '/),&
      units=(/'1','1','1'/),                           &
      origin=0.0, interval=real(ndisp*dt) )
    
    call HistoryPut('x',x_X(FileXMin:FileXMax))               ! 変数出力
    call HistoryPut('z',z_Z(FileZMin:FileZMax))               ! 変数出力

    call HistoryAddVariable( &                                ! 変数定義
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
