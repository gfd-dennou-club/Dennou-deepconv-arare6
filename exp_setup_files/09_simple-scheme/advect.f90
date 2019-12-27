! Sample program 
! 
!
program advect

!!!
!!! モジュールの設定. 
!!!

  ! gtool5 モジュール (I/O) 外部ライブラリ.
  ! gtool5 : http://www.gfd-dennou.org/arch/gtool/
  !
  use dc_types,          only: STRING, DP   ! DP は double precision の意味
  use gtool_history,     only: HistoryCreate, HistoryAddVariable, HistoryPut, HistoryClose

  ! 初期設定モジュール (src/setup/ 以下のファイル)
  !
  use mpi_wrapper,     only : MPIWrapperInit,                          &
    &                         MPIWrapperFinalize
  use argset,          only : argset_init
  use gridset,         only : gridset_init,                            &
    &                         imin, imax, jmin, jmax, kmin, kmax,      &
    &                         nx, ny, nz, ncmax
  use timeset,         only : timeset_init,                            &
    &                         TimesetDelTimeHalf, TimesetProgress,     &
    &                         TimeA, TimeN, TimeB,                     &
    &                         Nstep, NstepShort, NstepLong,            &
    &                         NstepOutput, FlagInitialRun,             &
    &                         DelTimeLong
  use axesset,         only : axesset_init, x_X, y_Y, z_Z, XMax, ZMax
  use namelist_util,   only : NmlutilInit
  
  ! 下請けモジュール
  !
  use setmargin,       only : SetMargin_init
  use setmargin,       only : SetMargin_xyz
  use differentiate_center4
  use average
  
!!!
!!! 暗黙の型宣言禁止
!!! 
  implicit none

!!!
!!! 変数宣言 
!!!  
  ! 予報変数
  !
  real(DP), allocatable :: xyz_ZetaAl(:,:,:)    ! 次の時刻の値
  real(DP), allocatable :: xyz_ZetaNl(:,:,:)    ! 現在の値
  real(DP), allocatable :: xyz_ZetaBl(:,:,:)    ! 前の時刻の値
  real(DP), allocatable :: xyz_DZetaDtNl(:,:,:) ! tendency
  
  ! 物理パラメタ
  ! 設定ファイル (.conf) では指定できない項目を陽に指定.  
  !
!  real(DP), parameter :: nu    = 1.0d-1  ! 拡散係数
  real(DP), parameter :: VelX0 = 1.0d0   ! 移流速度 X 方向
  real(DP), parameter :: VelZ0 = 0.0d0   ! 移流速度 Z 方向
  real(DP), parameter :: sigma = 0.1     ! 初期分布の大きさ
  real(DP), parameter :: DelMax = 1.0d0  ! 山の高さ
  real(DP), parameter :: Xc = 5.0d2      ! 山の中心位置 (割合) (0 < Xc < 1)
  real(DP), parameter :: Zc = 5.0d2      ! 山の中心位置 (割合) (0 < Zc < 1)
  real(DP), parameter :: Xr = 3.0d2      ! 山の幅
  real(DP), parameter :: Zr = 3.0d2      ! 山の幅

  ! I/O
  !
  character(STRING) :: cfgfile ! NAMELIST ファイル名 ; NAMELIST fine name

  ! do ループ変数 ; do loop variable
  !
  integer :: it = 0
  integer :: i, k


!!!
!!! 初期化 (main/arare.f90 の MainInit サブルーチンの簡単化)
!!! 
  call MainInit
  
!!!
!!! 初期値設定.
!!!

  ! 中心に山を与える
  !
  do k = kmin, kmax
     do i = imin, imax
        xyz_ZetaNl(i,:,k) =                                    & 
             &    DelMax                                       &
             &  * dexp(                                        &
             &       - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1  &
             &       - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1  &
             &    )        
    end do
  end do

!  write(*,*) xyz_ZetaNl(95:105,1,100)

  ! 境界条件
  !
  call SetMargin_xyz( xyz_ZetaNl )
  xyz_ZetaBl = xyz_ZetaNl

  
  !初期値の出力
  !
  call HistoryPut('zeta', xyz_ZetaNl(1:nx,1:ny,1:nz))

!!!
!!! 時間積分
!!!
  do while ( Nstep <= NstepLong )

     ! tendency (傾き) の計算.
     ! 計算の安定性のためには数値拡散項が必要だが今は無視.
     ! 長時間計算できないが今はよしとする.
     ! デフォルトでは VelZ0 = 0 としている.
     ! 
     xyz_DZetaDtNl =                                     &
          &  - VelX0 * xyz_pyz(pyz_dx_xyz(xyz_ZetaNl))   &
          &  - VelZ0 * xyz_xyr(xyr_dz_xyz(xyz_ZetaNl))   


     ! 温位の積分 (leap-frog)
     !
     xyz_ZetaAl = xyz_ZetaBl + 2.0d0 * DelTimeLong * xyz_DZetaDtNl

     write(*,*) xyz_ZetaAl(95:105,1,100)
     
     ! Asselin のタイムフィルターを利用
     !
     call AsselinTimeFilter
     
     ! 境界条件 ; Boundary condition
     !
     call SetMargin_xyz( xyz_ZetaAl )  
          
     ! ヒストリファイル出力. 
     !
     if ( it == NstepOutput ) then
        call HistoryPut('zeta', xyz_ZetaNl(1:nx,1:ny,1:nz))
        it = 0
     end if
     
     ! ループを回すための処理
     !
     xyz_ZetaNl = xyz_ZetaAl
     xyz_ZetaBl = xyz_ZetaNl
     it = it + 1
     
     ! 時刻の進行
     !
     call TimesetProgress

  end do

!!!
!!! 終了処理
!!!

  ! 出力ファイルのクローズ
  !
  call HistoryClose

  ! MPI END
  !
  call MPIWrapperFinalize
  
  stop

contains

!!!
!!! 初期化プロセス (src/main/arare.f90 の内部サブルーチンの簡略化)  
!!!
  subroutine MainInit
  
    ! MPI
    !
    call MPIWrapperInit

    ! NAMELIST ファイル名の読み込み
    !
    call argset_init( cfgfile ) !(out)
    
    ! NAMELIST ファイル名のモジュールへの登録
    !
    call NmlutilInit( cfgfile ) !(in)
    
    ! 時間積分の初期化
    !
    call timeset_init
    
    ! 格子点情報の初期化
    !
    call gridset_init
    
    ! 軸の情報の初期化                    
    !! 注) axesset と gridset をマージしても構わない.
    !
    call axesset_init
    
    ! I/O ファイル名の初期化 
    !! 注) 内部サブルーチン化.
    !
    call fileset_init
    
    ! マージンの設定の初期化
    !
    call SetMargin_Init
    
    ! 内部変数の初期化
    !
    call VariableAllocate
    
    ! t=0 から数値積分を行う場合は, 最初の一回目はオイラー法を利用. 
    !
    if ( FlagInitialRun ) call TimesetDelTimeHalf

  end subroutine MainInit
  
!!!
!!! 配列の allocate
!!!
  subroutine VariableAllocate
    !
    !初期化として, 配列を定義し, 値としてゼロを代入する.
    !
    
    !暗黙の型宣言禁止
    implicit none

    !配列割り当て
    allocate( xyz_ZetaBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ZetaNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ZetaAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DZetaDtNl(imin:imax,jmin:jmax,kmin:kmax) )

    !ゼロを代入
    xyz_ZetaBl    = 0.0d0
    xyz_ZetaNl    = 0.0d0
    xyz_ZetaAl    = 0.0d0
    xyz_DZetaDtNl = 0.0d0
    
  end subroutine VariableAllocate

!!!
!!! 時間フィルター. リープフロッグで解く時に必要. 
!!!
  
  subroutine AsselinTimeFilter
    !
    ! 時間フィルター; Asselin のタイムフィルターを利用    
    !   t = 0.0 の場合には tfil = 0.0d0, それ以外は tfil = 1.0d-1
    !   (t = 0 の時はオイラー法で, それ以外はリープフロッグ法で積分するため)
    !
    use TimeSet, only: tfil

    ! 暗黙の型宣言禁止    
    !
    implicit none

    ! 作業変数
    !
    real(DP) :: tfil2

    tfil2 = 1.0d0 - 2.0d0 * tfil 
    xyz_ZetaNl = tfil * ( xyz_ZetaBl + xyz_ZetaAl ) + tfil2 * xyz_ZetaNl

  end subroutine AsselinTimeFilter
  

!!!
!!! src/io 以下で行っていることの簡略化. gtool5 を利用するだけ. 
!!! deepconv の初期値・リスタートファイルは HistoryPut,
!!! deepconv の結果は HistoryAutoPut を使って出力している.
!!!
!!! gtool5 : http://www.gfd-dennou.org/arch/gtool/
!!!  
  
  ! I/O の初期化. 
  ! HistoryPut, HistoryAddVarible を利用. 
  !
  subroutine fileset_init
    
    call HistoryCreate( &                                  ! ヒストリー作成
      file='advect.nc', title='2D diffusion model',                             &
      source='Sample program of deepconv/arare6',                                &
      institution='GFD_Dennou Club deepconv project',                            &
      dims=(/'x','y','z','t'/), dimsizes=(/nx,ny,nz,0/),                         &
      longnames=(/'X-coordinate','Y-coordinate','Z-coordinate','time        '/), &
      units=(/'m','m','m','1'/),                                                 &
      origin=0.0, interval=real(NstepOutput * DelTimeLong) )
    
    call HistoryPut('x',x_X(1:nx))               ! 変数出力
    call HistoryPut('y',y_Y(1:ny))               ! 変数出力
    call HistoryPut('z',z_Z(1:nz))               ! 変数出力

    call HistoryAddVariable( &                   ! 変数定義
         varname='zeta', dims=(/'x','y','z','t'/), & 
         longname='variable', units='1', xtype='float')       !出力時は float に.

  end subroutine fileset_init

end program advect
