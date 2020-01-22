! Sample program 
! 
!
program soundwave

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
    &                         nx, ny, nz, ncmax, FlagCalc3D
  use timeset,         only : timeset_init,                            &
    &                         TimesetDelTimeHalf, TimesetProgress,     &
    &                         TimeA, TimeN, TimeB,                     &
    &                         Nstep, NstepShort, NstepLong,            &
    &                         NstepOutput, FlagInitialRun,             &
    &                         DelTimeLong, DelTimeShort
  use axesset,         only : axesset_init, x_X, y_Y, z_Z, XMax, ZMax, &
    &                         dx, dy, dz
  use constants,       only : constants_init, CpDry
  use namelist_util,   only : NmlutilInit
  
  ! 下請けモジュール
  !
  use setmargin,       only : SetMargin_init, & 
    &                         SetMargin_pyz, SetMargin_xqz, SetMargin_xyr, SetMargin_xyz
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
  real(DP), allocatable :: pyz_VelXAs(:,:,:)
  real(DP), allocatable :: pyz_VelXNs(:,:,:)
  real(DP), allocatable :: pyz_VelXAl(:,:,:)
  real(DP), allocatable :: pyz_VelXNl(:,:,:)
  real(DP), allocatable :: pyz_VelXBl(:,:,:)
  real(DP), allocatable :: xqz_VelYAs(:,:,:)
  real(DP), allocatable :: xqz_VelYNs(:,:,:)
  real(DP), allocatable :: xqz_VelYAl(:,:,:)
  real(DP), allocatable :: xqz_VelYNl(:,:,:)
  real(DP), allocatable :: xqz_VelYBl(:,:,:)
  real(DP), allocatable :: xyr_VelZAs(:,:,:)
  real(DP), allocatable :: xyr_VelZNs(:,:,:)
  real(DP), allocatable :: xyr_VelZAl(:,:,:)
  real(DP), allocatable :: xyr_VelZNl(:,:,:)
  real(DP), allocatable :: xyr_VelZBl(:,:,:)
  real(DP), allocatable :: xyz_ExnerAs(:,:,:)
  real(DP), allocatable :: xyz_ExnerNs(:,:,:)
  real(DP), allocatable :: xyz_ExnerAl(:,:,:)
  real(DP), allocatable :: xyz_ExnerNl(:,:,:)
  real(DP), allocatable :: xyz_ExnerBl(:,:,:)
  real(DP), allocatable :: xyz_PTempAs(:,:,:)
  real(DP), allocatable :: xyz_PTempNs(:,:,:)
  real(DP), allocatable :: xyz_PTempAl(:,:,:)
  real(DP), allocatable :: xyz_PTempNl(:,:,:)
  real(DP), allocatable :: xyz_PTempBl(:,:,:)

  real(DP), allocatable :: pyz_DVelXDtNs(:,:,:)
  real(DP), allocatable :: pyz_DVelXDtNl(:,:,:)
  real(DP), allocatable :: xqz_DVelYDtNs(:,:,:)
  real(DP), allocatable :: xqz_DVelYDtNl(:,:,:)
  real(DP), allocatable :: xyr_DVelZDtNs(:,:,:)
  real(DP), allocatable :: xyr_DVelZDtNl(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNs(:,:,:)
  real(DP), allocatable :: xyz_DExnerDtNl(:,:,:)
  real(DP), allocatable :: xyz_DPTempDtNs(:,:,:)
  real(DP), allocatable :: xyz_DPTempDtNl(:,:,:)

  real(DP), allocatable :: xyz_VelDivNs(:,:,:)

  real(DP), parameter :: AlphaSound = 5.0d-2
  real(DP)            :: AlphaH
  real(DP)            :: AlphaV  

  real(DP), parameter :: xyz_PTempBZ = 300.0d0      !将来的に配列に
  real(DP), parameter :: xyz_DensBZ = 1.1627d0      !将来的に配列に
  real(DP), parameter :: xyz_VelSoundBZ = 346.96d0  !将来的に配列に

  real(DP), parameter :: xyr_PTempBZ = 300.0d0      !将来的に配列に
  real(DP), parameter :: xyr_DensBZ = 1.1627d0      !将来的に配列に  
  
  ! 物理パラメタ
  ! 設定ファイル (.conf) では指定できない項目を陽に指定.  
  ! 
!  real(DP), parameter :: nu    = 1.0d-1  ! 拡散係数
  real(DP), parameter :: VelX0 = 1.0d0   ! 移流速度 X 方向
  real(DP), parameter :: VelZ0 = 0.0d0   ! 移流速度 Z 方向
  real(DP), parameter :: DelMax = 5.0d-2 ! 山の高さ
  real(DP), parameter :: Xc = 5.0d3      ! 山の中心位置
  real(DP), parameter :: Zc = 5.0d3      ! 山の中心位置
  real(DP), parameter :: Xr = 4.0d2      ! 山の幅
  real(DP), parameter :: Zr = 4.0d2      ! 山の幅

  ! I/O
  !
  character(STRING) :: cfgfile ! NAMELIST ファイル名 ; NAMELIST fine name

  ! do ループ変数 ; do loop variable
  !
  integer :: it = 0
  integer :: i, k, tau


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
        xyz_ExnerNl(i,:,k) =                                   & 
             &    DelMax                                       &
             &  * dexp(                                        &
             &       - ( (x_X(i) - Xc) / Xr )**2.0d0 * 5.0d-1  &
             &       - ( (z_Z(k) - Zc) / Zr )**2.0d0 * 5.0d-1  &
             &    )        
    end do
  end do

  ! 境界条件
  !
  call SetMargin_xyz( xyz_ExnerNl )

  ! 初期化
  !
  xyz_ExnerBl = xyz_ExnerNl

!  write(*,*) xyz_ExnerNl(1,1,:)
  

!!!
!!! 音波減衰項
!!!
  if ( FlagCalc3D ) then
     AlphaH = AlphaSound * ( Min( dx * dx, dy * dy) ) / DelTimeShort
     AlphaV = AlphaSound * ( Min( dx * dx, dy * dy, dz * dz) ) / DelTimeShort
  else
     AlphaH = AlphaSound * ( dx * dx ) / DelTimeShort
     AlphaV = AlphaSound * ( Min( dx * dx, dz * dz) ) / DelTimeShort
  end if

!  write(*,*) AlphaH
!  write(*,*) AlphaV
 
!!!
!!! 時間積分
!!!
  do while ( Nstep <= NstepLong )

     ! tendency (傾き) の計算.
     ! 計算の安定性のためには数値拡散項が必要だが今は無視.
     ! 長時間計算できないが今はよしとする.
     ! デフォルトでは VelZ0 = 0 としている.
     ! 
!     xyz_DZetaDtNl =                                     &
!          &  - VelX0 * xyz_pyz(pyz_dx_xyz(xyz_ZetaNl))   &
!          &  - VelZ0 * xyz_xyr(xyr_dz_xyz(xyz_ZetaNl))   
         
     ! 温位の積分 (leap-frog)
     !
!     xyz_ZetaAl = xyz_ZetaBl + 2.0d0 * DelTimeLong * xyz_DZetaDtNl
          
     ! 短い時間ステップの初期値作成.
     ! Initial values set up for time integration with short time step.
     !
     pyz_VelXNs  = pyz_VelXBl
     xqz_VelYNs  = xqz_VelYBl
     xyr_VelZNs  = xyr_VelZBl
     xyz_ExnerNs = xyz_ExnerBl
     xyz_PTempNs = xyz_PTempBl
     
     ! 短い時間ステップの時間積分. オイラー法を利用.
     ! Time integration with short time step.
     !
     Euler: do tau = 1, NstepShort
        
        ! 速度の発散
        ! divergence of velocity
        !
        xyz_VelDivNs =                       &
             &    xyz_dx_pyz( pyz_VelXNs )   &
             &  + xyz_dy_xqz( xqz_VelYNs )   &
             &  + xyz_dz_xyr( xyr_VelZNs )
       
        ! 圧力項の計算. 陽解放 (HE-VE)
        ! pressure terms 
        !
        pyz_DVelXDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * pyz_dx_xyz( xyz_ExnerNs ) &
             &   + AlphaH * pyz_dx_xyz( xyz_VelDivNs )

        xqz_DVelYDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * xqz_dy_xyz( xyz_ExnerNs ) &
             &   + AlphaH * xqz_dy_xyz( xyz_VelDivNs )

        xyr_DVelZDtNs =                                            &
             &   - CpDry * xyz_PTempBZ * xyr_dz_xyz( xyz_ExnerNs ) &
             &   + AlphaV * xyr_dz_xyz( xyz_VelDivNs )

        xyz_DExnerDtNs =                                                    &
             & - xyz_VelSoundBZ * xyz_VelSoundBZ / (CpDry * xyz_PTempBZ )   &
             &   * (  + xyz_dx_pyz( pyz_VelXNs )                            &
             &        + xyz_dy_xqz( xqz_VelYNs )                            &
             &        + xyz_dz_xyr( xyr_VelZNs )                            & 
             !             &        + xyz_dz_xyr( xyr_DensBZ * xyr_PTempBZ * xyr_VelZNs ) &
!             &           / (xyz_DensBZ * xyz_PTempBZ )                      &
             &     ) 
        
        ! 積分
        ! time integration
        !
        pyz_VelXAs = pyz_VelXNs  + DelTimeShort * ( pyz_DVelXDtNl  + pyz_DVelXDtNs  )
        xqz_VelYAs = xqz_VelYNs  + DelTimeShort * ( xqz_DVelYDtNl  + xqz_DVelYDtNs  )
        xyr_VelZAs = xyr_VelZNs  + DelTimeShort * ( xyr_DVelZDtNl  + xyr_DVelZDtNs  )
        xyz_ExnerAs= xyz_ExnerNs + DelTimeShort * ( xyz_DExnerDtNl + xyz_DExnerDtNs )

        ! 境界条件
        ! boundary condition
        !
        call SetMargin_pyz( pyz_VelXAs  )
        call SetMargin_xqz( xqz_VelYAs  )
        call SetMargin_xyr( xyr_VelZAs  )
        call SetMargin_xyz( xyz_ExnerAs )
        
        ! 短い時間ステップのループを回すための処置
        ! Renew prognostic variables for next short time step integration.
        !
        xyz_ExnerNs  = xyz_ExnerAs
        pyz_VelXNs   = pyz_VelXAs
        xqz_VelYNs   = xqz_VelYAs
        xyr_VelZNs   = xyr_VelZAs
        
     end do Euler
     
     ! 最終的な短い時間ステップでの値を長い時間ステップでの値とみなす
     ! Renew prognostic variables for next long time step integration.
     !
     xyz_ExnerAl  = xyz_ExnerAs
     pyz_VelXAl   = pyz_VelXAs
     xqz_VelYAl   = xqz_VelYAs
     xyr_VelZAl   = xyr_VelZAs
         
     ! Asselin のタイムフィルターを利用
     !
     call AsselinTimeFilter
     
     ! 境界条件 ; Boundary condition
     !
     call SetMargin_xyz( xyz_PTempAl )

     ! ヒストリファイル出力. 
     !
     if ( it == NstepOutput ) then
        call HistoryPut('VelX',  pyz_VelXNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('VelY',  xqz_VelYNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('VelZ',  xyr_VelZNl( 1:nx,1:ny,1:nz ))
        call HistoryPut('Exner', xyz_ExnerNl(1:nx,1:ny,1:nz ))
        call HistoryPut('PTemp', xyz_PTempNl(1:nx,1:ny,1:nz ))
        it = 0
     end if
     
     ! ループを回すための処理
     !
     pyz_VelXNl  = pyz_VelXAl
     pyz_VelXBl  = pyz_VelXNl
     xqz_VelYNl  = xqz_VelYAl
     xqz_VelYBl  = xqz_VelYNl
     xyr_VelZNl  = xyr_VelZAl
     xyr_VelZBl  = xyr_VelZNl
     xyz_ExnerNl = xyz_ExnerAl
     xyz_ExnerBl = xyz_ExnerNl
     xyz_PTempNl = xyz_PTempAl
     xyz_PTempBl = xyz_PTempNl
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

    ! 定数の情報の初期化
    !
    call constants_init
    
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
    allocate( pyz_VelXAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_VelXBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_VelYBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_VelZBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_ExnerBl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempAl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_PTempBl(imin:imax,jmin:jmax,kmin:kmax) )
    
    allocate( pyz_DVelXDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( pyz_DVelXDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_DVelYDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xqz_DVelYDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DVelZDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyr_DVelZDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DExnerDtNl(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DPTempDtNs(imin:imax,jmin:jmax,kmin:kmax) )
    allocate( xyz_DPTempDtNl(imin:imax,jmin:jmax,kmin:kmax) )

    allocate( xyz_VelDivNs(imin:imax,jmin:jmax,kmin:kmax) )
    
    !ゼロを代入
    pyz_VelXAs = 0.0d0; pyz_VelXNs = 0.0d0; pyz_VelXAl = 0.0d0; pyz_VelXNl = 0.0d0; pyz_VelXBl = 0.0d0
    xqz_VelYAs = 0.0d0; xqz_VelYNs = 0.0d0; xqz_VelYAl = 0.0d0; xqz_VelYNl = 0.0d0; xqz_VelYBl = 0.0d0
    xyr_VelZAs = 0.0d0; xyr_VelZNs = 0.0d0; xyr_VelZAl = 0.0d0; xyr_VelZNl = 0.0d0; xyr_VelZBl = 0.0d0
    xyz_ExnerAs= 0.0d0; xyz_ExnerNs= 0.0d0; xyz_ExnerAl= 0.0d0; xyz_ExnerNl= 0.0d0; xyz_ExnerBl= 0.0d0
    xyz_PTempAs= 0.0d0; xyz_PTempNs= 0.0d0; xyz_PTempAl= 0.0d0; xyz_PTempNl= 0.0d0; xyz_PTempBl= 0.0d0
    
    pyz_DVelXDtNs = 0.0d0;    pyz_DVelXDtNl = 0.0d0
    xqz_DVelYDtNs = 0.0d0;    xqz_DVelYDtNl = 0.0d0
    xyr_DVelZDtNs = 0.0d0;    xyr_DVelZDtNl = 0.0d0
    xyz_DExnerDtNs = 0.0d0;   xyz_DExnerDtNl = 0.0d0
    xyz_DPTempDtNs = 0.0d0;   xyz_DPTempDtNl = 0.0d0

    xyz_VelDivNs = 0.0d0
     
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
    xyz_PTempNl = tfil * ( xyz_PTempBl + xyz_PTempAl ) + tfil2 * xyz_PTempNl

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
      file='soundwave.nc', title='2D diffusion model',                           &
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
         varname='VelX', dims=(/'x','y','z','t'/), & 
         longname='Velocity (x-direction)', units='m/s', xtype='float')       !出力時は float に.
    call HistoryAddVariable( &                   ! 変数定義
         varname='VelY', dims=(/'x','y','z','t'/), & 
         longname='Velocity (y-direction)', units='m/s', xtype='float')       !出力時は float に.
    call HistoryAddVariable( &                   ! 変数定義
         varname='VelZ', dims=(/'x','y','z','t'/), & 
         longname='Velocity (z-direction)', units='m/s', xtype='float')       !出力時は float に.
    call HistoryAddVariable( &                   ! 変数定義
         varname='Exner', dims=(/'x','y','z','t'/), & 
         longname='Exner function', units='1', xtype='float')       !出力時は float に.
    call HistoryAddVariable( &                   ! 変数定義
         varname='PTemp', dims=(/'x','y','z','t'/), & 
         longname='Potential temperature', units='K', xtype='float')       !出力時は float に.

  end subroutine fileset_init

end program soundwave
