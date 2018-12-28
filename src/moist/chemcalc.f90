!= 化学関連の諸量を計算するためのモジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: chemcalc.f90,v 1.12 2014/07/08 01:05:32 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006-2014. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module ChemCalc
  !
  != 化学関連の諸量を計算するためのモジュール. 
  !
  ! AMP と Antoine の飽和蒸気圧式を用いて以下を求める. 
  ! デフォルトでは AMP 式を使うようにしてある. 
  !  * 飽和蒸気圧
  !  * 飽和蒸気圧の温度微分
  !  * 潜熱
  !

  !モジュール呼び出し
  use dc_types,   only: DP              !精度指定
  use ChemData,   only: ChemData_SpcNum !データベース上の化学種数

  ! 暗黙の型宣言禁止
  !
  implicit none

  ! 変数の定義
  !
  real(DP), save, public  :: ReactHeatNH4SH       !NH4SH 生成反応熱 [J/K kg]
  real(DP), save, public  :: ReactHeatNH4SHPerMol !NH4SH 生成反応熱 [J/K mol]

  integer,  save, private :: a_kmin(ChemData_SpcNum)  !物質毎に決めた配列の下限
  integer,  save, private :: a_kmax(ChemData_SpcNum)  !物質毎に決めた配列の下限
  real(DP), save, private :: a_SwAmp(ChemData_SpcNum) !スイッチ. AMP を使う場合は 1.0, そうでなければ 0.0  
  real(DP), save, private :: a_SwAnt(ChemData_SpcNum) !スイッチ. Antoine を使う場合は 1.0, そうでなければ 0.0
  real(DP), save, private :: a_antA(ChemData_SpcNum)  !Antoine の蒸気圧式の A 係数
  real(DP), save, private :: a_antB(ChemData_SpcNum)  !Antoine の蒸気圧式の B 係数
  real(DP), save, private :: a_antC(ChemData_SpcNum)  !Antoine の蒸気圧式の C 係数
  real(DP), save, private :: a_antU(ChemData_SpcNum)  !Antoine の蒸気圧式の単位換算のための係数
  real(DP), save, private :: a_ampA(ChemData_SpcNum)  !AMP 式の蒸気圧式の A 係数
  real(DP), save, private :: a_ampB(ChemData_SpcNum)  !AMP 式の蒸気圧式の B 係数
  real(DP), save, private :: a_ampC(ChemData_SpcNum)  !AMP 式の蒸気圧式の C 係数
  real(DP), save, private :: a_ampD(ChemData_SpcNum)  !AMP 式の蒸気圧式の D 係数
  real(DP), save, private :: a_ampE(ChemData_SpcNum)  !AMP 式の蒸気圧式の E 係数
  real(DP), save, private :: a_MolWt(ChemData_SpcNum) !分子量

  ! 公開するサブルーチンに public 属性を付ける
  !
  public ChemCalc_Init                 !初期化ルーチン
  public MolWt                         !分子量
  public GasR                          !気体定数
  public CpRef, CpPerMolRef, CvRef     !定圧比熱, 定積比熱
  public SvapPress, xyz_SvapPress      !飽和蒸気圧 [Pa]
  public xyz_LatentHeat                !潜熱 [J/K kg]
  public LatentHeatPerMol              !潜熱 [J/K mol]
  public xyz_DQMixSatDPTemp
  public xyz_DelQMixNH4SH
  public DelMolFrNH4SH
 
contains
  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ChemCalc_Init
    !
    !初期化ルーチン
    !

    !モジュール呼び出し
    use dc_types,   only: DP                               !精度指定
    use dc_message, only: MessageNotify                    !メッセージ表示
    use ChemData,   only: ChemData_init,                 & !初期化
      &                   ChemData_SpcNum,               & !データベース上の化学種数
      &                   ChemData_SvapPress_AntoineA,   & !Antoine 式の A 係数
      &                   ChemData_SvapPress_AntoineB,   & !Antoine 式の B 係数
      &                   ChemData_SvapPress_AntoineC,   & !Antoine 式の C 係数
      &                   ChemData_SvapPress_AntoineUnit,& !単位変換用係数
      &                   ChemData_SvapPress_AMPA,       & !AMP 式の A 係数
      &                   ChemData_SvapPress_AMPB,       & !AMP 式の B 係数
      &                   ChemData_SvapPress_AMPC,       & !AMP 式の C 係数
      &                   ChemData_SvapPress_AMPD,       & !AMP 式の D 係数
      &                   ChemData_SvapPress_AMPE,       & !AMP 式の E 係数
      &                   ChemData_MolWt,                & !分子量
      &                   GasRUniv,                      & !気体定数
      &                   ChemData_SpcSymbol,            & !分子名
      &                   ChemData_OneSpcID                !化学種の ID 検索     
    use gridset,    only: nz,                            & ! 物理領域の大きさ
      &                   kmin, kmax                       ! 配列の Z 方向の上限・下限
    use constants,  only: PressSfc                         !下部境界での圧力 [Pa]
    use basicset,   only: xyz_TempBZ,                    & !温度の基本場
      &                   xyz_PressBZ                      !圧力の基本場
    use axesset,    only: z_Z                              !z 軸
    use namelist_util, &
      &             only: namelist_filename
    use dc_iounit,  only: FileOpen

    !暗黙の型宣言禁止
    implicit none

    !内部変数
    character(20)      :: Name
    integer            :: id
    integer            :: k
    integer            :: unit, ierr
    real(DP)           :: Temp
    real(DP)           :: Press
    real(DP),parameter :: Temp0C = 273.15d0
    real(DP)           :: logsvap
    real(DP)           :: HeightUp = 0.0d0
    real(DP)           :: HeightDown = 0.0d0

    !NAMELIST の定義
    NAMELIST /chemcalc_nml/ HeightUp, HeightDown
 
    !ファイルオープン. 情報取得. 
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=chemcalc_nml, iostat=ierr, err=99)
    close(unit)
99  call MessageNotify( "M", "ChemCalc_Init", "No information of chemcalc_nml in config file; use default values")

    !-----------------------------------------------------------
    ! 初期化
    !

    ! データベースの初期化
    call chemdata_init

    !Antoine の飽和蒸気圧式の係数
    a_antA = ChemData_SvapPress_AntoineA
    a_antB = ChemData_SvapPress_AntoineB
    a_antC = ChemData_SvapPress_AntoineC
    a_antU = ChemData_SvapPress_AntoineUnit

    !AMP 式の飽和蒸気圧式の係数
    a_ampA = ChemData_SvapPress_AMPA
    a_ampB = ChemData_SvapPress_AMPB
    a_ampC = ChemData_SvapPress_AMPC
    a_ampD = ChemData_SvapPress_AMPD
    a_ampE = ChemData_SvapPress_AMPE

    !分子量
    a_MolWt = ChemData_MolWt
    
    !NH4SH の反応熱の初期化
    !  NH4SH 1kg に対する反応熱にする.
    Name = 'NH4SH-s'
    id   = ChemData_OneSpcID( Name )  
    
    ReactHeatNH4SHPerMol  = GasRUniv * 10834.0d0
    ReactHeatNH4SH = GasRUniv * 10834.0d0 / MolWt( id )

    !--------------------------------------------------------
    ! 物質によって, AMP を使うか Antoine を使うか決める. 
    ! AMP の係数がゼロならば Antoine を使うことにする. 
    !     

    do ID = 1, ChemData_SpcNum
      if ( a_ampA(ID) /= 0.0d0 ) then 

        ! AMP の係数が与えられている場合
        !
        a_SwAmp(ID) = 1.0d0
        a_SwAnt(ID) = 0.0d0

      elseif ( a_antA(ID) /= 0.0d0 ) then 

        ! Antoine の係数だけ与えられている場合
        !
        a_SwAmp(ID) = 0.0d0
        a_SwAnt(ID) = 1.0d0

      else

        ! 気体の場合
        !
        a_SwAmp(ID) = 0.0d0
        a_SwAnt(ID) = 0.0d0

      end if
    end do

    ! 初期化 (物理領域の上限を与える) 
    a_kmax = nz

    ! 初期化 (物理領域の下限を与える) 
    a_kmin = 1


    !!--------------------------------------------------------
    !! TempBZ, PressBZ が allocated されて値が確定している場合は, 
    !! 計算をサボるための以下の処置を行う
    !!
    if ( allocated( xyz_TempBZ) ) then 

    !--------------------------------------------------------
    ! 計算を適当にサボるための処置 (1) 
    !
    ! 飽和蒸気圧が十分小さければ, 計算する必要はないとして良いのだろうか?
    ! 木星計算の場合, 対流圏上部では水の蒸気圧はほとんどゼロだが, 移流
    ! によって雨が圏界面付近まで持ち上がり, 落下しながら蒸発する.
    ! 先天的に高度を指定するのは難しそうなので, 高度 (HeightUp) が指定された
    ! 場合には, それより上空の計算は行わないということに. 

    ! HeightUp が設定されている場合に処理を行う. 
    ! 物質毎に違うということはありえないので, 物質に対するループは回さない. 
    !
    if ( HeightUp > 0.0d0 ) then 

      do k = kmin, kmax
        if ( z_Z(k) > HeightUp ) then 
          a_kmax = k
          exit
        end if
      end do
     
    end if

    !--------------------------------------------------------
    ! 計算を適当にサボるための処置 (2) 
    ! 
    ! 飽和蒸気圧を計算するための配列添字の下限を決める. 
    ! HeightDown が設定されている場合にはそれを優先し, 
    ! HeightDown が負の場合には以下の手続きで下限設定する. 
    !
    ! * 飽和蒸気圧が下部境界での圧力を上回るはずがないということを基準に決める.
    !   * どの物質に対しても Antoine の係数は与えられていることが前提. 
    ! 
    ! Fujitsu Fortran では, exp(logsvap) [logsvap > 700] でエラーが出る. 
    ! HeightUp が正であっても, logsvap > 700 となる高度を下限とする. 

    do ID = 1, ChemData_SpcNum
      ! 凝結物の場合に計算を行う (気体の場合は a_antA = 0.0)
      !
      if ( a_antA(ID) /= 0.0d0 ) then 

        ! 鉛直方向は上空からループを回す.        
        do k = nz, 1, -1

          ! 基本場の温度に対して飽和蒸気圧の log を計算
          !
          Temp  = xyz_TempBZ(1,1,k)
          Press = xyz_PressBZ(1,1,k)
          logsvap =                                    &
            &       (                                  &
            &           a_antA(ID)                     &
            &         - a_antB(ID)                     &
            &           / (a_antC(ID) + Temp - Temp0C) &
            &        ) * dlog(10.0d0)                  &
            &        + a_antU(ID)

!          write(*,*) ChemData_SpcSymbol(ID), k, dexp( logsvap ), PressSfc, Press

          ! logsvap > 700 となれば添字を保管
          ! 
          if ( logsvap > 700 ) then
            a_kmin(ID) = k
            exit

          ! 下限が指定された場合
          !
          elseif ( z_Z(k) <= HeightDown ) then 
            a_kmin(ID) = k
            exit

          ! 飽和蒸気圧より決める場合
          !
!!          elseif( HeightDown < 0.0d0 .AND. logsvap >= dlog( PressSfc ) ) then 
          elseif( logsvap >= dlog( PressSfc ) ) then 
            a_kmin(ID) = k
            exit

          end if
        end do
        
      end if
    end do
    end if

    !--------------------------------------------------------
    ! 値の確認
    !
    call MessageNotify( "M", &
      & "ChemCalc_Init", "ReactHeatNH4SH = %f", d=(/ReactHeatNH4SH/) )
    id   = ChemData_OneSpcID( Name )  
    call MessageNotify( "M", &
      & "ChemCalc_Init", "NH4SH MolWt = %f", d=(/MolWt(id)/) )
    
    do k = 1, ChemData_SpcNum
      call MessageNotify( "M", "ChemCalc_Init", &
        & "%c : a_kmin= %d, a_kmax= %d, s_SwAmp= %f, s_SwAnt= %f", &
        & i=(/a_kmin(k), a_kmax(k)/), d=(/a_SwAMP(k), a_SwAnt(k)/), c1=trim(ChemData_SpcSymbol(k)))
    end do

  end subroutine ChemCalc_Init

!!!
!!! 飽和蒸気圧, 潜熱, etc. の基本関数. 
!!! 化学種の ID と温度に対して値を返す
!!!  

!!!==========================================================================
  function CpRef(ID)
    !
    !引数で与えられた化学種に対して, 標準状態での単位質量当たりの定圧比熱を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use chemdata,   only: ChemData_CpRef !標準状態での単位質量当たりの比熱

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: CpRef        !標準状態での単位質量当たりの比熱
    integer, intent(in) :: ID           !化学種の ID

    
    !データベースから情報取得
    CpRef = ChemData_CpRef(ID)

  end function CpRef


!!!==========================================================================
  function CpPerMolRef(ID)
    !
    !引数で与えられた化学種に対して, 標準状態での単位モル当たりの定圧比熱を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP                   !精度指定
    use chemdata,   only: ChemData_CpPerMolRef !標準状態での単位モル当たりの比熱

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: CpPerMolRef  !標準状態での単位モル当たりの比熱
    integer, intent(in) :: ID           !化学種の ID

    
    !データベースから情報取得
    CpPerMolRef = ChemData_CpPerMolRef(ID)

  end function CpPerMolRef


!!!==========================================================================
  function CvRef(ID)
    !
    !引数で与えられた化学種に対して, 標準状態での単位質量当たりの定圧比熱を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use chemdata,   only: ChemData_CvRef !標準状態での単位質量当たりの比熱

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: CvRef       !標準状態での単位質量当たりの比熱
    integer, intent(in) :: ID          !化学種の ID

    
    !データベースから情報取得
    CvRef = ChemData_CvRef(ID)

  end function CvRef


!!!==========================================================================
  function MolWt(ID)
    !
    !引数で与えられた化学種に対して, 分子量を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use chemdata,   only: ChemData_MolWt !分子量

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: MolWt         !分子量
    integer, intent(in) :: ID            !化学種の ID

    
    !データベースから情報取得
    MolWt = ChemData_MolWt(ID)

  end function MolWt


!!!==========================================================================
  function GasR(ID)
    !
    !引数で与えられた化学種に対して, 気体定数を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use chemdata,   only: ChemData_GasR  !気体定数 [J/K kg]

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: GasR          !分子量
    integer, intent(in) :: ID            !化学種の ID
    
    
    !データベースから情報取得
    GasR = ChemData_GasR(ID)

  end function GasR


!!!
!!! 空間 3 次元の関数群
!!!

!!!==========================================================================
  function xyz_SvapPress( ID, xyz_Temp )
    !
    != 引数で与えられた化学種と温度に対して, 飽和蒸気圧を計算. 
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use gridset,    only: nx, ny,      & !物理領域の大きさ
      &                   imin, imax,  & !配列の X 方向の上限・下限
      &                   jmin, jmax,  & !配列の Y 方向の上限・下限
      &                   kmin, kmax     !配列の Z 方向の上限・下限
    use constants,  only: PressSfc       !下部境界での圧力       [Pa]

    !暗黙の型宣言禁止
    implicit none
    
    ! 入出力変数  
    real(DP)            :: xyz_SvapPress(imin:imax,jmin:jmax,kmin:kmax) 
                                                          !飽和蒸気圧
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)  
                                                          !温度
    integer, intent(in) :: ID                             !化学種の ID
  
    !内部変数
    real(DP)            :: LogSvapPress
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! 初期化
    ! * 飽和蒸気圧は十分大きい値にしておく.
    !
    xyz_SvapPress = PressSfc * 100.0d0

    ! 飽和蒸気圧の計算
    ! a_SwAmp, a_SwAnt を用いることで, 選択された計算方法を用いる.
    !
    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx
          
          ! 飽和蒸気圧の log を計算
          !
          LogSvapPress =                                           &
            &      (                                               &
            &         a_ampA(ID) / xyz_Temp(i,j,k)                 &
            &       + a_ampB(ID)                                   &
            &       + a_ampC(ID) * dlog( xyz_Temp(i,j,k) )         &
            &       + a_ampD(ID) * xyz_Temp(i,j,k)                 &
            &       + a_ampE(ID) * ( xyz_temp(i,j,k) ** 2 )        &
            &       + dlog(1.0d-1)                                 &
            &      ) * a_SwAmp(ID)                                 &
            &    + (                                               &
            &       + (                                            &
            &          + a_antA(ID)                                &
            &          - a_antB(ID)                                &
            &            / (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) &
            &         ) * dlog(10.0d0)                             &
            &       + a_antU(ID)                                   &
            &      ) * a_SwAnt(ID)
          
          !飽和蒸気圧を計算
          !
          xyz_SvapPress(i,j,k) =  dexp( LogSvapPress )

        end do
      end do
    end do

  end function xyz_SvapPress  

!!!==========================================================================
  function xyz_LatentHeat(ID, xyz_Temp)
    !
    != 飽和蒸気圧より潜熱を計算する. 
    !
    
    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use gridset,    only: nx, ny,      & !物理領域の大きさ
      &                   imin, imax,  & !配列の X 方向の上限・下限
      &                   jmin, jmax,  & !配列の Y 方向の上限・下限
      &                   kmin, kmax     !配列の Z 方向の上限・下限

    !暗黙の型宣言禁止
    implicit none

    !入出力変数
    real(DP)            :: xyz_LatentHeat(imin:imax,jmin:jmax,kmin:kmax)
                                                            !潜熱[J/K kg]
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
                                                    !温度[K]
    integer, intent(in) :: ID                       !化学種の ID

    !内部関数
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: GasRUniv = 8.314d0
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! 初期化
    !
    xyz_LatentHeat = 0.0d0

    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx

          ! 飽和蒸気圧の温度微分
          ! a_SwAmp, a_SwAnt を用いることで, 選択された計算方法を用いる.
          !
          DLogSvapPressDTemp =                                             &
            &    (                                                         &
            &     - a_ampA(ID) / (xyz_Temp(i,j,k) ** 2.0d0)                &
            &     + a_ampC(ID) / xyz_Temp(i,j,k)                           &
            &     + a_ampD(ID)                                             &
            &     + a_ampE(ID) * 2.0d0 * xyz_Temp(i,j,k)                   &
            &    ) * a_SwAmp(ID)                                           &
            &  + (                                                         &
            &     + a_antB(ID) * dlog(10.0d0)                              &
            &       / ( (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) ** 2.0d0 ) &
            &    ) * a_SwAnt(ID)
          
          xyz_LatentHeat(i,j,k) =                                          &
            & DLogSvapPressDTemp * GasRUniv * (xyz_Temp(i,j,k) ** 2.0d0)   &
            &  / a_MolWt(ID)

        end do
      end do
    end do
    
  end function xyz_LatentHeat


!!!==========================================================================
  function SvapPress(ID, Temp)
    !
    != 引数で与えられた化学種と温度に対して, 飽和蒸気圧を計算. 
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定

    !暗黙の型宣言禁止
    implicit none
    
    !入出力変数  
    real(DP)            :: SvapPress   !飽和蒸気圧
    real(DP),intent(in) :: Temp        !温度 [K]
    integer, intent(in) :: ID          !化学種の ID

    !内部変数
    real(DP)            :: LogSvapPress
    real(DP), parameter :: Temp0C = 273.15d0

    ! 飽和蒸気圧の log を計算
    ! 対数が大きくなりすぎないようにする. 
    ! Fujitsu Fortran Compiler では 700 より大きい数の exp を取ると警告が出る.
    !
    LogSvapPress =                               &
      & min(                                     &
      &      (                                   &
      &         a_ampA(ID) / Temp                &
      &       + a_ampB(ID)                       &
      &       + a_ampC(ID) * dlog( Temp )        &
      &       + a_ampD(ID) * Temp                &
      &       + a_ampE(ID) * ( Temp ** 2 )       &
      &       + dlog(1.0d-1)                     &
      &      ) * a_SwAmp(ID)                     &
      &    + (                                   &
      &        (                                 &
      &         + a_antA(ID)                     &
      &         - a_antB(ID)                     &
      &           / (a_antC(ID) + Temp - Temp0C) &
      &        ) * dlog(10.0d0)                  &
      &       + a_antU(ID)                       &
      &      ) * a_SwAnt(ID),                    &
      &   700.0d0                                &
      & )
          
    !飽和蒸気圧を計算
    !
    SvapPress =  dexp( LogSvapPress )

  end function SvapPress


!!!==========================================================================
  function LatentHeatPerMol(ID, Temp)
    !
    != 引数で与えられた化学種と温度に対して, 潜熱 [J/K/mol] を計算
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定

    !暗黙の型宣言禁止
    implicit none

    !入出力変数
    real(DP)            :: LatentHeatPerMol   !潜熱
    real(DP),intent(in) :: Temp               !温度
    integer, intent(in) :: ID                 !化学種名
    
    !内部変数
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: GasRUniv = 8.314d0  !普遍気体定数
    real(DP),parameter  :: Temp0C = 273.15d0

    ! 飽和蒸気圧の温度微分
    ! a_SwAmp, a_SwAnt を用いることで, 選択された計算方法を用いる.
    !
    DLogSvapPressDTemp =                                  &
      &    (                                              &
      &     - a_ampA(ID) / (Temp ** 2.0d0)                &
      &     + a_ampC(ID) / Temp                           &
      &     + a_ampD(ID)                                  &
      &     + a_ampE(ID) * 2.0d0 * Temp                   &
      &    ) * a_SwAmp(ID)                                &
      &  + (                                              &
      &     + a_antB(ID) * dlog(10.0d0)                   &
      &       / ( (a_antC(ID) + Temp - Temp0C) ** 2.0d0 ) &
      &    ) * a_SwAnt(ID)
          
    ! 潜熱の計算
    !
    LatentHeatPerMol =                                    &
      & DLogSvapPressDTemp * GasRUniv * (Temp ** 2.0d0)   

  end function LatentHeatPerMol

!!!-----------------------------------------------------------------------!!!
  function xyz_DQMixSatDPTemp(ID, MolWt, xyz_Temp, xyz_Exner)
    !
    !飽和蒸気圧の θ 微分を行う
    !実際には, dq/dp * dp/dT * dT/dθ を実行. (但し p は飽和蒸気圧)
    !
    ! * dq/dp =  Mv / (Md * p_all) 
    !   (q = p * Mv / (Md * p_all) )
    ! * dT/dθ= \pi  (T = \pi \theta)
    !
    
    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use gridset,    only: nx, ny,      & !物理領域の大きさ
      &                   imin, imax,  & !配列の X 方向の上限・下限
      &                   jmin, jmax,  & !配列の Y 方向の上限・下限
      &                   kmin, kmax     !配列の Z 方向の上限・下限
    use constants,  only: PressBasis,  & !温位の標準圧力         [Pa]
      &                   CpDry,       & !乾燥成分の平均定圧比熱 [J/K kg]
      &                   MolWtDry,    & !乾燥成分の平均分子量   [kg/mol]
      &                   GasRDry        !乾燥成分の気体定数     [J/K kg]

    !暗黙の型宣言禁止
    implicit none 
    
    !入出力変数
    integer, intent(in) :: ID
    real(DP),intent(in) :: MolWt
    real(DP),intent(in) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
                                            !温度(擾乱 + 基本場)
    real(DP),intent(in) :: xyz_Exner(imin:imax,jmin:jmax,kmin:kmax)
                                            !エクスナー関数(擾乱 + 基本場)
    real(DP)            :: xyz_DQMixSatDPTemp(imin:imax,jmin:jmax,kmin:kmax)
                           
    !内部変数
    real(DP)            :: xyz_Press(imin:imax,jmin:jmax,kmin:kmax)
                                            !圧力(擾乱 + 基本場)
    real(DP)            :: DSvapPressDTemp
                                            !飽和蒸気圧の温度微分 [Pa/K]
    real(DP)            :: LogSvapPress
    real(DP)            :: DLogSvapPressDTemp
    real(DP),parameter  :: Temp0C = 273.15d0
    integer             :: i, j, k

    ! 初期化
    !
    xyz_DQMixSatDPTemp = 0.0d0
    xyz_Press = PressBasis * (xyz_Exner ** (CpDry / GasRDry))

    ! 飽和蒸気圧の温度微分
    !
    do k = a_kmin(ID), a_kmax(ID)
      do j = 1, ny
        do i = 1, nx
          ! 飽和蒸気圧の log を計算
          !
          LogSvapPress =                                           &
            &      (                                               &
            &         a_ampA(ID) / xyz_Temp(i,j,k)                 &
            &       + a_ampB(ID)                                   &
            &       + a_ampC(ID) * dlog( xyz_Temp(i,j,k) )         &
            &       + a_ampD(ID) * xyz_Temp(i,j,k)                 &
            &       + a_ampE(ID) * ( xyz_temp(i,j,k) ** 2 )        &
            &       + dlog(1.0d-1)                                 &
            &      ) * a_SwAmp(ID)                                 &
            &    + (                                               &
            &       + (                                            &
            &          + a_antA(ID)                                &
            &          - a_antB(ID)                                &
            &            / (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) &
            &         ) * dlog(10.0d0)                             &
            &       + a_antU(ID)                                   &
            &      ) * a_SwAnt(ID)

          ! 飽和蒸気圧の温度微分
          !
          DLogSvapPressDTemp =                                             &
            &    (                                                         &
            &     - a_ampA(ID) / (xyz_Temp(i,j,k) ** 2.0d0)                &
            &     + a_ampC(ID) / xyz_Temp(i,j,k)                           &
            &     + a_ampD(ID)                                             &
            &     + a_ampE(ID) * 2.0d0 * xyz_Temp(i,j,k)                   &
            &    ) * a_SwAmp(ID)                                           &
            &  + (                                                         &
            &     + a_antB(ID) * dlog(10.0d0)                              &
            &       / ( (a_antC(ID) + xyz_Temp(i,j,k) - Temp0C) ** 2.0d0 ) &
            &    ) * a_SwAnt(ID)
      
          DSvapPressDTemp = DLogSvapPressDTemp * dexp( LogSvapPress ) 

          xyz_DQMixSatDPTemp(i,j,k) =                               &
            &   MolWt / ( MolWtDry * xyz_Press(i,j,k) )             &
            &   * DSvapPressDTemp * xyz_Exner(i,j,k)   
          
        end do
      end do
    end do
    
  end function xyz_DQMixSatDPTemp


!!!-----------------------------------------------------------------------!!!
  function xyz_DelQMixNH4SH(xyz_TempAll, xyz_PressAll, xyz_PressDry, &
    &                       xyz_QMixNH3, xyz_QMixH2S, &
    &                       MolWtNH3, MolWtH2S)
    !
    ! NH4SH 生成反応に伴う, NH4SH の生成量(混合比)を求める
    !

    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定
    use gridset,    only: imin, imax,  & ! 配列の X 方向の上限・下限
      &                   jmin, jmax,  & ! 配列の Y 方向の上限・下限
      &                   kmin, kmax     ! 配列の Z 方向の上限・下限
    use constants,  only: MolWtDry

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: xyz_TempAll(imin:imax,jmin:jmax,kmin:kmax)
                                         !温度
    real(DP),intent(in) :: xyz_PressAll(imin:imax,jmin:jmax,kmin:kmax)
                                         !圧力
    real(DP),intent(in) :: xyz_PressDry(imin:imax,jmin:jmax,kmin:kmax)
                                         !圧力
    real(DP),intent(in) :: xyz_QMixNH3(imin:imax,jmin:jmax,kmin:kmax)
                                         !NH3 の混合比
    real(DP),intent(in) :: xyz_QMixH2S(imin:imax,jmin:jmax,kmin:kmax)
                                         !H2S の混合比
    real(DP),intent(in) :: MolWtNH3      !NH3 の分子量
    real(DP),intent(in) :: MolWtH2S      !H2S の分子量

    real(DP) :: xyz_DelQMixNH4SH(imin:imax,jmin:jmax,kmin:kmax)
                                         !NH4SH の混合比
    real(DP) :: xyz_EquivConst(imin:imax,jmin:jmax,kmin:kmax)
                                         !圧平衡定数
    real(DP) :: xyzf_PPress(imin:imax,jmin:jmax,kmin:kmax,2)
                                         !作業配列(分圧)
    real(DP) :: xyz_Solution(imin:imax,jmin:jmax,kmin:kmax)
                                         !作業配列(方程式の解)

    !初期化
!    xyz_DelQMixNH4SH = 0.0d0
    
    !アンモニアと硫化水素の分圧. 
    xyzf_PPress(:,:,:,1) = xyz_QMixNH3 * xyz_PressAll * MolWtDry / MolWtNH3 
    xyzf_PPress(:,:,:,2) = xyz_QMixH2S * xyz_PressAll * MolWtDry / MolWtH2S 

    !圧平衡定数
    xyz_EquivConst = 61.781d0 - 10834.0d0 / xyz_TempAll - dlog(1.0d2)

    !気圧変化を求める. 
    !  (P_NH3 - X) * (P_H2S - X) = exp(Kp)
    !  DelX^2 - (P_NH3 + P_H2S) * DelX + P_NH3 * P_H2S - exp( Kp ) = 0
    !  という二次方程式を求める必要があるが, (P_NH3 - X) > 0 と
    !  (P_H2S - X) > 0 を満たすためには, 解の公式のうちが負の方が選択される.
    !
    xyz_Solution  =                                                       &
      & (                                                                 &
      &     sum(xyzf_PPress, 4)                                           &
      &   - dsqrt( (xyzf_PPress(:,:,:,1) - xyzf_PPress(:,:,:,2)) ** 2.0d0 &
      &            + 4.0d0 * dexp( min( 700.0d0, xyz_EquivConst ) ) )     &
      &  ) * 5.0d-1

    !生成量を求める
    xyz_DelQMixNH4SH = xyz_Solution * ( MolWtNH3 + MolWtH2S ) &
      &                   / ( xyz_PressDry * MolWtDry )

  end function xyz_DelQMixNH4SH
  

!!!-----------------------------------------------------------------------!!!
  function DelMolFrNH4SH(TempAll, PressAll, MolFrNH3, MolFrH2S, Humidity)
    !
    ! NH4SH 生成反応に伴う H2S と NH3 のモル比の減少分を求める
    !
    
    !モジュール呼び出し
    use dc_types,   only: DP             !精度指定

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP),intent(in) :: TempAll       !温度
    real(DP),intent(in) :: PressAll      !圧力
    real(DP),intent(in) :: MolFrNH3      !NH3 のモル比
    real(DP),intent(in) :: MolFrH2S      !H2S のモル比
    real(DP),intent(in) :: Humidity      !飽和比
    real(DP)            :: DelMolFrNH4SH !NH4SH 生成に伴うモル比変化
    real(DP)            :: EquivConst    !圧平衡定数
    real(DP)            :: PPress(2)     !作業配列(分圧)
    real(DP)            :: Solution      !作業配列(方程式の解)

    !------------------------------------------------------------
    !NH4SH の平衡条件
    !------------------------------------------------------------
    !アンモニアと硫化水素の分圧
    PPress(1) = MolFrNH3 * PressAll
    PPress(2) = MolFrH2S * PressAll

    !圧平衡定数
    EquivConst = 61.781d0 - 10834.0d0 / TempAll - dlog(1.0d2) - 2.0d0 * dlog( Humidity )
    
    !気圧変化を二次方程式の解として求める. 
    Solution = 5.0d-1 * (sum(PPress)                                        &
      &        - dsqrt( (PPress(1) - PPress(2))**2.0d0                      &
      &                    + 4.0d0 * dexp( min( 700.0d0, EquivConst ))) )
    
    !NH4SH の生成量. 
    DelMolFrNH4SH = Solution / PressAll

  end function DelMolFrNH4SH

    
end module ChemCalc
