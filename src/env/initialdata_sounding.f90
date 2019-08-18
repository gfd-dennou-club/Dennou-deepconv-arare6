! サウンディングファイルから初期場を決めるためのサブルーチン
!
! Authors::   SUGIYAMA Koichiro, ODAKA Masatsugu
! Version::   $Id: initialdata_sounding.f90,v 1.8 2014/07/08 00:59:09 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module initialdata_sounding
  !
  ! サウンディングファイルから初期場を決めるためのサブルーチン
  ! サウンディングファイルはテキストファイルで書かれている. 
  ! 将来的には netCDF に変更する予定.
  !

  !モジュール読み込み
  use dc_types, only: DP, STRING
  use gtool_history, only : GT_HISTORY
  
  !暗黙の型宣言禁止
  implicit none

  !共通変数
  real(DP), save, private :: r_tmpAlt(10000)
  real(DP), save, private :: r_tmpTemp(10000)
  real(DP), save, private :: r_tmpPTemp(10000)
  real(DP), save, private :: r_tmpVelX(10000)
  real(DP), save, private :: r_tmpVelY(10000)
  real(DP), save, private :: TempTr = 0.0d0
  real(DP), save, private :: AltTr  = 0.0d0
  real(DP), save, private :: DelAlt = 4.0d3
  integer,  save, private :: NumRec = 0

  integer,  save, private :: AltCol   = 0
  integer,  save, private :: TempCol  = 0
  integer,  save, private :: PTempCol = 0
  integer,  save, private :: VelXCol  = 0
  integer,  save, private :: VelYCol  = 0

  integer , save, private :: NumVIRA = 23
  real(DP), save, private :: r_viraTemp(23)
  real(DP), save, private :: r_viraAlt(23)

  type(GT_HISTORY),  save, private :: rstat
  character(STRING), save, private :: OutputFile  = "BasicZ.Sounding.nc"

  !公開変数
  public  initialdata_sounding_init
  public  initialdata_sounding_basic
  public  initialdata_sounding_wind

contains

!!!------------------------------------------------------------------------------!!!
  subroutine initialdata_sounding_init
    !
    !ファイルからサウンディングデータを読み込む
    !

    !モジュール読み込み
    use dc_types,      only: DP
    use dc_iounit,     only: FileOpen      
    use namelist_util, only: namelist_filename
    
    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    integer             :: unit         !設定ファイル用装置番号  
    integer             :: io
    character(30)       :: SoundingFile    
    
    integer, parameter  :: maxch=12
    character(len=100)  :: buf, eachcol(maxch)
    integer             :: num, MaxCol
    
    !設定ファイルから読み込む出力ファイル情報
    !
    NAMELIST /initialdata_sounding_nml/ SoundingFile, AltCol, TempCol, PTempCol, VelXCol, VelYCol, TempTr, AltTr, DelAlt    

    !設定ファイルから読み込み
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_sounding_nml)
    close(unit)

    ! 初期化
    !
    io = 0
    NumRec = 0
    MaxCol = max( AltCol, TempCol )
    r_tmpAlt  = 0.0d0; r_tmpTemp = 0.0d0
    r_tmpVelX = 0.0d0; r_tmpVelY = 0.0d0; r_tmpPTemp = 0.0d0; 

    ! ファイルのオープン
    !
    open (17, file=SoundingFile, status='old')
    
    ! ファイル呼び出し
    !
    do while ( io == 0 ) 
      ! 1 行分読み出し
      !
      read (17, '(a)', IOSTAT=io) buf
      
!!      write(*,*) buf

      ! 行をカンマ区切りで分割
      !
      call devidecsv( buf, eachcol, maxch, num )

      ! num の値が小さいものはヘッダとみなす. 
      !
      if (num >= MaxCol) then 
        ! 行数の計算
        !
        NumRec = NumRec + 1

        ! 値の代入
        !
        if (AltCol > 0)   read( eachcol(AltCol)(1:len_trim(eachcol(AltCol))), *)     r_tmpAlt(NumRec) 
        if (TempCol > 0)  read( eachcol(TempCol)(1:len_trim(eachcol(TempCol))), *)   r_tmpTemp(NumRec)   
        if (PTempCol > 0) read( eachcol(PTempCol)(1:len_trim(eachcol(PTempCol))), *) r_tmpPTemp(NumRec)   
        if (VelXCol > 0)  read( eachcol(VelXCol)(1:len_trim(eachcol(VelXCol))), *)   r_tmpVelX(NumRec)   
        if (VelYCol > 0)  read( eachcol(VelYCol)(1:len_trim(eachcol(VelYCol))), *)   r_tmpVelY(NumRec) 

      end if

    end do

!!    write(*,*) r_tmpAlt(1:100)
!!    write(*,*) r_tmpTemp(1:100)

    ! ファイルのクローズ
    !
    close (17)    

    ! VIRA データ
    ! 本来はファイルから与えるべきだがテストとして. 
    !
    r_viraTemp = (/241.0d0, 235.4d0, 229.8d0, 224.1d0, 218.6d0, &
      &            212.1d0, 205.3d0, 197.1d0, 189.9d0, 183.8d0, &
      &            178.2d0, 173.6d0, 169.4d0, 167.2d0, 167.2d0, &
      &            169.2d0, 172.0d0, 175.4d0, 178.0d0, 189.0d0, &
      &            206.0d0, 225.0d0, 246.5d0 /)
    r_viraAlt = (/ 66000.0d0,  68000.0d0,  70000.0d0,  72000.0d0,  74000.0d0, &
      &            76000.0d0,  78000.0d0,  80000.0d0,  82000.0d0,  84000.0d0, &
      &            86000.0d0,  88000.0d0,  90000.0d0,  92000.0d0,  94000.0d0, &
      &            96000.0d0,  98000.0d0, 100000.0d0, 110000.0d0, 120000.0d0, &
      &           130000.0d0, 140000.0d0, 150000.0d0 /)


  end subroutine initialdata_sounding_init


!!!------------------------------------------------------------------------------!!!
  subroutine  initialdata_sounding_basic( z_Temp, z_Press )
    !
    ! サウンディングデータから基本場を作成する
    !

    !モジュール読み込み
    use dc_types,  only: DP
    use gridset,   only: kmin, kmax,       &!配列サイズ (Z 方向)
      &                  nz                 !物理領域の大きさ (Z 方向)
    use axesset,   only: r_Z, z_Z,         &!高度
      &                  dz                 !格子間隔 (Z 方向)
    use constants, only: Grav,             &
      &                  GasRDry,          &
      &                  PressSfc, PressBasis, &
      &                  CpDry
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out):: z_Press(kmin:kmax)           !圧力
    real(DP), intent(out):: z_Temp(kmin:kmax)            !温度
    real(DP)             :: r_Press(kmin:kmax)           !圧力
    real(DP)             :: r_Temp(kmin:kmax)            !温度
    real(DP)             :: z_PTemp(kmin:kmax)           !温位
    real(DP)             :: r_PTemp(kmin:kmax)           !温位
    real(DP)             :: z_Stab(kmin:kmax)            !安定度
    real(DP)             :: Dens                         !密度
    real(DP)             :: z_DTempDZ(kmin:kmax)
    real(DP)             :: DTempDZ
    real(DP)             :: ratio
    integer              :: i, k, k1, k2, kmax2

    ! 初期化
    !
    z_PTemp = 0.0d0
    r_PTemp = 0.0d0
    z_Temp  = 0.0d0
    r_Temp  = 0.0d0
    z_Press = 0.0d0
    r_Press = 0.0d0
    z_Stab  = 0.0d0

    !!
    !! 線形補完
    !!
    do k = kmin, kmax
      do i = 1, NumRec - 1 
        if ( r_tmpAlt(i) == r_Z(k) ) then 
          if (PTempCol > 0) r_PTemp(k) = r_tmpPTemp(i)
          if (TempCol  > 0) r_Temp(k)  = r_tmpTemp(i) 
        elseif ( r_tmpAlt(i) < r_Z(k) .AND. r_Z(k) <= r_tmpAlt(i+1) ) then 
          ratio = ( r_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          if (PTempCol > 0) r_PTemp(k) = r_tmpPTemp(i) + ( r_tmpPTemp(i+1) - r_tmpPTemp(i) ) * ratio
          if (TempCol  > 0) r_Temp(k)  = r_tmpTemp(i)  + ( r_tmpTemp(i+1)  - r_tmpTemp(i)  ) * ratio
        end if
      end do
    end do

    do k = kmin, kmax
      do i = 1, NumRec - 1 
        if ( r_tmpAlt(i) == z_Z(k) ) then 
          if (PTempCol > 0) z_PTemp(k) = r_tmpPTemp(i)
          if (TempCol  > 0) z_Temp(k)  = r_tmpTemp(i) 
        elseif ( r_tmpAlt(i) < z_Z(k) .AND. z_Z(k) <= r_tmpAlt(i+1) ) then 
          ratio = ( z_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          if (PTempCol > 0) z_PTemp(k) = r_tmpPTemp(i) + ( r_tmpPTemp(i+1) - r_tmpPTemp(i) ) * ratio
          if (TempCol  > 0) z_Temp(k)  = r_tmpTemp(i)  + ( r_tmpTemp(i+1)  - r_tmpTemp(i)  ) * ratio
        end if
      end do
    end do

!!    do k = kmin+1, kmax
!!      if (PTempCol > 0) z_PTemp(k) = ( r_PTemp(k-1) + r_PTemp(k) ) * 5.0d-1
!!      if (TempCol  > 0) z_Temp(k)  = ( r_Temp(k-1)  + r_Temp(k)  ) * 5.0d-1
!!    end do

    !! データの入っている配列サイズを得る
    !!
    kmax2 = nz
    do k = 1, nz      
      if ( max( r_Temp(k), r_PTemp(k) ) == 0.0d0 ) then 
        kmax2 = k - 1
        exit
      end if
    end do

    !! 静水圧平衡から基本場の圧力を作り直す
    !! 高度・温度 => 気圧
    !!
    if ( TempCol > 0 ) then 
      r_Press(0) = PressSfc
      z_Press(1) = r_Press(0) * exp( - Grav * dz * 0.5d0 / GasRDry / r_Temp(0) )

      do k = 1, kmax2
        z_Press(k+1) = z_Press(k) * exp( - Grav * dz / GasRDry / r_Temp(k) )
      end do

      ! 念の為に温位を計算しておく
      !
      do k = 1, kmax2
        z_PTemp(k) = z_Temp(k) * (PressBasis / z_Press(k)) ** (GasRDry / CpDry) 
      end do
    end if

    
    !! 静水圧平衡から基本場の圧力と温度を作り直す
    !! 高度・温位 => 温度・気圧
    !!
    if ( PTempCol > 0 ) then 
      r_Press(0) = PressSfc
      r_Temp(0)  = ( PressBasis / r_Press(0) ) ** ( - GasRDry / CpDry ) * r_PTemp(0)
      
      do k = 1, kmax2
        Dens = r_Press(k-1) * ( (r_Press(k-1) / PressBasis )**(- GasRDry / CpDry)) / GasRDry / r_PTemp(k-1)
        z_Press(k) = r_Press(k-1) + ( - Dens * Grav ) * ( dz * 0.5d0 )
        z_Temp(k)  = ( PressBasis / z_Press(k) ) ** ( - GasRDry / CpDry ) * z_PTemp(k)

        Dens = z_Press(k) * ( (z_Press(k) / PressBasis )**(- GasRDry / CpDry)) / GasRDry / z_PTemp(k)
        r_Press(k) = z_Press(k) + ( - Dens * Grav ) * ( dz * 0.5d0 )
        r_Temp(k)  = ( PressBasis / r_Press(k) ) ** ( - GasRDry / CpDry ) * r_PTemp(k)       
      end do
    end if

    !! サウンディングファイルで与えた高度より上空の処理
    !! VIRA データを使う場合
    !!
    if ( kmax2 < nz ) then 

      do k = kmax2+1, kmax

        if ( r_viraAlt(1) >= r_Z(k) ) then 
          ratio = ( r_Z(k) - r_Z(kmax2) ) / ( r_viraAlt(1) - r_Z(kmax2) ) 
          r_Temp(k)  = r_Temp(kmax2)  + ( r_viraTemp(1)  - r_Temp(kmax2)  ) * ratio
        else
          do i = 1, NumVIRA-1
            if ( r_viraAlt(i) == r_Z(k) ) then 
              r_Temp(k)  = r_viraTemp(i) 
            elseif ( r_viraAlt(i) < r_Z(k) .AND. r_Z(k) <= r_viraAlt(i+1) ) then 
              ratio = ( r_Z(k) - r_viraAlt(i) ) / ( r_viraAlt(i+1) - r_viraAlt(i) ) 
              r_Temp(k)  = r_viraTemp(i)  + ( r_viraTemp(i+1)  - r_viraTemp(i)  ) * ratio
            end if
          end do
        end if
        
        z_Temp(k)  = ( r_Temp(k-1) + r_Temp(k) ) * 5.0d-1    
        z_Press(k) = z_Press(k-1) * exp( - Grav * dz / GasRDry / r_Temp(k-1) )
      end do
    end if

!!    do k = 1, nz
!!      write(*,*) k, z_Z(k), z_Temp(k)
!!      write(*,*) k, r_Z(k), r_Temp(k)
!!    end do


    !! 安定度の計算
    !!
    do k = 1, nz
      z_Stab(k) = Grav / z_Temp(k)                  &
        &    * (                                    &
        &       + (r_Temp(k) - r_Temp(k-1)) / dz    &
        &       + Grav / CpDry                      &
        &      )   
    end do
    
    !! サウンディングデータのファイル出力
    !!
    call CheckOutput


!!!
!!! 圏界面が存在する場合には, そこでの温度勾配を緩やかにする.
!!!
    if ( AltTr > z_Z(1) .AND. AltTr < z_Z(nz) ) then

      ! 現在の温度分布の勾配を取る
      !
      do k = 1, kmax
        z_DTempDZ(k) = (z_Temp(k) - z_Temp(k-1)) / dz
      end do
      
      ! 対流圏界面より上の扱い. 指定された高度の温度減率を使い続ける. 
      !
      k1 = minloc( z_Z, 1, z_Z > AltTr ) 
      
      ! 温度減率
      !
      DTempDZ = z_DTempDZ(k1-1)
      k2 = int( DelAlt / dz )    !整数型へ変換
      
      do k = k1, kmax
        
        ! 温度源率をゼロに近づける. 
        !
        DTempDZ = min( -1.0d-14, DTempDZ - DTempDZ / k2 * ( k - k1 ) )
        
        !基本場の温度を決める
        z_Temp(k) = z_Temp(k-1) + DTempDZ * dz
        
        !圧力を静水圧平衡から計算
        z_Press(k) =                                      &
          &  z_Press(k-1) * ( ( z_Temp(k-1) / z_Temp(k) ) &
          &    ** (Grav / ( DTempDZ * GasRDry ) ) )
        
      end do

    end if

  contains    

    function CalcCp(Temp)

      real(DP) :: Temp
      real(DP) :: CalcCp
      
      real(DP), parameter  :: CpA =  19.796d0
      real(DP), parameter  :: CpB =  7.344d-2
      real(DP), parameter  :: CpC = -5.602d-5
      real(DP), parameter  :: CpD =  1.715d-8

      CalcCp = ( CpA + CpB * Temp + CpC * Temp**2 + CpD * Temp**3 ) / 44.0d-3

    end function CalcCp

    
    subroutine CheckOutput
      
      use gtool_history,  only : HistoryCreate, HistoryClose, &
        &                        HistoryPut, HistoryAddVariable
      use fileset,        only : filetitle,      &!データの表題
        &                        filesource,     &!データを作成する手順
        &                        FileInstitution  !最終変更者・組織
      use mpi_wrapper,    only : myrank

      real(DP) :: z_Stab2(kmin:kmax)            !安定度
      real(DP) :: Cp

      ! CPU のランクが 0 の場合にのみファイル出力
      !
      if (myrank /= 0) return
      
      ! 初期化
      !
      z_Stab2 = 0.0d0

      !! 主プログラムの計算によって, z_PTemp, z_Temp, z_Press は求まっている. 
      ! サウンディングファイルで温度が与えられた場合には, 比熱の温度依存性を
      ! 考慮して安定度を計算する. 
      !
      if (TempCol > 0) then 
        do k = 1, nz
          Cp = CalcCp( z_Temp(k) )          

          z_Stab2(k) = Grav / z_Temp(k)                 &
            &    * (                                    &
            &       + (r_Temp(k) - r_Temp(k-1)) / dz    &
            &       + Grav / Cp                         &
            &      )   
        end do
      end if

      !-------------------------------------------------------------    
      ! ヒストリー作成
      !-------------------------------------------------------------  
      call HistoryCreate(                              &
        & file = Outputfile,                           &
        & title = filetitle,                           &
        & source = filesource,                         &
        & institution = FileInstitution,               &
        & dims=(/'z'/),                                &
        & dimsizes=(/nz/),                             &
        & longnames=(/'Z-coordinate'/),                &
        & units=(/'m'/),                               &
        & xtypes=(/'float'/),                          &
        & history=rstat, quiet=.true. )
      
      call HistoryAddVariable(                         &
        & varname='Stab', dims=(/'z'/),                &
        & longname='Stability',                        &
        & units='1', xtype='double', history=rstat ) 
      
      call HistoryAddVariable(                         &
        & varname='Stab2', dims=(/'z'/),               &
        & longname='Stability with Cp(z)',             &
        & units='1', xtype='double', history=rstat ) 
      
      call HistoryAddVariable(                         &
        & varname='Temp', dims=(/'z'/),                &
        & longname='Temperature',                      &
        & units='1', xtype='double', history=rstat ) 
      
      call HistoryAddVariable(                         &
        & varname='PTemp', dims=(/'z'/),               &
        & longname='Temperature',                      &
        & units='1', xtype='double', history=rstat ) 
      
      call HistoryPut('z', z_Z(1:nz), rstat )

      call HistoryPut( 'Stab',   z_Stab(1:nz),  rstat )
      call HistoryPut( 'Stab2',  z_Stab2(1:nz), rstat )
      call HistoryPut( 'Temp',   z_Temp(1:nz),  rstat )
      call HistoryPut( 'PTemp' , z_PTemp(1:nz), rstat )

      call HistoryClose( rstat )
      
    end subroutine CheckOutput
    
  end subroutine Initialdata_sounding_basic


  subroutine initialdata_sounding_wind(pyz_VelX, xqz_VelY)

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: imin, imax,       &!配列サイズ (X 方向)
      &                 jmin, jmax,       &!配列サイズ (Y 方向)
      &                 kmin, kmax         !配列サイズ (Z 方向)
    use axesset,  only: r_Z                !高度
      
    implicit none

    real(DP), intent(out) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: pyr_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: xqr_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: ratio
    integer               :: i, k
    
    !初期化
    pyz_VelX = 0.0d0
    pyr_VelX = 0.0d0
    xqz_VelY = 0.0d0
    xqr_VelY = 0.0d0

    ! 線形補完
    !
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_tmpAlt(i) <= r_Z(k) .AND. r_Z(k) < r_tmpAlt(i+1) ) then 
          ratio = ( r_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          pyr_VelX(:,:,k) = r_tmpVelX(i) + ( r_tmpVelX(i+1) - r_tmpVelX(i) ) * ratio
          xqr_VelY(:,:,k) = r_tmpVelY(i) + ( r_tmpVelY(i+1) - r_tmpVelY(i) ) * ratio
        end if
      end do
    end do

    ! r => z の変換
    ! 
    do k = kmin+1, kmax
      pyz_VelX(:,:,k) = ( pyr_VelX(:,:,k-1) + pyr_VelX(:,:,k) ) * 5.0d-1
      xqz_VelY(:,:,k) = ( xqr_VelY(:,:,k-1) + xqr_VelY(:,:,k) ) * 5.0d-1
    end do

  end subroutine initialdata_sounding_wind
  

  subroutine devidecsv( buf, eachcol, maxch, num )
    !
    ! csv 形式のファイルから必要な情報を読み出す
    !

    !暗黙の型宣言禁止
    implicit none

    !変数定義
    integer          :: maxch, num
    character(len=*) :: buf, eachcol(maxch)
    integer          :: i, j, prev, now, now1, ll
    logical          :: quote
    
    quote = .FALSE.

    ll = len_trim(buf)

    prev=0; now=1; i=1
    do j=1, ll
      if( buf(j:j) == '"' ) then
        quote = .NOT. quote
      end if
      if( quote .eqv. .FALSE. ) then
        if(buf(j:j) == ',') then
          now = j
          now1 = now -1
          if( now1 < prev+1 ) then
            eachcol(i) = ''
          else
            eachcol(i) = buf(prev+1:now1)
          end if
          i=i+1
          if( i > maxch ) then
            write(0,*) 'maxch is too small'
            stop
          end if
          prev=now
        end if
      end if
    end do
    
    if( prev < ll ) then
      eachcol(i) =  buf(prev+1:ll)
      num = i
    else
      num=i-1
    end if
    
  end subroutine devidecsv



end module initialdata_sounding
