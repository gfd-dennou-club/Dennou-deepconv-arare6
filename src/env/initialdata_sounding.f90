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
  use dc_types, only: DP
  
  !暗黙の型宣言禁止
  implicit none

  !共通変数
  real(DP), save, private :: r_tmpAlt(10000)
  real(DP), save, private :: r_tmpTemp(10000)
  real(DP), save, private :: r_tmpPress(10000)
  real(DP), save, private :: r_tmpPTemp(10000)
  real(DP), save, private :: r_tmpVelX(10000)
  real(DP), save, private :: r_tmpVelY(10000)
  real(DP), save, private :: r_tmpQrad(10000)
  real(DP), save, private :: TempTr = 0.0d0
  real(DP), save, private :: AltTr  = 0.0d0
  real(DP), save, private :: DelAlt = 4.0d3
  integer,  save, private :: NumRec = 0

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
    integer             :: AltCol = 0
    integer             :: TempCol = 0
    integer             :: PressCol = 0
    integer             :: VelXCol = 0
    integer             :: VelYCol = 0
    integer             :: unit         !設定ファイル用装置番号  
    integer             :: io
    character(30)       :: SoundingFile    
    
    integer, parameter  :: maxch=12
    character(len=100)  :: buf, eachcol(maxch)
    integer             :: num, MaxCol
    
    !設定ファイルから読み込む出力ファイル情報
    !
    NAMELIST /initialdata_sounding_nml/ SoundingFile, AltCol, TempCol, PressCol, VelXCol, VelYCol, TempTr, AltTr, DelAlt    

    !設定ファイルから読み込み
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=initialdata_sounding_nml)
    close(unit)

    ! 初期化
    !
    io = 0
    NumRec = 0
    MaxCol = max( AltCol, max( TempCol, PressCol ) )
    r_tmpAlt  = 0.0d0; r_tmpTemp = 0.0d0; r_tmpPress = 0.0d0; 
    r_tmpVelX = 0.0d0; r_tmpVelY = 0.0d0

    ! ファイルのオープン
    !
    open (17, file=SoundingFile, status='old')
    
    ! ファイル呼び出し
    !
    do while ( io == 0 ) 
      ! 1 行分読み出し
      !
      read (17, '(a)', IOSTAT=io) buf
      
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
        if (AltCol > 0)   read( eachcol(AltCol)(1:len_trim(eachcol(AltCol))), *)    r_tmpAlt(NumRec) 
        if (TempCol > 0)  read( eachcol(TempCol)(1:len_trim(eachcol(TempCol))), *)   r_tmpTemp(NumRec)   
        if (PressCol > 0) read( eachcol(PressCol)(1:len_trim(eachcol(PressCol))), *) r_tmpPress(NumRec) 
        if (VelXCol > 0)  read( eachcol(VelXCol)(1:len_trim(eachcol(VelXCol))), *)   r_tmpVelX(NumRec)   
        if (VelYCol > 0)  read( eachcol(VelYCol)(1:len_trim(eachcol(VelYCol))), *)   r_tmpVelY(NumRec) 

      end if

    end do

    ! ファイルのクローズ
    !
    close (17)    

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
      &                  GasRDry
    
    !暗黙の型宣言禁止
    implicit none

    !変数定義
    real(DP), intent(out):: z_Press(kmin:kmax)           !圧力
    real(DP), intent(out):: z_Temp(kmin:kmax)            !温度
    real(DP)             :: r_Press(kmin:kmax)           !圧力
    real(DP)             :: r_Temp(kmin:kmax)            !温度
    real(DP)             :: z_DTempDZ(kmin:kmax)
    real(DP)             :: DTempDZ
    real(DP)             :: ratio
    integer              :: i, k, k1, k2
    logical              :: flag

    ! 初期化
    !
    z_Temp  = 0.0d0
    r_Temp  = 0.0d0
    z_Press = 0.0d0
    r_Press = 0.0d0

    flag = .false. 

    ! 線形補完
    !
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_tmpAlt(i) <= r_Z(k) .AND. r_Z(k) < r_tmpAlt(i+1) ) then 
          ratio = ( r_Z(k) - r_tmpAlt(i) ) / ( r_tmpAlt(i+1) - r_tmpAlt(i) ) 
          r_Temp(k)  = r_tmpTemp(i)  + ( r_tmpTemp(i+1) - r_tmpTemp(i) )   * ratio
          r_Press(k) = r_tmpPress(i) + ( r_tmpPress(i+1) - r_tmpPress(i) ) * ratio
        end if
      end do
    end do

    ! データ読み込み  
    ! うまく半格子に合う場合もあるので, その場合はデータを優先. 
    !    
    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
          z_Temp(k) = r_tmpTemp(i)
          z_Press(k) = r_tmpPress(i)
        end if
      end do
    end do
    
    ! r => z の変換
    ! 
    if (.NOT. flag) then 
      do k = kmin+1, kmax
        z_Temp(k)  = ( r_Temp(k-1)  + r_Temp(k)  ) / 2.0d0
        z_Press(k) = ( r_Press(k-1) + r_Press(k) ) / 2.0d0
      end do
    end if
  
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

  end subroutine Initialdata_sounding_basic


  subroutine initialdata_sounding_wind(pyz_VelX, xqz_VelY)

    !モジュール読み込み
    use dc_types, only: DP
    use gridset,  only: imin, imax,       &!配列サイズ (X 方向)
      &                 jmin, jmax,       &!配列サイズ (Y 方向)
      &                 kmin, kmax         !配列サイズ (Z 方向)
    use axesset,  only: r_Z, z_Z           !高度
      
    implicit none

    real(DP), intent(out) :: pyz_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(out) :: xqz_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: pyr_VelX(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: xqr_VelY(imin:imax,jmin:jmax,kmin:kmax)
    real(DP)              :: ratio
    integer               :: i, k
    logical               :: flag
    
    !初期化
    pyz_VelX = 0.0d0
    pyr_VelX = 0.0d0
    xqz_VelY = 0.0d0
    xqr_VelY = 0.0d0

    flag = .false. 

    ! データ読み込み  
    !        
!    do k = kmin, kmax
!      do i = 1, NumRec
!        if ( r_Z(k) == r_tmpAlt(i) ) then 
!          pyr_VelX(:,:,k) = r_tmpVelX(i)
!          xqr_VelY(:,:,k) = r_tmpVelY(i)
!        end if
!      end do
!    end do

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


    ! データ読み込み  
    ! うまく半格子に合う場合もあるので, その場合はデータを優先. 
    !        
    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
          pyz_VelX(:,:,k) = r_tmpVelX(i)
          xqz_VelY(:,:,k) = r_tmpVelY(i)
        end if
      end do
    end do

    ! r => z の変換
    ! 
    if (.NOT. flag) then 
      do k = kmin+1, kmax
        pyz_VelX(:,:,k) = ( pyr_VelX(:,:,k-1) + pyr_VelX(:,:,k) ) * 5.0d-1
        xqz_VelY(:,:,k) = ( xqr_VelY(:,:,k-1) + xqr_VelY(:,:,k) ) * 5.0d-1
      end do
    end if

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
