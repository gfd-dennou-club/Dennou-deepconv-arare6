!= 簡単放射: サウンディングファイルから加熱率を与える
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: radiation_sounding.f90,v 1.2 2014/03/04 04:49:42 sugiyama Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module Radiation_Sounding
  !
  ! 簡単放射: サウンディングファイルから加熱率を与える

  ! モジュール読み込み
  !
  use dc_types, only: DP
  
  ! 暗黙の型宣言禁止
  !
  implicit none

  ! private 属性
  !
  private

  !変数定義
  !
  real(DP), save, allocatable, public :: xyz_DPTempDtRadSndg(:,:,:)   
                                           !放射加熱項
  real(DP), save, allocatable, public :: xyz_ExnerRadSndg(:,:,:)  
                                           !放射加熱項
  real(DP), save :: FactorDExnerDtRad = 1.0d0

  public  radiation_sounding_init
  public  radiation_sounding_forcing

contains

!!!----------------------------------------------------------------------!!!
  subroutine radiation_sounding_init
    !
    ! サウンディングファイルから加熱率を読み込む
    !
    
    ! モジュール読み込み
    !
    use dc_types,          only : DP
    use dc_iounit,         only : FileOpen
!    use dc_message,        only : MessageNotify
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax        !z 方向の配列の上限
    use axesset,           only : z_Z, r_Z    !Z 座標軸
    use constants,         only : DayTime     ! 1 日の長さ [s]
    use basicset,          only : xyz_ExnerBZ !エクスナー関数の基本場
    use DExnerDt,          only : xyz_DExnerDt_xyz
    use namelist_util,     only : namelist_filename
    
    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    integer             :: AltCol = 0       !「高度」の列番号 (サウンディングファイル内)
    integer             :: SWaveCol = 0     !「短波放射による加熱率」の列番号 (サウンディングファイル内)
    integer             :: LWaveCol = 0     !「長波放射による加熱率」の列番号 (サウンディングファイル内)
    integer             :: unit             !設定ファイル用装置番号  
    integer             :: k
    integer             :: io
    character(30)       :: SoundingFile     ! サウンディングファイル
    
    integer, parameter  :: maxch=12
    character(len=200)  :: buf, eachcol(maxch)
    integer             :: i, num, MaxCol

    real(DP)            :: r_tmpQradSW(10000)
    real(DP)            :: r_tmpQradLW(10000)
    real(DP)            :: r_tmpAlt(10000)
    real(DP)            :: r_tmpQrad(10000)
    integer             :: NumRec = 0

    real(DP), allocatable :: xyr_DPTempDtRadSndg(:,:,:)   

    logical             :: flag

   
    !設定ファイルから読み込む出力ファイル情報
    !
    NAMELIST /radiation_sounding_nml/ &
      & SoundingFile, AltCol, SWaveCol, LWaveCol, FactorDExnerDtRad
              
    !設定ファイルから出力ファイルに記載する情報を読み込む
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=radiation_sounding_nml)
    close(unit)
      
    ! 初期化
    !
    io = 0
    NumRec = 0
    MaxCol = max( AltCol, max( SWaveCol, LWaveCol ) )
!    write(*,*) "MaxCol", MaxCol
    r_tmpAlt  = 0.0d0; r_tmpQrad = 0.0d0; r_tmpQradSW = 0.0d0; r_tmpQradLW = 0.0d0

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
      
      ! 確認
      !
!      write(*,*) num
!      do i=1, num
!        write(*,*) i, eachcol(i)(1:len_trim(eachcol(i)))
!      end do
      
      ! num の値が小さいものはヘッダとみなす. 
      !
      if (num >= MaxCol) then 
        ! 行数の計算
        !
        NumRec = NumRec + 1        
        
        ! 値の代入
        !
        if (AltCol > 0)   read( eachcol(AltCol)(1:len_trim(eachcol(AltCol))), *)   r_tmpAlt(NumRec) 
        if (SWaveCol > 0) read( eachcol(SWaveCol)(1:len_trim(eachcol(SWaveCol))), *) r_tmpQradSW(NumRec) 
        if (LWaveCol > 0) read( eachcol(LWaveCol)(1:len_trim(eachcol(LWaveCol))), *) r_tmpQradLW(NumRec) 

        r_tmpQrad(NumRec) = r_tmpQradSW(NumRec) + r_tmpQradLW(NumRec)

      end if

!      write(*,*) r_tmpAlt(NumRec), r_tmpQradSW(NumRec), r_tmpQradLW(NumRec), r_tmpQrad(NumRec) 

    end do

    ! ファイルのクローズ
    !
    close (17)

    !初期化
    !
    allocate( xyz_DPTempDtRadSndg(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyr_DPTempDtRadSndg(imin:imax, jmin:jmax, kmin:kmax) )
    allocate( xyz_ExnerRadSndg(imin:imax, jmin:jmax, kmin:kmax) )

    xyz_DPTempDtRadSndg = 0.0d0
    xyr_DPTempDtRadSndg = 0.0d0
    xyz_ExnerRadSndg = 0.0d0

    flag = .false.
    
    do k = kmin, kmax
      do i = 1, NumRec
        if ( r_Z(k) == r_tmpAlt(i) ) then           
          xyr_DPTempDtRadSndg(:,:,k) = r_tmpQrad(i)         
        end if
      end do
    end do

    do k = kmin, kmax
      do i = 1, NumRec
        if ( z_Z(k) == r_tmpAlt(i) ) then 
          flag = .true. 
!          write(*,*) r_tmpAlt(i), r_tmpQrad(i), '(', r_tmpQradSW(i), '+', r_tmpQradLW(i), ')'
          xyz_DPTempDtRadSndg(:,:,k) = r_tmpQrad(i)
        end if
      end do
    end do

    if (.NOT. flag) then 
      do k = kmin+1, kmax
        xyz_DPTempDtRadSndg(:,:,k) = ( xyr_DPTempDtRadSndg(:,:,k-1) + xyr_DPTempDtRadSndg(:,:,k) ) * 5.0d-1
      end do
    end if


    ! 単位換算
    !
    xyz_DPTempDtRadSndg = xyz_DPTempDtRadSndg / DayTime / xyz_ExnerBZ

    ! 圧力の tendency
    !
    xyz_ExnerRadSndg = xyz_DExnerDt_xyz(xyz_DPTempDtRadSndg) * FactorDExnerDtRad

    ! ヒストリデータ定義
    ! 
    call HistoryAutoAddVariable(  &
      & varname='DPTempDtRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of potential temperature', &
      & units='K.s-1',    &
      & xtype='float')

    call HistoryAutoAddVariable(  &
      & varname='ExnerRad', &
      & dims=(/'x','y','z','t'/),     &
      & longname='Radiation term of Exner function', &
      & units='K.s-1',    &
      & xtype='float')

  end subroutine radiation_sounding_init


  subroutine Radiation_Sounding_forcing(   &
    & xyz_DPTempDt, xyz_DExnerDt           & !(inout)
    & )

    ! モジュール読み込み
    !
    use dc_types,          only : DP
    use gtool_historyauto, only : HistoryAutoPut
    use timeset,           only : TimeN
    use gridset,           only : imin,     & !x 方向の配列の下限
      &                           imax,     & !x 方向の配列の上限
      &                           jmin,     & !y 方向の配列の下限
      &                           jmax,     & !y 方向の配列の上限
      &                           kmin,     & !z 方向の配列の下限
      &                           kmax,     & !z 方向の配列の上限
      &                           nx, ny, nz
    
    ! 暗黙の型宣言禁止
    !
    implicit none

    ! 入出力変数
    !
    real(DP), intent(inout) :: xyz_DPTempDt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(inout) :: xyz_DExnerDt(imin:imax,jmin:jmax,kmin:kmax)

    ! tendency の更新
    !
    xyz_DPTempDt = xyz_DPTempDt + xyz_DPTempDtRadSndg
    xyz_DExnerDt = xyz_DExnerDt + xyz_ExnerRadSndg
    
    ! ファイル出力
    !    
    call HistoryAutoPut(TimeN, 'DPTempDtRad', xyz_DPTempDtRadSndg(1:nx,1:ny,1:nz))
    call HistoryAutoPut(TimeN, 'ExnerRad', xyz_ExnerRadSndg(1:nx,1:ny,1:nz))
    
  end subroutine Radiation_Sounding_forcing


  subroutine devidecsv( buf, eachcol, maxch, num )
    !
    ! サウンディングファイルを読み込むためのルーチン

    implicit none

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
  
end module Radiation_Sounding
