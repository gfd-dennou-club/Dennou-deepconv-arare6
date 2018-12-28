!= NAMELIST ファイル入力に関するユーティリティ
!
!= Utilities for NAMELIST file input
!
! Authors::   Yoshiyuki O. Takahashi, Yasuhiro MORIKAWA, SUGIYAMA Ko-ichiro
! Version::   $Id: namelist_util.f90,v 1.2 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module namelist_util
  !
  != NAMELIST ファイル入力に関するユーティリティ
  !
  != Utilities for NAMELIST file input
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  !== Variables List
  !
  ! namelist_filename :: NAMELIST ファイルの名称. 
  ! MaxNmlArySize     :: NAMELIST から読み込む配列の最大サイズ. 
  ! ------------      :: ------------
  ! namelist_filename :: NAMELIST file name
  ! MaxNmlArySize     :: Maximum size of arrays loaded from NAMELIST
  !
  !== Procedures List
  !
  ! NmlutilInit      :: NAMELIST ファイル名の設定
  ! NmlutilMsg       :: NAMELIST ファイル入力に関するメッセージ表示
  ! NmlutilAryValid  :: NAMELIST ファイルから読み込んだ配列の有効性をチェック
  ! ------------  :: ------------
  ! NmlutilInit      :: Settings of NAMELIST file name
  ! NmlutilMsg       :: Print messages about NAMELIST file input
  ! NmlutilAryValid  :: Check validation of arrays loaded from NAMELIST file

  ! モジュール引用 ; USE statements
  !

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types,    only: STRING   ! 文字列. Strings. 

  ! メッセージ出力
  ! Message output
  !
  use dc_message, only: MessageNotify

  ! 宣言文 ; Declaration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public:: NmlutilInit, NmlutilMsg, NmlutilAryValid

  ! 公開変数
  ! Public variables
  !
  logical, save, public:: namelist_util_inited = .false.
                              ! 初期設定フラグ. 
                              ! Initialization flag
  character(STRING), save, public:: namelist_filename = ''
                              ! NAMELIST ファイルの名称. 
                              ! NAMELIST file name
  integer, parameter, public:: MaxNmlArySize = 256
                              ! NAMELIST から読み込む配列の最大サイズ. 
                              ! Maximum size of arrays loaded from NAMELIST

  ! 非公開変数
  ! Private variables
  !

  !  NAMELIST 変数群
  !  NAMELIST group name
  !
!!$  namelist /namelist_util_nml/ MaxNmlArySize

  character(*), parameter:: module_name = 'namelist_util'
                              ! モジュールの名称. 
                              ! Module name
  character(*), parameter:: version = &
    & '$Name:  $' // &
    & '$Id: namelist_util.f90,v 1.2 2014/01/20 08:12:41 sugiyama Exp $'
                              ! モジュールのバージョン
                              ! Module version

contains

  subroutine NmlutilInit(  &
    & namelist_filename_in & ! (in)
    & )
    !
    ! namelist_util モジュールの初期設定を行います. 
    !
    ! Initialize "namelist_util" module. 
    !

    ! モジュール引用 ; USE statements
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 宣言文 ; Declaration statements
    !
    implicit none

    character(*), intent(in ) :: namelist_filename_in
                              ! NAMELIST ファイルの名称. 
                              ! NAMELIST file name

!!$    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
!!$                              ! Unit number for NAMELIST file open
!!$    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
!!$                              ! IOSTAT of NAMELIST read

    ! 実行文 ; Executable statement
    !

    if ( namelist_util_inited ) return
    namelist_filename = namelist_filename_in

    ! NAMELIST の読み込み
    ! NAMELIST is input
    !
!!$    call FileOpen( unit_nml, &          ! (out)
!!$      & namelist_filename, mode = 'r' ) ! (in)
!!$
!!$    rewind( unit_nml )
!!$    read( unit_nml, &               ! (in)
!!$      & nml = namelist_util_nml, &  ! (out)
!!$      & iostat = iostat_nml )       ! (out)
!!$    close( unit_nml )
!!$
!!$    call NmlutilMsg( iostat_nml, module_name ) ! (in)

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  MaxNmlArySize = %d', i = (/ MaxNmlArySize /) )
    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    ! 初期化フラグ
    !
    namelist_util_inited = .true.

  end subroutine NmlutilInit


  subroutine NmlutilMsg( iostat, name & ! (in)
    & )
    !
    ! NAMELIST ファイル入力ステータスから
    ! 適切なメッセージを表示します. 
    !
    ! Appropriate messages are output from 
    ! status of NAMELIST file loading. 
    !

    implicit none
    integer, intent(in):: iostat
                              ! NAMELIST 読み込み時のステータス. 
                              ! Status of NAMELIST loading
    character(*), intent(in):: name
                              ! モジュールの名称. 
                              ! Module name

    ! 実行文 ; Executable statement
    !

    if ( iostat == 0 ) then
      call MessageNotify( 'M', name, &
        & 'NAMELIST group "%c" is loaded from "%c".', &
        & c1 = trim(name) // '_nml', &
        & c2 = trim(namelist_filename) )
    else
      call MessageNotify( 'W', name, &
        & 'NAMELIST group "%c" is not found in "%c" (iostat=%d).', &
        & c1 = trim(name) // '_nml', &
        & c2 = trim(namelist_filename), &
        & i = (/iostat/) )
    end if
  end subroutine NmlutilMsg

  subroutine NmlutilAryValid( name, &               ! (in)
    & array, array_name, need_num, need_num_name, & ! (in)
    & valid_limit &                                 ! (in) optional
    )
    !
    ! NAMELIST から読み込んだ配列型データの妥当性を
    ! チェックします. 
    ! 
    ! デフォルトでは, 正の値を有効と扱います. 
    ! 無効であると検証された場合には, エラーを発生させます. 
    ! 
    ! Check validation of array data loaded from NAMELIST. 
    !
    ! By defaut, positive values are treated as valid values.
    ! If invalidation is checked, an error is occurred. 
    ! 

    ! モジュール引用 ; USE statements
    !

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: DP      ! 倍精度実数型. Double precision. 

    ! 宣言文 ; Declaration statements
    !
    implicit none
    character(*), intent(in):: name
                              ! このサブルーチンを呼び出すモジュールの名称. 
                              ! Module name calling this subroutine
    real(DP), intent(in):: array(:)
                              ! 検証すべき配列データ. 
                              ! Checked array data 
    character(*), intent(in):: array_name
                              ! 検証すべき配列データの名前. 
                              ! Name of checked array data 
    integer, intent(in):: need_num
                              ! 必要なデータ数. 
                              ! 0 未満の数を与えるとエラーを生じます. 
                              ! 
                              ! Number of needed data. 
                              ! If number less than 0, an error is occurred. 
    character(*), intent(in):: need_num_name
                              ! 必要なデータ数を示す変数の名前. 
                              ! Name of a variable that indicates number of needed data
    real(DP), intent(in), optional:: valid_limit
                              ! 有効下限値 (デフォルトは 0.0). 
                              ! Lower limit of validation (defalt is 0.0) 

    ! 作業変数
    ! Work variables
    !
    real(DP):: valid_limit_work
                              ! 有効下限値 (デフォルトは 0.0). 
                              ! Lower limit of validation (defalt is 0.0) 
    integer:: valid_count     ! 配列データの有効値の数. 
                              ! Number of valid values in an array
    integer:: size_array      ! 配列データのサイズ
                              ! Size of array data

    ! 実行文 ; Executable statement
    !

    ! need_num が負ではないことをチェック
    ! Check that "need_num" is not negative 
    !
    if ( need_num < 0 ) then
      call MessageNotify( 'E', name, '%c=<%d> must not be negative.', &
        &                 c1 = trim(need_num_name), i = (/ need_num /) )
    end if

    ! array のサイズが十分であることをチェック
    ! Check that size of "array" is enough
    !
    size_array = size(array)
    if ( need_num > size_array ) then
      call MessageNotify( 'E', name, &
        &  'Maximum size=<%d> of "%c" is too smaller than %c=<%d>. ' // &
        &  'Please search for a statement "MaxNmlArySize = %d" in ' // &
        &  '"namelist_util.f90", and change it into "MaxNmlArySize = %d".', &
!!$        &  'Please add a following phrase to NAMELIST file. ' // &
!!$        &  ' "&namelist_util_nml  MaxNmlArySize=%d /"', &
        &  i = (/ size_array, need_num, MaxNmlArySize, need_num /), &
        &  c1 = trim(array_name), c2 = trim(need_num_name) )
    end if

    ! array をチェック
    ! Check "array"
    !
    if ( need_num > 0 ) then
      valid_limit_work = 0.0_DP
      if ( present( valid_limit ) ) valid_limit_work = valid_limit

      if ( any( array(1:need_num) < valid_limit_work ) ) then
        valid_count = count( .not. ( array(1:need_num) < valid_limit_work ) )
        if ( valid_count > 0 ) then
          call MessageNotify( 'E', name, &
            &   'Number of valid data of %c=<%*f> is %d. ' // &
            &   'Valid data is %c=<%d> necessary.', &
            &   c1 = trim( array_name ), c2 = trim( need_num_name ), &
            &   d = array(1:valid_count), n = (/ valid_count /), &
            &   i = (/ valid_count, need_num /) )
        else
          call MessageNotify( 'E', name, &
            &   'Valid data of %c is nothing. ' // &
            &   'Valid data is %c=<%d> necessary.', &
            &   c1 = trim( array_name ), c2 = trim( need_num_name ), &
            &   i = (/ need_num /) )
        end if
      end if
    end if

  end subroutine NmlutilAryValid

end module namelist_util
