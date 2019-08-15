!= コマンドライン引数の解釈を行うモジュール
!
! Authors::   ODAKA Masatsugu, SUGIYAMA Ko-ichiro
! Version::   $Id: argset.f90,v 1.4 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!

module argset
  !
  != コマンドライン引数の解釈を行うモジュール
  !
  != Procedures List
  !
  ! argset_init :: 初期化ルーチン
  !

  !モジュール引用
  !
  use dc_types, only : STRING
  use dc_args,  only : ARGS
  
  !暗黙の型宣言禁止
  !
  implicit none

  !公開変数
  !-----   コマンドライン引数解析用変数    -----
  type(ARGS), save, private        :: arg             ! コマンドライン引数情報
  logical, save, private           :: OPT_namelist    ! コマンドライン引数用論理変数
  character(STRING), save, private :: VAL_namelist    ! コマンドライン引数の値

contains

  subroutine argset_init(namelist_filename)
    !
    !コマンドライン引数を解釈して, 引数に与えられた NAMELIST ファイル名
    !を返す. 
    !

    !モジュール引用
    !
    use dc_types,    only : STRING
    use dc_string,   only : StoA
    use dc_args,     only : ARGS, Open, Debug, Help, Strict, Close, Option
    use dc_message,  only : MessageNotify
    
    !暗黙の型宣言禁止
    implicit none

    !入力変数
    character(STRING), intent(out)     :: namelist_filename

    !NAMELIST ファイルの取得
    !  gtool5 ライブラリの dc_args モジュールを利用.
    !  指定可能なオプションは以下の通り.
    !
    !  -N=VAL, --namelist=VAL
    !    specify Namelist file (default is 'arare.conf'). 
    !
    !  -D=VAL, --debug=VAL
    !    call dc_trace#SetDebug (display a lot of messages for debug).
    !    VAL is unit number (default is standard output)
    !
    !  -h=VAL, -H=VAL, --help=VAL
    !    display this help and exit. VAL is unit number (default is
    !    standard output)
    !
    call Open(arg)
    call Option(arg, StoA('-N', '--namelist'), &
      &         OPT_namelist, VAL_namelist, &
      &         help="specify Namelist file (default is 'arare.conf')." )
                      ! "-N/--namelist" オプションの設定
    call Debug(arg)   ! デバッグオプションの自動設定
    call Help(arg)    ! ヘルプオプションの自動設定
    call Strict(arg)  ! 無効なオプション指定時に警告を表示

    !"-N/-namelist" オプションの解釈
    !  与えられていない場合はデフォルト値 (arare.conf) を 
    !  NAMLIST ファイル名とする.
    !
    if (OPT_namelist) then
      call MessageNotify( "M", "main", &
        &                 "Namelist file is '%c'", c1=trim(VAL_namelist) )
      namelist_filename=trim(VAL_namelist)
    else
      call MessageNotify( "W", "main", &
        &                 "Namelist file is not specified." )
      call MessageNotify( "M", "main", &
        &                 "Use default Namelist file (arare.conf)." )
      namelist_filename="arare.conf"
    end if

    call Close(arg)

    ! 確認
    !
    call MessageNotify( "M", "argset_init",   &
      &                 "NAMELIST FILE = %c", c1=namelist_filename )

  end subroutine argset_init

end module argset
