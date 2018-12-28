!= ファイル名に関する設定を行うためのモジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: fileset.f90,v 1.10 2014/01/20 08:12:41 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module fileset
  !
  != ファイル名に関する設定を行うためのモジュール
  !
  !引数に与えられた NAMELIST ファイルから, I/O ファイル名を取得し, 
  !保管するための変数参照型モジュール
  !

  !モジュール読み込み
  !
  use dc_types,      only: STRING        

  !暗黙の型宣言禁止
  !
  implicit none

  ! デフォルトの属性
  !
  private

  !公開変数
  !
  character(STRING), save, public :: FileTitle = 'cloud moist convection experiment'
                              ! 実験名.
                              ! Title of experiment
  character(STRING), save, public :: FileSource = 'deepconv/arare5 (http://www.gfd-dennou.org/library/deepconv)'
                              ! データファイル作成プログラム名. 
                              ! Source of data file
  character(STRING), save, public :: FileInstitution = 'GFD Dennou Club (http://www.gfd-dennou.org)'
                              ! データファイル作成者/グループ.
                              ! Institution or person that changes data files for the last time


  character(STRING), save, public :: ExpTitle          !データの表題
  character(STRING), save, public :: ExpSrc            !データを作成する手順
  character(STRING), save, public :: ExpInst           !最終変更者・組織

  ! 初期化ルーチンの公開
  !
  public fileset_init

contains

  subroutine fileset_init
    !
    !設定ファイルから出力ファイルに記載する情報を読み込む
    !
    
    !モジュール読み込み
    use dc_iounit,     only: FileOpen      
    use dc_message,    only: MessageNotify 
    use namelist_util, only: namelist_filename
    
    !暗黙の型宣言禁止
    implicit none
    
    !内部変数
    integer                       :: unit     !設定ファイル用装置番号
     
    !設定ファイルから読み込む出力ファイル情報
    !
    NAMELIST /fileset_nml/ FileTitle, FileSource, FileInstitution
    
    !設定ファイルから出力ファイルに記載する情報を読み込む
    !
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=fileset_nml)
    close(unit)
    
    ! 読み込んだ情報を出力
    !
    call MessageNotify( "M", &
      & "fileset_init", "FileTitle = %c",    c1=trim(FileTitle) )
    call MessageNotify( "M", &
      & "fileset_init", "FileSource = %c",   c1=trim(FileSource) )
    call MessageNotify( "M", &
      & "fileset_init", "FileInstitution = %c", c1=trim(FileInstitution) )

    !----------------------------------------------------
    ! arare4 用の変数
    !
    ExpTitle= FileTitle
    ExpSrc  = FileSource 
    ExpInst = FileInstitution

  end subroutine fileset_init

end module fileset
