!= 物理・数学定数設定
!
!= Physical and mathematical constants settings
!
! Authors::   Yasuhiro MORIKAWA, Yoshiyuki O. Takahashi
! Version::   $Id: constants0.f90,v 1.3 2012/07/18 07:30:40 odakker Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module constants0
  !
  != 物理・数学定数設定
  !
  != Physical and mathematical constants settings
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! 物理・数学定数の設定および保管を行います. 
  ! デフォルト値は地球大気を想定した値が設定されています. 
  !
  ! Physical and mathematical constants are set and stored. 
  ! By default, values on atmosphere of earth are set. 
  !
  !== Procedures List
  !
  ! Constants0Init :: 物理定数の設定
  ! ------------  :: ------------
  ! Constants0Init :: Settings of physical constants
  !
  !== NAMELIST
  !
  ! N/A
  !


  ! モジュール引用 ; USE statements
  !

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types, only: DP     ! 倍精度実数型. Double precision. 

  ! 宣言文 ; Declaration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public:: Constants0Init

  ! 公開変数
  ! Public variables
  !
  logical, save, public:: constants0_inited = .false.
                              ! 初期設定フラグ. 
                              ! Initialization flag

  real(DP), parameter, public:: PI = 3.1415926535897932_DP
                              ! $ \pi $ .
                              ! 円周率.  Circular constant
  real(DP), parameter, public:: GasRUniv = 8.314_DP
                              ! $ R^{*} $ [J K-1 mol-1].
                              ! 普遍気体定数.  Universal gas constant
  real(DP), parameter, public:: StB = 5.67e-8_DP
                              ! $ \sigma_{SB} $ . 
                              ! ステファンボルツマン定数. 
                              ! Stefan-Boltzmann constant
  real(DP), parameter, public :: FKarm = 0.4d0 
                              ! カルマン定数
                              ! Karmann constant

  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'constants0'
                              ! モジュールの名称. 
                              ! Module name
  character(*), parameter:: version = &
    & '$Name:  $' // &
    & '$Id: constants0.f90,v 1.3 2012/07/18 07:30:40 odakker Exp $'
                              ! モジュールのバージョン
                              ! Module version

contains

  subroutine Constants0Init
    !
    ! constants0 モジュールの初期化を行います. 
    !
    ! "constants0" module is initialized. 
    !

    ! モジュール引用 ; USE statements
    !

    ! メッセージ出力
    ! Message output
    !
    use dc_message, only: MessageNotify

    ! 宣言文 ; Declaration statements
    !
    implicit none

    ! 実行文 ; Executable statement
    !

    if ( constants0_inited ) return


    ! 印字 ; Print
    !
!    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
!    call MessageNotify( 'M', module_name, '  PI               = %f', d = (/ PI               /) )
!    call MessageNotify( 'M', module_name, '  GasRUniv         = %f', d = (/ GasRUniv         /) )
!    call MessageNotify( 'M', module_name, '  StB              = %f', d = (/ StB              /) )
!    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    constants0_inited = .true.

  end subroutine Constants0Init

end module constants0
