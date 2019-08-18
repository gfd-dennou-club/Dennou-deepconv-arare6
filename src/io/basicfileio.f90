!= Module BasicFileIO
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: basicfileio.f90,v 1.17 2014/07/08 00:59:55 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module BasicFileIO
  !
  ! 基本場の情報を netCDF ファイルに出力するためのルーチン
  !

  !モジュール読み込み
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING

  !暗黙の型宣言禁止
  implicit none
  
  !関数を public に指定
  public BasicfileIO_Output

  type(GT_HISTORY),  save, private :: rstat
  character(STRING), save, private :: OutputFile  = "BasicZ.nc"
  
contains

  subroutine Basicfileio_Output
    !
    !基本場の情報を書き出し
    !

    !モジュール呼び出し
    !
    use gtool_history,  only : HistoryCreate, HistoryClose, &
      &                        HistoryPut, HistoryAddVariable
    use mpi_wrapper,    only : myrank
    use axesset,        only : z_Z              !Z 座標軸(スカラー格子点)
    use gridset,        only : nz,             &!物理領域の大きさ
      &                        ncmax            !凝縮成分の数
    use fileset,        only : filetitle,      &!データの表題
      &                        filesource,     &!データを作成する手順
      &                        FileInstitution  !最終変更者・組織
    use basicset,       only : xyz_DensBZ,        &
      &                        xyz_ExnerBZ,       &
      &                        xyz_PTempBZ,       &
      &                        xyz_VPTempBZ,      &
      &                        xyz_VelSoundBZ,    &
      &                        xyz_TempBZ,        &
      &                        xyz_PressBZ,       &
      &                        xyzf_QMixBZ,       &
      &                        xyz_EffMolWtBZ
    use composition,    only : SpcWetSymbol
    
    !暗黙の型宣言禁止
    !
    implicit none

    !変数定義
    !
    integer  :: s   

    ! ノードの番号が 0 の場合にはファイル出力する. 
    !
    if (myrank /= 0) return

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
      & units=(/'m  '/),                             &
      & xtypes=(/'double'/),                         &
      & history=rstat, quiet=.true. )
    
    !無次元圧力の基本場
    call HistoryAddVariable(                             &
      & varname='ExnerBZ', dims=(/'z'/),                 &
      & longname='nondimensional pressure', units='1',   &
      & xtype='double', history=rstat )
    
    !温位の基本場
    call HistoryAddVariable(                             &
      & varname='PTempBZ', dims=(/'z'/),                 &
      & longname='potential temperature',                &
      & units='K', xtype='double', history=rstat ) 

    !温位の基本場
    call HistoryAddVariable(                             &
      & varname='VPTempBZ', dims=(/'z'/),                &
      & longname='virtual potential temperature',        &
      & units='K', xtype='double', history=rstat ) 
    
    !密度の基本場
    call HistoryAddVariable(                             &
      & varname='DensBZ', dims=(/'z'/),                  &
      & longname='density',                              &
      & units='Kg.m-3', xtype='double', history=rstat )
    
    !音波速度の基本場
    call HistoryAddVariable(                             &
      & varname='VelSoundBZ', dims=(/'z'/),              &
      & longname='sound velocity',                       &
      & units='m.s-1', xtype='double', history=rstat )

    !温度の基本場
    call HistoryAddVariable(                             & 
      & varname='TempBZ', dims=(/'z'/),                  &
      & longname='Temperature of basic state',           &
      & units='K', xtype='double', history=rstat ) 
    
    !圧力の基本場
    call HistoryAddVariable(                             &
      & varname='PressBZ', dims=(/'z'/),                 &
      & longname='Pressure of basic state',              &
      & units='Pa', xtype='double', history=rstat ) 
    
    !分子量効果
    call HistoryAddVariable(                             &
      & varname='EffMolWtBZ', dims=(/'z'/),              &
      & longname='Effect of Mole Weight',                &
      & units='1', xtype='double', history=rstat ) 

    do s = 1, ncmax

      !混合比の基本場
      call HistoryAddVariable(                           &
        & varname=trim(SpcWetSymbol(s))//'BZ',           &
        & dims=(/'z'/),                                  &
        & longname=trim(SpcWetSymbol(s))//               &
        &   ' Mixing Ratio of basic state',              &
        & units='kg.kg-1', xtype='double', history=rstat ) 
      
    end do

    !-------------------------------------------------------------  
    ! 変数出力
    !-------------------------------------------------------------
    call HistoryPut('z', z_Z(1:nz), rstat )
    
    call HistoryPut( 'DensBZ',     xyz_DensBZ(1,1,1:nz),     rstat )
    call HistoryPut( 'ExnerBZ',    xyz_ExnerBZ(1,1,1:nz),    rstat )
    call HistoryPut( 'PTempBZ',    xyz_PTempBZ(1,1,1:nz),    rstat )
    call HistoryPut( 'VPTempBZ',   xyz_VPTempBZ(1,1,1:nz),   rstat )
    call HistoryPut( 'VelSoundBZ', xyz_VelSoundBZ(1,1,1:nz), rstat )
    call HistoryPut( 'TempBZ',     xyz_TempBZ(1,1,1:nz),     rstat )
    call HistoryPut( 'PressBZ',    xyz_PressBZ(1,1,1:nz),    rstat )
    call HistoryPut( 'EffMolWtBZ', xyz_EffMolWtBZ(1,1,1:nz), rstat )
    
    do s = 1, ncmax
      call HistoryPut( trim(SpcWetSymbol(s))//'BZ', xyzf_QMixBZ(1,1,1:nz,s), rstat )
    end do
    
    !-------------------------------------------------------------
    ! ファイルのクローズ
    !
    call HistoryClose( rstat )
    
  end subroutine Basicfileio_Output
          
end module BasicFileIO
