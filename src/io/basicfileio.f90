!= Module BasicFileIO
!
! Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
! Version::   $Id: basicfileio.f90,v 1.17 2014/07/08 00:59:55 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2007. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module BasicFileIO
  !
  ! $B4pK\>l$N>pJs$r(B netCDF $B%U%!%$%k$K=PNO$9$k$?$a$N%k!<%A%s(B
  !

  !$B%b%8%e!<%kFI$_9~$_(B
  use gtool_history, only : GT_HISTORY
  use dc_types,      only : STRING

  !$B0EL[$N7?@k8@6X;_(B
  implicit none
  
  !$B4X?t$r(B public $B$K;XDj(B
  public BasicfileIO_Output

  type(GT_HISTORY),  save, private :: rstat
  character(STRING), save, private :: OutputFile  = "BasicZ.nc"
  
contains

  subroutine Basicfileio_Output
    !
    !$B4pK\>l$N>pJs$r=q$-=P$7(B
    !

    !$B%b%8%e!<%k8F$S=P$7(B
    !
    use gtool_history,  only : HistoryCreate, HistoryClose, &
      &                        HistoryPut, HistoryAddVariable
    use mpi_wrapper,    only : FLAG_LIB_MPI
    use axesset,        only : x_X,            &!X $B:BI8<4(B($B%9%+%i!<3J;RE@(B)
      &                        y_Y,            &!Y $B:BI8<4(B($B%9%+%i!<3J;RE@(B)
      &                        z_Z              !Z $B:BI8<4(B($B%9%+%i!<3J;RE@(B)
    use gridset,        only : nx, ny, nz,     &!$BJ*M}NN0h$NBg$-$5(B
      &                        ncmax            !$B6E=L@.J,$N?t(B
    use fileset,        only : filetitle,      &!$B%G!<%?$NI=Bj(B
      &                        filesource,     &!$B%G!<%?$r:n@.$9$k<j=g(B
      &                        FileInstitution  !$B:G=*JQ99<T!&AH?%(B
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
    
    !$B0EL[$N7?@k8@6X;_(B
    !
    implicit none

    !$BJQ?tDj5A(B
    !
    integer  :: s   

    !-------------------------------------------------------------    
    ! $B%R%9%H%j!<:n@.(B
    !-------------------------------------------------------------  
    call HistoryCreate(                              &
      & file = Outputfile,                           &
      & title = filetitle,                           &
      & source = filesource,                         &
      & institution = FileInstitution,               &
      & dims=(/'x','y','z'/),                        &
      & dimsizes=(/nx, ny, nz/),                     &
      & longnames=(/'X-coordinate',                  &
      &             'Y-coordinate',                  &
      &             'Z-coordinate'/),                &
      & units=(/'m  ','m  ','m  '/),                 &
      & xtypes=(/'double', 'double', 'double'/),     &
      & flag_mpi_split = FLAG_LIB_MPI,               &
      & history=rstat, quiet=.true. )
    
    !$BL5<!8505NO$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='ExnerBZ', dims=(/'x','y','z'/),         &
      & longname='nondimensional pressure', units='1',   &
      & xtype='double', history=rstat )
    
    !$B290L$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='PTempBZ', dims=(/'x','y','z'/),         &
      & longname='potential temperature',                &
      & units='K', xtype='double', history=rstat ) 

    !$B290L$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='VPTempBZ', dims=(/'x','y','z'/),        &
      & longname='virtual potential temperature',        &
      & units='K', xtype='double', history=rstat ) 
    
    !$BL)EY$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='DensBZ', dims=(/'x','y','z'/),          &
      & longname='density',                              &
      & units='Kg.m-3', xtype='double', history=rstat )
    
    !$B2;GHB.EY$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='VelSoundBZ', dims=(/'x','y','z'/),      &
      & longname='sound velocity',                       &
      & units='m.s-1', xtype='double', history=rstat )

    !$B29EY$N4pK\>l(B
    call HistoryAddVariable(                             & 
      & varname='TempBZ', dims=(/'x','y','z'/),          &
      & longname='Temperature of basic state',           &
      & units='K', xtype='double', history=rstat ) 
    
    !$B05NO$N4pK\>l(B
    call HistoryAddVariable(                             &
      & varname='PressBZ', dims=(/'x','y','z'/),         &
      & longname='Pressure of basic state',              &
      & units='Pa', xtype='double', history=rstat ) 
    
    !$BJ,;RNL8z2L(B
    call HistoryAddVariable(                             &
      & varname='EffMolWtBZ', dims=(/'x','y','z'/),      &
      & longname='Effect of Mole Weight',                &
      & units='1', xtype='double', history=rstat ) 

    do s = 1, ncmax

      !$B:.9gHf$N4pK\>l(B
      call HistoryAddVariable(                           &
        & varname=trim(SpcWetSymbol(s))//'BZ',           &
        & dims=(/'x','y','z'/),                          &
        & longname=trim(SpcWetSymbol(s))//               &
        &   ' Mixing Ratio of basic state',              &
        & units='kg.kg-1', xtype='double', history=rstat ) 
      
    end do

    !-------------------------------------------------------------  
    ! $BJQ?t=PNO(B
    !-------------------------------------------------------------
    call HistoryPut('x', x_X(1:nx), rstat )
    call HistoryPut('y', y_Y(1:ny), rstat )
    call HistoryPut('z', z_Z(1:nz), rstat )
    
    call HistoryPut( 'DensBZ',     xyz_DensBZ(1:nx,1:ny,1:nz),     rstat )
    call HistoryPut( 'ExnerBZ',    xyz_ExnerBZ(1:nx,1:ny,1:nz),    rstat )
    call HistoryPut( 'PTempBZ',    xyz_PTempBZ(1:nx,1:ny,1:nz),    rstat )
    call HistoryPut( 'VPTempBZ',   xyz_VPTempBZ(1:nx,1:ny,1:nz),   rstat )
    call HistoryPut( 'VelSoundBZ', xyz_VelSoundBZ(1:nx,1:ny,1:nz), rstat )
    call HistoryPut( 'TempBZ',     xyz_TempBZ(1:nx,1:ny,1:nz),     rstat )
    call HistoryPut( 'PressBZ',    xyz_PressBZ(1:nx,1:ny,1:nz),    rstat )
    call HistoryPut( 'EffMolWtBZ', xyz_EffMolWtBZ(1:nx,1:ny,1:nz), rstat )
    
    do s = 1, ncmax
      call HistoryPut( trim(SpcWetSymbol(s))//'BZ', xyzf_QMixBZ(1:nx,1:ny,1:nz,s), rstat )
    end do
    
    !-------------------------------------------------------------
    ! $B%U%!%$%k$N%/%m!<%:(B
    !
    call HistoryClose( rstat )
    
  end subroutine Basicfileio_Output
          
end module BasicFileIO
