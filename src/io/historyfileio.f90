!= Module HistoryFileIO
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: historyfileio.f90,v 1.9 2014/03/04 05:55:04 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]


module HistoryFileIO
  !
  !$B%U%!%$%k=PNO(B. $BD9$$;~4V%9%F%C%W$NCM$r=PNO(B.
  !

  ! $B<oJL7?%Q%i%a%?(B
  ! Kind type parameter
  !
  use dc_types, only: STRING     ! $BJ8;zNs(B. Strings. 

  !$B0EL[$N7?@k8@6X;_(B
  implicit none
  
  !$BB0@-$N;XDj(B
  private

  !$B8x3+<jB3$-(B
  public HistoryFileIO_init
  public HistoryFileIO_finalize

  ! $BHs8x3+JQ?t(B
  ! Private variables
  !
  character(STRING), parameter, private :: module_name = 'historyfileio'
                              ! $B%b%8%e!<%k$NL>>N(B. 
                              ! Module name
  character(STRING), parameter, private :: version = &
    & '$Name:  $' // &
    & '$Id: historyfileio.f90,v 1.9 2014/03/04 05:55:04 sugiyama Exp $'
                              ! $B%b%8%e!<%k$N%P!<%8%g%s(B
                              ! Module version

contains 

!!!------------------------------------------------------------------------
  subroutine HistoryFileIO_init
    !
    ! history_file_io $B%b%8%e!<%k$N=i4|2=$r9T$$$^$9(B. 
    !--
    ! NAMELIST#history_file_io_nml $B$NFI$_9~$_$O$3$N<jB3$-$G9T$o$l$^$9(B. 
    !++
    !
    ! "history_file_io" module is initialized. 
    !--
    ! "NAMELIST#history_file_io_nml" is loaded in this procedure. 
    !++
    !

    ! $B%b%8%e!<%k0zMQ(B ; USE statements
    !

    ! gtool5 netCDF $B%G!<%?$NF~=PNO%$%s%?!<%U%'!<%9(B ($BBg5,LO%b%G%kMQ(B)
    ! Interface of Input/Output of gtool5 netCDF data (For large models)
    !
    use gtool_historyauto, only: HistoryAutoCreate,  &
      &                          HistoryAutoAddAttr, &
      &                          HistoryAutoPutAxis, &
      &                          HistoryAutoAddVariable

    ! $B%U%!%$%kF~=PNOJd=u(B
    ! File I/O support
    !
    use dc_iounit, only: FileOpen
    
    ! $B<oJL7?%Q%i%a%?(B
    ! Kind type parameter
    !
    use dc_types,      only: DP, &              ! $BG\@:EY<B?t7?(B. Double precision. 
      &                      STRING             ! $BJ8;zNs(B.       Strings. 
    use mpi_wrapper,   only: FLAG_LIB_MPI
    use namelist_util, only: namelist_filename    
    use axesset,       only: x_X, y_Y, z_Z
    use gridset,       only: nx, ny, nz, ncmax
    use fileset,       only: filetitle,        &!$B%G!<%?$NI=Bj(B
      &                      filesource,       &!$B%G!<%?$r:n@.$9$k<j=g(B
      &                      FileInstitution    !$B:G=*JQ99<T!&AH?%(B
    use timeset,       only: Restarttime, EndTime
    use composition,   only: SpcWetSymbol, GasNum
    
    ! $B@k8@J8(B ; Declaration statements
    !
    implicit none
    
    !$BJQ?tDj5A(B
    real(DP), parameter       :: TimeDisp = 1.0e5 !$B=PNO4V3V$N%G%U%)%k%HCM(B
    integer :: l, s

    !-----------------------------------------------------------
    ! $B%R%9%H%j!<:n@.(B
    !-----------------------------------------------------------
    call HistoryAutoCreate(                             &
      & title = FileTitle,                              &
      & source = FileSource,                            &
      & institution = FileInstitution,                  &
      & dims=(/'x','y','z','t'/),                       &
      & dimsizes=(/nx, ny, nz, 0/),                     &
      & longnames=(/'x-coordinate',                     &
      &             'y-coordinate',                     &
      &             'z-coordinate',                     &
      &             'time        '/),                   &
      & units=(/'m  ','m  ','m  ','sec'/),              &
      & xtypes=(/'double','double','double','double'/), &
      & origin   = Restarttime,                         & 
      & terminus = EndTime,                             &
      & interval = TimeDisp,                            &       
      & flag_mpi_split = FLAG_LIB_MPI,                  &
      & namelist_filename = namelist_filename)  

    call HistoryAutoAddAttr( &
      & varname = 'x', attrname = 'standard_name', &   ! (in)
      & value = 'x-coordinate' )                       ! (in)
    call HistoryAutoAddAttr( &
      & varname = 'y', attrname = 'standard_name', &   ! (in)
      & value = 'y-coordinate' )                       ! (in)
    call HistoryAutoAddAttr( &
      & varname = 'z', attrname = 'standard_name', &   ! (in)
      & value = 'z-coordinate' )                       ! (in)
    
    call HistoryAutoPutAxis('x', x_X(1:nx))
    call HistoryAutoPutAxis('y', y_Y(1:ny))
    call HistoryAutoPutAxis('z', z_Z(1:nz))

    ! $B0u;z(B ; Print
    !
!    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
!    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    !$BL5<!8505NO$N>qMp(B
    call HistoryAutoAddVariable(                           &
      & varname='Exner',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='disturbunce of nondimensional pressure', &
      & units=' ',                                         &
      & xtype='float' )

    !$BL5<!8505NO(B
    call HistoryAutoAddVariable(                           &
      & varname='ExnerAll',                                &
      & dims=(/'x','y','z','t'/),                          &
      & longname='nondimensional pressure',                &
      & units=' ',                                         &
      & xtype='float' )
    
    call HistoryAutoAddVariable(                           &
      & varname='PTemp',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='disturbunce of potential temperature',   &
      & units='K',                                         &
      & xtype='float' )

    !$B290L$N>qMp(B
    call HistoryAutoAddVariable(                           &
      & varname='PTempAll',                                &
      & dims=(/'x','y','z','t'/),                          &
      & longname='potential temperature',                  &
      & units='K',                                         &
      & xtype='float' )

    !$B?eJ?B.EY(B
    call HistoryAutoAddVariable(                           &
      & varname='VelX',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='zonal velocity',                         &
      & units='m.s-1',                                     &
      & xtype='float' )

    !$B?eJ?B.EY(B
    call HistoryAutoAddVariable(                           &
      & varname='VelY',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='meridional velocity',                    &
      & units='m.s-1',                                     &
      & xtype='float' )

    !$B1tD>B.EY(B
    call HistoryAutoAddVariable(                           &
      & varname='VelZ',                                    &
      & dims=(/'x','y','z','t'/),                          &
      & longname='vertical velocity',                      &
      & units='m.s-1',                                     &
      & xtype='float' )

    !$B12G4@-78?t(B($B1?F0NL(B)
    call HistoryAutoAddVariable(                           &
      & varname='Km',                                      &
      & dims=(/'x','y','z','t'/),                          &
      & longname='turbulet diffusion coefficient',         &
      & units='m2.s-1',                                    &
      & xtype='float' )
  
    !$B12G4@-78?t(B($BG.(B)
    call HistoryAutoAddVariable(                           &
      & varname='Kh',                                      &
      & dims=(/'x','y','z','t'/),                          &
      & longname='turbulet diffusion coefficient for heat',&
      & units='m2.s-1',                                    &
      & xtype='float')

    !$B1@L)EY(B
    call HistoryAutoAddVariable(                           &
      & varname='CDens',                                   &
      & dims=(/'x','y','z','t'/),                          &
      & longname='Cloud density',                          &
      & units='kg.m-3',                                    &
      & xtype='float')
    
    ! $B:.9gHf(B
    do l = 1, ncmax
      call HistoryAutoAddVariable(                         &
        & varname=trim(SpcWetSymbol(l)),                   &
        & dims=(/'x','y','z','t'/),                        &
        & longname=trim(SpcWetSymbol(l))//' Mixing Ratio', &
        & units='kg.kg-1',                                 &
        & xtype='float')
    end do

    do l = 1, GasNum      
      call HistoryAutoAddVariable(                         &
        & varname=trim(SpcWetSymbol(l))//'All',            &
        & dims=(/'x','y','z','t'/),                        &
        & longname=trim(SpcWetSymbol(l))//' Mixing Ratio', &
        & units='kg.kg-1',                                 &
        & xtype='float')
    end do


    !----------------------------------------------------------------
    ! Mixing Ratio time change
    !----------------------------------------------------------------
    do s = 1, ncmax
    
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(s))//'DtFill1',             &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Filling Negative term 1 of '               &
        &           //trim(SpcWetSymbol(s))//' mixing ratio',  &
        & units='kg.kg-1.s-1',                                 &
        & xtype='float' )
    
      call HistoryAutoAddVariable(                             &
        & varname='D'//trim(SpcWetSymbol(s))//'DtFill2',             &
        & dims=(/'x','y','z','t'/),                            &
        & longname='Filling Negative term 2 of '               &
        &           //trim(SpcWetSymbol(s))//' mixing ratio',  &
        & units='kg.kg-1.s-1',                                 &
        & xtype='float' )
    
    end do
   
  end subroutine HistoryFileIO_init


  subroutine HistoryFileIO_finalize
    !
    ! $B%R%9%H%j%G!<%?%U%!%$%k=PNO$N=*N;=hM}$r9T$$$^$9(B. 
    !
    ! Terminate history data files output. 

    ! $B%b%8%e!<%k0zMQ(B ; USE statements
    !
    use gtool_historyauto, only: HistoryAutoClose

    ! $B@k8@J8(B ; Declaration statements
    !
    implicit none

    ! $B:n6HJQ?t(B
    ! Work variables
    !

    ! $B<B9TJ8(B ; Executable statement
    !

    call HistoryAutoClose

  end subroutine HistoryFileIO_finalize
  
end module HistoryFileIO


