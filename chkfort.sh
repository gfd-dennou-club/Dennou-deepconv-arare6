#!/bin/sh
#
#= Script of checking Fortran Compilers.
#
# Authors::   Eizi TOYODA, Yasuhiro Morikawa
# Version::   $Id: chkfort.sh,v 1.1 2011/04/26 04:08:10 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2000-2006. All rights reserved.
# License::   See COPYRIGHT[link:COPYRIGHT]
#
#= Overview
#
# This script checks your fortran compilers, and return F90MODTYPE,
# SYSFFLAGS, SYSLDFLAGS, SYSLDLIBS and so on.
#

report_me() {
	echo Error: $*
	echo please report it to GFD-Dennou Club.
	exit 1
}

help() {
	cat <<'END_OF_HELP'
chkfort.sh - Fortran (at least 90) compiler auto-detector

 --- input ---
    FC         - fortran compiler command
    F90MODTYPE - see below
 --- input (sometimes ignored) ---
    SYSFFLAGS  - fortran compiler flags required always
    SYSLDFLAGS - fortran library flags required always
    SYSLDLIBS  - fortran libraries required always
 --- output ---
    output is written on 'chkfort.cfg' (or \$OUT if set)
    FC          - fortran compiler command
    SYSFFLAGS   - fortran compiler flags required always
    SYSLDFLAGS  - fortran library flags required always
    SYSLDLIBS   - fortran libraries required always
    MOD         - suffix of module dependency files (usually .mod)
    MODPATHFLAG - flag for module files search path (usually -I)
    F90MODTYPE  - one of (std.mod HP.mod fqs.mod intel.d hitachi.f90)
    IFCCEM      - PUC file editor, set non-null if F90MODTYPE=intel.d
    
END_OF_HELP
	exit 1
}

#
# Step 1: initialization; read chkfort.cfg if exists
#

output=${OUT:-chkfort.cfg}

killcache=no
for arg in "$@"; do
	case $arg in
	-reinit|-nocache)
		killcache=yes
		;;
	-help)
		help
		;;
	-h)
		help
		;;
	esac
done
[ X"$killcache" = X"yes" ] && [ -f $output ] && rm -f $output

if [ -f $output ] && [ X"$FC" = X"" ]; then
	echo previous configuration result \`$output\' used.
	. ./$output
fi

#
# Step 2: determine FC;  search if not specified.
#

if [ X"$FC" = X"" ] && [ X"$F90MODTYPE" = X"" ]; then
    if [ X"$FC" = X"" ]; then
	FC=none
	for fc in frt f95 f90 xlf cf95 cf90 lf95 lf90 g95 ifort ifc
	do
	    fcp=`which $fc 2> /dev/null` || continue
	    [ -x $fcp ] || continue
	    FC=$fc
	    break
	done
	if [ X"$FC" = X"none" ]; then
	    echo Fortran Compiler not found. If you have, specify with \$FC.
	fi
    fi

    cat <<END_OF_TEST_PROGRAM > conftest.f90
program conftest
print *, "Hello, World!"
end program conftest
END_OF_TEST_PROGRAM

    if $FC conftest.f90 > /dev/null 2>&1 && [ -x a.out ] ; then
	echo your fortran compiler is \`$FC\' and it works.
    else
	cat <<END_OF_TEST_F77 > conftest.f
      PROGRAM CONFTEST
      PRINT *, 'HELLO, WORLD'
      END PROGRAM
END_OF_TEST_F77
	if $FC conftest.f > /dev/null 2>&1 && [ -x a.out ] ; then
	    echo $FC is only FORTRAN 77 compiler. This is no good.
	else
	    echo $FC is NOT a fortran compiler
	fi
	rm -f conftest.* a.out work.pc*
	exit 1
    fi
    rm -f conftest.* a.out
else
    echo your fortran compiler is \`$FC\'.
fi

#
# Step 3: check module handling type
#

: ${F90MODTYPE:=undef}

if [ X"$F90MODTYPE" = X"undef" ]; then
	:
	echo I will examine how $FC pass module info to another files.
	rm -f conftes?.* CONFTES?.* work.pc*
	cat > conftes1.f90 <<EOF
module conftesa
logical:: b = .false.
end module conftesa
EOF
	$FC -c conftes1.f90 > /dev/null 2>&1
	if [ ! -f conftes1.o ] && [ X"$FC" = X"xlf90" ]; then
		$FC -qsuffix=f=f90 -c conftes1.f90 > /dev/null 2>&1
	fi
	if [ ! -f conftes1.o ]; then
		report_me $FC conftes1 failed.
	fi
	if [ -f conftes1.d ]; then
		F90MODTYPE=intel.d
	elif [ -f CONFTESA.mod ]; then
		F90MODTYPE=HP.mod
	elif [ -f conftesa.mod ]; then
		F90MODTYPE=std.mod
	else
		cat > conftes2.f90 <<EOF
program conftes2
use conftesa, only: b
b = .true.
end
EOF
		ln conftes1.f90 conftesa.f90
		if $FC -c conftes2.f90 > /dev/null 2>&1 && [ -f conftes2.o ]
		then
			F90MODTYPE=hitachi.f90
		elif $FC -c -Am conftes1.f90 && $FC -c -Am conftes2.f90 > /dev/null 2>&1
		then
			F90MODTYPE=fqs.mod
		else
			report_me unknown module system
		fi
	fi
	rm -f conftes?.* CONFTES?.* work.pc*
fi

#
# Step 4: check for F90MODTYPE
#

case "${F90MODTYPE:-undef}" in
intel.d)
	MOD=.d
	: ${IFCCEM:=ifccem}
	;;
HP.mod)
	MOD=.mod
	;;
std.mod)
	MOD=.mod
	;;
hitachi.f90)
	MOD=.f90
	;;
fqs.mod)
	case "$SYSFFLAGS" in
		*-Am*) : ;;
		*) SYSFFLAGS="-Am $SYSFFLAGS" ;;
	esac
	MOD=.mod
	;;
esac

#
# Step 5: Set default MODPATHFLAG
#
MODPATHFLAG=-I


#
# Step 6: special considerations to some environments
#

: ${SYSTYPE:=`uname -s`-`uname -r`-`uname -m`}

case $SYSTYPE in
HI-OSF*-*)
	SYSFFLAGS="$SYSFFLAGS -i,E,L,EU"
        ;;
HP-UX*)
	SYSLDLIBS="$SYSLDLIBS -lU77"
	;;
SunOS-5*)
	libpath=`which $FC`
	libpath=`dirname $libpath`
	libpath=`dirname $libpath`/lib
	SYSLDFLAGS="$SYSLDFLAGS -R$libpath"
	;;
esac

if [ X"$IFCCEM" != X"" ]; then
	IFCCEM=`which $IFCCEM`
	if [ X"$IFCCEM" = X"" ]; then
		echo IFCCEM command not found.
		exit 1
	fi
fi

#
# Step 7: special considerations to HITACHI Optimizing FORTRAN90
#         Compiler on SR11000
#
cat <<END_OF_TEST_PROGRAM > conftes3.f90
subroutine CONFTESB(c)
  integer, intent(inout) :: c
  inout = inout*2
end subroutine CONFTESB
END_OF_TEST_PROGRAM

cat <<END_OF_TEST_PROGRAM > conftes4.f90
program conftes4
  integer :: d
  d = 1
  call conftesb(d)
end
END_OF_TEST_PROGRAM

if $FC -c -i,L conftes3.f90 > /dev/null 2>&1 && [ -f conftes3.o ] ; then
    if ! $FC conftes4.f90 conftes3.o > /dev/null 2>&1 || [ ! -x a.out ] ; then
	if $FC -i,L conftes4.f90 conftes3.o > /dev/null 2>&1 || [ -x a.out ] ; then
	    echo $FC is HITACHI Optimizing FORTRAN90 Compiler
	    SYSFFLAGS="$SYSFFLAGS -nohugeary -i,L -parallel -allocinline -hf95"
	fi
    fi
fi
rm -f conftes?.* a.out work.pc*

#
# Step 8: special considerations to Absoft Pro Fortran
#         Compiler on Macintosh OS X
#
cat <<END_OF_TEST_PROGRAM > conftes5.f90
subroutine CONFTES5(c)
  integer, intent(inout) :: c
  inout = inout*2
end subroutine CONFTES5
END_OF_TEST_PROGRAM

cat <<END_OF_TEST_PROGRAM > conftes6.f90
program conftes6
  integer :: d
  d = 1
  call conftes5(d)
end
END_OF_TEST_PROGRAM

if $FC -c -YALL_NAMES=LCS conftes5.f90 > /dev/null 2>&1 && [ -f conftes5.o ] ; then
    if ! $FC conftes6.f90 conftes5.o > /dev/null 2>&1 || [ ! -x a.out ] ; then
	if $FC -YALL_NAMES=LCS conftes6.f90 conftes5.o > /dev/null 2>&1 || [ -x a.out ] ; then
	    echo $FC is Absoft Pro Fortran 90/95 Compiler
	    SYSFFLAGS="$SYSFFLAGS -YALL_NAMES=LCS -YEXT_SFX=_"
	    SYSLDLIBS="$SYSLDLIBS -lU77"
	    MODPATHFLAG=-p
	fi
    fi
fi
rm -f conftes?.* a.out work.pc*

#
# Last step: output
#

cat <<END_OF_REPORT > $output
 FC="$FC"			; export FC
 SYSFFLAGS="$SYSFFLAGS"		; export SYSFFLAGS
 SYSLDFLAGS="$SYSLDFLAGS"	; export SYSLDFLAGS
 SYSLDLIBS="$SYSLDLIBS"		; export SYSLDLIBS
 MOD="$MOD"			; export MOD
 MODPATHFLAG="$MODPATHFLAG"     ; export MODPATHFLAG
 F90MODTYPE="$F90MODTYPE"	; export F90MODTYPE
 IFCCEM="$IFCCEM"		; export IFCCEM
 SYSTYPE="$SYSTYPE"		; export SYSTYPE
END_OF_REPORT
echo my guess about the fortran compiler is written onto $output.
exit 0
