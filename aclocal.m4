dnl------------------------------------------------------------------
dnl      Copyright (C) GFD Dennou Club, 2005. All rights reserved.
dnl------------------------------------------------------------------
dnl
dnl= Definitions of Functions for configure
dnl
dnl Authors:: Masatsugu Odaka
dnl Version:: $Id: aclocal.m4,v 1.1 2011/04/26 04:08:09 sugiyama Exp $
dnl Tag Name:: $Name:  $
dnl Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
dnl License::   See COPYRIGHT[link:COPYRIGHT]
dnl
dnl == History
dnl 
dnl * 2006/09/28 (ODAKA Masatsugu): Initial release, copy from gt4f90io library
dnl
dnl == DC_ARG_WITH(withname, description, varname, ifnot)
dnl
dnl
AC_DEFUN(DC_ARG_WITH, [
	AC_ARG_WITH($1, [  --with-$1=ARG: $2], [
		$3=$withval
	], [
		AC_CACHE_CHECK([$2], $3, [$4])
	])
])

dnl
dnl == DC_ARG_ENABLE(feature, description, varname, ifnot)
dnl
AC_DEFUN(DC_ARG_ENABLE, [
	AC_ARG_ENABLE($1, [  --enable-$1: $2], [
		$3=$enableval
	], [
		AC_CACHE_CHECK([$2], $3, [$4])
	])
])

dnl
dnl == Set DEEPCONVDIR (current working directory)
dnl
dnl  usage: DC_SET_DEEPCONVDIR
dnl
AC_DEFUN(DC_SET_DEEPCONVDIR,
[
dnl     AC_PATH_PROG(PWD, pwd)
    case "${PWD-}" in
        '')
            AC_MSG_WARN(*** WARN *** Please set environment variable GT4DIR by yourself)
            DEEPCONVDIR=
            ;;
        *)
            DEEPCONVDIR=${PWD-}
    esac
    AC_SUBST(DEEPCONVDIR)
])

dnl
dnl == Check libfile and set LIBDIR and  LIBNAME
dnl
dnl  Check existence of "libfile" file, and set LIBDIR, LIBNAME
dnl  from "libfile".
dnl
dnl  usage: DC_SET_LIBDIR_LIBNAME(libfile, LIBDIR, LIBNAME)
dnl
AC_DEFUN(DC_SET_LIBDIR_LIBNAME, [
	if test ! -f $1 ; then
		AC_MSG_ERROR(specified library file \"$1\" is not exist)
	fi

	$2=`dirname $1`
	$3=`basename $1 .a | sed 's/^lib//'`
])

dnl
dnl == Modify INSTALL (if "./install-sh", set absolute path)
dnl
dnl  usage: DC_MOD_INSTALL
dnl
AC_DEFUN(DC_MOD_INSTALL,
[
dnl     AC_PATH_PROG(PWD, pwd)
    case "${INSTALL}" in
        './'*)
            case "${PWD-}" in
                '')
                    AC_MSG_WARN(*** WARN *** Please set environment variable INSTALL with absolute path by yourself)
                    ;;
                *)
                    INSTALL=${PWD}/$INSTALL
                    AC_SUBST(INSTALL)
                    AC_MSG_NOTICE(*** MSG *** environment variable INSTALL is reconfigured with absolute path)
            esac
            ;;
        *) ;;
    esac
])
