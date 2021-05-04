# ZLIB_CHECK_CONFIG ([DEFAULT-ACTION] [MINIMUM-VERSION],
#                    [ACTION-IF-YES], [ACTION-IF-NO])
#
AC_DEFUN([ZLIB_CHECK_CONFIG],
[
  AC_ARG_WITH(zlib,
	      AC_HELP_STRING([--with-zlib=DIR],[look for zlib in DIR]),
	      [_zlib_with=$withval],[_zlib_with="no"])

  ZLIB_ROOT=""
  if test "$_zlib_with" != "no"
  then
     if test -f "$_zlib_with/include/zlib.h"
     then
         ZLIB_ROOT=$_zlib_with
     fi
  fi

  # Check if it's a working library
  zlib_ok=no
  if test "$ZLIB_ROOT" != ""
  then
    _cppflags=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS -I${ZLIB_ROOT}/include"
    _ldflags=$LDFLAGS
    LDFLAGS="$LFDLAGS -L${ZLIB_ROOT}/lib"
    AC_LANG_SAVE
    AC_LANG_C
    AC_CHECK_LIB(z, inflateEnd,
	[AC_CHECK_HEADER(zlib.h, zlib_ok=yes, zlib_ok=no)])
    AC_LANG_RESTORE
    if test "$zlib_ok" != "yes"
    then
        # Backout and whinge
        CPPFLAGS=$_cppflags
        LDFLAGS=$_ldflags
        AC_MSG_WARN("--with-zlib specified, but non functioning")
    fi

  else
    # Maybe it works "out of the box"?
    AC_CHECK_LIB(z, inflateEnd,
	[AC_CHECK_HEADER(zlib.h, zlib_ok=yes, zlib_ok=no)])
  fi

  if test "$zlib_ok" = "yes"
  then
      AC_DEFINE(HAVE_ZLIB, 1,
         [Define to 1 if you have a functional libz.])
      if test "$ZLIB_ROOT" != ""
      then
          LIBZ="-L${ZLIB_ROOT}/lib -lz"
      else
          LIBZ=-lz
      fi
      AC_SUBST(LIBZ)
  else
    AC_MSG_WARN("No functioning zlib found")
  fi

  # Not sure how many of these are needed, but it's belt-and-braces mode
  AH_TEMPLATE([HAVE_ZLIB], [Define if zlib is installed])
  AM_CONDITIONAL(HAVE_ZLIB, test "$zlib_ok" = "yes")
])

