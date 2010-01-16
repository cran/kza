AC_PREREQ(2.50)
AC_INIT([kza], [1.02], [brian.close@gmail.com])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
 echo "could not determine R_HOME"
 exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

AC_PROG_CC

AC_CHECK_LIB(fftw3, fftw_execute, [], [], [-lm])

# Checks for header files.
AC_HEADER_STDC

AC_CHECK_HEADERS([fftw3.h], [HAVE_FFTW=TRUE; AC_DEFINE(HAVE_FFTW)], [HAVE_FFTW=FALSE; AC_MSG_WARN([Can't find fftw3 header])])

AC_SUBST(LIBS)
AC_SUBST(HAVE_FFTW)

AC_CONFIG_HEADERS([src/config.h])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
