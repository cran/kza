AC_PREREQ(2.50)
AC_INIT([kza], [1.01], [brian.close@gmail.com])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
 echo "could not determine R_HOME"
 exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

# Checks for programs.
AC_PROG_CC


# Checks for libraries.
AC_CHECK_LIB([fftw3_junk], [fftw_plan_dft], [HAVE_FFTW=TRUE], [HAVE_FFTW=FALSE])
AC_CHECK_LIB([fftw3_junk], [fftw_plan_dft],,)

# Checks for header files.
AC_HEADER_STDC

AC_CHECK_HEADERS([fftw3.h], [HAVE_FFTW_H=TRUE],
                 [HAVE_FFTW_H=FALSE; AC_MSG_WARN([Can't find fftw3 header])])

AC_SUBST(LIBS)
AC_SUBST(HAVE_FFTW)

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars R/kzft.R])
AC_OUTPUT
