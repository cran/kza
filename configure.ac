
AC_INIT([src/kza.c])

AC_PATH_PROG([FFTW_CONFIG], [ftw-config])
## If fftw-config was found, let's use it
if test "${FFTW_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   FFTW_CFLAGS=`${FFTW_CONFIG} --cflags`
   FFTW_LIBS=`${FFTW_CONFIG} --libs`
else
   AC_MSG_ERROR([fftw-config not found, is FFTW installed?])
fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(FFTW_CFLAGS) 
AC_SUBST(FFTW_LIBS)
AC_OUTPUT(src/Makevars)

