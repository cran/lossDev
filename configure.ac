AC_INIT([lossDev], [0.8.0])
#AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)


#AC_PROG_CXX

AC_PREFIX_PROGRAM(jags)

#Set optional include path with option --with-jags-include or
#environment variable JAGS_INCLUDE

AC_ARG_WITH([jags-include],
            AC_HELP_STRING([--with-jags-include=INCLUDE_PATH],
                           [the location of JAGS header files]),
            [jags_include_path=$withval])
if test [ -n "$jags_include_path" ] ; then
   JAGS_INCLUDE=${jags_include_path}
else
   if test [ -z "${JAGS_INCLUDE}" ] ; then
      if test [ "$prefix" == "NONE" ]; then
         AC_MSG_ERROR("JAGS headers not found. Use --with-jags-include configure option")
      else
         JAGS_INCLUDE=${prefix}/include/JAGS
      fi
   fi
fi

# There is no consensus on where 64- versus 32-bit libraries should be
# installed.  We take a hint from R, assuming that if R is installed
# in, say, /usr/LIBnn/R then we can find the jags library in
# ${prefix}/LIBnn.
#
# This is a reasonable assumption if both JAGS and R have been
# properly packaged, but may not be true if either has been manually
# installed. The user can set the environment variable LIBnn to
# override this.

if test -z "${LIBnn}"; then 
   if test -z "${R_HOME}"; then
      echo "R_HOME is not defined"
      exit 1
   fi
   LIBnn=`${R_HOME}/bin/R CMD config LIBnn`
fi

#Set optional library path with option --with-jags-lib or
#environment variable JAGS_LIBS

AC_ARG_WITH([jags-lib],
            AC_HELP_STRING([--with-jags-lib=LIB_PATH],
                           [the location of the JAGS library]),
            [jags_lib_path=$withval])
if test [ -n "$jags_lib_path" ] ; then
   JAGS_LIB=$jags_lib_path
else
   if test [ -z "${JAGS_LIB}" ] ; then
      if test [ "${prefix}" == "NONE" ] ; then 
         AC_MSG_ERROR("Location of JAGS library not defined. Use --with-jags-lib configure option")
      else
         JAGS_LIB=${prefix}/${LIBnn}
      fi
   fi
fi



dnl fortran stuff
ACX_F77_WRAPPERS
ACX_F77_LIBRARY_LDFLAGS
ACX_LAPACK([],[
	AC_MSG_ERROR("You need to install the LAPACK library")
])

AC_SUBST(JAGS_INCLUDE)
AC_SUBST(JAGS_LIB)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LIBS)
AC_SUBST(FLIBS)
#AC_SUBST(CXXFLAGS)

AC_OUTPUT(src/Makevars)
