
# AC_OGSI_MIPSOLVER()
# ------------------------------------------------------------------
# Selects the MIP solver used by OsiGlpk. Newer versions of glpk provide
# lpx_intopt. For older versions, we'll fall back on lpx_integer. All that's
# done here is to check for lpx_integer and define the compilation symbol
# GLPK_HAS_INTOPT. At least as far back as glpk 4.4, lpx_intopt is present.
# ------------------------------------------------------------------
AC_DEFUN([AC_OGSI_MIPSOLVER],
[ AC_MSG_CHECKING([for glpk mip solver ])
  osiglpk_candidates="_glp_lpx_intopt glp_lpx_intopt _glp_lpx_integer glp_lpx_integer"
  coin_save_LIBS="$LIBS"
  LIBS="$ADDLIBS $LIBS"
  osiglpk_mipsolver=
  AC_LANG_PUSH(C)
  for fnm in $osiglpk_candidates ; do
#   AC_MSG_CHECKING([whether symbol $fnm is available with glpk])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[void $fnm() ;]],[[$fnm()]])],
        [osiglpk_mipsolver=$fnm ; break], [])
  done
  AC_LANG_POP(C)
  LIBS="$coin_save_LIBS"
  if test x"$osiglpk_mipsolver" = x ; then
    AC_MSG_RESULT([none!])
    AC_MSG_ERROR([No glpk MIP solver! You really should upgrade glpk.])
  else
    case $osiglpk_mipsolver in
      *intopt)
	AC_MSG_RESULT([lpx_intopt])
	AC_DEFINE(GLPK_HAS_INTOPT,[1],
	    [Define to 1 if GLPK has the advanced B&B solver lpx_intopt])
	break ;;
      *integer)
	AC_MSG_RESULT([lpx_integer])
	break ;;
    esac
  fi
])


