
# AC_ODSI_PARANOIA(dflt)
# ------------------------------------------------------------------
# Processes the paranoia option.
# ------------------------------------------------------------------
AC_DEFUN([AC_ODSI_PARANOIA],
[ AC_ARG_ENABLE([osidylp-paranoia],
      AS_HELP_STRING([--enable-osidylp-paranoia],
          [Control osidylp's paranoid checks. 0 is off, 1 is index checks,
	   2 is consistency checks (expensive) (default=$1)]),
      [osidylp_paranoia=$enableval],
      [osidylp_paranoia=$1])
  if test "$osidylp_paranoia" = "yes"; then
    osidylp_paranoia=1
  elif test "$osidylp_paranoia" = "no"; then
    osidylp_paranoia=0
  fi
  AC_DEFINE_UNQUOTED([ODSI_PARANOIA],[$osidylp_paranoia],
	[Control OsiDylp's paranoid checks. Legal values: 0 - off; 1 - normal;
	 2 - consistency (expensive)])
  case $osidylp_paranoia in
    0)
      AC_MSG_NOTICE([OsiDylp paranoid checks disabled.])
      ;;
    1)
      AC_MSG_NOTICE([OsiDylp paranoid checks at normal level.])
      ;;
    *)
      AC_MSG_NOTICE([OsiDylp paranoid checks at level $osidylp_paranoia.])
      ;;
  esac
])

# AC_ODSI_STATISTICS(dflt)
# ------------------------------------------------------------------
# Processes the statistics option.
# ------------------------------------------------------------------
AC_DEFUN([AC_ODSI_STATISTICS],
[ AC_ARG_ENABLE([osidylp-stats],
      AS_HELP_STRING([--enable-osidylp-stats],
          [Enable OsiDylp's support for dylp's statistics collection
	   features (default=$1)]),
      [osidylp_stats=$enableval],
      [osidylp_stats=$1])
  if test "$osidylp_stats" = "yes"; then
    AC_DEFINE([ODSI_STATISTICS],[1],
	[Define this variable to enable support for dylp's statistics
	 collection features.])
    AC_MSG_NOTICE([OsiDylp support for dylp statistics collection enabled.])
  else
    AC_MSG_NOTICE([OsiDylp support for dylp statistics collection disabled.])
  fi
])

# AC_ODSI_INFO(dflt)
# ------------------------------------------------------------------
# Processes the information printing (info) option.
# ------------------------------------------------------------------
AC_DEFUN([AC_ODSI_INFO],
[ AC_ARG_ENABLE([osidylp-info],
      AS_HELP_STRING([--enable-osidylp-info],
          [Enable OsiDylp's informational printing (output will depend on log
	   level) (default=$1)]),
      [osidylp_info=$enableval],
      [osidylp_info=$1])
  if test "$osidylp_info" = "yes"; then
    AC_DEFINE([ODSI_INFOMSGS],[1],
	[Define this variable to enable OsiDylp's informational printing
	 features.])
    AC_MSG_NOTICE([OsiDylp informational printing enabled.])
  else
    AC_MSG_NOTICE([OsiDylp informational printing disabled.])
  fi
])
