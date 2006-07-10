
# AC_OSI_CANONICAL(solver)

# Given a solver name (e.g., clp, all letters in lower case), set four variables:
# osi_lc_solver: solver in lower case (e.g., clp)
# osi_mc_solver: solver in mixed case (e.g., Clp)
# osi_uc_solver: solver in upper case (e.g., CLP)
# osi_exists_solver: yes/no according to whether solver is present or not,
#	based on presence/absence of underlying project.

# These variants are useful in various contexts: lc for configure variables,
# mc for OsiXxxSolverInterface, uc for COIN_HAS_XXX.

AC_DEFUN([AC_OSI_CANONICAL],
[ AC_REQUIRE([AC_COIN_HAS_PROJECT])

# Convert to canonical values. Unfortunately, there's no easy way to do this in
# sh and/or sed. Assuming the existence of `tr' seems risky --- highly probable
# in any environment that supports configure, but not guaranteed.  If this list
# starts to get any longer, tr deserves a harder look. Not to mention that m4
# eats square braces for all three meals, hence the ugly quadrigraphs: @<:@ for
# open square bracket, @:>@ for close square bracket.  Of course, that last
# sentence will be completely opaque to people reading the configure script,
# 'cause the quadrigraphs will be replaced with brackets, even in a comment!

  case $1 in
    clp)
      osi_lc_solver=clp ;
      osi_mc_solver=Clp ;
      osi_uc_solver=CLP ;
      osi_exists_solver=${coin_has_clp-"unavailable"}
      ;;
    dylp)
      osi_lc_solver=dylp ;
      osi_mc_solver=Dylp
      osi_uc_solver=DYLP
      osi_exists_solver=${coin_has_dylp-"unavailable"}
      ;;
    cpx)
      osi_lc_solver=cpx ;
      osi_mc_solver=Cpx
      osi_uc_solver=CPX
      osi_exists_solver=${coin_has_cpx-"unavailable"}
      ;;
    fmp)
      osi_lc_solver=fmp ;
      osi_mc_solver=Fmp
      osi_uc_solver=FMP
      osi_exists_solver=${coin_has_fmp-"unavailable"}
      ;;
    glpk)
      osi_lc_solver=glpk ;
      osi_mc_solver=Glpk
      osi_uc_solver=GLPK
      osi_exists_solver=${coin_has_glpk-"unavailable"}
      ;;
    msk)
      osi_lc_solver=msk ;
      osi_mc_solver=Msk
      osi_uc_solver=MSK
      osi_exists_solver=${coin_has_msk-"unavailable"}
      ;;
    osl)
      osi_lc_solver=osl ;
      osi_mc_solver=Osl
      osi_mc_solver=OSL
      osi_exists_solver=${coin_has_osl-"unavailable"}
      ;;
    spx)
      osi_lc_solver=spx ;
      osi_mc_solver=Spx
      osi_uc_solver=SPX
      osi_exists_solver=${coin_has_spx-"unavailable"}
      ;;
    sym)
      osi_lc_solver=sym ;
      osi_mc_solver=Sym
      osi_uc_solver=SYM
      osi_exists_solver=${coin_has_sym-"unavailable"}
      ;;
    vol)
      osi_lc_solver=vol ;
      osi_mc_solver=Vol
      osi_uc_solver=VOL
      osi_exists_solver=${coin_has_vol-"unavailable"}
      ;;
    xpr)
      osi_lc_solver=xpr ;
      osi_mc_solver=Xpr
      osi_uc_solver=XPR
      osi_exists_solver=${coin_has_xpr-"unavailable"}
      ;;
    cbc)
      osi_lc_solver=cbc ;
      osi_mc_solver=Cbc
      osi_uc_solver=CBC
      osi_exists_solver=${coin_has_cbc-"unavailable"}
      ;;
    *)
      osi_lc_solver=clp ;
      osi_mc_solver=Clp ;
      osi_uc_solver=CLP ;
      osi_exists_solver=${coin_has_clp-"unavailable"}
      AC_MSG_WARN([Unrecognised solver $1; defaulting to $osi_lc_solver.])
      ;;
  esac

# Now that we have something in standard form, turn the result on solver
# availability into something easy to test.

  if test $osi_exists_solver = unavailable || \
     test $osi_exists_solver = skipping || \
     test $osi_exists_solver = false ; then
    osi_exists_solver=no
  else
    osi_exists_solver=yes
  fi
])
