#ifdef _WIN32
#  ifdef DLL_EXPORT
#    ifdef OSILIB_BUILD  /* build of Osi DLL */
#      define OSILIB_EXPORT __declspec(dllexport)
#    else  /* use of Osi DLL */
#      define OSILIB_EXPORT __declspec(dllimport)
#    endif
#  endif
#  ifdef OSILIB_DLLIMPORT  /* alternative to indicate using Osi DLL */
#    define OSILIB_EXPORT __declspec(dllimport)
#  endif
#endif

#ifndef OSILIB_EXPORT
#  define OSILIB_EXPORT
#endif
