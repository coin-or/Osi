#ifdef _WIN32
#  ifdef OSILIB_BUILD
#    ifdef DLL_EXPORT  /* build of Osi DLL */
#      define OSILIB_EXPORT __declspec(dllexport)
#    else  /* build of Osi LIB */
#      define OSILIB_EXPORT
#    endif
#  elif defined(OSILIB_DLLIMPORT)  /* using Osi DLL */
#    define OSILIB_EXPORT __declspec(dllimport)
#  else   /* using Osi LIB */
#    define OSILIB_EXPORT
#  endif
#else  /* Unix */
#    define OSILIB_EXPORT
#endif
