
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_osi_default.h */
#ifndef OSILIB_EXPORT
#ifdef _WIN32
/* assuming we build an Osi DLL */
#define OSILIB_EXPORT __declspec(dllexport)
#else
#define OSILIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_osi_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to the debug sanity check level (0 is no test) */
#define COIN_OSI_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_OSI_VERBOSITY 0

/* Define to 1 if the CoinUtils package is used.
 * Don't undef this unless you really know what you're doing.
 */
#define COIN_HAS_COINUTILS 1

/* Define to 1 if the Cplex package is used */
/* #define COIN_HAS_CPLEX 1 */

/* Define to 1 if the Glpk package is used */
/* #define COIN_HAS_GLPK 1 */

/* Define to 1 if the Gurobi package is used */
/* #define COIN_HAS_GUROBI 1 */

/* Define to 1 if the Mosek package is used */
/* #define COIN_HAS_MOSEK 1 */

/* Define to 1 if the SoPlex package is used */
/* #define COIN_HAS_SOPLEX 1 */

/* Define to 1 if the Xpress package is used */
/* #define COIN_HAS_XPRESS 1 */
