/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Include file for the configuration of OsiCommonTest.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header file is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these files and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h
 * files.
 */

#ifndef __OSIUNITTESTSCONFIG_H__
#define __OSIUNITTESTSCONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef OSICOMMONTEST_BUILD

#ifdef DLL_EXPORT
#define OSICOMMONTESTLIB_EXPORT __declspec(dllexport)
#endif

#else
#include "config_osicommontest.h"
#endif

#else /* HAVE_CONFIG_H */


#ifndef OSICOMMONTESTLIB_EXPORT
#ifdef _WIN32
#ifdef OSILIB_BUILD
#define OSICOMMONTESTLIB_EXPORT __declspec(dllexport)
#else
#define OSICOMMONTESTLIB_EXPORT __declspec(dllimport)
#endif
#endif
#endif


#endif /* HAVE_CONFIG_H */

#ifndef OSICOMMONTESTLIB_EXPORT
#define OSICOMMONTESTLIB_EXPORT
#endif

#endif /*__OSICONFIG_H__*/