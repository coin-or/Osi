# Osi

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Osi?sort=semver)](https://github.com/coin-or/Osi/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme).
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation script._

Osi (*O*pen *S*olver *I*nterface) provides an abstract base class to a generic linear programming (LP) solver, along with derived classes for specific solvers.
Many applications may be able to use the Osi to insulate themselves from a specific LP solver.
That is, programs written to the OSI standard may be linked to any solver with an OSI interface and should produce correct results.
The OSI has been significantly extended compared to its first incarnation.
Currently, the OSI supports linear programming solvers and has rudimentary support for integer programming.
Among others the following operations are supported:
 * creating the LP formulation;
 * directly modifying the formulation by adding rows/columns;
 * modifying the formulation by adding cutting planes provided by [CGL](https://www.github.com/coin-or/Cgl);
 * solving the formulation (and resolving after modifications);
 * extracting solution information;
 * invoking the underlying solver's branch-and-bound component.

The following is a list of derived Osi classes:

|Solver|Derived Class|Note|
|------|-------------|----|
|[Cbc](https://www.github.com/coin-or/Cbc)|OsiCbc| unmaintained | 
|[Clp](https://www.github.com/coin-or/Clp)|OsiClp| |
|[CPLEX](https://www.ibm.com/analytics/cplex-optimizer)|OsiCpx| |
|[DyLP](https://www.github.com/coin-or/DyLP)|OsiDylp| |
|[GLPK](http://www.gnu.org/software/glpk/glpk.html)|OsiGlpk| Glpk |
|[Gurobi](http://www.gurobi.com)|OsiGrb| |
|[HiGHS](https://www.github.com/coin-or/HiGHS)|OsiHiGHS| under development |
|[MOSEK](http://www.mosek.com)|OsiMsk| |
|[SoPlex](http://soplex.zib.de)|OsiSpx| SoPlex < 4.0 |
|[SYMPHONY](https://www.github.com/coin-or/SYMPHONY)|OsiSym| |
|[Vol](https://www.github.com/coin-or/Vol)|OsiVol| |
|[XPRESS-MP](https://www.fico.com/en/products/fico-xpress-optimization)|OsiXpr| |

Each solver interface is in a separate directory of Osi or distributed
with the solver itself.

Within COIN-OR, Osi is used by [Cgl](https://www.github.com/coin-or/Cgl), [Cbc](https://www.github.com/coin-or/Cbc), and [Bcp](https://www.github.com/coin-or/Bcp), among others.

The main project managers are Lou Hafer (@LouHafer) and Matt Saltzmann (@mjsaltzman).

An incomplete list of recent changes to Osi are found in the [CHANGELOG](Osi/CHANGELOG)


Osi is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org)

The Osi website is https://github.com/coin-or/Osi.

## CITE

[![DOI](https://zenodo.org/badge/173476455.svg)](https://zenodo.org/badge/latestdoi/173476455)

## CURRENT BUILD STATUS

[![Build Status](https://travis-ci.org/coin-or/Osi.svg?branch=master)](https://travis-ci.org/coin-or/Osi)

[![Build status](https://ci.appveyor.com/api/projects/status/frpvf6totmchmjmv/branch/master?svg=true)](https://ci.appveyor.com/project/tkralphs/osi-x2d8y/branch/master)

## DOWNLOAD

Binaries for most platforms are available as part of [Cbc](https://bintray.com/coin-or/download/Cbc). 

 * *Linux*: On Debian/Ubuntu, Osi is available in the package `coinor-osi` and can be installed with apt. On Fedora, Osi is available in the package `coin-or-Osi`.
 * *Windows*: The easiest way to get Osi on Windows is to download from *[Bintray](https://bintray.com/coin-or/download/Cbc)*.
 * *Mac OS X*: The easiest way to get Cbc on Mac OS X is through [Homebrew](https://brew.sh).
   * `brew tap coin-or-tools/coinor`
   * `brew install coin-or-tools/coinor/osi`

Due to license incompatibilities, pre-compiled binaries lack some functionality.
If binaries are not available for your platform for the latest version and you would like to request them to be built and posted, feel free to let us know on the mailing list.

*Source code* can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Osi from the
 [releases](https://github.com/coin-or/Osi/releases) page.
 * Cloning this repository from [Github](https://github.com/coin-or/Osi) or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html).

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew`

To build CoinUtils from source, obtain the `coinbrew` script, do
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Osi@master
./coinbrew build Osi
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Osi
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxygen-docs` 

in the build directory. If Osi was built via `coinbrew`, then the build
directory will be `./build/Osi/master` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Osi/Doxygen).

## Project Links

 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [mailing list](http://list.coin-or.org/mailman/listinfo/osi)
 * [Report a bug](https://github.com/coin-or/Osi/issues/new)
 * [Doxygen-generated html documentation](http://www.coin-or.org/Doxygen/Osi/hierarchy.html)
 * [OSI2 Discussion](https://projects.coin-or.org/Osi2/wiki/Osi2Discussion)
 * The most recent tutorial on OSI can be accessed from the [page on presentations from the 2004 CORS/INFORMS Joint Meeting in Banff](http://www.coin-or.org/Presentations/CORSINFORMSWorkshop04/index.html).
 * [The COIN-OR Open Solver Interface: Technology Overview](http://www.coin-or.org/Presentations/CORS2004-OSI.pdf): An overview of the COIN-OR OSI and design issues for a next-generation version given at CORS/INFORMS 2004 by Matthew Saltzman.

-------

## Dynamically loading commercial solver libraries

### At build time

It is possible to create an osi build that supports cplex, gurobi and xpress even if you don't have (yet) any of these solvers on your machine using [lazylpsolverlibs](https://code.google.com/p/lazylpsolverlibs/). To do so, follow these steps:

 1. Install lazylpsolverlibs (follow the instructions of the [lazylpsolverlibs wiki](https://code.google.com/p/lazylpsolverlibs/wiki/HowToSetup))
 2. Use the following command line to configure Osi:
```
./configure --with-cplex-incdir="$(pkg-config --variable=includedir lazycplex)/lazylpsolverlibs/ilcplex" \
            --with-cplex-lib="$(pkg-config --libs lazycplex)" \ 
            --with-gurobi-incdir="$(pkg-config --variable=includedir lazygurobi)/lazylpsolverlibs" \
            --with-gurobi-lib="$(pkg-config --libs lazygurobi)" \
            --with-xpress-incdir="$(pkg-config --variable=includedir lazyxprs)/lazylpsolverlibs" \
            --with-xpress-lib="$(pkg-config --libs lazyxprs)"
```
 3. Then follow the normal installation process (make, make install)

### At run time

Your build should now support cplex, gurobi and xpress, which means that if you install one of these solvers, osi will be able to use it.
At run time, you just need to point one of the environment variables LAZYLPSOLVERLIBS_GUROBI_LIB, LAZYLPSOLVERLIBS_CPLEX_LIB or LAZYLPSOLVERLIBS_XPRS_LIB to the full path of the corresponding solver library.
For example:
```
export LAZYLPSOLVERLIBS_CPLEX_LIB=/usr/ilog/cplex121/bin/x86_debian4.0_4.1/libcplex121.so
```

### Troubleshooting

If pkg-config reports errors during the configure step, try modifying the PKG_CONFIG_PATH variable. Most likely, you need to do:
```
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
```

