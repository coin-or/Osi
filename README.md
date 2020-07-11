# Osi

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Osi?sort=semver)](https://github.com/coin-or/Osi/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](https://github.com/coin-or/coinbrew/tree/master/scripts/generate_readme).
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation script._

The COIN-OR Open Solver Interface is a uniform API for interacting
with callable solver libraries.
It supports linear programming solvers as well as the ability to \"finish off\"
a mixed-integer problem calling the solver library\'s MIP solver.
A list of supported solvers appears at the bottom of the page.

Osi is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/eclipse-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org)

The Osi website is https://github.com/coin-or/Osi.

## CITE

[![DOI](https://zenodo.org/badge/173476455.svg)](https://zenodo.org/badge/latestdoi/173476455)

## CURRENT BUILD STATUS

[![Build Status](https://travis-ci.org/coin-or/Osi.svg?branch=master)](https://travis-ci.org/coin-or/Osi)

[![Build status](https://ci.appveyor.com/api/projects/status/frpvf6totmchmjmv/branch/master?svg=true)](https://ci.appveyor.com/project/tkralphs/osi-x2d8y/branch/master)

## DOWNLOAD

Binaries for most platforms are available as part of [Cbc](https://bintray.com/coin-or/download/Cbc). 

 * *Linux*: On Debian/Ubuntu, CoinUtils is available in the package `coinor-osi` and can be installed with apt. On Fedora, Osi is available in the package `coin-or-Osi`.
 * *Windows*: The easiest way to get Osi on Windows is to download from *[Bintray](https://bintray.com/coin-or/download/Cbc)*.
 * *Mac OS X*: The easiest way to get Cbc on Mac OS X is through [Homebrew](https://brew.sh).
   * `brew tap coin-or-tools/coinor`
   * `brew install osi`

Due to license incompatibilities, pre-compiled binaries lack some functionality.
If binaries are not available for your platform for the latest version and you would like to request them to be built and posted, feel free to let us know on the mailing list.

*Source code* can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Osi from the
 [releases](https://github.com/coin-or/Osi/releases) page.
 * Cloning the repository from [Github](https://github.com/coin-or/Osi) or using the 
`coinbrew` script (recommended).  

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
./coinbrew fetch Osi@stable/2.10
./coinbrew build Osi
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

Obtain the source code, e.g., by cloning the git repo https://github.com/coin-or/Osi
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
directory will be `./build/Osi/version` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Osi/Doxygen).

## Project Links

 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [Mailing list](http://list.coin-or.org/mailman/listinfo/osi)
 * [Report a bug](https://github.com/coin-or/Osi/issues/new)
 * [Doxygen-generated html documentation](http://coin-or.github.io/Osi/Doxygen)

