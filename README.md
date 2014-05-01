VFGEN
-----

A *Vector Field File Generator* for differential equation solvers and
other computational tools.  Given a specification of an ordinary or delay
differential equation, VFGEN can generate source code that defines the
equations for a wide assortment of solvers and tools.

* Web page: `www.warrenweckesser.net/vfgen`
* License:  GPL-2
* Author:   Warren Weckesser, `www.warrenweckesser.net`
* Email:    vfgen DOT help AT gmail com  (fix the obvious minor obfuscation)

This file is in the top directory for the source code distribution of the
program VFGEN.  This README gives the steps necessary to build VFGEN from
the source code.  Before building from source, check the web page
www.warrenweckesser.net/vfgen; a compiled binary may be available for your
system.

Before you can build VFGEN, you must have installed these libraries:

* GiNaC    (`http://www.ginac.de`)
* Mini-XML (`http://www.minixml.org`)

*Note*: Version 2.5 of the Mini-XML library has a spelling error in one of
its configuration files. If you install Mini-XML version 2.5 from source,
the spelling error causes the build of VFGEN to fail. To fix this, reinstall
version 2.5 of Mini-XML  after changing `@PHTREAD_LIBS@` to `@PTHREAD_LIBS@`
in the file `mxml.pc.in`.

If you are using Linux, binary packages (.deb and .rpm) are available for
these libraries.  Check your standard repositories.  Be sure you install
the "development" versions of the packages.  In Debian-based systems, these
usually end with -dev.  Note that GiNaC depends on the arbitrary precision
numerical library CLN, so you will also have to install CLN.

Also be sure you have pkg-config (`www.freedesktop.org/software/pkgconfig`)
installed.  Most Linux systems have pkg-config available as a binary package.

There are three methods that you can use to build the program: CMake,
configure (i.e. the autotools script), or plain old Make.  These are discussed
below.


** CMake **

CMake (`http://www.cmake.org/`) is a platform-independent build tool. To build
VFGEN with CMake, you must have CMake installed, and you must have the GiNaC
(`http://www.ginac.de/`) and Mini-XML (`http://www.easysw.com/~mike/mxml/`)
libraries installed. (If you are using a debian-based Linux distribution,
packages for these libraries are available.)  If these libraries are installed
in either /usr or /usr/local, the following commands (run in the top directory)
will build VFGEN:

    $ cd cmake_build
    $ cmake ../src
    $ make

If the libraries are installed in some other location, define the environment
variables `CMAKE_INCLUDE_PATH` and `CMAKE_LIBRARY_PATH` before running cmake.

By default, the command

    $ make install

will install the executable file in `/usr/local/bin`.  To change the
installation directory, set the environment variable `CMAKE_INSTALL_PREFIX`
before running cmake.  Then `make install` will install the executable file
in `$(CMAKE_INSTALL_PREFIX)/bin`.


** configure (i.e. autotools) **

VFGEN also comes with a "configure" script.  The simplest way to use this
script is the following sequence of commands in the top directory:

    $ ./configure
    $ make
    $ make install

This will install the vfgen executable in `/usr/local/bin`.  You may change the
installation directory with the `--prefix` option to the `./configure` command,
e.g.

    $ ./configure --prefix=/opt
    $ make
    $ make install

This will put the vfgen executable in `/opt/bin`.

The configure script has several other options.  They may be listed with the
command

    $ ./configure --help


** Plain old Makefile **

If you are using Linux (or some other Unix-like system), you have the GNU C++
compiler installed (g++), you have the program pkg-config installed, and you
have installed GiNaC and Mini-XML, you can use the file `Makefile.vfgen` to
build the program:

    $ cd src
    $ make -f Makefile.vfgen

GiNaC version 1.3.4, and Mini-XML version 2.3 (and later versions) provide the
appropriate files to work with pkg-config; I'm not sure if older versions do.

The file `Makefile.vfgen` is a *very* simple Makefile; for example, there is
no `install` target.  Besides the default `vfgen` target, the only other
target of interest is `clean`; `make clean` will remove all the `.o` files.


Building from version control source
------------------------------------

The code under version control does not include the `configure` script.
To create the `configure` script:

    $ aclocal
    $ automake --add-missing
    $ autoconf
