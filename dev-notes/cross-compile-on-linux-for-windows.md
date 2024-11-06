Notes on cross-compiling VFGEN on Linux for Windows
===================================================

The dependencies must be cross-compiled before VFGEN can be built.

mxml
----

This was using mxml-3.3.tar.gz.

Before building, remove the definitions of the macro 'FORTIFY_SOURCE'
(or whatever it is called) from the files `configure` and `configure.ac`.
Then run

```
./configure --host=x86_64-w64-mingw32 --prefix=$HOME/local-win --disable-shared --disable-threads
make
make install
```

cln
---

Use the git source `git://www.ginac.de/cln.git/`, not a released tar file.
In the git repo, run

```
./autogen.sh
./configure --host=x86_64-w64-mingw32 --prefix=$HOME/local-win --disable-shared
make
make install
```

ginac
-----

Using `ginac-1.8.3.tar.bz2`.

```
export CLN_CFLAGS=-I$HOME/local-win/include
export CLN_LIBS="-L$HOME/local-win/lib -lcln"
./configure --host=x86_64-w64-mingw32 --prefix=$HOME/local-win --disable-shared
make
make install
```


vfgen
-----

In the vfgen git repo, run the commands that create `configure`
(given in README.md), then run `./configure` and `make` as follows:

```
export MXML_CFLAGS=-I$HOME/local-win/include
export MXML_LIBS="-L$HOME/local-win/lib -lmxml"
export GINAC_CFLAGS=-I$HOME/local-win/include
export GINAC_LIBS="-L$HOME/local-win/lib -lginac -lcln"
./configure --host=x86_64-w64-mingw32 --prefix=$HOME/local-win
make LDFLAGS=-static
```

The executable `vfgen.exe` is in the `src/` subdirectory:

```
$ wine ./src/vfgen.exe help
VFGEN (Version:2.6.0.dev6)
Use: vfgen command  vector-field-file.vf
or:  vfgen command:option=value,...,option=value vector-field-file.vf
or:  vfgen help command
where command is one of:
    adolc, auto, boostodeint, check, cvode7, dde23, dde_solver, ddebiftool,
    delay2ode, evf, gsl, help, javamath, javascript, julia, latex,
    lsoda, matcont, matlab, octave, pygsl, r, radau5, scilab,
    scipy, taylor, xml, xpp
```
