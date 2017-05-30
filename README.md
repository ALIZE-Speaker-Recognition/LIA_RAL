<img src="http://alize.univ-avignon.fr/images/LIA_RAL.png" alt="The LIA_RAL logo" height="198" >

# LIA_RAL

*This package is part of the ALIZÉ project: <http://alize.univ-avignon.fr>*



Welcome to LIA_RAL!
-------------------

LIA_RAL is a set of utilities for Automatic Speaker Recognition developed at LIA, providing a high-level access to the [ALIZÉ platform](http://alize.univ-avignon.fr). All these utilities are based on the core `ALIZE` library, which is required in order to compile LIA_RAL.

LIA_RAL is an open project. Feel free to contact us to propose to work on it! 

In order to use Support Vector Machines, the `libsvm` library is included in this package. For more information on `libsvm`, please refer to:
<http://www.csie.ntu.edu.tw/~cjlin/libsvm/>

The `Eigen` library is also included into this package. For more information on `Eigen`, please refer to:
<http://eigen.tuxfamily.org>


FAQ
---

#### What is needed to build and install the LIA_RAL library?

`aclocal`, `autoconf`, `automake` and `libtool` are required.
This package also requires the core `ALIZE` library.


#### How to compile the LIA_RAL tools under Linux, Mac OS and Cygwin

Follow these four steps:

1. Install the `ALIZE` library.
2. Run, in this order:
	- `aclocal`
	- `automake`
	- `autoconf`
3. Then run `./configure`.
   By default, the `ALIZE` library is searched for in `../alize-core`. It may be in a folder with a different name, depending on how you downloaded the library, or if you have decided to install it in a non-default location. If so, you can specify the absolute path by using the `--with-alize=ABSOLUTE_PATH` option.
4. Finally, run `make`.

**Note:** At step 2, you may need to use `automake --add-missing` if the file `compile` cannot be found.


#### How to compile the LIA_RAL tools under Windows with Visual C++

Use the `LIA_RAL.sln` solution file. If your version of Visual Studio is newer than 2010, you just have to convert it using Visual studio tools.


#### I have problems compiling LIA_RAL

Please read this technote:
<http://alize.univ-avignon.fr/mediawiki/index.php/Main_Page>


#### How to generate documentation

In a shell, in the main folder of LIA_RAL:

	doxygen doxygen.cfg

Of course, you need to have [`doxygen`](http://www.doxygen.org) installed.


#### How to compile LIA_RAL with multithreading enabled

Run `./configure` with the option `--enable-MT`.
 

#### How to link with Lapack

Install `LapackE`  then run `./configure` with the option `--enable-LAPACK`.


#### Linking with the SPro library for the speaker detection server

LIA_RAL includes a server for remote speaker detection (in `LIA_SpkDet/RemoteSpkDet`). The server can be sent features, in order to train models and run tests. It can also optionally be sent audio, which it will parameterize. In order to do this, the server must be compiled with a link to the SPro 4 library.

After downloading SPro 4.0.1 (<http://www.irisa.fr/metiss/guig/spro/download.html>), you must apply the patch that is found in the file `LIA_SpkDet/RemoteSpkDet/spro4_fft.c.patch`:

    patch path/to/spro/fft.c path/to/LIA_RAL/LIA_SpkDet/RemoteSpkDet/spro4_fft.c.patch

Then compile SPro according to the instructions in the package. If you are only going to use SPro for this purpose, there is no need to `make install` at the end, the compiled version can just stay in the distribution directory. If you install SPro somewhere else (by default, it installs in `/usr/local/`), then you need to copy the `system.h` file from the SPro distribution directory into the SPro `include` installation directory, since this file is needed for the compilation of `SpkDetServer` but not copied there by the SPro installation process.

Finally, use the `--with-spro` option when running `configure` for LIA_RAL. By default, it will look for SPro in `../spro4`. You can specify a different path if required (`--with-spro=path`).

Along with the speaker detection server, an example of a client is included. It is compiled with the rest of LIA_RAL, but it does not link with the `ALIZE` and `LIA_RAL`libraries and can be compiled completely independently.


------------------------------------------------------------

For more information please refer to:
<http://alize.univ-avignon.fr/mediawiki/index.php/Main_Page>

