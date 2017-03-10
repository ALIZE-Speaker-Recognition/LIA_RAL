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
2. Run aclocal, automake and autoconf.
3. Then run `./configure`.
   By default, the `ALIZE` library is searched for in `../alize-core`. It may be in a folder with a different name, depending on how you downloaded the library, or if you have decided to install it in a non-default location. If so, you can specify the absolute path by using the `--with-alize=ABSOUTE_PATH` option.
4. Finally, run `make`.


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


------------------------------------------------------------

For more information please refer to:
<http://alize.univ-avignon.fr/mediawiki/index.php/Main_Page>

