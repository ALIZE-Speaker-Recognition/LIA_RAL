<img src="http://alize.univ-avignon.fr/images/LIA_RAL.png" alt="The LIA_RAL logo" height="198" >

# LIA_RAL

*This package is part of the ALIZÉ project: <http://alize.univ-avignon.fr>*



Welcome to LIA_RAL!
-------------------

LIA_RAL is a set of utilities for Automatic Speaker Recognition developed at LIA, providing a high-level access to the [ALIZÉ platform](http://alize.univ-avignon.fr). All these utilities are based on the core `ALIZE` library, which is required in order to compile LIA_RAL.

LIA_RAL is an open project. Feel free to contact us to propose to work on it! 

In order to use Support Vector Machines, the library `libsvm` is included in this package. For more information on `libsvm`, please refer to:
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
 

#### Feature extraction

LIA_RAL provides tools for normalization of input features, but does not handle feature extraction from the audio signal.
A separate toolkit is required for this task.

ALIZÉ/LIA_RAL is most commonly used in conjunction with the free speech signal processing toolkit SPro, developed by Guillaume Gravier at IRISA:
<https://gforge.inria.fr/projects/spro/>

⚠️**Warning:** Note that only the revisions 155 and up of SPro are fully compatible with 64-bit CPUs. However, at the time of this writing, these versions of SPro are only available through Subversion, and the direct download link given on the website above points to an older revision of SPro 5 which includes a bug leading to corrupted feature files when compiled for 64 bit systems.
If you want to be sure to get the right version of SPro for use with ALIZÉ, you can download it from ALIZÉ's website: <http://alize.univ-avignon.fr/spro-5.0-157.tar.gz>.


#### Linking with the SPro library for SimpleSpkDetSystem

LIA_RAL includes a class (`LIA_SpkDet/SimpleSpkDetSystem`) which exposes a high-level API for a speaker verification/identification system, completely masking the details of how to implement such a system with ALIZÉ. This API targets application developers who want an easy way to embed speaker recognition in their applications, without having to learn all the details of speaker recognition.

The system can be fed features, in order to train models and run tests. But in its intended use as part of an application, it will be directly given an audio signal, which it will parameterize. In order to do this, the system must be compiled with a link to the SPro library, which will handle parameterization.

After downloading SPro 5 (please read the note above, under [Feature extraction]), compile it according to the instructions in the package.

If you are only going to use SPro for the purpose of linking LIA_RAL with its library, there is no need to `make install` after compiling it. The compiled version can just stay in the distribution directory. If you install SPro somewhere else (by default, it installs in `/usr/local/`), then you need to copy the `system.h` file from the SPro distribution directory into the SPro `include` installation directory, since this file is needed for the compilation of `SpkDetServer` but not copied there by the SPro installation process.

Finally, use the `--with-spro` option when running `configure` for LIA_RAL. By default, it will look for SPro in `../spro`. You can specify a different path if required (`--with-spro=path`).

The high-level speaker recognition API can be used with a local instance, in C++ (using `LIA_SpkDet/SimpleSpkDetSystem`) or in Java for Android applications (using the code found in the `Android-ALIZÉ`repository, which adds a JNI/Java layer over this class).

It is also usable over the network in a client/server mode (using the code in `LIA_SpkDet/RemoteSpkDet`). Along with a speaker detection server, an example of a client software is included. The client is compiled with the rest of LIA_RAL, but it does not link with the `ALIZE` and `LIA_RAL`libraries and can be compiled completely independently.


------------------------------------------------------------

For more information please refer to:
<http://alize.univ-avignon.fr/mediawiki/index.php/Main_Page>

