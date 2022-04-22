# Overview

This is a C++ library for computing obstructions to symplectic embeddings using a certain combinatorial criterion developed by Michael Hutchings in his paper ["Beyond ECH Capacities"](https://arxiv.org/abs/1409.1352). These obstructions are derived from embedded contact homology (ECH) and are stronger than the obstructions obtained from ECH capacities. There are several known cases where symplectic folding can be used to show that the obstructions produced by Hutchings' criterion are optimal.

I wrote this code in 2016 as part of a Research Experience for Undergraduates (REU) program at Columbia University. This research ultimately resulted in a [paper](https://arxiv.org/abs/1610.00566) (co-authored with Jo Nelson, and published under my former legal name) that extends some of Hutchings' results in "Beyond ECH Capacities."
[Note: this research was partially supported by NSF grants DMS-1206667, DMS-0970108, and a Graduate Research Fellowship.]

This version of the library is written in C++. The original version was written in Haskell; however, some of the combinatorics involved are very computationally intensive, so I rewrote it in C++, which resulted in a huge performance improvement. If you're interested in doing computations for research with this library, I recommend this version. However, if you just want to test out some of these combinatorial computations and see how they work, the Haskell version should work just fine. If you're interested, you can found the Haskell version [here](https://github.com/michristianson/ech-capacities-haskell).



# Comparison with the Haskell Version
As mentioned above, this project is a rewriting in C++ of a Haskell library called "CIP." This rewriting made the code orders of magnitude faster. As an informal benchmark, I tried running the "generate less thans" function (which is the most computationally intensive part of the code) in both programs with the same parameters:

- ctd1 = Polydisk(2.43, 1)
- ctd2 = Ellipsoid(3.205, 3.205)
- cg = 30e(1,-1)

In words, this means computing all convex generators less than 30e(1,-1) with respect to Polydisk(2.43, 1) and Ellipsoid(3.205, 3.205). Both programs were compiled with maximal optimization (for GHC, the -O2 flag, and for g++, the -O3 flag), and each was run several times each to account for differences in cache, background processes, etc. The results were:

Haskell version: fluctuated between about 14.5 and 46.5 seconds.
C++ version: fluctuated between 0.13 and 0.43 seconds.

This is of course a very informal test, but I think it makes it clear that the C++ version is significantly faster. (Note: the Haskell version, fortunately, is able to use thm114Heuristics here. Without using any heuristics, the runtime of the Haskell version would probably be at least 10x longer. The C++ version, by contrast, applies similar heuristics automatically whenever they apply.)

In addition to the huge speed increase, the C++ version is generally more polished than the Haskell version. This is probably in part because I'm more experienced with object-oriented programming than with functional programming and in part because it's always easier to write things well the second time around. A couple examples of such improvements are:

- The Haskell version uses as its main "less than" heuristic an upper-bound on the y-values of potential less thans. The C++ version uses this and an upper-bound on x-values as well.
- The Haskell version requires explicitly specifying the y-value upper-bound; if this is not done, then no heuristic will be used. The C++ version computes its upper-bound heuristics automatically, using polymorphism of the relevant objects, and a (relatively weak but still useful) default bound is applied even in the most general cases.
- The Haskell version originally memoized "less than" computations using a relational database with an ODBC driver. (This memoization has been removed from the current version because of portability concerns, but it could be added back in; see the readme of that version for details.) By contrast, the C++ version allows memoization using flat files, which are simple and portable. 

In spite of these benefits, there is one major disadvantage to CTD (the C++ version) over CIP (the Haskell version):  the functions for actually obstructing symplectic embeddings (i.e. everything in the Obstruct module of CIP) was never implemented in CTD. This is because the data obtained from the "less than" calculations (which could now be done for much larger convex generators than before, due to the speed improvements) were already enough to show me very interesting numerical patterns. I ended up spending the rest of my time on this project doing mathematical research guided by just that data, so I never needed the full obstruction functions. These would not be particularly difficult to implement; the bulk of the work, by far, is the "generate less thans" procedure, and that is already done. I would love to revist this project some day, finish the obstruction functions, and maybe make a more full UI; unfortunately, I can't make any promises about if or when that will happen. If you're interested in doing this, please send me a pull request!



# Installation

Since this code was intended for personal use, it was never packaged for distribution in any nice way. However, it is a relatively small and portable library, so it isn't too hard to get it to work.

## Dependencies

The only dependence for this library is boost. This is unfortunately a very large dependency for a very small library; on the bright side, it is very common and readily available. See [this page](https://www.boost.org/doc/libs/1_79_0/more/getting_started/index.html) for instructions on how to install boost.

On Linux, boost will most likely work out-of-the-box if you just install the right package from the main repository of your package manager (apt, pacman, yum, rpm, etc.) If you choose to install it this way, you can look for a package named something like "libboost-all" to install all of boost, or you can look for piecemeal packages like "libboost-math" or "libboost-filesystem" to install only the boost libraries needed for this library. The specific boost libraries that this project depends on are:

- bind
- filesystem 
- integer 
- spirit

## Compiling from Source

Once boost is installed, you just need to compile the code in this repository using the compiler of your choice. A makefile is provided in this repository, so if you have gcc (as many Linux systems do by default), you can just use the following commands on the command line:

1. cd local/directory/for/this/library
2. git clone https://github.com/michristianson/ech-capacities-cpp
3. make

On Windows, the easiest solution may be to download (or git clone) this repository,
then import all the files into Visual Studio and compile from there.

However you choose to compile the source code, here are a couple important details to keep in mind. (These are already taken care of by the provided makefile, so if you're using make, you don't need to worry about them.)

- Boost's filesystem library is one of the few boost libraries that is not header-only, so it must be built separately. Boost's installation instructions describe how to build this library (specifically, see [this section](https://www.boost.org/doc/libs/1_79_0/more/getting_started/unix-variants.html#prepare-to-use-a-boost-library-binary) of those instructions). If the library has already been built during installation (for instance, this was true for me after installing boost using my native package manager on Linux), then you just need to make sure to link it during compilation. For instance, on gcc/g++, this means using the flag -lboost_filesystem.
- Since this is a relatively computationally intensive library, you should probably turn on compiler optimzations. For gcc/g++, this means -O3 (or at least -O2; but in my tests, compilation was well under a minute in both cases, and -O3 did provide a noticeable speedup).



# Overview of Project Files
- ConvexToricDomains.cpp: contains main(), a few useful I/O functions, and a couple of functions that can be used for basic computations. More specifically, the following computations are provided:
  - testLTs(): runs "less than" computations for a large range of different inputs to test functionality and prints them to stdout. This function generates quite a lot of output (over 39,000 lines of text, which came out to a 2.4MB text file on my system), so it's probably best to pipe output to a file if you want to run this. 
  - computeGivenLTs(): asks for parameters for a polydisk, a ball, and a convex generator, then computes all less thans for the convex generator with respect to the polydisk and the ball. Results are printed to stdout.
  Currently, main() is set to just run computeGivenLTs() once.

- data.hpp and data.cpp: the equivalent of Internals.hs from the Haskell version. These files contain definitions of the main data types as well as some basic operations on them. The main types of interest here are:
  - CGen: a class which represents a convex generator.
  - CTD: a class which represents a convex toric domain.
  Both classes have various methods for computing basic values (e.g. cgen.I() computes the index I(cgen), cgen.x() computes x(cgen), ctd.actionOf(cgen) computes the action of cgen with respect to CTD, etc.)

- lessthans.hpp and lessthans.cpp: the equivalent of LessThans.hs from the Haskell version. These files contain the most computationally difficult work of the application, which is determining the CGen's that are less than or equal to a given CGen. For this, a few important classes are supplied:
  - LTGen: uses a generator pattern to produce "less than" results one-by-one. Users are intended to use the functions start(), atEnd(), next(),
    and end() to iterate through all "less than" results. This is a template
    class with one template parameter, which specifies a Memoizer policy.
  - NoMemoizer and FSMemoizer: implementations of the Memoizer policy for LTGen. NoMemoizer causes no memoization to be done, while FSMemoizer
    memoizes using flat files in the directory

    current_working_directory/memos/FSMemoizer/.

    Memos will be arranged in subdirectories of the above directory based on the parameters of the "less than" search. (Note: the FS in the name
    FSMemoizer stands for "filesystem.")

- obstruct.hpp and obstruct.cpp: the equivalent of Obstruct.hs from the Haskell version. Intended to be used for implementing functions to obstruct symplectic embeddings. As of yet, however, these have not been implemented.



# Documentation

There is no formal documentation of this library. However, I've left very extensive comments in the code, and especially in the header files; I think these are probably more thorough than the documentation of the average project. If you're looking for documentation, I recommend first reading the "Overview of Project Files" section above, then reading through the header files data.hpp and lessthans.hpp. Looking at the code in ConvexToricDomains.cpp (especially the code in computeGivenLTs()) may also be a good idea, as this will give an example of how to use this library for "less than" computations.


# Licensing

This library is provided under the MIT license. See LICENSE.txt for details. The gist of it is: feel free to do anything you want with this code, provided you retain the information in LICENSE.txt whenever you make copies of it.
