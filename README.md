# NGSlib: A C++ library for reading high-throughtput sequencing data base on htslib.

## Installation

You should install htslib first, before you compile ngslib.


### How to install htslib


1. Shift to htscodecs directory and run the following commands: 

```bash

$ cd htslib/htscodecs
$ autoreconf -i
$ ./configure
$ make

```



2. Go back to the upper directory and install main htslib by running the commands below:

```bash

$ cd htslib

$ autoreconf -i
$ ./configure
$ make

```




