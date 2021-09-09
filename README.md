# NGSlib: A C++ library for reading high-throughtput sequencing data base on htslib.

## Installation

You should install htslib first, before you compile ngslib.


### How to install htslib

1. Download NGSlib from github:

```bash
$ git clone --recursive git@github.com:ShujiaHuang/ngslib.git 
```

> WARNING: Please try several times if fail to clone the data causing by 
> the network problem.


2. Shift to htscodecs directory and run the following commands: 

```bash

$ cd htslib/htscodecs
$ autoreconf -i
$ ./configure
$ make

```



3. Go back to the upper directory and install main htslib by running the commands below:

```bash

$ cd htslib

$ autoreconf -i
$ ./configure
$ make

```




