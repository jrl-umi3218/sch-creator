sch-creator
===========

[![Build Status](https://travis-ci.org/jrl-umi3218/sch-creator.svg?branch=master)](https://travis-ci.org/jrl-umi3218/sch-creator)
[![Coverage Status](https://coveralls.io/repos/jrl-umi3218/sch-creator/badge.png)](https://coveralls.io/r/jrl-umi3218/sch-creator)

Documentation
-------------

This package contains scripts and methods allowing to build convex surfaces 
from vrml files or blender objects.
It also convert those files into qhull volumes.

Execution
---------

Several scripts are installed in the `${CMAKE_INSTALL_PREFIX}/bin` directory
Those scripts allow to convert vrml files into cloud files (simplified geometry
files).

### Building from vrml files

FIXME: write this.

### Building from blender files:

- install this package (see instructions below)
- open blender
- click on «File/User Preferences...» menu
- click on «Addons» tab
- click «Import-Export» sub section
- enable «Import-Export: qconvex cloud format»
- click on «Save As Default» button
- close blend

### Converting

- run `script/blender_2cloud.sh output_directory file1 file2 ... fileN`
- output_directory must be fill with .qc files
- `cloud2qhull.sh qc_directory output_directory` convert .qc files into convex hull
- `cloud2sch qc_directory output_directory [R,r]` convert .qc files into sch hull

Dependency
----------

* sch-core
* blender (optional)
* qhull (optional)

Installation
------------

Building using cmake:
- create a new build directory
- run `cmake ..` in this directory 
  (you probably want to define `CMAKE_INSTALL_PREFIX`)
- `make`
- `make install`

Indentation
-----------

    astyle --style=allman --lineend=linux --indent=spaces=2
