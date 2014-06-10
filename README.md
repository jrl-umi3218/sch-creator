sch-creator
===========

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

FIXME: write this.

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
