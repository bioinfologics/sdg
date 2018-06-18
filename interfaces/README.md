# Python interface

## Building the interface
To build the python interface along with the library and the rest of the executables,
make sure to execute the cmake configuration with the BUILD_PYTHON_INTERFACE variable set to ON.

`cmake -DBUILD_PYTHON_INTERFACE=ON`

This will generate the `_bsg.so` and `bsg.py` files linked to the `libbsg` c++ library.

## Using the python interface

To use the interface in a python interpreter, make sure to make the `_bsg.so, libbsg and bsg.py`  files
available to the interpreter either by `LD_LIBRARY_PATH` or `DYLIB_LIBRARY_PATH` (depending on your OS) or
by copying the files to the working directory of the interpreter.

The library should be available and can be loaded using:

```python
import bsg
```

Documentation for the functions is available to the user by doing:

```python
bsg.SequenceGraph?
```

as done usually within a python environment.

### Notes
Make sure the python interpreter available whilst compiling the interface is the same
as the one you are trying to use the library with, having different interpreters might
mean that the library won't be loaded in the best case or will crash the 
interpreter in the worst, this is particularly true for OSX environments using brew
which can lead to having multiple python versions coexisting.

