# Python interface

## Building the interface
To build the python interface along with the library and the rest of the executables,
make sure to execute the cmake configuration with the BUILD_PYTHON_INTERFACE variable set to ON.

`cmake -DBUILD_PYTHON_INTERFACE=ON`

This will generate the `_pysdg.so` and `pysdg.py` files under the `pysdg` directory linked to the `libsdg` c++ library.

## Using the python interface

To use the interface in a python interpreter, make sure to make the `_pysdg.so, libsdg and pysdg.py`  files
available to the interpreter either by `LD_LIBRARY_PATH` or `DYLIB_LIBRARY_PATH` (depending on your OS) or
by copying the files to the working directory of the interpreter.

Another option is to append the location of the compiled `bsg` objects to the python interpreter path using, e.g:

```python
import sys
sys.path.append("~/git_sources/sdg/build/pysdg")
```

The library should be available and can be loaded using, this command prints the library version:

```python
import pysdg as sdg
```

Documentation for the functions is available to the user by doing:

```python
sdg.SequenceDistanceGraph?
```

as done usually within a python environment.

### Notes
Make sure the python interpreter available whilst compiling the interface is the same
as the one you are trying to use the library with, having different interpreters might
mean that the library won't be loaded in the best case or will crash the 
interpreter in the worst, this is particularly true for OSX environments using brew
which can lead to having multiple python versions coexisting.

In order to redirect the library output to the interpreter, make sure you have the 
`wurlitzer` package installed. To install it, use `pip install wurlitzer`