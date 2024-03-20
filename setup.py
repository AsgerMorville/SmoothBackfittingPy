from distutils.core import setup
from distutils.extension import Extension

USE_CYTHON = False

if USE_CYTHON:
    from Cython.Build import cythonize
    sourcefiles = ['SmoothBackfittingPy.pyx']
    extensions = [Extension(name="SmoothBackfittingPy", sources=sourcefiles, extra_compile_args=["-std=c++20"],extra_link_args=["-std=c++20"])]
    setup(
        ext_modules=cythonize(extensions,compiler_directives={'language_level' : "3"}),
        include_dirs = ['SmoothBackfittingCpp/include', 'SmoothBackfittingCpp/extern']
        )
else:
    sourcefiles = ['SmoothBackfittingPy.cpp']
    setup(
        ext_modules = [Extension(name="SmoothBackfittingPy", sources=sourcefiles)],
        include_dirs = ['SmoothBackfittingCpp/include', 'SmoothBackfittingCpp/extern']
        )
