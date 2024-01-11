from distutils.core import setup
from distutils.extension import Extension


#
#sourcefiles = ['SmoothBackfittingPy.pyx']

#extensions = [Extension(name="SmoothBackfittingPy", sources=sourcefiles, extra_compile_args=["-std=c++20"],extra_link_args=["-std=c++20"])]

#setup(
#    ext_modules=cythonize(extensions,compiler_directives={'language_level' : "3"}),
#    include_dirs = ['./include', './extern']    
#)



USE_CYTHON = True

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



#setup(
#    name='TestName',
#    ext_modules=cythonize("sumCython.pyx"),
#    include_dirs=[numpy.get_include()]
#)