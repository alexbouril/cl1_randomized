# python setup.py build_ext --inplace
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
package = "CL1R"
file = "fib"
ext_modules=[
    Extension("%s"%(file),    # location of the resulting .so
             ["%s.pyx"%(file)],) ]

setup(name='package',
      packages=find_packages(),
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules,
     )