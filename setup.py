from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name="cysolar",
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("cysolar", ["cysolar.pyx"], 
                             include_dirs=["."],
                            depends=["solar_constants.h"])]
)

