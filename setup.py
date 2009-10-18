from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name="cysolar",
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("cysolar", ["cysolar.pyx"], 
                             include_dirs=["."],
                            depends=["solar_constants.h"])],
    license="BSD",

      keywords      = 'solar cython',
      author        = 'Brent Pedersen',
      author_email  = 'bpederse@gmail.com',
      url   = 'http://github.com/brentp/cysolar/',
      long_description = open('README.rst').read(),
      classifiers   = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        ],

)

