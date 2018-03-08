from distutils.core import setup, Extension
import numpy

module = Extension('rkernel',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = [numpy.get_include(), 'python/', 'lib/'],
                    # language = "c++",
                    libraries = [],
                    library_dirs = [],
                    extra_compile_args=['-O0'],
                    sources = ['python/rkernel.cpp'])

setup (name = 'RKERNEL',
       version = '1.0',
       description = 'Fast kernels for kernel methods',
       author = 'Remi Lespinet',
       author_email = 'remi.lespinet@gmail.com',
       url = '',
       long_description = '''
Native implementation of some kernels for the kernel method course (challenge)
''',
       ext_modules = [module])
