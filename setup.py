from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

extensions=[
			Extension('multidim_integration.symbolic_definite_integration',
              sources = ['symbolic_definite_integration.pyx']),
            Extension('multidim_integration.numerical_integration',
              sources = ['numerical_integration.pyx'])
]

ext_modules = cythonize(
               extensions, 
               compiler_directives={'language_level' : "2"}   # or "2" or "3str"
             ) 
setup(name='multidim_integration',cmdclass={'build_ext':build_ext},ext_modules=ext_modules)
