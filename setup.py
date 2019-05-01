import platform
from setuptools import setup, Extension
# from distutils.core import setup
# from distutils.extension import Extension
# from Cython.Build import cythonize

compile_extra_args = []
link_extra_args = []

if platform.system() == "Windows":
    compile_extra_args = ["/std:c++latest", "/EHsc"]
elif platform.system() == "Darwin":
    compile_extra_args = ['-std=c++11', "-mmacosx-version-min=10.9"]
    link_extra_args = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

with open("README.md", "r") as fh:
	long_description = fh.read()

setup(
	name = 'matsim',
	version = '0.2.3',
	author = "Tomas Barta",
	author_email = "tomas.barta@fgu.cas.cz",
	description = "Package primarily for simulation of MAT model",
	long_description = long_description,
	# ext_modules = [
	# 	Extension('matsim'),
	# 	sources=['matsim.pyx'],
	# 	extra_compile_args=compile_extra_args,
	# 	extra_link_args
	# ],
	# ext_modules=cythonize("matsim/matsim.pyx"),
	ext_modules = [Extension(
		"matsim",
		sources=["matsim/matsim.cpp"],
	    extra_compile_args = compile_extra_args,
        extra_link_args = link_extra_args
	)],
	package_data = {'matsim': ["matsim/bio.h", "matsim/simulation.h",
	    "matsim/bio.cpp","matsim/simulation.cpp"]}
)