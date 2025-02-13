import os
import sys
import pathlib

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig

class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])

ename='specex._libspecex'
class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.parent.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config
        ]
        # example of build args
        build_args = [
            '--config', config,
            '--', '-j8'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(str(cwd))

        for ext in self.extensions:
            dest_path = pathlib.Path(self.get_ext_fullpath(ext.name)).resolve()
            print('dest_path:',dest_path)

from desiutil.setup import get_version
#
# Begin setup
#
pname = 'specex'
pkg_version = get_version(pname)


#
# Warning about old features.
#
VERSION_HELP = """
Note: Generating version strings is no longer done using 'python setup.py version'. Instead
you will need to run:

    desi_update_version [-t TAG] desiutil

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    pytest

If you don't already have pytest installed, you can install it with:

    pip install pytest
"""

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py {0}'. Instead you will need to run:

    sphinx-build -W --keep-going -b html doc doc/_build/html

If you don't already have Sphinx installed, you can install it with:

    pip install Sphinx
"""

message = {'test': TEST_HELP,
           'version': VERSION_HELP,
           'build_docs': DOCS_HELP.format('build_docs'),
           'build_sphinx': DOCS_HELP.format('build_sphinx'), }

for m in message:
    if m in sys.argv:
        print(message[m])
        sys.exit(1)

packages=find_packages('py')
setup(
    name=pname,
    description='DESI PSF Fitting',
    author='DESI Collaboration',
    author_email='desi-data@desi.lbl.gov',
    packages=packages,
    package_dir={'': 'py'},
    license='BSD',
    zip_safe=False,
    version=pkg_version,
    ext_modules=[CMakeExtension(ename)],
    cmdclass={
        'build_ext': build_ext,
    },
    package_data={
        pname: ["data/*"]
    }
)
