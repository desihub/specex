import os
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
            
from desiutil.setup import DesiTest, DesiVersion, get_version
#
# Begin setup
#
pname = 'specex'
pkg_version = get_version(pname)

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
        'version': DesiVersion
    },
    package_data={
        pname: ["data/*"]
    }
)
