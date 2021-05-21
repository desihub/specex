#!/bin/bash
build_dev=1
compare=1
dirac=1
tag=

do_dev=1
do_master=0

pull_master=1
build_master=0

loc_dev=1

mstcode='mst'
devcode='dev'

mstbranch='master'
devbranch='rm_harpio'

if [ $dirac -gt 0 ] ; then
    devtestdir=/global/home/users/malvarez/dev/dev-test
    inputroot=$devtestdir
else
    devtestdir=/global/common/software/desi/users/malvarez/dev-test
    inputroot=/global/cfs/cdirs/desi/spectro/redux/blanc/preproc/20201216/00068217/
fi
cd $devtestdir

# activate desi environment
if [ $dirac -gt 0 ] ; then
    source dirac_setup.sh
    export CXX=`which c++` # on dirac force C++ compiler to version after dirac_setup.sh
else
    source /global/common/software/desi/desi_environment.sh
    source ${DESICONDA}/etc/profile.d/conda.sh
    conda activate
    export SPECEXDATA=/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/code/specex/0.6.7/data
    # clear default system specex and desispec
    module unload specex
    module unload desispec
fi

specex_system=/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/code/specex/master/lib/python3.8/site-packages/specex-0.6.7.dev617-py3.8-linux-x86_64.egg

desi_compute_psf_args="--input-image /global/cfs/cdirs/desi/spectro/redux/blanc/preproc/20201216/00068217/preproc-b1-00068217.fits --input-psf /global/cfs/cdirs/desi/spectro/redux/blanc/exposures/20201216/00068217/shifted-input-psf-b1-00068217.fits --debug" 

if [ $do_master -gt 0 ]; then
    # clean
    rm -f *mst*.fits mst_version.log
    
    # change current desispec to master
    if [ $dirac -eq 0 ] ; then
	module load desispec
    fi
    
    # change specex to master
    if [ $pull_master -eq 0 ]; then
	specex_loc=$specex_system
    else
	specex_loc=$PWD/specex-$mstcode/code/lib/python3.8/site-packages/specex-0.7.0.dev637-py3.8-linux-x86_64.egg
    fi
    export PYTHONPATH=$specex_loc:$PYTHONPATH
    if [ $build_master -gt 0 ]; then
	rm -rf specex-$mstcode    
	git clone --single-branch --branch $mstbranch \
	    https://github.com/desihub/specex \
	    specex-$mstcode
	cd specex-$mstcode
	rm -rf CMakeF* CMakeC* build code
	python setup.py install --prefix ./code -v 
	cd ../
    fi
    # do desispec-mst and specex-mst
    code=mst
    echo $specex_loc >> mst_version.log    
    echo `ls -rlt $specex_loc` >> mst_version.log
    echo `which desi_compute_psf` >> mst_version.log
    (time python test_specex.py psf-$code.fits $inputroot) |& tee psf-b1-$code.log
    
fi

if [ $do_dev -gt 0 ]; then
    # clean
    rm -f *dev*.fits dev_version.log

    if [ $build_dev -gt 0 ]; then
	rm -rf specex-$devcode    
	if [ $loc_dev -gt 0 ]; then
	    # get local branch
	    cp -r ../specex-$devcode ./	
	else
	    # get current specex dev branch
	    git clone --single-branch --branch $devbranch \
		https://github.com/desihub/specex \
		specex-$devcode
	fi    

	# compile specex dev 
	cd specex-$devcode; mkdir code
	rm -rf CMakeF* CMakeC* build code
	python setup.py install --prefix ./code -v 
	cd ../

    fi
    
    # change current desispec to master
    if [ $dirac -eq 0 ] ; then
	module load desispec
    fi
    
    # change specex to dev branch
    specex_loc=$PWD/specex-$devcode/code/lib/python3.8/site-packages/specex-0.7.0.dev637-py3.8-linux-x86_64.egg
    export PYTHONPATH=$specex_loc:$PYTHONPATH
    
    # do specex-dev
    #./linkpython.sh
    echo $specex_loc >> dev_version.log    
    echo `ls -rlt $specex_loc` >> dev_version.log
    echo `which desi_compute_psf` >> dev_version.log
    code=dev
    (time python test_specex.py psf-$code.fits $inputroot) |& tee psf-b1-$code\_$tag.log
    
fi

if [ $compare -gt 0 ]; then
    diff psf-mst.fits psf-dev.fits
    python compare_psf.py psf-mst.fits psf-dev.fits
    grep user *log
fi
