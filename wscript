# -*- mode: python; -*- 


import os
import sys
import commands
import string
import os.path as op

import Options
import Configure
import wafutils 


APPNAME  = 'specex'
VERSION  = '0.0.0'
top   = '.'
out   = wafutils.get_out_name()
DESCRIPTION = "spectral extraction"
requirements = [("harp","0.0",True)]
debug = True
optimize = 2
openmp = False

def harp_pkgconfig(pkgpath):

    if not os.path.isdir(pkgpath) :
        print "ERROR in harp_pkgconfig: cannot access %s"%pkgpath
        sys.exit(12)
    
    name="harp"
    description="description of harp"
    version="0.0"
    pkgconfig_file="%s/harp-%s.pc"%(pkgpath,version)
    
    if os.path.isfile(pkgconfig_file) :
        backup_file=pkgconfig_file+".back"
        print "WARNING in harp_pkgconfig: file %s exists"%pkgconfig_file
        print "                           save it in file %s"%backup_file
        if os.path.isfile(backup_file) :
            os.unlink(backup_file)
        os.link(pkgconfig_file,backup_file)
        
    harpconfpath=os.popen("which %s" % "harpconfig").read().strip() 
    if len(harpconfpath)==0 :
        print "ERROR in harp_pkgconfig: need 'harpconfig' in path"
        sys.exit(12)

    libdir=os.popen("harpconfig  --ldflags").read().strip().replace("-L","")
    prefix=libdir.replace("/lib","")
    libs=os.popen("harpconfig  --link").read().strip()
    cppflags=os.popen("harpconfig --cppflags").read().strip()
    includedir=string.split(str(cppflags)," ")[0].replace("-I","")

    file=open(pkgconfig_file,"w")
    file.write("prefix=%s\n"%prefix)
    file.write("libdir=%s\n"%libdir)
    file.write("includedir=%s\n"%includedir)
    file.write("\n")
    file.write("Name: %s\n"%name)
    file.write("Description: %s\n"%description)
    file.write("Version: %s\n"%version)
    file.write("Libs: %s\n"%libs)
    file.write("Requires:\n")
    file.write("Cflags: %s\n"%cppflags)
    file.close()
    print "wrote",pkgconfig_file


def options(opt):    
    opt.load('compiler_cc')
    opt.load('compiler_cxx')
    #opt.load('compiler_fc')
    opt.add_option('--debug', 
                   action='store', 
                   default=True, 
                   dest='debug', 
                   help='True or False, turn on/off the -g / -DNDEBUG option, (True by default)')
    opt.add_option('--optimize', 
                   action='store', 
                   default=optimize, 
                   dest='optimize', 
                   help='0, 1, 2 or 3, default is %d, i.e. -O%d'%(optimize,optimize)) 
    opt.add_option('--enable-openmp', 
                   action='store', 
                   default=openmp, 
                   dest='openmp', 
                   help='enable OpenMP support [default=%s]'%str(openmp)) 
    opt.add_option("--static",
                   action="store_true",
                   dest='static',
                   help="build static library")
    opt.add_option('--autogen-harp-pkgconfig',
                   action='store', 
                   default=False, 
                   dest='harp_pkgconfig_dir', 
                   help='auto generate a pkgconfig file for HARP in this path')

def configure(conf):    
    
    # c compiler 
    conf.load( 'compiler_c' )
    conf.env['CCFLAGS'] = ['-fPIC', '-DPIC']
    
    # c++ compiler 
    conf.load( 'compiler_cxx' )
    conf.env['CXXFLAGS'] = ['-fPIC', '-DPIC','-Wuninitialized','-Wunused-value','-Wunused-variable']
# option -Wmaybe-uninitialized not recognized by all compilers, I don't want Wall because very many unused typdef in boost

    conf.env['PKG_INCDIR'] = op.join('include', '%s-%s' % (APPNAME,VERSION))

    print "optimize=",conf.options.optimize
    opti=string.atoi(str(conf.options.optimize))
    if opti > 0 :
        conf.env['CFLAGS'].append('-O%d' % opti )
        conf.env['CXXFLAGS'].append('-O%d' % opti )

    if conf.options.openmp :
        conf.env['CCFLAGS'].append('-fopenmp')
        conf.env['CXXFLAGS'].append('-fopenmp')
        conf.env['LINKFLAGS'].append('-fopenmp')
        conf.env['LINKFLAGS_HARP'].append('-fopenmp')

    print "debug=",conf.options.debug
    if conf.options.debug == "True" :
        conf.env['CFLAGS'].append('-g')
        conf.env['CXXFLAGS'].append('-g')
    else :
        conf.env['CFLAGS'].append('-DNDEBUG')
        conf.env['CXXFLAGS'].append('-DNDEBUG')
        conf.env['CFLAGS'].append('-DBOOST_UBLAS_NDEBUG')
        conf.env['CXXFLAGS'].append('-DBOOST_UBLAS_NDEBUG')
    if conf.options.harp_pkgconfig_dir :
        harp_pkgconfig(conf.options.harp_pkgconfig_dir)
    
    if conf.options.static:
        conf.env['static'] = True
    else:
        conv.env['static'] = False

    #conf.check_cc(lib='z', msg='Checking for zlib')    
    conf.check_packages(requirements)    

    conf.write_config_header("config.h")

    print conf.env

def load_pkg_config(conf):
    """
    Load pkg-config (used to read the package metadata)
    """
    
    # first, we need pkg-config 
    if not conf.env['PKG_CONFIG'] or \
       not conf.env['PKG_CONFIG_PATH']:
        try:
            conf.find_program('pkg-config', var='PKG_CONFIG')
            pkgcpath = op.join( conf.env.PREFIX, 'lib', 'pkgconfig')
            if os.environ.has_key('PKG_CONFIG_PATH'):
                os.environ['PKG_CONFIG_PATH'] = pkgcpath + ":" + os.environ['PKG_CONFIG_PATH']
            else:
                os.environ['PKG_CONFIG_PATH'] = pkgcpath
        except Configure.ConfigurationError:
            conf.fatal('pkg-config not found')
    

@Configure.conftest
def check_packages(conf, pkg_list):
    """
    - Check whether the packages passed in argument 
      can be found by pkg-config. 
    - Parse the cflags and libs and store the information 
      in the uselib_store.
    """
    load_pkg_config(conf)
    # check for the packages passed in arguments 
    for pkg in pkg_list:
        try:
            pkg_name, pkg_version, pkg_mandatory = pkg
            conf.check_cfg(args='--cflags --libs', 
                           package = pkg_name + '-' + pkg_version, 
                           mandatory=pkg_mandatory, 
                           uselib_store=pkg_name.replace('-', '_').upper())
        except Configure.ConfigurationError:
            conf.fatal('unable to locate %s-%s (mandatory)' % (pkg_name, pkg_version))


def build(bld):
    #bld.add_subdirs( ['src/library', 'src/apps', 'src/tests'])
    bld.add_subdirs( ['src/library', 'src/plugin', 'src/apps', 'src/tests'])
    gen_pkgconfig(bld)


def install_headers(bld, headers):
    """
    Just hide the install header commands 
    """
    install_dir = op.join('$PREFIX', 'include', '%s-%s' % (APPNAME,VERSION))
    bld.install_files(install_dir, headers)


def gen_pkgconfig(bld):

    # solution directly from T. Nagy (see email in google-groups)
    from waflib.TaskGen import feature, before, after
    @feature('subst')
    @before('process_subst')
    def read_libs(self):
        self.env.append_value('ALL_LIBS', " ")
        for g in self.bld.groups:
            for tg in g:
                try:
                    if 'cshlib' in tg.features or \
                       'cxxshlib' in tg.features or \
                       'fcshlib' in tg.features:
                        # or 'cxxshlib'/'fcshlib'/'dshlib' in tg.features
                        self.env.append_value('ALL_LIBS', "-l" + tg.name)                    
                except:
                    pass
    
    requirements = " "
        
    obj = bld(features = 'subst', 
              target = '%s-%s.pc' % (APPNAME, VERSION),
              source = 'pc.in', 
              install_path =  '${PREFIX}/lib/pkgconfig/',
              PREFIX = bld.env['PREFIX'], 
              APPNAME   = APPNAME,
              DESCRIPTION = DESCRIPTION,
              VERSION = VERSION,
              REQUIREMENTS = requirements)

