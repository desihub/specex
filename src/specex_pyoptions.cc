#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <specex_pyoptions.h>

using namespace std;

#ifdef GETOPT

void PrintHelp()
{
  return;
  std::cout <<
    "--help,-h             display usage information\n" 
    "--arc,-a              preprocessed fits image file name (mandatory)\n" 
    "--in-psf              input psf file name (fits or xml, can contain only the traces, mandatory)\n"  
    "--out-psf             output psf fits file name (mandatory)\n"
    "--lamp-lines          input lamp lines file name (mandatory)\n" 
    "--only-trace-fit      only fit the trace coordinates (default is trace,sigma,psf)\n"
    "--no-trace-fit        do not fit the trace coordinates (default is trace,sigma,psf)\n"
    "--no-sigma-fit        do not fit the gaussian sigma (default is trace,sigma,psf)\n"
    "--no-psf-fit          do not fit the psf (default is trace,sigma,psf)\n"
    "--flux-hdu            flux hdu in input arc fits, for unusual data format\n"
    "--ivar-hdu            ivar hdu in input arc fits, for unusual data format\n"
    "--mask-hdu            mask hdu in input arc fits, for unusual data format\n"
    "--header-hdu          header hdu in input arc fits, for unusual data format\n"
    "--xcoord-hdu          hdu of x trace legendre polynomial of wavelength (default is extension XTRACE)\n" 
    "--ycoord-hdu          hdu of y trace legendre polynomial of wavelength (default is extension YTRACE)\n" 
    "--first-bundle        first fiber bundle to fit\n"
    "--last-bundle         last fiber bundle to fit\n"
    "--single-bundle       fit a single bundle with all fibers\n"
    "--first-fiber         first fiber (must be in bundle)\n"
    "--last-fiber          last fiber (must be in bundle)\n"
    "--half-size-x         half size of PSF stamp (full size is 2*half_size+1)\n"
    "--half-size-y         half size of PSF stamp (full size is 2*half_size+1)\n"
    "--psfmodel            PSF model, default is GAUSSHERMITE\n"
    "--positions           fit positions of each spot individually after global fit for debugging\n"
    "--verbose,v           turn on verbose mode (deprecated, true by default)\n" 
    "--quiet               no info message, only warning\n" 
    "--debug               turn on debug mode\n" 
    "--trace-prior-deg     force equal trace coeff in bundle starting at this degree\n"
    "--lamplines           lamp lines ASCII file name (def. is $SPECEXDATA/specex_linelist_desi.txt)\n" 
    "--core                dump core files when an exception is thrown\n"
    "--gauss-hermite-deg   degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)\n"
    "--gauss-hermite-deg2  degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)\n"
    "--legendre-deg-wave   degree of Legendre polynomials along wavelength (can be reduced if missing data)\n"
    "--legendre-deg-x      degree of Legendre polynomials along x_ccd (can be reduced if missing data)\n"
    "--trace-deg-wave      degree of Legendre polynomials along wavelength for fit of traces\n"
    "--trace-deg-x         degree of Legendre polynomials along x_ccd for fit of traces\n"
    "--psf-error           psf fractional uncertainty (default is 0.01, for weights in the fit)\n"
    "--psf-core-wscale     scale up the weight of pixels in 5x5 PSF core\n"
    "--broken-fibers       broken fibers (comma separated list)\n"

    "--variance-model      refit at the end with a model of the variance to avoid Poisson noise bias\n"    
    "--out-psf-xml         output psf xml file name\n"
    "--out-spots           output spots file name\n"  
    "--prior               gaussian prior on a param : 'name' value error\n"  
    "--tmp_results         write tmp results\n"  
#ifdef EXTERNAL_TAIL
    "--fit-psf-tails       unable fit of psf tails\n"
#endif
#ifdef CONTINUUM
    "--fit-continuum       unable fit of continuum\n"
#endif
    "--nlines              max # emission lines used (uses an algorithm to select best ones \n"
    "                      based on S/N and line coverage\n";

  return;
  
}

void specex::PyOptions::loadmap(map<string,vector<int>>& optmap, const string argstring, int req){
  optmap.insert(pair<string,vector<int>>({argstring,vector<int>({req, noptions++})}));
  return;
}

int specex::PyOptions::argint(map<string,vector<int>>& optmap, const string argstring){
  return optmap.find(argstring)->second[1];
}

int specex::PyOptions::parse(int argc, char *argv[] )
{

  using namespace std;
  map<string, vector<int>> optmap;
  map<string, vector<int>>::iterator it;
  
  int noptsmax=500;
  option long_opts[noptsmax];

  noptions=0;
  loadmap(optmap, "help",               required_argument);
  loadmap(optmap, "arc",                required_argument);
  loadmap(optmap, "help",               required_argument);
  loadmap(optmap, "verbose",            optional_argument);
  loadmap(optmap, "in-psf",             required_argument);
  loadmap(optmap, "out-psf",            required_argument);
  loadmap(optmap, "lamp-lines",         required_argument);
  loadmap(optmap, "only-trace-fit",     optional_argument);
  loadmap(optmap, "no-trace-fit",       optional_argument);
  loadmap(optmap, "no-sigma-fit",       optional_argument);
  loadmap(optmap, "no-psf-fit",         optional_argument);
  loadmap(optmap, "flux-hdu",           required_argument);
  loadmap(optmap, "ivar-hdu",           required_argument);
  loadmap(optmap, "mask-hdu",           required_argument);
  loadmap(optmap, "header-hdu",         required_argument);
  loadmap(optmap, "xcoord-hdu",         required_argument);
  loadmap(optmap, "ycoord-hdu",         required_argument);
  loadmap(optmap, "first-bundle",       required_argument);
  loadmap(optmap, "last-bundle",        required_argument);
  loadmap(optmap, "single-bundle",      optional_argument);
  loadmap(optmap, "first-fiber",        required_argument);
  loadmap(optmap, "last-fiber",         required_argument);
  loadmap(optmap, "half-size-x",        required_argument);
  loadmap(optmap, "half-size-y",        required_argument);
  loadmap(optmap, "psfmodel",           required_argument);
  loadmap(optmap, "positions",          optional_argument);
  loadmap(optmap, "quiet",              optional_argument);
  loadmap(optmap, "debug",              optional_argument);
  loadmap(optmap, "trace-prior-deg",    required_argument);
  loadmap(optmap, "lamplines",          required_argument);
  loadmap(optmap, "core",               optional_argument);
  loadmap(optmap, "gauss-hermite-deg",  required_argument);
  loadmap(optmap, "gauss-hermite-deg2", required_argument);
  loadmap(optmap, "legendre-deg-wave",  required_argument);
  loadmap(optmap, "legendre-deg-x",     required_argument);
  loadmap(optmap, "trace-deg-wave",     required_argument);
  loadmap(optmap, "trace-deg-x",        required_argument);
  loadmap(optmap, "psf-error",          required_argument);
  loadmap(optmap, "psf-core-wscale",    required_argument);
  loadmap(optmap, "broken-fibers",      required_argument);
  loadmap(optmap, "variance-model",     optional_argument);
  loadmap(optmap, "out-psf-xml",        required_argument);
  loadmap(optmap, "out-spots",          required_argument);
  loadmap(optmap, "prior",              required_argument);
  loadmap(optmap, "tmp_results",        optional_argument);
#ifdef EXTERNAL_TAIL
  loadmap(optmap, "fit-psf-tails",      optional_argument);
#endif
#ifdef CONTINUUM
  loadmap(optmap, "fit-continuum",      optional_argument);
#endif
  loadmap(optmap, "nlines",             required_argument);
 
  /*
  const option long_opts[] = {
    {"help",               required_argument, nullptr, 'h'},
    {"arc",                required_argument, nullptr, 'arc'},
    {"verbose",            optional_argument, nullptr, 'ver'},
    {"in-psf",             required_argument, nullptr, 'ip'},
    {"out-psf",            required_argument, nullptr, 'op'},
    {"lamp-lines",         required_argument, nullptr, 'll'},
    {"only-trace-fit",     optional_argument, nullptr, 'otf'},
    {"no-trace-fit",       optional_argument, nullptr, 'ntf'},
    {"no-sigma-fit",       optional_argument, nullptr, 'nsf'},
    {"no-psf-fit",         optional_argument, nullptr, 'npf'},
    {"flux-hdu",           optional_argument, nullptr, 'fh'},
    {"ivar-hdu",           optional_argument, nullptr, 'ih'},
    {"mask-hdu",           optional_argument, nullptr, 'mh'},
    {"header-hdu",         optional_argument, nullptr, 'hh'},
    {"xcoord-hdu",         optional_argument, nullptr, 'xh'},
    {"ycoord-hdu ",        optional_argument, nullptr, 'yh'},
    {"first-bundle",       optional_argument, nullptr, 'fb'},
    {"last-bundle",        optional_argument, nullptr, 'lb'},
    {"single-bundle",      optional_argument, nullptr, 'sb'},
    {"first-fiber",        optional_argument, nullptr, 'ff'}, 
    {"last-fiber",         optional_argument, nullptr, 'lf'},
    {"half-size-x",        optional_argument, nullptr, 'hsx'},
    {"half-size-y",        optional_argument, nullptr, 'hsy'},
    {"psfmodel",           optional_argument, nullptr, 'pm'},
    {"positions",          optional_argument, nullptr, 'pos'},
    {"quiet",              optional_argument, nullptr, 'q'},
    {"debug",              optional_argument, nullptr, 'd'},
    {"trace-prior-deg",    optional_argument, nullptr, 'tpd'},
    {"lamplines",          optional_argument, nullptr, 'lls'},
    {"core",               optional_argument, nullptr, 'c'},
    {"gauss-hermite-deg",  optional_argument, nullptr, 'ghd'},
    {"gauss-hermite-deg2", optional_argument, nullptr, 'ghd2'},
    {"legendre-deg-wave",  optional_argument, nullptr, 'ldw'},
    {"legendre-deg-x",     optional_argument, nullptr, 'ldx'},
    {"trace-deg-wave",     optional_argument, nullptr, 'tdw'},
    {"trace-deg-x",        optional_argument, nullptr, 'tdx'},
    {"psf-error",          optional_argument, nullptr, 'pe'},
    {"psf-core-wscale",    optional_argument, nullptr, 'pcw'},
    {"broken-fibers",      optional_argument, nullptr, 'bf'},
    {"variance-model",     optional_argument, nullptr, 'vm'},
    {"out-psf-xml",        optional_argument, nullptr, 'opx'},
    {"out-spots",          optional_argument, nullptr, 'os'},
    {"prior",              optional_argument, nullptr, 'p'},
    {"tmp_results",        optional_argument, nullptr, 'tr'},
#ifdef EXTERNAL_TAIL
    {"fit-psf-tails",      optional_argument, nullptr, 'fpt'},
#endif
#ifdef CONTINUUM
    {"fit-continuum",      optional_argument, nullptr, 'fc'},
#endif
    {"nlines",             optional_argument, nullptr, 'nl'}
  };
  */
  
  int i = 0;
  for (it = optmap.begin(); it !=optmap.end(); it++){
    long_opts[i] = {it->first.c_str(), it->second[0], nullptr, it->second[1]}; i++;
  }
  for (int i=noptions; i<noptsmax; i++) long_opts[i] = {"NULL", 0, nullptr, -1};

  option long_optse[] = {{"",optional_argument,nullptr,'Z'}};
  const char* const short_opts = "ha:v";
  
  specex_set_verbose(true);

  while (true)
    {
      const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

      if (opt == -1){
	break;
      } else if ( opt == argint(optmap, "help")){
      } else if (opt == argint(optmap, "arc") || opt == 'a'){
	arc_image_filename = string(optarg); 
      } else if (opt == argint(optmap, "in-psf")){
	input_psf_filename = string(optarg);
      } else if (opt == argint(optmap, "out-psf")){
	output_fits_filename = string(optarg);	  
      } else if (opt == argint(optmap, "lamp-lines")){
	lamp_lines_filename = string(optarg);	  
      } else if (opt == argint(optmap, "only-trace-fit")){
	fit_traces = true;
	fit_sigmas = true;
	fit_thepsf = true;	
      } else if (opt == argint(optmap, "no-trace-fit")){
	fit_traces = false;	  
      } else if (opt == argint(optmap, "no-sigma-fit")){
	fit_sigmas = false;
      } else if (opt == argint(optmap, "no-psf-fit")){
	fit_thepsf = false;
      } else if (opt == argint(optmap, "flux-hdu")){
	flux_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "ivar-hdu")){
	ivar_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "mask-hdu")){
	mask_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "header-hdu")){
	header_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "xcoord-hdu")){
	xtrace_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "ycoord-hdu")){
	ytrace_hdu = stoi(optarg);
      } else if (opt == argint(optmap, "first-bundle")){
	first_fiber_bundle = stoi(optarg);
      } else if (opt == argint(optmap, "last-bundle")){
	last_fiber_bundle = stoi(optarg);
      } else if (opt == argint(optmap, "single-bundle")){
	single_bundle = true;
      } else if (opt == argint(optmap, "first-fiber")){
	first_fiber = stoi(optarg);
      } else if (opt == argint(optmap, "last-fiber")){
	last_fiber = stoi(optarg);
      } else if (opt == argint(optmap, "half-size-x")){
	half_size_x = stoi(optarg);
	half_size_x_def = false;		  
      } else if (opt == argint(optmap, "half-size-y")){
	half_size_y = stoi(optarg); 
	half_size_y_def = false;
      } else if (opt == argint(optmap, "psfmodel")){
	psf_model = stoi(optarg);
      } else if (opt == argint(optmap, "positions")){
	fit_individual_spots_position = true;
      } else if (opt == argint(optmap, "quiet")){
	specex_set_verbose(false);
      } else if (opt == argint(optmap, "debug")){
	specex_set_debug(true);
      } else if (opt == argint(optmap, "trace-prior-deg")){
	trace_prior_deg = stoi(optarg);
      } else if (opt == argint(optmap, "lamplines")){
	lamp_lines_filename = string(optarg);
      } else if (opt == argint(optmap, "core")){
	specex_set_dump_core(true);
      } else if (opt == argint(optmap, "gauss-hermite-deg")){
	gauss_hermite_deg = stoi(optarg);
	gauss_hermite_deg_def = false;
      } else if (opt == argint(optmap, "gauss-hermite-deg2")){
	gauss_hermite_deg2 = stoi(optarg);
	gauss_hermite_deg2_def = false;
      } else if (opt == argint(optmap, "legendre-deg-x")){
	legendre_deg_x = stoi(optarg);
	legendre_deg_x_def = false;
      } else if (opt == argint(optmap, "legendre-deg-wave")){
	legendre_deg_wave = stoi(optarg);
	legendre_deg_wave_def = false;
      } else if (opt == argint(optmap, "trace-deg-x")){
	trace_deg_x = stoi(optarg);
	trace_deg_x_def = false;
      } else if (opt == argint(optmap, "trace-deg-wave")){
	trace_deg_wave = stoi(optarg);
	trace_deg_wave_def = false;
      } else if (opt == argint(optmap, "psf-error")){
	psf_error = stod(optarg);
      } else if (opt == argint(optmap, "psf-core-wscale")){
	psf_core_wscale = stod(optarg);
      } else if (opt == argint(optmap, "broken-fibers")){
	broken_fibers_string = string(optarg);
#ifdef EXTERNAL_TAIL
      } else if (opt == argint(optmap, "fit-psf-tails")){
	fit_psf_tails = true;
#endif
#ifdef CONTINUUM
      } else if (opt == argint(optmap, "fit-continuum")){
	fit_continuum = true;
#endif
      } else if (opt == argint(optmap, "variance-model")){
	use_variance_model = true;
      } else if (opt == argint(optmap, "out-psf-xml")){
	output_xml_filename = string(optarg);
      } else if (opt == argint(optmap, "out-spots")){
	output_spots_filename = string(optarg);
      } else if (opt == argint(optmap, "prior")){
	string sentence = string(optarg);
	istringstream iss(sentence);
	copy(istream_iterator<string>(iss),
	     istream_iterator<string>(),
	     back_inserter(argurment_priors));	  
      } else if (opt == argint(optmap, "tmp_results")){
	write_tmp_results = true;
      } else if (opt == argint(optmap, "nlines")){
	max_number_of_lines = stoi(optarg);
      }
    }

  return EXIT_SUCCESS;
  
}

#else

int specex::PyOptions::parse(int argc, char *argv[] ){

  // reading arguments
  // --------------------------------------------

  desc.add_options()
    ( "help,h", "display usage information" )
    ( "arc,a", popts::value<string>( &arc_image_filename ), "arc preprocessed fits image file name (mandatory)" )
    ( "in-psf", popts::value<string>( &input_psf_filename ), " input psf file name (fits or xml, can contain only the traces, mandatory)")  
    ( "out-psf", popts::value<string>( &output_fits_filename ), " output psf fits file name (mandatory)")
    ( "lamp-lines", popts::value<string>(&lamp_lines_filename ), " input lamp lines file name (mandatory)") 
    ( "only-trace-fit", "only fit the trace coordinates (default is trace,sigma,psf)")
    ( "no-trace-fit", "do not fit the trace coordinates (default is trace,sigma,psf)")
    ( "no-sigma-fit", "do not fit the gaussian sigma (default is trace,sigma,psf)")
    ( "no-psf-fit", "do not fit the psf (default is trace,sigma,psf)")
    ( "flux-hdu", popts::value<int>( &flux_hdu )->default_value(1), " flux hdu in input arc fits, for unusual data format")
    ( "ivar-hdu", popts::value<int>( &ivar_hdu )->default_value(2), " ivar hdu in input arc fits, for unusual data format")
    ( "mask-hdu", popts::value<int>( &mask_hdu )->default_value(3), " mask hdu in input arc fits, for unusual data format")
    ( "header-hdu", popts::value<int>( &header_hdu )->default_value(1), " header hdu in input arc fits, for unusual data format")
    ( "xcoord-hdu", popts::value<int>( &xtrace_hdu )->default_value(-1), "hdu of x trace legendre polynomial of wavelength (default is extension XTRACE)" )
    ( "ycoord-hdu", popts::value<int>( &ytrace_hdu )->default_value(-1), "hdu of y trace legendre polynomial of wavelength (default is extension YTRACE)" )
    ( "first-bundle", popts::value<int>( &first_fiber_bundle )->default_value(0), "first fiber bundle to fit")
    ( "last-bundle", popts::value<int>( &last_fiber_bundle )->default_value(0), "last fiber bundle to fit")
    ( "single-bundle", "fit a single bundle with all fibers")
    ( "first-fiber", popts::value<int>( &first_fiber )->default_value(0), "first fiber (must be in bundle)")
    ( "last-fiber", popts::value<int>( &last_fiber )->default_value(100000), "last fiber (must be in bundle)")
    ( "half-size-x", popts::value<int>( &half_size_x )->default_value(8), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half-size-y", popts::value<int>( &half_size_y )->default_value(5), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model )->default_value("GAUSSHERMITE"), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode (deprecated, true by default)" )
    ( "quiet", "no info message, only warning" )
    ( "debug", "turn on debug mode" )
    ( "trace-prior-deg", popts::value<int>( &trace_prior_deg )->default_value(0) , "force equal trace coeff in bundle starting at this degree")
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $SPECEXDATA/specex_linelist_desi.txt)" )
    ( "core", "dump core files when an exception is thrown" )
    ( "gauss-hermite-deg",  popts::value<int>( &gauss_hermite_deg )->default_value(6), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ("gauss-hermite-deg2",  popts::value<int>( &gauss_hermite_deg2 )->default_value(2), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)")
    //( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre-deg-wave",  popts::value<int>( &legendre_deg_wave )->default_value(3), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre-deg-x",  popts::value<int>( &legendre_deg_x )->default_value(1), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "trace-deg-wave",  popts::value<int>( &trace_deg_wave )->default_value(6), "degree of Legendre polynomials along wavelength for fit of traces")
    ( "trace-deg-x",  popts::value<int>( &trace_deg_x )->default_value(6), "degree of Legendre polynomials along x_ccd for fit of traces")
    ( "psf-error",  popts::value<double>( &psf_error )->default_value(0), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
    ( "psf-core-wscale",  popts::value<double>( &psf_core_wscale ), "scale up the weight of pixels in 5x5 PSF core")
    ( "broken-fibers", popts::value<string>( &broken_fibers_string ), "broken fibers (comma separated list)")    
#ifdef EXTERNAL_TAIL
    ( "fit-psf-tails", "unable fit of psf tails")
#endif
#ifdef CONTINUUM
    ( "fit-continuum", "unable fit of continuum")
#endif
    ( "variance-model", "refit at the end with a model of the variance to avoid Poisson noise bias")
    
    ( "out-psf-xml", popts::value<string>( &output_xml_filename ), " output psf xml file name")
   
    ( "out-spots", popts::value<string>( &output_spots_filename ), " output spots file name")  
    ( "prior", popts::value< vector<string> >( &argurment_priors )->multitoken(), " gaussian prior on a param : 'name' value error")  
    ( "tmp_results", " write tmp results")  
    ( "nlines", popts::value<int>( &max_number_of_lines )->default_value(200)," maximum number of emission lines used in fit (uses a algorithm to select best ones based on S/N and line coverage")  
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  try {    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
  }catch(std::exception& e) {
    cerr << "error in arguments : " << e.what() << endl;
    cerr << "try --help for options" << endl;
    return EXIT_FAILURE;
  }
  
  if ( ( argc < 2 ) || vm.count( "help" ) ) {
      cerr << endl;
      cerr << desc << endl;      
      return EXIT_FAILURE;
  }
  if ( ! vm.count( "arc" ) ) {
    cerr << endl;
    cerr << "missing --arc , try --help for options" << endl;      
    return EXIT_FAILURE;
  }
  if ( ! vm.count( "in-psf" ) ) {
    cerr << endl;
    cerr << "missing --in-psf , try --help for options" << endl;      
    return EXIT_FAILURE;
  }    
  if(lamp_lines_filename == "") {
    cerr << endl;
    cerr << "missing lamp_lines_filename either define env. variable SPECEXDATA or use option --lamplines" << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
  
}
#endif
