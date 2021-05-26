#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <specex_pyoptions.h>

using namespace std;

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

  // reset the index of the option to be read from argv in case get_opt is invoked again
  optind = 1;
  
  return EXIT_SUCCESS;
  
}
