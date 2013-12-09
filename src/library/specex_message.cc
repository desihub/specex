#include "specex_message.h"

static bool static_specex_verbose = false; 
static bool static_specex_dump_core = false; 

void specex_set_verbose(bool yesorno) { static_specex_verbose=yesorno;}
void specex_set_dump_core(bool yesorno) { static_specex_dump_core=yesorno;}
bool specex_dump_core() {return static_specex_dump_core;}
bool specex_verbose() {return static_specex_verbose;}


void specex_info(const std::string& mess) {
  if(specex_verbose()) {
    std::cout << "INFO " << mess << std::endl;
  }
}
void specex_warning(const std::string& mess) {
  std::cerr << "WARNING " << mess << std::endl;
}

void specex_error(const std::string& mess) {
  std::cerr << "ERROR " << mess << std::endl;
  if(specex_dump_core()) {
    std::cerr << "dumping core ..." << std::endl;
    abort();
  }
}
