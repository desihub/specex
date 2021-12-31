#include "specex_message.h"

static bool static_specex_debug = false;
static bool static_specex_verbose = false;
static bool static_specex_dump_core = false;
static std::string static_message_prefix = "";

void specex_set_message_prefix(const std::string &mess) { static_message_prefix = " "+mess;}
void specex_set_debug(bool yesorno) { static_specex_debug=yesorno;}
void specex_set_verbose(bool yesorno) { static_specex_verbose=yesorno;}
void specex_set_dump_core(bool yesorno) { static_specex_dump_core=yesorno;}
bool specex_dump_core() {return static_specex_dump_core;}
bool specex_is_verbose() {return static_specex_verbose;}
bool specex_is_debug() {return static_specex_debug;}

void specex_debug(const std::string& mess) {
  if(specex_is_debug()) {
    std::cout << "DEBUG" << static_message_prefix << " " << mess << std::endl;
  }
}
void specex_info(const std::string& mess) {
  if(specex_is_verbose()) {
    std::cout << "INFO" << static_message_prefix << " " << mess << std::endl;
  }
}
void specex_warning(const std::string& mess) {
  std::cerr << "WARNING" << static_message_prefix << " " << mess << std::endl;
}

void specex_error(const std::string& mess) {
  std::cerr << "ERROR" << static_message_prefix << " " << mess << std::endl;
  if(specex_dump_core()) {
    std::cerr << "dumping core ..." << std::endl;
    abort();
  }
}
