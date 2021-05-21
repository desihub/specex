#ifndef SPECEX_MESSAGE__H
#define SPECEX_MESSAGE__H

#include <iostream>
#include <string>
#include <sstream>

#include <unhrp.h>

void specex_set_debug(bool yesorno);
void specex_set_verbose(bool yesorno);
void specex_set_dump_core(bool yesorno);
bool specex_is_verbose();
bool specex_is_debug();
bool specex_dump_core();

void specex_debug(const std::string& mess);
void specex_info(const std::string& mess);
void specex_warning(const std::string& mess);
void specex_error(const std::string& mess);

#define SPECEX_DEBUG(mess) {std::stringstream ss; ss << mess; specex_debug(ss.str()); }
#define SPECEX_INFO(mess) {std::stringstream ss; ss << mess; specex_info(ss.str()); }
#define SPECEX_WARNING(mess) {std::stringstream ss; ss << mess; specex_warning(ss.str()); }
#define SPECEX_ERROR(mess) {std::stringstream ss; ss << mess; if(specex_dump_core()) {ss << " (at line " << __LINE__ << " of file " << __FILE__ << ")"; specex_error(ss.str());} else {HARP_THROW(ss.str().c_str()); }}

#endif
