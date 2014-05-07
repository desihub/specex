#ifndef SPECEX_MESSAGE__H
#define SPECEX_MESSAGE__H

#include <iostream>
#include <string>
#include <sstream>

#ifdef USE_MPI
#  include <harp_mpi.hpp>
#else
#  include <harp.hpp>
#endif

void specex_set_verbose(bool yesorno);
void specex_set_dump_core(bool yesorno);
bool specex_verbose();
bool specex_dump_core();

void specex_info(const std::string& mess);
void specex_warning(const std::string& mess);
void specex_error(const std::string& mess);

#define SPECEX_INFO(mess) {std::stringstream ss; ss << mess; specex_info(ss.str()); }
#define SPECEX_WARNING(mess) {std::stringstream ss; ss << mess << " (at line " << __LINE__ << " of file " << __FILE__ << ")"; specex_warning(ss.str()); }
#define SPECEX_ERROR(mess) {std::stringstream ss; ss << mess; if(specex_dump_core()) {ss << " (at line " << __LINE__ << " of file " << __FILE__ << ")"; specex_error(ss.str());} else {HARP_THROW(ss.str().c_str()); }}

#endif
