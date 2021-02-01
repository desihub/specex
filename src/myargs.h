// Licensed under a 3-clause BSD style license - see LICENSE.rst

#define MYARGS_INT 38

#ifndef MYARGS_H
#define MYARGS_H

#include <iostream>
#include <cstdint>
#include <ctime>

#include <string>
#include <vector>
#include <map>
#include <memory>

namespace specex {

// This class does myargs

class MyArgs : public std::enable_shared_from_this <MyArgs> {

    public :

        typedef std::shared_ptr <MyArgs> pshr;

	MyArgs(int id);
        //void get_args(int argc, char **argv );
	void get_args(int argc, char **argv);
	void print_args();
	int ids;
	int argc;
	std::vector<std::string> argv; 

};

}
#endif
