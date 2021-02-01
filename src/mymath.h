// Licensed under a 3-clause BSD style license - see LICENSE.rst

#define MYMATH_INT 39

#ifndef MYMATH_H
#define MYMATH_H

#include <cstdint>
#include <ctime>

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>

namespace specex {

// This class does mymath

class MyMath : public std::enable_shared_from_this <MyMath> {

    public :

        typedef std::shared_ptr <MyMath> pshr;

        MyMath(int ids);
        int add(int, int);
	int id;

};

}
#endif
