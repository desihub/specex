#include <specex_stamp.h>

std::ostream& operator << (std::ostream &stream, const specex::Stamp &stamp) { 
  stamp.write(stream); 
  return stream;
}
