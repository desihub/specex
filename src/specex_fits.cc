#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

#include <specex_unbls.h>

#include <specex_image_data.h>
#include <specex_fits.h>
#include <specex_message.h>

using namespace std;

static vector<string> DecomposeString(const string &Source,vector<string> &Tokens)
{
  vector<string> SubStrings;

  string source_copy = Source;
  
  string::size_type start = 0;

  string::size_type tokpos = source_copy.npos;
  for(size_t t=0;t<Tokens.size();t++) {
    string::size_type pos = source_copy.find_first_of(Tokens[t],start);
    if(pos<tokpos) tokpos=pos;
  }

  while (tokpos != (string::size_type) source_copy.npos)
    {
      string sub = source_copy.substr(start,tokpos-start);
      if (sub.length() > 0) SubStrings.push_back(sub);
      start = tokpos+1;
      
      tokpos = source_copy.npos;
      for(size_t t=0;t<Tokens.size();t++) {
	string::size_type pos = source_copy.find_first_of(Tokens[t],start);
	if(pos<tokpos) tokpos=pos;
      }
    }
  
  // get the last element
  string sublast = source_copy.substr(start, source_copy.length());
  if (sublast.length() > 0) SubStrings.push_back(sublast);
  
  return SubStrings;
}

std::vector<int> specex::FitsTable::decode_dimension(const string& tdim) const {
  
  vector<string> tokens; 
  tokens.push_back(" "); 
  tokens.push_back(","); 
  tokens.push_back("("); 
  tokens.push_back(")"); 
  vector<string> substrings = DecomposeString(tdim,tokens);
  std::vector<int> dimension;
  for(size_t s=0;s<substrings.size();s++) {
    dimension.push_back(atoi(substrings[s].c_str()));
  }

  return dimension;
}

string specex::FitsTable::encode_dimension(const std::vector<int>& dimension) const {
  if(dimension.size()==0) return "";
  std::stringstream tdim;
  tdim << "( ";
  for(size_t i=0;i<dimension.size();i++) {
    tdim << dimension[i];
    if(i<dimension.size()-1) tdim << ", ";
  }
  tdim << ")";
  return string(tdim.str());
}

bool specex::FitsColumnDescription::IsDouble() const {
  if(format.empty()) return false;
  if(format=="UNKNOWN") return false;
  char last = format[format.size()-1];
  return (last=='D' || last=='E');
}
int specex::FitsColumnDescription::SizeOfVectorOfDouble() const {
  if(!IsDouble()) return 0;
  if(format.size()==1) return 1;
  return atoi(format.substr(0,format.size()-1).c_str());
}

specex::FitsTable::FitsTable(){}

void specex::FitsTable::AddColumnDescription(const string& ttype, const string& tform,  const string& tdim, const string& tunit) {
  FitsColumnDescription column;
  column.unit=tunit;
  column.format=tform;
  column.dimension=decode_dimension(tdim);
  column.col = int(columns.size());
  columns[ttype]=column;
}
