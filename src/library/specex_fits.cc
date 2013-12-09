#include <iostream>
#include <cstdio>
#include <string>

#include <harp.hpp>

#include <specex_image_data.h>
#include <specex_fits.h>
#include <specex_message.h>

using namespace std;


void specex::write_new_fits_image(std::string const & path, size_t n_rows, size_t n_cols, const harp::vector_double& data) {
  
  if(data.size() != (n_rows*n_cols))
    SPECEX_ERROR("Incompatible number of columns x rows and vector size");

  fitsfile * fp;  
  harp::fits::create ( fp, path );
  harp::fits::img_append < double > ( fp, n_rows, n_cols );
  harp::fits::img_write ( fp, data );
  harp::fits::close ( fp );
  
  
  
}

void specex::write_new_fits_image(std::string const & path, const image_data& img) {
  specex::write_new_fits_image(path, img.n_rows(), img.n_cols() , img.data);
  SPECEX_INFO("wrote image in file " << path);
}

void specex::write_new_fits_image(std::string const & path, const harp::matrix_double& mat) {
  specex::write_new_fits_image(path, mat.size1(), mat.size2() , mat.data()); 
  SPECEX_INFO("wrote matrix in file " << path);
}
  
void specex::read_fits_image(std::string const & path, int hdu, image_data& img) {
  
  fitsfile * fp; 
  harp::fits::open_read ( fp, path);
  harp::fits::img_seek ( fp, hdu);
  
  size_t nrows,ncols;
  
  harp::fits::img_dims ( fp, nrows, ncols );
  
  img.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row

  harp::fits::img_read ( fp, img.data );
  harp::fits::close ( fp ); 
  SPECEX_INFO("read one image in '" << path << "' hdu=" << hdu);
}

void specex::read_fits_image(std::string const & path, image_data& img) {
  read_fits_image(path,1,img); 
  
}

void specex::read_fits_images(std::string const & path, image_data& img_in_hdu1, image_data& img_in_hdu2) {
  fitsfile * fp; 
  harp::fits::open_read ( fp, path);
  harp::fits::img_seek ( fp, 1);
  size_t nrows,ncols;
  harp::fits::img_dims ( fp, nrows, ncols );
  img_in_hdu1.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
  harp::fits::img_read ( fp, img_in_hdu1.data );
  harp::fits::img_seek ( fp, 2);
  harp::fits::img_dims ( fp, nrows, ncols );
  img_in_hdu2.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
  harp::fits::img_read ( fp, img_in_hdu2.data );
  harp::fits::close ( fp ); 
  SPECEX_INFO("read 2 images in '" << path << "'");
}


// the following should be better written and placed in HARP


bool specex::FitsColumnDescription::IsString() const {
  if(format.empty()) return false;
  if(format=="UNKNOWN") return false;
  char last = format[format.size()-1];
  return (last=='A');
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



#define CHECKERROR if(status) {fits_report_error(stdout, status); SPECEX_ERROR("fits error");}

specex::FitsTable::FitsTable() : fptr(0) {
}

specex::FitsTable::FitsTable(const string& filename, int hdu_number, bool verbose) : fptr(0) {
  Read(filename,hdu_number);
}

void specex::FitsTable::Read(const string& filename, int hdu_number, bool verbose)  {
  int status = 0;
  fits_open_file(&fptr, filename.c_str() , READONLY, &status);
  CHECKERROR;
  fits_movrel_hdu(fptr, hdu_number, NULL, &status);
  CHECKERROR;
  
  // checking HDU TYPE
  int hdutype;
  fits_get_hdu_type(fptr, &hdutype, &status);
  CHECKERROR;
  if(hdutype != BINARY_TBL && hdutype != ASCII_TBL ) {
    SPECEX_ERROR("FitsTable::Read file " << filename << " hdu " << hdu_number << " is not a table");
  }
  
  long nrows;
  int ncols;
  fits_get_num_rows(fptr, &nrows, &status); CHECKERROR; 
  fits_get_num_cols(fptr, &ncols, &status); CHECKERROR;

  // Columns Descriptors :
  
  {
    char key[12];
    for (int c=0; c<ncols; c++) {
      
      specex::FitsColumnDescription column;
      
      sprintf(key,"TTYPE%d", c+1);
      column.name = StringKeyValue(key);

      sprintf(key,"TUNIT%d", c+1);
      if(HasKey(key))
	column.unit = StringKeyValue(key);
      else
	column.unit = "UNKNOWN";

      sprintf(key,"TFORM%d", c+1);
      if(HasKey(key))
	column.format = StringKeyValue(key);
      else
	column.format = "UNKNOWN";
      
      sprintf(key,"TDIM%d", c+1);
      if(HasKey(key)) {
	column.dimension=StringKeyValue(key);
      } else {
	column.dimension="(1)";
      }

      columns.push_back(column);
    }
  }
  
  
  for(size_t c=0;c<columns.size();c++) {
    const specex::FitsColumnDescription &column = columns[c];
    SPECEX_INFO("FitsColumnDescription col=" << c << " " << column.name << " " << column.format << " " << column.dimension);
  }
  
  
  char** strptr  = new char*[ncols]; 
  strptr[0] = new char[256];

  data.resize(nrows);
  
  // loop on rows
  for(int r=0; r<nrows; r++) {
    //{std::vector<FitsTableEntry> toto; data.push_back(toto);}
    std::vector<specex::FitsTableEntry>& entries = data[r];
    entries.resize(columns.size());
    
    for(size_t c=0;c<columns.size();c++) {
      const specex::FitsColumnDescription &col = columns[c];
      
      //{FitsTableEntry toto; entries.push_back(toto);}
      specex::FitsTableEntry& entry = entries[c];
      
      
      if(col.IsString()) {
	//cout << "col " << c << " is a string" << endl;

	char nullval[16] = "none";
	
	int anynul;
	fits_read_col(fptr, TSTRING, c+1, long(r+1), 1,1, nullval, strptr, &anynul, &status); 
	CHECKERROR;
	//cout << strptr[0] << endl;
	entry.string_val = strptr[0]; // can we pass directly the address ?
	
	if(r==0) SPECEX_INFO("first FitsTableEntry col=" << c << " : '" << entry.string_val << "'");
	
      } else if(col.IsDouble()) {
	int nvals = col.SizeOfVectorOfDouble();
	
	double *values = new double[nvals];
	double nullval = 0.;
	
	int anynul;
	//fits_read_col(fptr, TDOUBLE, c+1, long(r+1), 1,nvals, &nullval, entry.double_vals.NonConstData(), &anynul, &status); 
	fits_read_col(fptr, TDOUBLE, c+1, long(r+1), 1,nvals, &nullval, values, &anynul, &status); 
	CHECKERROR;
	
	entry.double_vals.resize(nvals);
	for(int i=0;i<nvals;i++)
	  entry.double_vals(i)=values[i];
	delete [] values;

	if(r==0) {
	  SPECEX_INFO("first FitsTableEntry col=" << c << " : " << entry.double_vals(0));
	  if(entry.double_vals.size()>1) SPECEX_INFO(entry.double_vals(1) << ",...");
	}
      } else {
	SPECEX_ERROR("data type not implemented in FitsTable");
      }
    }
    
  }
  delete[] strptr[0];
  delete[] strptr;
}


bool specex::FitsTable::HasKey(const string& key) const {
  int status = 0;
  char a_C_string[80];
  fits_read_key(fptr, TSTRING, key.c_str(), a_C_string, NULL, &status);
  return (status != KEY_NO_EXIST);
}

double specex::FitsTable::DoubleKeyValue(const string& key) const  {
  double val;
  int status = 0;
  fits_read_key(fptr, TDOUBLE, key.c_str(), &val, NULL, &status);
  CHECKERROR;
  return val;
}

string specex::FitsTable::StringKeyValue(const string& key) const  {
  int status = 0;
  char a_C_string[80];
  fits_read_key(fptr, TSTRING, key.c_str(), a_C_string, NULL, &status);
  CHECKERROR;
  return string(a_C_string);
}

int specex::FitsTable::IntKeyValue(const string& key)  const {
  int val;
  int status = 0;
  fits_read_key(fptr, TINT, key.c_str(), &val, NULL, &status);
  CHECKERROR;
  return val;
}

int specex::FitsTable::DumpKeys(ostream& stream) const {
  int nkeys = 0;
  int status = 0;
  fits_get_hdrspace(fptr, &nkeys, NULL, &status); CHECKERROR;  
  for (int ii = 1; ii <= nkeys; ii++) {
    char card[81];
    fits_read_record(fptr, ii, card, &status);
    stream << "@ FITS_" << card << endl;
  }
  return 0;
}
