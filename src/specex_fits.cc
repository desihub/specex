#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

//#include <unhrp.h>
#include <specex_image_data.h>
#include <specex_fits.h>
#include <specex_message.h>

using namespace std;


int specex::find_hdu( fitsfile *fp, const std::string& extname, const std::string& alternate_extname) {
  //SPECEX_DEBUG("Looking for extension '" << extname << "'");
  
  int nhdus = harp::fits::nhdus(fp);
  map<string,int> extname_in_hdus;
  for(int i=0;i<nhdus;i++) {
    int hdu=i+1; // starts at 1
    int status = 0;
    fits_movabs_hdu ( fp, hdu, NULL, &status );
    
    string extname_in_hdu="";
    try {
      harp::fits::key_read (fp,"EXTNAME",extname_in_hdu);
    }catch(...) { 
      if(i!=0) SPECEX_WARNING("could not read EXTNAME in hdu " << hdu);
      continue;
    }
    extname_in_hdus[extname_in_hdu]=hdu;
  }
  map<string,int>::iterator it=extname_in_hdus.find(extname);
  if(it != extname_in_hdus.end()) {
    SPECEX_DEBUG("found '" << extname << "' in hdu "<< it->second);
    return it->second;
  }
  if(alternate_extname != "") {
    it=extname_in_hdus.find(alternate_extname);
    if(it != extname_in_hdus.end()) {
      SPECEX_DEBUG("found '" << alternate_extname << "' in hdu "<< it->second);
      return it->second;
    }
  }

  if(alternate_extname != "") {
    SPECEX_WARNING("Didn't find extname '" << extname << "' nor '" << alternate_extname << "' in file");
  } else {
    SPECEX_WARNING("Didn't find extname '" << extname << "' in file"); 
  }
  return -1;
}


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
    // cout << "dim='" << substrings[s] << "'" << endl;
    dimension.push_back(atoi(substrings[s].c_str()));
  }

  /*
  SPECEX_INFO("decode_dimension '" << tdim << "'");
  for(size_t d = 0 ; d < dimension.size(); d++) {
    SPECEX_INFO("d[" << d << "]=" << dimension[d]);
  }
  */
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


////////////////////////////////////////



void specex::write_new_fits_image(std::string const & path, size_t n_rows, size_t n_cols, const unhrp::vector_double& data) {
  
  if(data.size() != (n_rows*n_cols))
    SPECEX_ERROR("Incompatible number of columns x rows and vector size");

  fitsfile * fp;  
  harp::fits::create ( fp, path );
  harp::fits::img_append < double > ( fp, n_rows, n_cols );
  harp::fits::img_write ( fp, data, false );
  harp::fits::close ( fp );
  
  
  
}

void specex::write_new_fits_image(std::string const & path, const image_data& img) {
  specex::write_new_fits_image(path, img.n_rows(), img.n_cols() , img.data);
  SPECEX_INFO("wrote image in file " << path);
}

void specex::write_new_fits_image(std::string const & path, const unhrp::matrix_double& mat) {
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
  
  int bitpix;
  harp::fits::key_read(fp ,"BITPIX", bitpix);
  int status = 0;
  int compressed = fits_is_compressed_image(fp,&status);
  harp::fits::check ( status );
  
  harp::fits::img_read ( fp, img.data , false); // calls  fits_read_pix  , works with compressed BITPIX=8 for cfitsio 337
    
  harp::fits::close ( fp ); 
  SPECEX_INFO("read one image in '" << path << "' hdu=" << hdu << " size=" << nrows << "x" << ncols << " bitpix=" << bitpix << " compressed=" << compressed);
  
  
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
  harp::fits::img_read ( fp, img_in_hdu1.data, false );
  harp::fits::img_seek ( fp, 2);
  harp::fits::img_dims ( fp, nrows, ncols );
  img_in_hdu2.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
  harp::fits::img_read ( fp, img_in_hdu2.data, false );
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
bool specex::FitsColumnDescription::IsInt() const {
  if(format.empty()) return false;
  if(format=="UNKNOWN") return false;
  char last = format[format.size()-1];
  return (last=='I' || last=='U' || last=='J'|| last=='V' || last=='K');
}
int specex::FitsColumnDescription::SizeOfVectorOfDouble() const {
  if(!IsDouble()) return 0;
  if(format.size()==1) return 1;
  return atoi(format.substr(0,format.size()-1).c_str());
}
int specex::FitsColumnDescription::SizeOfVectorOfInt() const {
  if(!IsInt()) return 0;
  if(format.size()==1) return 1;
  return atoi(format.substr(0,format.size()-1).c_str());
}



#define CHECKERROR if(status) {fits_report_error(stdout, status); SPECEX_ERROR("fits error");}

specex::FitsTable::FitsTable() : fptr(0) {
}

specex::FitsTable::FitsTable(const string& filename, int hdu_number, bool verbose) : fptr(0) {
  Read(filename,hdu_number);
}


void specex::FitsTable::AddColumnDescription(const string& ttype, const string& tform,  const string& tdim, const string& tunit) {
  FitsColumnDescription column;
  column.unit=tunit;
  column.format=tform;
  column.dimension=decode_dimension(tdim);
  column.col = int(columns.size());
  columns[ttype]=column;
}


bool specex::FitsTable::Write(fitsfile *fp) const {
  
  int status = 0;
  int nrows  = data.size();
  int ncols  = columns.size();
  
  char **ttype = new char*[ncols]; for(int c=0;c<ncols;c++) ttype[c] = new char[8];
  char **tform = new char*[ncols]; for(int c=0;c<ncols;c++) tform[c] = new char[8];
  char **tunit = new char*[ncols]; for(int c=0;c<ncols;c++) tunit[c] = new char[8];
  char **tdim = new char*[ncols]; for(int c=0;c<ncols;c++) tdim[c] = new char[8];
  
  for(std::map<std::string,FitsColumnDescription>::const_iterator it = columns.begin(); it!=columns.end(); ++it) {
    sprintf(ttype[it->second.col],"%s",it->first.c_str());
    sprintf(tform[it->second.col],"%s",it->second.format.c_str());
    sprintf(tunit[it->second.col],"%s",it->second.unit.c_str());
    sprintf(tdim[it->second.col],"%s",encode_dimension(it->second.dimension).c_str());  
  }
  
  
  char extname[] = "";
  int tfields = ncols;
  
  fits_create_tbl( fp, BINARY_TBL, nrows, tfields, ttype, tform,
		   tunit, extname, &status);
  CHECKERROR;
  
  for(int c=0;c<ncols;c++) {
    if(strlen(tdim[c])>0) {
      char key[8];
      sprintf(key,"TDIM%d",c+1);
      harp::fits::key_write(fp,key,tdim[c],"dimension");
    }
  }
 
  // now we need to write the data .........

    // Automatic data type conversion is performed for numerical data types (only) if the data type of the column (defined by the TFORMn keyword) differs from the data type of the array in the calling routine. ASCII and binary tables support the following data type values: TSTRING, TBYTE, TSBYTE, TSHORT, TUSHORT, TINT, TUINT, TLONG, TLONGLONG, TULONG, TFLOAT, or TDOUBLE. Binary tables also support TLOGICAL (internally mapped to the `char' data type), TCOMPLEX, and TDBLCOMPLEX. 

    
    int firstrow  = 1;  // first row in table to write 
    int firstelem = 1; // first element in row  (ignored in ASCII tables)
    
    for(std::map<std::string,FitsColumnDescription>::const_iterator it = columns.begin(); it!=columns.end(); ++it) {
      
      
      const FitsColumnDescription& column = it->second;
      
      if(column.IsString()) {
	char ** string_vals= new char*[nrows];
	for(int r=0;r<nrows;r++) {
	  const char *input_val = data[r][column.col].string_val.c_str();
	  string_vals[r] = new char[strlen(input_val)];
	  strcpy(string_vals[r],input_val);
	}
	fits_write_col(fp, TSTRING, column.col+1, firstrow, firstelem, nrows, string_vals, &status);
	CHECKERROR;

	// delete
	for(int r=0;r<nrows;r++) {
	  delete [] string_vals[r];
	}
	delete [] string_vals;
	

      }else if(column.IsDouble()) {
	
	int nd = column.SizeOfVectorOfDouble();
	
	double * double_vals = new double[nrows*nd];
	
	for(int r=0;r<nrows;r++) {
	  
	  if(nd !=  int(data[r][column.col].double_vals.size())) SPECEX_ERROR("specex::FitsTable::Write inconsistent data size " << nd << " != " << int(data[r][column.col].double_vals.size()));
	  
	  for(int d=0;d<nd;d++) {
	    
	    double_vals[r*nd+d] = data[r][column.col].double_vals(d); // ordering checked with pyfits
	    //double_vals[r+d*nrows] = data[r][column.col].double_vals(d);
	  
	  }
	}
	//SPECEX_INFO("specex::FitsTable::Write column " << column.col << " with " << nrows*nd << " data");
	fits_write_col(fp, TDOUBLE, column.col+1, firstrow, firstelem, nrows*nd, double_vals, &status);
	CHECKERROR;
	
	delete [] double_vals;

      }else if(column.IsInt()) {
	
	int nd = column.SizeOfVectorOfInt();
	
	int * int_vals = new int[nrows*nd];
	
	for(int r=0;r<nrows;r++) {
	  
	  if(nd !=  int(data[r][column.col].int_vals.size())) SPECEX_ERROR("specex::FitsTable::Write inconsistent data size " << nd << " != " << int(data[r][column.col].int_vals.size()));
	  
	  for(int d=0;d<nd;d++) {
	    
	    int_vals[r*nd+d] = data[r][column.col].int_vals(d); // ordering checked with pyfits
	    
	  }
	}
	//SPECEX_INFO("specex::FitsTable::Write column " << column.col << " with " << nrows*nd << " data");
	fits_write_col(fp, TINT, column.col+1, firstrow, firstelem, nrows*nd, int_vals, &status);
	CHECKERROR;
	
	delete [] int_vals;

      }else{
	SPECEX_ERROR("specex::FitsTable::Write format not yet implemented");
      }
    }
    return 0;
}
void specex::FitsTable::Read(const string& filename, int hdu_number, bool verbose)  {
  int status = 0;
  fits_open_file(&fptr, filename.c_str() , READONLY, &status);
  CHECKERROR;
  fits_movrel_hdu(fptr, hdu_number, NULL, &status);
  CHECKERROR;
  Read(fptr,verbose);
}

void specex::FitsTable::Read(fitsfile *fp, bool verbose)  {
  fptr = fp;
  int status = 0;
  
  SPECEX_DEBUG("Starting specex::FitsTable::Read");

  // checking HDU TYPE
  int hdutype;
  fits_get_hdu_type(fptr, &hdutype, &status);
  CHECKERROR;
  if(hdutype != BINARY_TBL && hdutype != ASCII_TBL ) {
    SPECEX_ERROR("FitsTable::Read : hdu is not a table");
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
      column.col = c;
      
      sprintf(key,"TTYPE%d", c+1);
      string name = StringKeyValue(key);

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
	column.dimension = decode_dimension(StringKeyValue(key));
      } else {
	column.dimension.clear();
	column.dimension.push_back(1);
      }
      
      columns[name]=column;
    }
  }
  
  
  //for(size_t c=0;c<columns.size();c++) {
  for(ColumnDescriptionConstIterator it=columns.begin(); it != columns.end(); ++it) {
    SPECEX_INFO("FitsColumnDescription : " << it->first << " " << it->second.format);
  }
  
  
  char** strptr  = new char*[ncols]; 
  strptr[0] = new char[256];
  
  data.resize(nrows);
  
  SPECEX_DEBUG("specex::FitsTable::Read start reading rows");
  
  // loop on rows
  for(int r=0; r<nrows; r++) {
    
    std::vector<specex::FitsTableEntry>& entries = data[r];
    entries.resize(columns.size());
    
    for(ColumnDescriptionConstIterator it=columns.begin(); it != columns.end(); ++it) {
      const specex::FitsColumnDescription &col = it->second;
      int c = col.col;
      specex::FitsTableEntry& entry = entries[c];
      
      if(col.IsString()) {
	//cout << "col " << c << " is a string" << endl;

	char nullval[16] = "none";
	
	int anynul;
	fits_read_col(fptr, TSTRING, c+1, long(r+1), 1,1, nullval, strptr, &anynul, &status); 
	CHECKERROR;
	//cout << strptr[0] << endl;
	entry.string_val = strptr[0]; // can we pass directly the address ?
	
	if(r==0) SPECEX_INFO("first FitsTableEntry col=" << it->first << " : '" << entry.string_val << "'");
	
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
	  SPECEX_INFO("first FitsTableEntry col=" << it->first << " : " << entry.double_vals(0));
	}
      } else if(col.IsInt()) {
	int nvals = col.SizeOfVectorOfInt();
	
	int *values = new int[nvals];
	int nullval = 0;
	
	int anynul;
	//fits_read_col(fptr, TDOUBLE, c+1, long(r+1), 1,nvals, &nullval, entry.double_vals.NonConstData(), &anynul, &status); 
	fits_read_col(fptr, TINT, c+1, long(r+1), 1,nvals, &nullval, values, &anynul, &status); 
	CHECKERROR;
	
	entry.int_vals.resize(nvals);
	for(int i=0;i<nvals;i++)
	  entry.int_vals(i)=values[i];
	delete [] values;

	if(r==0) {
	  SPECEX_INFO("first FitsTableEntry col=" << it->first << " : " << entry.int_vals(0));
	}
      } else {
	SPECEX_ERROR("data type not implemented in FitsTable");
      }
    }
    
  }
  SPECEX_DEBUG("specex::FitsTable::Read done reading rows");
  delete[] strptr[0];
  delete[] strptr;
}


bool specex::FitsTable::HasKey(const string& key) const {
  int status = 0;
  char a_C_string[80];
  fits_read_key(fptr, TSTRING, const_cast<char*>(key.c_str()), a_C_string, NULL, &status);
  return (status != KEY_NO_EXIST);
}

double specex::FitsTable::DoubleKeyValue(const string& key) const  {
  double val;
  int status = 0;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>(key.c_str()), &val, NULL, &status);
  CHECKERROR;
  return val;
}

string specex::FitsTable::StringKeyValue(const string& key) const  {
  int status = 0;
  char a_C_string[80];
  fits_read_key(fptr, TSTRING, const_cast<char*>(key.c_str()), a_C_string, NULL, &status);
  CHECKERROR;
  return string(a_C_string);
}

int specex::FitsTable::IntKeyValue(const string& key)  const {
  int val;
  int status = 0;
  fits_read_key(fptr, TINT, const_cast<char*>(key.c_str()), &val, NULL, &status);
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
