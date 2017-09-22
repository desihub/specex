#ifndef SPECEX_FITS__H
#define SPECEX_FITS__H

#include <boost/numeric/ublas/io.hpp>

namespace specex {

  class image_data;

  void write_new_fits_image(std::string const & path, size_t n_rows, size_t n_cols, const harp::vector_double& data);
  void write_new_fits_image(std::string const & path, const image_data& img);
  void write_new_fits_image(std::string const & path, const harp::matrix_double& mat);
  
  void read_fits_image(std::string const & path, int hdu, image_data& img);
  void read_fits_image(std::string const & path, image_data& img_in_hdu1);
  void read_fits_images(std::string const & path, image_data& img_in_hdu1, image_data& img_in_hdu2);


  int find_hdu( fitsfile *fp, const std::string& extname, const std::string& alternate_extname="");
  
  
  // the following shoudl be somewhere in harp : it's about reading fits tables
  
  class FitsColumnDescription {
  public :
    int col;
    //std::string name;
    std::string unit;
    std::string format;
    std::vector<int> dimension;
    
    bool IsString() const;
    bool IsDouble() const;
    bool IsInt() const;
    int SizeOfVectorOfDouble() const;
    int SizeOfVectorOfInt() const;
    
  };
  
  class FitsTableEntry {
  public :
    std::string string_val; // can be empty
    boost::numeric::ublas::vector<double> double_vals; // can be empty
    boost::numeric::ublas::vector<int>    int_vals; // can be empty
  };
  
  class FitsTable {
    
  private :
    fitsfile *fptr;
    
  public :
  
    std::map<std::string,FitsColumnDescription> columns;
    std::vector< std::vector<FitsTableEntry> > data;
     

    FitsTable();
    FitsTable(const std::string& filename, int hdu_number, bool verbose = false);
    void Read(const std::string& filename, int hdu_number, bool verbose = false);
    void Read(fitsfile *fp, bool verbose = false);
    bool HasKey(const std::string& key) const;
    bool Write(fitsfile *fp) const;
    int  IntKeyValue(const std::string& key) const;
    double DoubleKeyValue(const std::string& key) const;
    std::string StringKeyValue(const std::string& key) const;
    
    std::vector<int> decode_dimension(const std::string& tdim) const;
    std::string encode_dimension(const std::vector<int>& dimension) const;

    int DumpKeys(std::ostream& stream) const;
    void AddColumnDescription(const std::string& ttype, const std::string& tform,  const std::string& tdim="", const std::string& tunit="");
  };

  typedef std::map<std::string,FitsColumnDescription>::const_iterator ColumnDescriptionConstIterator;
  typedef std::map<std::string,FitsColumnDescription>::iterator ColumnDescriptionIterator;
  
}

#endif
