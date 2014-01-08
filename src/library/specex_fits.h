#ifndef SPECEX_FITS__H
#define SPECEX_FITS__H

namespace specex {

  class image_data;

  void write_new_fits_image(std::string const & path, size_t n_rows, size_t n_cols, const harp::vector_double& data);
  void write_new_fits_image(std::string const & path, const image_data& img);
  void write_new_fits_image(std::string const & path, const harp::matrix_double& mat);
  
  void read_fits_image(std::string const & path, int hdu, image_data& img);
  void read_fits_image(std::string const & path, image_data& img_in_hdu1);
  void read_fits_images(std::string const & path, image_data& img_in_hdu1, image_data& img_in_hdu2);



  
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
    int SizeOfVectorOfDouble() const;
    
  };
  
  class FitsTableEntry {
  public :
    std::string string_val; // can be empty
    harp::vector_double double_vals; // can be empty
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
    bool HasKey(const std::string& key) const;
    
    double DoubleKeyValue(const std::string& key) const;
    std::string StringKeyValue(const std::string& key) const;
    int IntKeyValue(const std::string& key) const;
    
    int DumpKeys(std::ostream& stream) const;
    
  
  };

  typedef std::map<std::string,FitsColumnDescription>::const_iterator ColumnDescriptionConstIterator;
  typedef std::map<std::string,FitsColumnDescription>::iterator ColumnDescriptionIterator;
  
}

#endif
