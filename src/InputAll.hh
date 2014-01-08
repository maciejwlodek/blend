#ifndef INPUTALL_HEADER
#define INPUTALL_HEADER

#include "keywords.hh"

namespace phaser_io {
  class InputAll :
    public HKLIN,
    public XDSIN,
    public SCAIN,
    public MERGED,
    public HKLOUT,
    public HKLREF,
    public XMLOUT,
    public XYZIN,
    public CHIRAL,
    public RESO,
    public ISIG,
    public LABI,
    public LABR,
    public ORIG,
    public LAUE,
    public SPAC,
    public REIN,
    public TOLE,
    public ACCE,
    public CENT,
    public PARTIALS,
    public SYSABS,
    public NAME,
    public COPY,
    public CELL,
    public TEST,
    public NOTEST,
    public ASSU,
    public NOASSU,
    public EXCLUDE,
    public CHOOSE,
    public NEIGHBOUR,
    public SETTING,
    public OUTPUT,
    public MULTIPLY,
    public SCORE,
    public ALLOW,
    public WAVELENGTH,
    public BLANK,
    public WILDFILE,
    public POLARISATION
  {
  public:
    InputAll(Preprocessor&); 
    InputAll();
    // Entry to do nothing if "online"  ie if stdin not connected to file)
    InputAll(const bool& OnLine);
    InputAll(const bool& OnLine, Output& output);  // ... with echo to output
    ~InputAll();
    void Analyse();

    // true is input has been read (ie not online)
    bool Valid() const {return !(OnLine_);}

  private:
    bool OnLine_;

    void init(const bool& OnLine, const bool& echo, Output& output);

  };
    
} // phaser_io
  
#endif

