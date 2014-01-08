#include "InputAll.hh"

namespace phaser_io {


  InputAll::InputAll() : OnLine_(false)
  {
    Preprocessor capture(0,0);
    parseCCP4(capture);
  }


  InputAll::InputAll(Preprocessor& capture)  : OnLine_(false)
  { parseCCP4(capture); }

  // Entry to do nothing if "online" (ie if stdin not connected to file)
  InputAll::InputAll(const bool& OnLine)
  {
    Output output;  // dummy
    init(OnLine, false, output);
  }
  // Entry to do nothing if "online" (ie if stdin not connected to file)
  InputAll::InputAll(const bool& OnLine, Output& output)
  {
    init(OnLine, true, output);
  }
  void InputAll::init(const bool& OnLine, const bool& echo, Output& output)
  {
    OnLine_ = OnLine;
      if (! OnLine)
	{
	  Preprocessor capture(0,0);
	  if (echo) {
	    output.logTab(0,phaser_io::LOGFILE, ">>>>> Input command lines <<<<<\n\n");
	    output.logKeywords(phaser_io::LOGFILE, capture.Echo());
	    output.logTab(0,phaser_io::LOGFILE, ">>>>>     End of input    <<<<<\n\n");
	  }
	  parseCCP4(capture);
	}
  }

  InputAll::~InputAll() {}

  //performs the analysis between correlated keyword Input
  void InputAll::Analyse()
  {
    HKLIN::analyse();
    XDSIN::analyse();
    HKLOUT::analyse();
    HKLREF::analyse();
    XMLOUT::analyse();
    CHIRAL::analyse();
    RESO::analyse();
    ISIG::analyse();
    LABI::analyse();
    LABR::analyse();
    ORIG::analyse();
    LAUE::analyse();
    SPAC::analyse();
    REIN::analyse();
    TOLE::analyse();
    ACCE::analyse();
    CENT::analyse();
    PARTIALS::analyse();
    SYSABS::analyse();
    NAME::analyse();
    COPY::analyse();
    TEST::analyse();
    NOTEST::analyse();
    ASSU::analyse();
    NOASSU::analyse();
    EXCLUDE::analyse();
    CHOOSE::analyse();
    NEIGHBOUR::analyse();
    SETTING::analyse();
    OUTPUT::analyse();
    MULTIPLY::analyse();
    SCORE::analyse();
    ALLOW::analyse();
    WAVELENGTH::analyse();
    BLANK::analyse();
    WILDFILE::analyse();
    POLARISATION::analyse();

    if (RESO::ValidRESO())
      ISIG::setISIG(-1.0);  // Set MinIsig flag to indicate that
			    // resolution limits have been set
  }
} // phaser_io
