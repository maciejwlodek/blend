// keywords.cpp

#include "keywords.hh"
#include "pointgroup.hh"
#include "latsym.hh"
#include "scala_util.hh"
#include "jiffy.hh"
#include "expandfilename.hh"
#include "string_util.hh"
#include "lattice.hh"
#include "hkl_symmetry.hh"

namespace phaser_io {

//--------------------------------------------------------------
  bool SetLabel(const std::string& lkey,
		const std::string& label,
		int& nlab,
		std::string& Ilabel,
		std::string& sigIlabel)
    // return false if error
  {
    if ((lkey == "F" || lkey == "I") ||
	(lkey == "" && nlab == 0)) {
      // store I|F label
      Ilabel = label;
    } else if ((lkey == "SIGF" || lkey == "SIGI") ||
	       (lkey == "" && nlab == 1)) {
      // store sigI|F label
      sigIlabel = label;
    } else {
      return false;
    }
    nlab += 1;
    return true;
  }
  //--------------------------------------------------------------
  std::string CheckGroup(const std::string& Group)
  // Strip quotes & spaces from group name Group, check that it corresponds to
  // to a valid group (fail if not) & return group name
  {
    std::string name = Group;
    // strip out quotes (& spaces)
    name = StringUtil::Strip(name,'\'');
    name = StringUtil::Strip(name,'"');
    scala::SpaceGroup SG(name);
    if (SG.Null()) {
      throw SyntaxError("","ERROR in input command, invalid group name "+name);
    }
    /*
      try {
      std::cout << "Name1: " << name << "\n"; //^
      CCtbxSym::SpaceGroupName(name);
      //      name = CCtbxSym::SpaceGroupName(name);
      std::cout << "Name2: " << name << "\n"; //^
      CCtbxSym::PointGroup LG(name);
      }
      catch (std::exception const& err) {
      throw SyntaxError("ERROR in input command, invalid group name ",name);
      }
    */
    return name;
  }
  //--------------------------------------------------------------
  CHIRAL::CHIRAL() : CCP4base(), InputBase()
{
  Add_Key("CHIRAL");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  Chiral = scala::UNKNOWN;
}

  //--------------------------------------------------------------
Token_value CHIRAL::parse(std::istringstream& input_stream)
  // Syntax  CHIRAL  CHIRAL | NONCHIRAL | CENTROSYMMETRIC
{
  Chiral = scala::NONCHIRAL;
  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NAME))
	{
	  if (keyIs("CHIR"))
	    {Chiral = scala::CHIRAL;}
	  else if (keyIs("NONC"))
	    {Chiral = scala::NONCHIRAL;}
	  else if (keyIs("CENT")) 
	    {Chiral = scala::CENTROSYMMETRIC;}
	  else 
	    {
	      throw SyntaxError(keywords,"unrecognised keyword");
	    }
	}
    }
  return skip_line(input_stream);
}

//--------------------------------------------------------------
Chirality CHIRAL::getCHIRAL() const
{
  return Chiral;
}

  //--------------------------------------------------------------
void CHIRAL::analyse(void)
{
}
//--------------------------------------------------------------
RESO::RESO() : CCP4base(), InputBase()
{
  Add_Key("RESO");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  Valid = false;
}
//--------------------------------------------------------------
Token_value RESO::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  RESOlution  <high> [<low>]  either order
  //  RESOlution  HIGH <high>
  //  RESOlution  LOW  <low>
  double high = input_resolution_range.min();
  double  low = input_resolution_range.max();
  bool onefound = false;

  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NAME))
	{
	  if (keyIs("HIGH")) high = get1num(input_stream);
	  else if (keyIs("LOW")) low = get1num(input_stream);
	}
      else
	{
	  if (onefound)
	    {
	      low = number_value;
	    }
	  else
	    {
	      high = number_value;
	      onefound = true;
	    }
	}
    }
  setRESO(ResoRange(low,high));
  return skip_line(input_stream);
}
//--------------------------------------------------------------
void RESO::setRESO(const ResoRange& reso_range)
{
  input_resolution_range = reso_range;
  Valid = true;
}
//--------------------------------------------------------------
ResoRange RESO::getRESO() const {return input_resolution_range;}
//--------------------------------------------------------------

static float MinIsigRatioDefault = 6.0;  // Default value

ISIG::ISIG() : CCP4base(), InputBase()
{
  Add_Key("ISIG");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  MinIsigRatio = MinIsigRatioDefault;
}
//--------------------------------------------------------------
Token_value ISIG::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  ISIGLIMIT <minimum<I>/<sigmaI>>

  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  setISIG(number_value);
	}
    }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
void ISIG::setISIG(const float& MinIsig)
{
  MinIsigRatio = MinIsig;
}
//--------------------------------------------------------------
float ISIG::getISIG() const {return MinIsigRatio;}
//--------------------------------------------------------------
LABI::LABI() : CCP4base(), InputBase()
{
  Add_Key("LABI");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value LABI::parse(std::istringstream& input_stream)
{
  // Syntax: LABIN [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  int nlab = 0;
  int proglab = -1; // Initially no "program label" read
  std::string label;
  std::string lkey;

  while (get_token(input_stream) != ENDLINE) {
    if (proglab >= 0) {
      // Have potential program label, is next field "="?
      if (curr_tok == ASSIGN) {
	lkey = label;
	label = "";
	proglab = +1;
      } else {
	// No, store as actual label
	if (proglab > 0) label = string_value;
	if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
	  throw SyntaxError(keywords,"Too many column labels");
	if (proglab > 0) {
	  // "=" was read, just had label
	  proglab = -1;
	} else {
	  // current string is next key or label
	  proglab = 0;
	  label = string_value;
	  lkey = "";
	}
      }
    } else {
      label = string_value;
      lkey = "";
      proglab = 0;
    }
  }
  if (proglab == 0)
    if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
      throw SyntaxError(keywords,"Too many column labels");
  if (sigFIlabel != "" && FIlabel == "") {
      throw SyntaxError(keywords,"Can't specify just SIG");
  }

  return skip_line(input_stream);
}

//--------------------------------------------------------------
LABR::LABR() : CCP4base(), InputBase()
{
  Add_Key("LABR");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value LABR::parse(std::istringstream& input_stream)
{
  // Syntax: LABREF [F|I = ] <F|Ilabel> [[SIGF|I = ] <sigF|Ilabel>]
  int nlab = 0;
  int proglab = -1; // Initially no "program label" read
  std::string label;
  std::string lkey;

  while (get_token(input_stream) != ENDLINE)
    {
      if (proglab >= 0)
	{
	  // Have potential program label, is next field "="?
	  if (curr_tok == ASSIGN)
	    {
	      lkey = label;
	      label = "";
	      proglab = +1;
	    }
	  else
	    // No, store as actual label
	    {
	      if (proglab > 0) label = string_value;
	      if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
		throw SyntaxError(keywords,"Too many column labels");
	      if (proglab > 0) 
		// "=" was read, just had label
		{
		  proglab = -1;
		}
	      else
		// current string is next key or label
		{
		  proglab = 0;
		  label = string_value;
		  lkey = "";
		}
	    }
	}
      else
	{
	  label = string_value;
	  lkey = "";
	  proglab = 0;
	}
    }
  if (proglab == 0)
    if (! SetLabel(lkey, label, nlab, FIlabel, sigFIlabel))
      throw SyntaxError(keywords,"Too many column labels");
  if (sigFIlabel != "" && FIlabel == "") {
      throw SyntaxError(keywords,"Can't specify just SIG");
  }

  return skip_line(input_stream);
}

//--------------------------------------------------------------
ORIG::ORIG() : CCP4base(), InputBase(), OriginalLat(false)
{
  Add_Key("ORIG");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value ORIG::parse(std::istringstream& input_stream)
{
  OriginalLat = true;
  return skip_line(input_stream);
}
//--------------------------------------------------------------
LAUE::LAUE() : CCP4base(), InputBase(),  Lauegroup("")
{
  Add_Key("LAUE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value LAUE::parse(std::istringstream& input_stream)
{
  // Check for string HKLIN or ALL, force to full-length uppercase
  Lauegroup = getLine(input_stream);
  std::string thiskey(phaser_io::stoup(Lauegroup));

  if (thiskey.find("HKLIN") < thiskey.size())
    {Lauegroup = "HKLIN";}
  else if (thiskey.find("ALL") < thiskey.size())
    {Lauegroup = "ALL";}
  else
    {
      // Check for valid name
      Lauegroup = CheckGroup(Lauegroup);
    }
  return ENDLINE;
}
//--------------------------------------------------------------
SPAC::SPAC() : CCP4base(), InputBase(),  SPGname("")
{
  Add_Key("SPAC");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value SPAC::parse(std::istringstream& input_stream)
{
  // Check for string HKLIN, force to full-length uppercase
  SPGname = getLine(input_stream);
  std::string thiskey(phaser_io::stoup(SPGname));

  if (thiskey.find("HKLIN") < thiskey.size())
    {SPGname = "HKLIN";}
  else
    {
      SPGname = StringUtil::Strip(SPGname,'\''); // strip out quotes & space
      SPGname = StringUtil::Strip(SPGname,'"');
      scala::SpaceGroup SG(SPGname);
      if (SG.Null()) {
	throw SyntaxError("ERROR in SPACEGROUP command: ", SPGname);
      }
      SPGname = SG.Symbol_hm();

      ///      SPGname = CCtbxSym::FixCentreLatSymbol(SPGname, 'R');
	///      SPGname = CCtbxSym::SpaceGroupName(SPGname, 'H', false); // not reference setting
	///      try {
	///	CCtbxSym::PointGroup SG(SPGname);
	///      }
	///      catch (std::exception const& err) {
	///	throw SyntaxError("ERROR in SPACEGROUP command: ", SPGname);
	///      }
    }
  return ENDLINE;
}
//--------------------------------------------------------------
  REIN::REIN() : CCP4base(), InputBase(), Reindex(), set(false), left(false)
{
  Add_Key("REIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value REIN::parse(std::istringstream& input_stream)
// Syntax REINDEX [LEFTHAND] <reindex operator>
//        REINDEX MATRIX H11 H21 H31 H12 H22 H32 H13 H23 H33
//          ie matrix in the same order as the XDS REIDX command
{
  std::string Op = "";
  bool matrix = false;
  std::vector<double> xm(9);  // XDS matrix
  int i = 0;
  while (get_token(input_stream) != ENDLINE)  {

    if (keyIs("LEFT"))  {
      left = true;
    } else if (keyIs("MATRIX"))  {
	matrix = true;
    } else {
      if (matrix) {
	if (tokenIs(1,NUMBER)) {
	  xm[i++] = number_value;
	} else {
	  throw SyntaxError
	    (keywords,
	     "ERROR in command REINDEX: non-number in matrix");
	}
      } else {
	Op += string_value;
      }
    }
  }
    if (matrix) {
      if (i != 9) {
	throw SyntaxError
	  (keywords,
	   "ERROR in command REINDEX: 9 numbers expected in matrix");
      }
      // Matrix for ReindexOp post-multiplies hkl row vector
      Mat33<double> H(xm[0], xm[3], xm[6],
		      xm[1], xm[4], xm[7],
		      xm[2], xm[5], xm[8]);
      Reindex = scala::ReindexOp(H);
    } else {
      Reindex = scala::ReindexOp(Op);
    }
  if (!Reindex.Positive() && !left)
    throw SyntaxError(keywords,"ILLEGAL reindex operator inverts hand");

  set = true;

  return ENDLINE;
}
//--------------------------------------------------------------

static float ToleranceDefault = 2.0;  // Default value

TOLE::TOLE() : CCP4base(), InputBase()
{
  Add_Key("TOLE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  Tolerance = ToleranceDefault;
}
//--------------------------------------------------------------
Token_value TOLE::parse(std::istringstream& input_stream)
{
  // Syntax: TOLERANCE  <LatticeTolerance>
  //   tolerance in degrees for test lattice to fit input cell

  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  setTOLERANCE(number_value);
	}
    }
  return skip_line(input_stream);
}
//--------------------------------------------------------------

static float AcceptanceDefault = 0.2;  // Default value

ACCE::ACCE() : CCP4base(), InputBase()
{
  Add_Key("ACCE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  acceptlimit = AcceptanceDefault;
}
//--------------------------------------------------------------
Token_value ACCE::parse(std::istringstream& input_stream)
{
  // Syntax: ACCEPT <AcceptanceFraction> 
  int i = 0;
  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  if (i == 0)
	    {
	      if (number_value > 0.0 && number_value <=1.0)
		setACCEPTANCELIMIT(number_value);
	      else
		throw InputError
		  ("ERROR in command ACCEPT",
		   " value must be between 0 & 1.0, not"+
		   phaser_io::dtos(number_value));
	    }
	  i++;
	}
    }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
CENT::CENT() : centre(false)
{
  Add_Key("CENT");
  gridmax = clipper::Vec3<int>(2,2,2);
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value CENT::parse(std::istringstream& input_stream)
{
  centre = true;
  int i = 0;

  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  if (i < 3 && std::abs(number_value) < 10)
	    {
	      gridmax[i] = Nint(number_value);
	    }
	  else
	    {throw SyntaxError
		(keywords,
		 "ERROR in command CENTRE: three numbers must be given, < 10");
	    }
	  i++;
	}
    }
  if (i != 0 && i != 3)
    {throw SyntaxError
	(keywords,
	 "ERROR in command CENTRE: either no or three numbers must be given, < 10");
    }
  return skip_line(input_stream);
}

//--------------------------------------------------------------
PARTIALS::PARTIALS()
{
  Add_Key("PARTIALS");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  // Default values
  minfraclim = 0.95;
  maxfraclim = 1.05;
  minscalefrac = 0.5;
  check = true;
  maxgap = 0;
}
//--------------------------------------------------------------
Token_value PARTIALS::parse(std::istringstream& input_stream)
// Syntax: PARTials [TEST <lower_limit> <upper_limit>] [CORRECT <minimum_fraction>]
//  [[NO]CHECK] [[NO]GAP <gap limit>]
{
  bool limits = false;
  int nlim = 0;
  bool scalelim = false;
  int nsclim = 0;
  bool gapval = false;

  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NAME))
	{
	  limits = false;
	  scalelim = false;
	  gapval = false;
	  if (keyIs("TEST"))
	    {
	      limits = true;
	    }  
	  else if (keyIs("CORRECT"))
	    {
	      scalelim = true;
	    }  
	  else if (keyIs("CHECK"))
	    {
	      check = true;
	    }  
	  else if (keyIs("NOCHECK"))
	    {
	      check = false;
	    }
	  else if (keyIs("GAP"))
	    {
	      gapval = true;
	      maxgap = 1;
	    }
	  else if (keyIs("NOGAP"))
	    {
	      maxgap = 0;
	    }
	  else
	    {throw SyntaxError
		(keywords, "unrecognised keyword");}
	}
      else if (tokenIs(1,NUMBER))
	{
	  if (limits)
	    {
	      if (nlim == 0)
		{minfraclim = number_value;}
	      else if (nlim == 1)
		{maxfraclim = number_value;}
	      else
		{throw SyntaxError
		    (keywords, "ERROR in TEST two numbers must be given");}
	      nlim++;
	    }
	  else if (scalelim)
	    {
	      if (nsclim == 0)
		{minscalefrac = number_value;}
	      else
		{throw SyntaxError
		    (keywords,
		     "ERROR in CORRECT one numbers must be given");}
	      nsclim++;
	    }
	  else if (gapval)
	    {
	      maxgap = Nint(number_value);
	      if (maxgap < 0 || maxgap > 5)
		{
		  throw SyntaxError(keywords,"unreasonable MaxGap ");
		}
	    }
	  else
	    {throw SyntaxError
		(keywords, "unexpected number");}
	}
    }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
SYSABS::SYSABS() : CCP4base(), InputBase()
{
  Add_Key("SYST");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  check = true;
}
//--------------------------------------------------------------
Token_value SYSABS::parse(std::istringstream& input_stream)
{
  check = getBoolean(input_stream);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
HKLIN::HKLIN() : CCP4base(), InputBase()
{
  Add_Key("HKLIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  wild = false;
  fileseries = 0;
}
//--------------------------------------------------------------
void HKLIN::storeNames(const scala::ExpandFileName& filenames)
{
  std::vector<FileNameTime> fnames = filenames.NameTimes(); // all names & times
  if (fnames.size() > 0) {
    fileseries++;
    for (size_t i=0;i<fnames.size();i++) {
      names.push_back(fnames[i]);
      pxdnames.push_back(NAME::getNAME());
      seriesnumber.push_back(fileseries);  // store file series number
    }
  }
}
//--------------------------------------------------------------
Token_value HKLIN::parse(std::istringstream& input_stream)
{
  std::string fname = getLine(input_stream);
  // If WILDFILE command has been given, see if the filenames given
  // can be treated as if it had wild-cards from Mosflm
  if (WILDFILE::WildFile()) {
    fname = MakeWildTemplate(fname, WILDFILE::Template());
  }
  // Expand file name if it contains wild-cards "*" or "?"
  scala::ExpandFileName filenames(fname, "mtz");
  //  wild set true if wild-cards in name
  wild = wild || filenames.Wild();
  storeNames(filenames);  // store names etc
  seriesnames.push_back(fname);  // raw file name
  return ENDLINE;
}
//--------------------------------------------------------------
std::string HKLIN::MakeWildTemplate(const std::string& fname,
				    const std::string& ftemplate)
// return fname modified to include wildcard characters "?"
// Match fname to template from end,
//  replace characters equivalent to "#" in template by "?" in fname
{
  // Check from end of template, looking for '#'
  std::string newname = fname;
  int i;
  int j = -1;
  int n = 0;
  for (i=int(ftemplate.size())-1;i>=0;i--) {
    if (ftemplate[i] == '#') {
      if (j < 0) {
	j = i;
      }
      n++; // count # characters
    }
  }
  if (n == 0) {
    Message::message(Message_fatal("No '#' characters in file template\n  "+
				   ftemplate));
  }
  // Store the character after the last '#'
  char marker = ' ';
  if (j < int(ftemplate.size())-1) {marker = ftemplate[j+1];}
  // Now search the filename backwards for the marker
  for (i=fname.size()-1;i>=0;i--) {
    if (fname[i] == marker) break;
  }
  // marker at point i, replace previous n characters with ?
  for (j=i-1;j>=i-n;j--) {
    // Check that this character is a digit
    if (!std::isdigit(newname[j])) {
      Message::message(Message_fatal(
	"'#' characters in file template don't match digit in filename\n  "+
	fname));
    }
    newname[j] = '?';
  }
  return newname;
}
//--------------------------------------------------------------
int HKLIN::AddHKLINfilename(const std::vector<std::string>& Names)
// Stores file names from command line
// returns number of files added
{
  if (Names.size() <= 0) return 0;
  //  wild set true if > 1 file, assumed to come from series
  wild = (Names.size() > 1);
  scala::ExpandFileName filenames(Names);  // get times
  storeNames(filenames);  // store names etc
  seriesnames.push_back("on command line");  // raw file name
  return Names.size();
}
//--------------------------------------------------------------
std::vector<std::string>  HKLIN::getHKLINnames() const
{
  std::vector<std::string> fnames;
  for (size_t i=0;i<names.size();i++) {
    fnames.push_back(names[i].FileName());
  }
  return fnames;
}
//--------------------------------------------------------------
std::vector<std::string> HKLIN::RejectOutofSequenceHKLINFiles()
// Strip out files in all series which are out of time order
// ie where "later" files in order have earlier modification times
// returns list of filenames rejected
{
  std::vector<std::string> rejects;
  int nseries = seriesnames.size();
  int nfiles = names.size();
  if (nfiles <= 0 ) return rejects;
  for (int iser=0;iser<nseries;++iser) { // loop file series
    // find files in this series (series numbered from 1)
    int i1 = -1;
    int i2 = -1;
    for (int ifile=0;ifile<nfiles;++ifile) {
      if (seriesnumber[ifile] == iser+1) {
	if (i1 < 0) {
	  i1 = ifile; // start of series
	}
      } else if (i1 >= 0 && i2 < 0) {
	i2 = ifile;  // end of series
	break;
      }
    }
    if (i2 < 0) i2 = nfiles-1;
    // Series runs from files i1 to i2
    if ((i2-i1) < 2) continue;  // only 1 file in series
    // Loop files in series
    time_t mtime = names[i1].ModTime();  // time for 1st file
    for (int ifile=i1+1;ifile<=i2;++ifile) {
      //^
      //^      std::cout << "ROOSF " << names[ifile].FileName()
      //^		<<  " time " << names[ifile].ModTime() << "\n";
      if (names[ifile].ModTime() < mtime) {
	// this file is older than the previous one,
	// reject all the rest in the series & list them
	std::vector<std::string> r = Reject(ifile, i2);
	rejects.insert(rejects.end(), r.begin(), r.end());
	break;
      }
      mtime = names[ifile].ModTime();
    }
  } // end loop series
  return rejects;
}
//--------------------------------------------------------------
std::vector<std::string> HKLIN::Reject(const int& i1, const int& i2)
// reject & remove from list files with indices i1 to i2 (inclusive)
// return list of filenames rejected
{
  std::vector<std::string> rejects;
  for (int ifile=i1;ifile<=i2;++ifile) {
    rejects.push_back(names[ifile].FileName());
  }
  names.erase(names.begin()+i1, names.begin()+i2+1);
  pxdnames.erase(pxdnames.begin()+i1, pxdnames.begin()+i2+1);
  seriesnumber.erase(seriesnumber.begin()+i1, seriesnumber.begin()+i2+1);
  return rejects;
}
//--------------------------------------------------------------
XDSIN::XDSIN() : CCP4base(), InputBase()
{
  Add_Key("XDSIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value XDSIN::parse(std::istringstream& input_stream)
{
  std::string fname = StringUtil::Trim(getLine(input_stream));
  names.push_back(fname);
  pxdnames.push_back(NAME::getNAME());
  return ENDLINE;
}
//--------------------------------------------------------------
int XDSIN::AddXDSINfilename(const std::string& Name)
// Stores file names from command line, series 1
// returns number of files added
{
  if (Name == "") return 0;
  names.push_back(Name);
  pxdnames.push_back(NAME::getNAME());
  return 1;
}
//--------------------------------------------------------------
SCAIN::SCAIN() : CCP4base(), InputBase()
{
  Add_Key("SCAIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value SCAIN::parse(std::istringstream& input_stream)
{
  std::string fname = StringUtil::Trim(getLine(input_stream));
  names.push_back(fname);
  pxdnames.push_back(NAME::getNAME());
  return ENDLINE;
}
//--------------------------------------------------------------
int SCAIN::AddSCAINfilename(const std::string& Name)
// Stores file names from command line, series 1
// returns number of files added
{
  if (Name == "") return 0;
  names.push_back(Name);
  pxdnames.push_back(NAME::getNAME());
  return 1;
}
//--------------------------------------------------------------
MERGED::MERGED() : CCP4base(), InputBase()
{
  Add_Key("MERGED");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
  spacegroupName = "";
}
//--------------------------------------------------------------
Token_value MERGED::parse(std::istringstream& input_stream)
{
  spacegroupName = StringUtil::Trim(getLine(input_stream));
  if (spacegroupName == "") {spacegroupName = "UNKNOWN";}
  return ENDLINE;
}
//--------------------------------------------------------------
HKLOUT::HKLOUT() : CCP4base(), InputBase()
{
  Add_Key("HKLOUT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value HKLOUT::parse(std::istringstream& input_stream)
{
  name = getLine(input_stream);
  scala::AddFileExtension(name, "mtz");
  return ENDLINE;
}
//--------------------------------------------------------------
HKLREF::HKLREF() : CCP4base(), InputBase()
{
  Add_Key("HKLREF");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value HKLREF::parse(std::istringstream& input_stream)
{
  name = StringUtil::Trim(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
XMLOUT::XMLOUT() : CCP4base(), InputBase()
{
  Add_Key("XMLOUT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value XMLOUT::parse(std::istringstream& input_stream)
{
  name = StringUtil::Trim(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
XYZIN::XYZIN() : CCP4base(), InputBase()
{
  Add_Key("XYZIN");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value XYZIN::parse(std::istringstream& input_stream)
{
  name = StringUtil::Trim(getLine(input_stream));
  return ENDLINE;
}
//--------------------------------------------------------------
NAME::NAME() : CCP4base(), InputBase()
{
  Add_Key("NAME");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);  
}
//--------------------------------------------------------------
Token_value NAME::parse(std::istringstream& input_stream)
{
  // Syntax: NAME PROJECT <Pname> CRYSTAL <Xname> DATASET <Dname>
  std::string pname = "";
  std::string xname = "";
  std::string dname = "";
  int which = -1;
  while (get_token(input_stream) != ENDLINE) {
    if (keyIs("PROJECT")) {
      which = 0;
    } else if (keyIs("CRYSTAL")) {
      which = 1;
    } else if (keyIs("DATASET")) {
      which = 2;
    } else {
      // store string
      if (which == 0) {pname = string_value;}
      else if (which == 1) {xname = string_value;}
      else if (which == 2) {dname = string_value;}
      else {
	Message::message(Message_fatal("NAME: Unrecognised subkey "+string_value));
      }
    }
  }
  pxdname = scala::PxdName(pname, xname, dname);
  return ENDLINE;
}
//--------------------------------------------------------------
scala::PxdName NAME::pxdname;
//--------------------------------------------------------------
void NAME::setNAME(const scala::PxdName& PXDname) 
// Only reset blank fields
  {
    if (pxdname.pname() == "") {
      pxdname.pname() = PXDname.pname();
    }
    if (pxdname.xname() == "") {
      pxdname.xname() = PXDname.xname();
    }
    if (pxdname.dname() == "") {
      pxdname.dname() = PXDname.dname();
    }
  }
//--------------------------------------------------------------
COPY::COPY() : CCP4base(), InputBase()
{
  Add_Key("COPY");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  copyflag = false;
}
//--------------------------------------------------------------
Token_value COPY::parse(std::istringstream& input_stream)
{
  // Syntax: COPY
  copyflag = true;
  return ENDLINE;
}
//--------------------------------------------------------------
CELL::CELL() : CCP4base(), InputBase()
{
  Add_Key("CELL");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value CELL::parse(std::istringstream& input_stream)
{
  int i = 0;
  std::vector<Dtype> vcell(6);
  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  vcell[i++] = number_value;
	}
    }
  if (i != 6)
    {throw SyntaxError
	(keywords,
	 "ERROR in command CELL: six numbers must be given");
    }
  cell = scala::Scell(vcell);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
bool TEST::testfirstfile = false;         // default no test
//--------------------------------------------------------------
TEST::TEST() : CCP4base(), InputBase()
{
  Add_Key("TEST");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value TEST::parse(std::istringstream& input_stream)
{
  setTESTFIRSTFILE(true);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
NOTEST::NOTEST() : CCP4base(), InputBase()
{
  Add_Key("NOTEST");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value NOTEST::parse(std::istringstream& input_stream)
{
  TEST::setTESTFIRSTFILE(false);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
int ASSU::assumesameindexing = -1;  // = -1 unset, = 0 false, = +1 true
//--------------------------------------------------------------
ASSU::ASSU() : CCP4base(), InputBase()
{
  Add_Key("ASSU");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value ASSU::parse(std::istringstream& input_stream)
{
  bool check = true;
  while (get_token(input_stream) != ENDLINE)
    {
      if (keyIs("ON")) {check = true;}
      else if (keyIs("OFF")) {check = false;}
      else {
	throw SyntaxError(keywords,"flag should be ON or OFF");
      }
    }
  setASSUMESAMEINDEXING(check);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
void ASSU::setASSUMESAMEINDEXING(const bool& AssumeSame)
{
  if (AssumeSame) {
    assumesameindexing = +1;
  } else {
    assumesameindexing =  0;
  }
}
//--------------------------------------------------------------

NOASSU::NOASSU() : CCP4base(), InputBase()
{
  // Syntax: DONOTASSUMESAMEINDEXING
  //   Do not assume all HKLIN files have the same indexing scheme
  Add_Key("DONOTASSU");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value NOASSU::parse(std::istringstream& input_stream)
{
  ASSU::setASSUMESAMEINDEXING(false);
  return skip_line(input_stream);
}
//--------------------------------------------------------------
EXCLUDE::EXCLUDE() : CCP4base(), InputBase()
{
  Add_Key("EXCLUDE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
}
//--------------------------------------------------------------
Token_value EXCLUDE::parse(std::istringstream& input_stream)
// Syntax:
//   EXCLUDE DATASET <datasetname> | <crystalname>/<datasetname>
//     exclude dataset
//   EXCLUDE BATCH [FILE|SERIES  <Jfile>]  <b1> <b2> <b3> ... | <b1> TO <b2>
//     exclude batch list or range: if FILE key present, b1 etc
//     refer to original (file) batch numbers (SERIES is the equivalent for a
//     series of files defined with a wild-card)
{
  int Select = 0; // 0 nothing; +1 batch selection; -1 dataset selection
  int fileSeries = 0; // >0 if FILE|SERIES keywords given
  bool brange = false;
  std::vector<int> bnum;

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME) && (keyIs("FILE") || keyIs("SERIES"))) {
	// Keyword FILE or SERIES (synonymous)
	fileSeries = Nint(get1num(input_stream));
	if (fileSeries < 0) {
	  throw SyntaxError(keywords,"FILE | SERIES value must be > 0");
	}
    } else if (Select == 0) {
      if (tokenIs(1,NAME)) {
	if (keyIs("BATCH")) {
	  Select = +1;
	  continue;
	} else if (keyIs("DATASET")) {
	  Select = -1;
	  continue;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      }
    } else if (Select == +1) {
      // Batch selection
      if (tokenIs(1,NAME)) {
	if (keyIs("TO")) {
	  // Range
	  if (brange) {
	    throw SyntaxError
	      (keywords,"only one batch range allowed per EXCLUDE command");
	  }
	  brange = true;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      } else {
	// gather numbers
	bnum.push_back(Nint(number_value));
      }
    } else if (Select == -1) {
      throw SyntaxError(keywords,"EXCLUDE DATASET option not yet implemented");
    }
  }
  // Line finished, store results
  if (Select == 0) {
    throw SyntaxError(keywords,"subkeyword BATCH must be given");
    //    throw SyntaxError(keywords,"subkeyword BATCH or DATASET must be given");
  }
  if (brange) {
    if (bnum.size() != 2) {
      throw SyntaxError(keywords,"batch range must be given as 'n1 TO n2'");
    }
    batchexclude.AddRange(bnum[0], bnum[1], fileSeries);
  } else {
    // List
    for (size_t i=0;i<bnum.size();i++) {
      batchexclude.AddBatch(bnum[i], fileSeries);
    }
  }  
  return skip_line(input_stream);
}
//--------------------------------------------------------------
CHOOSE::CHOOSE() : CCP4base(), InputBase(),
		   lauegroup(""), spacegroup("")
{
  Add_Key("CHOOSE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  solution_number = 0;
}
//--------------------------------------------------------------
Token_value CHOOSE::parse(std::istringstream& input_stream)
// Syntax:
//  CHOOSE [SOLUTION] <N>
//  CHOOSE LAUEGROUP <group>
//  CHOOSE SPACEGROUP <group>
//  choose solution from Laue group or space group list
{
  //  bool solution = true;
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (keyIs("SOLUTION")) {
	continue;
      } else if (keyIs("LAUEGROUP")) {
	// read & check group name
	lauegroup = CheckGroup(getLine(input_stream));
	break;
      } else if (keyIs("SPACEGROUP")) {
	// read & check group name
	spacegroup = CheckGroup(getLine(input_stream));
	break;
      }
    } else if (tokenIs(1,NUMBER)) {
      solution_number = Nint(number_value);
      if (solution_number <= 0) {
	throw SyntaxError(keywords,"solution number must be > 0");
      }
      break;
    }
  }
  return ENDLINE;
}
//--------------------------------------------------------------
NEIGHBOUR::NEIGHBOUR() : CCP4base(), InputBase()
{
  Add_Key("NEIGHBOUR");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  neighbourfraction = -1.0;
}
//--------------------------------------------------------------
Token_value NEIGHBOUR::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  NEIGHBOUR <neighbourfraction>
  // Set fraction of neighbouring intensities to subtract from
  // axial intensities in systematic absence analysis
  int i = 0;
  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NUMBER))
	{
	  if (i == 0)
	    {
	      if (number_value >= 0.0 && number_value <= 0.2)
		neighbourfraction = number_value;
	      else
		throw InputError
		  ("ERROR in command NEIGHBOUR",
		   " value must be between 0 & 0.2, not"+
		   phaser_io::dtos(number_value));
	    }
	  i++;
	}
    }
  return skip_line(input_stream);
}
  //--------------------------------------------------------------
SETTING::SETTING() : CCP4base(), InputBase()
{
  Add_Key("SETTING");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  ITCsetting = 0; //  SETTING CELL-BASED default
}
  //--------------------------------------------------------------
Token_value SETTING::parse(std::istringstream& input_stream)
  // Syntax:
  //   SETTING CELL-BASED | SYMMETRY-BASED | C2
{
  ITCsetting = true; //  SETTING CELL-BASED default
  while (get_token(input_stream) != ENDLINE)
    {
      if (tokenIs(1,NAME))
	{
	  if (keyIs("CELL")) {
	    ITCsetting = 0;
	  } else if (keyIs("SYMM")) {
	    ITCsetting = +1;
	  } else if (keyIs("C2")) {
	    ITCsetting = -1;
	  } else {
	      throw SyntaxError(keywords,"unrecognised keyword");
	  }
	}
    }
  return skip_line(input_stream);
}
  //--------------------------------------------------------------
OUTPUT::OUTPUT() : CCP4base(), InputBase()
{
  Add_Key("OUTPUT");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  filename_summed = "";
}
  //--------------------------------------------------------------
Token_value OUTPUT::parse(std::istringstream& input_stream)
  // Syntax:
  //   OUTPUT [SUMMEDLIST <filename_summed>]
  //     SUMMEDLIST    write summed partials to formatted  file filename_summed
{
  bool sumlist = false;
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NAME)) {
      if (keyIs("SUMM")) {
	sumlist = true;
      } else {
	if (sumlist) {
	  filename_summed = string_value;
	} else {
	  throw SyntaxError(keywords,"unrecognised keyword");
	}
      }
    }
  }
  if (sumlist && filename_summed == "") {
    filename_summed = "SUMMEDLIST";
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
MULTIPLY::MULTIPLY() : CCP4base(), InputBase()
{
  Add_Key("MULTIPLY");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  scale = 1.0;
}
  //--------------------------------------------------------------
Token_value MULTIPLY::parse(std::istringstream& input_stream)
  // MULTIPLY  scale factor to apply on input
  // Syntax:
  //   MULTIPLY <input_scale>
{
  while (get_token(input_stream) != ENDLINE) {
    if (tokenIs(1,NUMBER)) {
      scale = number_value;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
SCORE::SCORE() : CCP4base(), InputBase()
{
  Add_Key("SCORE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  maxmult = 25;
}
//--------------------------------------------------------------
Token_value SCORE::parse(std::istringstream& input_stream)
//
// SCORE  parameters for score functions
//   MAXMULTIPLICITY   maximum multiplicity to use in scoring
//   
// Syntax:
//   SCORE MAXMULTIPLICITY <maximum_multiplicity>  [default 25]
{
  while (get_token(input_stream) != ENDLINE)    {
    if (tokenIs(1,NAME))      {
      if (keyIs("MAXMULTIPLICITY")) maxmult = Nint(get1num(input_stream));
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
ALLOW::ALLOW() : CCP4base(), InputBase()
{
  Add_Key("ALLOW");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);
  allowoutofsequencefiles = false;
}
//--------------------------------------------------------------
Token_value ALLOW::parse(std::istringstream& input_stream)
//
// Syntax:
//  ALLOW  OUTOFSEQUENCEFILES
//        allow files in a series to be out of time sequence
{
  while (get_token(input_stream) != ENDLINE)    {
    if (tokenIs(1,NAME))      {
      if (keyIs("OUTOFSEQ")) {
	allowoutofsequencefiles = true;
      }
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
WAVELENGTH::WAVELENGTH() : CCP4base(), InputBase()
{
  Add_Key("WAVELENGTH");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  wavelength = -1.0;;
}
//--------------------------------------------------------------
Token_value WAVELENGTH::parse(std::istringstream& input_stream)
{
  // Syntax:
  //  WAVELENGTH  <wavelength>

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NUMBER)) {
      wavelength = number_value;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
BLANK::BLANK() : CCP4base(), InputBase()
{
  Add_Key("BLANK");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  nullResolutionfraction = -1.0;
  nullNegativeReject = -1.0;
}
//--------------------------------------------------------------
Token_value BLANK::parse(std::istringstream& input_stream)
{
  // Read parameters for detecting blank images
  // Syntax:
  //  BLANK [NONE] [RESO[FRACTION] <fraction>] [NEGATIVE <negativefraction>]
  //    <fraction>   blank image test will be done out to this fraction
  //                 of the maximum resolution, default 0.5
  //    <negativefraction> fraction of reflections with negative intensity
  //                 above which the image is considered to be blank
  //                 Default 0.3
  // Defaults
  nullResolutionfraction = 0.5; // fraction of maximum resolution to use
  nullNegativeReject = 0.3;  // threshold on proportion of negative reflections

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (keyIs("NONE")) {nullResolutionfraction = -1.0;
      } else {
	if (keyIs("RESO")) nullResolutionfraction = get1num(input_stream);
	if (keyIs("NEGATIVE")) nullNegativeReject = get1num(input_stream);
      }
    }
  }
  if (nullResolutionfraction > 0.0) {
    if (nullResolutionfraction < 0.1 || nullResolutionfraction > 1.0001) {
      throw SyntaxError("ERROR in BLANK command, Resolutionfraction must be between 0.1 and 1.0",
			"");
    }}
  if (nullNegativeReject > 0.0) {
    if (nullNegativeReject < 0.01 || nullNegativeReject > 0.501) {
      throw SyntaxError("ERROR in BLANK command, negative fraction must be between 0.01 and 0.5",
			"");
    }}
  return skip_line(input_stream);
}
//--------------------------------------------------------------
bool WILDFILE::wildfile;
std::string WILDFILE::ftemplate;
//--------------------------------------------------------------
WILDFILE::WILDFILE() : CCP4base(), InputBase()
{
  Add_Key("WILDFILE");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  wildfile = false;
  ftemplate = "*_###.*";
}
//--------------------------------------------------------------
Token_value WILDFILE::parse(std::istringstream& input_stream)
{
  wildfile = true;
  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      ftemplate = string_value;
    }
  }
  return skip_line(input_stream);
}
//--------------------------------------------------------------
POLARISATION::POLARISATION() : CCP4base(), InputBase()
{
  Add_Key("POLARISATION");
  //Add to CCP4base;
  inputPtr iPtr(this);
  possible_fns.push_back(iPtr);

  polarisationfactor = 0.0;
  set = false;
}
//--------------------------------------------------------------
Token_value POLARISATION::parse(std::istringstream& input_stream)
{
  bool xds = true; // by default read the value as for XDS
		   //  FRACTION_OF_POLARIZATION, ie unpolarised = 0.5
  polarisationfactor = 0.99;  // XDS-style default
  set = true;

  while (get_token(input_stream) != ENDLINE)  {
    if (tokenIs(1,NAME)) {
      if (keyIs("XDS")) {
	xds = true;
      } else if (keyIs("MOSFLM")) {
	// "Mosflm" convention (Kahn et al), fraction in range 0 to +1
	// unpolarised = 0.0
	// This convention is used in the rest of the program
	xds = false;
      }
    } else if (tokenIs(1,NUMBER)) {
      polarisationfactor = number_value;
    }
  }
  // We want to store the factor in the range 0 to +1.0
  // so if the value was read in XDS style, convert
  if (xds) {
    if (std::abs(polarisationfactor) < 0.5) {
      polarisationfactor = 0.0; // unpolarised
    } else {
      bool neg = (polarisationfactor < 0.0);
      polarisationfactor = (std::abs(polarisationfactor) - 0.5) * 2.0;
      if (neg) polarisationfactor = -polarisationfactor;
    }
  }
  return skip_line(input_stream);
}
} // phaser_io
