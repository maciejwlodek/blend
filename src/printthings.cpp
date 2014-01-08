// printthings.cpp

// Routines for printing stuff

#include "printthings.hh"
#include "version.hh"
#include "string_util.hh"
#include "hkl_unmerge.hh"

using phaser_io::LOGFILE;
using phaser_io::LXML;
using phaser_io::RESULT;

std::string MkXML(const std::string& text, const std::string& tag)
{
  return "<"+tag+">"+text+"</"+tag+">";
}

//--------------------------------------------------------------
void dump_list(const hkl_unmerge_list& ref_list, phaser_io::Output& output)
// Retrieve all reflections and print them
{
  for (int jref = 0; jref < ref_list.num_reflections(); jref++)
    // loop all reflections jref
    {
      reflection this_refl = ref_list.get_reflection(jref);
      output.logTab(0,LOGFILE, "New reflection: " + this_refl.hkl().format() +
                       " Nobs " +
                       clipper::String(this_refl.num_observations()) + "\n");
      for (int lobs = 0; lobs < this_refl.num_observations(); lobs++)
        {
          observation this_obs = this_refl.get_observation(lobs);
          output.logTab(0,LOGFILE, "New observation: "
                           + this_obs.hkl_original().format() + "  "
                           + " Nparts " +
                           clipper::String(this_obs.num_parts()) + "\n");
          for (int kpart = 0; kpart < this_obs.num_parts(); kpart++)
            {
              observation_part this_part = this_obs.get_part(kpart);
              output.logTab(0,LOGFILE, clipper::String(this_part.I_sigI().I())  + "  " 
                               + clipper::String(this_part.I_sigI().sigI())  + "  " 
                               + clipper::String(this_part.isym()) + "  " 
                               + clipper::String(this_part.batch())  + "\n");
            }
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void PrintTitle( phaser_io::Output& output)
{
  // Make version line
  int sp = 52; // number of spaces in line
  std::string version = PROGRAM_VERSION;
  int n1 = (sp - version.size())/2;
  int n2 = sp - version.size() - n1;
  std::string pad1(n1,' ');
  std::string pad2(n2,' ');
  std::string line = "        *"+pad1+version+pad2+"*\n";

  output.logTabPrintf(0,LOGFILE,
          "\n        ******************************************************\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                     POINTLESS                      *\n");
  output.logTab(0,LOGFILE, line);
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *   Determine Laue group from unmerged intensities   *\n");
  output.logTabPrintf(0,LOGFILE,
          "        *     Phil Evans MRC LMB, Cambridge                  *\n");
  output.logTabPrintf(0,LOGFILE,
          "        * Uses cctbx routines by Ralf Grosse-Kunstleve et al.*\n");
  output.logTabPrintf(0,LOGFILE,
          "        *                                                    *\n");
  output.logTabPrintf(0,LOGFILE,
          "        ******************************************************\n\n");
}
//--------------------------------------------------------------
void PrintSymElementScores(const hkl_unmerge_list& hkl_list,
                           std::vector<SetScores>& scores,
                           std::vector<SCsignificance>& CCsig,
			   const RefListType& isMerged,
                           phaser_io::Output& output)
  //  isMerged  > 0 if input file is merged
  //            < 0 if input file is unmerged
{
  output.logTabPrintf(0,LOGFILE,"\n*******************************************\n");
  output.logTabPrintf(0,LOGFILE,
		       "\nAnalysing rotational symmetry in lattice group %s\n",
		       CCtbxSym::PointGroup(hkl_list.symmetry().symbol_xHM()).RefLGname().c_str());
  output.logTab(0,LOGFILE,"----------------------------------------------");
  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN-->\n");
  if (isMerged == MERGED) {
    // merged file
    output.logTab(0,LOGFILE,
		  "\nInput file is merged so symmetry elements marked '---' are implicitly present\n");
  }
  output.logTabPrintf(0,LOGFILE,"\nScores for each symmetry element\n");
  output.logTabPrintf(0,LOGFILE,"\nNelmt  Lklhd  Z-cc    CC        N  Rmeas    Symmetry & operator (in Lattice Cell)\n\n");

  output.logTab(0,LXML,"<ElementScores>");

  for (int k=0;k<  hkl_list.symmetry().Nelement();k++) {
    // true if this element is implicit in merged symmetry
    bool implicitmerged = CCsig[k].Dummy();
    output.logTab(1,LXML,"<Element>");
    float Zcc  = CCsig[k].Z();
    //      float Zmsd = MSDsig[k].Z();
    float lklhd = CCsig[k].Likelihood();
    int nd = CCsig[k].Nsample();
    if (implicitmerged) nd = 0;
    output.logTabPrintf(0,LOGFILE,"%3d%8.3f%7.2f%7.2f%8d%7.3f",
			k+1, lklhd,
			Zcc, CCsig[k].SC(),
			nd,
			scores[k].OverallRfactor().result().val);
    output.logTabPrintf(1,LXML,
			"<number>%3d </number>\n   <Likelihood>%9.3f</Likelihood>\n   <ZCC>%7.2f </ZCC> <CC>%7.2f </CC> <NCC>%8d </NCC> <R>%7.3f </R>",
			k+1, lklhd,
			Zcc, CCsig[k].SC(),
			CCsig[k].Nsample(),
			scores[k].OverallRfactor().result().val);
    std::string  stars = "     ";
    if (!  hkl_list.symmetry().IsElementIdent(k)) {
      // No stars for identity
      //	  if (Zcc > 3.0) stars = " *   ";
      //	  if (Zcc > 5.0) stars = " **  ";
      //	  if (Zcc > 7.0) stars = " *** ";
      if (lklhd > 0.5) stars = " *   ";
      if (lklhd > 0.7) stars = " **  ";
      if (lklhd > 0.9) stars = " *** ";
      if (implicitmerged) stars = " --- ";  // dummy for merged
    }
    output.logTabPrintf(0,LOGFILE,"%s", stars.c_str());
    output.logTabPrintf(0,LOGFILE,"%s\n",
			hkl_list.symmetry().format_element(k).c_str());
    output.logTabPrintf(1,LXML,"%s\n",
			hkl_list.symmetry().XML_element(k).c_str());
    output.logTab(1,LXML,"</Element>");
  }
  output.logTab(0,LXML,"</ElementScores>");
  output.logTab(0,LOGFILE, "\n<!--SUMMARY_END-->\n");
}
//--------------------------------------------------------------
void PrintSubGroupScores(const std::vector<PGscore>& SGscores,
			 const RefListType& isMerged,
			 const double& overallcompleteness,
			 phaser_io::Output& output)
  // Print results for all subgroups
  // 0.5.1 don't print rms scores, replace by combined element Z-CC
  //
  // On entry:
  //  SGscores   for each subgroup, PG information and scores
  //             ranked by net Z-score
  //  isMerged  MERGED or UNMERGED
  //  overallcompleteness   overall completeness fraction for merged files only
  //
{
  output.logTabPrintf(0,LOGFILE,"\n\nScores for all possible Laue groups which are sub-groups");
  output.logTabPrintf(0,LOGFILE," of lattice group\n");
  output.logTabPrintf(0,LOGFILE,    "--------------------------------------------------------");
  output.logTabPrintf(0,LOGFILE,"-----------------\n");
  output.logTabPrintf(0,LOGFILE,"\nNote that correlation coefficients are from intensities approximately normalised\n");
  output.logTabPrintf(0,LOGFILE,"by resolution, so will be worse than the usual values\n");
  output.logTabPrintf(0,LOGFILE,"Rmeas is the multiplicity weighted R-factor\n");
  output.logTabPrintf(0,LOGFILE,"\nLklhd is a likelihood measure, a probability used in the ranking of space groups\n");
  output.logTabPrintf(0,LOGFILE,"\nZ-scores are from combined scores for");
  output.logTabPrintf(0,LOGFILE," all symmetry elements\nin the sub-group (Z+)");
  output.logTabPrintf(0,LOGFILE," or not in sub-group (Z-)\n");
  output.logTabPrintf(0,LOGFILE,"\n    NetZ = Z+ - Z-\n");
  output.logTabPrintf(0,LOGFILE,"\nNet Z-scores are calculated for correlation coefficients (cc)\n");
  output.logTabPrintf(0,LOGFILE,"The point-group Z-scores Zc are calculated  \n");
  output.logTabPrintf(0,LOGFILE,"    as the Zcc-scores recalculated for all symmetry elements for or against,\n");
  //  output.logTabPrintf(0,LOGFILE,"The point-group Z-scores are calculated either \n");
  //  output.logTabPrintf(0,LOGFILE,"    Zc: as the Zcc-scores recalculated for all symmetry elements for or against,\n");
  //  output.logTabPrintf(0,LOGFILE," or Za: as the RMS Zcc-scores over the symmetry elements for and against\n");
  output.logTabPrintf(0,LOGFILE,
      "\nCC- and R- are the correlation coefficients and R-factors for symmetry elements not in the group\n");
  output.logTabPrintf(0,LOGFILE,"\nDelta is maximum angular difference (degrees) between original cell\n");
  output.logTabPrintf(0,LOGFILE,"and cell with symmetry constraints imposed\n");
  output.logTabPrintf(0,LOGFILE,"\nThe reindex operator converts original index scheme into the conventional\n");
  output.logTabPrintf(0,LOGFILE,"scheme for sub-group\n");
  output.logTabPrintf(0,LOGFILE,
		       "\nAccepted Laue groups are marked '>'\nThe HKLIN Laue group is marked '=' if accepted, '-' if rejected\n");
  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN-->\n");
  if (isMerged == MERGED) {
    // merged file
    output.logTab(0,LOGFILE,
		  "\nInput file is merged so subgroups of the input symmetry are omitted\n");
    output.logTab(1,LOGFILE,
		  "Only groups of the same or higher symmetry are kept");
    if (overallcompleteness < 0.5) {
      output.logTabPrintf(1,LOGFILE,
       "\nWARNING: very low completeness %6.2f in merged file, may not be able to determine Laue group\n",
			  overallcompleteness);
    }
  }

  output.logTabPrintf(0,LOGFILE,
      "\n\n   Laue Group        Lklhd   NetZc  Zc+   Zc-    CC    CC-  Rmeas   R-  Delta ReindexOperator\n\n");

  output.logTab(0,LXML,"<LaueGroupScoreList>");

  for (size_t k=0; k<SGscores.size();k++)
    {
      char accept = ' ';
      std::string AcceptStatus = "NotAccepted";
      if (SGscores[k].Accepted()) 
	{
	  accept = '>';
	  AcceptStatus = "Accepted";
	}
      if (SGscores[k].Original()) {
	if (SGscores[k].Accepted()) {
	  accept = '=';
	  AcceptStatus = "Original";
	} else {
	  accept = '-';
	  AcceptStatus = "OriginalRejected";
	}
      }
      output.logTab(1,LXML,"<LaueGroupScore Accept=\""+AcceptStatus+"\">");

      //      float Zc = SGscores[k].CC_Zgain();
      //float Zc = SGscores[k].CC_Zagain();
      //r      float Zr = SGscores[k].MSD_Zgain();
      std::string  stars = "    ";
      //      if (Zc > 3.0) stars = " *  ";
      //      if (Zc > 5.0) stars = " ** ";
      //      if (Zc > 7.0) stars = " ***";

      double L = SGscores[k].Likelihood();
      if (L > 0.4) stars = " *  ";
      if (L > 0.6) stars = " ** ";
      if (L > 0.8) stars = " ***";


      output.logTabPrintf(0,LOGFILE,"%1c%2d %10s ",
			   accept,  k+1, SGscores[k].RefLGname().c_str());

      // Lattice type
      std::string ct = CrystalType(SGscores[k].crystal_system(),SGscores[k].GetLatType()).format();


      output.logTabPrintf(1,LXML,"<number>%2d</number> <LaueGroupName>%10s</LaueGroupName>\n",
			   k+1, SGscores[k].RefLGname().c_str());
      output.logTabPrintf(1,LXML,"<LatticeType>%2s</LatticeType>\n", ct.c_str());

      output.logTabPrintf(0,LOGFILE,"%4s", stars.c_str());



      output.logTabPrintf(0,LOGFILE,"%7.3f",
			  SGscores[k].Likelihood());
      output.logTabPrintf(1,LXML,"<Likelihood>%9.3f</Likelihood>",
			   SGscores[k].Likelihood());

      output.logTabPrintf(0,LOGFILE," %6.2f%6.2f%6.2f ",
			   SGscores[k].CC_Zgain(),
			   SGscores[k].CC_Zfor(),
			   SGscores[k].CC_Zagainst());
      output.logTabPrintf(1,LXML,"<NetZCC>%6.2f</NetZCC><ZCC_plus>%6.2f</ZCC_plus><ZCC_minus>%6.2f</ZCC_minus>\n",
			   SGscores[k].CC_Zgain(),
			   SGscores[k].CC_Zfor(),
			   SGscores[k].CC_Zagainst());
      /* 
      output.logTabPrintf(0,LOGFILE,"%6.2f%6.2f%6.2f ",
			   SGscores[k].CC_Zagain(),
			   SGscores[k].CC_Zafor(),
			   SGscores[k].CC_Zaagainst());
      output.logTabPrintf(1,LXML,
			   "<NetZCCa>%6.2f</NetZCCa><ZCCa_plus>%6.2f</ZCCa_plus><ZCCa_minus>%6.2f</ZCCa_minus>\n",
			   SGscores[k].CC_Zagain(),
			   SGscores[k].CC_Zafor(),
			   SGscores[k].CC_Zaagainst());
      */
      //r      output.logTabPrintf(0,LOGFILE,"  %6.2f%6.2f%6.2f",
      //r             Zr,
      //r             SGscores[k].MSD_Zfor(),
      //rr             SGscores[k].MSD_Zagainst());

      output.logTabPrintf(0,LOGFILE,"%6.2f%6.2f %6.2f%6.2f",
			   SGscores[k].OverallCC().result().val,
			   SGscores[k].CC_against().SC(),
			   SGscores[k].OverallRfactor().result().val,
			   SGscores[k].Rfactor_against().result().val);

      output.logTabPrintf(0,LXML,"<CC>%6.2f </CC><R>%6.2f </R>",
			   SGscores[k].OverallCC().result().val,
			   SGscores[k].OverallRfactor().result().val);

      output.logTabPrintf(0,LOGFILE," %5.1f",SGscores[k].Delta());
      output.logTabPrintf(0,LXML,"<CellDelta>%5.1f</CellDelta>\n",SGscores[k].Delta());
      
      // Change of basis
      // Ch_b operator for input->LatticeSymmetry
      //      std::vector<double> vop=SGscores[k].RefSGreindex();
      output.logTabPrintf(0,LOGFILE," %s",SGscores[k].RefSGreindexFormat().c_str());
      output.logTab(1,LXML,SGscores[k].RefSGreindex().as_hkl_XML()+"\n");
      output.logTab(1,LXML,SGscores[k].RefSGreindex().as_XML());

      output.logTabPrintf(0,LOGFILE,"\n");
      output.logTab(1,LXML,"</LaueGroupScore>");
    }

  output.logTabPrintf(0,LOGFILE,"\n");
  output.logTab(0,LOGFILE, "\n<!--SUMMARY_END-->\n");
  output.logTab(0,LXML,"</LaueGroupScoreList>");

}
//--------------------------------------------------------------
void PrintSubGroupInfo(const PGscore& SGscores, const bool& PrintHeader,
                       const Chirality& chiral, phaser_io::Output& output)
  // Print various information about subgroup
{
  if (PrintHeader)
    {
      output.logTab(0,LOGFILE, "'New cell' is the cell reindexed into this Laue group\n");
      output.logTab(0,LOGFILE, " in a 'standard' setting, unconstrained to group\n");
      output.logTab(0,LOGFILE, "\nFor 'New cell', 'Deviation' is maximum angular deviation (degrees)\n");
      output.logTab(0,LOGFILE, " from cell with Laue group symmetry imposed\n");
      output.logTab(0,LOGFILE, "\nAlternative indexing schemes which lead to identical or similar\n");
      output.logTab(0,LOGFILE, " cells are grouped on continuation lines if they are symmetry-related\n");
      output.logTab(0,LOGFILE, "'Same cell' and 'Other cell' list alternative indexings:\n");
      output.logTab(0,LOGFILE, " 'Other cell' lists the cell difference defined as mean square\n");
      output.logTab(0,LOGFILE, " difference in basis vectors (columns of orthogonalisation matrix)\n");
    }

  if (SGscores.Accepted())
    {
      output.logTabPrintf(0,LOGFILE,"\n>>>>> %10s ",SGscores.RefLGname().c_str());
      output.logTabPrintf(0,LOGFILE," NetZc =%6.2f", SGscores.CC_Zgain());
      output.logTabPrintf(0,LOGFILE,"  %s\n",SGscores.RefSGreindexFormat().c_str());
      printf
	("               a     b     c  alpha  beta gamma Deviation ReindexOperator\n");
      output.logTabPrintf(0,LOGFILE,"   New cell");
      std::vector<double> tcell = SGscores.TransformedCell();
      for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%6.1f",tcell[i]);
      output.logTabPrintf(0,LOGFILE," %7.2f",SGscores.Delta());
      output.logTabPrintf(0,LOGFILE,"\n");
      SGscores.PrintAlternativeCells(output, false, 3.0, true);
      output.logTabPrintf(0,LOGFILE,"\n");
      SGscores.PrintAllSpacegroups(output, chiral);
    }
}
//--------------------------------------------------------------
void PrintAllSubGroupInfo(const std::vector<PGscore>& SGscores,
		       const ScoreAccept& Accept,
                       const Chirality& chiral, phaser_io::Output& output)
  // Print various information about subgroup
  // Only print groups with score acceptable, & mark these in SGscores
{
  if (Accept.IfSet())
    {
      output.logTabPrintf(0,LOGFILE,
     "\n\nMore information about possible Laue groups with scores above %5.2f\n\n",
			   Accept.Threshold());
    }
  else
    {
      output.logTab(0,LOGFILE, "\n\nMore information about specified Laue group\n\n");
    }
  output.logTab(0,LOGFILE, "'New cell' is the cell reindexed into this Laue group\n");
  output.logTab(0,LOGFILE, " in a 'standard' setting, unconstrained to group\n");
  output.logTab(0,LOGFILE, "\nFor 'New cell', 'Deviation' is maximum angular deviation (degrees)\n");
  output.logTab(0,LOGFILE, " from cell with Laue group symmetry imposed\n");
  output.logTab(0,LOGFILE, "\nAlternative indexing schemes which lead to identical or similar\n");
  output.logTab(0,LOGFILE, " cells are grouped on continuation lines if they are symmetry-related\n");
  output.logTab(0,LOGFILE, "'Same cell' and 'Other cell' list alternative indexings:\n");
  output.logTab(0,LOGFILE, " 'Other cell' lists the cell difference defined as mean square\n");
  output.logTab(0,LOGFILE, " difference in basis vectors (columns of orthogonalisation matrix)\n");

  for (size_t k=0; k<SGscores.size();k++)
    {
      if (SGscores[k].Accepted())
	{
	  output.logTabPrintf(0,LOGFILE,"\n>>>>>%3d %10s ",k+1,SGscores[k].RefLGname().c_str());
	  output.logTabPrintf(0,LOGFILE," NetZc =%6.2f", SGscores[k].CC_Zgain());
	  output.logTabPrintf(0,LOGFILE,"  %s\n",SGscores[k].RefSGreindexFormat().c_str());
	  printf
	    ("               a     b     c  alpha  beta gamma Deviation ReindexOperator\n");
	  output.logTabPrintf(0,LOGFILE,"   New cell");
	  std::vector<double> tcell = SGscores[k].TransformedCell();
	  for (int i=0;i<6;i++) output.logTabPrintf(0,LOGFILE,"%6.1f",tcell[i]);
	  output.logTabPrintf(0,LOGFILE," %7.2f",SGscores[k].Delta());
	  output.logTabPrintf(0,LOGFILE,"\n");
	  SGscores[k].PrintAlternativeCells(output, false, 3.0, true);
	  output.logTabPrintf(0,LOGFILE,"\n");
	  SGscores[k].PrintAllSpacegroups(output, chiral);
	}
    }
}
//--------------------------------------------------------------
void PrintIndexScores
(const std::vector<ReindexScore>& Scores, const bool& result, phaser_io::Output& output)
// Print scores for alternative indexing schemes
// result == true if to be flagged as Result block
{
  //  output.logTab(0,LOGFILE,
  //                 "\nPossible alternative indexing schemes ranked by correlation coefficient (CC)\n");
  //  output.logTab(0,LOGFILE,
  //                 "\nRfactor is calculated from normalised intensities (E^2) so is larger than normal\n");
  //  output.logTab(1,LOGFILE,
  //		 "RMSdeviation is root mean square deviation between base vectors of the\n reference and test cells\n");
  //  output.logTab(0,LOGFILE,
  //		"\n                                                               <!--SUMMARY_BEGIN-->");

  int longop = 30;  // Minimum length
  for (size_t i=0;i<Scores.size();i++)  {
    // Find longest reindex operator
    longop = Max(longop, int(Scores[i].as_hkl().size()));
  }

  phaser_io::outStream OUTSTREAM = LOGFILE;

  if (result) OUTSTREAM = RESULT;
  //*/*/*  if (result) output.logTab(0,LOGFILE,"\n$TEXT:Result:$$ $$\n\n");

  output.logTab(0,OUTSTREAM,
		StringUtil::CentreString("Alternative reindexing", longop)+
    "     Lklhd      CC     R(E^2)    Number Cell_deviation\n");


  output.logTab(0,LXML,"<IndexScores>");
  output.logTabPrintf(1,LXML,
		       "<ScoreCount> %4d</ScoreCount>\n", Scores.size());


  for (size_t i=0;i<Scores.size();i++) {
    output.logTabPrintf(0,OUTSTREAM,"  %s",
    StringUtil::CentreString(Scores[i].as_hkl(),longop).c_str());
    output.logTabPrintf(0,OUTSTREAM, "%8.3f %8.3f %8.3f %9d  %8.2f\n",
			Scores[i].Likelihood(),
			Scores[i].OverallCC().result().val,
			Scores[i].OverallRfactor().result().val,
			Scores[i].OverallCC().result().count,
			Scores[i].Deviation());

    output.logTab(1,LXML,"<Index>");
    output.logTabPrintf(1,LXML,
"  <number>%3d</number> <Likelihood>%7.3f</Likelihood> <CC>%7.2f </CC>  <NCC>%8d </NCC> <R>%7.3f </R>\n",
			i+1,
			Scores[i].Likelihood(),
			Scores[i].OverallCC().result().val,
			Scores[i].OverallCC().result().count,
			Scores[i].OverallRfactor().result().val);
    output.logTabPrintf(1,LXML,"  %s\n",
			Scores[i].as_hkl_XML().c_str());
    output.logTab(1,LXML,"  "+Scores[i].as_XML());
    output.logTab(1,LXML,"</Index>");
  }

  //*/*/*/*  if (result) output.logTab(0,LOGFILE,"$$");
  //  output.logTab(0,LOGFILE,
  //		"\n                                                                 <!--SUMMARY_END-->\n");
  output.logTab(0,LXML,"</IndexScores>");
}
//--------------------------------------------------------------
void PrintAlternativeIndexing
(const std::vector<scala::ReindexOp>& ReidxList, const double& max_diff,
 phaser_io::Output& output)
{
  output.logTab(0,LOGFILE, "\nPossible alternative indexing schemes");
  output.logTab(1,LOGFILE,
                 "Operators labelled \"exact\" are exact by symmetry");
  output.logTab(1,LOGFILE,
                 "For inexact options, deviations are from original cell");
  output.logTab(2,LOGFILE,
                 "(root mean square deviation between base vectors)\n"); 
  output.logTabPrintf(1,LOGFILE,
		 "Maximum accepted RMS deviation between test and reference cells (TOLERANCE) =%5.1f\n\n",
		 max_diff);

  for (size_t i=0;i<ReidxList.size();i++)
    {
      output.logTabPrintf(0,LOGFILE, "%38s", ReidxList[i].as_hkl().c_str());
      if (ReidxList[i].Strict())
        {
          output.logTab(0,LOGFILE, "    exact");
        }
      else
        {
          output.logTabPrintf(0,LOGFILE, "    root mean square deviation %8.3f\n",
                               ReidxList[i].Deviation());
        }
    }
  output.logTab(0,LOGFILE,"");
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void PrintZoneScores(std::vector<Zone>& SZones, 
		     phaser_io::Output& output)
{
  if (SZones.size() <= 0)    {
    output.logTab(0,LOGFILE,"No systematic absence zones to test in this lattice group\n");
    return;
  }
  output.logTab(0,LOGFILE,
		 "\nEach 'zone' (axis or plane) in which some reflections may be systematically absent\n");
  output.logTab(0,LOGFILE,
		 "are scored by Fourier analysis of I'/sigma(I). 'PeakHeight' is the value\n");
  output.logTab(0,LOGFILE,
		 "in Fourier space at the relevent point (eg at 1/2 for a 2(1) axis)\n");
  output.logTab(0,LOGFILE,
		 "relative to the origin. This has an ideal value of 1.0 if the corresponding\n");
  output.logTab(0,LOGFILE,
     "symmetry element is present. Zone directions (a,b,c) shown here are in the\nlattice group frame\n");
  output.logTab(0,LOGFILE,
		 "\n'Probability' is an estimate of how likely the element is to be  present\n");
  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN-->\n");
  output.logTab(0,LOGFILE,
 "\n         Zone                Number PeakHeight  SD  Probability  ReflectionCondition");
	
  output.logTab(0,LXML,"<ZoneScoreList>");

  hkl_symmetry symm;
  bool prunedData = false;

  for (size_t iz=0;iz<SZones.size();iz++) {
    if (symm != SZones[iz].LGsymm()) {
      symm = SZones[iz].LGsymm();
      output.logTab(0,LOGFILE,
		    "\nZones for Laue group "+symm.symbol_xHM());
    }
    std::string zoneType = "screw axis";
    if (!SZones[iz].Axis()) {zoneType = "glide plane";}
    if (SZones[iz].Valid()) {
      std::vector<double> P = SZones[iz].p();
      for (int i=1;i<SZones[iz].NgridPoints();i++) { // Omit origin
	if (SZones[iz].ValidPoint(i)) {
	  if (!SZones[iz].PrunedData(i)) {
	  std::string  stars = "   ";
	  if (P[i] > 0.7) stars = "  *";
	  if (P[i] > 0.8) stars = " **";
	  if (P[i] > 0.9) stars = "***";
	  if (P[i] >= 0.0) {
	    output.logTabPrintf(0,LOGFILE,
				"%2d%24s %8d%8.3f%8.3f   %3s%6.3f   %s\n", iz+1,
				SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
				SZones[iz].Nobs(),
				SZones[iz].FourierVal()[i],
				SZones[iz].GetSD(i),
				stars.c_str(),
				P[i], 
				SZones[iz].FormatConditionLGFrame(i).c_str());
	  } else {
	    output.logTabPrintf(0,LOGFILE,
				"%2d%24s %8d%8.3f%8.3f   %3s  Null   %s\n", iz+1,
			      SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
			      SZones[iz].Nobs(),
			      SZones[iz].FourierVal()[i],
			      SZones[iz].GetSD(i),
			      stars.c_str(),
			      SZones[iz].FormatConditionLGFrame(i).c_str());
	  }
	  output.logTab(1,LXML,"<Zone>");
	  output.logTabPrintf(1,LXML,
			      "<Number>%2d</Number><ZoneType>%22s</ZoneType><Nobs>%8d</Nobs><PeakHeight>%8.3f</PeakHeight><SDPkHt>%8.3f</SDPkHt><Prob>%8.3f</Prob><Condition>%s</Condition>\n",
			      iz+1,
			      SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
			      SZones[iz].Nobs(),
			      SZones[iz].FourierVal()[i],
			      SZones[iz].GetSD(i),
			      P[i], 
			      SZones[iz].FormatConditionLGFrame(i).c_str());
	  output.logTab(1,LXML,"</Zone>");
	} else {
	    // ValidPoint(i) && PrunedData(i)
	    output.logTabPrintf(0,LOGFILE,
				"%2d%24s    Systematically missing data, %s may be present\n", iz+1,
				SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
				zoneType.c_str());
	    output.logTab(1,LXML,"<Zone>");
	    output.logTabPrintf(1,LXML,
				"<Number>%2d</Number><ZoneType>%22s</ZoneType><Nobs>%8d</Nobs><PeakHeight>%8.3f</PeakHeight><SDPkHt>%8.3f</SDPkHt><Prob>%8.3f</Prob><Condition>%s</Condition>\n",
				iz+1,
				SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
				0, 0.0, 0.0, 0.0,
				SZones[iz].FormatConditionLGFrame(i).c_str());
	    output.logTab(1,LXML,"</Zone>");
	    prunedData = true;
	  }
	} else {
	  // Invalid point
	}
      }  // end loop points
    } else {
      // Invalid zone
      for (int i=1;i<SZones[iz].NgridPoints();i++) { // Omit origin
	output.logTabPrintf(0,LOGFILE,
			    "%2d%24s    No observations\n", iz+1,
			    SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str());
	output.logTab(1,LXML,"<Zone>");
	output.logTabPrintf(1,LXML,
			    "<Number>%2d</Number><ZoneType>%22s</ZoneType><Nobs>%8d</Nobs><PeakHeight>%8.3f</PeakHeight><SDPkHt>%8.3f</SDPkHt><Prob>%8.3f</Prob><Condition>%s</Condition>\n",
			    iz+1,
			    SZones[iz].formatLGFrame(SZones[iz].Ngrid()[i]).c_str(),
			    0, 0.0, 0.0, 0.0,
			    SZones[iz].FormatConditionLGFrame(i).c_str());
	output.logTab(1,LXML,"</Zone>");
      }
    }
  }  // end loop zones
  if (prunedData) {
    output.logTab(0,LOGFILE,
        "\n\nWARNING! one or more zones have data systematically missing, see reflection list,");
    output.logTab(0,LOGFILE,
      "   thus we cannot determine if the reflections are truly systematically absent\n\n");
    output.logTab(0,RESULT,
        "\n\nWARNING! one or more zones have data systematically missing from the input file");
    output.logTab(0,RESULT,
      "  thus we cannot determine if reflections are truly systematically absent\n\n");
  }

  output.logTab(0,LOGFILE, "\n<!--SUMMARY_END-->\n");
  output.logTab(0,LXML,"</ZoneScoreList>");
}
//--------------------------------------------------------------
void PrintPossibleSpaceGroups(const std::vector<scala::PossibleSpaceGroup>& groups,
			      const double& Pcutoff,
			      const Chirality& chiral,
			      const int& settingOption,
			      phaser_io::Output& output)
// On entry:
//  settingOption     = 0 CELL-BASED = +1 SYMMETRY-BASED = -1 C2
//                    <= 0 to keep primitive orthorhombic spacegroups in the Lauegroup setting
//                              ie a < b < c
//                    == 0 to allow C2 as possibly I2
{			      
  // Maximum score
  double Pmax = groups[0].SysAbsProb();
  double Plimit = Pmax * Pcutoff;
  int naccept = 0;
  for (size_t i=0;i<groups.size();i++)
    {
      if (groups[i].SysAbsProb() > Plimit) naccept++;
    }
  if (naccept == 0)
    {
      output.logTab(0,LOGFILE,
        "\n\nWARNING! WARNING! WARNING! Insufficient information to choose spacegroup\n");
      output.logTab(0,LOGFILE,"Printing all possible groups\n\n");
      Plimit = -10000.;
    }

  output.logTab(0,LOGFILE,"\n\nPossible spacegroups:");
  output.logTab(0,LOGFILE,  "--------------------\n");
  output.logTab(0,LOGFILE,"Indistinguishable space groups are grouped together on successive lines\n");
  output.logTab(0,LOGFILE,
		 "\n'Reindex' is the operator to convert from the input hklin frame to the standard spacegroup frame.\n");
  if (chiral == CHIRAL) 
    output.logTab(0,LOGFILE,
		   "\n'TotProb' is a total probability estimate (unnormalised)\n");
  output.logTab(0,LOGFILE,
		 "\n'SysAbsProb' is an estimate of the probability of the space group based on\n");
  output.logTab(0,LOGFILE,
		 "the observed systematic absences.\n");
  output.logTab(0,LOGFILE,
		 "\n'Conditions' are the reflection conditions (absences)\n");


  if (chiral != CHIRAL) 
    {
      output.logTab(0,LOGFILE,
		     "'TotProb' is a total probability estimate (unnormalised) including the probability\n");
      output.logTab(0,LOGFILE,
		     "of the crystal being centrosymmetric from the <|E^2-1|> statistic.\n");
      output.logTab(0,LOGFILE,
		     "Chiral space groups are marked '*' and centrosymmetric ones 'O'\n");
    }

  output.logTab(0,LOGFILE,"\n\n   Spacegroup         TotProb SysAbsProb     Reindex         Conditions\n\n");

  double prob_last = -1000000.;

  output.logTab(0,LXML,"\n<SpacegroupList>");

  bool PrimitiveOrthSymmetryBased = false;
  bool C2setting = false;

  for (size_t i=0;i<groups.size();i++)
    {
      if (groups[i].SysAbsProb() > Plimit)
	{
	  // set flag true for keeping Laue group setting for P, I, F orthorhombic groups, or I2
	  bool LGsetting = true;  // conventional setting (a<b<c) and allow I2
	  if (groups[i].IsPIForthorhombic() && settingOption > 0) {
	    LGsetting = false; // oP && SYMMETRY-BASED
	    PrimitiveOrthSymmetryBased = true;
	  }
	  if (groups[i].IsI2() && settingOption != 0) {
	    LGsetting = false;  // mC/I && C2
	    C2setting = true;
	  }

	  output.logTab(1,LXML,"\n<Spacegroup>");
	  if (i > 0 && std::abs(prob_last - groups[i].Prob()) > 0.001)
	    output.logTabPrintf(0,LOGFILE,"    ..........\n");
	  prob_last = groups[i].Prob();

	  //	  std::string sname = "<"+groups[i].Name()+">";
	  std::string rname = "<"+groups[i].Name(LGsetting)+">";
	  std::string ch = " ";
	  if (groups[i].Chiral() == CHIRAL && chiral != CHIRAL) ch = "*";
	  if (groups[i].Chiral() == CENTROSYMMETRIC) ch = "O";

	  output.logTabPrintf(0,LOGFILE,"%14s (%3d) %1s",
			       rname.c_str(),
			       groups[i].Sgnumber(),
			       ch.c_str());
	  output.logTabPrintf(1,LXML,
			       "<SpacegroupName>%14s</SpacegroupName><SGnumber>%3d</SGnumber>",
			       groups[i].Name(LGsetting).c_str(),
			       groups[i].Sgnumber());
	  if (chiral != CHIRAL)
	    {output.logTab(1,LXML,Chiral_as_string(groups[i].Chiral()));}

	  // total reindex operator H =
	  //     H(original->Lauegroup) * H(LaueGroup->SpaceGroup)
	  ReindexOp TotalReindex = groups[i].Reindex(LGsetting);
	  
	  output.logTabPrintf(0,LOGFILE,"%7.3f%7.3f",
			       groups[i].Prob(),
			       groups[i].SysAbsProb());
	  output.logTabPrintf(1,LXML,
			       "<TotalProb>%7.3f</TotalProb><SysAbsProb>%7.3f</SysAbsProb>",
			       groups[i].Prob(),
			       groups[i].SysAbsProb());

	  if (TotalReindex.IsIdentity())
	    {output.logTabPrintf(0,LOGFILE,"                        ");}
	  else
	    {output.logTabPrintf(0,LOGFILE,"%24s",
		 StringUtil::CentreString(TotalReindex.as_hkl(),24).c_str());}

	  output.logTab(1,LXML, "\n   "+TotalReindex.as_hkl_XML());
	  output.logTab(1,LXML, "   "+TotalReindex.as_XML());

	  output.logTabPrintf(0,LOGFILE," %s",
			       groups[i].Condition(LGsetting).c_str());
	  output.logTabPrintf(1,LXML,"<Condition>%s</Condition>",
			       groups[i].Condition(LGsetting).c_str());

	  std::vector<int> ZoneList = groups[i].ZoneList();  // list of zone numbers
	  if (ZoneList.size() > 0)
	    {
	      std::string zl;
	      for (size_t iz=0;iz<ZoneList.size();iz++)
		{
		  zl += clipper::String(ZoneList[iz]);
		  if (iz < ZoneList.size()-1) zl += ",";
		}
	      if (ZoneList.size() == 1) {
		output.logTabPrintf(0,LOGFILE,
			 " (zone %s)\n", StringUtil::Strip(zl).c_str());
	      } else {
		output.logTabPrintf(0,LOGFILE," (zones %s)\n",
				    StringUtil::Strip(zl).c_str());
	      }
	      output.logTabPrintf(1,LXML,"<ZoneNumbers>%s</ZoneNumbers>\n",
				  StringUtil::Strip(zl).c_str());
	    } else {
	    output.logTabPrintf(0,LOGFILE,"\n");
	  }
	  //^
	  //	  std::cout << "Reindex, total, orig " << TotalReindex.as_hkl() << " "
	  //		    << groups[i].ReindexOrig().as_hkl() << "\n";
	  //^-
	  //^!
	  //	  std::cout << "Original setting " << groups[i].Name() << " "
	  //		    << groups[i].ReindexOrig().as_hkl()
	  //		    << " " << groups[i].ConditionLG() << "\n" ;
	  //^!-

	  output.logTab(1,LXML,"</Spacegroup>");
	}
    }
  output.logTab(0,LXML,"</SpacegroupList>\n");
  if (PrimitiveOrthSymmetryBased) {
    output.logTab(0,LOGFILE,
		  "\nPrimitive orthorhombic space groups have been put in the 'reference' (SYMMETRY-BASED) setting");
  }
  if (C2setting) {
    output.logTab(0,LOGFILE,
		  "\nCentred monoclinic space groups have been kept as C-centred");
  }

}
//--------------------------------------------------------------
void OutputZoneTable(const Zone& SZone,
		     phaser_io::Output& output)
// Write out zone data as table suitable for loggraph
// Averaged values already stored in Zone objects
{
  if (SZone.AverageIsigI().size())
    {
      if (SZone.Axis())
	// Axis
	{
	  output.logTab(0,LOGFILE,"\n$TABLE: Axial reflections, axis "+
			 SZone.Direction()+" (lattice frame):\n"+
			 "$GRAPHS: I/sigI vs. index :N: 1,4,5 :"+
			 ": I vs. index :N: 1,2 :$$\n\n"+
			 "    Index         I       sigI     I/sigI   I'/sigI  $$ $$\n\n");
	  output.logTab(0,LXML,"\n<AxialReflections latticeaxis=\""+
			 SZone.Direction()+"\">");

	  //	  int n = 0;

	  ASSERT (SZone.AverageIsigI().size() == SZone.AdjustedIsigI().size());

	  for (size_t i=0;i<SZone.AdjustedIsigI().size();i++)
	    {
	      IsigI Isig = SZone.AverageIsigI()[i];
	      IsigI IsigAdj = SZone.AdjustedIsigI()[i];
	      if (Isig.sigI() > 0.0) {
		output.logTabPrintf(0,LOGFILE,
				    "%8d %10.0f %10.0f %9.2f %9.2f\n",
				    i, Isig.I(), Isig.sigI(), 
				    Isig.I()/Isig.sigI(), IsigAdj.I()/IsigAdj.sigI());
		output.logTabPrintf(1,LXML,
				    "<AxialObservation>\n      <index>%8d</index><I>%10.0f</I><sigI>%10.0f</sigI><I_ov_sigI>%9.2f</I_ov_sigI><Iadj_ov_sigI>%9.2f</Iadj_ov_sigI></AxialObservation>\n",
				    i, Isig.I(), Isig.sigI(), 
				    Isig.I()/Isig.sigI(), IsigAdj.I()/IsigAdj.sigI());
	      }
	    }
	  output.logTab(0,LOGFILE, "$$\n");
	  output.logTab(0,LXML,"</AxialReflections>");
	}
    }
}
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
			phaser_io::Output& output)
{
  std::vector<int> runOffsets;  // Dummy
  PrintFileInfoToXML(StreamName, FileName, cell, SpaceGroupName, runOffsets, output);
}
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
			const std::vector<int>& runOffsets,
			phaser_io::Output& output)
{
  if (output.doXmlout())
    {
      output.logTab(0,LXML,
		     "<ReflectionFile stream=\""+StreamName+
		     "\" name=\""+FileName+"\">\n");
      output.logTab(0,LXML,
		     cell.xml());
      output.logTab(0,LXML,
		     "<SpacegroupName> "+SpaceGroupName+"</SpacegroupName>");
      if (runOffsets.size() > 0) {
	output.logTab(1,LXML,"<BatchNumberIncrements>");
	for (size_t i=0;i<runOffsets.size();i++) {
	  output.logTabPrintf(2,LXML,"<Run> %3d </Run>\n", i+1);
	  output.logTabPrintf(2,LXML,"<BatchOffset> %8d </BatchOffset>\n",
			      runOffsets[i]);
	}
	output.logTab(1,LXML,"</BatchNumberIncrements>");
      }
      output.logTab(0,LXML,
		    "</ReflectionFile>");
    }
}
//--------------------------------------------------------------
void PrintNormalisationResult(const Normalise& NormRes,
			      const ResoRange& ResRange, const double& MinIsigRatio,
			      phaser_io::Output& output)
{
  std::vector<BfactorModel> BFcorr = NormRes.BfactorCorr();  // Correction factors
  int Nruns = BFcorr.size();
  if (Nruns > 1) {
    output.logTab(0,LOGFILE,"\nIntensity normalisation for each run:");
    for (int ir=0;ir<Nruns;ir++) {
      std::vector<double> p = BFcorr[ir].Params();
      if (BFcorr[ir].Nparam() > 2) {
	output.logTabPrintf(1,LOGFILE,
			    "Run %3d: B-factor = %6.1f  + %8.4f * time  (final B %5.1f)\n",
			    ir+1, p[1], p[2], p[1] + p[2]*BFcorr[ir].DTime());
      }
      else {
	output.logTabPrintf(1,LOGFILE,
			    "Run %3d: B-factor = %6.1f\n", ir+1, p[1]);
      }
    }
  }
  else {
    std::vector<double> p = BFcorr[0].Params();
    if (BFcorr[0].Nparam() > 2) {
      output.logTabPrintf(1,LOGFILE,
			  "\nIntensity normalisation: B-factor = %6.1f  + %8.4f * time  (final B %5.1f)\n",
			  p[1], p[2],  p[1] + p[2]*BFcorr[0].DTime());
    }
    else {
      output.logTabPrintf(1,LOGFILE,
			  "\nIntensity normalisation: B-factor = %6.1f\n", p[1]);
    }
  }
  if (MinIsigRatio > 0.0) {
    output.logTabPrintf(0,LOGFILE, "\nResolution range reset to %8.2f to %8.2f\n",
			ResRange.ResLow(), ResRange.ResHigh());
    output.logTabPrintf(1,LOGFILE, "using I/sigmaI cutoff %6.1f\n", MinIsigRatio);
  }
}
//--------------------------------------------------------------
void PrintReindexSummary(const RefListType& reflisttype,
			 const IO_files& AllFiles,
			 phaser_io::Output& output)
{
  int if1;
  //*/  output.logTab(0,LOGFILE,"\n<!--SUMMARY_BEGIN-->\n");
  std::string outstring;
  
  if (reflisttype == NONE) {
    outstring += FormatOutput::logTab(0,"\nAlternative indexing relative to first file(s):");
    if1 = 1;
  } else {
    outstring += FormatOutput::logTab(0,"\nAlternative indexing relative to reference file "+
		  BaseFileName(AllFiles.Filename("HKLREF"),false));
    if1 = 0;
  }
  int maxlength = -1;
  for (int i=if1;i<AllFiles.NumberOfHKLINnames();i++) {
    if (AllFiles.IsHKLINreindex(i)) {
      maxlength = Max(maxlength,
		      int(AllFiles.HKLINreindex(i).as_hkl().length()));
    }
  }    // maximum length of reindex operator
  maxlength = Max(20, maxlength+4);

  outstring += FormatOutput::logTab(0,
		      "\n    "+StringUtil::CentreString
		      ("Reindex operator",maxlength)+
		      "       CC      Lklhd  Confidence     File name");
  std::string xmlstring = "<BestReindex>\n";

  ReflectionFileType type;
  PxdName pxd;
  int cs = -1;
  for (int i=if1;i<AllFiles.NumberOfHKLINnames();i++) {
    if (cs < 1 || cs != AllFiles.HKLIN_fileSeries(i)) {
      cs = AllFiles.HKLIN_fileSeries(i);
      if (AllFiles.IsHKLINreindex(i)) {
	outstring += FormatOutput::logTabPrintf(0, "%3d %s", i+1,
	  StringUtil::CentreString(AllFiles.HKLINreindex(i).as_hkl(),maxlength).c_str());
	outstring += FormatOutput::logTabPrintf(0, "  %8.3f",
			    AllFiles.HKLIN_CC(i).result().val);
	outstring += FormatOutput::logTabPrintf(0, "  %8.3f  %8.3f",
			AllFiles.HKLIN_Data(i).likelihood,
			AllFiles.HKLIN_Data(i).confidence);
	outstring += FormatOutput::logTabPrintf(0, "    %s\n",
			BaseFileName(AllFiles.HKLINfilename(i), false).c_str());
	// XML
	xmlstring += "   "+StringUtil::MakeXMLtag("Filenumber",StringUtil::itos(i+1,3))+"\n";
	xmlstring += "   "+AllFiles.HKLINreindex(i).as_hkl_XML()+"\n";
	xmlstring += "   "+AllFiles.HKLINreindex(i).as_XML()+"\n";
	xmlstring += "   "+
	  StringUtil::MakeXMLtag("CC",
				StringUtil::ftos(AllFiles.HKLIN_CC(i).result().val,8,3))+"\n";
	xmlstring += "   "+
	  StringUtil::MakeXMLtag("Likelihood",
				 StringUtil::ftos(AllFiles.HKLIN_Data(i).likelihood,8,3))+"\n";
	xmlstring += "   "+
	  StringUtil::MakeXMLtag("Confidence",
				 StringUtil::ftos(AllFiles.HKLIN_Data(i).confidence,8,3))+"\n";
	xmlstring += "   "+
	  StringUtil::MakeXMLtag("Filename",BaseFileName(AllFiles.HKLINfilename(i), false))+"\n";
      }
    }
  }
  output.logTab(0,RESULT,outstring+"\n");

  xmlstring += "</BestReindex>\n";
  output.logTab(0,LXML, xmlstring);
  //*/  output.logTab(0,LOGFILE,outstring);
  //*/  output.logTab(0,LOGFILE, "\n<!--SUMMARY_END-->\n");
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void PrintUnmergedHeaderStuff(const scala::hkl_unmerge_list& hkl_list,
			      phaser_io::Output& output,
			      const int& verbose)
  // Optional summary printing
  {
    if (verbose > 0)  {
      if (verbose == 1 || !hkl_list.IsReady()) {
	output.logTab(0,LOGFILE,
		      "      ResolutionRange    NobsParts  Nbatches  Ndatasets");
	output.logTabPrintf(0,LOGFILE,
			    "     %8.2f %6.2f  %10d%7d%10d\n",
			    hkl_list.ResRange().ResLow(), hkl_list.ResRange().ResHigh(),
			    hkl_list.num_parts(), hkl_list.num_batches(),
			    hkl_list.num_datasets());
	
      } else {
	//        output.logTabPrintf(0,LOGFILE,
	//			    "\nSummary of reflection list\n");
	output.logTabPrintf(0,LOGFILE,
			    "   Resolution range accepted: %8.2f    %8.2f\n",
			    hkl_list.ResRange().ResLow(),
			    hkl_list.ResRange().ResHigh());
	
	output.logTabPrintf(0,LOGFILE,
			    "\n   Number of reflections      =    %10d\n",
			    hkl_list.num_reflections_valid());
	output.logTabPrintf(0,LOGFILE,
			    "   Number of observations     =    %10d\n",
			    hkl_list.num_observations());
	output.logTabPrintf(0,LOGFILE,
			    "   Number of parts            =    %10d\n",
			    hkl_list.num_parts());
	output.logTabPrintf(0,LOGFILE,
			    "   Number of batches in file  =    %10d\n",
			    hkl_list.num_batches());
	if (hkl_list.num_accepted_batches() != hkl_list.num_batches()) {
	  output.logTabPrintf(0,LOGFILE,
			      "   Number of accepted batches =    %10d\n",
			      hkl_list.num_accepted_batches());
	}
	output.logTabPrintf(0,LOGFILE,
			    "   Number of datasets         =    %10d\n",
			    hkl_list.num_datasets());
      }
      if (verbose <= 3) {
	int ndatasets = hkl_list.num_datasets();
	std::vector<Batch> batches = hkl_list.Batches();
	int nbatches = batches.size();
	// 1st, last batch
	std::vector<std::pair<int,int> >  rejectedbatches;
	// Make list of rejected batch ranges
	int b1 = -1;
	for (int ib=0;ib<nbatches;ib++) {
	  if (!batches[ib].Accepted()) {
	    if (b1 < 0) {b1 = ib;}
	  } else {
	    if (b1 >= 0) {
	      // store range
	      if (b1 == ib-1) {
		// just one
		rejectedbatches.push_back
		  (std::pair<int,int>(batches[b1].num(),0));
	      } else {
		rejectedbatches.push_back
		  (std::pair<int,int>(batches[b1].num(),batches[ib-1].num()));
	      }
	      b1 = -1;
	    }
	  }
	}
	int nrb = rejectedbatches.size();
	std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
	std::vector<Run> runlist = hkl_list.RunList();
	for (int k=0; k<ndatasets; k++) {  // loop datasets
	  output.logTab(0,LOGFILE,datasets[k].pxdname().formatPrint());
	  if (verbose == 3) {
	    for (size_t i=0;i<runlist.size();i++) {
	      if (runlist[i].DatasetIndex() == k) {
		output.logTab(0,LOGFILE,runlist[i].formatPrintBrief(datasets));
		std::string rejlist;
		for (int k=0;k<nrb;k++) { 
		  if (runlist[i].IsInList(rejectedbatches[k].first)) {
		    // Rejected batches in this run
		    if (rejlist.size() > 0) rejlist += ", ";
		    if (rejectedbatches[k].second == 0) {
		      rejlist +=
		StringUtil::Strip(clipper::String(rejectedbatches[k].first,6));
		    } else {
		      rejlist +=
	      	StringUtil::Strip(clipper::String(rejectedbatches[k].first,6)+
	       	      "-"+clipper::String(rejectedbatches[k].second,6));
		    }
		  }
		}
		if (rejlist.size() > 0) {
		  output.logTab(3,LOGFILE,"Excluded batches: "+rejlist);}
	      }
	    }
	    // Check for different unit cells within dataset
	    if (datasets[k].NumberofCells() > 1) {
	      output.logTab(3,LOGFILE,
		       "Dataset contains multiple unit cells and wavelengths\n");
	      output.logTab(0,LOGFILE,datasets[k].formatAllCells());
	      if (datasets[k].WorstDeviation() >  0.5*hkl_list.ResRange().ResHigh()) {
		output.logTabPrintf(0,LOGFILE,
    "\n**** WARNING: cells may have unacceptable differences at resolution %6.2fA\n",
				    hkl_list.ResRange().ResHigh());
		output.logTabPrintf(0,LOGFILE,
				    "****  worst deviation %7.3fA\n\n",
				    datasets[k].WorstDeviation());
		output.logTabPrintf(0,RESULT,
    "\n\nWARNING! in dataset %s\n different cells may have unacceptable differences at the maximum resolution %6.2fA\n",
				    datasets[k].pxdname().format().c_str(),
				    hkl_list.ResRange().ResHigh());
		output.logTabPrintf(0,RESULT,
				    "   worst deviation %7.3fA\n\n",
				    datasets[k].WorstDeviation());
	      }
	    }
	  } // verbose == 3
	}   // end loop datasets
      } else if (verbose > 3) {
	std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
	for (size_t k=0; k<datasets.size(); k++) {
	  output.logTab(0,LOGFILE,
			datasets[k].formatPrint());
	}
      }
      output.logTabPrintf(1,LOGFILE, "Average unit cell: ");
      output.logTab(0,LOGFILE,hkl_list.Cell().formatPrint());
      output.logTabPrintf(1,LOGFILE,"");
      
      if (verbose > 3) {
	std::vector<Run> runlist = hkl_list.RunList();
	std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
	for (size_t i=0;i<runlist.size();i++) {
	  output.logTab(0,LOGFILE,runlist[i].formatPrint(datasets));}
      }
    }

    // XML things
    if (verbose > 2) {
      std::string hklstream = "HKLIN";
      // Is there more than one unique file number?
      std::vector<Run> runlist = hkl_list.RunList();
      int fn = runlist[0].FileNumber();
      bool OneFile = true;
      if (runlist.size() > 1) {
	for (size_t i=1;i<runlist.size();i++) {
	  if (runlist[i].FileNumber() != fn) {OneFile = false;}
	}
      }
      output.logTab(0,LXML,"<ReflectionData>");
      output.logTabPrintf(1,LXML,
			  "<NumberReflections>  %10d </NumberReflections>\n",
			  hkl_list.num_reflections_valid());
      output.logTabPrintf(1,LXML,
			  "<NumberObservations> %10d </NumberObservations>\n",
			  hkl_list.num_observations());
      output.logTabPrintf(1,LXML,
			  "<NumberParts>        %10d </NumberParts>\n",
			  hkl_list.num_parts());
      output.logTabPrintf(1,LXML,
			  "<NumberBatches>      %10d </NumberBatches>\n", 
			  hkl_list.num_batches());
      output.logTabPrintf(1,LXML,
			  "<NumberDatasets>     %10d </NumberDatasets>\n",
			  hkl_list.num_datasets());

      int ndatasets = hkl_list.num_datasets();
      std::vector<Xdataset> datasets = hkl_list.AllXdatasets();
      for (int k=0; k<ndatasets; k++) {
	output.logTabPrintf(1,LXML, "<Dataset  name=\"%s\">\n",
			    datasets[k].pxdname().format().c_str());
	for (size_t i=0;i<runlist.size();i++) {
	  if (runlist[i].DatasetIndex() == k) {
	    output.logTabPrintf(2,LXML,"<Run> <number> %3d </number>\n",i+1);
	    output.logTabPrintf(2,LXML,
				"<BatchRange> %8d %8d </BatchRange>\n",
				runlist[i].BatchRange().first, runlist[i].BatchRange().second);
	    output.logTabPrintf(2,LXML,"<BatchOffset> %8d </BatchOffset>\n",
				runlist[i].BatchNumberOffset());
	    if (!OneFile) {
	      hklstream =
		StringUtil::Strip("HKLIN"+clipper::String(runlist[i].FileNumber()));
	    }
	    output.logTab(2,LXML,"<FileStream> "+hklstream+" </FileStream>");
	    output.logTabPrintf(2,LXML,"</Run>\n",i+1);
	  }
	}
	output.logTabPrintf(1,LXML, "</Dataset>\n");
      }
      output.logTab(0,LXML,"</ReflectionData>");
    }
  }
//--------------------------------------------------------------
  //--------------------------------------------------------------
  void PrintError(const std::string& ErrorMessage,
		  phaser_io::Output& output)
  {
    output.logTab(0,LOGFILE, "\n**** ERROR ****\n");
    output.logTab(1,LOGFILE, ErrorMessage);
    output.logTab(0,LOGFILE,   "**** ERROR ****\n\n");
  }
  //--------------------------------------------------------------
  void PrintWarning(const std::string& ErrorMessage,
		  phaser_io::Output& output)
  {
    output.logTab(0,LOGFILE, "\n**** WARNING ****\n");
    output.logTab(1,LOGFILE, ErrorMessage);
    output.logTab(0,LOGFILE,   "**** WARNING ****\n\n");
  }
