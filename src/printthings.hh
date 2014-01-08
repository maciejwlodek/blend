// printthings.hh

#ifndef PRINTTHINGS_HEADER
#define PRINTTHINGS_HEADER

// Routines for printing stuff

#include "scala.hh"
#include "io_files.hh"
#include "readallhklinfiles.hh"

using namespace scala;

//--------------------------------------------------------------
//--------------------------------------------------------------
void PrintXML(const std::string text,
	      phaser_io::Output& output);
//--------------------------------------------------------------
void dump_list(const hkl_unmerge_list& ref_list, phaser_io::Output& output);
// Retrieve all reflections and print them
//--------------------------------------------------------------
void PrintTitle( phaser_io::Output& output);
//--------------------------------------------------------------
//  On entry: 
//  isMerged  MERGED or UNMERGED
void PrintSymElementScores(const hkl_unmerge_list& hkl_list,
                           std::vector<SetScores>& scores,
                           std::vector<SCsignificance>& CCsig,
			   const RefListType& isMerged,
                           phaser_io::Output& output);
//--------------------------------------------------------------
void PrintSubGroupScores(const std::vector<PGscore>& SGscores, 
			 const RefListType& isMerged,
			 const double& overallcompleteness,
			 phaser_io::Output& output);
  // Print results for all subgroups
  // 0.5.1 don't print rms scores, replace by combined element Z-CC
  //
  // On entry:
  //  SGscores   for each subgroup, PG information and scores
  //             ranked by net Z-score
  //  isMerged  MERGED or UNMERGED
  //  overallcompleteness   overall completeness fraction for merged files only
  //
//--------------------------------------------------------------
void PrintSubGroupInfo(const PGscore& SGscores, const bool& PrintHeader,
                       const Chirality& chiral, phaser_io::Output& output);
// Print various information about subgroup
//--------------------------------------------------------------
void PrintAllSubGroupInfo(const std::vector<PGscore>& SGscores,
		       const ScoreAccept& Accept,
                       const Chirality& chiral, phaser_io::Output& output);
  // Print various information about subgroup
  // Only print groups with score acceptable
//--------------------------------------------------------------
void PrintIndexScores
(const std::vector<ReindexScore>& Scores,  const bool& result,
 phaser_io::Output& output);
  // Print scores for alternative indexing schemes
//--------------------------------------------------------------
void PrintAlternativeIndexing
(const std::vector<scala::ReindexOp>& ReidxList, const double& max_diff,
 phaser_io::Output& output);
//--------------------------------------------------------------
void PrintZoneScores(std::vector<Zone>& SZones, 
		     phaser_io::Output& output);
//--------------------------------------------------------------
void PrintPossibleSpaceGroups(const std::vector<scala::PossibleSpaceGroup>& groups,
			      const double& Pcutoff,
			      const Chirality& chiral,
			      const int& settingOption,
			      phaser_io::Output& output);
//--------------------------------------------------------------
void OutputZoneTable(const Zone& SZone,
		     phaser_io::Output& output);
// Write out zone data as table suitable for loggraph
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
			const std::vector<int>& runOffsets,
			phaser_io::Output& output);
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
			phaser_io::Output& output);
//--------------------------------------------------------------
void PrintNormalisationResult(const Normalise& NormRes,
			      const ResoRange& ResRange, const double& MinIsigRatio,
			      phaser_io::Output& output);
//--------------------------------------------------------------
void PrintReindexSummary(const RefListType& reflisttype,
			 const IO_files& Allfiles,
			 phaser_io::Output& output);
//--------------------------------------------------------------
void PrintUnmergedHeaderStuff(const scala::hkl_unmerge_list& hkl_list,
			      phaser_io::Output& output,
			      const int& verbose);
//--------------------------------------------------------------
//--------------------------------------------------------------
void PrintError(const std::string& ErrorMessage,
		phaser_io::Output& output);
//--------------------------------------------------------------
void PrintWarning(const std::string& ErrorMessage,
		  phaser_io::Output& output);
#endif
