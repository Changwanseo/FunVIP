/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.2: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses.

    readAl v1.2: a tool for automated alignment conversion among different
                 formats.

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
                       [scapella, tgabaldon]@crg.es

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "alignment.h"
#include "rwAlignment.cpp"
#include "autAlignment.cpp"

alignment::alignment(void) {

  /* Init variable value to 0 */
  sequenNumber = 0;
  residNumber =  0;
  isAligned = false;

  iformat = 0;
  oformat = 0;

  dataType = 0;

  ghWindow = 0;
  shWindow = 0;

  /* New Info */
  oldAlignment  = false;
  saveResidues  = NULL;
  saveSequences = NULL;

  /* And pointer value to NULL */
  inputFileName =   NULL;
  alignmentInfo =   NULL;

  alignmentMatrix = NULL;
  residuesNumber  = NULL;
  sequenNames =     NULL;
  seqInfo =         NULL;

  sgaps = 	    NULL;
  scons =           NULL;
  seqMatrix = 	    NULL;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| alignment::alignment( int _sequenNumber, ..., ... ) 	|
|  	Class constructor                               |
|  							|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment::alignment(char *inputFileName_, char *alignmentInfo_, char **Matrix, char **Names, char **seqInfo_, int _sequenNumber, 
                     int _residNumber, int _iformat, int _oformat, int _dataType, int isAligned_, int _OldSequences, int _OldResidues,
                     int *residuesNumber_, int *saveResidues_, int *saveSequences_, int _ghWindow, int _shWindow) {
  int i, j, length;

  /* Assign the parameter values to the variables */
  sequenNumber = _sequenNumber;
  residNumber =  _residNumber;

  iformat = _iformat;
  oformat = _oformat;

  dataType = _dataType;

  ghWindow = _ghWindow;
  shWindow = _shWindow;

  isAligned = isAligned_;

  if(inputFileName_ != NULL) {
    length = strlen(inputFileName_);
    inputFileName = new char[length + 1];
    strcpy(inputFileName, inputFileName_);
  } else { inputFileName = NULL; }

  if(alignmentInfo_ != NULL) {
    length = strlen(alignmentInfo_);
    alignmentInfo = new char[length + 1];
    strcpy(alignmentInfo, alignmentInfo_);
  } else { alignmentInfo = NULL; }

  sequenNames = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {	
    length = strlen(Names[i]);
    sequenNames[i] = new char[length + 1];
    strcpy(sequenNames[i], Names[i]);
  }

  alignmentMatrix = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    alignmentMatrix[i] = new char[residNumber + 1];
    strcpy(alignmentMatrix[i], Matrix[i]);
  }

  residuesNumber = new int[sequenNumber];
  if((isAligned) || (residuesNumber_ != NULL)) {
    for(i = 0; i < sequenNumber; i++)
      residuesNumber[i] = residNumber;
  } else {
    for(i = 0; i < sequenNumber; i++)
      residuesNumber[i] = residuesNumber_[i];
  }

  if(seqInfo_ != NULL) {
    seqInfo = new char*[sequenNumber];

    for(i = 0; i < sequenNumber; i++) {
      length = strlen(seqInfo_[i]);
      seqInfo[i] = new char[length + 1];
      strcpy(seqInfo[i], seqInfo_[i]);
    }
  } else { seqInfo = NULL; }

  sgaps =     NULL;
  scons =     NULL;
  seqMatrix = NULL;

  oldAlignment = true;
  saveResidues  = NULL;
  saveSequences = NULL;

  if(saveResidues_ != NULL) {
    saveResidues = new int[residNumber];
    for(i = 0, j = 0; i < _OldResidues; i++)
      if(saveResidues_[i] != -1) {
        saveResidues[j] = saveResidues_[i];
        j++;
      }
  }

  if(saveSequences_ != NULL) {
    saveSequences = new int[sequenNumber];
    for(i = 0, j = 0; i < _OldSequences; i++)
      if(saveSequences_[i] != -1) {
        saveSequences[j] = saveSequences_[i];
        j++;
      }
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| alignment::~alignment()			 	|
|  	Class destructor				|
|  							|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  
alignment::~alignment(void) {

  if(alignmentMatrix != NULL)
    freeAlignment();
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::freeAlignment() 	                |
|  	Attribute memory deletion                       |
|  							|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::freeAlignment(void) {
  int i;

  for(i = 0; i < sequenNumber; i++)
    delete[] alignmentMatrix[i];
  delete [] alignmentMatrix;

  for(i = 0; i < sequenNumber; i++)
    delete[] sequenNames[i];
  delete[] sequenNames;

  if(residuesNumber != NULL)
    delete [] residuesNumber;

  if(inputFileName != NULL)
    delete[] inputFileName;

  if(alignmentInfo != NULL) 
    delete[] alignmentInfo;

  if(seqInfo != NULL) {
    for(i = 0; i < sequenNumber; i++)
      delete[] seqInfo[i];
    delete[] seqInfo;
  }

  if(sgaps != NULL)
    delete sgaps;

  if(scons != NULL)
    delete scons;

  if(seqMatrix != NULL)
    delete seqMatrix;

  if(saveResidues != NULL)
    delete[] saveResidues;

  if(saveSequences != NULL)
    delete[] saveSequences;

  oldAlignment = false;

  sequenNumber = 0;
  residNumber =  0;

  isAligned = false;

  iformat =  0;
  oformat =  0;

  dataType = 0;

  ghWindow = 0;
  shWindow = 0;

  inputFileName = NULL;
  alignmentInfo = NULL;

  seqInfo =       NULL;

  alignmentMatrix = NULL;
  sequenNames    = NULL;
  residuesNumber  = NULL;

  sgaps =     NULL;
  scons =     NULL;
  seqMatrix = NULL;

  saveResidues  = NULL;
  saveSequences = NULL;
}


alignment &alignment::operator=(const alignment &old) {

   int i;

   if(this != &old) {

     /* Assign the parameter values to the variables */
    sequenNumber = old.sequenNumber;
    residNumber =  old.residNumber;

    isAligned =  old.isAligned;

    iformat = old.iformat;
    oformat = old.oformat;

    dataType = old.dataType;

    ghWindow = old.ghWindow;
    shWindow = old.shWindow;

    delete [] inputFileName;
    if(old.inputFileName) {
      inputFileName = new char[strlen(old.inputFileName) + 1];
      strcpy(inputFileName, old.inputFileName);
    } else inputFileName = NULL;

    delete [] alignmentInfo;
    if(old.alignmentInfo) {
      alignmentInfo = new char[strlen(old.alignmentInfo) + 1];
      strcpy(alignmentInfo, old.alignmentInfo);
    } else alignmentInfo = NULL;

    delete [] sequenNames;
    if(old.sequenNames) {
      sequenNames = new char*[sequenNumber];
      for(i = 0; i < sequenNumber; i++) {
        sequenNames[i] = new char[strlen(old.sequenNames[i]) + 1];
        strcpy(sequenNames[i], old.sequenNames[i]);
      }
    } else sequenNames = NULL;

    delete [] alignmentMatrix;
    if(old.alignmentMatrix) {
      alignmentMatrix = new char*[sequenNumber];
      for(i = 0; i < sequenNumber; i++) {
        alignmentMatrix[i] = new char[residNumber + 1];
        strcpy(alignmentMatrix[i], old.alignmentMatrix[i]);
      }
    } else alignmentMatrix = NULL;

    delete [] residuesNumber;
    if(old.residuesNumber) {
      residuesNumber = new int[sequenNumber];
      for(i = 0; i < sequenNumber; i++)
        residuesNumber[i] = old.residuesNumber[i];
    } residuesNumber = NULL;

    delete [] seqInfo;
    if(old.seqInfo) {
      seqInfo = new char*[sequenNumber];
      for(i = 0; i < sequenNumber; i++) {
        seqInfo[i] = new char[strlen(old.seqInfo[i]) + 1];
        strcpy(seqInfo[i], old.seqInfo[i]);
      }
    } else seqInfo = NULL;

    delete sgaps;
    sgaps = NULL;

    delete scons;
    scons = NULL;

    delete seqMatrix;
    seqMatrix = old.seqMatrix;

    oldAlignment = true;

    delete [] saveResidues;
    if(old.saveResidues) {
      saveResidues = new int[residNumber];
      for(i = 0; i < residNumber; i++)
        saveResidues[i] = old.saveResidues[i];
    } else saveResidues = NULL;

    delete [] saveSequences;
    if(old.saveSequences) {
      saveSequences = new int[sequenNumber];
      for(i = 0; i < sequenNumber; i++)
        saveSequences[i] = old.saveSequences[i];
    } else saveSequences = NULL;
  }

  return *this;
}

bool alignment::loadAlignment(char *alignmentFile) {

  iformat = formatInputAlignment(alignmentFile);
  oformat = iformat;

  switch(iformat) {

    case 1:
      return loadClustalAlignment(alignmentFile);

    case 3:
      return loadNBRF_PirAlignment(alignmentFile); 

    case 8:
      return loadFastaAlignment(alignmentFile);

    case 11:
      return loadPhylip3_2Alignment(alignmentFile);

    case 12:
      return loadPhylipAlignment(alignmentFile);

    case 17:
      return loadNexusAlignment(alignmentFile);

    case 21:
      return loadMegaInterleavedAlignment(alignmentFile);

    case 22:
      return loadMegaNonInterleavedAlignment(alignmentFile);

    default:
      return false;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| int alignment::formatInputFile()			|
|							|
|							|
|							|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int alignment::formatInputFile(void) {

  return iformat;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| int alignment::formatInputFile()			|
|							|
|							|
|							|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int alignment::typeInputFile(void) {

  return SINGLE;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| bool alignment::printAlignment()			|
|  	This method prints the alignment to the standard|
|	output in PHYLIP2 format                        |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::printAlignment(void){

  if(alignmentMatrix == NULL) return false;

  switch(oformat) {

    case 1: alignmentClustalToFile(cout); break;

    case 3: alignmentNBRF_PirToFile(cout); break;

    case 8: alignmentFastaToFile(cout); break;

    case 11: alignmentPhylip3_2ToFile(cout); break;

    case 12: alignmentPhylipToFile(cout); break;

    case 17: alignmentNexusToFile(cout); break;

    case 21:
    case 22:
      alignmentMegaToFile(cout); 
      break;

    default: return false;

  }

  return true;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| bool alignment::saveAlignment(char *destFile)		|
|  	This method puts the alignment on the file      |
|  	in PHYLIP2 format				|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::saveAlignment(char *destFile) {

  ofstream file;

  if(alignmentMatrix == NULL) 
    return false;

  /* File open and correct open check */
  file.open(destFile);
  if(!file) return false;

  /* Alignment saving */
  switch(oformat) {

    case 1: alignmentClustalToFile(file); break;

    case 3: alignmentNBRF_PirToFile(file); break;

    case 8: alignmentFastaToFile(file); break;

    case 11: alignmentPhylip3_2ToFile(file); break;

    case 12: alignmentPhylipToFile(file); break;

    case 17: alignmentNexusToFile(file); break;

    case 21:
    case 22: 
      alignmentMegaToFile(file);
      break;

    default: return false;

  }
  /* Destination file closing */
  file.close();

  /* All was OK, return true */
  return true;

}

/* ****************************************************************************************************************** */
/* ******************************************** New Code ************************************************************ */
/* ********************************************************************************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanGaps(float baseLine, float gapsPct, bool complementary = false)                          |
|                                                                                                                      |
|       This method cleans the alignment using gap's method. It takes two inputs parameters, these parameters are:     |
|       - baseline: Columns' percentage of the alignment that we have to conserve.                                     |
|       - gapsPct: Gaps' percentage per column that we will have as maximum in the new alignment.                      |
|       The method calls other methods to select the correct threshold to clean the alignment.                         |
|       This method allow to get the complementary alignemnt. This complementary alignment have the columns rejects for|
|       the selection's method.                                                                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanGaps(float baseLine, float gapsPct, bool complementary) {

  alignment *ret;
  double cut;

  /* If gap's statistics are not calculated, we calculate them */
  if(calculateGapStats() != true) 
    return NULL;

  /* Obtain the cut point using the given parameters */
  cut = sgaps -> calcCutPoint(baseLine, gapsPct);

  /* Once we have the cut value proposed, we call the appropiate method to clean the alignment and, then, generate 
     the new alignment. */
  ret = cleanByCutValue(cut, baseLine, sgaps -> getGapsWindow(), complementary);

  /* Return a reference of the new alignment */
  return ret;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanConservation(float baseLine, float conservationPct, bool complementary)                  |
|                                                                                                                      |
|       This method cleans the alignment using conservation's method. It takes two inputs parameters, these parameters |
|       are:                                                                                                           |
|       - baseline: Columns' percentage of the alignment that we have to conserve.                                     |
|       - conservationPct: Conservation's value (0 - 1) per column that we will have as minimum in the new alignment.  |
|       The method calls other method to select the correct threshold to clean the alignment.                          |
|       This method allow to get the complementary alignemnt. This complementary alignment have the columns rejects for|
|       the selection's method.                                                                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanConservation(float baseLine, float conservationPct, bool complementary) {

  alignment *ret;
  float cut;

  /* If conservation's statistics are not calculated, we calculate them */
  if(calculateConservationStats() != true)
    return NULL;

  /* Calculate the cut point using the given parameters */
  cut = scons -> calcCutPoint(baseLine, conservationPct);

  /* Once we have the cut value proposed, we call the appropiate method to clean the alignment and, then, generate 
     the new alignment. */
  ret = cleanByCutValue(cut, baseLine, scons -> getMdkwVector(), complementary);

  /* Return a reference of the new alignment */
  return ret;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::clean(float baseLine, float gapsPct, float conservationPct, bool complementary)               |
|                                                                                                                      |
|       This method cleans the alignment using a combination between the conservation's method and gap's method. It    |
|       takes three inputs parameters, these parameters are:                                                           |
|       - baseline: Columns' percentage of the alignment that we have to conserve.                                     |
|       - gapsPct: Gaps' percentage per column that we will have as maximum in the new alignment.                      |
|       - conservationPct: Conservation's value (0 - 1) per column that we will have as minimum in the new alignment.  |
|       The method calls other methods to select the thresholds to clean the alignment. Finally, the method calls to   |
|       the cleaning's method that will combine the two thresholds to select one of them or, in somecases, select the  |
|       baseline value as the columns' number to conserve in the new alignment.                                        |
|       This method allow to get the complementary alignemnt. This complementary alignment have the columns rejects for|
|       the cleaning's method.                                                                                         |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::clean(float baseLine, float GapsPct, float conservationPct, bool complementary) {

  alignment *ret;
  float cutCons;
  double cutGaps;

  /* If gap's statistics are not calculated, we calculate them */
  if(calculateGapStats() != true)
    return NULL;

  /* If conservation's statistics are not calculated, we calculate them */
  if(calculateConservationStats() != true)
    return NULL;

  /* Calculate the two cut points using the given parameters */
  cutGaps = sgaps->calcCutPoint(baseLine, GapsPct);
  cutCons = scons->calcCutPoint(baseLine, conservationPct);

  /* Clean the alingment using the two cut values, the gapsWindow and MDK_Windows vectors and the baseline value */
  ret = cleanByCutValue(cutGaps, sgaps -> getGapsWindow(), baseLine, cutCons, scons -> getMdkwVector(), complementary);

  /* Return a reference of the clean alignment object */
  return ret;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary)       |
|                                                                                                                      |
|       This method cleans the alignment using three inputs parameters, these parameters are:                          |
|       - cutpoint: Lower limit (0-1) of comparefile value admits in the new alignment.                                |
|       - baseline: Columns' percentage of the alignment that we have to conserve.                                     |
|       - vectValues: A vector with alignment's comparefile values.                                                    |
|       The method calls other method to select the correct threshold to clean the alignment.                          |
|       This method allow to get the complementary alignemnt. This complementary alignment have the columns rejects for|
|       the selection's method.                                                                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary) {

  alignment *ret;
  float cut, *vectAux;

  /* Allocate memory */
  vectAux = new float[residNumber];

  /* Sort a copy of the vectValues vector, and take the value at 100% - baseline position. */
  utils::copyVect((float *) vectValues, vectAux, residNumber);
  utils::quicksort(vectAux, 0, residNumber-1);
  cut = vectAux[(int) ((float)(residNumber - 1) * (100.0 - baseLine)/100.0)];

  /* We have to decide which is the smallest value between cutpoint value 
     and value from minimum percentage threshold */
  cut = cutpoint < cut ? cutpoint : cut;

  /* Deallocate memory */
  delete [] vectAux;

  /* Clean the selected alignment using the input parameters. */
  ret = cleanByCutValue(cut, baseLine, vectValues, complementary);

  /* Return a refernce of the new alignment */
  return ret;

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanSpuriousSeq(float, float)                                                                |
|                                                                                                                      |
|       This method selects and removes the sequences that have a overlap's value less than minimum overlap's cut point|
|       For this purpose, the methods calls a other method to compute the overlap's values for each sequences. To      |
|       computes this value, the method called uses the overlap column value to decide if a column for a sequence      |
|       given is greater or equal to this value or not.
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanSpuriousSeq(float overlapColumn, float minimumOverlap) {

  float *overlapVector;
  alignment *newAlig;

  /* Allocate local memory */
  overlapVector = new float[sequenNumber];

  /* Compute the overlap's vector using the overlap column's value */
  if(!calculateSpuriousVector(overlapColumn, overlapVector))
    return NULL;

  /* Select and remove the sequences with a overlap less than threshold's overlap then create the new alignemnt */
  newAlig = cleanOverlapSeq(minimumOverlap, overlapVector);

  /* Deallocate local memory */
  delete [] overlapVector;

  /* Return a reference of the clean alignment object */
  return newAlig;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanMixSlope(bool complementarity)                                                           |
|                                                                                                                      |
|       This method cleans the alignmnet using a automatic method. The automatic method compute the rate between the   |
|       first slope ratio between alignment's length and gaps percentage and the "second" slope ratio between          |
|       alignment's length and gaps percentage. As the all methods, this method will be able to compute the            |
|       complementary alignment.                                                                                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanMixSlope(bool complementarity) {

  alignment *ret;
  int cut;

  /* If gap statistics are not calculated, we calculate them */
  if(calculateGapStats() != true)
    return NULL;

  /* We get the cut point using a automatic method for this purpose. */
  cut = sgaps -> calcCutPointMixSlope();

  /* Using the cut point calculates in last steps, we clean the alignment and generate a new alignment. */
  ret = cleanByCutValue(cut, 0, sgaps->getGapsWindow(), complementarity);

  /* Returns the new alignment. */
  return ret;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::clean2ndSlope(bool complementarity)                                                           |
|                                                                                                                      |
|       This method cleans the alignmnet using a automatic method. The automatic method compute the "second" slope     |
|       ratio between alignment's length and gaps percentage. Moreover, this method will be able to compute the        |
|       complementary alignment.                                                                                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::clean2ndSlope(bool complementarity) {

  alignment *ret;
  int cut;

  /* If gap statistics are not calculated, we calculate them */
  if(calculateGapStats() != true) 
    return NULL;

  /* We get the cut point using a automatic method for this purpose. */
  cut = sgaps -> calcCutPoint2ndSlope();

  /* Using the cut point calculates in last steps, we clean the alignment and generate a new alignment. */
  ret = cleanByCutValue(cut, 0, sgaps->getGapsWindow(), complementarity);

  /* Returns the new alignment. */
  return ret;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanNoAllGaps(bool complementarity)                                                          |
|                                                                                                                      |
|       This method removes all columns with only gaps. Moreover, this method will be able to compute the complementary|
|       alignment.                                                                                                     |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanNoAllGaps(bool complementarity) {

  alignment *ret;

  /* If gap statistics are not calculated, we calculate them */
  if(calculateGapStats() != true) 
    return NULL;

  /* We want to conserve the columns with gaps' number less or equal than sequences' number - 1  */
  ret = cleanByCutValue((sequenNumber - 1), 0, sgaps->getGapsWindow(), complementarity);

  /* Returns the new alignment. */
  return ret;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanCombMethods(bool complementarity)                                                        |
|                                                                                                                      |
|       This method computes a cut point using for this purpose information related to gaps' statistics and            |
|       conservation's statistics. First, rejects some columns using the strict method to clening gaps then computes   |
|       a similarity's value using the columns that have not been rejected for this purpose. The method calls a        |
|       cleaning method using this value as input. This method also will be able to get the complementary alignment.   |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanCombMethods(bool complementarity, bool variable) {

  float simCut, first20Point, last80Point, *simil, *vectAux;
  int i, j, acm, gapCut, *positions, *gaps;
  double inic, fin, vlr;

  /* If gaps' statistics are not calculated, we calculate them */
  if(calculateGapStats() != true)
    return NULL;

  /* Computes the gap cut point using a automatic method and at the same time, we get the gaps'
     values of the alignment. */
  gapCut = sgaps -> calcCutPoint2ndSlope();
  gaps = sgaps -> getGapsWindow();

  /* If conservation's statistics are not calculated, we calculate them */
  if(calculateConservationStats() != true)
    return NULL;

  /* Computes the conservation's value for each column in the alignment. At the same time, the method get the
     vector with those values. */
  scons -> calculateVectors(alignmentMatrix, sgaps -> getGapsWindow());
  simil = scons -> getMdkwVector();

  /* Allocate local memory and initializate it to -1 */
  positions = new int[residNumber];
  utils::initlVect(positions, residNumber, -1);

  /* The method only selects columns with gaps' number less or equal than the gap's cut point. Counts the number of
     columns that have been selected */
  for(i = 0, acm = 0; i < residNumber; i++) {
    if(gaps[i] <= gapCut) {
      positions[i] = i;
      acm++;
    }
  }

  /* Allocate local memory and save the similarity's values for the columns that have been selected */
  vectAux = new float[acm];
  for(i = 0, j = 0; i < residNumber; i++)
    if(positions[i] != -1)
      vectAux[j++] = simil[i];

  /* Sort the conservation's value vector. */
  utils::quicksort(vectAux, 0, acm-1);

  /* ... and search for the vector points at the 20 and 80% of length. */
  first20Point = 0;
  last80Point  = 0;

  for(i = acm - 1, j = 1; i >= 0; i--, j++) {
    if((((float) j/acm) * 100.0) <= 20.0) 
      first20Point = vectAux[i];
    if((((float) j/acm) * 100.0) <= 80.0) 
      last80Point = vectAux[i];
  }

  /* Computes the logaritmic's values for those points. Finally the method computes the similarity cut point 
     using these values. */
  inic = log10(first20Point);
  fin  = log10(last80Point);
  vlr  = ((inic - fin) / 10) + fin;
  simCut = (float) pow(10, vlr);

  /* Clean the alignment and generate a new alignment object using the gap's cut and the similarity's cut values */
  alignment *ret = cleanStrictPlus(gapCut, sgaps -> getGapsWindow(), simCut, scons -> getMdkwVector(), complementarity, variable);

  /* Deallocate local memory */
  delete [] vectAux;
  delete [] positions;

  /* Return a reference of the new alignment */
  return ret;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanByCutValue(int, float, const int *, bool)                                                |
|                                                                                                                      |
|       This method uses a given cut point value over a gapswindow vector to select the columns in the alignment to    |
|       conserve. At the same time, the method uses the baseline value to add some columns to the new alignment.       |
|       This method can be used to compute the complementary alignment, in other words, it can be used to get the      |
|       removed columns from the original alignment.                                                                   |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanByCutValue(double cut, float baseLine, const int *gInCol, bool complementary) {

  int i, j, k, jn, oth, newResidNumber, *vectAux;
  alignment *newAlig;
  char **matrixAux;

  /* Select the columns with a gaps value less or equal than the cut point. */
  for(i = 0, newResidNumber = 0; i < residNumber; i++)
    if(gInCol[i] <= cut) newResidNumber++;
    else saveResidues[i] = -1;

  /* Compute, if it's necessary, the column's number needed to achieve the new alignment column's number fixed by 
     conservation value. */
  oth = utils::roundInt((((baseLine/100.0) - (float) newResidNumber/residNumber)) * residNumber);

  if(oth > 0) {
    newResidNumber += oth;

    /* Allocate memory */
    vectAux = new int[residNumber];

    /* Sort a copy of the gInCol vector, and take the value of the column that marks the % baseline */
    utils::copyVect((int *) gInCol, vectAux, residNumber);
    utils::quicksort(vectAux, 0, residNumber-1);
    cut = vectAux[(int) ((float)(residNumber - 1) * (baseLine)/100.0)];

    /* Deallocate memory */
    delete [] vectAux;
  }

  /* Fixed the initial size of blocks as 0.5% of alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* We start in the alignment middle then we move on right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

      /* if block's size is greater or equal than the fixed size then we save all columns that
         have not been saved previously. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if(gInCol[jn] <= cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;

      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* if block's size is greater or equal than the fixed size then we save all columns that
         have not been saved previously. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if(gInCol[jn] <= cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
    }
  }

  /* Once we've selected the columns, if the complementary flag is true, we will have to change the selected and
    non-selected columns. */
  if(complementary == true) {
    newResidNumber = residNumber - newResidNumber;
    for(i = 0; i < residNumber; i++) {
      if(saveResidues[i] == -1) saveResidues[i] = i;
      else saveResidues[i] = -1;
    }
  }

  /* We allocate memory to save the columns selected. */
  matrixAux = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i] = new char[newResidNumber+1];

  /* Copy the columns from the original alignment to the new alignment, only if the column has been selected. */
  for(i = 0, k = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = 0; j < sequenNumber; j++)
        matrixAux[j][k] = alignmentMatrix[j][i];
      k++;
    }

  /* Concatenate the character '\0' to the last position of all rows of the new alignment to facilitate character 
     vector manipulation and printing functions. */
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i][newResidNumber] = '\0';

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, sequenNames, seqInfo, sequenNumber, newResidNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);

  /* Deallocated auxiliar memory */
  for(i = 0; i < sequenNumber; i++)
    delete[] matrixAux[i]; 
  delete[] matrixAux;

  /* Return the new alignment reference */
  return newAlig;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanByCutValue(float, float, const float *, bool)                                            |
|                                                                                                                      |
|       This method uses a given cut point value over a conservationwindow vector to select the columns in the         |
|       alignment to conserve. At the same time, the method uses the baseline value to add some columns to the new     |
|       alignment. This method can be used to compute the complementary alignment, in other words, it can be used to   |
|       get the removed columns from the original alignment. Like ohers methods, this method can obtain the            |
|       complementary alignment.                                                                                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanByCutValue(float cut, float baseLine, const float *ValueVect, bool complementary) {

  int i, j, k, jn, oth, newResidNumber;
  alignment *newAlig;
  char **matrixAux;

  /* Select the columns with a conservation's value greater than the cut point. */
  for(i = 0, newResidNumber = 0; i < residNumber; i++)
    if(ValueVect[i] > cut) newResidNumber++;
    else saveResidues[i] = -1;

  /* Compute, if it's necessary, the column's number needed to achieve the new alignment column's number fixed by 
     conservation value. */
  oth = utils::roundInt((((baseLine/100.0) - (float) newResidNumber/residNumber)) * residNumber);
  if(oth > 0) newResidNumber += oth;

  /* Fixed the initial size of blocks as 0.5% of alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* We start in the alignment middle then we move on right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

      /* if block's size is greater or equal than the fixed size then we save all columns that
         have not been saved previously. */
      /* Here, we only accept column with a conservation's value equal to conservation cut point. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if(ValueVect[jn] == cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;

      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* if block's size is greater or equal than the fixed size then we select the column and save the block's size
         for the next iteraction. it's obvius that we decrease the column's number needed to finish. */
     /* Here, we only accept column with a conservation's value equal to conservation cut point. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if(ValueVect[jn] == cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
    }
  }

  /* Once we've selected the columns, if the complementary flag is true, we will have to change the selected and
    non-selected columns. */
  if(complementary == true) {
    newResidNumber = residNumber - newResidNumber;
    for(i = 0; i < residNumber; i++) {
      if(saveResidues[i] == -1) saveResidues[i] = i;
      else saveResidues[i] = -1;
    }
  }

  /* We allocate memory to save the columns selected. */
  matrixAux = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i] = new char[newResidNumber+1];

  /* Copy the columns from the original alignment to the new alignment, only if the column has been selected. */
  for(i = 0, k = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = 0; j < sequenNumber; j++)
        matrixAux[j][k] = alignmentMatrix[j][i];
      k++;
    }

  /* Concatenate the character '\0' to the last position of all rows of the new alignment to facilitate character 
     vector manipulation and printing functions. */
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i][newResidNumber] = '\0';

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, sequenNames, seqInfo, sequenNumber, newResidNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);

  /* Deallocated auxiliar memory */
  for(i = 0; i < sequenNumber; i++)
    delete[] matrixAux[i]; 
  delete[] matrixAux;

  /* Return the new alignment reference */
  return newAlig;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanByCutValue(int, const int *, float, float, const float *, bool)                          |
|                                                                                                                      |
|       This method cleans the alignment using two approaches. The method selects the columns with a conservation's    |
|       value greater than conservation's cut point AND with a gaps' number less or equal than gap's cut point. If it's|
|       necessary add some columns to achieve the minimum baseline, the method selects columns with a conservation's   |
|       value greater or equal than a new conservation's cut point OR with a gaps' number less or equal than a new     |
|       gaps' cut point. Like ohers methods, this method can obtain the complementary alignment.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanByCutValue(double cutGaps, const int *gInCol, float baseLine, float cutCons, const float *MDK_Win, bool complementary) {

  int i, j, k, oth, jn, newResidNumber, blGaps, *vectAuxGaps;
  float blCons, *vectAuxCons;
  alignment *newAlig;
  char **matrixAux;

  /* Select the columns with a conservation's value greater than the conservation cut point AND less or equal than the gap cut point. */
  for(i = 0, newResidNumber = 0; i < residNumber; i++)
    if((MDK_Win[i] > cutCons) && (gInCol[i] <= cutGaps)) newResidNumber++;
    else saveResidues[i] = -1;

  /* Compute, if it's necessary, the column's number needed to achieve the new alignment column's number fixed by 
     conservation value. */
  oth = utils::roundInt((((baseLine/100.0) - (float) newResidNumber/residNumber)) * residNumber);

  /* If it's needed to add new columns, we compute the news thresholds */
  if(oth > 0) {
    newResidNumber += oth;

    /* Allocate memory */
    vectAuxCons = new float[residNumber];
    vectAuxGaps = new int[residNumber];

    /* Sort a copy of the MDK_Win vector and of the gInCol vector, and take the value of the column that marks the % baseline */
    utils::copyVect((float *) MDK_Win, vectAuxCons, residNumber);
    utils::copyVect((int *) gInCol, vectAuxGaps, residNumber);

    utils::quicksort(vectAuxCons, 0, residNumber-1);
    utils::quicksort(vectAuxGaps, 0, residNumber-1);

    blCons = vectAuxCons[(int) ((float)(residNumber - 1) * (100.0 - baseLine)/100.0)];
    blGaps = vectAuxGaps[(int) ((float)(residNumber - 1) * (baseLine)/100.0)];

    /* Deallocate memory */
    delete [] vectAuxCons;
    delete [] vectAuxGaps;
  }

  /* Fixed the initial size of blocks as 0.5% of alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* We start in the alignment middle then we move on right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

      /* if block's size is greater or equal than the fixed size then we select the column and save the block's size
         for the next iteraction. it's obvius that we decrease the column's number needed to finish. */
      /* Here, we accept column with a conservation's value greater or equal than the conservation cut point OR less or
         equal than the gap cut point. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;

      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* if block's size is greater or equal than the fixed size then we select the column and save the block's size
         for the next iteraction. it's obvius that we decrease the column's number needed to finish. */
      /* Here, we accept column with a conservation's value greater or equal than the conservation cut point OR less or
         equal than the gap cut point. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
    }
  }

  /* Once we've selected the columns, if the complementary flag is true, we will have to change the selected and
    non-selected columns. */
  if(complementary == true) {
    newResidNumber = residNumber - newResidNumber;
    for(i = 0; i < residNumber; i++) {
      if(saveResidues[i] == -1) saveResidues[i] = i;
      else saveResidues[i] = -1;
    }
  }

  /* We allocate memory to save the columns selected. */
  matrixAux = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i] = new char[newResidNumber+1];

  /* Copy the columns from the original alignment to the new alignment, only if the column has been selected. */
  for(i = 0, k = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = 0; j < sequenNumber; j++)
        matrixAux[j][k] = alignmentMatrix[j][i];
      k++;
    }

  /* Concatenate the character '\0' to the last position of all rows of the new alignment to facilitate character 
     vector manipulation and printing functions. */
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i][newResidNumber] = '\0';

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, sequenNames, seqInfo, sequenNumber, newResidNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);

  /* Deallocated auxiliar memory */
  for(i = 0; i < sequenNumber; i++)
    delete[] matrixAux[i]; 
  delete[] matrixAux;

  /* Return the new alignment reference */
  return newAlig;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanStrictPlus(int, const int *, float, const float *, bool, bool)                           |
|                                                                                                                      |
|       This method cleans the alignment using two approaches. The method selects the columns with a conservation's    |
|       value greater than conservation's cut point AND with a gaps' number less or equal than gap's cut point. If it's|
|       necessary add some columns to achieve the minimum baseline, the method selects columns with a conservation's   |
|       value greater or equal than a new conservation's cut point OR with a gaps' number less or equal than a new     |
|       gaps' cut point. Like ohers methods, this method can obtain the complementary alignment.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanStrictPlus(int gapCut, const int *gInCol, float simCut, const float *MDK_W, bool complementary, bool variable) {

  int i, k, j = 0, num, newResidNumber, lenBlock;
  alignment *newAlig;
  char **matrixAux;

  /* Rejects the columns with gaps' number greater than the gap's cut point. */
  for(i = 0; i < residNumber; i++)
    if(gInCol[i] > gapCut)
      saveResidues[i] = -1;

  /* Rejects the columns with conservation'value less than the conservation's cut point. */
  for(i = 0; i < residNumber; i++)
    if(MDK_W[i] < simCut)
      saveResidues[i] = -1;

  /* Search for columns that has been rejected and that has 3 or more adjacent columns selected */
  /* For the second column in the alignment */
  if((saveResidues[0] != -1) && (saveResidues[2] != -1) && (saveResidues[3] != -1)) 
    saveResidues[1] = 1;
  else
    saveResidues[1] = -1;

  /* For the penultimate column in the alignment */
  if((saveResidues[residNumber-1] != -1) && (saveResidues[residNumber-3] != -1) && (saveResidues[residNumber-4] != -1))
    saveResidues[(residNumber - 2)] = (residNumber - 2);
  else
    saveResidues[(residNumber - 2)] = -1;

  /* For the rest of columns in the alignment */
  for(i = 2, num = 0; i < (residNumber - 2); i++, num = 0)
    if(saveResidues[i] == -1) {
      if(saveResidues[(i - 2)] != -1) num++;
      if(saveResidues[(i - 1)] != -1) num++;
      if(saveResidues[(i + 1)] != -1) num++;
      if(saveResidues[(i + 2)] != -1) num++;
      if(num >= 3) saveResidues[i] = i;
    }


  /* Fix the block's size based on variable flag. The Block's size can be fixed to 5 or can be variable between a
     Minimum Block's size 3 and a Maximum Block's size: 12 depend on percentage of alignment's length. */
  if(!variable)
    lenBlock = 5;
  else {
    lenBlock = utils::roundInt(residNumber * 0.01) > 3 ? utils::roundInt(residNumber * 0.01) : 3;
    lenBlock = lenBlock < 12 ? lenBlock : 12;
  }

  /* The method searchs for columns' blocks with LONGBLOCK or greater size */
  for(i = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = (i + 1); ((j < residNumber) && (saveResidues[j] != -1)); j++) ;
      if((j - i) < lenBlock)
        for(k = i; k < j; k++)
          saveResidues[k] = -1;
      i = j;
    }

  /* Finally, the method computes the new alignment columns' number */
  for(i = 0, newResidNumber = 0; i < residNumber; i++) 
    if(saveResidues[i] != -1)
      newResidNumber++;

  /* Once we've selected the columns, if the complementary flag is true, we will have to change the selected and
    non-selected columns. */
  if(complementary == true) {
    newResidNumber = residNumber - newResidNumber;
    for(i = 0; i < residNumber; i++) {
      if(saveResidues[i] == -1) saveResidues[i] = i;
      else saveResidues[i] = -1;
    }
  }

  /* We allocate memory to save the columns selected. */
  matrixAux = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i] = new char[newResidNumber+1];

  /* Copy the columns from the original alignment to the new alignment, only if the column has been selected. */
  for(i = 0, k = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = 0; j < sequenNumber; j++)
        matrixAux[j][k] = alignmentMatrix[j][i];
      k++;
    }

  /* Concatenate the character '\0' to the last position of all rows of the new alignment to facilitate character 
     vector manipulation and printing functions. */
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i][newResidNumber] = '\0';

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, sequenNames, seqInfo, sequenNumber, newResidNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);

  /* Deallocated auxiliar memory */
  for(i = 0; i < sequenNumber; i++)
    delete[] matrixAux[i]; 
  delete[] matrixAux;

  /* Return the new alignment reference */
  return newAlig;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  alignment *alignment::cleanOverlapSeq(float, float *)                                                               |
|                                                                                                                      |
|       This method cleans the alignment using the overlaps' values for each sequences in the alignment.               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

alignment *alignment::cleanOverlapSeq(float minimumOverlap, float *overlapSeq) {

  char **matrixAux, **newSpeciesNames;
  int i, j, lenNames, newSequences;
  alignment *newAlig;

  /* Computes the new sequences' number. For this purpose, selects the sequences with a overlap's value equal or
     greater than the minimum overlap value. At the same time, computes the sequence's name length. */
  for(i = 0, newSequences = 0, lenNames = 0; i < sequenNumber; i++) {
    if(overlapSeq[i] >= minimumOverlap) {
      newSequences++;
      lenNames = lenNames > (int) strlen(sequenNames[i]) ? lenNames : (int) strlen(sequenNames[i]);
    } else {
      saveSequences[i] = -1;
    }
  }

  /* Allocate memory for the sequences selected */
  matrixAux = new char*[newSequences];
  newSpeciesNames = new char*[newSequences];

  /* Initializate the memory allocated */
  for(i = 0; i < newSequences; i++) {
    matrixAux[i] = new char[residNumber+1];
    matrixAux[i][residNumber] = '\0';

    newSpeciesNames[i] = new char[lenNames+1];
    newSpeciesNames[i][lenNames] = '\0';
  }

  /* Copy to new structures the information that have been selected previously. */
  for(i = 0, j = 0; i < sequenNumber; i++)
    if(overlapSeq[i] >= minimumOverlap) {
       strcpy(newSpeciesNames[j], sequenNames[i]);
       strcpy(matrixAux[j], alignmentMatrix[i]);
       j++;
    }

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, newSpeciesNames, seqInfo, newSequences, residNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);

  /* Deallocate the local memory */
  for(i = 0; i < newSequences; i++) {
    delete [] matrixAux[i];
    delete [] newSpeciesNames[i];
  }

  delete [] matrixAux;
  delete [] newSpeciesNames;

  /* Return the new alignment reference */
  return newAlig;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool alignment::calculateSpuriousVector(float, float *)                                                             |
|                                                                                                                      |
|       This method computes the spurious's or overlap's value for each alignment's sequences using an input parameter |
|       as threshold to decides if a column for a sequence given can be considered as a hit or not.                    |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::calculateSpuriousVector(float overlap, float *spuriousVector) {

  int i, j, k, seqValue, ovrlap, hit;
  float floatOverlap;
  char indet;

  floatOverlap = overlap * float(sequenNumber-1);
  ovrlap = int(overlap * (sequenNumber-1));

  if(floatOverlap > float(ovrlap))
    ovrlap++;

  /* If the spurious vectos is NULL, returns false. */
  if(spuriousVector == NULL)
    return false;

  if(getTypeAlignment() == AAType)
    indet = 'X';
  else
    indet = 'N';

  /* For each alignment's sequence, computes its overlap */
  for(i = 0, seqValue = 0; i < sequenNumber; i++, seqValue = 0) {

    /* For each alignment's column, computes the overlap between the selected sequence and the other ones */
    for(j = 0; j < residNumber; j++) {

      /* For sequences are before the sequence selected */
      for(k = 0, hit = 0; k < i; k++) {

        /* If the element of sequence selected is the same that the element of sequence considered, computes a hit */
        if(alignmentMatrix[i][j] == alignmentMatrix[k][j]) 
          hit++;

        /* If the element of sequence selected isn't a 'X' nor 'N' (indetermination) or a '-' (gap) and the element of
           sequence considered isn't a  a 'X' nor 'N' (indetermination) or a '-' (gap), computes a hit */
        else if((alignmentMatrix[i][j] != indet) && (alignmentMatrix[i][j] != '-')
          && (alignmentMatrix[k][j] != indet) && (alignmentMatrix[k][j] != '-'))
          hit++;
      }

      /* For sequences are after the sequence selected */
      for(k = (i + 1); k < sequenNumber; k++) {

        /* If the element of sequence selected is the same that the element of sequence considered, computes a hit */
        if(alignmentMatrix[i][j] == alignmentMatrix[k][j]) 
          hit++;

        /* If the element of sequence selected isn't a 'X' nor 'N' (indetermination) or a '-' (gap) and the element of
           sequence considered isn't a  a 'X' nor 'N' (indetermination) or a '-' (gap), computes a hit */
        else if((alignmentMatrix[i][j] != indet) && (alignmentMatrix[i][j] != '-')
          && (alignmentMatrix[k][j] != indet) && (alignmentMatrix[k][j] != '-'))
          hit++;
      }
      /* Finally, if the hit's number divided by sequences' number minus one is greater or equal than overlap's value,
         computes a column's hit. */     
      if(hit >= ovrlap)
        seqValue++;
    }

    /* For each alignment's sequence, computes its spurious's  or overlap's value as the column's hits -for that 
       sequence- divided by column's number. */
    spuriousVector[i] = ((float) seqValue / residNumber);
  }
  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool alignment::calculateGapStats(void)                                                                             |
|                                                                                                                      |
|       This method computes the gap's statistics for the alignment if it haven't been calculated yet                  |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::calculateGapStats(void) {

  /* If alignment matrix is not created, return false */
  if(alignmentMatrix == NULL)
    return false;

  /* If sgaps object is not created, we create them and calculate the statistics */
  if(sgaps == NULL) {
    sgaps = new statisticsGaps(alignmentMatrix, sequenNumber, residNumber, dataType);
    sgaps -> applyWindow(ghWindow);
  }

  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::printStatisticsGapsColumns(void)                                                                    |
|                                                                                                                      |
|       This method prints the gaps's number of each column in the alignment                                           |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printStatisticsGapsColumns(void) {

  /* We check if the gaps' statistics have been calculated */
  if(calculateGapStats())
    /* then prints the information */
    sgaps -> printGapsColumns();
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::printStatisticsGapsTotal(void)                                                                      |
|                                                                                                                      |
|       This method prints the the accumulated columns' number for each gaps' number in the alignment.                 |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printStatisticsGapsTotal(void) {

  /* We check if the gaps' statistics have been calculated */
  if(calculateGapStats())
    /* then prints the information */
    sgaps -> printGapsAcl();

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool alignment::setSimilarityMatrix(similarityMatrix *)                                                             |
|                                                                                                                      |
|       This method associates the similarity matrix pointer to conservation's statistics object.                      |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::setSimilarityMatrix(similarityMatrix *sm) {

  /* If scons object is not created, we create them */
  if(scons == NULL) 
    scons = new statisticsConservation(alignmentMatrix, sequenNumber, residNumber, dataType);

  /* If the similarity matrix is wrong, return false */
  if(!scons -> setSimilarityMatrix(sm)) 
    return false;

  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool alignment::calculateConservationStats(void)                                                                    |
|                                                                                                                      |
|       This method computes the conservation's statistics for the alignment if it hasn't been calculated yet          |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool alignment::calculateConservationStats(void) {

  /* If gap statistics are not calculated, we calculate them */
  if(calculateGapStats() != true)
    return false;

  /* If scons object is not created, return false */
  if(scons == NULL) 
    return false;

  if(scons -> isSimMatrixDef() != true)
    return false;

  if(!scons -> calculateVectors(alignmentMatrix, sgaps->getGapsWindow()))
    return false;

  if(scons->isDefinedWindow())
    return true;

  else return scons->applyWindow(shWindow);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::printStatisticsConservationColumns(void)                                                            |
|                                                                                                                      |
|       This method prints the conservation's value of each column in the alignment                                    |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printStatisticsConservationColumns(void) {

  /* We check if the conservation's statistics have been calculated */
  if(calculateConservationStats())
    /* then prints the information */
    scons -> printConservationColumns();
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::printStatisticsConservationTotal(void)                                                              |
|                                                                                                                      |
|       This method prints the accumulated columns' number for each conservation's value in the alignment              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printStatisticsConservationTotal(void) {

  /* We check if the conservation's statistics have been calculated */
  if(calculateConservationStats())
    /* then prints the information */
    scons -> printConservationAcl();
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation::printCorrespondence(void)                                                              |
|                                                                                                                      |
|       This method prints the correspondence between new alignment's columns and old alignment's columns.             |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printCorrespondence(void) {

  int i;

  cout << endl;
  for(i = 0; i < residNumber - 1; i++)
    cout << saveResidues[i] << ", ";
  cout << saveResidues[i] << endl << endl;

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  int alignment::getNumSpecies(void)                                                                                  |
|                                                                                                                      |
|       This method returns the alignment's sequences number                                                           |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int alignment::getNumSpecies(void) {

  return sequenNumber;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  int alignment::getNumAminos(void)                                                                                   |
|                                                                                                                      |
|       This method returns the alignment's columns number                                                             |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int alignment::getNumAminos(void) {

  return residNumber;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::setWindowsSize(int, int)                                                                            |
|                                                                                                                      |
|       This method sets the values to gap's and similarity's half windows size                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::setWindowsSize(int ghWindow_, int shWindow_) {

  ghWindow = ghWindow_;
  shWindow = shWindow_;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void alignment::setOutputFormat(int)                                                                                |
|                                                                                                                      |
|       This method allows to changes the default output format.                                                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::setOutputFormat(int format_) {

  oformat = format_;

}


/* ****************************************************************************************************************** */
/* ****************************************************************************************************************** */

/* ****************************************************************************************************************** */
/* ******************************************** New Code ************************************************************ */
/* ********************************************************************************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::getSpecies(char **Names)               |
|       Method that returns the alignment's sequenNumber name |
|                                                        |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::getSpecies(char **Names) {

  for(int i = 0; i < sequenNumber; i++)
    strcpy(Names[i], sequenNames[i]);

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::sequenMatrix(void)                       |
|       Method that builds an alignment's sequences matrix |
|                                                          |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::sequenMatrix(void) {

  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(alignmentMatrix, sequenNumber, residNumber);

}

void alignment::destroySequenMatrix(void) {

  if(seqMatrix != NULL)
    delete seqMatrix;
  seqMatrix = NULL;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::printSequenMatrix(void)                |
|       Method that prints the alignment's sequence      |
|       matrix                                           |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::printSequenMatrix(void) {

  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(alignmentMatrix, sequenNumber, residNumber);
  
  seqMatrix -> printMatrix();
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::getColumSeqMatrix(int, int *)          |
|       Method that gets a "column"-column from the      |
|       alignment's sequence matrix                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::getColumnSeqMatrix(int column, int *columnSeqMatrix) {

  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(alignmentMatrix, sequenNumber, residNumber);

  seqMatrix -> getColumn(column, columnSeqMatrix);
	
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::getColumSeqMatrix(int, int, int *)     |
|       Method that gets a column from the aligment's    |
|       sequence matrix with the same value that "value" |
|       at matrix's position (row, i)                    |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::getColumnSeqMatrix(int value, int row, int *columnSeqMatrix) {

  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(alignmentMatrix, sequenNumber, residNumber);

  seqMatrix -> getColumn(value, row, columnSeqMatrix);

}

/* ************************************************************************* */
void alignment::getColumn(int columnNumber, char *column) {
  int i;

  column[sequenNumber] = '\0';
  for(i = 0; i < sequenNumber; i++) 
    column[i] =  alignmentMatrix[i][columnNumber];
}

/* ****************************************************************************************************************** */
/* ****************************************************************************************************************** */


void alignment::checkTypeAlignment(void) {

  int i, j, k, l, hitDNA, hitRNA, gDNA, gRNA;
  char listRNA[6] = "AGCUN";
  char listDNA[6] = "AGCTN";

  /* For each sequences, this method locks at the 100 letters (excluding gaps). If 95% or more of those letters are 
     valid nucleotides, then the files is treated as nucleotides. The method also recognizes between ADN and ARN. */
  for(i = 0, gDNA = 0, gRNA = 0; i < sequenNumber; i++) {

    /* Looks at the 100 letters (excluding gaps) while doesn's get the sequence's end */
    for(j = 0, k = 0, hitDNA = 0, hitRNA = 0; ((j < residNumber) && (k  < 100)); j++)
      if(alignmentMatrix[i][j] != '-') {
        k++;

        /* Recognizes between DNA and RNA. */
        for(l = 0; l < (int) strlen(listDNA); l++)
          if(listDNA[l] == alignmentMatrix[i][j])
            hitDNA++;

        for(l = 0; l < (int) strlen(listRNA); l++)
          if(listRNA[l] == alignmentMatrix[i][j])
            hitRNA++;
      }

    /* If for an alignment's sequences the nucleotides don't achieve the threshold, then the method finish and fix
       the alignment's datatype as AminoAcids. */
    if((((float) hitDNA/k) < 0.95) && (((float) hitRNA/k) < 0.95)) {
      dataType = AAType;
      return;
    }

    /* Computes the greater value between DNA's nucleotides and RNA's nucleotides */
    else if(hitRNA > hitDNA)
      gRNA++;

    else
      gDNA++;
 }
  /* Return the datatype with the greater value */
  if(gRNA > gDNA)
    dataType = RNAType;

  else
    dataType = DNAType;
}

int alignment::getTypeAlignment(void) {

  if(dataType == 0)
    checkTypeAlignment();

  return dataType;
}

alignment *alignment::removeColumns(int *columns, int size, bool complementary) {

  alignment *ret;

  ret = cleanByColumnsNumber(columns, size, complementary);

  return ret;
}

alignment *alignment::cleanByColumnsNumber(int *columns, int size, bool complementary) {

  int i, j, k, delAminos, newResidNumber;
  alignment *newAlig;
  char **matrixAux;

  for(i = 0, delAminos = 0; i < size; i += 2) {
    for(j = columns[i]; j <= columns[i+1]; j++) {
      saveResidues[j] = -1;
      delAminos++;
    }
  }

  /* The method computes the new alignment columns' number */
  newResidNumber = residNumber - delAminos;

  /* Once we've selected the columns, if the complementary flag is true, we will have to change the selected and
    non-selected columns. */
  if(complementary == true) {
    newResidNumber = residNumber - newResidNumber;
    for(i = 0; i < residNumber; i++) {
      if(saveResidues[i] == -1)
        saveResidues[i] = i;
      else
        saveResidues[i] = -1;
    }
  }

  /* We allocate memory to save the columns selected. */
  matrixAux = new char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i] = new char[newResidNumber+1];

  /* Copy the columns from the original alignment to the new alignment, only if the column has been selected. */
  for(i = 0, k = 0; i < residNumber; i++)
    if(saveResidues[i] != -1) {
      for(j = 0; j < sequenNumber; j++)
        matrixAux[j][k] = alignmentMatrix[j][i];
      k++;
    }

  /* Concatenate the character '\0' to the last position of all rows of the new alignment to facilitate character 
     vector manipulation and printing functions. */
  for(i = 0; i < sequenNumber; i++)
    matrixAux[i][newResidNumber] = '\0';

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(inputFileName, alignmentInfo, matrixAux, sequenNames, seqInfo, sequenNumber, newResidNumber, iformat, oformat,
                          dataType, isAligned, sequenNumber, residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow);
 
  /* Deallocated auxiliar memory */
  for(i = 0; i < sequenNumber; i++)
    delete[] matrixAux[i]; 
  delete[] matrixAux;

  /* Return the new alignment reference */
  return newAlig;
}

int *alignment::getCorrespResidues(void) {
  return saveResidues;
}

int *alignment::getCorrespSequences(void) {
  return saveSequences;
}

bool alignment::isFileAligned(void) {
  return isAligned;
}

