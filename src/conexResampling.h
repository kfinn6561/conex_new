#ifndef _include_conex_resampling_h_
#define _include_conex_resampling_h_

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include <conexConfig.h>
#include <ResamplingMode.h>

#include <fstream>
class TFile;
class TTree;


#ifdef CONEX_EXTENSIONS

// RU: to change the CrossSection by a given factor
extern double gFactorCrossSection;
extern double gFactorExtraMeson;
extern double gFactorResampling;
extern double gResamplingThreshold;


extern resample::EResampleMode gResamplingMode;


//KF: Declare Classicalization variables
extern bool gClassicalizationFlag; //flag to see if classicalization occurred
extern double gClassicalizationThreshold; //classicalization energy threshold
extern double gClassicalizationFraction; //fraction of energy in classicalization event
extern double gNscaling;//classicalization number scaling
extern double gClassicalonMass;//Mass of classical quanta
extern double gClasigma;//Clasicalization cross section
extern double gSigma0;//regular cross section
extern bool gClassicalizationOff;

double ProbDistDat[1001];
bool ReadData=true;
// ------------------------------------------------
//  function definitions
//

class CommonBlockCONEX;
class CommonBlockSIBYLL;
class CommonBlockSIBYLL_LAB;

extern "C" {

  void addclassicalization_(CommonBlockCONEX& blockPtr, double* pfive, int& primaryId, double& Mtarg);//KF
  
  void resampleconex_(CommonBlockCONEX& p, double& x, int& primaryId);
  void resamplesibyll_(CommonBlockSIBYLL& p, double& x, int& primaryId,
		       float& sqs, double& plab);
  void resamplesibylllab_(CommonBlockCONEX& p, CommonBlockSIBYLL_LAB& s,
			  double& x, int& primaryId);
  //				    float& sqs, double& plab);}
  
  void modifiercx_(double& factMod, const double& energy, const int& pid);
  void sibyllf19_(double& factModPP, const double& factModPair);
  void classicalcx_(double& factMod, const double& energy, const int& pid, const double& sigma);

  void reducee_(double& Elab, double& ECM, double& Engy, double& Pmod, double& Ekin, double& Mproj, double& Mtarg);//KF
  void checkclassicalization_(double& dz);//KF

}

extern "C" {
  struct xschadron {
    double xsamproj;
    double xsamtarg;
    double xsypjtl;
    double xsyhaha;
    double xspnullx;
  };
}


// ---------------------------------------------
//  debugging
//
#ifdef DEBUG_RESAMPLING
class TClonesArray;
std::string gROOTfilename = "";
// extern "C" {void cxutlob5_(double& yboost, double& x1, double& x2, double& x3, double& x4,  double& x5);}
TFile* fResamplingDebug = 0;
TTree* tResamplingDebug = 0;
int resampling_primary = 0;
double resampling_Elab = 0;
double resampling_factor = 0;
double resampling_yboost = 0;
int resampling_inucleon = 0;
int resampling_nnucleon = 0;
int resampling_cms = false;
TClonesArray* resampling_before = 0;
TClonesArray* resampling_after = 0;
int resampling_before_n = 0;
int resampling_after_n = 0;
#endif // DEBUG_RESAMPLING
#endif



#endif
