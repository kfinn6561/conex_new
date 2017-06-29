#include <ConexDynamicInterface.h>
#include <CxOutputPrefix.h>
#include <conexConfig.h>

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <dlfcn.h>

#include <cmath>//KF: for pow

using namespace std;

class TTree;

#ifdef CONEX_EXTENSIONS
// the full relative path is needed here to allow compilation
// of ConexShowerGenerator module of Auger Offline
#include <../resampler/ResamplingMode.h> 

// RU: to change the CrossSection by a given factor, etc.
double gFactorCrossSection = 1;
double gFactorExtraMeson = 1;
double gFactorResampling = 1;
double gResamplingThreshold = 15-9;

//KF: Classicalization
double gClassicalizationThreshold= 15.;//GeV
bool gClassicalizationFlag=false;
double gClassicalizationFraction=0.0;
double gFixedFraction=0.9;
double gNscaling=1.0;
double gClassicalonMass=0.17;
double gClasigma=0.;
double gSigma0=1.;
double gForwardThreshold=0.1;
bool gClassicalizationOff=false;
bool gFinalState=true;
bool gForward=false;

resample::EResampleMode gResamplingMode = resample::eNone;

int gParticleListMode = 0;
std::string gParticleListModeFile = "unknown.dat";

#endif // extensions



ConexDynamicInterface::ConexDynamicInterface(const EHEModel model)
: fModel(model),
  fConexLib(0) {
  initFunctionPointers(model);
}


ConexDynamicInterface::~ConexDynamicInterface()
{
  if (fConexLib)
    dlclose(fConexLib);
}


std::string
ConexDynamicInterface::GetHEModelName()
const
{
  switch(fModel) {
  case eSibyll23:
    return "Sibyll2.3";
    break;
  case eQGSJet01:
    return "QGSJet01";
    break;
  case eQGSJetII:
    return "QGSJetII-04";
    break;
  case eEposLHC:
    return "EposLHC";
    break;
  default:
    return "unknonw";
  }
}


/*************************************************************
 *
 *************************************************************/
void
ConexDynamicInterface::initFunctionPointers(EHEModel model)
{
  const string conex_prefix = ""_CONEX_PREFIX;
  const string conex_system = ""_CONEX_SYSTEM;

  string sharedLibExt = ".so";
  if ( conex_system == "Darwin" )
    sharedLibExt = ".dynlib";

  ostringstream conexLibNameS;
  conexLibNameS << conex_prefix << "/lib/" << conex_system << "/";
  if ( model == eQGSJet01)
    conexLibNameS << "libCONEXqgsjet" << sharedLibExt;
  if ( model == eSibyll23)
    conexLibNameS << "libCONEXsibyll" << sharedLibExt;
  else if ( model == eQGSJetII )
    conexLibNameS << "libCONEXqgsjetII" << sharedLibExt;
  else if ( model == eEposLHC )
    conexLibNameS << "libCONEXepos" << sharedLibExt;

  const string conexLibName = conexLibNameS.str();

  cout << fgOutputPrefix << "loading library"
       << "\n" << setw(fgOutputPrefix.size()) << "" << "\'"
       << conexLibName << "\'" << endl;

  fConexLib = dlopen(conexLibName.c_str(), RTLD_LAZY);

  if ( ! fConexLib ) {
    ostringstream errMsg;
    errMsg << "\n can not open CONEX shared library " << conexLibName
	   << " at path \'" << conex_prefix << "/lib/" << conex_system << "\'\n\n"
	   << " Dynamic-link error:\n \"" << dlerror() << "\"\n";

    cout << errMsg.str() << endl;
    exit(1);
  }

  bool failed = false;

  InitConex =  (void(*)(int&, int*, EHEModel&, int&, 
#ifdef CONEX_EXTENSIONS
			int&,
#endif
			const char*, int)) dlsym(fConexLib, "initconex_");
  if (!InitConex) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function initconex_"
	 << endl;
    failed = true;
  }

  RunConex  = (void (*)(int*, double&, double&, double&, double&, int&)) dlsym(fConexLib, "runconex_");
  if (!RunConex) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function runconex_" << endl;
    failed = true;
  }

  GetHeaderVariables = (void (*)(int&, float&, float&, int&, float&, int&, float&,  
                                 int&, float&, int& , float&, float&, float&,
                                 float&, float&, int&, float&)) dlsym(fConexLib, "get_data_card1_");
  if (!GetHeaderVariables) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_data_card1_"
	 << endl;
    failed = true;
  }

#if __MC3D__
  GetHeaderVariables3D = (void (*)(float&, float&, float&, float&, float&,
			     float&, float&, int&, float&, float&)) dlsym(fConexLib, "get_data_card2_");
  if (!GetHeaderVariables3D) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_data_card2_"
	 << endl;
    failed = true;
  }
  
  Get3DOutput = (void (*)(int&, float&, float&, int&, float&, float&, float&)) dlsym(fConexLib, "get_shower_lat_");
  if (!Get3DOutput) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_lat_"
	 << endl;
    failed = true;
  }
  
#endif

  GetShowerData =  (void (*)(int&, int&, int&, float&, float&,
                             float&, float&, float&)) dlsym(fConexLib, "get_shower_data_");
  if (!GetShowerData) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_data_"
	 << endl;
    failed = true;
  }

  GetdEdXProfile = (void (*)(int&, int&, float&, float&)) dlsym(fConexLib, "get_shower_edep_");
  if (!GetdEdXProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_edep_"
	 << endl;
    failed = true;
  }

  GetMuonProfile = (void (*)(int&, int&, float&, float&))dlsym(fConexLib, "get_shower_muon_");
  if (!GetMuonProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_muon_"
	 << endl;
    failed = true;
  }

  GetGammaProfile = (void (*)(int&, int&, float&)) dlsym(fConexLib, "get_shower_gamma_");
  if (!GetGammaProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_gamma_"
	 << endl;
    failed = true;
  }

  GetElectronProfile = (void (*)(int &, int&, float&)) dlsym(fConexLib, "get_shower_electron_");
  if (!GetElectronProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_electron_"
	 << endl;
    failed = true;
  }

  GetHadronProfile = (void (*)(int &, int&, float&)) dlsym(fConexLib, "get_shower_hadron_");
  if (!GetHadronProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_hadron_"
	 << endl;
    failed = true;
  }

  GetNumberOfDepthBins = (int (*)()) dlsym(fConexLib, "get_number_of_depth_bins_");
  if (!GetNumberOfDepthBins) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_number_of_depth_bins_"
	 << endl;
    failed = true;
  }

  ConexRandom = (double (*)(double&)) dlsym(fConexLib, "drangen_");
  if (!ConexRandom) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function drangen_" << endl;
    failed = true;
  }

  InteractionLength = (double(*)(int&,double&,double&)) dlsym(fConexLib, "rlam_");
  if (!InteractionLength) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function rlam_"
	 << endl;
    failed = true;
  }

  ParticleMass = (void(*)(int&,double&)) dlsym(fConexLib, "cxidmass_");
  if (!ParticleMass) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function cxidmass_"
	 << endl;
    failed = true;
  }

#ifdef __SIBYLL21__
  DoSignucIniIni = (void(*)(const int& init)) dlsym(fConexLib, "do_signuc_ini_ini_");
  if (!DoSignucIniIni) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function do_signuc_ini_ini_"
	 << endl;
    failed = true;
  }
#endif

  
#ifdef CONEX_EXTENSIONS
  initListTree = (void(*)(const int, const int)) dlsym(fConexLib, "initListTree");
  if (!initListTree) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function initListTree"
	 << endl;
    failed = true;
  }
  
  finishListTree = (void(*)(const int)) dlsym(fConexLib, "finishListTree");
  if (!finishListTree) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function finishListTree"
	 << endl;
    failed = true;
  }
  
  presetSeed = (void(*)(int seed[3])) dlsym(fConexLib, "presetSeed");
  if (!presetSeed) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function presetSeed"
	 << endl;
    failed = true;
  }
#endif

  
  if (failed) {
    cout << fgOutputPrefix << " ConexDynamicInterface failed!" << endl;
    exit(1);
  }
}

void* 
ConexDynamicInterface::GetSymbol(const string& name) 
  const
{
  void* ptr = (void*) dlsym(fConexLib, name.c_str());
  if (!ptr) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX symbol: \'" 
	 << name << "\";" 
	 << endl;
  }
  return ptr;
}
