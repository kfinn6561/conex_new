#ifndef _include_ConexDynamicInterface_h_
#define _include_ConexDynamicInterface_h_

#include <conexConfig.h>
#include <conexHEModels.h>

#include <string>

class ConexDynamicInterface {

  ConexDynamicInterface();
  ConexDynamicInterface(const ConexDynamicInterface&);
  ConexDynamicInterface& operator=(const ConexDynamicInterface&);

public:
  ConexDynamicInterface(const EHEModel model);
  ~ConexDynamicInterface();

  EHEModel GetHEModel() const { return fModel; }
  std::string GetHEModelName() const;

  // generic library symbol fetcher
  void* GetSymbol(const std::string& name) const;

  // function pointers to CONEX subroutines
  void (*InitConex)(int&, int*, EHEModel&, int&, 
#ifdef CONEX_EXTENSIONS
		    int&,
#endif
		    const char*, int);
  void (*RunConex)(int*, double&, double&, double&, double&, int&);

  void (*GetHeaderVariables)(int&, float&, float&, int&, float&, int&, float&,  
                             int&, float&, int& , float&, float&, float&,
                             float&, float&, int&, float&);
#ifdef __MC3D__
  void (*GetHeaderVariables3D)(float&, float&, float&, float&, float&,
                               float&, float&, int&, float&, float&);
  void (*Get3DOutput)(int&, float&, float&, int&, float&, float&, float&);
#endif

  void (*GetShowerData)(int&, int&, int&, float&, float&,
			float&, float&, float&);
  void (*GetdEdXProfile)(int&, int&, float&, float&);
  void (*GetMuonProfile)(int&, int&, float&, float&);
  void (*GetGammaProfile)(int&, int&, float&);
  void (*GetElectronProfile)(int &, int&, float&);
  void (*GetHadronProfile)(int &, int&, float&);
  int (*GetNumberOfDepthBins)();

  double (*ConexRandom)(double& dummy);

  double (*InteractionLength)(int&, double&, double&);
  void (*ParticleMass)(int&, double&);

#ifdef __SIBYLL21__
  // prepare call to initialize (modified) nucleus-air cross sections (VERY SLOW)
  void (*DoSignucIniIni)(const int& init);
#endif      

#ifdef CONEX_EXTENSIONS
  void (*initListTree)(const int i, const int mode);
  void (*finishListTree)(const int mode); 
  void (*presetSeed)(int seed[3]);
#endif
  
private:
  void initFunctionPointers(EHEModel model);


 protected:
  EHEModel fModel;

 private:
  void* fConexLib;
};

#endif
