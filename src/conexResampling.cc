#include <conexResampling.h>
#include <conexConfig.h>

#include <Modifier.h>

#include <Verbosity.h>
#include <CommonBlockWrapperCONEX.h>
#include <CommonBlockParticleCONEX.h>
#include <CommonBlockWrapperSIBYLL.h>
#include <CommonBlockParticleSIBYLL.h>
#include <CommonBlockWrapperSIBYLL_LAB.h>
#include <resample.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <map>
#include <sstream>

#include <TDirectory.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>

#include <TRandom3.h>//KF: for random numbers
#include <LambertW.h>//KF: for the Lambert W function
#include <fstream>//KF: to read in data

using namespace std;
using namespace resample;


#ifdef CONEX_EXTENSIONS

void 
resampleconex_(CommonBlockCONEX& blockPtr, double& Elab, int& primaryId) 
{
  if (ruVerbosity>2) {
    cout << " resampleconex_ Elab=" << Elab << endl;
  }
  
  /*
    for(int i=0; i<10; i++) {
    cout << " " << *((double*)(&blockPtr)+i);
    }
    cout << endl;
    
    for(int i=-10; i<10; i++) {
    cout << " " << i << "=" << *(&(blockPtr.number)+i);
    }
    cout << endl;
  */
  
  CommonBlockWrapperCONEX pblock(blockPtr, false);
  
  /*
    const double factorResampling = (Elab>10.e6 ? 
    1.0 + (gFactorResampling-1.0) * (log10(Elab)-6.)/4. :
    1.0);
  */
  const double factorResampling = Modifier(gFactorResampling, 
					   gResamplingThreshold,
					   Elab);   
  
#ifdef DEBUG_RESAMPLING
  resampling_primary = primaryId;
  resampling_Elab = Elab;
  resampling_factor = factorResampling;
  resampling_yboost = xschadron_.xsyhaha;
  resampling_cms = pblock.IsCenterOfMassSystem();
  resampling_nnucleon = 1;
  resampling_inucleon = 0;
  
  resampling_before_n = 0;
  resampling_after_n = 0;
  resampling_before->Clear();
  resampling_after->Clear();
  
  const int nBefore = pblock.Size();
  for (int iParticle=0; iParticle<nBefore; ++iParticle) {
    const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
    if (p.IsGood()) {  // These are real particles, all the rest is crap
      const int id = p.GetId();
      const string name = p.GetName();
      double px = p.GetPx();
      double py = p.GetPy();
      double pz = p.GetPz();
      double E = p.GetEnergy();
      double M = p.GetMass();
      new((*resampling_before)[resampling_before_n++]) particle(id, name, px, py, pz, E, M);
    }
  } // loop particles
#endif // DEBUG_RESAMPLING

      
   // ## DEBUG ############################################
   if (ruVerbosity>3) {
     cout << " CONEX-STACK BEFORE RESAMPLING " << endl;
     const int n = pblock.Size();
     for (int iParticle=0; iParticle<n; ++iParticle) {
       CommonBlockParticleCONEX p(&pblock.GetData(), iParticle);//, false);
       if (p.IsGood()) { // These are real particles, all the rest is crap
	 cout << " index=" << setw(3) << iParticle
	      << " id=" << setw(5) << p.GetId() 
	      << " (" << setw(11) << p.GetName() << ")"
	      << " status=" << setw(2) << p.GetStatus()
	      << " E=" << setw(12) << p.GetMomentum(3)
	      << " M=" << setw(12) << setprecision(6) << p.GetMomentum(4)
	      << " p=" << setprecision(2) 
	      << sqrt(pow(p.GetPx(),2)+pow(p.GetPy(),2)+pow(p.GetPz(),2))
	   //<< " " << p.GetPx() << " " << p.GetPy() << " " << p.GetPz()
	   //<< " " << setprecision(4) << sqrt((p.GetMomentum(3)+p.GetMomentum(4))*(p.GetMomentum(3)-p.GetMomentum(4)))
	      << endl;
       }
     } // loop particles
     cout << " -------- END CONEX-STACK ------- " << endl;
   }
   // ## end-DEBUG ############################################

     
   // ## DEBUG ############################################
   if (ruVerbosity>9) {
     int nEntries = blockPtr.GetNumberEntries();
     cout << " ---> resample> nSecondaries=" << nEntries 
	  << " Elab=" << Elab
	  << " primaryId=" << primaryId
	  << " factorResampling=" << factorResampling 
	  << endl;
   }
   // ## end-DEBUG ############################################
   
   if (factorResampling==1.0) {

     // ## DEBUG ############################################
     if (ruVerbosity>1) {
       cout << "resample-CONEX: No need for resampling. Skipping." << endl;
     }
     // ## end-DEBUG ############################################

   } else {
     
     // ## DEBUG ############################################
     if (ruVerbosity>=0) {
       cout << " ********* do resampling (CONEX) > nStack=" << pblock.Size() 
	    << " primaryId=" << primaryId
	    << " Elab=" << Elab
	    << " factor=" << factorResampling
	    << " ********* "
	    << endl;
     }
     // ## end-DEBUG ############################################
     
     resample::SecondaryResampling(gResamplingMode, factorResampling, primaryId, pblock);
     
     // ## DEBUG ############################################
     if (ruVerbosity>3) {
       cout << " CONEX-STACK AFTER RESAMPLING " << endl;
       const int n = pblock.Size();
       for (int iParticle=0; iParticle<n; ++iParticle) {
	 CommonBlockParticleCONEX p(&pblock.GetData(), iParticle);//, false);
	 if (p.IsGood()) { // These are real particles, all the rest is crap
	   cout << " index=" << setw(3) << iParticle
		<< " id=" << setw(5) << p.GetId()
		<< " (" << setw(11) << p.GetName() << ")"
		<< " status=" << setw(2) << p.GetStatus()
		<< " E=" << setw(12) << p.GetEnergy()
		<< " M=" << setw(12) << setprecision(6) << p.GetMass()
		<< " p=" << setprecision(2) 
		<< sqrt(pow(p.GetPx(),2)+pow(p.GetPy(),2)+pow(p.GetPz(),2))
		<< endl;
	 }
       } // loop particles
       cout << " -------- END CONEX-STACK ------- " << endl;
     }
     // ## end-DEBUG ############################################
     
#ifdef DEBUG_RESAMPLING
     const int nAfter = pblock.Size();
     for (int iParticle=0; iParticle<nAfter; ++iParticle) {
       const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
       if (p.IsGood()) {  // These are real particles, all the rest is crap
	 const int id = p.GetId();
	 const string name = p.GetName();
	 double px = p.GetPx();
	 double py = p.GetPy();
	 double pz = p.GetPz();
	 double E = p.GetEnergy();
	 double M = p.GetMass();
	 new((*resampling_after)[resampling_after_n++]) particle(id, name, px, py, pz, E, M);
       }
     } // loop particles
     // FILL TREE
     tResamplingDebug->Fill();
#endif // DEBUG_RESAMPLING
 
   }
   
}

//KF function to reduce energy if classicalization event
void
reducee_(double& Elab, double& ECM, double& Engy, double& Pmod, double& Ekin, double& Mproj, double& Mtarg)
{
  //cout << "Entered ClassicalizationReduceEnergy"<<endl;//KF:debug
  //cout << "starting energy is " << Elab << endl;//KF:debug
  if (gClassicalizationFlag){
    //cout << "Classicalization event, reducing energy" << endl;//KF:debug


    /*
    ECM=ECM*(1-gClassicalizationFraction);
    Engy=ECM;
       
    Elab=(ECM*ECM-Mtarg*Mtarg-Mproj*Mproj)/(2*Mtarg);
    Pmod=sqrt(max(0.0,Elab*Elab-Mproj*Mproj));
    Ekin=Elab-Mtarg;//Make other variables fit the reduction in CM energy
    */
    double N=(ECM*gClassicalizationFraction/gClassicalonMass)*utl::LambertW<0>((ECM*gClassicalizationFraction*gClassicalonMass)/(gClassicalizationThreshold*gClassicalizationThreshold))*gNscaling;
    
    
    if (N>2){//N<2 should have been caught earlier, but just in case
    Elab=Elab*(1-gClassicalizationFraction);
    Mtarg=Mtarg*(1-gClassicalizationFraction);
    Mproj=Mproj*(1-gClassicalizationFraction);

    ECM=sqrt(max(0.0,2*Elab*Mtarg+Mtarg*Mtarg+Mproj*Mproj));
    Engy=ECM;
    Pmod=sqrt(max(0.0,Elab*Elab-Mproj*Mproj)); 
    }else{
      cout<<"ERROR: classicalized state would contain "<<N<<" particles. Reverting to regular collision"<<endl;
      gClassicalizationFlag=false;
    }


    //NOTE: MAY NOT NEED TO DO THIS SINCE NOT ACTUALLY A REDUCED ENERGY EVENT WITH SAME MASS
  }



  
}

/*
  Direct resampling of SIBYLL secondaries (CM-System!). Not used for
  CONEX resampling ...
*/
void resamplesibyll_(CommonBlockSIBYLL& blockPtr, double& Elab, int& primaryId,
		     float& SQS, double& Mprojectile) { // FOR DEBUG

  if (ruVerbosity>2) {
    cout << " resamplesibyll_ Elab=" << Elab << endl;
  }
  
  CommonBlockWrapperSIBYLL pblock(blockPtr, true);
  
  
  // TODO: check whether we do have to use Elab/A for this type of resampling ...
  const double factorResampling = Modifier(gFactorResampling,
					   gResamplingThreshold,
					   Elab);
  
  
#ifdef DEBUG_RESAMPLING
  resampling_primary = primaryId;
  resampling_Elab = Elab;
  resampling_factor = factorResampling;
  resampling_yboost = xschadron_.xsyhaha;
  resampling_cms = pblock.IsCenterOfMassSystem();
  resampling_nnucleon = 1;
  resampling_inucleon = 0;
  
  resampling_before_n = 0;
  resampling_after_n = 0;
  resampling_before->Clear();
  resampling_after->Clear();
  
  const int nBefore = pblock.Size();
  for (int iParticle=0; iParticle<nBefore; ++iParticle) {
    const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
    if (p.IsGood()) {  // These are real particles, all the rest is crap
      const int id = p.GetId();
      const string name = p.GetName();
      double px = p.GetPx();
      double py = p.GetPy();
      double pz = p.GetPz();
      double E = p.GetEnergy();
      double M = p.GetMass();
      new((*resampling_before)[resampling_before_n++]) particle(id, name, px, py, pz, E, M);
    }
  } // loop particles
#endif // DEBUG_RESAMPLING
  
  
  
  // ## DEBUG ############################################
  if (ruVerbosity>3) {
    cout << " SIBYLL-STACK BEFORE RESAMPLING " << endl;
    //const double gamma = Elab/SQS;
    //const double betagamma = sqrt(Elab*Elab-Mprojectile*Mprojectile)/SQS;
    //cout << " gamma=" << gamma << " betagamma=" << betagamma << endl;
    const int n = pblock.Size();
    for (int iParticle=0; iParticle<n; ++iParticle) {
      CommonBlockParticleSIBYLL p(&pblock.GetData(), iParticle);//, false);
      if (p.IsGood()) { // These are real particles, all the rest is crap
	
	const double E = p.GetEnergy();
	const double P = sqrt(pow(p.GetPx(),2)+pow(p.GetPy(),2)+pow(p.GetPz(),2));
	//const double Elab = gamma*E + betagamma*P;
	
	cout << " index=" << setw(3) << iParticle
	     << " id=" << setw(5) << p.GetId()
	     << " (" << setw(11) << p.GetName() << ")"
	  //<< " status=" << setw(2) << p.GetStatus()
	     << " E_cm=" << setw(12) << E // p.GetMomentum(3)
	     << " E_lab=" << setw(12) << Elab
	     << " m=" << setw(12) << setprecision(6) << p.GetMass()
	     << " p=" << setprecision(2) << P
	     << endl;
      }
    } // loop particles
    cout << " -------- END SIBYLL-STACK ------- " << endl;
  }
  // ## end-DEBUG ############################################
  
  
  // ## DEBUG ############################################
  if (ruVerbosity>9) {
    int nEntries = blockPtr.GetNumberEntries();
    cout << " ---> resample> nSecondaries=" << nEntries 
	 << " Elab=" << Elab
	 << " primaryId=" << primaryId
	 << " factorResampling=" << factorResampling 
	 << endl;
  }
  // ## end-DEBUG ############################################
  
  if (factorResampling==1.0) {
    
    // ## DEBUG ############################################
    if (ruVerbosity>1) {
      cout << "resample-SIBYLL: No need for resampling. Skipping." << endl;
    }
    // ## end-DEBUG ############################################
    
  } else {
    
    // ## DEBUG ############################################
    if (ruVerbosity>=0) {
      cout << " ********* do resampling (SIBYLL) > nStack=" << pblock.Size() 
	   << " primaryId=" << primaryId
	   << " Elab=" << Elab
	   << " factor=" << factorResampling
	   << " ********* "
	   << endl;
    }
    // ## end-DEBUG ############################################
    
    resample::SecondaryResampling(gResamplingMode, factorResampling, primaryId, pblock);
    
    // ## DEBUG ############################################
    if (ruVerbosity>3) {
      cout << " SIBYLL-STACK AFTER RESAMPLING " << endl;
      const int n = pblock.Size();
      for (int iParticle=0; iParticle<n; ++iParticle) {
	const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
	if (p.IsGood()) { // These are real particles, all the rest is crap
	  cout << " index=" << setw(3) << iParticle
	       << " id=" << setw(5) << p.GetId()
	       << " (" << setw(11) << p.GetName() << ")"
	    //<< " status=" << setw(2) << p.GetStatus()
	       << " E=" << setw(12) << p.GetEnergy()
	       << " m=" << setw(12) << setprecision(6) << p.GetMass()
	    //<< " p=" << setprecision(2) 
	    //<< sqrt(pow(p.GetMomentum(0),2)+pow(p.GetMomentum(1),2)+pow(p.GetMomentum(2),2))
	       << endl;
	}
      } // loop particles
      cout << " -------- END SIBYLL-STACK ------- " << endl;
    }
    // ## end-DEBUG ############################################
    
#ifdef DEBUG_RESAMPLING
    const int nAfter = pblock.Size();
    for (int iParticle=0; iParticle<nAfter; ++iParticle) {
      const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
      if (p.IsGood()) {  // These are real particles, all the rest is crap
	const int id = p.GetId();
	const string name = p.GetName();
	double px = p.GetPx();
	double py = p.GetPy();
	double pz = p.GetPz();
	double E = p.GetEnergy();
	double M = p.GetMass();
	new((*resampling_after)[resampling_after_n++]) particle(id, name, px, py, pz, E, M);
      }
    } // loop particles
    // FILL TREE
    tResamplingDebug->Fill();
#endif // DEBUG_RESAMPLING
    
  }   
  
} // end SIBYLL resampling


/*
  Resample SIBYLL interaction, after transformation into LAB system
  The SIBYLL mapping of individual nucleon-nucleon interactions is
  preserved, and thus allows a unified handling of p-air and A-air
  interactions.
*/
void resamplesibylllab_(CommonBlockCONEX& conexPtr, 
			CommonBlockSIBYLL_LAB& sibyllPtr,
			double& Elab, int& primaryId) {
  //			 float& SQS, double& Mprojectile) { // FOR DEBUG

  if (ruVerbosity>2) {
    cout << " resamplesibylllab_ Elab=" << Elab << " id=" << primaryId
	 << " nint=" << sibyllPtr.GetNumberOfInteractions()
	 << endl;
  }
  
  const int nInteraction = sibyllPtr.GetNumberOfInteractions();
  const double factorResampling = Modifier(gFactorResampling,
					   gResamplingThreshold,
					   Elab);
  
  if (factorResampling==1.0) {
    
    // ## DEBUG ############################################
    if (ruVerbosity>1) {
      cout << "resample-SIBYLL-LAB: No need for resampling. Skipping." << endl;
    }
    // ## end-DEBUG ############################################
    
  } else { // factorResampling==1.0
    
    for (int iInteraction=0; 
	 iInteraction<nInteraction;
	 ++iInteraction) {
      
      if (ruVerbosity>1) {
	cout << " *********  Resampling nucleus-air interaction, nucleon number " << iInteraction+1
	     << " of " << sibyllPtr.GetNumberOfInteractions()
	     << " ********* "
	     << endl;
      }
      
      CommonBlockWrapperSIBYLL_LAB pblock(conexPtr,
					  sibyllPtr,
					  iInteraction,
					  false);     
      
#ifdef DEBUG_RESAMPLING
      
      resampling_primary = primaryId;
      resampling_Elab = Elab;
      resampling_factor = factorResampling;
      resampling_yboost = xschadron_.xsyhaha;
      resampling_cms = pblock.IsCenterOfMassSystem();
      resampling_nnucleon = sibyllPtr.GetNumberOfInteractions();
      resampling_inucleon = iInteraction;
      
      resampling_before_n = 0;
      resampling_after_n = 0;
      resampling_before->Clear();
      resampling_after->Clear();
      
      const int nBefore = pblock.Size();
      for (int iParticle=0; iParticle<nBefore; ++iParticle) {
	const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
	if (p.IsGood()) {  // These are real particles, all the rest is crap
	  const int id = p.GetId();
	  const string name = p.GetName();
	  double px = p.GetPx();
	  double py = p.GetPy();
	  double pz = p.GetPz();
	  double E = p.GetEnergy();
	  double M = p.GetMass();
	  new((*resampling_before)[resampling_before_n++]) particle(id, name, px, py, pz, E, M);
	}
      } // loop particles
#endif // DEBUG_RESAMPLING
      
      
       // ## DEBUG ############################################
      if (ruVerbosity>3) {
	cout << " STACK BEFORE RESAMPLING " << endl;
	const int n = pblock.Size();
	for (int iParticle=0; iParticle<n; ++iParticle) {
	  const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
	   if (p.IsGood()) { // These are real particles, all the rest is crap
	     const double E = p.GetEnergy();
	     //const double P = sqrt(pow(p.GetMomentum(0),2)+pow(p.GetMomentum(1),2)+pow(p.GetMomentum(2),2));
	     
	     cout << " index=" << setw(3) << iParticle
		  << " id=" << setw(5) << p.GetId()
		  << " (" << setw(11) << p.GetName() << ")"
	       //<< " status=" << setw(2) << p.GetStatus()
		  << " E=" << setw(12) << E // p.GetMomentum(3)
	       //<< " E_lab=" << setw(12) << Elab
		  << " m=" << setw(12) << setprecision(6) << p.GetMass()
	       //<< " p=" << setprecision(2) << P
		  << endl;
	   }
	} // loop particles
	cout << " -------- END STACK ------- " << endl;
      }
      // ## end-DEBUG ############################################
      
      // ## DEBUG ############################################
      if (ruVerbosity>9) {
	const int nEntriesCONEX = conexPtr.GetNumberEntries();
	const int nEntriesNUC = sibyllPtr.GetNumberEntries(iInteraction);
	cout << " ---> resample> nSecondaries=" << nEntriesNUC 
	     << " (conex-stack: " << nEntriesCONEX << ")"
	      << " Elab=" << Elab
	     << " primaryId=" << primaryId
	     << " factorResampling=" << factorResampling 
	     << endl;
      }
      // ## end-DEBUG ############################################
      
      // ## DEBUG ############################################
      if (ruVerbosity>=0) {
	cout << " ********* do resampling (SIBYLL-LAB) > nStack=" << pblock.Size() 
	     << " primaryId=" << primaryId
	     << " Elab=" << Elab
	     << " factor=" << factorResampling
	      << " ********* "
	     << endl;
      }
      // ## end-DEBUG ############################################
      
      resample::SecondaryResampling(gResamplingMode, factorResampling, primaryId, pblock);
      
      // ## DEBUG ############################################
      if (ruVerbosity>3) {
	 cout << " STACK AFTER RESAMPLING " << endl;
	 const int n = pblock.Size();
	 for (int iParticle=0; iParticle<n; ++iParticle) {
	   const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
	   if (p.IsGood()) { // These are real particles, all the rest is crap
	     cout << " index=" << setw(3) << iParticle
		  << " p=" << setw(5) << p.GetId()
		  << " (" << setw(11) << p.GetName() << ")"
		  << " E=" << setw(12) << p.GetEnergy()
		  << " m=" << setw(12) << setprecision(6) << p.GetMass()
	       //<< " p=" << setprecision(2) 
	       // << sqrt(pow(p.GetMomentum(0),2)+pow(p.GetMomentum(1),2)+pow(p.GetMomentum(2),2))
		  << endl;
	   }
	 } // loop particles
	 cout << " -------- END STACK ------- " << endl;
      }
      // ## end-DEBUG ############################################
      
#ifdef DEBUG_RESAMPLING
      const int nAfter = pblock.Size();
      for (int iParticle=0; iParticle<nAfter; ++iParticle) {
	const ParticleBlockEntry& p = pblock.GetEntry(iParticle);
	if (p.IsGood()) {  // These are real particles, all the rest is crap
	  const int id = p.GetId();
	  const string name = p.GetName();
	  double px = p.GetPx();
	  double py = p.GetPy();
	  double pz = p.GetPz();
	  double E = p.GetEnergy();
	  double M = p.GetMass();
	  new((*resampling_after)[resampling_after_n++]) particle(id, name, px, py, pz, E, M);
	}
      } // loop particles
      // FILL TREE
      tResamplingDebug->Fill();
#endif // DEBUG_RESAMPLING
      
    } // loop interactions (nucleon-air)
    
  } // if factorResampling!=1.0
  
} // end resamplesibylllab_ 
 

void
modifiercx_(double& factMod, const double& energy, const int& pid)
{
  double f19 = gFactorCrossSection;
  if (pid>=2 && pid<=6) { // pions, kaons
    f19 = f19 + (f19-1) * (gFactorExtraMeson-1);
    // f19 = f19 * gFactorExtraMeson - (gFactorExtraMeson-1);
    }
  factMod = Modifier(f19,
		     gResamplingThreshold,
		     energy);
  //cout << " gFactorCrossSection=" << gFactorCrossSection << ", gFactorExtraMeson=" << gFactorExtraMeson  << ", e=" << energy << ", f=" << factMod << endl;
}


void
addclassicalization_(CommonBlockCONEX& blockPtr, double* pfive, int& primaryId, double& Mtarg)
{
  CommonBlockWrapperCONEX pblock(blockPtr, false);
  //cout<<"pfive: "<<pfive[0]<<", "<<pfive[1]<<", "<<pfive[2]<<", "<<pfive[3]<<", "<<pfive[4]<<endl;//KF:debug
  //cout<<"pmass: "<<sqrt(pfive[3]*pfive[3]-pfive[1]*pfive[1]-pfive[2]*pfive[2]-pfive[0]*pfive[0])<<endl;//KF:debug
  if (gClassicalizationFlag)
    {
      //cout<<"resampling classicalization event\n";//KF:debug
      resample::ClassicalizationResampling(pfive,primaryId, pblock, Mtarg);//KF: resamplesibyll
    }else{
    //cout<<"checked for classicalization but found none\n";//KF:debug
  }
    
  /*KF: This was to deal with Sybill. I think it will be more trouble than it's worth to add it.
  if (gClassicalizationFlag)
    {
      cout<<"resampling classicalization event\n";//KF:debug

      for (int iInteraction=0;
	   iInteraction<nInteraction;
	   ++iInteraction) {

	CommonBlockWrapperSIBYLL_LAB pblock(conexPtr,
					    sibyllPtr,
					    iInteraction,
					    false);
	cout<<"Elab="<<Elab<<", primaryId="<<primaryId<<endl;//KF:debug
	resample::ClassicalizationResampling(Elab,primaryId, pblock);//KF: resamplesibylllab
      }//KF:loop over interactions with nucleus THIS MAY REQUIRE MORE THOUGHT
    }else{
    cout<<"K checked for classicalization but found none\n";
    }*/

}

void
ReadProbData()
{
  //KF: read in probability data
  std::ifstream probfile("overlap_data.dat");
  std::string temp;
  if (probfile.is_open())
    {
      int i=0;
      while (!probfile.eof())
	{
	  getline(probfile,temp);
	  ProbDistDat[i]=atof(temp.c_str());
	  //cout<<"data "<<i<<": "<<temp<<endl;
	  i++;
	}
    }
  ReadData=false;//set flag to prevent data being read in more than once
}


//KF: Calculate the fractional energy that classicalises
double
GetFraction(double ranNo)
{
  if (ReadData)
    {
      ReadProbData();
    }
  int i=(int)(ranNo*1000);//Assumes the probability data is sampled at 1001 evenly spaced out places between 0 and 1. This is quicker, but less robust
  double x1=i*0.001;
  double x2=(i+1)*0.001;
  double y1=ProbDistDat[i];
  double y2=ProbDistDat[i+1];
  return ((y1-y2)/(x1-x2))*(ranNo-x1)+y1;//Linear interpolation
}


//KF: decide whether this is a classicalization event and alter the cross section accordingly. This sets gClassicalizationFlag
void
classicalcx_(double& factMod, const double& energy, const int& pid, const double& sigma)
{
  if (gClassicalizationOff){//look for gClassicalizationOff flag and do nothing if raised
    factMod=1.0;
    gClassicalizationFlag=false;//This should be enough to prevent any other code performing classicalization
    return;
  }
  
  const double mtarg=0.94;
  const double mproj=0.94;//KF: assume both projectile and target are protons mass=0.94 hardcoded, may want to update if important
  //gClassicalizationFraction=gRandom->Uniform();//choose fraction of energy to classicalize. TODO this distribution may need to change. currently uniform
  gClassicalizationFraction=GetFraction(gRandom->Uniform());//using overlap of two spheres
  gClassicalizationFraction=0.9;//const fraction
  const double comEnergy=gClassicalizationFraction*sqrt(2*mtarg*energy+mtarg*mtarg+mproj*mproj);//KF: assume target is proton mass=0.94 hardcoded, may want to update if important
  double N=(comEnergy/gClassicalonMass)*utl::LambertW<0>((comEnergy*gClassicalonMass)/(gClassicalizationThreshold*gClassicalizationThreshold))*gNscaling;
  //N=3;//rm_cl

  if (N<2){
    //below classicalization threshold
    gClasigma=0.;
  }else{
    gClasigma = M_PI*pow((1./gClassicalonMass)*utl::LambertW<0>((comEnergy*gClassicalonMass)/(gClassicalizationThreshold*gClassicalizationThreshold)),2.0)*0.3893793;//number converts between 1/Gev^2 and mb
  }
  
  gSigma0=sigma;
  
  factMod=(sigma+gClasigma)/sigma;
  
  }

void
checkclassicalzation_(double& dz)
{
  if (gClassicalizationOff){//No classicalization
    gClassicalizationFlag=false;
    return;
  }
  //cout<<std::scientific;//KF:debug
  //cout<<"\n\nentered checkclassicalization"<<endl;//KF:debug
  //cout<<"avog: "<<avog<<endl;//KF:debug
    double sigEffective=(gClasigma-gSigma0)//Effective cross section in cm^2/g (avog already has conversion from mb to cm^2)
    double pClassicalize=1./((gSigma0/gClasigma)*exp(sigEffective*dz)+1);//Probability to classicalize
    //cout<<"Probability to classicalize: "<<pClassicalize<<endl;//KF:debug

    double r=gRandom->Uniform();

    if (r<pClassicalize){
      gClassicalizationFlag=true;
    }else{
      gClassicalizationFlag=false;
    }
    
}


void 
sibyllf19_(double& factModPP, const double& factModPair) 
{
  // output of: https://devel-ik.fzk.de/svn/users/ralf/trunk/glauber/SibyllPAir2PP_f19
  static const int nData = 73;
  static double dataX[nData]={0.0288897,0.031692,0.0347569,0.0381071,0.0417671,0.0457627,0.0501217,0.0548732,0.0600482,0.0656791,0.0717999,0.0784458,0.0856532,0.0934595,0.101902,0.11102,0.120851,0.131431,0.142797,0.154984,0.168023,0.181942,0.196767,0.212518,0.229213,0.246863,0.265474,0.28505,0.305589,0.327091,0.34955,0.372966,0.397342,0.422689,0.449029,0.476398,0.504849,0.534456,0.565316,0.597549,0.631302,0.666746,0.704078,0.743521,0.785322,0.829755,0.877116,0.927732,0.981955,1.04017,1.10281,1.17032,1.24321,1.32203,1.4074,1.49998,1.60052,1.70982,1.82879,1.95842,2.09981,2.25415,2.42279,2.60719,2.80896,3.02988,3.27192,3.53722,3.82816,4.14735,4.49766,4.88225,5.30458};
  static double dataY[nData]={0.01,0.011,0.0121,0.01331,0.014641,0.0161051,0.0177156,0.0194872,0.0214359,0.0235795,0.0259374,0.0285312,0.0313843,0.0345227,0.037975,0.0417725,0.0459497,0.0505447,0.0555992,0.0611591,0.067275,0.0740025,0.0814027,0.089543,0.0984973,0.108347,0.119182,0.1311,0.14421,0.158631,0.174494,0.191943,0.211138,0.232252,0.255477,0.281024,0.309127,0.340039,0.374043,0.411448,0.452593,0.497852,0.547637,0.602401,0.662641,0.728905,0.801795,0.881975,0.970172,1.06719,1.17391,1.2913,1.42043,1.56247,1.71872,1.89059,2.07965,2.28762,2.51638,2.76801,3.04482,3.3493,3.68423,4.05265,4.45792,4.90371,5.39408,5.93349,6.52683,7.17952,7.89747,8.68722,9.55594};
  static TGraph* gFactPair = new TGraph(nData, dataX, dataY);
  if (factModPair<dataX[0] || factModPair>dataX[nData-1]) {
    cerr << "conex.cc::sibyllf19 ERROR factModPair=" << factModPair << " out of range " << dataX[0] << " - " << dataX[nData-1] << endl;
    exit(4);
  }
  factModPP = gFactPair->Eval(factModPair);
  //cout << " sibyllf19 = " << factModPair << " " << factModPP << endl;
}

#endif // CONEX_EXTENSIONS


