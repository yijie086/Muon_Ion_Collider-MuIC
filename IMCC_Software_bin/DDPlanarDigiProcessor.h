#ifndef DDPlanarDigiProcessor_h
#define DDPlanarDigiProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
using namespace lcio ;
using namespace marlin ;

namespace EVENT {
  class SimTrackerHit;
}


/** ======= DDPlanarDigiProcessor ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions. 
 * The geometry of the surface is retreived from DDRec::Surface associated to the hit via cellID.
 * 
 * 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
 * (default value false)
 * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
 * (default value false) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
 * (default value lcio::ILDDetID::VXD) <br>
 * <br>
 * 
 * @author F.Gaede CERN/DESY, S. Aplin DESY
 * @date Dec 2014
 */
class DDPlanarDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new DDPlanarDigiProcessor ; }
  
  
  DDPlanarDigiProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  

  
protected:
  
  std::string _inColName ;
  
  std::string _outColName ;
  std::string _outRelColName ;
 
  std::string _subDetName ;
  
  int _nRun ;
  int _nEvt ;
  
  FloatVec _resU ;
  FloatVec _resV ;
  FloatVec _resT ;
  
  bool _isStrip;
  
  gsl_rng* _rng ;
  
  const dd4hep::rec::SurfaceMap* _map ;

  bool _forceHitsOntoSurface  ;
  double _minEnergy ;

  bool _useTimeWindow ;
  bool _correctTimesForPropagation ;
  FloatVec _timeWindow_min ;
  FloatVec _timeWindow_max ;

  std::vector<TH1F*> _h ;
  std::vector<TH2F*> _h2d;
  std::vector<TH3F*> _h3d;
} ;

#endif



