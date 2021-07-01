// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include<map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class CosmicAnalysis : public edm::EDAnalyzer {
public:
  explicit CosmicAnalysis(const edm::ParameterSet&);
  ~CosmicAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::Service<TFileService> fs;
  
  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 

  int nEvents;
  int nGEMRecHits;
  int nCSCSegments;
  int nME11Segments;
  int nME11Segments6;
  int nInGEMBounds;
  TTree *t_csc;

  float csc_r, csc_x, csc_y, csc_z, csc_pr, csc_px, csc_py, csc_pz, csc_pth, chi2;
  float gem_x, gem_y, gem_z;
  float dr,in;
  int re, ch;
};

CosmicAnalysis::CosmicAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nGEMRecHits(0),
  nCSCSegments(0),
  nME11Segments(0),
  nME11Segments6(0),
  nInGEMBounds(0)
{
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  cscSegments_ = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  // vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());
  
  t_csc = fs->make<TTree>("csc", "csc");
  t_csc->Branch("r", &csc_r, "r/F");
  t_csc->Branch("x", &csc_x, "x/F");
  t_csc->Branch("y", &csc_y, "y/F");
  t_csc->Branch("z", &csc_z, "z/F");
  t_csc->Branch("pr", &csc_pr, "pr/F");
  t_csc->Branch("pth", &csc_pth, "pth/F");
  t_csc->Branch("px", &csc_px, "px/F");
  t_csc->Branch("py", &csc_py, "py/F");
  t_csc->Branch("pz", &csc_pz, "pz/F");
  t_csc->Branch("chi2", &chi2, "chi2/F");
  t_csc->Branch("gem_x", &gem_x, "gem_x/F");
  t_csc->Branch("gem_y", &gem_y, "gem_y/F");
  t_csc->Branch("gem_z", &gem_z, "gem_z/F");
  t_csc->Branch("dr", &dr, "dr/F");
  t_csc->Branch("in", &in, "in/F");
  t_csc->Branch("re", &re, "re/I");
  t_csc->Branch("ch", &ch, "ch/I");
}

CosmicAnalysis::~CosmicAnalysis()
{
  std::cout << "::: GEM Slice Test Results :::" << std::endl;
  std::cout << ": From " << nEvents << " events" << std::endl;
  std::cout << std::endl;
  std::cout << ": # GEMs " << nGEMRecHits << "" << std::endl;
  std::cout << ": # CSCS " << nCSCSegments << "" << std::endl;
  std::cout << ": # ME11 " << nME11Segments << "" << std::endl;
  std::cout << ": # ME11 w 6 hits " << nME11Segments6 << "" << std::endl;
  std::cout << ": # ME11 w 6 hits in GEM bounds " << nInGEMBounds << "" << std::endl;
  // std::cout << " # Muons " << nMuonTotal << std::endl;
  // std::cout << " # FidMu " << nGEMFiducialMuon << std::endl;
  // std::cout << " # GEMMu " << nGEMTrackWithMuon << std::endl;
  std::cout << std::endl;
}

void
CosmicAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nEvents++;

  edm::ESHandle<GEMGeometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const GEMGeometry* GEMGeometry_ = &*hGeom;
  edm::ESHandle<CSCGeometry> hCSCGeom;
  iSetup.get<MuonGeometryRecord>().get(hCSCGeom);
  const CSCGeometry* CSCGeometry_ = &*hCSCGeom;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  // iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagator_);
  // iSetup.get<IdealMagneticFieldRecord>().get(bField_); 
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegments_, cscSegments);

  // std::cout << "NCSC : " << cscSegments->size() << std::endl;
  
  // edm::Handle<reco::VertexCollection> vertexCollection;
  // iEvent.getByToken( vertexCollection_, vertexCollection );
  // if(vertexCollection.isValid()) {
  //   vertexCollection->size();
  //   //    std::cout << "vertex->size() " << vertexCollection->size() <<std::endl;
  // }

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);
  for (auto& gem : *gemRecHits) {
    nGEMRecHits++;
  }

  for (auto& csc : *cscSegments) {
    nCSCSegments++;
    if (csc.cscDetId().station()==1) {
      nME11Segments++;
      if (csc.recHits().size() == 6) {
        nME11Segments6++;
        chi2 = csc.chi2()/csc.degreesOfFreedom();
        auto chamber = CSCGeometry_->chamber(csc.cscDetId());
        auto pos = chamber->toGlobal(csc.localPosition());
        csc_r = pos.mag();
        csc_x = pos.x();
        csc_y = pos.y();
        csc_z = pos.z();
        auto dir = chamber->toGlobal(csc.localDirection());
        re = csc.cscDetId().zendcap();
        ch = csc.cscDetId().chamber();
        csc_pr = dir.mag();
        csc_pth = dir.theta();
        csc_px = dir.x();
        csc_py = dir.y();
        csc_pz = dir.z();
        //   constexpr GEMDetId(int region, int ring, int station, int layer, int chamber, int ieta)

        // Take eta partition as the chamber plane is slightly offset
        // from the readout plane. All parts should be on the same z
        // plane though
        GEMDetId gemDetId(csc.cscDetId().zendcap(),1,1,2,csc.cscDetId().chamber(),1);
        auto gem_part = GEMGeometry_->etaPartition(gemDetId);
        gem_z = gem_part->toGlobal(LocalPoint(0,0)).z();

        // find intercept of straight line from csc to gem chamber
        float s = (gem_z - pos.z()) / dir.z();
        gem_x = pos.x() + s * dir.x();
        gem_y = pos.y() + s * dir.y();
        // check if it propagates inside the bounds of one of the eta partitions
        in=0;
        for (int i = 1; i <= 8; ++i) {
          GEMDetId gemDetId(csc.cscDetId().zendcap(),1,1,2,csc.cscDetId().chamber(),i);
          auto gem_part = GEMGeometry_->etaPartition(gemDetId);
          auto loc = gem_part->toLocal(GlobalPoint(gem_x,gem_y,gem_z));
          if (gem_part->surface().bounds().inside(loc)) {
            in=1;
            nInGEMBounds++;
          }
        }

        dr=999.;
        for (int i = 1; i <= 8; ++i) {
          GEMDetId gemDetId(csc.cscDetId().zendcap(),1,1,2,csc.cscDetId().chamber(),i);
          auto gem_part = GEMGeometry_->etaPartition(gemDetId);
          auto recHitsRange = gemRecHits->get(gemDetId); 
          auto gemRecHit = recHitsRange.first;
          for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
            auto hit_pos = gem_part->toGlobal(hit->localPosition());
            float this_dr = sqrt(pow(gem_x-hit_pos.x(),2)+pow(gem_y-hit_pos.y(),2)+pow(gem_z-hit_pos.z(),2));
            if (this_dr < dr) dr = this_dr;
          }
        }
        
        t_csc->Fill();
      }
    }
  }
}

void CosmicAnalysis::beginJob(){}
void CosmicAnalysis::endJob(){}

void CosmicAnalysis::beginRun(Run const& run, EventSetup const&){
}
void CosmicAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(CosmicAnalysis);
