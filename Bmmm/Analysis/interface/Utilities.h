
#include "FWCore/Framework/interface/ESHandle.h"


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


class Utilities {

    private:

    public:
          float  computeTrkMuonIsolation(const pat::Muon& muon, 
			  const pat::Electron& the_electron,
              //edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_,
              std::vector<pat::PackedCandidate> pfCandHandle_,
			  unsigned int primaryVertexIndex,
			  float minPt=0.5, float dR=0.5,
			  std::vector<const pat::PackedCandidate*> ignoreTracks = 
			  std::vector<const pat::PackedCandidate*>());
};
