#include "Bmmm/Analysis/interface/Utilities.h"

float Utilities::computeTrkMuonIsolation(const pat::Muon& the_muon, 
                      const pat::Electron& the_electron, 
                      std::vector<pat::PackedCandidate> pfCandidates,
					  unsigned int primaryVertexIndex,
					  float minPt, float dR,
					  std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  float sumPt(0);
  for (const auto& pfCand: pfCandidates){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(the_muon, pfCand) < 0.01 || deltaR(the_electron, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if (deltaR(the_muon, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return the_muon.pt()/(the_muon.pt()+sumPt);
}

