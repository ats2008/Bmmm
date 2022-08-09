import numpy as np
import types
from scipy import stats
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from itertools import product, combinations

import ROOT
ROOT.gSystem.Load('libBmmmAnalysis')
from ROOT import KVFitter   # VertexDistance3D is contained here, dirt trick!!
from ROOT import Utilities
cpp_code = """
// Function definition
reco::Track* getTrackFromRef(reco::GsfTrack *trkRf) { return new reco::Track(*trkRf); }

"""
ROOT.gInterpreter.ProcessLine(cpp_code)

# make these available everywhere in here
global vtxfit
vtxfit = KVFitter()

global utilities
utilities= Utilities()

global tofit
tofit = ROOT.std.vector('reco::Track')()

def candTrack(self):
    return self.track_candidate

class BToEMuCandidate():
    '''
    E-Mu candidate.
    MISSING: use the post fit muon momenta
    '''
    def __init__(self, mu , ele , vertices, beamspot):
        # sort by pt
        self.mu =  mu
        self.ele = ele
        self.ele.track_candidate = ROOT.reco.Track(ele.gsfTrack().get())

        self.ele.candTrack =  types.MethodType( candTrack , self.ele)
        # check that the muon track covariance matrix is pos-def
        self.mu.cov = self.convert_cov(self.mu.bestTrack().covariance())
        self.ele.cov = self.convert_cov(self.ele.candTrack().covariance())
        self.mu.is_cov_pos_def = self.is_pos_def(self.mu.cov)
        self.ele.is_cov_pos_def = self.is_pos_def(self.ele.cov)
        # choose as PV the one that's closest to the leading muon in the dz parameter
        self.pv = sorted( [vtx for vtx in vertices], key = lambda vtx : abs( self.mu.bestTrack().dz(vtx.position() ) ) )[0]
        dZs = [ abs( self.mu.bestTrack().dz(vtx.position() ) ) for vtx in vertices ]
        self.pvIndex = np.argsort(dZs)[0]
        # create a Vertex type of object from the bs coordinates at the z of the chosen PV
        bs_point = ROOT.reco.Vertex.Point(
            beamspot.x(self.pv.z()),
            beamspot.y(self.pv.z()),
            beamspot.z0(),
        )

        bs_error = beamspot.covariance3D()
        chi2 = 0.
        ndof = 0.
        self.bs = ROOT.reco.Vertex(bs_point, bs_error, chi2, ndof, 3) # size? say 3? does it matter?

        # we'll fit a vertex out of the three muons, shall we? 
        # ideally this can be triggered on demand, and just build a skinny candidate to 
        # check simple things, such as mass etc
        #tofitTrkRefs.clear()

        tofit.clear()
        tofit.push_back(self.mu.bestTrack())
        tofit.push_back(self.ele.candTrack())
        
        self.vtx = vtxfit.Fit(tofit)
        #self.vtx = ele
        #self.vtx.isValid()=False
        self.vtx.chi2 = self.vtx.normalisedChiSquared() if self.vtx.isValid() else np.nan
        self.vtx.prob = (1. - stats.chi2.cdf(self.vtx.chi2, 1)) if self.vtx.isValid() else np.nan 

        # now compute some displacement related quantities, here in the transverse plane.
        # later can add 3D quantities
        self.lxy = ROOT.VertexDistanceXY().distance(self.bs, self.vtx.vertexState()) if self.vtx.isValid() else np.nan

        vect_lxy = ROOT.Math.DisplacementVector3D('ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag')( 
                    self.vtx.position().x() - self.bs.position().x(),
                    self.vtx.position().y() - self.bs.position().y(),
                    0. ) if self.vtx.isValid() else np.nan

        vect_pt = ROOT.Math.DisplacementVector3D('ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag')( 
                    self.px(),
                    self.py(),
                    0. ) if self.vtx.isValid() else np.nan

        self.vtx.cos = vect_pt.Dot(vect_lxy) / (vect_pt.R() * vect_lxy.R()) if (self.vtx.isValid() and vect_lxy.R() > 0.) else np.nan
        
        self.pv_to_sv = ROOT.Math.XYZVector(
                            (self.vtx.position().x() - self.pv.position().x()), 
                            (self.vtx.position().y() - self.pv.position().y()),
                            (self.vtx.position().z() - self.pv.position().z())
                        ) if self.vtx.isValid() else np.nan
        self.Bdirection  = self.pv_to_sv/np.sqrt(self.pv_to_sv.Mag2()) if self.vtx.isValid() else np.nan                  
        self.Bdir_eta    = self.Bdirection.eta() if self.vtx.isValid() else np.nan                                
        self.Bdir_phi    = self.Bdirection.phi() if self.vtx.isValid() else np.nan                                
        self.emu_p4_par  = self.p4().Vect().Dot(self.Bdirection) if self.vtx.isValid() else np.nan                   
        self.emu_p4_perp = np.sqrt(self.p4().Vect().Mag2() - self.emu_p4_par*self.emu_p4_par) if self.vtx.isValid() else np.nan
        self.mcorr       = np.sqrt(self.p4().mass()*self.p4().mass() + self.emu_p4_perp*self.emu_p4_perp) + self.emu_p4_perp if self.vtx.isValid() else np.nan


    def computeTrkMuonIsolation( self,pfHandle ):
        return utilities.computeTrkMuonIsolation(self.mu , self.ele , pfHandle , self.pvIndex ) 

    def convert_cov(self, m):
        return np.array([[m(i,j) for j in range(m.kCols)] for i in range(m.kRows)])

    def is_pos_def(self, x):
        '''
        https://stackoverflow.com/questions/16266720/find-out-if-matrix-is-positive-definite-with-numpy
        '''
        return np.all(np.linalg.eigvals(x) > 0)

    def p4(self):
        return self.mu.p4() + self.ele.p4(1)
    def pt(self):
        return self.p4().pt()
    def eta(self):
        return self.p4().eta()
    def phi(self):
        return self.p4().phi()
    def mass(self):
        return self.p4().mass()
    def mass_corrected(self):    
        return self.mcorr
    def energy(self):
        return self.p4().energy()
    def px(self):
        return self.p4().px()
    def py(self):
        return self.p4().py()
    def pz(self):
        return self.p4().pz()
    def charge(self):
        return self.mu.charge() + self.ele.charge()
    def r(self):
        '''
        Cone radius parameter: max distance between the 4-mu candidate direction and one of the muons
        '''
        return max([deltaR(self.p4(), self.mu) ,deltaR(self.p4(), self.ele)])

    def max_dr(self):
        '''
        Max distance between pairwise muons
        '''
        return self.drMuE()
    def drMuE(self):
        return deltaR(self.mu, self.ele)
    def mass12(self):
        return self.p4.mass()
    def __str__(self):
        to_return = [
            'cand mass %.2f pt %.2f eta %.2f phi %.2f' %(self.mass(), self.pt(), self.eta(), self.phi()),
            'cand vtx prob %2f vtx chi2 %.2f lxy %.4f lxy sig %.2f cos %.2f' %(self.vtx.prob, self.vtx.chi2, self.lxy.value(), self.lxy.significance(), self.vtx.cos),
            '\t mu pt %.2f eta %.2f phi %.2f' %(self.mu.pt(), self.mu.eta(), self.mu.phi()),
            '\t ele pt %.2f eta %.2f phi %.2f' %(self.ele.pt(), self.ele.eta(), self.ele.phi()),
        ]
        return '\n'.join(to_return)
