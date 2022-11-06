set ExecutionPath {
	ParticlePropagator
	
	GenJetFinder
	
	ChargedHadronTrackingEfficiency
	ElectronTrackingEfficiency
	MuonTrackingEfficiency
	ForwardMuonEfficiency
	
	ChargedHadronSmearing
	ElectronSmearing
	MuonSmearing
	ForwardMuonSmearing
	
	TrackMerger
	
	ECal
	HCal
	
	Calorimeter
	
	PIDSystems
	
	EFlowMerger
	
	MissingET
	
	FastJetFinder
	JetEnergyScale
	
	BTagging
	
	PhotonEfficiency
	PhotonIsolation
	
	ElectronFilter
	ElectronEfficiency
	ElectronIsolation
	
	MuonEfficiency
	MuonIsolation
	
	UniqueObjectFinder
	
	TreeWriter
}

module ParticlePropagator ParticlePropagator {
	set InputArray Delphes/stableParticles
	set OutputArray stableParticles
	set ChargedHadronOutputArray chargedHadrons
	set ElectronOutputArray electrons
	
	set Radius 1.5
	set HalfLength 2.10
	set Bz 3.0
}

module FastJetFinder GenJetFinder {
	set InputArray  ParticlePropagator/stableParticles
	set OutputArray jets
	
	# algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
	set JetAlgorithm 6
	set ParameterR 1.0
	set JetPTMin 4.0
}

set CommonTrackingEfficiency {
		(eta > -5.0 && eta < 2.4) * (pt > 0.250) * (1.0) +
		(eta > 2.4) * (0.0) +
		(eta < -5.0) * (0.00)+
		0.0
}

module Efficiency ChargedHadronTrackingEfficiency {
	set InputArray ParticlePropagator/chargedHadrons
	set OutputArray chargedHadrons
	set EfficiencyFormula $CommonTrackingEfficiency
}
module Efficiency ElectronTrackingEfficiency {
	set InputArray ParticlePropagator/electrons
	set OutputArray electrons
	set EfficiencyFormula $CommonTrackingEfficiency
}
module Efficiency MuonTrackingEfficiency {
	set InputArray ParticlePropagator/muons
	set OutputArray muons
	set EfficiencyFormula $CommonTrackingEfficiency
}
module Efficiency ForwardMuonEfficiency {
		set InputArray ParticlePropagator/muons
		set OutputArray muons
	
		set EfficiencyFormula {
			(pt < 1.0) * (0.0) +
			(eta < -5.0 && eta > -8.0) * (pt >= 0.5) * (0.80) +
			(eta < -8.0 ) * (0.00)
		}
}


set CommonTrackingResolution {
		(eta > -5.0 && eta < 2.4) * (pt > 0.250) * sqrt(0.01^2 + (pt*cosh(eta))^2*(1e-3)^2 )
}

module MomentumSmearing ChargedHadronSmearing {
	set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
	set OutputArray chargedHadrons
	set ResolutionFormula $CommonTrackingResolution
}
module MomentumSmearing ElectronSmearing {
	set InputArray ElectronTrackingEfficiency/electrons
	set OutputArray electrons
	set ResolutionFormula $CommonTrackingResolution
}
module MomentumSmearing MuonSmearing {
	set InputArray MuonTrackingEfficiency/muons
	set OutputArray muons
	set ResolutionFormula $CommonTrackingResolution
}
module MomentumSmearing ForwardMuonSmearing {
		set InputArray ForwardMuonEfficiency/muons
		set OutputArray fmuons
		set ResolutionFormula {
			(eta < -5.0 && eta > -8.0) * (0.1)
	}
}

module Merger TrackMerger {
# add InputArray InputArray
	add InputArray ChargedHadronSmearing/chargedHadrons
	add InputArray ElectronSmearing/electrons
	add InputArray MuonSmearing/muons
	set OutputArray tracks
}

module SimpleCalorimeter ECal {
	set ParticleInputArray ParticlePropagator/stableParticles
	set TrackInputArray TrackMerger/tracks
	
	set TowerOutputArray ecalTowers
	set EFlowTrackOutputArray eflowTracks
	set EFlowTowerOutputArray eflowPhotons
	
	set IsEcal true
	set EnergyMin 0.050
	#does not seem possible to set minimum dependent on eta as spec in the YR.
	
	
	set EnergySignificanceMin 1.0
	
	set SmearTowerCenter true
	
	set pi [expr {acos(-1)}]
	
	# lists of the edges of each tower in eta and phi
	# each list starts with the lower edge of the first tower
	# the list ends with the higher edged of the last tower
	
	# Granularity is not discussed in EIC detector handbook.
	##BARREL
	#assume 0.1 x 0.1 (real cell size will be smaller, so this is to represent some cluster)
	
		set PhiBins {}
		for {set i -135} {$i <=135} {incr i} {
			add PhiBins [expr {$i * $pi/135.0}]
		}
		for {set i -211} {$i <=101} {incr i} {
			set eta [expr {$i * 0.0236}]
			add EtaPhiBins $eta $PhiBins
		}
	
	
	add EnergyFraction {0} {0.0}
	# energy fractions for e, gamma and pi0
	add EnergyFraction {11} {1.0}
	add EnergyFraction {22} {1.0}
	add EnergyFraction {111} {1.0}
	# energy fractions for muon, neutrinos and neutralinos
	add EnergyFraction {12} {0.0}
	add EnergyFraction {13} {0.0}
	add EnergyFraction {14} {0.0}
	add EnergyFraction {16} {0.0}
	add EnergyFraction {1000022} {0.0}
	add EnergyFraction {1000023} {0.0}
	add EnergyFraction {1000025} {0.0}
	add EnergyFraction {1000035} {0.0}
	add EnergyFraction {1000045} {0.0}
	# energy fractions for K0short and Lambda
	add EnergyFraction {310} {0.3}
	add EnergyFraction {3122} {0.3}
	
	set ResolutionFormula {
		(eta <= 2.4 && eta > -5.0) *	 sqrt((0.1)^2 * energy + 0.02^2 * energy^2)
	}
	
}

module SimpleCalorimeter HCal {
	set ParticleInputArray ParticlePropagator/stableParticles
	set TrackInputArray ECal/eflowTracks
	
	set TowerOutputArray hcalTowers
	set EFlowTrackOutputArray eflowTracks
	set EFlowTowerOutputArray eflowNeutralHadrons
	
	set IsEcal false
	
	##Assumes noise 100 MeV per tower.
	set EnergyMin 0.5
	set EnergySignificanceMin 1.0
	
	set SmearTowerCenter true
	
	set pi [expr {acos(-1)}]
	
		set PhiBins {}
		for {set i -135} {$i <=135} {incr i} {
			add PhiBins [expr {$i * $pi/135.0}]
		}
		for {set i -211} {$i <=101} {incr i} {
			set eta [expr {$i * 0.0236}]
			add EtaPhiBins $eta $PhiBins
		}
	
	
	add EnergyFraction {0} {1.0}
	# energy fractions for e, gamma and pi0
	add EnergyFraction {11} {0.0}
	add EnergyFraction {22} {0.0}
	add EnergyFraction {111} {0.0}
	# energy fractions for muon, neutrinos and neutralinos
	add EnergyFraction {12} {0.0}
	add EnergyFraction {13} {0.0}
	add EnergyFraction {14} {0.0}
	add EnergyFraction {16} {0.0}
	add EnergyFraction {1000022} {0.0}
	add EnergyFraction {1000023} {0.0}
	add EnergyFraction {1000025} {0.0}
	add EnergyFraction {1000035} {0.0}
	add EnergyFraction {1000045} {0.0}
	# energy fractions for K0short and Lambda
	add EnergyFraction {310} {0.7}
	add EnergyFraction {3122} {0.7}
	
	# set HCalResolutionFormula {resolution formula as a function of eta and energy}
	set ResolutionFormula {
		(eta <= 2.4 && eta > -5.0) *	 sqrt((0.5)^2 * energy + 0.1^2 * energy^2)
	}
	
}

module Merger Calorimeter {
# add InputArray InputArray
	add InputArray ECal/ecalTowers
	add InputArray HCal/hcalTowers
	add InputArray MuonSmearing/muons
	set OutputArray towers
}

module IdentificationMap PIDSystems {
		set InputArray HCal/eflowTracks
		set OutputArray tracks
		
		add EfficiencyFormula {-11} {-11} {1.0}
		add EfficiencyFormula {211} {211} {1.0}
		# pi/K/proton identification, assuming 3 sigma separation between all species! (all pair-wise separations are 3 sigma)
		# Kaons must have p > 0.135 GeV/c
		# Pions must have p > 0.100 GeV/c
		add EfficiencyFormula {321} {321} {1.0}
		add EfficiencyFormula {211} {211} {1.0}
		add EfficiencyFormula {2212} {2212} {1.0}
		# Everything else with no PID system coverage is 100% identified as pion
		add EfficiencyFormula {0} {0} {0.00}
}

module Merger EFlowMerger {
		# add InputArray InputArray
		add InputArray HCal/eflowTracks
		add InputArray ECal/eflowPhotons
		add InputArray HCal/eflowNeutralHadrons
		set OutputArray eflow
}

module Merger MissingET {
		# add InputArray InputArray
		add InputArray EFlowMerger/eflow
		set MomentumOutputArray momentum
}

module FastJetFinder FastJetFinder {
	set InputArray EFlowMerger/eflow
	set OutputArray jets
	
	# algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
	set JetAlgorithm 6
	set ParameterR 1.0
	set JetPTMin 20.0
}

module BTagging BTagging {
	set JetInputArray FastJetFinder/jets
	set PartonInputArray Delphes/partons
	
	set BitNumber 0
	
	# based on ATL-PHYS-PUB-2015-022
	
	# default efficiency formula (misidentification rate)
	add EfficiencyFormula {0} {0.002+7.3e-06*pt}
	
	# efficiency formula for c-jets (misidentification rate)
	add EfficiencyFormula {4} {0.20*tanh(0.02*pt)*(1/(1+0.0034*pt))}
	
	# efficiency formula for b-jets
	add EfficiencyFormula {5} {0.80*tanh(0.003*pt)*(30/(1+0.086*pt))}
}

module EnergyScale JetEnergyScale {
	set InputArray FastJetFinder/jets
	set OutputArray jets
	
	# scale formula for jets
	set ScaleFormula {1.00}
}

module Efficiency PhotonEfficiency {
		set InputArray ECal/eflowPhotons
		set OutputArray photons
		set EfficiencyFormula {1}
}
module Isolation PhotonIsolation {
		set CandidateInputArray PhotonEfficiency/photons
		set IsolationInputArray EFlowMerger/eflow
		set OutputArray photons
		
		set DeltaRMax 0.1
		set PTMin 0.5
		set PTRatioMax 0.2
}

module PdgCodeFilter ElectronFilter {
		set InputArray HCal/eflowTracks
		set OutputArray electrons
		set Invert true
		add PdgCode {11}
		add PdgCode {-11}
}
module Efficiency ElectronEfficiency {
	set InputArray ElectronFilter/electrons
	set OutputArray electrons
	set EfficiencyFormula {1}
}
module Isolation ElectronIsolation {
	set CandidateInputArray ElectronEfficiency/electrons
	set IsolationInputArray EFlowMerger/eflow
	set OutputArray electrons
	
	set DeltaRMax 0.1
	set PTMin 0.5
	set PTRatioMax 0.2
}

module Efficiency MuonEfficiency {
		set InputArray MuonSmearing/muons
		set OutputArray muons

		set EfficiencyFormula {1}
}
module Isolation MuonIsolation {
		set CandidateInputArray MuonEfficiency/muons
		set IsolationInputArray EFlowMerger/eflow
		set OutputArray muons
		
		set DeltaRMax 0.1
		set PTMin 0.5
		set PTRatioMax 0.2
}

module UniqueObjectFinder UniqueObjectFinder {
		add InputArray PhotonIsolation/photons photons
		add InputArray ElectronIsolation/electrons electrons
		add InputArray MuonIsolation/muons muons
		add InputArray EFlowMerger/eflow eflow
		add InputArray JetEnergyScale/jets jets
}

module TreeWriter TreeWriter {
	add Branch Calorimeter/towers Tower Tower
	
	add Branch TrackMerger/tracks Track Track
	add Branch ForwardMuonSmearing/fmuons forwardmuon Track
	
	add Branch HCal/hcalTowers HCalTower Tower
	add Branch ECal/ecalTowers ECalTower Tower
	
	add Branch UniqueObjectFinder/photons Photon Photon
	add Branch UniqueObjectFinder/electrons Electron Electron
	add Branch UniqueObjectFinder/muons Muon Muon
	
	add Branch FastJetFinder/jets Jet Jet
	add Branch GenJetFinder/jets GenJet Jet
	
	add Branch PIDSystems/tracks PIDSystemsTrack Track
	add Branch MissingET/momentum MissingET MissingET
}
