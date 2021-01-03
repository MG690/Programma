#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"

int main() {
TFile *file = new TFile("data.root", "recreate");

/*gROOT->LoadMacro("ParticleType.cpp++");
gROOT->LoadMacro("ResonanceType.cpp++");
gROOT->LoadMacro("Particle.cpp++");*/

int nGen = 100000;
char p_plus[] = "pion + ";
char p_minus[] = "pion - ";
char K_plus[] = "Kaone + ";
char K_minus[] = "Kaone- ";
char P_plus[] = "Proton + ";
char P_minus[] = "Proton- ";
char K_star[] = "Resonance ";
Particle::AddParticleType(p_plus, 0.13957, 1, 0);
Particle::AddParticleType(p_minus, 0.13957, -1, 0);
Particle::AddParticleType(K_plus, 0.49367, 1, 0);
Particle::AddParticleType(K_minus, 0.49367, -1, 0);
Particle::AddParticleType(P_plus, 0.93827, 1, 0);
Particle::AddParticleType(P_minus, 0.93827, -1, 0);
Particle::AddParticleType(K_star, 0.89166, 0, 0.050);

	// Histogram of the distribution of the types of generated particles (h1)
	TH1F *h1 = new TH1F("h1 ", "distribution of the types of generated particles ", 7, 0, 7);
	// Histogram of the distribution of azimutal angles (h2)
	TH1F *h2 = new TH1F("h2 ", "distribution of azimutal angles ", 1000, 0, 2 * M_PI);
	// Histogram of the distribution of polar angles (h3)
	TH1F *h3 = new TH1F("h3 ", "distribution of polar angles ", 500, 0, 2 * M_PI);
	// Histogram of the distribution of the momentum (h4)
	TH1F *h4 = new TH1F("h4 ", "distribution of the momentum ", 700, 0, 7);
	// Histogram of the distribution of the trasverse momentum (h5)
	TH1F *h5 = new TH1F("h5 ", "distribution of the trasverse momentum ", 600, 0, 6);
	// Histogram of the distribution of the energy (h6)
	TH1F *h6 = new TH1F("h6 ", "distribution of the energy ", 600, 0, 6);
	// Histogram of the distribution of the invariant mass generated in every event (h7)
	TH1F *h7 = new TH1F("h7 ", "distribution of the invariant mass ", 100, 0, 5);
	// Histogram of the invariant mass generated in every event of particle with opposite sign (h8)
	TH1F *h8 = new TH1F("h8 ", "distribution of the invariant mass of particles with discordant sign ", 100, 0, 5);
	// Histogram of the invariant mass generated in every event of particle with the same sign (h9)
	TH1F *h9 = new TH1F("h9 ", "distribution of the invariant mass with concordant sign ", 100, 0, 5);
 	// Histogram of the invariant mass generated in every event with combination of p+/K- o p-/K+ (h10)
 	TH1F *h10 = new TH1F("h10 ", "distribution of the invariant mass of p+ / K- , p- / K+ ", 100, 0, 5);
 	// Histogram of the invariant mass generated in every event with combination of p+/K+ o p-/K- (h11)
 	TH1F *h11 = new TH1F("h11 ", "distribution of the invariant mass of p+ / K+ , p- / K- ", 100, 0, 5);
	// Histogram of the invariant mass generated in every event which derive from decayment of K* (h12)
	TH1F *h12 = new TH1F("h12 ", "distribution of the invariant mass of the resonances ", 30, 0, 2);

gRandom -> SetSeed();

for(int i = 0; i < nGen; ++i) {
		
	if (i == 0) { std::cout << "running... " << '\n';}
	int generatedParticles = 100;
	int decayedParticles = 0;
	int N = 120;
	Particle particles[N];		
		
	for (int n = 0; n < generatedParticles; ++n) { // Particle particle {};
		double phi = gRandom->Uniform(2 * M_PI);
		h2 -> Fill(phi);
		double theta = gRandom->Uniform(M_PI);
		h3 -> Fill(theta);
		double TotalMomentum = gRandom->Exp(1);
		h4 -> Fill(TotalMomentum);
		double fPx = TotalMomentum*sin(theta)*cos(phi);
		double fPy = TotalMomentum*sin(theta)*sin(phi);
		double fPz = TotalMomentum*cos(theta);
		particles[n].SetP(fPx, fPy, fPz);
		double trasverseImpulse = sqrt(fPx * fPx + fPy * fPy);
		h5 -> Fill(trasverseImpulse);
		double energy = particles[n].GetParticleTotalEnergy();
		h6 -> Fill(energy);
		
		double x = gRandom->Uniform(1);
		
			if (x < 0.4) { particles[n].SetIndex(0);
						   h1 -> Fill(0); }
			else if (x < 0.8) { particles[n].SetIndex(1);
			       				h1 -> Fill(1); }
			else if (x < 0.85) { particles[n].SetIndex(2);
						     	h1 -> Fill(2); }
			else if (x < 0.9) { particles[n].SetIndex(3);
						     	h1 -> Fill(3); }
			else if (x < 0.945) { particles[n].SetIndex(4);
							 	  h1 -> Fill(4); }
			else if (x < 0.99) { particles[n].SetIndex(5);
							 	 h1 -> Fill(5); }
			else {
			particles[n].SetIndex(6);
			Particle pion;
			Particle Kaon;
			if (x < 0.995) {
				pion.SetIndex(1);
				Kaon.SetIndex(3);
			}
			else {
				pion.SetIndex(0);
				Kaon.SetIndex(2);
			}
			
			particles[n].Decay2Body(pion , Kaon);
			particles[generatedParticles+decayedParticles] = pion;
			particles[generatedParticles+decayedParticles+1] = Kaon;
			decayedParticles += 2;		
	        } 
		
		
	  for (int j = 0; j < generatedParticles + decayedParticles - 1; j++) {

	  	for (int k = j + 1; k < generatedParticles + decayedParticles; k++) {
	  	  int jIndex = particles[j].GetIndex();
	  	  int kIndex = particles[k].GetIndex();
	  		
	  	  // filling h7
	  		h7 -> Fill(particles[j].InvMass(particles[k]));
	  	  // filling h8	
	  		if ((jIndex == 0 || jIndex == 2 || jIndex ==  4) 
	  			&& (kIndex == 1 || kIndex == 3 || kIndex == 5)) {
	  		h8 -> Fill(particles[j].InvMass(particles[k]));
	  		} 
	  	  // filling h9
	  	    if (((jIndex == 0 || jIndex == 2 || jIndex ==  4) 
	  	        && (kIndex == 0 || kIndex == 2 || kIndex ==  4))
	  	       || ((jIndex == 1 || jIndex == 3 || jIndex ==  5) 
	  	        && (kIndex == 1 || kIndex == 3 || kIndex ==  5)))
	  	       
	  	         {
	  	    h9 -> Fill(particles[j].InvMass(particles[k]));
	  	    }
	  	  // filling h10
	  	    if ((jIndex == 0 && kIndex == 3) || 
	  	    	(jIndex == 1 && kIndex == 2)) {
	  	    h10 -> Fill(particles[j].InvMass(particles[k]));
	  	    }
	  	  // filling h11
	  	    if ((jIndex == 0 && kIndex == 2 ) || 
	  	    	(jIndex == 1 && kIndex == 3)) {
	  	    h11 -> Fill(particles[j].InvMass(particles[k]));
	  	    }
	     }
	  }
   }
		if(i % 1000 == 0) {std::cout << i/1000 << "%" << '\n';}
}

	TCanvas* c1 = new TCanvas();
	TCanvas* c2 = new TCanvas();
	c1 -> Divide(3,2);
	c2 -> Divide(3,2);
	
	c1 -> cd(1);
	h1 -> Draw();
	c1 -> cd(2);
	h2 -> Draw();
	c1 -> cd(3);
	h3 -> Draw();
	c1 -> cd(4);
	h4 -> Draw();
	c1 -> cd(5);
	h5 -> Draw();
	c1 -> cd(6);
	h6 -> Draw();
	c2 -> cd(1);
	h7 -> Draw();
	c2 -> cd(2);
	h8 -> Draw();
	c2 -> cd(3);
	h9 -> Draw();
	c2 -> cd(4);
	h10 -> Draw();
	c2 -> cd(5);
	h11 -> Draw();
	c2 -> cd(6);
	h12 -> Draw();
	
	c1 -> Print("c1.pdf ", "pdf ");
	c2 -> Print("c2.pdf ", "pdf ");
	
	file->Write();
    file->Close();
 
}

	

 
  
  
  
  
  
  
  
  
  
  
  
       
