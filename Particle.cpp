#include "Particle.hpp"
#include "ResonanceType.hpp"
#include "ParticleType.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
	
	int Particle::fNParticleType = 0;
	ParticleType* Particle::fParticleType[fMaxNumParticleType];
	// Getters
	const char* Particle::GetParticleName() const {
		return fParticleType[fIParticle_]->ParticleType::GetName();
	}
	double Particle::GetParticleMass() const {
		return fParticleType[fIParticle_]->ParticleType::GetMass();
	}
	double Particle::GetPx() const { return Px_; }
    double Particle::GetPy() const { return Py_; }
    double Particle::GetPz() const { return Pz_; } 
	int Particle::GetIndex() { return fIParticle_; }
	// Setters
	void Particle::SetPx(double fPx) { Px_ = fPx; }
	void Particle::SetPy(double fPy) { Py_ = fPy; }
	void Particle::SetPz(double fPz) { Pz_ = fPz; }
	void Particle::SetIndex(int fIndex) { 
		if (fIndex >= fNParticleType) { std::cout << "The particle doesn't exist " << '\n';}
		fIParticle_ = fIndex;
		}
	void Particle::SetName(char* particleName) {
		int index = FindParticle(particleName);
		fIParticle_ = index ;
	}
	// Private Member
	int Particle::FindParticle(char* Name) {
			for(int i = 0 ; i < fNParticleType ; ++i) { // se c'è corrispondenza fra il nome in input e quello dell'oggetto puntato ritorna l'indice
				  if (fParticleType[i]->ParticleType::GetName() == Name) {
				  		return i;
				  }                                         
			} 
			
			return fNParticleType;
	    }
	// Public Member
	void Particle::AddParticleType(char* lName, double fMass, int fCharge, double fWidth) {
		if (FindParticle(lName) < fNParticleType) {
			std::cout << "Particle " << lName << "already exists " << '\n';
		} else if (FindParticle(lName) == fMaxNumParticleType) {
				std::cout << "Array is already full " << '\n';
			} else if (fWidth == 0) {
			  		fParticleType[fNParticleType] = new ParticleType(lName, fMass, fCharge);
			  		++fNParticleType;	
				} else {
					fParticleType[fNParticleType] = new ResonanceType(lName, fMass, fCharge, fWidth);
					++fNParticleType;			
			} 
	} 
	
	void Particle::Printf() {
			std::cout << "Index: " << fIParticle_ << '\n';
			std::cout << "Name: " << Particle::GetParticleName() << '\n';
			std::cout << "Px: " << Px_ << '\n';
			std::cout << "Py: " << Py_ << '\n';
			std::cout << "Pz: " << Pz_ << '\n';
		}
		
	void Particle::PrintArray() {
			for(int i = 0; i != fNParticleType; i++) {
				fParticleType[i]->Print();
			}
	}
	
	double Particle::GetTotalMomentum() const {
			return pow(Px_ , 2) + pow(Py_ , 2) + pow(Pz_ , 2) ;
	}
	
	double Particle::GetParticleTotalEnergy() const {
			return sqrt((pow(Particle::GetParticleMass() , 2))+GetTotalMomentum() );	
	}
	
	double Particle::InvMass(Particle const& p) const {
			return sqrt(pow((p.GetParticleTotalEnergy()+GetParticleTotalEnergy()) , 2)+
						pow((p.GetTotalMomentum()+GetTotalMomentum()) , 2));
	}
	
	void Particle::SetP(double Px, double Py, double Pz) {
			Px_ = Px;
			Py_ = Py;
			Pz_ = Pz;
	}
	
	int Particle::Decay2Body(Particle &dau1 , Particle &dau2) {
  if(GetParticleMass() == 0.0){
    printf("Decayment cannot be performed if mass is zero\n");
    return 1;
  }
  
  double massMot = GetParticleMass();
  double massDau1 = dau1.GetParticleMass();
  double massDau2 = dau2.GetParticleMass();

  if(fIParticle_ > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleType[fIParticle_]->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*
  				(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(Px_*Px_ + Py_*Py_ + Pz_*Pz_ + massMot*massMot);

  double bx = Px_/energy;
  double by = Py_/energy;
  double bz = Pz_/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

  double energy = GetParticleTotalEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*Px_ + by*Py_ + bz*Pz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  Px_ += gamma2*bp*bx + gamma*bx*energy;
  Py_ += gamma2*bp*by + gamma*by*energy;
  Pz_ += gamma2*bp*bz + gamma*bz*energy;
}


	
	
	
	
	
	

		
		
		
		
		
		
		
		
		
		
	
	
