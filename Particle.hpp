#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <iostream>

class Particle {
public:
	Particle(char* Name, double Px, double Py, double Pz): Px_{Px}, Py_{Py}, Pz_{Pz} {
	fIParticle_ = Particle::FindParticle(Name);
    if (fIParticle_ == fNParticleType) std::cout << "Particle " << Name << "is not in the array " << '\n';
	};
	Particle() = default;
	const char* GetParticleName() const;
	double GetParticleMass() const;
	double GetPx() const; 
	double GetPy() const;
	double GetPz() const;
	int GetIndex();
	void SetPx(double fPx);
	void SetPy(double fPy);
	void SetPz(double fPz);
	void SetIndex(int fIndex);// fornirne una in overload col nome della particella
	void SetName(char* particleName);
	static void AddParticleType(char* lName, double fMass, int fCharge, double fWidth = 0);
	void Printf();
	void PrintArray();
	double GetTotalMomentum() const;
	double GetParticleTotalEnergy() const;
	double InvMass(Particle const& p) const;
	void SetP(double Px,double Py,double Pz);
	int Decay2Body(Particle &dau1,Particle &dau2);
	void Boost(double bx, double by, double bz);
private:
	static const int fMaxNumParticleType = 10; // size dell'array
	static ParticleType* fParticleType[fMaxNumParticleType]; // array di puntatori
	static int fNParticleType; // numero di particelle nell'array
	int fIParticle_ = 0; // indice della particella nell'array
	char* Name_;
	double Px_ = 0;
	double Py_ = 0;
	double Pz_ = 0;
	static int FindParticle(char* fParticleName);
};
#endif
