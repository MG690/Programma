#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP
#include <iostream>

class ParticleType {
public:
	ParticleType(char* Name, double Mass, int Charge): Name_{Name}, Mass_{Mass}, Charge_{Charge} {};
	const char* GetName() const ;
	double GetMass() const ;
	int GetCharge() const ;
	virtual void Print() const ;
	virtual double GetWidth() const ;
		
private:
	char const* Name_;
	double const Mass_;
	int const Charge_;
};

#endif
