#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "ParticleType.hpp"
class ResonanceType: public ParticleType {
public:
	ResonanceType(char* Name, double Mass, int Charge, double Width);
	double GetWidth() const override;
	void Print() const override;
	
private:
	double const Width_;
};

#endif
