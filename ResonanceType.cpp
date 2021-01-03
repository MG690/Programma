#include "ResonanceType.hpp"
#include "ParticleType.hpp"
#include <iostream>

ResonanceType::ResonanceType(char* Name, double Mass, int Charge, double Width): ParticleType(Name, Mass, Charge), Width_{Width} {};

double ResonanceType::GetWidth() const {
	return Width_; 
	}

void ResonanceType::Print() const {
	
	ParticleType::Print();
	std::cout << "Width: " << Width_ << '\n';
	
}



