#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include<iostream>

void ParticleType::Print() const {
	std::cout << "Name: " << Name_ << '\n';
	std::cout << "Mass: " << Mass_ << '\n';
	std::cout << "Charge: " << Charge_ << '\n';
}

const char* ParticleType::GetName() const { 
	return Name_; 
	}
	
double ParticleType::GetMass() const {
	return Mass_; 
	}
	
int ParticleType::GetCharge() const {
	return Charge_; 
	}

double ParticleType::GetWidth() const {
	return GetWidth();
}
