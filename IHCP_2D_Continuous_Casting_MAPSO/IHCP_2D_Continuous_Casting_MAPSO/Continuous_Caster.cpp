#include "Continuous_Caster.h"
#include <iostream>
using namespace std;

Continuous_Caster::Continuous_Caster(const int section, const int coolsection, const int moldsection, const float* ccml)
{
	this->section = section;
	this->coolsection = coolsection;
	this->moldsection = moldsection;
	this->ccml = new float[section];
	for (int i = 0; i < section; i++)
		this->ccml[i] = ccml[i];
}

Continuous_Caster::~Continuous_Caster()
{
	delete [] ccml;
}

void Continuous_Caster::print()
{
	cout << "section = " << section << ", " << "coolsection = " << coolsection << " ";
	cout << "moldsection = " << moldsection << ", " << endl;
	for (int i = 0; i < section; i++)
		cout << ccml[i] << ", ";
	cout << endl;
}