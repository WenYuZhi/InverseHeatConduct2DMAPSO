#include "steel.h"
#include <iostream>
#include <map>
#include <string>
using namespace std;

Steel::Steel(map<string, float> * components)
{
	this->components = components;
	temperature_l = 1536.6f - 88.0f * components->at("C") - 8.0f * components->at("Si") - 5.0f * components->at("Mn");
	temperature_s = 1527.0f - 187.0f * components->at("C") - 700.0f * components->at("S") - 500.0f * components->at("P") - 20.5f * components->at("Si") - 6.5f * components->at("Si");

}

void Steel::print()
{
	map<string, float>::iterator iter;
	for (iter = components->begin(); iter != components->end(); ++iter)
		cout << iter->first << ": " << iter->second << "  ";
	cout << "liquid temperature = " << temperature_l << ", " << "solid temperature = " << temperature_s << endl;
}

float Steel::compute_fs(float temperature_point)
{
	static float fs;
	if (temperature_point >= temperature_l)
		fs = 0.0f;
	if (temperature_point > temperature_s && temperature_point < temperature_l)
		fs = (temperature_l - temperature_point) / (temperature_l - temperature_s);
	if (temperature_point <= temperature_s)
		fs = 1.0f;
	return fs;
}

map<string, float> * Steel::compute_physicial_parameters(float temperature_point)
{
	float L = 268000.0f;
	float fs = this->compute_fs(temperature_point);
	static float ce, k, pho;
	
	if (temperature_point >= temperature_l)
	{
		ce = 540.0f;
		k = 50.0f;
		pho = 7250.0f;
	}
	if (temperature_point > temperature_s && temperature_point < temperature_l)
	{
		ce = 540 + L / (temperature_l - temperature_s);
		k = 25.0f * fs + 50.0f * (1.0f - fs);
		pho = 7250.0f;
	}
	if (temperature_point <= temperature_s)
	{
		ce = 540.0f;
		k = 28.0f;
		pho = 7250.0f;
	}

	static map<string, float> steel_parameters = { { "ce", ce },{ "pho", pho },{ "k", k } };
	return &steel_parameters;
}

	
	