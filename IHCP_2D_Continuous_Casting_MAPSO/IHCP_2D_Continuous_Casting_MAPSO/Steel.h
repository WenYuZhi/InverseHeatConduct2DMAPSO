#ifndef STEEL_H
#define STEEL_H
#include <string>
#include <map>
#include <stdlib.h> 
using namespace std;
class Steel
{
    public:
	    float pho;
	    float ce;
	    float lamda;
		float temperature_l, temperature_s;
		map<string, float> * components;
		Steel(map<string, float> *);
		float compute_fs(float);
		map<string, float> * compute_physicial_parameters(float);
		void print();
};
#endif