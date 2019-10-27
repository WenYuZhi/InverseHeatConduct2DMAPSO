#ifndef TEMPERATURE_H
#define TEMPERATURE_H
#include "Continuous_Caster.h"
#include "Mesh.h"
#include "Steel.h"
class Temperature
{
    private:
	    float vcast, cast_temperature;
	    float h;
	    float* T_New;
	    float* T_Last;
	    float* T_Surface;
		float* h_mold;
	    bool disout;
	    int nx, ny, tnpts;
	    float dx, dy, tao, tf, lx, ly;
		Steel* steel;
		Continuous_Caster* caster;
    public:
	    int tstep, n_dim, measurednum;
		float objvalues;
		float* calculation_temperature, *measuredpoistion, *measuredtemperature;
	    Temperature(Mesh*, Steel*, Continuous_Caster*, float*, const float, const float);
	    ~Temperature();
		void setmeasure(const float*, const float*, int);
	    void differencecalculation();
	    void boundarycondition(float*);
	    void initcondition();
		void solve(float *);
	    void print();
	    float computeobjvalues();
};
#endif