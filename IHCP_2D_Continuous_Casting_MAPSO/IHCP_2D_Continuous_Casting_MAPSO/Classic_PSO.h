#ifndef CLASSIC_PSO_H
#define CLASSIC_PSO_H
#include "Continuous_Caster.h"
#include "Mesh.h"
#include "Steel.h"
#include <Eigen/Dense>
#include "Temperature.h"
using namespace Eigen;
class Classic_PSO
{
    private:
	    float c1, c2, omga;
	    float vmax, vmin;
	    int n_agent;
	    
    public:
		float gbest;
		float* obj_values;
		int iteration_times;
		MatrixXf gbest_position;
		MatrixXf lbest;
		MatrixXf lbest_position;
		MatrixXf position;
		MatrixXf v;
		Temperature* opt_problem;
		Classic_PSO(int n_agent, float c1, float c2, float omga, float vmax, float vmin, float, float, Temperature*);
		void computeobjvalues();
		void updategbest();
		void updatelbest();
		void updatev();
		void updateposition();
		void print();
		~Classic_PSO();
};
#endif