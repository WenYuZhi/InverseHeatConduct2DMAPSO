#pragma once
#ifndef MA_PSO_H
#define MA_PSO_H
#include <Eigen/Dense>
#include "Temperature.h"
#include <vector>
using namespace Eigen;
class MA_PSO
{
    private:
	    float c1, c2, omga;
	    float vmax, vmin;
	    int n_agent, width, width_sl;

    public:
		bool flag;
	    float gbest;
	    float* obj_values;
		float* obj_values_sl;
		int* best_neighbor_index, *best_neighbor_index_sl;
	    int iteration_times;
		vector<int>* neighbor_index;
		vector<int>* neighbor_index_sl;
	    MatrixXf gbest_position;
		MatrixXf lbest_position;
		MatrixXf lbest;
	    MatrixXf position;
		MatrixXf position_copy;
	    MatrixXf v;
		MatrixXf position_sl;
		MatrixXf position_sl_copy;
	    Temperature* opt_problem;
		MA_PSO(int width, int sl_width, float c1, float c2, float omga, float vmax, float vmin, float, float, Temperature*);
	    void computeobjvalues();
		void getneighbor();
		void getbestneighbor();
		void competition();
	    void updategbest();
		void updatelbest();
	    void updatev();
	    void updateposition();
		void getselflearninggrid();
		void computeselflearningobjvalues();
		void getneighborselflearning();
		void getbestneighborselflearning();
		void competitionselflearning();
	    void print();
		void print_neighbor();
	   ~MA_PSO();
};
#endif
