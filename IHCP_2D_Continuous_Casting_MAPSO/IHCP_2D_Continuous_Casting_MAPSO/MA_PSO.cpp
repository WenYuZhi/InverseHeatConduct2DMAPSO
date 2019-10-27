#include "MA_PSO.h"
#include "algorithm"
#include <iostream>
#include <vector>
#include <map>
#include <string>
using namespace std;
MA_PSO::MA_PSO(int width, int width_sl, float c1, float c2, float omga, float vmax, float vmin, float pmax, float pmin, Temperature* opt_problem)
{
	flag = true;
	n_agent = width * width;
	this->c1 = c1;
	this->c2 = c2;
	this->omga = omga;
	this->vmax = vmax;
	this->vmin = vmin;
	this->width = width;
	this->width_sl = width_sl;
	this->opt_problem = opt_problem;
	iteration_times = 0;
	obj_values = new float[n_agent];
	obj_values_sl = new float[width_sl * width_sl];
	best_neighbor_index = new int[n_agent];
	best_neighbor_index_sl = new int[width_sl * width_sl];
	neighbor_index = new vector<int>[n_agent];
	neighbor_index_sl = new vector<int>[width_sl * width_sl];
	for (int i = 0; i < n_agent; i++)
		obj_values[i] = 0.0f;
	gbest = 1000000000000.0f;
	gbest_position = MatrixXf::Zero(1, this->opt_problem->n_dim);
	lbest = MatrixXf::Constant(n_agent, 1, 1000000000000000.0f);
	lbest_position = MatrixXf::Zero(n_agent, this->opt_problem->n_dim);
	position = (pmax - pmin) * MatrixXf::Random(n_agent, this->opt_problem->n_dim) + MatrixXf::Constant(n_agent, this->opt_problem->n_dim, pmin);
	position_copy = position;
	position_sl = MatrixXf::Zero(width_sl * width_sl, this->opt_problem->n_dim);
	position_sl_copy = MatrixXf::Zero(width_sl * width_sl, this->opt_problem->n_dim);
	v = (vmax - vmin) * MatrixXf::Random(n_agent, this->opt_problem->n_dim) + MatrixXf::Constant(n_agent, this->opt_problem->n_dim, vmin);
}

void MA_PSO::computeobjvalues()
{
	float *h_inv = new float[opt_problem->n_dim];
	for (int i = 0; i < n_agent; i++)
	{
		for (int j = 0; j < opt_problem->n_dim; j++)
			h_inv[j] = position(i, j);
		opt_problem->solve(h_inv);
		obj_values[i] = opt_problem->computeobjvalues();
	}
	delete[] h_inv;
}

void MA_PSO::getneighbor()
{
	map<string, int> neighbor_action = { { "N", -width },{ "E", 1 },{ "S", width },{ "W", -1 } };
	for (int index = 0; index < n_agent; index++)
	  for (map<string, int>::iterator iter = neighbor_action.begin(); iter != neighbor_action.end(); ++iter)
		if ((index < width && iter->first == "N") || (index > width * (width - 1) - 1 && iter->first == "S") || ((index + 1) % (width) == 0 && iter->first == "E") || (index % width == 0 && iter->first == "W"))
		{
			;
		}
		else
		    neighbor_index[index].push_back(index  + iter->second);
}

void MA_PSO::getbestneighbor()
{
	for (int i = 0; i < n_agent; i++)
	{
		best_neighbor_index[i] = neighbor_index[i][0];
		for (int j = 0; j < neighbor_index[i].size(); j++)
		{
			if (obj_values[best_neighbor_index[i]] > obj_values[neighbor_index[i][j]])
				best_neighbor_index[i] = neighbor_index[i][j];
		}
	}       
}

void MA_PSO::competition()
{
	MatrixXf r = MatrixXf::Random(1, 1);
	position_copy = position;
	for (int i = 0; i < n_agent; i++)
		if (obj_values[i] > obj_values[best_neighbor_index[i]])
			position_copy.row(i) = position.row(best_neighbor_index[i]) + r(0,0) * (position.row(i) - position.row(best_neighbor_index[i]));
	position = position_copy;
}


void MA_PSO::updategbest()
{
	if (flag == true) {
	    int index = min_element(obj_values, obj_values + n_agent) - obj_values;
	    if (gbest > obj_values[index])
	    {
		    gbest = obj_values[index];
		    gbest_position = position.row(index);
	    }
    }

	else {
		int index = min_element(obj_values_sl, obj_values_sl + width_sl * width_sl) - obj_values_sl;
		if (gbest > obj_values_sl[index])
		{
			gbest = obj_values_sl[index];
			gbest_position = position_sl.row(index);
		}
	}
	flag = !flag;
}

void MA_PSO::updatelbest()
{
	for (int i = 0; i < n_agent; i++)
	{
		if (lbest(i) > obj_values[i])
		{
			lbest(i) = obj_values[i];
			lbest_position.row(i) = position.row(i);
		}
	}
}

void MA_PSO::updatev()
{
	MatrixXf r = MatrixXf::Random(1, 2);
	for (int i = 0; i < n_agent; i++)
	{
		v.row(i) = omga * v.row(i) + c1 * r(0, 0) * (gbest_position - position.row(i)) ;
	}
	for (int i = 0; i < n_agent; i++)
		for (int j = 0; j < opt_problem->n_dim; j++)
		{
			if (v(i, j) < vmin)
				v(i, j) = vmin;
			if (v(i, j) > vmax)
				v(i, j) = vmax;
		}
}

void MA_PSO::updateposition()
{
	position += v;
	iteration_times += 1;
}

void MA_PSO::getselflearninggrid()
{
	float sr = 0.1f;
	for (int i = 0; i < width_sl * width_sl; i++)
		position_sl.row(i) = gbest_position * (2 * sr * MatrixXf::Random(1, opt_problem->n_dim) + (1 - sr) *  MatrixXf::Constant(1, opt_problem->n_dim, 1.0f));
}

void MA_PSO::computeselflearningobjvalues()
{
	float *h_inv = new float[opt_problem->n_dim];
	for (int i = 0; i < width_sl * width_sl; i++)
	{
		for (int j = 0; j < opt_problem->n_dim; j++)
			h_inv[j] = position_sl(i, j);
		opt_problem->solve(h_inv);
		obj_values_sl[i] = opt_problem->computeobjvalues();
	}
	delete[] h_inv;
}

void MA_PSO::getneighborselflearning()
{
	map<string, int> neighbor_action = { { "N", -width_sl },{ "E", 1 },{ "S", width_sl },{ "W", -1 } };
	for (int index = 0; index < width_sl * width_sl; index++)
		for (map<string, int>::iterator iter = neighbor_action.begin(); iter != neighbor_action.end(); ++iter)
		{
			if ((index < width_sl && iter->first == "N") || (index > width_sl * (width_sl - 1) - 1 && iter->first == "S") || ((index + 1) % (width_sl) == 0 && iter->first == "E") || (index % width_sl == 0 && iter->first == "W"))
			{
				;
			}
			else
				neighbor_index_sl[index].push_back(index + iter->second);
		}
}

void MA_PSO::getbestneighborselflearning()
{
	for (int i = 0; i < width_sl * width_sl; i++)
	{
		best_neighbor_index_sl[i] = neighbor_index_sl[i][0];
		for (int j = 0; j < neighbor_index_sl[i].size(); j++)
		{
			if (obj_values[best_neighbor_index_sl[i]] > obj_values[neighbor_index_sl[i][j]])
				best_neighbor_index_sl[i] = neighbor_index_sl[i][j];
		}
	}
}

void MA_PSO::competitionselflearning()
{
	MatrixXf r = MatrixXf::Random(1, 1);
	position_sl_copy = position_sl;
	for (int i = 0; i < width_sl * width_sl; i++)
		if (obj_values_sl[i] > obj_values_sl[best_neighbor_index_sl[i]])
			position_sl_copy.row(i) = position_sl.row(best_neighbor_index_sl[i]) + r(0, 0) * (position_sl.row(i) - position_sl.row(best_neighbor_index_sl[i]));
	position_sl = position_sl_copy;
}


void MA_PSO::print()
{
	cout << "iteration time: " << iteration_times << "  " << "gbest: " << gbest << endl;
}

void MA_PSO::print_neighbor()
{
	for (int i = 0; i < n_agent; i++)
	{
		for (int j = 0; j < neighbor_index[i].size(); j++)
			cout << neighbor_index[i][j] << ", ";
		cout << endl;
	}

	for (int i = 0; i < n_agent; i++)
		cout << obj_values[i] << ",";
	cout << endl;

	for (int i = 0; i < n_agent; i++)
		cout << best_neighbor_index[i] << ", ";
	cout << endl;

}

MA_PSO::~MA_PSO()
{
	delete[] obj_values;
}