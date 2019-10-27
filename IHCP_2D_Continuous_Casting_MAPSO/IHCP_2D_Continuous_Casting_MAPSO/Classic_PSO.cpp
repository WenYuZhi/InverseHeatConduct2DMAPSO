#include "Classic_PSO.h"
#include "algorithm"
#include <iostream>
using namespace std;
Classic_PSO::Classic_PSO(int n_agent, float c1, float c2, float omga, float vmax, float vmin, float pmax, float pmin, Temperature* opt_problem)
{
	this->n_agent = n_agent;
	this->c1 = c1;
	this->c2 = c2;
	this->omga = omga;
	this->vmax = vmax;
	this->vmin = vmin;
	this->opt_problem = opt_problem;
	iteration_times = 0;
	obj_values = new float[n_agent];
	for (int i = 0; i < n_agent; i++)
		obj_values[i] = 0.0f;
	gbest = 1000000000000.0f;
	gbest_position = MatrixXf::Zero(1, this->opt_problem->n_dim);
	lbest = MatrixXf::Constant(n_agent, 1, 1000000000000000.0f);
	lbest_position = MatrixXf::Zero(n_agent, this->opt_problem->n_dim);
	position = (pmax - pmin) * MatrixXf::Random(n_agent, this->opt_problem->n_dim) + MatrixXf::Constant(n_agent, this->opt_problem->n_dim, pmin);
	v = (vmax - vmin) * MatrixXf::Random(n_agent, this->opt_problem->n_dim) + MatrixXf::Constant(n_agent, this->opt_problem->n_dim, vmin);
}

void Classic_PSO::computeobjvalues()
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

void Classic_PSO::updategbest()
{
	int index = min_element(obj_values, obj_values + n_agent) - obj_values;
	if (gbest > obj_values[index])
	{
		gbest = obj_values[index];
		gbest_position = position.row(index);
	}
}

void Classic_PSO::updatelbest()
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

void Classic_PSO::updatev()
{
	MatrixXf r = MatrixXf::Random(1, 2);
	for (int i = 0; i < n_agent; i++)
	{
		v.row(i) = omga * v.row(i) + c1 * r(0, 0) * (gbest_position - position.row(i)) + c2 * r(0, 1) * (lbest_position.row(i) - position.row(i));
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

void Classic_PSO::updateposition()
{
	position += v;
	iteration_times += 1;
}

void Classic_PSO::print()
{
	cout << "iteration time: " << iteration_times << "  " << "gbest: " << gbest << ", ";
	cout << "min obj: " << *min_element(obj_values, obj_values + n_agent) << ", ";
	cout << "max obj: " << *max_element(obj_values, obj_values + n_agent) << endl;
}

Classic_PSO::~Classic_PSO()
{
	delete[] obj_values;
}