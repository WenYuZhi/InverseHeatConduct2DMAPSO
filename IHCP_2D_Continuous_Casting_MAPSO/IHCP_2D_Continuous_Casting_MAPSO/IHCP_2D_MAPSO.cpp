#include <iostream>
#include "Continuous_Caster.h"
#include "Mesh.h"
#include "Steel.h"
#include "Temperature.h"
#include "Classic_PSO.h"
#include <string>
#include <map>
#include <stdlib.h> 
#include <Eigen/Dense>
#include <time.h> 
#include <vector>
#include "MA_PSO.h"
using namespace Eigen;
using namespace std;
int main()
{
	const int section = 12, coolsection = 8, moldsection = 4, measurednum = 8;
	const float ccml[section + 1] = { 0.0f,0.2f,0.4f,0.6f,0.8f,1.0925f,2.27f,4.29f,5.831f,9.6065f,13.6090f,19.87014f,28.599f };
	const float measuredpoistion[measurednum] = { 0.9463f, 1.6812f, 3.28f, 5.0605f, 7.7188f, 11.6077f, 16.7395f, 24.235f };
	const float measuredtemperature[measurednum] = {1010.0f, 960.0f, 932.0f, 905.0f, 883.0f, 862.0f, 852.0f, 815.0f };
	float h_init[coolsection] = { 873.16f,765.05f,524.32f,392.83f,328.94f,251.64f,146.16f,120.96f };
	float h_mold[moldsection] = { 1380.0f,1170.0f,980.0f,800.0f };
	const float lx = 0.25f, ly = 0.25f, vcast = 0.03666666f, cast_temperature = 1558.0f;
	const float tf = ccml[section] / vcast;
	const int nx = 15, ny = 15, tnpts = 3001;
	srand((unsigned int)time(0));

	map<string, float> steel_components = { {"C", 0.12f},{ "Mn", 0.35f },{ "Si", 0.2f },{ "S", 0.04f } ,{ "P", 0.03f } ,{ "Cr", 0.0f } ,{ "Ni", 0.0f } ,{ "Cu", 0.0f } };

	Continuous_Caster caster = Continuous_Caster(section, coolsection, moldsection, ccml);
	caster.print();
	Mesh mesh = Mesh(nx, ny, tnpts, lx, ly, tf);
	mesh.print();
	Steel steel = Steel(&steel_components);
	steel.print();

	Temperature SteelTemperature = Temperature(&mesh, &steel, &caster, h_mold, vcast, cast_temperature);
	SteelTemperature.solve(h_init);
	SteelTemperature.setmeasure(measuredpoistion, measuredtemperature, measurednum);
	float objvalues = SteelTemperature.computeobjvalues();
	cout << "objvalues = " << objvalues << endl;

	int n_agent = 10;
	float c1 = 1.5f, c2 = 1.5f, omga = 0.6f, vmax = 10.0f, vmin = -10.0f, pmax = 300.0f, pmin = 100.0f;

	Classic_PSO ClassicPSOalgorithm = Classic_PSO(n_agent, c1, c2, omga, vmax, vmin, pmax, pmin, &SteelTemperature);

	int width = 10;
	int sl_width = 3;
	MA_PSO MAPSOalgorithm = MA_PSO(width, sl_width, c1, c2, omga, vmax, vmin, pmax, pmin, &SteelTemperature);

	clock_t start;
	int max_iteration_times = 200;

	start = clock();
	for (int i = 0; i < max_iteration_times; i++)
	{
		MAPSOalgorithm.computeobjvalues();
		MAPSOalgorithm.getneighbor();
		MAPSOalgorithm.getbestneighbor();
		//MAPSOalgorithm.print_neighbor();
		MAPSOalgorithm.competition();

		MAPSOalgorithm.updategbest();
		MAPSOalgorithm.updatelbest();
		MAPSOalgorithm.updatev();
		MAPSOalgorithm.updateposition();

		MAPSOalgorithm.getselflearninggrid();
		MAPSOalgorithm.computeselflearningobjvalues();
		MAPSOalgorithm.getneighborselflearning();
		MAPSOalgorithm.getbestneighborselflearning();
		MAPSOalgorithm.competitionselflearning();

		MAPSOalgorithm.updategbest();
		MAPSOalgorithm.print();
	}

	cout << "running time = " << clock() - start << endl;
	max_iteration_times = 0;
	for (int i = 0; i < max_iteration_times; i++) 
	{
		ClassicPSOalgorithm.computeobjvalues();
		ClassicPSOalgorithm.updategbest();
		ClassicPSOalgorithm.updatelbest();
		ClassicPSOalgorithm.updatev();
		ClassicPSOalgorithm.updateposition();
	}
}