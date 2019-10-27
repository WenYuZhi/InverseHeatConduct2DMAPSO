#include "Temperature.h"
#include "Mesh.h"
#include <iostream>
using namespace std;

Temperature::Temperature(Mesh* mesh, Steel* steel, Continuous_Caster* caster, float* h_mold, const float vcast, const float cast_temperature)
{
	nx = mesh->nx;
	ny = mesh->ny;
	tnpts = mesh->tnpts;
	tf = mesh->tf;
	lx = mesh->lx;
	ly = mesh->ly;
	dx = mesh->dx;
	dy = mesh->dy;
	tao = mesh->tao;
	T_New = new float[nx * ny];
	T_Last = new float[nx * ny];
	T_Surface = new float[tnpts];
	this->h_mold = h_mold;
	this->vcast = vcast;
	this->cast_temperature = cast_temperature;
	tstep = 0;
	disout = 0;
	this->steel = steel;
	this->caster = caster;
	this->n_dim = caster->coolsection;
	objvalues = 0.0f;
}

Temperature::~Temperature()
{
	delete [] T_New;
	delete [] T_Last;
	delete [] T_Surface;
}

void Temperature::setmeasure(const float* measuredpoistion, const float* measuredtemperature, int measurednum)
{
	this->measuredpoistion = new float[measurednum];
	this->measuredtemperature = new float[measurednum];
	for (int i = 0; i < measurednum; i++){
		this->measuredpoistion[i] = measuredpoistion[i];
		this->measuredtemperature[i] = measuredtemperature[i];
	}
	this->measurednum = measurednum;
}

void Temperature::differencecalculation()
{
	float a, Tw = 30.0, T_Up, T_Down, T_Right, T_Left, T_Middle;
	map<string, float>* steel_parameters;
	if (disout == 0)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				steel_parameters = steel->compute_physicial_parameters(T_Last[i * ny + j]);
				a = steel_parameters->at("k") * tao / (steel_parameters->at("pho") * steel_parameters->at("ce"));
				if (i == 0 && j != 0 && j != ny - 1)  //1
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j] - 2 * dx * h * (T_Last[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j != 0 && j != ny - 1)//2
				{
					T_Up = T_Last[(i - 1)*ny + j] - 2 * dx*h*(T_Last[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == 0 && i != 0 && i != nx - 1)//3
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1] - 2 * dy*h*(T_Last[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == ny - 1 && i != 0 && i != nx - 1)//4
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1] - 2 * dy*h*(T_Last[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == 0)//5
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == ny - 1)//6
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == 0)//7
				{
					T_Up = T_Last[(i - 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == ny - 1)//8
				{
					T_Up = T_Last[(i - 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else//9
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
			}
		}
		T_Surface[tstep] = T_New[int((ny - 1) / 2)];
	}
	else
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				steel_parameters = steel->compute_physicial_parameters(T_New[i * ny + j]);
				a = steel_parameters->at("k") * tao / (steel_parameters->at("pho") * steel_parameters->at("ce"));
				if (i == 0 && j != 0 && j != ny - 1)  //1
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j] - 2 * dx*h*(T_New[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j != 0 && j != ny - 1)//2
				{
					T_Up = T_New[(i - 1)*ny + j] - 2 * dx*h*(T_New[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == 0 && i != 0 && i != nx - 1)//3
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1] - 2 * dy*h*(T_New[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == ny - 1 && i != 0 && i != nx - 1)//4
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j - 1] - 2 * dy*h*(T_New[i*ny + j] - Tw) / steel_parameters->at("k");
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == 0)//5
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == ny - 1)//6
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j];
					T_Right = T_New[i*ny + j - 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == 0)//7
				{
					T_Up = T_New[(i - 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == ny - 1)//8
				{
					T_Up = T_New[(i - 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j - 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else//9
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
			}
		}
		T_Surface[tstep] = T_Last[int((ny - 1) / 2)];
	}
	disout = !disout;
	tstep++;
}

void Temperature::boundarycondition(float *h_init)
{
	float zposition = fabs(tstep * tao * vcast);
	for (int i = 0; i < caster->section; i++)
	{
		if (zposition >= *(caster->ccml + i) && zposition <= *(caster->ccml + i + 1) && i < caster->moldsection)
			h = *(h_mold + i);
		if (zposition >= *(caster->ccml + i) && zposition <= *(caster->ccml + i + 1) && i >= caster->moldsection)
			h = *(h_init + i - caster->moldsection);
	}	
}

void Temperature::initcondition()
{
	tstep = 0;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			T_Last[ny * i + j] = cast_temperature;
			T_New[ny * i + j] = cast_temperature;
		}
	disout = 0;
}

void Temperature::solve(float *h)
{
	this->initcondition();
	while (tstep < tnpts) 
	{
		this->boundarycondition(h);
		this->differencecalculation();
	}
}

float Temperature::computeobjvalues()
{
	calculation_temperature = new float[measurednum];
	objvalues = 0.0;
	for (int i = 0; i < measurednum; i++)
	{
		for (int j = 0; j < tnpts; j++)
		{
			if ((fabs(j * vcast * tao - *(measuredpoistion + i))) <= fabs(tao * vcast / 2.0f))
				calculation_temperature[i] = T_Surface[j];
		}
		objvalues += (calculation_temperature[i] - measuredtemperature[i]) * (calculation_temperature[i] - measuredtemperature[i]);
	}
	
	delete [] calculation_temperature;
	objvalues = objvalues / measurednum;
	return objvalues;
}

void Temperature::print()
{
	cout << "the length of steel billets = " << lx << ", ";
	cout << "the width of steel billets = " << ly << ", ";
	cout << "casting speed = " << vcast << ", " << endl;
	cout << "dx = " << dx << ", ";
	cout << "dy = " << dy << ", ";
	cout << "time step = " << tao << ", " << endl;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			cout << T_New[ny * i + j] << ", ";
	cout << endl;
}

