#include "Mesh.h"
#include <iostream>
using namespace std;
Mesh::Mesh(const int nx, const int ny, const int tnpts, const float lx, const float ly, const float tf)
{
	this->nx = nx;
	this->ny = ny;
	this->tnpts = tnpts;
	this->lx = lx;
	this->ly = ly;
	this->tf = tf;
	this->dx = this->lx / (this->nx - 1);
	this->dy = this->ly / (this->ny - 1);
	this->tao = this->tf / (this->tnpts - 1);
}

void Mesh::print()
{
	cout << "nx = " << nx << ", " << "ny = " << ny << ", " << "tnpts = " << tnpts << endl;
	cout << "lx = " << lx << ", " << "ly = " << ly << ", " << "tf = " << tf << endl;
	cout << "dx = " << dx << ", " << "dy = " << dy << ", " << "tao = " << tao << endl;
}

Mesh::~Mesh()
{

}