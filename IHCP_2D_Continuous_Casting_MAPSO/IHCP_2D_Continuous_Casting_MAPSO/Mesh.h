#ifndef MESH_H
#define MESH_H

class Mesh
{
    public:
	    int nx, ny, tnpts;
	    float dx, dy, tao;
	    float lx, ly, tf;
	    Mesh(const int nx, const int ny, const int tnpts, const float lx, const float ly, const float tf);
		~Mesh();
		void print();
};
#endif
