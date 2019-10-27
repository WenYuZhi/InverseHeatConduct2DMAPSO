#ifndef CONTINUOUS_CASTER_H
#define CONTINUOUS_CASTER_H
class Continuous_Caster
{
    public:
	    int section, coolsection, moldsection;
	    float *ccml;
	    Continuous_Caster(const int, const int, const int, const float*);
	    ~Continuous_Caster();
	    void print();
};

#endif