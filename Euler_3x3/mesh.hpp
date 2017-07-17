#ifndef MESH_HPP
#define MESH_HPP

#include<vector>
#include<cmath>

using namespace std;

class mesh
{
private:
	double xa;
	double xb;
	double Dx;
	int N;
public:
	mesh();
	mesh(double,double,double);
	mesh(double,double,int);
	inline double get_xa() {return xa;};
	inline double get_xb() {return xb;};
	inline double get_Dx() {return Dx;};
	inline double get_xi(int i) {return xa+Dx*(i+0.5);}
	inline int get_N() {return N;}
	inline void set_xa(double A) {xa = A;};
	inline void set_xb(double B) {xb = B;};
	inline void set_Dx(double delta) {Dx = delta; N = int(round((xb-xa)/Dx)); };
	inline void set_N(int n) {N=n; Dx=(xb-xa)/N;}
};

#endif