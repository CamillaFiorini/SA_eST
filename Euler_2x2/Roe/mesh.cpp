#include"mesh.hpp"

mesh::mesh () : xa(0), xb(1), Dx(1e-1), N(10) {};

mesh::mesh (double a,double b, double delta) : xa(a), xb(b), Dx(delta), N(int((xb-xa)/Dx)) {};
mesh::mesh(double a, double b, int n) : xa(a), xb(b), Dx((b-a)/n), N(n) {};
