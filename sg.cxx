#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* f1, cmplx* f0,
          const double alpha, const double lambda,
          const double dx, const double dt, const int N,
		  const double xmin, const double kappa);
//-----------------------------------
int main(){

//+++++++++++++++USERINTERFACE+++++++++++++++++++//
	const int Nx = 300;
	const double xmin = -40.0;
  const double xmax = 40.0;
	const double Tend = 11.0*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt =  Tend/400;//Crank Nicolson is unconditionally stable so choice is not that important
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	
	const double lambda = 10.0;
  const double omega = 0.2;
//+++++++++++++++++++++++++++++++++++++++++++++++//

	//auxillary values
	const double kappa = pow(omega, 2.0);
	const double alpha = pow(kappa,0.25);
  

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx]; // contains solution for timeinterval n
	cmplx* psi1 = new cmplx[Nx]; //contains computed solution for timeinterval n+1

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);


	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
			
			//iterate one step in time
			step(psi1, psi0, alpha, lambda, dx, dt, Nx, xmin, kappa);
			t+=dt;
			
			// swap psi0 <-> psi1
			cmplx* temp;
			temp = psi0;
			psi0 = psi1;
			psi1 = temp;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
	cout << "t = " << t << endl;
  
	delete[] psi0;
	delete[] psi1;

	return 0;
}
//-----------------------------------
void step(cmplx* f1, cmplx* f0,
          const double alpha, const double lambda,
          const double dx, const double dt, const int N, const double xmin, const double kappa)
{
	cmplx a = cmplx(0.0, -dt/(4*dx*dx) );
	cmplx* d = new cmplx[N];
	
	for (int i=0; i<N; i++){
		double x = xmin + i*dx;
		d[i] = cmplx(1.0, dt/(2.0*dx*dx)+dt*kappa*x*x/4.0);
	}
	
	//Forward Substitution: 
	//LHS of eq. sys: A becomes A= diag(dvec', 0) + diag(a, 1)
	//RHS of eq. sys: A* becomes A* = diag(dvex'*, 0) + diag(a*,1)
	//after forloop computes dvec' and write it into d[i]
	for (int i= 1; i<N; i++){
		d[i] = d[i] - a*a/d[i-1];
	}
	
	//Backward Substitution:
	f1[N-1]= conj(d[N-1])*f0[N-1]/d[N-1];
	for (int i= N-2; i>=0; i--){
		f1[i] = ( conj(a)*f0[i+1]+ conj(d[i])*f0[i]- a*f1[i+1] )/d[i];
	}
	
	delete[] d;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
