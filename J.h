/*
   *** HEADER FOR J's.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   ELIMINATING kp!! Final integral over Qi
   Full-RPA for the conductivity.
   T = 0K
   Different doping combinations for the graphene sheets.
   Calculates J(w) according to notes 15/feb/2019.
   I1(kp,w) is tabulated and interpolated (NR3)
   Combinations of Ef are: 
      Ef2=0.3eV    Ef1=0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0eV
      Ef2=0.5eV    Ef1=1.0eV   
*/ 
# pragma once

#include "sve/sve.h"
#include "sve/interp_linear.h"
#include "sve/adapt.h"


// w2.dat and IntPhiWo
int    nk = 1100;
int    nw = 900;

// w3.dat and IntW
//int nk = 510;
//int nw = 500;

double Ef1,Ef2,Vb,w,g0,eta0;
double Qi;
double En1,En2;
double En1m,En2m;
double En1p,En2p;
double kp01m,kp01p,kp02m,kp02p; 
double kp03m,kp03p,kp04m,kp04p;
double Qf01m,Qf01p,Qf02m,Qf02p; 
double Qf03m,Qf03p,Qf04m,Qf04p;
double A1m,A2m,A3m,A4m;
double A1p,A2p,A3p,A4p;
double B1m,B2m,B3m,B4m;
double B1p,B2p,B3p,B4p;
double infites=0.0;

int i,j,l;
double sum;
double var1 = 0.001;
double var2 = 0.999*pi;

VecDoub kpt(nk); // tabulated kp vector for interpolating
VecDoub I1t(nk); // tabulated I1 vector for interpolating
Linear_interp I1s(kpt,I1t); // interpolation constructor
Adapt apt(1.0e-6);

// --- Fermi-Dirac distribution at T=0K (Step function)
double FD(double En){
	if(En>0)       return 1.0;
	else if(En==0) return 0.5;
	else           return 0.0;
}

// --- kp0: poles of the deltas
// --- (1)
double kp001m(double Qi,double vph,double g0){ 
	double izq   = Qi*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   double der;
   if(inqrt<0) der = 0.0;
	else        der = sqrt(inqrt);
	
   double res = (izq - der);
   if(res<0) return 0.0;
	else      return res;
}

double kp001p(double Qi,double vph,double g0){ 
	double izq   = Qi*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   double der;
   if(inqrt<0) der = 0.0;
	else        der = sqrt(inqrt);
	
   double res = izq + der;
   if(res<0) return 0.0;
	else      return res;
}

// --- (2)
double kp002m(double Qi,double vph,double g0){ 
	double izq   = Qi*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   double der;
   if(inqrt<0) der = 0.0;
	else        der = sqrt(inqrt);
	
   double res = izq - der;
   if(res<0) return 0.0;
	else      return res;
}

double kp002p(double Qi,double vph,double g0){ 
	double izq   = Qi*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   double der;
   if(inqrt<0) der = 0.0;
	else        der = sqrt(inqrt);
	
   double res = izq + der;
   if(res<0) return 0.0;
	else      return res;
}

// --- (3)
double kp003m(double Qi,double vph,double g0){ 
	return kp001m(Qi,vph,g0);
}

double kp003p(double Qi,double vph,double g0){ 
	return kp001p(Qi,vph,g0);
}

// --- (4)
double kp004m(double Qi,double vph,double g0){ 	
	return kp002m(Qi,vph,g0);
}

double kp004p(double Qi,double vph,double g0){ 	
	return kp002p(Qi,vph,g0);
}



// --- A = [1+Qbi*Qbf/Qi/Qf]/|F'(kp0)| 
// ---     delta(F(Qi))=delta(Qi-kp0)/|F'(kp0)|
// --- (1)
double A01m(double Qi,double vph,double g0){ 
	double k01   = kp001m(Qi,vph,g0);
   double num   = 2.0*Qi - g0 - k01*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
  
}

double A01p(double Qi,double vph,double g0){ 
	double k01   = kp001p(Qi,vph,g0);
   double num   = 2.0*Qi - g0 - k01*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}


// --- (2)
double A02m(double Qi,double vph,double g0){ 
	double k02   = kp002m(Qi,vph,g0);
   double num   = - g0 - k02*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

double A02p(double Qi,double vph,double g0){ 
	double k02   = kp002p(Qi,vph,g0);
   double num   = - g0 - k02*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

// --- (3)
double A03m(double Qi,double vph,double g0){ 
	double k03   = kp003m(Qi,vph,g0);
   double num   = g0 - k03*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

double A03p(double Qi,double vph,double g0){ 
	double k03   = kp003p(Qi,vph,g0);
   double num   = g0 - k03*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

// --- (4)
double A04m(double Qi,double vph,double g0){ 
	double k04   = kp004m(Qi,vph,g0);
   double num   = g0 + 2.0*Qi - k04*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
   else         return num/sqrt(inqrt);
}

double A04p(double Qi,double vph,double g0){ 
	double k04   = kp004p(Qi,vph,g0);
   double num   = g0 + 2.0*Qi - k04*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
   else         return num/sqrt(inqrt);
}


// --- B = [1-Qbi*Qbf/Qi/Qf]/|F'(kp0)| 
// ---     delta(F(Qi))=delta(Qi-kp0)/|F'(kp0)|
// --- (1)
double B01m(double Qi,double vph,double g0){ 
	double k01   = kp001m(Qi,vph,g0);
   double num   = - g0 + k01*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
  
}

double B01p(double Qi,double vph,double g0){ 
	double k01   = kp001p(Qi,vph,g0);
   double num   = - g0 + k01*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}


// --- (2)
double B02m(double Qi,double vph,double g0){ 
	double k02   = kp002m(Qi,vph,g0);
   double num   = - g0 - 2.0*Qi + k02*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

double B02p(double Qi,double vph,double g0){ 
	double k02   = kp002p(Qi,vph,g0);
   double num   = - g0 - 2.0*Qi + k02*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

// --- (3)
double B03m(double Qi,double vph,double g0){ 
	double k03   = kp003m(Qi,vph,g0);
   double num   = g0 - 2.0*Qi + k03*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

double B03p(double Qi,double vph,double g0){ 
	double k03   = kp003p(Qi,vph,g0);
   double num   = g0 - 2.0*Qi + k03*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 - 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
	else         return num/sqrt(inqrt);
}

// --- (4)
double B04m(double Qi,double vph,double g0){ 
	double k04   = kp004m(Qi,vph,g0);
   double num   = g0 + k04*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
   else         return num/sqrt(inqrt);
}

double B04p(double Qi,double vph,double g0){ 
	double k04   = kp004p(Qi,vph,g0);
   double num   = g0 + k04*cos(vph);
   double inqrt = Qi*Qi*cos(vph)*cos(vph) + g0*g0 + 2.0*Qi*g0;
   if(inqrt<=0) return 0.0;
   else         return num/sqrt(inqrt);
}



// --- Qf0: Qf's evaluated at kp0
double Qf00(double Qi,double vph,double kp0){
   double inqrt = Qi*Qi + kp0*kp0 - 2.0*Qi*kp0*cos(vph);

   if(inqrt<=0) return 0.0;
	else         return sqrt(inqrt);
} 
