/*
   *** HEADER FOR Jee.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   ELIMINATING kp!! Final integral over Qi
   Full-RPA for the conductivity.
   T = 0K
   Electron (e) doped graphene to electron (e) doped graphene.
   Calculates Jee(w) according to notes 15/feb/2019.
   I1(kp,w) is tabulated and interpolated (NR3)
   Combinations of Ef are: 
      Ef2=0.3eV    Ef1=0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0eV
      Ef2=0.5eV    Ef1=1.0eV
	Notation: 
		c: conduction
		v: valence
		Iif_x: integrand for integral over x that corresponds to the tunnel from
		         i (c or v) to f (c or v)
*/ 

#include "J.h"

//double err=1.0e-21;  int nmax=6000;   // integration parameters

// --- Integral 1: cond-val
double Icv_vph(double vph){ // Integrand for  Int d_vph
   double mm,pp;
   double  I102m,I102p;
   kp02m = kp002m(Qi,vph,g0);
   kp02p = kp002p(Qi,vph,g0);
   Qf02m = Qf00(Qi,vph,kp02m);
   Qf02p = Qf00(Qi,vph,kp02p);
	A2m   = A02m(Qi,vph,g0);
   A2p   = A02p(Qi,vph,g0); 
   En2m  = Ef2 - vf*Qf02m;
   En2p  = Ef2 - vf*Qf02p;

   I102m = I1s.interp(kp02m);
   if     (w<vf*(-2.0*Qi-kp02m-eta0)) mm=0.0;
   else if(w>vf*(-2.0*Qi+kp02m-eta0)) mm=0.0;
   else  mm=kp02m*I102m*A2m*FD(En2m);

   I102p = I1s.interp(kp02p);
   if     (w<vf*(-2.0*Qi-kp02p-eta0)) pp=0.0;
   else if(w>vf*(-2.0*Qi+kp02p-eta0)) pp=0.0;
   else  pp=kp02p*I102p*A2p*FD(En2p);

   return mm+pp;
}

double Icv_Qi(double Qi_){ // integrand for  Int d_Qi
   Qi = Qi_;

   sum = 0.0;
   for(j=0;j<(nev-1);j++){
      sum += 0.5*(Icv_vph(phi[j+1])+Icv_vph(phi[j]))*(phi[j+1]-phi[j]);
   }

   return Qi*sum;
   //return Qi*Qi*qtrap(var1,var2,Icv_vph); // much slower
}



// --- Integral 2: cond-cond
double Icc_vph(double vph){ // integrand for  Int d_vph
   double  mm,pp;
   double  I101m,I101p;
   kp01m = kp001m(Qi,vph,g0);
   kp01p = kp001p(Qi,vph,g0);
   Qf01m = Qf00(Qi,vph,kp01m);
   Qf01p = Qf00(Qi,vph,kp01p);
	A1m   = A01m(Qi,vph,g0); 
   A1p   = A01p(Qi,vph,g0); 
	En1   = Ef1 - vf*Qi;
   En2m  = vf*Qf01m - Ef2;
   En2p  = vf*Qf01p - Ef2;
   
   I101m = I1s.interp(kp01m);
   if      (w<vf*(-kp01m-eta0)) mm=0.0;
   else if (w>vf*( kp01m-eta0)) mm=0.0;
   else    mm=kp01m*I101m*A1m*FD(En1)*FD(En2m);
   
   I101p = I1s.interp(kp01p);
   if      (w<vf*(-kp01p-eta0)) pp=0.0;
   else if (w>vf*( kp01p-eta0)) pp=0.0;
   else    pp=kp01p*I101p*A1p*FD(En1)*FD(En2p);

   return  mm+pp;
}

double Icc_Qi(double Qi_){ // Integrand for  Int d_Qi
   Qi = Qi_;

   sum = 0.0;
   for(j=0;j<(nev-1);j++){
      sum += 0.5*(Icc_vph(phi[j+1])+Icc_vph(phi[j]))*(phi[j+1]-phi[j]);
   }

   return Qi*sum;
   //return Qi*qtrap(var1,var2,Icc_vph); // much slower
}
