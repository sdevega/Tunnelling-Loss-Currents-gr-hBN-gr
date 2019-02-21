/*
   *** HEADER FOR Jhh.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   ELIMINATING kp!! Final integral over Qi
   Full-RPA for the conductivity.
   T = 0K
   Hole (h) doped graphene to hole (h) doped graphene.
   Calculates Jhh(w) according to notes 15/feb/2019.
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


// --- Integral 1: val-cond
double Ivc_vph(double vph){ // Integrand for  Int d_vph
   double  mm,pp;
   double  I102m,I102p;
   kp02m = kp002m(Qi,vph,g0);
   kp02p = kp002p(Qi,vph,g0);
	A2m   = A02m(Qi,vph,g0);
   A2p   = A02p(Qi,vph,g0); 
   En1   = vf*Qi - Ef1;
   
   I102m = I1s.interp(kp02m);
   if     (w<vf*(-2.0*Qi-kp02m-eta0)) mm=0.0;
   else if(w>vf*(-2.0*Qi+kp02m-eta0)) mm=0.0;
   else  mm=kp02m*I102m*A2m*FD(En1);

   I102p = I1s.interp(kp02p);
   if     (w<vf*(-2.0*Qi-kp02p-eta0)) pp=0.0;
   else if(w>vf*(-2.0*Qi+kp02p-eta0)) pp=0.0;
   else  pp=kp02p*I102p*A2p*FD(En1);

   return mm+pp;
}

double Ivc_Qi(double Qi_){ // integrand for  Int d_Qi
   Qi = Qi_;
   return Qi*apt.integrate(Ivc_vph,var1,var2);
}



// --- Integral 2: val-val
double Ivv_vph(double vph){ // integrand for  Int d_vph
   double  mm,pp;
   double  I104m,I104p;
   kp04m = kp004m(Qi,vph,g0);
   kp04p = kp004p(Qi,vph,g0);
   Qf04m = Qf00(Qi,vph,kp04m);
   Qf04p = Qf00(Qi,vph,kp04p);
	A4m   = A04m(Qi,vph,g0); 
   A4p   = A04p(Qi,vph,g0); 
	En1   = vf*Qi - Ef1;
   En2m  = Ef2 - vf*Qf04m;
   En2p  = Ef2 - vf*Qf04p;
   
   I104m = I1s.interp(kp04m);
   if      (w<vf*(-kp04m-eta0)) mm=0.0;
   else if (w>vf*( kp04m-eta0)) mm=0.0;
   else    mm=kp04m*I104m*A4m*FD(En1)*FD(En2m);
   
   I104p = I1s.interp(kp04p);
   if      (w<vf*(-kp04p-eta0)) pp=0.0;
   else if (w>vf*( kp04p-eta0)) pp=0.0;
   else    pp=kp04p*I104p*A4p*FD(En1)*FD(En2p);

   return  mm+pp;
}

double Ivv_Qi(double Qi_){ // Integrand for  Int d_Qi
   Qi = Qi_;
   return Qi*apt.integrate(Ivv_vph,var1,var2);
}
