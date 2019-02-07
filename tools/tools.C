#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "lineDist.h" 
#include "hYrnpt.h"
#include <time.h>
#include <iomanip>
#include <omp.h>

using namespace std;
extern float ran2(long *idum);
// translation move is done preferably along axis
// weighted with 5*cosine of angle between rod and driection of displacement
// higher weights can speed up equilibration but also often lead to jamming
// maxDisplace must be < 1. 

void hYr::displaceParticle(Rod &p) {
  
  double displace[3], stretchfactor, pR[3];
  double psi0, psi1, save;
  int ch=0;
  
  displace[0] = ran2(iran) - 0.5;
  displace[1] = ran2(iran) - 0.5;
  displace[2] = ran2(iran) - 0.5;
  stretchfactor = (fabs(displace[0]*p.M[0] +  displace[1]*p.M[1]
		    + displace[2]*p.M[2])*ARInv*5.0+1)*maxDisplace;
  
if(pbc==2){  
ch =  int((p.R[2]+stretchfactor*displace[2])*box->halfxInv[2]);
if(ch != 0){

// phi=pi/2:
	pR[0] = i2phi*ch*p.R[1];
	pR[1] = -i2phi*ch*p.R[0];
	pR[2] = p.R[2];

	psi0 = acos(-p.M[1]/sqrt(p.M[1]*p.M[1]+p.M[2]*p.M[2]));
	save = p.M[0];
	p.M[0] = i2phi*ch*p.M[1];
	p.M[1] = -i2phi*ch*save;
	psi1 = acos(-p.M[1]/sqrt(p.M[1]*p.M[1]+p.M[2]*p.M[2]));
	p.P = p.P+psi1-psi0;



}
else{
	for (int i2 =0; i2<3; i2++){
		pR[i2] = p.R[i2];
		}
}	
}
else{
	for (int i2 =0; i2<3; i2++){
		pR[i2] = p.R[i2];
		}
}	

  for (int i2 =0; i2<3; i2++){
    p.R[i2] = -int((pR[i2]+stretchfactor*displace[i2])*box->halfxInv[i2])
      *box->halfx[i2] + fmod((pR[i2]+stretchfactor*displace[i2]),box->halfx[i2]);
    
  } 
}


bool hYr::checkOverlapRodAll(Rod &p) {

  bool overlap = false;
  Cell *curCell;
  Piece *curPiece, *ppiece;
  int t=0;
  bool test;
  Rod cpr;
  for (int i=0; i<N_r; i++){
      Tested[i]=0;
  }
  for (int i=0; i<p.npieces; i++){
      ppiece = p.pieces[i];
      curCell = ppiece->cell;
     
      // rods in cell
      curPiece = curCell->firstPiece;
     
      while (curPiece){

	  if (curPiece->rod != &p){
	      test = false;
	      for (int k=0; k<t; k++){
		  if (Tested[k] == curPiece->rod){
		      test=true;
		      break;
		  }
	      }
	      if (!test){
		  overlap = checkOverlapRods(*(curPiece->rod), p);
		  if (overlap) return true;
		  Tested[t] = curPiece->rod;
		  t++;
	      }
	  }
	  curPiece = curPiece->next;
      }

     // neighbouring cells
     for(int j=0; j<26; j++){
	 curCell = ppiece->cell->neighbours[j];
	
         curPiece = curCell->firstPiece;
	 while (curPiece){

	     if (curPiece->rod != &p){

		 test = false;
		 for (int k=0; k<t; k++){
		     if (Tested[k] == curPiece->rod){
			 test=true;
			 break;
		     }
		 }
		 if (!test){
		     overlap = checkOverlapRods(*(curPiece->rod), p);
		     if (overlap) return true;
		     Tested[t] = curPiece->rod;
		     t++;               
		 }
	     }
	     curPiece = curPiece->next;
         }
     }
  } 
  return false;
}

  
  

bool hYr::checkOverlapRods(Rod &p1, Rod &p2) {

    double d[3], dm[3], p1R[3], p1M[3], distance;
    double save;
    double s0, t0;    // Dummy -- not used
    bool overlap = false;
	int ch=0;
  

if (pbc==2){
  ch = int((p1.R[2] - p2.R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*p1.R[1];
	p1R[1] = -i2phi*ch*p1.R[0];
	p1R[2] = p1.R[2];

	save = p1.M[0];
	p1M[0] = i2phi*ch*p1.M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = p1.M[2];
	

}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = p1.R[i2];
		p1M[i2] = p1.M[i2];
	}
}	
}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = p1.R[i2];
		p1M[i2] = p1.M[i2];
	}
}	

    for (int j=0; j<3; j++){
	d[j] = p1R[j] - p2.R[j];
	if (fabs(d[j]) > box->halfx[j]) 
	    d[j] -= copysign(box->x[j], d[j]);
    } 
    
    if ((d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) < AR12){
	
	for (int k=0; k<3; k++){
	    dm[k] = 0.5*(p1M[k] - p2.M[k]);
	    d[k] -= dm[k];
	}

//	distance = ld.getDistance(d, p1M, p2.M, &s0, &t0);
	distance = getDistance(d, p1M, p2.M, &s0, &t0);
	if (distance < 1) {
	    overlap = true;
	    
//	    cout << "Distance " << distance << endl;
//	    cout << "Parameters " << s0 << " " << t0 << endl;
	}
    }
  
    return overlap;
}


void hYr::randomVector(double v[3]) {
  double rana, ranb, ransq, factor;

  ransq=2;
  while(ransq >= 1){
    rana = 1-2*ran2(iran);
    ranb = 1-2*ran2(iran);
    ransq=rana*rana+ranb*ranb;
  }
  factor = 2*sqrt(1-ransq);
  v[0] = rana*factor;
  v[1] = ranb*factor;
  v[2] = 1-2*ransq;
  
}


void hYr::nematicOPconfig() {

  double xx = -0.5;
  double yy = -0.5;
  double zz = -0.5;
  double xy = 0;
  double xz = 0;
  double yz = 0;
  double a0,a1,a2,sq;
  double qq,rr,theta;
  double z1,z2,z3;
  double onethird = 1.0/3.0;


  /* orientational tensor elements */

  double tensor_norm = 1.5*Nrinv*ARInv*ARInv;
  for (int i=0; i<N_r; i++)
    {
      xx+= tensor_norm*rods[i].M[0]*rods[i].M[0];
      yy+= tensor_norm*rods[i].M[1]*rods[i].M[1];
      zz+= tensor_norm*rods[i].M[2]*rods[i].M[2];
      xy+= tensor_norm*rods[i].M[0]*rods[i].M[1];
      xz+= tensor_norm*rods[i].M[0]*rods[i].M[2];
      yz+= tensor_norm*rods[i].M[1]*rods[i].M[2];
    }


  /* solve cubic equation */
  a2 = -xx - yy - zz;
  a1 = xx*yy + xx*zz + yy*zz - xy*xy - xz*xz - yz*yz;
  a0 = xz*xz*yy + xx*yz*yz + xy*xy*zz - 2*xy*xz*yz - xx*yy*zz;
  qq = 3*a1-a2*a2; qq/=9.0;
  rr = 9*a2*a1-27*a0-2*a2*a2*a2; rr/=54.0;

  sq=sqrt(-qq); theta=acos( rr/(-qq*sq) ); 

  z2 = 2*sq*cos( (theta+4*M_PI)*onethird ) - onethird*a2;


  z3 = 2*sq*cos( theta*onethird ) - onethird*a2;
  z1 = 2*sq*cos( (theta+2*M_PI)*onethird ) - onethird*a2;

  if(z1>z3 || z2>z3) {
      cerr << "ERROR: eigenvalues not ordered" << endl;
      exit(8);
  }

   currOP = -2*z2;

}

void hYr::locnematicOPconfig() {

  double a0,a1,a2,sq;
  double qq,rr,theta;
  double z1,z2,z3;
  double onethird = 1.0/3.0;
  int i;
  double val;


for (int pNr=0; pNr<N_r; pNr++){
  double xx = -0.5;
  double yy = -0.5;
  double zz = -0.5;
  double xy = 0;
  double xz = 0;
  double yz = 0;
  int pN = pNr*N_l;

  /* orientational tensor elements */

  double tensor_norm = 1.5*ARInv*ARInv/nbn[pNr];
  for (int j=0; j<nbn[pNr]; j++)
    {
      i=nlist[pN+j];
      xx+= tensor_norm*rods[i].M[0]*rods[i].M[0];
      yy+= tensor_norm*rods[i].M[1]*rods[i].M[1];
      zz+= tensor_norm*rods[i].M[2]*rods[i].M[2];
      xy+= tensor_norm*rods[i].M[0]*rods[i].M[1];
      xz+= tensor_norm*rods[i].M[0]*rods[i].M[2];
      yz+= tensor_norm*rods[i].M[1]*rods[i].M[2];
    }


  /* solve cubic equation */
  a2 = -xx - yy - zz;
  a1 = xx*yy + xx*zz + yy*zz - xy*xy - xz*xz - yz*yz;
  a0 = xz*xz*yy + xx*yz*yz + xy*xy*zz - 2*xy*xz*yz - xx*yy*zz;
  qq = 3*a1-a2*a2; qq/=9.0;
  rr = 9*a2*a1-27*a0-2*a2*a2*a2; rr/=54.0;

  sq=sqrt(-qq); theta=acos( rr/(-qq*sq) ); 

  z2 = 2*sq*cos( (theta+4*M_PI)*onethird ) - onethird*a2;


  z3 = 2*sq*cos( theta*onethird ) - onethird*a2;
  z1 = 2*sq*cos( (theta+2*M_PI)*onethird ) - onethird*a2;

  if(z1>z3 || z2>z3) {
      cerr << "ERROR: eigenvalues not ordered" << endl;
      exit(8);
  }

   val = -2*z2;
   if (val <= 0.3 ) locOP[pNr] = 1;
   else locOP[pNr] = 0;
}

}


void hYr::usage(){

  cerr << "Usage is hYrExe [-n -f<ConfigFile>] -p<ParameterFile> \n";
  exit(8);

}

double hYr::integrate(double funcin[], int Nintin){

double integral = 0;

    for (int i=0; i<Nintin+1; i++){
	integral += funcin[i];
   }
  
    return integral;

}


double hYr::calculateEnergy(int pNr){

    double energy=0.0;

    for (int i=pNr; i<pNr+1; i++){
      int iN = i*N_l;

#pragma omp parallel num_threads(2)
{

//   lineDist ldl;
    double energypair=0.0;
    double v0,v1,v2,v[3],w[3],w2[3], sh[3],vr2[3];
    double delx[3], dist, d[3], dm[3];
    double psi0,psi1, p1R[3], p1M[3], p1P, save, vabs, s0,t0,p2P,p2R[3],p2M[3];

    double** rh1 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh1[l] = new double[3];
    }
    double** rh2 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh2[l] = new double[3];
    }

    double *func = new double[Nint+1]();
    double *func2 = new double[Nint+1]();
    int i2, ch=0;


#pragma omp for
   for (int n=0; n<nbn[i]; n++){
	i2=nlist[iN+n];
if(pbc==2){
ch = int((rods[i].R[2] - rods[i2].R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*rods[i].R[1];
	p1R[1] = -i2phi*ch*rods[i].R[0];
	p1R[2] = rods[i].R[2];

	save = rods[i].M[0];
	p1M[0] = i2phi*ch*rods[i].M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = rods[i].M[2];

	v0=rods[i].M[1]*rods[i2].M[2]-rods[i].M[2]*rods[i2].M[1];
	v1=rods[i].M[2]*rods[i2].M[0]-rods[i].M[0]*rods[i2].M[2];
	v2=rods[i].M[0]*rods[i2].M[1]-rods[i].M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;

	}

	w[0]=(rods[i].M[1]*v[2]-rods[i].M[2]*v[1])/L;
	w[1]=(rods[i].M[2]*v[0]-rods[i].M[0]*v[2])/L;
	w[2]=(rods[i].M[0]*v[1]-rods[i].M[1]*v[0])/L;

	psi0 = acos((w[0]*rods[i].M[1]-w[1]*rods[i].M[0])/L);

	v0=p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])/L;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])/L;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])/L;

	psi1 = acos((w[0]*p1M[1]-w[1]*p1M[0]))/L;

	p1P = rods[i].P+psi1-psi0;

}
else{
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
}	
}
else{
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
}	

    for (int k=0; k<3; k++){
        d[k] = p1R[k] - p2R[k];
        if (fabs(d[k]) > box->halfx[k])
            d[k] -= copysign(box->x[k], d[k]);
	dm[k] = 0.5*(p1M[k] - p2M[k]);
	d[k] -= dm[k];

    }

//        dist = ld.getDistance(d, p1M, p2M, &s0, &t0);
        dist = getDistance(d, p1M, p2M, &s0, &t0);

        if(dist<dpcof2){

	for (int m=0; m<3; m++){
		sh[m]=d[m]+dm[m]+p2R[m];

	}

	v0=-p1M[2];
	v1=0.0;
	v2=p1M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	v[0]=v0*vabs;
	v[1]=v1*vabs;
	v[2]=v2*vabs;

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])*Linv;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])*Linv;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])*Linv;

	v0=-p2M[2];
	v1=0.0;
	v2=p2M[0];

	if (v0==0&&v1==0&&v2==0) {
	  vr2[0]=0.0;
	  vr2[1]=0.0;
	  vr2[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	vr2[0]=v0*vabs;
	vr2[1]=v1*vabs;
	vr2[2]=v2*vabs;

	}

	w2[0]=(p2M[1]*vr2[2]-p2M[2]*vr2[1])*Linv;
	w2[1]=(p2M[2]*vr2[0]-p2M[0]*vr2[2])*Linv;
	w2[2]=(p2M[0]*vr2[1]-p2M[1]*vr2[0])*Linv;

	
	double tj=-0.5;
	for (int j=0; j<Nint+1; j++){
//	  double tj=0.5*(-1.0+2.0/Nint*j);
	  double tkj=tkwave*tj;
	  double c1=Dhalf*cos(tkj+p1P);
	  double c2=Dhalf*cos(tkj+p2P);
	  double s1=Dhalf*sin(tkj+p1P);
	  double s2=Dhalf*sin(tkj+p2P);
		for (int k=0; k<3; k++){
			rh1[j][k]=sh[k]+tj*p1M[k];
			rh1[j][k]+=(c1*v[k]+s1*w[k]);

			rh2[j][k]=p2R[k]+tj*p2M[k];
			rh2[j][k]+=(c2*vr2[k]+s2*w2[k]);
		}
	    tj += NintInv;
	}




	for (int j=0; j<Nint+1; j++){
	     double rh1j0=rh1[j][0];
	     double rh1j1=rh1[j][1];
	     double rh1j2=rh1[j][2];
		for (int k=0; k<Nint+1; k++){
		     double drx = rh2[k][0]-rh1j0;
		     double dry = rh2[k][1]-rh1j1;
		     double drz = rh2[k][2]-rh1j2;
		    dist = (drx*drx+dry*dry+drz*drz);

		if(dist < (cofka2)){
		  dist=sqrt(dist);
			func[k]=exp(-ka*dist);
			func[k]/=dist;
			func[k]-=shift;
		}
		else func[k]=0.0;


		}

//		func2[j]=integrate(func,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func[in];
   }
		func2[j]=integral;
  
	}

//	energypair=integrate(func2,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func2[in];
   }
		energypair=integral;
	
#pragma omp critical
{
	energy+=energypair;
}

    }
    }
delete[] func;
func = NULL;
delete[] func2;
func2 = NULL;

for (int i = 0; i < Nint+1; ++i)
    delete [] rh1[i];
delete [] rh1;
for (int i = 0; i < Nint+1; ++i)
    delete [] rh2[i];
delete [] rh2;

    }
}
    energy *= tzl;//energy*0.25*tempr*Z*Z*lam;
    return energy;

}



void hYr::calculateVirial(){

    double energy=0;
      for (int i=0; i<N_r; i++){
	int iN = i*N_l;
#pragma omp parallel num_threads(2)
{

//    lineDist ldl;
    double energypair=0.0;
    double v0,v1,v2,v[3],w[3],w2[3], sh[3],vr2[3];
    double delx[3], dist, d[3], dm[3];
    double psi0, psi1, p1R[3], p1M[3], p1P, save, vabs, s0,t0,p2P, p2R[3],p2M[3];

    double** rh1 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh1[l] = new double[3];
    }
    double** rh2 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh2[l] = new double[3];
    }

    double *func = new double[Nint+1]();
    double *func2 = new double[Nint+1]();
    int i2, ch=0;


#pragma omp for
   for (int n=0; n<nbn[i]; n++){
	i2=nlist[iN+n];
	if (i2 < i) continue;

if(pbc==2){
ch = int((rods[i].R[2] - rods[i2].R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*rods[i].R[1];
	p1R[1] = -i2phi*ch*rods[i].R[0];
	p1R[2] = rods[i].R[2];

	save = rods[i].M[0];
	p1M[0] = i2phi*ch*rods[i].M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = rods[i].M[2];

	v0=rods[i].M[1]*rods[i2].M[2]-rods[i].M[2]*rods[i2].M[1];
	v1=rods[i].M[2]*rods[i2].M[0]-rods[i].M[0]*rods[i2].M[2];
	v2=rods[i].M[0]*rods[i2].M[1]-rods[i].M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{
	  vabs=v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;
	}

	w[0]=(rods[i].M[1]*v[2]-rods[i].M[2]*v[1])/L;
	w[1]=(rods[i].M[2]*v[0]-rods[i].M[0]*v[2])/L;
	w[2]=(rods[i].M[0]*v[1]-rods[i].M[1]*v[0])/L;

	psi0 = acos((w[0]*rods[i].M[1]-w[1]*rods[i].M[0])/L);

	v0=p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0)  {
	    v[0]=0.0;
	    v[1]=0.0;
	    v[2]=0.0;
	}
	else{
	    vabs=v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;
	}
	
	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])/L;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])/L;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])/L;

	psi1 = acos((w[0]*p1M[1]-w[1]*p1M[0]))/L;

	p1P = rods[i].P+psi1-psi0;
}
else{
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
	
}	
}
else{
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
	
}	

    for (int k=0; k<3; k++){
        d[k] = p1R[k] - p2R[k];
        if (fabs(d[k]) > box->halfx[k])
            d[k] -= copysign(box->x[k], d[k]);
	dm[k] = 0.5*(p1M[k] - p2M[k]);
	d[k] -= dm[k];
    }

//        dist = ldl.getDistance(d, p1M, p2M, &s0, &t0);
        dist = getDistance(d, p1M, p2M, &s0, &t0);

        if(dist<dpcof2){


	for (int m=0; m<3; m++){
		sh[m]=d[m]+dm[m]+p2R[m];

	}

	v0=-p1M[2];
	v1=0.0;
	v2=p1M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	v[0]=v0*vabs;
	v[1]=v1*vabs;
	v[2]=v2*vabs;

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])*Linv;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])*Linv;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])*Linv;

	v0=-p2M[2];
	v1=0.0;
	v2=p2M[0];

	if (v0==0&&v1==0&&v2==0) {
	  vr2[0]=0.0;
	  vr2[1]=0.0;
	  vr2[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	vr2[0]=v0*vabs;
	vr2[1]=v1*vabs;
	vr2[2]=v2*vabs;
	}

	w2[0]=(p2M[1]*vr2[2]-p2M[2]*vr2[1])*Linv;
	w2[1]=(p2M[2]*vr2[0]-p2M[0]*vr2[2])*Linv;
	w2[2]=(p2M[0]*vr2[1]-p2M[1]*vr2[0])*Linv;

	double tj = -0.5;
	for (int j=0; j<Nint+1; j++){
//	  double tj=0.5*(-1.0+2.0/Nint*j);
	  double tkj=tkwave*tj;
	  double c1=Dhalf*cos(tkj+p1P);
	  double c2=Dhalf*cos(tkj+p2P);
	  double s1=Dhalf*sin(tkj+p1P);
	  double s2=Dhalf*sin(tkj+p2P);
		for (int k=0; k<3; k++){
			rh1[j][k]=sh[k]+tj*p1M[k];
			rh1[j][k]+=(c1*v[k]+s1*w[k]);

			rh2[j][k]=p2R[k]+tj*p2M[k];
			rh2[j][k]+=(c2*vr2[k]+s2*w2[k]);
		}
	    tj += NintInv;
	}

	for (int j=0; j<Nint+1; j++){
	    double rh1j0=rh1[j][0];
	    double rh1j1=rh1[j][1];
	    double rh1j2=rh1[j][2];
		for (int k=0; k<Nint+1; k++){
		    double drx = rh2[k][0]-rh1j0;
		    double dry = rh2[k][1]-rh1j1;
		    double drz = rh2[k][2]-rh1j2;
		    dist = (drx*drx+dry*dry+drz*drz);

		if(dist < cofka2){
		    dist=sqrt(dist);
			func[k]=exp(-ka*dist);
			func[k]/=dist;
//			func[k]*=(ka+1.0/dist1);		//between rods wrong: dist || force ?!
			func[k]*=(ka*dist+1.0);		//between point charges
		}
		else func[k]=0;
		}

//		func2[j]=integrate(func,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func[in];
   }
		func2[j]=integral;
  

	}

//	energypair=integrate(func2,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func2[in];
   }
		energypair=integral;
//	energy+=energypair*dist;	//between rods wrong: dist || force ?!
#pragma omp critical
	energy+=energypair;			//between point charges

    }
    }
delete[] func;
func = NULL;
delete[] func2;
func2 = NULL;

for (int i = 0; i < Nint+1; ++i)
    delete [] rh1[i];
delete [] rh1;
for (int i = 0; i < Nint+1; ++i)
    delete [] rh2[i];
delete [] rh2;

}
    }

    energy=(nt+energy*tzl*tempr/3.0)/box->V;
    Ecalc=energy;

}


void hYr::neighborlist(){

  double dist, d[3], p1R[3], p1M[3], dm[3];
  double s0, t0, save;
  Rod p1, p2;
  int ch=0;

	
  for (int i=0; i<N_r; i++){
    int iN = i*N_l;
    nbn[i] = 0;
    for (int j=0; j<N_l; j++){
	nlist[iN+j] = 0;
    }
  }
  
	for (int i=0; i<N_r-1; i++){
    int iN = i*N_l;
		for (int j=i+1; j<N_r; j++){

		for (int k=0; k<3; k++){
		p1.R[k] = rods[i].R[k];
		p2.R[k] = rods[j].R[k];
		p1.M[k] = rods[i].M[k];
		p2.M[k] = rods[j].M[k];
		}

if(pbc==2){
ch = int((p1.R[2] - p2.R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*p1.R[1];
	p1R[1] = -i2phi*ch*p1.R[0];
	p1R[2] = p1.R[2];

	save = p1.M[0];
	p1M[0] = i2phi*ch*p1.M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = p1.M[2];


}
else{
	for (int k =0; k<3; k++){
		p1R[k] = p1.R[k];
		p1M[k] = p1.M[k];
	}
}	
}
else{
	for (int k =0; k<3; k++){
		p1R[k] = p1.R[k];
		p1M[k] = p1.M[k];
	}
}	

    for (int k=0; k<3; k++){
	d[k] = p1R[k] - p2.R[k];
	if (fabs(d[k]) > box->halfx[k]) 
	    d[k] -= copysign(box->x[k], d[k]);
    } 
    
    if ((d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) < (AR+1+cofka+4*maxDisplace)*(AR+1+cofka+4*maxDisplace)){
	
	for (int k=0; k<3; k++){
	    dm[k] = 0.5*(p1M[k] - p2.M[k]);
	    d[k] -= dm[k];
	}
	
//	dist = ld.getDistance(d, p1M, p2.M, &s0, &t0);
	dist = getDistance(d, p1M, p2.M, &s0, &t0);

			if(dist<(dpcof+4*maxDisplace)*(dpcof+4*maxDisplace)){
				nlist[iN+nbn[i]]=j;

				nlist[j*N_l+nbn[j]]=i;

				nbn[i]+=1;

				nbn[j]+=1;
			}

	}

		}

	}
}


void hYr::compressBox(double oV, double nV, int dir){

  double noo = nV/oV;
  double snoo = sqrt(noo);
if(dir==1){
	for (int i=0; i<N_r; i++){
	for (int m=0; m<2; m++){
		rods[i].R[m]*=snoo;
	}
	}

	for (int m=0; m<2; m++){
		box->x[m]*=snoo;
	}
}
else{
	for (int i=0; i<N_r; i++){
	for (int m=2; m<3; m++){
		rods[i].R[m]*=noo;
	}
	}

	for (int m=2; m<3; m++){
		box->x[m]*=noo;
	}
}

	delete[] cells;
	cells = NULL;
	delete[] Tested;
	Tested = NULL;
	delete[] PBuffer;
	PBuffer = NULL;
    for (int i=0; i<N_r; i++){
	delete[] rods[i].pieces;
	rods[i].pieces = NULL;
    }
	

	setupList();
	ka=sqrt(4.0*M_PI*lam*(Z*N_r/box->V+2.0*cs));
	cof=pow((box->V)/(N_r)/M_PI/AR*4.0,1.0/2.0)*ka;
	cofka=cof/ka;
        dpcof=D+cofka;
	dpcof2=dpcof*dpcof;
	cofka2=cofka*cofka;

}

void hYr::compressBoxV(double oV, double nV, int dir){

  double noo = nV/oV;
  double snoo = sqrt(noo);
if(dir==1){
	for (int i=0; i<N_r; i++){
	for (int m=0; m<2; m++){
		rods[i].R[m]*=snoo;
	}
	rods[i].R[2]/=noo;
	}

	for (int m=0; m<2; m++){
		box->x[m]*=snoo;
	}
	box->x[2]/=noo;
}
else{
	for (int i=0; i<N_r; i++){
	for (int m=2; m<3; m++){
		rods[i].R[m]*=noo;
	}
	for (int m=0; m<2; m++){
		rods[i].R[m]/=snoo;
	}
	}

	for (int m=2; m<3; m++){
		box->x[m]*=noo;
	}
	for (int m=0; m<2; m++){
		box->x[m]/=snoo;
	}
}

	delete[] cells;
	cells = NULL;
	delete[] Tested;
	Tested = NULL;
	delete[] PBuffer;
	PBuffer = NULL;
    for (int i=0; i<N_r; i++){
	delete[] rods[i].pieces;
	rods[i].pieces = NULL;
    }
	

	setupList();
	ka=sqrt(4.0*M_PI*lam*(Z*N_r/box->V+2.0*cs));

}

// Measures the cluster size distribution and finds the size of the largest
// cluster, the number of contacts in the largest cluster, and the number of
// contacts in the whole system.

void hYr::measureClusterDist(){
    
    long ilargest;   // Index of largest cluster
    Rod *curPart;


    findClusters();

    ilargest = 0;
    nlargest = 0;
    allContacts = 0;

    for (int i=0; i<N_r+1; i++){
	clusterDist[i] = 0;
    }

    for (int i=0; i<N_r; i++){
	clusterDist[clusters[i].npart]++;
        if (clusters[i].npart > nlargest) {
           nlargest = clusters[i].npart;
           ilargest = i;
        }
        allContacts += rods[i].nneb;
    }

    clustContacts = 0;
    curPart = clusters[ilargest].firstParticle;
    while (curPart) {
       clustContacts += curPart->nneb;
       curPart = curPart->nextCl;
    }
}


void hYr::findClusters(){
    
    Cluster *curCluster;
    
    clearClusters();

    //Teilchen sind Nachbarn, wenn sie naeher als D+1/ka liegen und 
    
    for (int i=0; i<N_r; i++) {
        rods[i].nneb = 0;
	rods[i].insertToCluster(clusters[i]);
        // This makes every rod a one-particle cluster for now.
    }
    for (int i=0; i<N_r; i++){
	for (int j=i+1; j<N_r; j++){
	    if (Neighbours(rods[i], rods[j])){
		    rods[i].neighbours[rods[i].nneb] = &rods[j];
		    rods[i].nneb ++;
		    rods[j].neighbours[rods[j].nneb] = &rods[i];
		    rods[j].nneb ++;
                    // Add i and j mutually to each other's neighbour lists.
	    }
	}
    }

    // Cluster finden

    for (int i=0; i<N_r; i++){
	for (int j=0; j<rods[i].nneb; j++){
	    curCluster = rods[i].neighbours[j]->cluster;
	    if (rods[i].cluster != curCluster){
		while(curCluster->firstParticle){
		    curCluster->firstParticle->
			moveBetweenClusters(*curCluster, (* rods[i].cluster));
		}
	    }
	}
    }
	

}


// Function to detect percolation amongst the rods in the present configuration.
// MAM
int hYr::percolate() {

    double d[3], dm[3], p1R[3];
    double dCont;         // Squared distance between rods in contiguous cluster
    double s0, t0;        // Dummy -- not used
    long placed;          // Number of particles already placed in cluster
    int percolating;      // 1 = a percolating cluster exists
    Cluster *curCluster;  // Current cluster being tested for percolation
    Rod *curRod;          // Current rod being placed in cluster
    Rod *explore;         // Rod being exlpored for neighbours
    Rod *newest;          // Most recent rod to be added to copy of cluster
    Rod *rod1, *rod2;     // A pair of rods being checked for distance apart
	int ch=0;

    findClusters();

    percolating = 0;
    for (int i=0; i<N_r; i++) rods[i].done=false;

    // Test each cluster in turn for percolation
    for (curCluster=clusters; curCluster<clusters+N_r; curCluster++) {
       if (curCluster->npart < 2) continue;

       // Build up a contiguous copy of the cluster without the boundary conditions.

       // Place first particle at the origin.
       curRod = curCluster->firstParticle;
       curRod->copyR[0] = curRod->copyR[1] = curRod->copyR[2] = 0.0;
       curRod->done = true;
       placed = 1;
       newest = curRod;

             // Iteratively add neighbours of particles that have already been placed.
       explore = curRod;
       while (placed < curCluster->npart) {
          for (int i=0; i<explore->nneb; i++) {
             curRod = explore->neighbours[i];
             if (!curRod->done) {

                // Find nearest image of neighbour.

if(pbc==2){
ch = int((curRod->R[2] - explore->R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*curRod->R[1];
	p1R[1] = -i2phi*ch*curRod->R[0];
	p1R[2] = curRod->R[2];


}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = curRod->R[i2];
		
	}
}	
}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = curRod->R[i2];
		
	}
}	

                for (int k=0; k<3; k++){
                   d[k] = p1R[k] - explore->R[k];    
                   if (fabs(d[k]) > box->halfx[k]) d[k] -= copysign(box->x[k], d[k]);
                   curRod->copyR[k] = explore->copyR[k] + d[k];
                } 

                curRod->done = true;
                newest->nextPlace = curRod;
                newest = curRod;
                placed++;
             }
          }
          explore = explore->nextPlace;
       }

       // Now check whether any two rods that are far apart in the contiguous
       // cluster are neighbours when the periodic boundary conditions are
       // turned on.  If so, the distant rods are connected directly to
       // each other's periodic images as well as indirectly through the
       // contiguous cluster, and the cluster therefore percolates.  If no
       // such pair of rods exists in the cluster then the cluster does not
       // span the system.

       for (rod1=curCluster->firstParticle; rod1; rod1=rod1->nextCl) {

          // Loop over neighbours defined WITH boundary conditions
          for (int i=0; i<rod1->nneb; i++) {
             rod2 = rod1->neighbours[i];

             // Distance in contiguous cluster (no boundary conditions)
             for (int k=0; k<3; k++){
                dm[k]=0.5*(rod1->M[k] - rod2->M[k]);
                d[k] = rod1->copyR[k] - rod2->copyR[k] - dm[k];
             } 
//             dCont = ld.getDistance(d, rod1->M, rod2->M, &s0, &t0);
             dCont = getDistance(d, rod1->M, rod2->M, &s0, &t0);
             // The +0.1 is for numerical hygeine.  If the rods are not neighbours they
             // will be MUCH more than D+1/ka apart.
             if (dCont > (D+1.0/ka)*(D+1.0/ka)+0.1) {
                percolating = 1;
                goto escape;
             }

          }
       }
    }

    escape:
    return percolating;
}


// closest approach of axes
void hYr::measureDistances(){
    
    double distance, d[3], dm[3], p1R[3], p1M[3], save;
    double s0, t0;
    double centreSep, parSep, perpSep;  // Separation of centres of the rods
    double dot;        // Dot product of rod directions
    int bin, binpar, binperp;
	int ch=0;

    for (int i=0; i<nbins; i++){
	distances[i] = 0;
    }
    for (int i=0; i<N_r; i++){
	for (int j=i+1; j<N_r; j++){
            centreSep = 0.0;
	    parSep = 0.0;
	    perpSep = 0.0;

if(pbc==2){
ch = int((rods[i].R[2] - rods[j].R[2])*box->halfxInv[2]);
if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*rods[i].R[1];
	p1R[1] = -i2phi*ch*rods[i].R[0];
	p1R[2] = rods[i].R[2];

	save = rods[i].M[0];
	p1M[0] = i2phi*ch*rods[i].M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = rods[i].M[2];
	

}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = rods[i].R[i2];
		p1M[i2] = rods[i].M[i2];
	}
}	
}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = rods[i].R[i2];
		p1M[i2] = rods[i].M[i2];
	}
}	

	    for (int k=0; k<3; k++){
		d[k] = p1R[k] - rods[j].R[k];    
		if (fabs(d[k]) > box->halfx[k]) 
		    d[k] -= copysign(box->x[k], d[k]);
                centreSep += d[k]*d[k];
		parSep += d[k]*p1M[k];
		dm[k] = 0.5*(p1M[k] - rods[j].M[k]);
		d[k] -= dm[k];
	    }
	
//	    distance = ld.getDistance(d, p1M, rods[j].M, &s0, &t0);
	    distance = getDistance(d, p1M, rods[j].M, &s0, &t0);
	    distance = sqrt(distance);
	    bin = int(distance/binWidth);
	    if (bin <= nbins) distances[bin]++;
            if (distance < A) {
               s0 = 1.0-2.0*fabs(s0 - 0.5);
               if (s0 < 1.0e-8) {
                  endCont++;
               } else {
                  bin = (long)(s0*sbins);
                  sHist[bin]++;
               }
               t0 = 1.0-2.0*fabs(t0 - 0.5);
               if (t0 < 1.0e-8) {
                  endCont++;
               } else {
                  bin = (long)(t0*sbins);
                  sHist[bin]++;
               }
            }

	    parSep /= L;
	    binpar = int(fabs(parSep)/binWidth);
	    perpSep = sqrt(centreSep-parSep*parSep);
	    binperp = int(perpSep/binWidth);
            centreSep = sqrt(centreSep);
            bin = int(centreSep/binWidth);
            if (bin <= nbins) {
               gofr[bin]++;
               dot = (p1M[0]*rods[j].M[0] + p1M[1]*rods[j].M[1] +
                  p1M[2]*rods[j].M[2]) / AR2;
               g2ofr[bin] += 1.5*dot*dot - 0.5;
	       gofrpp[binpar][binperp]++;
	       g2ofrpp[binpar][binperp] += 1.5*dot*dot - 0.5;
	       
            }
	}
    } 
    
}


bool hYr::Neighbours(Rod &p1, Rod &p2) {

    bool neigh = false;
    double d[3], dm[3], d2, p1R[3], p1M[3], save;
    double s0, t0;    // Dummy -- not used
	int ch=0;

if(pbc==2){	
ch = int((p1.R[2] - p2.R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*p1.R[1];
	p1R[1] = -i2phi*ch*p1.R[0];
	p1R[2] = p1.R[2];

	save = p1.M[0];
	p1M[0] = i2phi*ch*p1.M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = p1.M[2];
	

}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = p1.R[i2];
		p1M[i2] = p1.M[i2];
	}
}	
}
else{
	for (int i2 =0; i2<3; i2++){
		p1R[i2] = p1.R[i2];
		p1M[i2] = p1.M[i2];
	}
}	

    for (int j=0; j<3; j++){
	d[j] = p1R[j] - p2.R[j];
	if (fabs(d[j]) > box->halfx[j]) 
	    d[j] -= copysign(box->x[j], d[j]);
    } 
    
    if ((d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) < (AR+1+1.0/ka)*(AR+1+1.0/ka)){
	
	for (int k=0; k<3; k++){
	    dm[k] = 0.5*(p1M[k] - p2.M[k]);
	    d[k] -= dm[k];
	}
	
//	d2 = ld.getDistance(d, p1M, p2.M, &s0, &t0);
	d2 = getDistance(d, p1M, p2.M, &s0, &t0);

	if (d2 < (D+1.0/ka)*(D+1.0/ka)) neigh = true;
    }
    return neigh;
}   


double hYr::twistBox(double twphi){

  double Enew;
  double rad, rad1, vphi, vphi1, psi0, psi1, oldV[3],oldW[3],newV[3],newW[3];
  double philoc,save,r;
  Rod *rodscp = new Rod[N_r];

	for (int i=0; i<N_r; i++){
	  rodscp[i] = rods[i];
		philoc = twphi/box->x[2]*rods[i].R[2]+0.5*twphi;
		double cpl=cos(philoc);
		double spl=sin(philoc);
		
		save = rods[i].R[0];
		rods[i].R[0] = rods[i].R[0]*cos(philoc)-rods[i].R[1]*sin(philoc);
		rods[i].R[1] = rods[i].R[1]*cos(philoc)+save*sin(philoc);
		
		psi0=rods[i].P;	
//v'.v~
oldV[0]=-rods[i].M[2];
oldV[1]=0.0;
oldV[2]=rods[i].M[0];
r=2;
      newV[(int(r)+1)%3] = cos(philoc)*oldV[(int(r)+1)%3]-sin(philoc)*oldV[(int(r)+2)%3]; 
      newV[(int(r)+2)%3] = cos(philoc)*oldV[(int(r)+2)%3]+sin(philoc)*oldV[(int(r)+1)%3];
      newV[int(r)] = oldV[int(r)];

//v'.w~
oldW[0]=rods[i].M[0]*rods[i].M[1];
oldW[1]=-(rods[i].M[0]*rods[i].M[0]+rods[i].M[2]*rods[i].M[2]);
oldW[2]=rods[i].M[1]*rods[i].M[2];
r=2;
      newW[(int(r)+1)%3] = cos(philoc)*oldW[(int(r)+1)%3]-sin(philoc)*oldW[(int(r)+2)%3]; 
      newW[(int(r)+2)%3] = cos(philoc)*oldW[(int(r)+2)%3]+sin(philoc)*oldW[(int(r)+1)%3];
      newW[int(r)] = oldW[int(r)];
		
		save = rods[i].M[0];
		rods[i].M[0] = rods[i].M[0]*cpl-rods[i].M[1]*spl;
		rods[i].M[1] = rods[i].M[1]*cpl+save*spl;
		
		if(fabs(rods[i].M[1]) < 1e-9) {
		  rods[i].M[1] = 0.0;
		}
		
psi1=cos(psi0)*(-rods[i].M[2]*newV[0]+rods[i].M[0]*newV[2])/sqrt((rods[i].M[0]*rods[i].M[0]+rods[i].M[2]*rods[i].M[2])*(oldV[0]*oldV[0]+oldV[2]*oldV[2]));

psi1+=sin(psi0)*(-rods[i].M[2]*newW[0]+rods[i].M[0]*newW[2])/L/sqrt((rods[i].M[0]*rods[i].M[0]+rods[i].M[2]*rods[i].M[2])*(oldV[0]*oldV[0]+oldV[2]*oldV[2]));

if(fabs(psi1)>1.0) psi1=copysign(1,psi1);
if(fabs(psi0)<M_PI) psi1=copysign(acos(psi1),psi0);
else psi1=copysign(2.0*M_PI-acos(psi1),psi0);

		rods[i].P = rods[i].P+psi1-psi0;
	}

//if (twphi == M_PI/2.0) writeOutConfig("twistconf");

	    Enew =0.0;
	for (int i=0; i<N_r; i++){
	  Enew += calculateEnergyTwist(i, twphi);
	}

	for (int i=0; i<N_r; i++){
	  rods[i] = rodscp[i];
	}

	delete[] rodscp;
	rodscp = NULL;
	return Enew;

}

/// Not corrected, because not used at the moment! ///
double hYr::calculateEnergyTwist(int pNr, double twphi){

    double energy=0.0;

    for (int i=pNr; i<pNr+1; i++){

#pragma omp parallel num_threads(2)
{

    double energypair=0.0;
    double v0,v1,v2,v[3],w[3],w2[3], sh[3],vr2[3];
    double delx[3], dist, d[3], dm[3], philoc;
    double p1R[3], p1M[3], p1P, vabs, s0,t0,p2P,p2R[3],p2M[3];

    double** rh1 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh1[l] = new double[3];
    }
    double** rh2 = new double*[Nint+1];
    for (int l = 0; l < Nint+1; ++l) {
    	rh2[l] = new double[3];
    }

    double *func = new double[Nint+1];
    double *func2 = new double[Nint+1];
    for (int k=0; k<Nint+1; k++){
	    func[k]=0;
	    func2[k]=0;
		}
    int i2, intphi;


#pragma omp for
   for (int n=0; n<nbn[i]; n++){

	i2=nlist[i*N_l+n];

	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
	}
		p1P = rods[i].P;

	philoc = (twphi/box->x[2])*p1R[2]+0.5*twphi;
	philoc = fmod(philoc+9.0/4.0*M_PI,2.0*M_PI);
	intphi = int(philoc/M_PI*2.0);
	philoc = philoc+(49.0/6.0*intphi-13.0/2.0*intphi*intphi+4.0/3.0*intphi*intphi*intphi-0.5)*0.5*M_PI;
    d[0] = p1R[0] - rods[i2].R[0];
    d[1] = p1R[1] - rods[i2].R[1];
    d[2] = p1R[2] - rods[i2].R[2];
	if (fabs(d[2]) > box->halfx[2]) {
	  continue;
//	  d[2] -= box->x[2];
//	  d[1] -= sin(philoc-twphi)*p1R[1];
//	  d[0] -= cos(philoc-twphi)*p1R[0];
	  }
        if (fabs(d[0]) > box->halfx[0]) {
	  d[1] -= sin(philoc)*box->x[1]*copysign(1,d[0]);
	  d[0] -= cos(philoc)*box->x[0]*copysign(1,d[0]);
	}
	if (fabs(d[1]) > box->halfx[1]) {
	  d[0] -= -sin(philoc)*box->x[0]*copysign(1,d[1]);
	  d[1] -= cos(philoc)*box->x[1]*copysign(1,d[1]);
	    if (fabs(d[0]) > box->halfx[0]) {
		d[1] -= sin(philoc)*box->x[1]*copysign(1,d[0]);
		d[0] -= cos(philoc)*box->x[0]*copysign(1,d[0]);
	    }
	} 



    for (int k=0; k<3; k++){
	dm[k] = 0.5*(p1M[k] - rods[i2].M[k]);
	d[k] -= dm[k];

    }
#pragma omp critical
//        dist = ld.getDistance(d, p1M, rods[i2].M, &s0, &t0);
	dist = getDistance(d, p1M, rods[i2].M, &s0, &t0);

        dist = sqrt(dist);

        if(dist<(D+cof/ka)){

	sh[0]=0.0;
	sh[1]=0.0;
	sh[2]=0.0;
	for (int m=0; m<3; m++){
	  sh[m] = d[m]+dm[m]-p1R[m]+rods[i2].R[m];
	  sh[m] /= -box->x[m];
//		delx[m]=p1R[m]-rods[i2].R[m];
	}
/*	if(fabs(delx[0]) > box->halfx[0]) {
	  sh[0] += cos(philoc);
	  sh[1] += sin(philoc);
	}
	if(fabs(delx[1]) > box->halfx[1]) {
          sh[1] += cos(philoc);
          sh[0] += -sin(philoc);
        } 
*/	


	v0=-p1M[2];//p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=0.0;//p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=p1M[0];//p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];
	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{
	  vabs = sqrt(v0*v0+v1*v1+v2*v2);
	v[0]=v0/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	v[1]=v1/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	v[2]=v2/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	}
	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])/L;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])/L;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])/L;

	v0=-rods[i2].M[2];//p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=0.0;//p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=rods[i2].M[0];//p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];
	if (v0==0&&v1==0&&v2==0) {
	  vr2[0]=0.0;
	  vr2[1]=0.0;
	  vr2[2]=0.0;
	}
	else{
	  vabs = sqrt(v0*v0+v1*v1+v2*v2);
	vr2[0]=v0/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	vr2[1]=v1/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	vr2[2]=v2/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	}
	w2[0]=(rods[i2].M[1]*vr2[2]-rods[i2].M[2]*vr2[1])/L;
	w2[1]=(rods[i2].M[2]*vr2[0]-rods[i2].M[0]*vr2[2])/L;
	w2[2]=(rods[i2].M[0]*vr2[1]-rods[i2].M[1]*vr2[0])/L;

	p1P=p1P;//acos((cos(p1P)*(-v[0]*p1M[2]+v[2]*p1M[0])+sin(p1P)*(v[0]*p1M[0]*p1M[1]-v[1]*(p1M[0]*p1M[0]+p1M[2]*p1M[2])+v[2]*p1M[1]*p1M[2]))/sqrt(p1M[0]*p1M[0]+p1M[2]*p1M[2]));
	p2P=rods[i2].P;//acos((cos(rods[i2].P)*(-v[0]*rods[i2].M[2]+v[2]*rods[i2].M[0])+sin(rods[i2].P)*(v[0]*rods[i2].M[0]*rods[i2].M[1]-v[1]*(rods[i2].M[0]*rods[i2].M[0]+rods[i2].M[2]*rods[i2].M[2])+v[2]*rods[i2].M[1]*rods[i2].M[2]))/sqrt(rods[i2].M[0]*rods[i2].M[0]+rods[i2].M[2]*rods[i2].M[2]));
//t = clock();

//#pragma omp critical
	for (int j=0; j<Nint+1; j++){
		for (int k=0; k<3; k++){
			rh1[j][k]=p1R[k]-sh[k]*box->x[k]+0.5*(-1.0+2.0/Nint*j)*p1M[k];
			rh1[j][k]+=0.5*D*(cos(kwave*(-1.0+2.0/Nint*j)+p1P)*v[k]+sin(kwave*(-1.0+2.0/Nint*j)+p1P)*w[k]);

			rh2[j][k]=rods[i2].R[k]+0.5*(-1.0+2.0/Nint*j)*rods[i2].M[k];
			rh2[j][k]+=0.5*D*(cos(kwave*(-1.0+2.0/Nint*j)+p2P)*vr2[k]+sin(kwave*(-1.0+2.0/Nint*j)+p2P)*w2[k]);
		}
	}

	/*v0=p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = sqrt(v0*v0+v1*v1+v2*v2);
	v[0]=v0/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	v[1]=v1/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));
	v[2]=v2/vabs;///sqrt(pow(v0,2)+pow(v1,2)+pow(v2,2));

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])/L;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])/L;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])/L;

	w2[0]=(rods[i2].M[1]*v[2]-rods[i2].M[2]*v[1])/L;
	w2[1]=(rods[i2].M[2]*v[0]-rods[i2].M[0]*v[2])/L;
	w2[2]=(rods[i2].M[0]*v[1]-rods[i2].M[1]*v[0])/L;

//t = clock();
	for (int j=0; j<Nint+1; j++){
		for (int k=0; k<3; k++){
			rh1[j][k]=p1R[k]-sh[k]*box->x[k]+0.5*(-1.0+2.0/Nint*j)*p1M[k];
			rh1[j][k]+=0.5*D*(cos(kwave*(-1.0+2.0/Nint*j)+p1P)*v[k]+sin(kwave*(-1.0+2.0/Nint*j)+p1P)*w[k]);

			rh2[j][k]=rods[i2].R[k]+0.5*(-1.0+2.0/Nint*j)*rods[i2].M[k];
			rh2[j][k]+=0.5*D*(cos(kwave*(-1.0+2.0/Nint*j)+rods[i2].P)*v[k]+sin(kwave*(-1.0+2.0/Nint*j)+rods[i2].P)*w2[k]);
		}
	}*/
//t = clock()-t;
//cout<<"seconds1 "<< setiosflags(ios::fixed)<<setprecision(20)<<((float)t)/CLOCKS_PER_SEC<<endl;
//t = clock();
	for (int j=0; j<Nint+1; j++){
//s = clock();
		for (int k=0; k<Nint+1; k++){

		 dist = sqrt(pow(rh2[k][0]-rh1[j][0],2)+pow(rh2[k][1]-rh1[j][1],2)+pow(rh2[k][2]-rh1[j][2],2));

		 if(dist < cof/ka){
			func[k]=exp(-ka*dist);
			func[k]/=dist;
			func[k]-=exp(-cof)/cof*ka;
		 }
		 else func[k]=0.0;
//cout<<"dist "<<dist<<endl;
//cout<<k<<" "<<func[k]<<" "<<sqrt(pow(rh2[k][0]-rh1[j][0],2)+pow(rh2[k][1]-rh1[j][1],2)+pow(rh2[k][2]-rh1[j][2],2))<<endl;

		}
//cout<<endl;
//s = clock()-s;
//cout<<"seconds3 "<< setiosflags(ios::fixed)<<setprecision(20)<<((float)s)/CLOCKS_PER_SEC<<endl;

//s = clock();
//		func2[j]=integrate(func,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func[in];
   }
		func2[j]=integral;
  
//s = clock()-s;
//cout<<"seconds4 "<< setiosflags(ios::fixed)<<setprecision(20)<<((float)s)/CLOCKS_PER_SEC<<endl;


	}
//t = clock()-t;
//cout<<"seconds2 "<< setiosflags(ios::fixed)<<setprecision(20)<<((float)t)/CLOCKS_PER_SEC<<endl;

//	energypair=integrate(func2,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func2[in];
   }
		energypair=integral;
	
#pragma omp critical
	energy+=energypair;
////cout<<"here"<<energypair<<endl;

    }
    }
delete[] func;
func = NULL;
delete[] func2;
func2 = NULL;

for (int i = 0; i < Nint+1; ++i)
    delete [] rh1[i];
delete [] rh1;
for (int i = 0; i < Nint+1; ++i)
    delete [] rh2[i];
delete [] rh2;

}
    }
//cout<<energy<<endl;
//cout<< tempr<<endl;
//cout<<Z<<endl;
//cout<<lam<<endl;
//cout<<energy<<" "<<tempr<<" "<<Z<<" "<<lam<<endl;
    energy=energy*0.25*tempr*Z*Z*lam;
    return energy;

}

double hYr::calculateEnergyPair(int pNr, int pNr2){
    double energy=0.0;
//    lineDist ldl;
    double energypair=0.0;
    double v0,v1,v2,v[3],w[3],w2[3], sh[3],vr2[3],r,oldV[3],oldW[3],newV[3],newW[3];
    double delx[3], dist, d[3], dm[3], philoct, val;
    double rad, rad1, vphi, vphi1, psi0, psi1, save;
    double p1R[3], p1M[3], p1P, p1philoc, vabs, s0,t0,p2P, p2R[3], p2M[3];
    int ch=0, intphi;

    double **rh2p = new double*[Nint+2];
    for (int l = 0; l < Nint+2; l++) {
    	rh2p[l] = new double[3];
    }
    double **rh1p = new double*[Nint+2];
    for (int l = 0; l < Nint+2; l++) {
    	rh1p[l] = new double[3];
    }

    double *func = new double[Nint+2]();
    double *func2 = new double[Nint+2]();
    int i2,i;
      i2 = pNr2;
      i=pNr;
            
if(pbc==2){      
ch = int((rods[i].R[2] - rods[i2].R[2])*box->halfxInv[2]);

if (ch !=0){

// phi=pi/2:
	p1R[0] = i2phi*ch*rods[i].R[1];
	p1R[1] = -i2phi*ch*rods[i].R[0];
	p1R[2] = rods[i].R[2];

	save = rods[i].M[0];
	p1M[0] = i2phi*ch*rods[i].M[1];
	p1M[1] = -i2phi*ch*save;
	p1M[2] = rods[i].M[2];

	v0=rods[i].M[1]*rods[i2].M[2]-rods[i].M[2]*rods[i2].M[1];
	v1=rods[i].M[2]*rods[i2].M[0]-rods[i].M[0]*rods[i2].M[2];
	v2=rods[i].M[0]*rods[i2].M[1]-rods[i].M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;

	}

	w[0]=(rods[i].M[1]*v[2]-rods[i].M[2]*v[1])/L;
	w[1]=(rods[i].M[2]*v[0]-rods[i].M[0]*v[2])/L;
	w[2]=(rods[i].M[0]*v[1]-rods[i].M[1]*v[0])/L;

	psi0 = acos((w[0]*rods[i].M[1]-w[1]*rods[i].M[0])/L);

	v0=p1M[1]*rods[i2].M[2]-p1M[2]*rods[i2].M[1];
	v1=p1M[2]*rods[i2].M[0]-p1M[0]*rods[i2].M[2];
	v2=p1M[0]*rods[i2].M[1]-p1M[1]*rods[i2].M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = v0*v0+v1*v1+v2*v2;
	v[0]=v0/vabs;
	v[1]=v1/vabs;
	v[2]=v2/vabs;

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])/L;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])/L;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])/L;

	psi1 = acos((w[0]*p1M[1]-w[1]*p1M[0]))/L;

	p1P = rods[i].P+psi1-psi0;

	for (int k =0; k<3; k++){
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p2P = rods[i2].P;
		
}
else {
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
		
}
}
else {
	for (int k =0; k<3; k++){
		p1R[k] = rods[i].R[k];
		p1M[k] = rods[i].M[k];
		p2R[k] = rods[i2].R[k];
		p2M[k] = rods[i2].M[k];
	}
		p1P = rods[i].P;
		p2P = rods[i2].P;
		
}


    for (int k=0; k<3; k++){
        d[k] = p1R[k] - p2R[k];
        if (fabs(d[k]) > box->halfx[k])
            d[k] -= copysign(box->x[k], d[k]);
	dm[k] = 0.5*(p1M[k] - p2M[k]);
	d[k] -= dm[k];

    }
//        dist = ldl.getDistance(d, p1M, p2M, &s0, &t0);
	dist = getDistance(d, p1M, p2M, &s0, &t0);
        if(dist<dpcof2){

	for (int m=0; m<3; m++){
		sh[m]=d[m]+dm[m]+p2R[m];

	}

	v0=-p1M[2];
	v1=0.0;;
	v2=p1M[0];

	if (v0==0&&v1==0&&v2==0) {
	  v[0]=0.0;
	  v[1]=0.0;
	  v[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	v[0]=v0*vabs;
	v[1]=v1*vabs;
	v[2]=v2*vabs;

	}

	w[0]=(p1M[1]*v[2]-p1M[2]*v[1])*Linv;
	w[1]=(p1M[2]*v[0]-p1M[0]*v[2])*Linv;
	w[2]=(p1M[0]*v[1]-p1M[1]*v[0])*Linv;

	v0=-p2M[2];
	v1=0.0;
	v2=p2M[0];

	if (v0==0&&v1==0&&v2==0) {
	  vr2[0]=0.0;
	  vr2[1]=0.0;
	  vr2[2]=0.0;
	}
	else{

	  vabs = 1.0/sqrt(v0*v0+v1*v1+v2*v2);
	vr2[0]=v0*vabs;
	vr2[1]=v1*vabs;
	vr2[2]=v2*vabs;

	}

	w2[0]=(p2M[1]*vr2[2]-p2M[2]*vr2[1])*Linv;
	w2[1]=(p2M[2]*vr2[0]-p2M[0]*vr2[2])*Linv;
	w2[2]=(p2M[0]*vr2[1]-p2M[1]*vr2[0])*Linv;

	double tj = -0.5;
	for (int j=0; j<Nint+1; j++){
//	  double tj=0.5*(-1.0+2.0/Nint*j);
	  double tkj=tkwave*tj;
	  double c1=Dhalf*cos(tkj+p1P);
	  double c2=Dhalf*cos(tkj+p2P);
	  double s1=Dhalf*sin(tkj+p1P);
	  double s2=Dhalf*sin(tkj+p2P);
		for (int k=0; k<3; k++){
			rh1p[j][k]=sh[k]+tj*p1M[k];
			rh1p[j][k]+=(c1*v[k]+s1*w[k]);

			rh2p[j][k]=p2R[k]+tj*p2M[k];
			rh2p[j][k]+=(c2*vr2[k]+s2*w2[k]);
		}
	    tj += NintInv;
	}




	for (int j=0; j<Nint+1; j++){
	     double rh1j0=rh1p[j][0];
	     double rh1j1=rh1p[j][1];
	     double rh1j2=rh1p[j][2];
		for (int k=0; k<Nint+1; k++){
		     double drx = rh2p[k][0]-rh1j0;
		     double dry = rh2p[k][1]-rh1j1;
		     double drz = rh2p[k][2]-rh1j2;
		    dist = (drx*drx+dry*dry+drz*drz);

		if(dist < (cofka2)){
		  dist=sqrt(dist);
			func[k]=exp(-ka*dist);
			func[k]/=dist;
			func[k]-=shift;
		}
		else func[k]=0.0;
		}

//		func2[j]=integrate(func,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func[in];
   }
		func2[j]=integral;
  
	}

//	energypair=integrate(func2,Nint);
double integral = 0;

    for (int in=0; in<Nint+1; in++){
	integral += func2[in];
   }
		energypair=integral;
	
	energy+=energypair;

    }
delete[] func;
func = NULL;
delete[] func2;
func2 = NULL;

for (int i = 0; i < Nint+2; i++)
    delete [] rh1p[i];
delete [] rh1p;
for (int i = 0; i < Nint+2; i++)
    delete [] rh2p[i];
delete [] rh2p;

    energy *= tzl;
    return energy;

}

