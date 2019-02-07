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


double hYr::getDistance(double d[3], double m0[3], double m1[3], double *s0, double *t0){

  for (int i=0; i<3; i++){
    M0[i]=m0[i];
    M1[i]=m1[i];
    DD[i] = d[i];
  }

  dinitialize();

  if(ddelta*ddelta < 1e-10) {
    parallel();
  }

  else {
    nonparallel();
  }

  *s0 = smin;
  *t0 = tmin;

  return distance;

}


void hYr::dinitialize(){

//  a = M0[0]*M0[0] + M0[1]*M0[1] + M0[2]*M0[2];
//  b = - M1[0]*M0[0] - M1[1]*M0[1] - M1[2]*M0[2];
//  c = M1[0]*M1[0] + M1[1]*M1[1] + M1[2]*M1[2];
//  d = M0[0]*D[0] + M0[1]*D[1] + M0[2]*D[2];
//  e = - M1[0]*D[0] - M1[1]*D[1] - M1[2]*D[2];
//  f = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
  a = M0[0]*M0[0] + M0[1]*M0[1] + M0[2]*M0[2];
  b = - M1[0]*M0[0] - M1[1]*M0[1] - M1[2]*M0[2];
  c = M1[0]*M1[0] + M1[1]*M1[1] + M1[2]*M1[2];
  d = M0[0]*DD[0] + M0[1]*DD[1] + M0[2]*DD[2];
  e = - M1[0]*DD[0] - M1[1]*DD[1] - M1[2]*DD[2];
  f = DD[0]*DD[0] + DD[1]*DD[1] + DD[2]*DD[2];


  ddelta = a*c - b*b;

}

void hYr::nonparallel() {

  smin = (b*e - c*d);
  tmin = (b*d - a*e);

  if (smin >= 0){
    if (smin <= ddelta) {
      if (tmin >= 0) { 
	if(tmin <= ddelta) quad0(); 
	else quad3();
      }
      else {
	quad7();
      }
    }
    else {
      if (tmin >= 0 ){
	if (tmin <=ddelta) quad1();
	else quad2();
      }
      else {
	quad8();
      }
    }
  }
  else {
    if (tmin >= 0 ){
      if(tmin <=ddelta) quad5();
      else quad4();
    } 
    else {
      quad6();
    }
  }
}


void hYr::parallel() {

  if (b > 0){
    if (d >= 0) {smin=0; tmin=0;}
    else if (-d <= a) {smin=-d/a; tmin=0;}
    else { 
      smin = 1;
      if ( -(a+d) >= b) tmin = 1;
      else tmin = -(a+d)/b;
    }
  }
  else {
    if ( -d >= a) {smin=1; tmin=0;}
    else if ( d <= 0 ) {smin=-d/a; tmin=0;}
    else {
      smin = 0;
      if ( d >= -b ) tmin = 1;
      else tmin = -d/b;
    }
  }

  calcdist();
}

void hYr::calcdist(){

  /*distance = sqrt(a*smin*smin + 2*b*smin*tmin + c*tmin*tmin + 
    2*d*smin + 2*e*tmin + f);*/
  distance = a*smin*smin + 2*b*smin*tmin + c*tmin*tmin + 
    2*d*smin + 2*e*tmin + f;

}

void hYr::quad0() {

  double invdet = 1/ddelta;
  smin *= invdet;
  tmin *= invdet;

  calcdist();

}

void hYr::quad1() {

  smin = 1;
  tmin = -(b+e)/c;

  if(tmin < 0) tmin = 0;
  if(tmin > 1) tmin = 1;

  calcdist();

}

void hYr::quad2() {

  double dQdshalf = a+b+d;
  double dQdthalf = b+c+e;
  
  if (dQdshalf > 0) {
    
    tmin = 1;
    smin = -(b+d)/a;
    if (smin < 0) smin = 0;

  }
  else {
    smin = 1;
    if (dQdthalf > 0){
      tmin = -(b+e)/c;
      if (tmin < 0) tmin = 0;
    }
    else tmin = 1;
  }

  calcdist();

}

void hYr::quad3() {

  tmin = 1;
  smin = -(b+d)/a;

  if(smin < 0) smin = 0;
  if(smin > 1) smin = 1;

  calcdist();

}

void hYr::quad4() {

  double dQdshalf = b+d;
  double dQdthalf = c+e;
  
  if (dQdshalf < 0) {
    
    tmin = 1;
    smin = -(b+d)/a;
    if (smin > 1) smin = 1;

  }
  else {
    smin = 0;
    if (dQdthalf > 0){
      tmin = -e/c;
      if (tmin < 0) tmin = 0;
    }
    else tmin = 1;
  }

  calcdist();

}

void hYr::quad5() {

  smin = 0;
  tmin = -e/c;

  if(tmin < 0) tmin = 0;
  if(tmin > 1) tmin = 1;

  calcdist();

}

void hYr::quad6() {

  double dQdshalf = d;
  double dQdthalf = e;
  
  if (dQdshalf < 0) {
    
    tmin = 0;
    smin = -d/a;
    if (smin > 1) smin = 1;

  }
  else {
    smin = 0;
    if (dQdthalf < 0){
      tmin = -e/c;
      if (tmin > 1) tmin = 1;
    }
    else tmin = 0;
  }

  calcdist();

}

void hYr::quad7() {

  tmin = 0;
  smin = -d/a;

  if(smin < 0) smin = 0;
  if(smin > 1) smin = 1;

  calcdist();

}

void hYr::quad8() {

  double dQdshalf = a+d;
  double dQdthalf = b+e;
  
  if (dQdshalf > 0) {
    
    tmin = 0;
    smin = -d/a;
    if (smin < 0) smin = 0;

  }
  else {
    smin = 1;
    if (dQdthalf < 0){
      tmin = -(b+e)/c;
      if (tmin > 1) tmin = 1;
    }
    else tmin = 0;
  }

  calcdist();

}












