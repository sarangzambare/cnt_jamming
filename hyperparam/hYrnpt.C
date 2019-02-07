#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include "lineDist.h"
#include "hYrnpt.h"
#include <omp.h>

#define return_data_tag 2002
#define send_data_tag 2001

using namespace std;

const double Pi = M_PI;

extern float ran2(long *idum);

void hYr::run(){

  if (newSetUp)setUpNew();
  else setUpFromFile();
cout<<en<<endl;
cout<<pbc<<endl;
cout<<ka<<endl;
cout<<cof<<endl;
  writeOutConfig("Start");
  equil();
  MC();
  writeOutConfig("Config");
  measureClusterDist();
  writeOutLargestCluster("Cluster");
  writeOutSimData("Stats");
  writeOutDistances("Dist");
  writeOutContacts("Contacts");
  
}


void hYr::workerFunc(int pNr){
double flip;
TLS = new hYrtls;


flip=ran2(iran);
		    if (flip < 1.0/3.0) {

		      TLS->accepted = particleMove(pNr);
		    }
		    else if (flip < 2.0/3.0) {

TLS->accepted = orientationalMove(pNr);
}
		    else {

TLS->accepted = axisrotMove(pNr);
}
		
		TLS->move=int(flip*3.0);


{
                if (TLS->move == 0) transMoves++;
                else if (TLS->move == 1) rotMoves++;
                else axrotMoves++;
                
                if (TLS->accepted){
                        if (TLS->move == 0) transAccept++;
                        else if (TLS->move == 1) rotAccept++;
                        else axrotAccept++;
                }
}

delete TLS;
TLS = NULL;
    
}

void hYr::equil(){

  double flip, pre;
  int eqInterval;
  if (eqEvalSteps > 0 )
      eqInterval = eqSteps/eqEvalSteps;
  else eqInterval = 0;
  int measuring = eqInterval;
  int snapInterval;
  if (eqSnapSteps > 0)
      snapInterval = eqSteps/eqSnapSteps;
  else snapInterval = 0;
  int snapshot = snapInterval;
  char outstring[16];
  char temp[16];
  char outdirt[200];
  int *check = new int[N_r];

  hYrtls TLS;

strcpy(outdirt,outdir);

ofstream out(strcat(outdirt, "/Outfile"), ios::app);


	  out << "# StepsEQ: S2 [0 0 0] tAcc rAcc arAcc vAcc [0] pre" << endl;

neighborlist();
if( en == 2 ) nematicOPconfig();
for (int pNr=0; pNr<N_r; pNr++){
    if(Z != 0) ener[pNr]=calculateEnergy(pNr);
    check[pNr]=0;
}
int Nm = 1;
int *mlist = new int[Nm];

  for (int steps=1; steps<=eqSteps; steps++){
//cout<<steps<<" "<<endl;
    measuring--;
      snapshot--;

	int m=0;
      for (int i=0; i<N_r+en-1; i++){      

	flip=ran2(iran);

	if(flip>1.0/(N_r+1) || en==1){


	  int pNr=int((N_r)*ran2(iran));

	if (Nm == 1) {
	  workerFunc(pNr);
	}
	else {

	  if (check[pNr] == 0 && m <Nm) {
		mlist[m]=pNr;
		for (int j=0; j<nbn[pNr]; j++){
			check[nlist[pNr*N_l+j]] = 1;
		}
		check[pNr] = 1;
		m++;
	    }
	    else {
		
	      {
		
	
		for(int an_id = 0; an_id < m; an_id++){
		    int PNr = mlist[an_id];
		    workerFunc(PNr);
		}
		}
		m=0;
		for (int j=0; j<N_r; j++){
//		    ener[j]=calculateEnergy(j);
			check[j]=0;
		}
		}
	}

	}
	else
	 volumeMove();
       
     }

      neighborlist();
      press += pressgrad;
if( en == 2 ) nematicOPconfig();
    
      if (measuring==0){
	  nematicOPconfig();
	  transRate = float(transAccept)/float(transMoves);
	  rotRate = float(rotAccept)/float(rotMoves);
 	  axrotRate = float(axrotAccept)/float(axrotMoves);
  	  volRate = float(volAccept)/float(volMoves);
	calculateVirial();
	pre = Ecalc;
         
	  transAccept = 0; transMoves=0;
	  rotAccept = 0; rotMoves = 0;
	  axrotAccept = 0; axrotMoves = 0;
	  volAccept = 0; volMoves = 0;
	  
//	  out << "# StepsEQ: S2 [0 0 0] tAcc rAcc arAcc vAcc [0]" << endl;
	  out << steps << " " << currOP << " 0 0 0 " 
	       << transRate << " " << rotRate << " " << axrotRate << " " 
	       << " " << volRate << " 0" <<" "<< pre << " " << float(N_r)*Vr/box->V <<" "<<  maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax <<   endl;
	  measuring = eqInterval;

	  if (transRate > 0.55) maxDisplace *= transRate/0.5;
	  else if (transRate < 0.45) maxDisplace *= transRate/0.5;
	  if (rotRate > 0.55) maxRot *= rotRate/0.5;
	  else if (rotRate < 0.45) maxRot *= rotRate/0.5;
	  if (axrotRate > 0.55) maxARot *= min(axrotRate/0.5,2*M_PI/maxARot);
	  else if (axrotRate < 0.45) maxARot *= min(axrotRate/0.5,2*M_PI/maxARot);
	  if (volRate != 0.0) {
	  if (volRate > 0.55) dVmax *= volRate/0.5;
	  else if (volRate < 0.45) dVmax *= volRate/0.5;
	  }
	  else dVmax *= 0.5;
//  cout <<  maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax << endl;
    }
    
    if(snapshot == 0){
	snapshot = snapInterval;
	strcpy(outstring, "CoEQ");
	sprintf(temp,"%.*d", 0, int(steps/snapInterval));
	strcat(outstring, temp);
	writeOutConfig(outstring);
	
    }
  }
  out << "# Equilibration ended maxDisplace maxRot maxARot dVmax: "<< maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax << endl;

delete[] check;
check = NULL;
delete[] mlist;
mlist = NULL;
}


void hYr::MC(){

  double flip;
  char outstring[16];
  char temp[16];
  int mcInterval;
  if (mcEvalSteps > 0 )
      mcInterval = mcSteps/mcEvalSteps;
  else mcInterval = 0;
  int measuring = mcInterval;
  int snapInterval;
  if (mcSnapSteps > 0)
      snapInterval = mcSteps/mcSnapSteps;
  else snapInterval = 0;
  int snapshot = snapInterval;
  double pre;
  char outdirt[200];
  int *check = new int[N_r];
  hYrtls TLS;

  if (twist == 1) {
    twistprob = new double[100];
    for (int i=0; i<100; i++) {
      twistprob[i] = 0.0;
    }
  }
  
  strcpy(outdirt,outdir);

  ofstream out(strcat(outdirt, "/Outfile"), ios::app);


  transAccept = 0; transMoves=0;
  rotAccept = 0; rotMoves = 0;
  axrotAccept = 0; axrotMoves = 0;
  volAccept = 0; volMoves = 0;
  transRate = 0; rotRate = 0; axrotRate = 0; volRate = 0;
  pre = 0;
  meantwAngle = 0;
  eqPitch = 0;
  meanOP = 0;

out << "# StepsMC: S2 Lfrac AllCont ClustCont tAcc rAcc arAcc vAcc Perc pre" << endl;

neighborlist();
if( en == 2 ) nematicOPconfig();
for (int pNr=0; pNr<N_r; pNr++){
      if(Z != 0) ener[pNr]=calculateEnergy(pNr);
    check[pNr] = 0;
}


int Nm = 1;
int *mlist = new int[Nm];

  for (int steps=1; steps<=mcSteps; steps++){
//cout<<steps<<" "<<endl;
    measuring--;
      snapshot--;

	int m=0;
      for (int i=0; i<N_r+en-1; i++){      

	flip=ran2(iran);

	if(flip>1.0/(N_r+1) || en==1){


	  int pNr=int((N_r)*ran2(iran));

	if (Nm == 1) {
	  workerFunc(pNr);
	}
	else {

	  if (check[pNr] == 0 && m <Nm) {
		mlist[m]=pNr;
		for (int j=0; j<nbn[pNr]; j++){
			check[nlist[pNr*N_l+j]] = 1;
		}
		check[pNr] = 1;
		m++;
	    }
	    else {
		
	      {
		
	
		for(int an_id = 0; an_id < m; an_id++){
		    int PNr = mlist[an_id];
		    workerFunc(PNr);
		}
		}
		m=0;
		for (int j=0; j<N_r; j++){
//		    ener[j]=calculateEnergy(j);
			check[j]=0;
		}
		}
	}
	}
	else
	 volumeMove();
       
     }

     neighborlist();
      press += pressgrad;
if( en == 2 ) nematicOPconfig();
    
     if (measuring==0) {
	 measuring = mcInterval;

	 transRate += float(transAccept)/float(transMoves);
	 rotRate += float(rotAccept)/float(rotMoves);
 	 axrotRate += float(axrotAccept)/float(axrotMoves);
  	 volRate += float(volAccept)/float(volMoves);
	calculateVirial();
	pre += Ecalc;
         
	 transAccept = 0; transMoves=0;
	 rotAccept = 0; rotMoves = 0;
	 axrotAccept = 0; axrotMoves = 0;
	 volAccept = 0; volMoves = 0;

	 if (twist == 1){
	 twAngle = twistPitch();
	 meantwAngle += twAngle;
	 eqPitch = (2*M_PI/meantwAngle*float(measureCount))*box->x[2];
//	 eqPitch += (2*M_PI/twAngle)*box->x[2];
//	 cout<<"tw "<<twAngle<<" "<<eqPitch<<endl;
	 }

	 nematicOPconfig();
	 meanOP += currOP;
	 nPerc += percolate();

	 measureClusterDist();
	 for (int i=0; i<=N_r; i++){
	     meanClusterDist[i] += clusterDist[i];
	 }
         clustFrac += ((double)nlargest/N_r);
         allAveCont += ((double)allContacts/N_r);
         clustAveCont += ((double)clustContacts/nlargest);

	 measureDistances();
	 for (int i=0; i<nbins; i++){
	     meanDistances[i] += distances[i];
	 }
	 
	 measureCount++;
	 
//	 out << "# StepsMC: S2 Lfrac AllCont ClustCont tAcc rAcc arAcc vAcc Perc pre" << endl;
	 out << steps + eqSteps << " " << meanOP/float(measureCount) << " " 
              << (double)nlargest/N_r << " "
              << (double)allContacts/N_r << " "
              << (double)clustContacts/nlargest << " "
	      << transRate/float(measureCount) << " " 
	      << rotRate/float(measureCount) << " "
	      << axrotRate/float(measureCount) << " "
	      << volRate/float(measureCount) << " "
	      << nPerc/float(measureCount) << " "
	      << pre/float(measureCount) << " "
		  << float(N_r)*Vr/box->V ;
	      if (twist == 1) out << " " << eqPitch << " " << meantwAngle/float(measureCount);
//	      if (twist == 1) out << " " << eqPitch/float(measureCount) << " " << meantwAngle/float(measureCount);
	      out << endl;

//	  if (transRate/float(measureCount) > 0.55) maxDisplace *= transRate/float(measureCount)/0.5;
//	  else if (transRate/float(measureCount) < 0.45) maxDisplace /= transRate/float(measureCount)/0.5;
//	  if (rotRate/float(measureCount) > 0.55) maxRot *= rotRate/float(measureCount)/0.5;
//	  else if (rotRate/float(measureCount) < 0.45) maxRot /= rotRate/float(measureCount)/0.5;
//	  if (axrotRate/float(measureCount) > 0.55) maxARot *= axrotRate/float(measureCount)/0.5;
//	  else if (axrotRate/float(measureCount) < 0.45) maxARot /= axrotRate/float(measureCount)/0.5;
//	  if (volRate/float(measureCount) > 0.55) dVmax *= volRate/float(measureCount)/0.5;
//	  else if (volRate/float(measureCount) < 0.45) dVmax /= volRate/float(measureCount)/0.5;

     }
     
     if(snapshot == 0) {
	 snapshot = snapInterval;

	 locnematicOPconfig();

	 strcpy(outstring, "CoMC");
	 sprintf(temp,"%.*d", 0, int(steps/snapInterval));
	 strcat(outstring, temp);
	 writeOutConfig(outstring);
	 strcpy(outstring, "CluMC");
	 sprintf(temp,"%.*d", 0, int(steps/snapInterval));
	 strcat(outstring, temp);
	 writeOutLargestCluster(outstring);
  writeOutSimData("Stats");
  writeOutDistances("Dist");
  writeOutContacts("Contacts");
     }
  }

//  out << "# Simulation ended maxDisplace maxRot maxARot dVmax: "<< maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax << endl;
out.close();

/*
//cout<<outdirt<<endl;
ofstream out2;
  strcpy(outdirt,outdir);
  out2.open(strcat(outdirt,"/TwistProb"), ios::out);
for (int i=0; i<100; i++) {
//  cout<<i<<endl;
if(i != 50)  out2<<i*2*M_PI/maxtwi/100.0-M_PI/maxtwi<<" "<<log(twistprob[i]/float(measureCount))<<endl;
}
*/

delete[] check;
check = NULL;
delete[] mlist;
mlist = NULL;
}


bool hYr::particleMove(int pNr){
  

    bool overlap = false, acc = false;
    double Rold[3];
    double Enew, Eold, r;
    double *enerpair = new double[N_r]();
    double Eoldtest, Eoldtest2, Etest2, Etest ;
    double *eners = new double[N_r]();

   Eold = ener[pNr];
   if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
   for (int i=0; i<nbn[pNr]; i++){
    enerpair[nlist[pNr*N_l+i]] = calculateEnergyPair(nlist[pNr*N_l+i],pNr);
   }
}
   }


   Eoldtest = 0.0;
    for (int pNr2=0; pNr2<N_r; pNr2++){
	Eoldtest += ener[pNr2];
    }

    Etest = 0.0;
	if(Z !=0 ){
	for (int pNr2=0; pNr2<N_r; pNr2++){
	eners[pNr2] = calculateEnergy(pNr2);
	//else eners[pNr] = 0.0;
	Etest += eners[pNr2];
	}
	}

    for (int i=0; i<3; i++){
	Rold[i] = rods[pNr].R[i];
    }

    {
    displaceParticle(rods[pNr]);

    rods[pNr].clearPieces();
    insertRodIntoCells(rods[pNr]);


    overlap = checkOverlapRodAll(rods[pNr]);
    }
    if (!overlap){

      if (Z!=0){
	Enew = calculateEnergy(pNr);
	}
      else Enew = 0.0;




	if (Enew<Eold || ran2(iran)<(exp(-(Enew-Eold)))){

	acc = true;
	if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
	for (int i=0; i<nbn[pNr]; i++){
	    ener[nlist[pNr*N_l+i]] = ener[nlist[pNr*N_l+i]]-enerpair[nlist[pNr*N_l+i]]+calculateEnergyPair(nlist[pNr*N_l+i],pNr);
//	    ener[nlist[pNr*N_l+i]] = calculateEnergy(nlist[pNr*N_l+i]);
	}
}
	}
	ener[pNr] = Enew;
	

	Eoldtest2 = 0.0;
    for (int pNr2=0; pNr2<N_r; pNr2++){
	Eoldtest2 += ener[pNr2];
    }

    Etest2 = 0.0;
	if(Z !=0 ){
	for (int pNr2=0; pNr2<N_r; pNr2++){
	eners[pNr2] = calculateEnergy(pNr2);
	//else eners[pNr] = 0.0;
	Etest2 += eners[pNr2];
	}
	}

	cout<<"Etest: "<<Etest<<"Etest2: "<<Etest2<<"Eoldtest: "<<Eoldtest<<"Eoldtest2: "<<Eoldtest2<<endl;
    }
	else {
		for (int i=0; i<3; i++){
	    		rods[pNr].R[i] = Rold[i];
		}	

		{
		rods[pNr].clearPieces();     
		insertRodIntoCells(rods[pNr]);
		}

	}
    }
    
    else {
	for (int i=0; i<3; i++){
	    rods[pNr].R[i] = Rold[i];
	}	

	{
	rods[pNr].clearPieces();     
	insertRodIntoCells(rods[pNr]);
	}

    }
delete[] enerpair;
enerpair = NULL;
delete[] eners;
eners = NULL;
return acc;    

    
}


bool hYr::orientationalMove(int pNr){
    double v[3];
    double newM[3], oldM[3],newV[3],newW[3],oldV[3],oldW[3];
    double newMAbsInv, newVabs, newWabs;
    bool overlap = false, acc = false;
    double *enerpair = new double[N_r]();

    randomVector(v);

    double Enew, Eold, r, psi0=0.0, psi1;
 

//if(en == 2)  nematicOPconfig();
    Eold = ener[pNr];
    if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
   for (int i=0; i<nbn[pNr]; i++){
    enerpair[nlist[pNr*N_l+i]] = calculateEnergyPair(nlist[pNr*N_l+i],pNr);
   }
}
    }
    for (int i=0;i<3;i++){
      oldM[i] = rods[pNr].M[i];
    }
if(en == 2) Eold += currOP*oldM[2]*oldM[2];

    for (int i=0; i<3; i++){
      newM[i] = maxRot*v[i]+ARInv*rods[pNr].M[i]; 
    }

    newMAbsInv = 1/sqrt(newM[0]*newM[0]+newM[1]*newM[1]
		       +newM[2]*newM[2]);


psi0=rods[pNr].P;
//v'.v~
double om=oldM[0]*oldM[0]+oldM[2]*oldM[2];
double so=1.0/sqrt(om);
double sn=1.0/sqrt(newM[0]*newM[0]+newM[2]*newM[2]);
oldV[0]=-oldM[2]*so;
oldV[1]=0.0;
oldV[2]=oldM[0]*so;
    for (int i=0; i<3; i++){
      newV[i] = maxRot*v[i]+oldV[i];
    }
    newVabs = sqrt(newV[0]*newV[0]+newV[1]*newV[1]+newV[2]*newV[2]);
    for (int i=0; i<3; i++){
      newV[i] /= newVabs;
    } 

psi1=cos(psi0)*(-newM[2]*newV[0]+newM[0]*newV[2])*sn;

//v'.w~
oldW[0]=oldM[0]*oldM[1]*Linv*so;
oldW[1]=-(om)*Linv*so;
oldW[2]=oldM[1]*oldM[2]*Linv*so;
    for (int i=0; i<3; i++){
      newW[i] = maxRot*v[i]+oldW[i];
    }
    newWabs = sqrt(newW[0]*newW[0]+newW[1]*newW[1]+newW[2]*newW[2]);
    for (int i=0; i<3; i++){
      newW[i] /= newWabs;
    } 

psi1+=sin(psi0)*(-newM[2]*newW[0]+newM[0]*newW[2])*sn;
//cout<<psi1<<endl;
if(fabs(psi1)>1.0) psi1=copysign(1,psi1);
if(fabs(psi0)<M_PI) psi1=copysign(acos(psi1),psi0);
else psi1=copysign(2.0*M_PI-acos(psi1),psi0);


    for (int i=0; i<3; i++){
      rods[pNr].M[i] = newM[i]*newMAbsInv*AR; 
    }
	rods[pNr].P = rods[pNr].P+psi1-psi0;
	{
    rods[pNr].clearPieces();
    insertRodIntoCells(rods[pNr]);

       overlap = checkOverlapRodAll(rods[pNr]);
	}
  
    if (!overlap){

      if (Z!=0) {
	Enew = calculateEnergy(pNr);
      }
      else Enew = 0.0;
if(en == 2) Enew += currOP*newM[2]*newMAbsInv*AR*newM[2]*newMAbsInv*AR;


	if (Enew<Eold || ran2(iran)<(exp(-(Enew-Eold)))){

      acc = true;
      if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
        for (int i=0; i<nbn[pNr]; i++){
	    ener[nlist[pNr*N_l+i]] = ener[nlist[pNr*N_l+i]]-enerpair[nlist[pNr*N_l+i]]+calculateEnergyPair(nlist[pNr*N_l+i],pNr);
            }
}
      }
	if(en == 2) Enew -= currOP*newM[2]*newMAbsInv*AR*newM[2]*newMAbsInv*AR;
	ener[pNr] = Enew;
	}
    
	else {
		for (int i=0; i<3; i++){
	    		rods[pNr].M[i] = oldM[i];
		}	
	rods[pNr].P = rods[pNr].P+psi0-psi1;

		rods[pNr].clearPieces();     


		insertRodIntoCells(rods[pNr]);

	}
    }
    else {
	for (int i=0; i<3; i++){
	    rods[pNr].M[i] = oldM[i];
	}
	rods[pNr].P = rods[pNr].P+psi0-psi1;

	{
	rods[pNr].clearPieces();
	insertRodIntoCells(rods[pNr]);
	}

    }

delete[] enerpair;
enerpair = NULL;
    return acc;
 }


bool hYr::axisrotMove(int pNr){
  
    double oldPsi, newPsi, u;
    double Enew, Eold, r;
	bool acc = false;
    double *enerpair = new double[N_r]();

    Eold = ener[pNr];
    if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
   for (int i=0; i<nbn[pNr]; i++){
    enerpair[nlist[pNr*N_l+i]] = calculateEnergyPair(nlist[pNr*N_l+i],pNr);
   }
}
    }

    u=ran2(iran)-0.5;
   

    oldPsi=rods[pNr].P;
    newPsi=fmod(oldPsi+u*maxARot,2.0*M_PI);

   rods[pNr].P=newPsi;

   if (Z!=0) {
	Enew = calculateEnergy(pNr);
   }
   else Enew = 0.0;


	if (Enew<Eold || ran2(iran)<(exp(-(Enew-Eold)))){

      acc = true;

	if(Z != 0){
#pragma omp parallel num_threads(2)
{
#pragma omp for
        for (int i=0; i<nbn[pNr]; i++){
	    ener[nlist[pNr*N_l+i]] = ener[nlist[pNr*N_l+i]]-enerpair[nlist[pNr*N_l+i]]+calculateEnergyPair(nlist[pNr*N_l+i],pNr);
              }
}
	}
	ener[pNr] = Enew;
	}
    
	else {
		for (int i=0; i<3; i++){
	    		rods[pNr].P = oldPsi;
		}	
   	}
delete[] enerpair;
enerpair = NULL;
    return acc;
}

void hYr::volumeMove(){
  
  	double Etest;
    double oldVol, newVol, u;
    double Enew, Eold, r;
    double *eners = new double[N_r]();
   bool overlap = false;
	int dir;

	Eold = 0.0;
    for (int pNr=0; pNr<N_r; pNr++){
	Eold += ener[pNr];
    }

    Etest = 0.0;
	if(Z !=0 ){
	for (int pNr=0; pNr<N_r; pNr++){
	eners[pNr] = calculateEnergy(pNr);
	//else eners[pNr] = 0.0;
	Etest += eners[pNr];
	}
	}

    u=ran2(iran)-0.5;

        oldVol=box->V;
	newVol=oldVol+u*dVmax;
	u=ran2(iran);
	if (u>0.5) dir=1;
	else dir=2;

if (en==2) compressBox(oldVol, newVol, dir);
else compressBoxV(oldVol, newVol, dir);
    volMoves++;

	
	for (int pNr=0; pNr<N_r; pNr++){

    overlap = checkOverlapRodAll(rods[pNr]);
  
	if (overlap){
if (en==2) compressBox(newVol, oldVol, dir);
else compressBoxV(newVol, oldVol, dir);
delete[] eners;
eners = NULL;
	return;
    	}

	}

//    if (!overlap){
	Enew = 0.0;
	if(Z !=0 ){
	for (int pNr=0; pNr<N_r; pNr++){
	eners[pNr] = calculateEnergy(pNr);
	//else eners[pNr] = 0.0;
	Enew += eners[pNr];
	}
	}
	//cout<<"Enew: "<<Enew<<" Eold: "<<Eold<<"Etest: "<<Etest<<endl;

	r=ran2(iran);

	if (en==2 && r<(exp(-(Enew-Eold)/2.0-press*(newVol-oldVol)/tempr+N_r*log(newVol/oldVol)))){

      volAccept ++;
      for (int pNr=0; pNr<N_r; pNr++){
	ener[pNr] = eners[pNr];
      }
	}
	else if (en==3 &&  r<(exp(-(Enew-Eold)/2.0))){

	volAccept ++;
	for (int pNr=0; pNr<N_r; pNr++){
	    ener[pNr] = eners[pNr];
	} 

	}
	else {
if (en==2) compressBox(newVol, oldVol, dir);
else compressBoxV(newVol, oldVol, dir);
		
    	}
//    }
//    else {
//	compressBox(newVol, oldVol);
//	rods[pNr].clearPieces();
//	insertRodIntoCells(rods[pNr]);
//    }

delete[] eners;
eners = NULL;


}

double hYr::twistPitch(){

  double twPhi, twEner, minPhi, minEner, r, Ener0;
//  neighborlistTwist();

  minPhi=0;
  minEner = 0.0;
  for (int pNr=0; pNr<N_r; pNr++) {
  minEner += calculateEnergyTwist(pNr,0.0);
  }
  Ener0 = minEner;
//  cout<<"minE "<<minEner<<endl;
  for (int i=0; i<100; i++){
    twPhi=-M_PI/maxtwi + i*2*M_PI/maxtwi/100.0;
    twEner = twistBox(twPhi);

    r = ran2(iran);
    if (r < (exp(-(twEner-Ener0)/tempr)) && i !=50 ) twistprob[i] += 1.0;

    if (twEner < minEner) {
	minEner = twEner;
	minPhi = twPhi;
	}
//  cout<<"minE "<<minPhi<<" "<<minEner<<" "<<twPhi<<" "<<twEner<<endl;
  }
//cout<<endl;
  return minPhi;

}


int main(int argc, char *argv[]){

   hYr HYr;
   hYrtls htls;

  if (argc < 3) {
    cerr << "Too few arguments \n";
    HYr.usage();
  }

  while ((argc > 1) && (argv[1][0] == '-')) {

    switch (argv[1][1]) {
    case 'n':              // new Setup
      HYr.newSetUp=true;
      break;
    case 'f':
      HYr.newSetUp=false;
      HYr.ConfigFile=&argv[1][2];
      break;
    case 'p' :
      HYr.ParameterFile=&argv[1][2];
      break;
    default:
      cerr << "Bad option " << argv[1] << '\n';
      HYr.usage();
    }

    ++argv;
    --argc;
  }

//  out << "# Parameterfile: " << HYr.ParameterFile << ", new SetUp: " <<
//    HYr.newSetUp << endl;



  HYr.run();





}
