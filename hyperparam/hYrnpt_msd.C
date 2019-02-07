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
//extern double tem;

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
//double accum = 0.0;
double flip;
TLS = new hYrtls;
double bfm[3];


//{
flip=ran2(iran);
//cout<<"fl "<<pNr<<" "<<flip<<endl;
//}
		    if (flip < 1.0/3.0) {
		    for (int m=0; m<3; m++) {
			bfm[m] = rods[pNr].R[m];
		      }

		      TLS->accepted = particleMove(pNr);
		      if (TLS->accepted){
			  for (int m=0; m<3; m++) {
			  delta[pNr][m] += rods[pNr].R[m]-bfm[m] - box->x[m]*int((rods[pNr].R[m]-bfm[m])/box->halfx[m]);
			  }
			}
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

/*
		    if (flip < 1.0/3.0) TLS.accepted = TLS.particleMove(pNr);
		    else if (flip < 2.0/3.0) TLS.accepted = TLS.orientationalMove(pNr);
		    else TLS.accepted = TLS.axisrotMove(pNr);
		
		TLS.move=int(flip*3.0);

                if (TLS.move == 0) transMoves++;
                else if (TLS.move == 1) rotMoves++;
                else axrotMoves++;
                
                if (TLS.accepted){
                        if (TLS.move == 0) transAccept++;
                        else if (TLS.move == 1) rotAccept++;
                        else axrotAccept++;
                }
*/

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
//double accum = 0.0;

  double ***deltap = new double**[N_r];
  for (int i=0; i<N_r; i++){
    deltap[i] = new double*[10];
    for (int j=0; j<10; j++){
      deltap[i][j] = new double[3]();
    }
  }
  long double *msd1 = new long double[10]();
  long double *msd1xy = new long double[10]();
  long double *msd1z = new long double[10]();

  long double *msd1s = new long double[10]();
  long double *msd1xys = new long double[10]();
  long double *msd1zs = new long double[10]();

  double ***deltap10 = new double**[N_r];
  for (int i=0; i<N_r; i++){
    deltap10[i] = new double*[10];
    for (int j=0; j<10; j++){
      deltap10[i][j] = new double[3]();
    }
  }
  long double *msd10 = new long double[10]();
  long double *msd10xy = new long double[10]();
  long double *msd10z = new long double[10]();

  long double *msd10s = new long double[10]();
  long double *msd10xys = new long double[10]();
  long double *msd10zs = new long double[10]();

  double ****deltaps = new double***[5];
  for (int k=0; k<5; k++){
    deltaps[k] = new double**[N_r];
  for (int i=0; i<N_r; i++){
    deltaps[k][i] = new double*[10];
    for (int j=0; j<10; j++){
      deltaps[k][i][j] = new double[3]();
    }
  }
  }
  long double **msds = new long double*[5];
  for (int k=0; k<5; k++){
    msds[k] = new long double[10]();
  }
  long double **msdsxy = new long double*[5];
  for (int k=0; k<5; k++){
    msdsxy[k] = new long double[10]();
  }
  long double **msdsz = new long double*[5];
  for (int k=0; k<5; k++){
    msdsz[k] = new long double[10]();
  }

  hYrtls TLS;

strcpy(outdirt,outdir);

ofstream out(strcat(outdirt, "/Outfile"), ios::app);


	  out << "# StepsEQ: S2 [0 0 0] tAcc rAcc arAcc vAcc [0] pre" << endl;

neighborlist();
if( en == 2 ) nematicOPconfig();
for (int pNr=0; pNr<N_r; pNr++){
    if(Z != 0) ener[pNr]=calculateEnergy(pNr);
//cout<<pNr<<" "<<ener[pNr]<<endl;
    check[pNr]=0;
}
//////cout<<"h1"<<endl;
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

	  flip=ran2(iran);

	  int pNr=int((N_r)*ran2(iran));

	if (Nm == 1) {
/*for (int j=0; j<N_r; j++){
  cout<<"pbf "<<j<<" "<<rods[j].P<<" ";
}
cout<<endl;
*/
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
//cout<<"done own move"<<endl;
		m=0;
		for (int j=0; j<N_r; j++){
//		    ener[j]=calculateEnergy(j);
			check[j]=0;
		}
		}
	}

//out<<flip<<endl;
	}
	else
	 volumeMove();
       
     }

  for (int pNr=0; pNr<N_r; pNr++){
    for (int m=0; m<3; m++){
      msd1[0] += delta[pNr][m]*delta[pNr][m];
      for (int k=1; k<min(10,steps); k++){
        msd1[k] += (deltap[pNr][k-1][m]+delta[pNr][m])*(deltap[pNr][k-1][m]+delta[pNr][m]);
      }
//if(pNr==0) cout<<"1 "<<(deltap[pNr][8][m]+delta[pNr][m])*(deltap[pNr][8][m]+delta[pNr][m])<<endl;
    }
    for (int m=0; m<2; m++){
      msd1xy[0] += delta[pNr][m]*delta[pNr][m];
      for (int k=1; k<min(10,steps); k++){
        msd1xy[k] += (deltap[pNr][k-1][m]+delta[pNr][m])*(deltap[pNr][k-1][m]+delta[pNr][m]);
      }
    }
    for (int m=2; m<3; m++){
      msd1z[0] += delta[pNr][m]*delta[pNr][m];
      for (int k=1; k<min(10,steps); k++){
        msd1z[k] += (deltap[pNr][k-1][m]+delta[pNr][m])*(deltap[pNr][k-1][m]+delta[pNr][m]);
      }
    }
    for (int m=0; m<3; m++){
      for (int k=min(10,steps)-1; k>0; k--){
	deltap[pNr][k][m] = deltap[pNr][k-1][m]+delta[pNr][m];
      }
      deltap[pNr][0][m] = delta[pNr][m];
      delta[pNr][m] = 0.0;
    }
}
//cout<<endl;

if( steps%10 == 0 ){
for (int pNr=0; pNr<N_r; pNr++){
    for (int m=0; m<3; m++){
      msd10[0] += deltap[pNr][9][m]*deltap[pNr][9][m];
//if(pNr==0) cout<<"10 "<<deltap[pNr][9][m]*deltap[pNr][9][m]<<endl;
      for (int k=1; k<min(10,(steps)/10); k++){
        msd10[k] += (deltaps[0][pNr][k-1][m]+deltap[pNr][9][m])*(deltaps[0][pNr][k-1][m]+deltap[pNr][9][m]);
      }
    }
    for (int m=0; m<2; m++){
      msd10xy[0] += deltap[pNr][9][m]*deltap[pNr][9][m];
      for (int k=1; k<min(10,(steps)/10); k++){
        msd10xy[k] += (deltaps[0][pNr][k-1][m]+deltap[pNr][9][m])*(deltaps[0][pNr][k-1][m]+deltap[pNr][9][m]);
      }
    }
    for (int m=2; m<3; m++){
      msd10z[0] += deltap[pNr][9][m]*deltap[pNr][9][m];
      for (int k=1; k<min(10,(steps)/10); k++){
        msd10z[k] += (deltaps[0][pNr][k-1][m]+deltap[pNr][9][m])*(deltaps[0][pNr][k-1][m]+deltap[pNr][9][m]);
      }
    }
    for (int m=0; m<3; m++){
      for (int k=min(10,(steps)/10)-1; k>0; k--){
	deltaps[0][pNr][k][m] = deltaps[0][pNr][k-1][m]+deltap[pNr][9][m];
      }
      deltaps[0][pNr][0][m] = deltap[pNr][9][m];
    }
}
}

for (int ks=1; ks<5; ks++){
if( steps%int(pow(10,ks+1)) == 0 ){
for (int pNr=0; pNr<N_r; pNr++){
    for (int m=0; m<3; m++){
      msds[ks][0] += deltaps[ks-1][pNr][9][m]*deltaps[ks-1][pNr][9][m];
//if(pNr==0) cout<<"10 "<<deltap[pNr][9][m]*deltap[pNr][9][m]<<endl;
      for (int k=1; k<min(10,(steps)/int(pow(10,ks+1))); k++){
        msds[ks][k] += (deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m])*(deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m]);
      }
    }
    for (int m=0; m<2; m++){
      msdsxy[ks][0] += deltaps[ks-1][pNr][9][m]*deltaps[ks-1][pNr][9][m];
      for (int k=1; k<min(10,(steps)/int(pow(10,ks+1))); k++){
        msdsxy[ks][k] += (deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m])*(deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m]);
      }
    }
    for (int m=2; m<3; m++){
      msdsz[ks][0] += deltaps[ks-1][pNr][9][m]*deltaps[ks-1][pNr][9][m];
      for (int k=1; k<min(10,(steps)/int(pow(10,ks+1))); k++){
        msdsz[ks][k] += (deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m])*(deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m]);
      }
    }
    for (int m=0; m<3; m++){
      for (int k=min(10,(steps)/int(pow(10,ks+1)))-1; k>0; k--){
	deltaps[ks][pNr][k][m] = deltaps[ks][pNr][k-1][m]+deltaps[ks-1][pNr][9][m];
      }
      deltaps[ks][pNr][0][m] = deltaps[ks-1][pNr][9][m];
    }
}
}
}


      neighborlist();
      press += pressgrad;
//cout<<press<<endl;
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
	       << transRate << " " << rotRate << " " <<axrotRate << " " 
	       << " " <<volRate<< " 0" <<" "<<pre<< endl;
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
//  cout <<  maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax << endl;
    }
    
    if(snapshot == 0){
	snapshot = snapInterval;
	strcpy(outstring, "CoEQ");
	sprintf(temp,"%.*d", 0, int(steps/snapInterval));
	strcat(outstring, temp);
	writeOutConfig(outstring);

  char outdirt[200];
  char temp[5];
  ofstream msdfile;
  strcpy(outdirt,outdir);
	sprintf(temp,"/%.*d", 0, int(steps/snapInterval));
	strcat(outdirt, temp);
  msdfile.open (strcat(outdirt,"msd1.txt"), ios::out);

  for (int k=0; k<10; k++) {
    msd1s[k]=(1-float(snapInterval)/float(steps))*msd1s[k]+float(snapInterval)/float(steps)*msd1[k]/(snapInterval-k)/N_r;
    msd1xys[k]=(1-float(snapInterval)/float(steps))*msd1xys[k]+float(snapInterval)/float(steps)*msd1xy[k]/(snapInterval-k)/N_r;
    msd1zs[k]=(1-float(snapInterval)/float(steps))*msd1zs[k]+float(snapInterval)/float(steps)*msd1z[k]/(snapInterval-k)/N_r;
    msdfile << k+1 <<" "<< msd1s[k] <<" "<< msd1xys[k] <<" "<< msd1zs[k] <<endl;
//if(k==9) cout<<msd1[k]/(snapInterval-k)/N_r<<" "<<msd1s[k]<<endl;
    msd1[k]=0.0;
    msd1xy[k]=0.0;
    msd1z[k]=0.0;
  }
  for (int k=0; k<10; k++) {
    msd10s[k]=(1-float(snapInterval)/float(steps))*msd10s[k]+float(snapInterval)/float(steps)*msd10[k]/(snapInterval/10-k)/N_r;
    msd10xys[k]=(1-float(snapInterval)/float(steps))*msd10xys[k]+float(snapInterval)/float(steps)*msd10xy[k]/(snapInterval/10-k)/N_r;
    msd10zs[k]=(1-float(snapInterval)/float(steps))*msd10zs[k]+float(snapInterval)/float(steps)*msd10z[k]/(snapInterval/10-k)/N_r;
    msdfile << (k+1)*10 <<" "<< msd10s[k] <<" "<< msd10xys[k] <<" "<< msd10zs[k] <<endl;
//if(k==9 || k==0) cout<<msd10[k]/(snapInterval/10-k)/N_r<<" "<<msd10s[k]<<endl;
    msd10[k]=0.0;
    msd10xy[k]=0.0;
    msd10z[k]=0.0;
  }
  for (int ks=1; ks<5; ks++){
  for (int k=0; k<10; k++) {
    msdfile << (k+1)*pow(10,ks+1) <<" "<< msds[ks][k]/(steps/pow(10,ks+1)-k)/N_r <<" "<< msdsxy[ks][k]/(steps/pow(10,ks+1)-k)/N_r <<" "<< msdsz[ks][k]/(steps/pow(10,ks+1)-k)/N_r <<endl;
  }
  }

  msdfile.close();
	
    }
  }
  out << "# Equilibration ended maxDisplace maxRot maxARot dVmax: "<< maxDisplace <<" "<<maxRot<<" "<< maxARot<<" "<< dVmax << endl;

//     printf( "%lf\n", accum );
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

	  flip=ran2(iran);

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
//cout<<"done own move"<<endl;
		m=0;
		for (int j=0; j<N_r; j++){
//		    ener[j]=calculateEnergy(j);
			check[j]=0;
		}
		}
	}
//out<<flip<<endl;
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
//cout<<"ecalc "<<Ecalc<<endl;
         
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
	      << pre/float(measureCount);
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
//int tid;

//cout<<"tid "<<boost::this_thread::get_id()<<endl;

 
//  extern int N_r;
//  int nt=N_r;
//cout<<"N_r m "<<N_r<<endl;
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

    for (int i=0; i<3; i++){
	Rold[i] = rods[pNr].R[i];
//cout<<"rod "<<pNr<<" "<<rods[pNr].R[0]<<" "<<rods[pNr].R[1]<<" "<<rods[pNr].R[2]<<endl;
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

	r=ran2(iran);
	if (r<(exp(-(Enew-Eold)))){

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
//cout<<"yes"<<endl;
//cout<<"jap2 "<<pNr<<endl;
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
return acc;    


//out<<pNr<<endl;
    
}


bool hYr::orientationalMove(int pNr){
    double v[3];
    double newM[3], oldM[3],newV[3],newW[3],oldV[3],oldW[3];
    double newMAbsInv, newVabs, newWabs;
    bool overlap = false, acc = false;
    double *enerpair = new double[N_r]();

    randomVector(v);

//cout<<"h15"<<endl;
    double Enew, Eold, r, psi0=0.0, psi1;
 

//if(en == 2)  nematicOPconfig();
//cout<<"cOP "<<currOP<<endl;
//cout<<pNr<<"or"<<endl;    
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
//cout<<"Eold"<<Eold<<endl;
    for (int i=0;i<3;i++){
      oldM[i] = rods[pNr].M[i];
    }
if(en == 2) Eold += currOP*oldM[2]*oldM[2];

    for (int i=0; i<3; i++){
      newM[i] = maxRot*v[i]+ARInv*rods[pNr].M[i]; 
    }

    newMAbsInv = 1/sqrt(newM[0]*newM[0]+newM[1]*newM[1]
		       +newM[2]*newM[2]);


/*if (rods[pNr].M[1] != 0.0 || rods[pNr].M[2] != 0.0)	psi0 = acos(-rods[pNr].M[1]/sqrt(rods[pNr].M[1]*rods[pNr].M[1]+rods[pNr].M[2]*rods[pNr].M[2]));
    for (int i=0; i<3; i++){
      rods[pNr].M[i] = newM[i]*newMAbsInv*AR; 
    }
	psi1 = acos(-rods[pNr].M[1]/sqrt(rods[pNr].M[1]*rods[pNr].M[1]+rods[pNr].M[2]*rods[pNr].M[2]));
	rods[pNr].P = rods[pNr].P+psi1-psi0;
	{
    rods[pNr].clearPieces();
    insertRodIntoCells(rods[pNr]);*/

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
//cout<<psi1<<" ";//<<cos(psi0)<<" "<<newM[0]*newM[0]+newM[2]*newM[2]<<" "<<oldM[0]*oldM[0]+oldM[2]*oldM[2]<<" "<<-newM[2]/L*newV[0]+newM[0]/L*newV[2]<<endl;
//cout<<-newM[2]/L<<" "<<newV[0]<<" "<<newM[0]/L<<" "<<newV[2]<<endl;
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


//if (rods[pNr].M[1] != 0.0 || rods[pNr].M[2] != 0.0)	psi0 = acos(-rods[pNr].M[1]/sqrt(rods[pNr].M[1]*rods[pNr].M[1]+rods[pNr].M[2]*rods[pNr].M[2]));
//else psi0 = 0;
    for (int i=0; i<3; i++){
      rods[pNr].M[i] = newM[i]*newMAbsInv*AR; 
    }
//	psi1 = acos(-rods[pNr].M[1]/sqrt(rods[pNr].M[1]*rods[pNr].M[1]+rods[pNr].M[2]*rods[pNr].M[2]));
	rods[pNr].P = rods[pNr].P+psi1-psi0;
//cout<<pNr<<" "<<psi0<<" "<<psi1<<endl;
	{
    rods[pNr].clearPieces();
    insertRodIntoCells(rods[pNr]);
//cout<<endl;
//cout<<"onm "<<oldM[0]<<" "<<oldM[1]<<" "<<oldM[2]<<" "<<rods[pNr].M[0]<<" "<<rods[pNr].M[1]<<" "<<rods[pNr].M[2]<<endl;

       overlap = checkOverlapRodAll(rods[pNr]);
	}
  
    if (!overlap){

      if (Z!=0) {
	Enew = calculateEnergy(pNr);
      }
      else Enew = 0.0;
if(en == 2) Enew += currOP*newM[2]*newMAbsInv*AR*newM[2]*newMAbsInv*AR;

	r=ran2(iran);

	if (r<(exp(-(Enew-Eold)))){

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
//cout<<"yes"<<endl;
//cout<<"jap3 "<<pNr<<endl;
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

//cout<<"N_r m "<<N_r<<endl;


//out<<pNr<<"ax"<<endl;
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
//cout<<"Eo"<<Eold<<endl;

    u=ran2(iran)-0.5;
   

    oldPsi=rods[pNr].P;
    newPsi=fmod(oldPsi+u*maxARot,2.0*M_PI);

   rods[pNr].P=newPsi;

   if (Z!=0) {
	Enew = calculateEnergy(pNr);
   }
   else Enew = 0.0;
//cout<<"En"<<Enew<<endl;

 	r=ran2(iran);

	if (r<(exp(-(Enew-Eold)))){

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
//cout<<"psi "<<pNr<<" "<<rods[pNr].P<<" "<<oldPsi<<" "<<newPsi<<" "<<maxARot<<endl;
delete[] enerpair;
enerpair = NULL;
    return acc;
//out<<"acc mov"<<axrotAccept<<axrotMoves<<endl;
}

void hYr::volumeMove(){
  
    double oldVol, newVol, u;
    double Enew, Eold, r;
    double *eners = new double[N_r]();
   bool overlap = false;
	int dir;

	Eold = 0.0;
    for (int pNr=0; pNr<N_r; pNr++){
	Eold += ener[pNr];
//cout<<ener[pNr]<<endl;
    }
//    etest=calculateEnergy(401);
//cout<<"Eold "<<Eold<<" "<<etest<<endl;
    u=ran2(iran)-0.5;

        oldVol=box->V;
	newVol=oldVol+u*dVmax;
	u=ran2(iran);
//cout<<newVol<<" "<<dVmax<<endl;	
	if (u>0.5) dir=1;
	else dir=2;

if (en==2) compressBox(oldVol, newVol, dir);
else compressBoxV(oldVol, newVol, dir);
    volMoves++;

	
	for (int pNr=0; pNr<N_r; pNr++){

//for (int k=0; k<pNr; k++){
//overlap = checkOverlapRods(rods[pNr], rods[k]);
//if(overlap) cout<<pNr<<" "<<k<<endl;
//}

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
//cout<<"Enew "<<Enew<<endl;
//for (int pNr=0; pNr<N_r; pNr++){
//  if (fabs(ener[pNr]-eners[pNr])>1) cout<<pNr<<" "<<ener[pNr]<<" "<<eners[pNr]<<" "<<nbn[pNr]<<" "<<nlist[401*N_l]<<" "<<nlist[401*N_l+1]<<" "<<nlist[401*N_l+2]<<" "<<nlist[401*N_l+3]<<" "<<nlist[401*N_l+4]<<" "<<nlist[401*N_l+5]<<endl;
//}
//cout<<oldVol<<" "<<Eold<<" "<<newVol<<" "<<Enew<<endl;
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
