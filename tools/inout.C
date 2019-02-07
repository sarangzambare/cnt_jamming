#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include "lineDist.h"
#include "hYrnpt.h"

using namespace std;

void hYr::writeOutConfig(const char name[16]){

  char temp[16];
  char allname[200];



  strcpy(allname, outdir);
  strcat(allname, "/");
  strcat(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 5, float(N_r)*Vr/box->V);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 1, AR);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 2, D+1.0/ka);
  strcat(allname, temp);

//cout<<allname<<endl;
  ofstream outfile(allname, ios::out);
  outfile.precision(14);
//out<<"Nr"<<N_r<<endl;

  outfile << "#NumberR " << N_r << endl;
  outfile << "#L/D " << AR << endl;
  outfile << "#DensityR " << float(N_r)*Vr/box->V << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;

  for (int i=0; i<N_r; i++){
    outfile << rods[i].R[0] << " " << rods[i].R[1]  << " "
	    << rods[i].R[2] << " " << rods[i].M[0] << " " 
	    << rods[i].M[1] << " " << rods[i].M[2] << " " << rods[i].P << " " << locOP[i] << endl;
  }
  

  outfile.close();
}

void hYr::writeOutLargestCluster(const char name[16]){

  char temp[16];
  char allname[200];
  int largest = 0;
  Rod *curPart;

  strcpy(allname, outdir);
  strcat(allname, "/");
  strcat(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 5, float(N_r)*Vr/(box->V));
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 1, AR);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 2, D+1.0/ka);
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);

  outfile << "#NumberR " << N_r << endl;
  outfile << "#L/D " << AR << endl;
  outfile << "#DensityR " << float(N_r)*Vr/(box->V) << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;

  for (int i=0; i<N_r; i++){
      if (clusters[i].npart > clusters[largest].npart)
	  largest = i;
  }
  outfile << "# Size of this cluster: " << clusters[largest].npart << endl;
  curPart = clusters[largest].firstParticle;
  while (curPart){
      outfile << curPart->R[0] <<" "<< curPart->R[1] <<" "<< curPart->R[2] 
	      << " " 
	      << curPart->M[0] <<" "<< curPart->M[1] <<" "<< curPart->M[2] 
	      << endl;
      curPart = curPart->nextCl;
  }

  outfile.close();
}


void hYr::writeOutSimData(const char name[16]){

  char temp[16];
  char allname[200];

  strcpy(allname, outdir);
  strcat(allname, "/");
  strcat(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 5, float(N_r)*Vr/(box->V));
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 1, AR);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 2, D+1.0/ka);
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);

  outfile << "# L/D, N_r: " << AR << " " << N_r << endl;
  outfile << "# TransAkzeptanzrate: " 
	  << transRate/float(measureCount) << endl;
  outfile << "# RotAkzeptanzrate: " 
	  << rotRate/float(measureCount) << endl;
  outfile << "# AxRotAkzeptanzrate: "
            << axrotRate/float(measureCount) << endl;
  outfile << "# VolAkzeptanzrate: "
            << volRate/float(measureCount) << endl;
  outfile << "# x, y, z " << box->x[0] << " " << box->x[1] 
	  << " " << box->x[2] << endl;
  outfile << "# average OP: " << meanOP/float(measureCount) << endl;
  outfile << "# average fraction of rods in largest cluster: "
          << clustFrac/float(measureCount) << endl;
  outfile << "# average number of contacts per rod in whole system: "
          << allAveCont/float(measureCount) << endl;
  outfile << "# average number of contacts per rod in largest cluster: "
          << clustAveCont/float(measureCount) << endl;
  outfile << "# percolation probability: " << nPerc/float(measureCount) 
	  << endl;
  outfile << "# clustersize P(clustersize)" 
	  << endl;
  for (int i=1; i<=N_r; i++){
      if (meanClusterDist[i] > 0)
      outfile <<  i <<  " " << meanClusterDist[i]/float(measureCount) 
	      << endl;
  }
  outfile.close();
}


void hYr::writeOutDistances(const char name[16]){

  char temp[16];
  char allname[200];
  double shell;
  double ideal; // Number of idea gas particles in g(r) shell

  shell = (4.0*M_PI/3.0) * binWidth*binWidth*binWidth *
          N_r / (box->x[0] * box->x[1] * box->x[2]);

  strcpy(allname, outdir);
  strcat(allname, "/");
  strcat(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 5, float(N_r)*Vr/(box->V));
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 1, AR);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 2, D+1.0/ka);
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);

  outfile << "# L/D, N_r: " << AR << " " << N_r << endl;
  outfile << "# TransAkzeptanzrate: " 
	  << transRate/float(measureCount) << endl;
  outfile << "# RotAkzeptanzrate: " 
	  << rotRate/float(measureCount) << endl;
  outfile << "# AxRotAkzeptanzrate: "
            << axrotRate/float(measureCount) << endl;
  outfile << "# VolAkzeptanzrate: "
            << volRate/float(measureCount) << endl;
  outfile << "# x, y, z " << box->x[0] << " " << box->x[1] 
	  << " " << box->x[2] << endl;
  outfile << "# average S2: " << meanOP/float(measureCount) << endl;
  outfile << "# distance P(distance) g(r) g2(r)" << endl;
  for (int i=0; i<nbins; i++){
          ideal = ((i+1)*(i+1)*(i+1) - i*i*i) * shell;
	  outfile << (i+0.5)*binWidth << " " 
		  << meanDistances[i]/(N_r/2.0)/float(measureCount)/ideal << " "
                  << gofr[i]/(N_r/2.0)/float(measureCount)/ideal << " ";
          if (gofr[i] > 0) {
             outfile << g2ofr[i] / float(gofr[i]);
          } else {
             outfile << "0.0";
          }
          outfile << endl;
  }
  outfile.close();

  strcat(allname,"pp");
  ofstream outfile2(allname, ios::out);
  for (int i=0; i<nbins; i++){
    for (int j=0; j<nbins; j++){
      outfile2 << (i+0.5)*binWidth << " " << (j+0.5)*binWidth << " ";
      if (gofrpp[i][j] > 0) {
	outfile2 << g2ofrpp[i][j] / float(gofrpp[i][j]);
	} 
      else {
	outfile2 << "0.0";
      }
        outfile2 << endl;
    }
    outfile2 << endl;
  }
  outfile2.close();

}


// Write out histogram of contact positions, measured between 0 (end of rod) and 0.5 (middle).

void hYr::writeOutContacts(const char name[16]){

  char temp[16];
  char allname[200];
  double norm;

  strcpy(allname, outdir);
  strcat(allname, "/");
  strcat(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 5, float(N_r)*Vr/(box->V));
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 1, AR);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*f", 2, D+1.0/ka);
  strcat(allname, temp);

  norm = (double)N_r*measureCount/(2.0*sbins);

  ofstream outfile(allname, ios::out);

  outfile << "# End contacts per rod: " << (double)endCont/(N_r*measureCount) << endl;
  outfile << "# smin  frequency" << endl;
  for (int i=0; i<sbins; i++) {
     outfile << (i+0.5)/(2.0*sbins) << " " << sHist[i]/norm << endl;;
  }
  outfile.close();
}
