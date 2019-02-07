#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include "lineDist.h"
#include "hYrnpt.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>

#define send_data_tag 2001
#define return_data_tag 2002

using namespace std;


extern float ran2(long *idum);

// reads simulation parameters in case option -n ("set up new
// simulation") is chosen and sets up a start configuration

void hYr::setUpNew(){

    ifstream infile;
    char line[200];
    char quantity[200];
    float number;
    int check;
    char outdirt[200];
    char temp[200];
    box = new Box;

    cof=2.0;

    infile.open(ParameterFile);
    if (infile.bad()) {
	cerr << "Error " << ParameterFile << " not found\n";
	exit(8);
    }
    while (infile.peek() != EOF) {
	infile.getline(line,sizeof(line));
	check = sscanf(line, "%s%f", quantity, &number);
	if (check != 2) {
	    cerr << "Error: wrong line format in " << ParameterFile << endl;
	    exit(8);
	}
	
	if (strcmp(quantity,"NumberRods")==0) {
	    Nin_r=int(number);
	}
	else if (strcmp(quantity,"Boxx")==0) box->x[0]=box->x[1]=number;
	else if (strcmp(quantity,"Boxz")==0) box->x[2]=number;
	else if (strcmp(quantity,"CutoffFactor")==0) cof=number;
	else if (strcmp(quantity,"L/D")==0) AR=number;
	else if (strcmp(quantity,"Seed")==0) seed=int(number); 
	else if (strcmp(quantity,"McSteps")==0) mcSteps=int(number);
	else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
	else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
	else if (strcmp(quantity,"EqSteps")==0) eqSteps=int(number); 
	else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
	else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);  
	else if (strcmp(quantity,"MaxDisplacement")==0) maxDisplace=number;
	else if (strcmp(quantity,"MaxRotation")==0) maxRot=number;
	else if (strcmp(quantity,"MaxAxRotation")==0) maxARot=number;
	else if (strcmp(quantity,"MaxVolumeCh")==0) dVmax=number;
	else if (strcmp(quantity,"Distance")==0) A=number;
	else if (strcmp(quantity,"NInt")==0) Nint=int(number);
	else if (strcmp(quantity,"Temperature")==0) tempr=number;
	else if (strcmp(quantity,"Charges")==0) Z=int(number); 
	else if (strcmp(quantity,"BjLength")==0) lam=number;
	else if (strcmp(quantity,"ScreenConst")==0) ka=number;  
	else if (strcmp(quantity,"Diameter")==0) D=number;
	else if (strcmp(quantity,"PitchFac")==0) pf=number;
	else if (strcmp(quantity,"kWaveCalc")==0) kw=int(number);
	else if (strcmp(quantity,"Ensemble")==0) en=int(number);  
	else if (strcmp(quantity,"Pressure")==0) press=number;
	else if (strcmp(quantity,"PBC")==0) pbc=int(number);  
 	else if (strcmp(quantity,"torsionAngle")==0) phi=number;
	else if (strcmp(quantity,"SaltConc")==0) cs=number;
	else if (strcmp(quantity,"Twist")==0) twist=int(number);
	else if (strcmp(quantity,"Pi/MaxTwist")==0) maxtwi=number;

    }
    
    infile.close();

//    cof=max(2,int(2.0*D*sqrt(4.0*M_PI*lam*(Z*(4*Nin_r*Nin_r)/(box->x[0]*box->x[1]*box->x[2])+2.0*cs)))+1);
    cof=pow((box->x[0]*box->x[1]*box->x[2])/(4*Nin_r*Nin_r)/M_PI/AR*4.0,1.0/2.0)*sqrt(4.0*M_PI*lam*(Z*(4*Nin_r*Nin_r)/(box->x[0]*box->x[1]*box->x[2])+2.0*cs));
//    if (en==1 && pbc==1) strcpy(outdirt, "NVT_pbc/");
//    else if (en==1 && pbc==2) strcpy(outdirt, "NVT_tpbc/");
//    else if (en==2 && pbc==1) strcpy(outdirt, "NpT_pbc/");
//    else if (en==2 && pbc==2) strcpy(outdirt, "NpT_tpbc/");
    strcpy(outdirt, "");
    strcat(outdirt, "N");
    sprintf(temp, "%i", 4*Nin_r*Nin_r);
    strcat(outdirt, temp);
    if (en==1 || en==3) {
	sprintf(temp, "%.*f",0, box->x[0]*box->x[1]*box->x[2]);
	strcat(outdirt, "_V");
	strcat(outdirt, temp);
    }
    else {
	sprintf(temp, "%.*f",4, press);
	strcat(outdirt, "_p");
	strcat(outdirt, temp);
    }
    strcat(outdirt, "_AR");
    sprintf(temp, "%.*f", 0, AR);
    strcat(outdirt, temp);
    strcat(outdirt, "_T");
    sprintf(temp, "%.*f", 1, tempr);
    strcat(outdirt, temp);
    strcat(outdirt, "_Z");
    sprintf(temp, "%.i", Z);
    strcat(outdirt, temp);
    strcat(outdirt, "_npc");
    sprintf(temp, "%.i", Nint+1);
    strcat(outdirt, temp);
    strcat(outdirt, "_l");
    sprintf(temp, "%.*f", 2, lam);
    strcat(outdirt, temp);
    strcat(outdirt, "_cs");
    sprintf(temp, "%.*f", 3, cs);
    strcat(outdirt, temp);
    strcat(outdirt, "_pf");
    sprintf(temp, "%.*f", 1, pf);
    strcat(outdirt, temp);
    if (pbc==2){
     	strcat(outdirt, "_phi");
    	sprintf(temp, "%.*f", 1, phi);
    	strcat(outdirt, temp);
    }

	struct stat st = {0};

strcpy(outdir,outdirt);
	if (stat(outdir, &st) == -1) {
   	 mkdir(outdir, 0700);
	}
 

    ofstream out(strcat(outdirt, "/Outfile"));

    if (maxDisplace >= 1){
      cerr << "Maxdisplace must be smaller than 1" << endl;
      exit(8);
    }
    
    N_r=Nin_r*Nin_r*4;
    
    rods = new Rod[N_r];
    
    out << "# Rods: " << N_r <<  endl;
    out << "# AR: " << AR << " Box " << box->x[0] << "x" << box->x[1]
	 << "x" << box->x[2] << endl; 
    out << "# Eq " << eqSteps << " Mc " << mcSteps 
	 << " EqEval " << eqEvalSteps << " McEval " << mcEvalSteps
	 << " EqSnaps " << eqSnapSteps << " Mcsnaps " << mcSnapSteps << endl;
    
    iran = &seed;
    
    packAAAExpandAndShift();
    
    initialize();
}


void hYr::packAAAExpandAndShift(){

  double zpos, ypos, xpos, dy, dz;
  int l;

  box->update();
  
  if (N_r > 0){
    if (box->x[0] < 2*(AR+1)){
      cerr << "Density too high, rods do not fit into box" << endl;
      exit(8);
    }
  }

  dy = box->x[1]/float(Nin_r);
  dz = box->x[2]/float(2*Nin_r);
  if (dy < 1 || dz < 1){
    cerr << "Nin_rods too high!" << endl;
    exit(8);
  }
  
  zpos = -0.5*box->x[2]+0.5;// + 0.5*(AR+1);
  ypos = -0.5*box->x[1] + 0.5;
  xpos = -0.5*box->halfx[0];
  
  for (int j=0; j<2*Nin_r; j++){
    for (int k=0; k<Nin_r; k++){
      
      l=Nin_r*j+k;	  
      rods[l].R[0] = xpos;
      rods[l].R[1] = ypos + k*dy;
      rods[l].R[2] = zpos;
      rods[l].M[0] = AR;
      rods[l].M[1] = 0;
      rods[l].M[2] = 0;
      rods[l].P=0;
    }
    zpos += dz;
  }
    
  zpos = -0.5*box->x[2]+0.5;// - 0.5*(AR+1);
  ypos = -0.5*box->x[1] + 0.5;
  xpos = 0.5*box->halfx[0];
  
  for (int j=0; j<2*Nin_r; j++){
    for (int k=0; k<Nin_r; k++){
      
      l=2*Nin_r*Nin_r+Nin_r*j+k;	  
      rods[l].R[0] = xpos;
      rods[l].R[1] = ypos + k*dy;
      rods[l].R[2] = zpos;
      rods[l].M[0] = AR;
      rods[l].M[1] = 0;
      rods[l].M[2] = 0;
      rods[l].P=0;
    }
    zpos += dz;
  }
  

  double shift;
//  double save;

  for (int j=0; j<2*Nin_r; j++){
    for (int k=0; k<Nin_r; k++){	
      l=Nin_r*j+k;
     
//      shift = ran2(iran)*box->halfx[2];
      shift = 0.5*ran2(iran)*(box->halfx[0]-(AR+1));

/*ch = int((rods[l].R[2]+shift)*box->halfxInv[2]);

if(pbc==2 && ch != 0){

// phi=pi/2:
	save = rods[l].R[0];
	rods[l].R[0] = int(2*phi)*ch*rods[l].R[1];
	rods[l].R[1] = -int(2*phi)*ch*save;

	save = rods[l].M[0];
	psi0 = acos(-rods[l].M[1]/sqrt(rods[l].M[1]*rods[l].M[1]+rods[l].M[2]*rods[l].M[2]));
	rods[l].M[0] = int(2*phi)*ch*rods[l].M[1];
	rods[l].M[1] = -int(2*phi)*ch*save;

	psi1 = acos(-rods[l].M[1]/sqrt(rods[l].M[1]*rods[l].M[1]+rods[l].M[2]*rods[l].M[2]));
	rods[l].P = rods[l].P+psi1-psi0;


}
*/

      rods[l].R[0] = -int((rods[l].R[0]+shift)*box->halfxInv[0])*box->halfx[0] 
      + fmod((rods[l].R[0]+shift),box->halfx[0]);

/*ch = int((rods[l+Nin_r*Nin_r].R[2]+shift)*box->halfxInv[2]);

if(pbc==2 && ch != 0){

// phi=pi/2:
	save = rods[l+Nin_r*Nin_r].R[0];
	rods[l+Nin_r*Nin_r].R[0] = int(2*phi)*ch*rods[l+Nin_r*Nin_r].R[1];
	rods[l+Nin_r*Nin_r].R[1] = -int(2*phi)*ch*save;

	save = rods[l+Nin_r*Nin_r].M[0];
	psi0 = acos(-rods[l+Nin_r*Nin_r].M[1]/sqrt(rods[l+Nin_r*Nin_r].M[1]*rods[l+Nin_r*Nin_r].M[1]+rods[l+Nin_r*Nin_r].M[2]*rods[l+Nin_r*Nin_r].M[2]));
	rods[l+Nin_r*Nin_r].M[0] = int(2*phi)*ch*rods[l+Nin_r*Nin_r].M[1];
	rods[l+Nin_r*Nin_r].M[1] = -int(2*phi)*ch*save;

	psi1 = acos(-rods[l+Nin_r*Nin_r].M[1]/sqrt(rods[l+Nin_r*Nin_r].M[1]*rods[l+Nin_r*Nin_r].M[1]+rods[l+Nin_r*Nin_r].M[2]*rods[l+Nin_r*Nin_r].M[2]));
	rods[l+Nin_r*Nin_r].P = rods[l+Nin_r*Nin_r].P+psi1-psi0;


}
*/
      rods[l+2*Nin_r*Nin_r].R[0] = -int((rods[l+ 2*Nin_r*Nin_r].R[0]+shift)*box->halfxInv[0])*box->halfx[0] 
      + fmod((rods[l+ 2*Nin_r*Nin_r].R[0]+shift),box->halfx[0]);
    
    }
  }

  return;
}


void hYr::initialize(){

    double lX;
//    int ierr;
    int my_id;
    //set counters
    
    transAccept = 0; transMoves=0;
    rotAccept = 0; rotMoves = 0;
    axrotAccept = 0; axrotMoves = 0;
    volAccept = 0; volMoves = 0;
    transRate = 0; rotRate = 0; axrotRate = 0; volRate = 0;
    measureCount = 0;
    
//    rhocp_r = 2/(sqrt(2.0) + AR*sqrt(3.0));

    Vr = M_PI*(0.25*AR + 1.0/6.0);

    AR2 =AR*AR;
    AR12 = (AR+1)*(AR+1);
    ARInv = 1.0/AR;
    hAR = 0.5*AR;
    A2 = A*A;
    ARA2 = (AR+A)*(AR+A);
    A12 = (A+1)*(A+1);
    
    Nrinv = 1/float(N_r);
    NN1hInv = 2.0/float(N_r)/float(N_r+1);

    //parameters

//    tempr=1.0;
//    Z=1;
//    lam=0.1;
//    ka=5.0;
//    D=1.0;
    L=AR*D;
    Linv=1.0/L;
    Dhalf=0.5*D;

    
//    pitch=pf*L;
    pitch=pf*D;
if(kw == 1) kwave=M_PI*L/pitch;
else kwave=0.0;
    tkwave=kwave*2.0;
    ka=sqrt(4.0*M_PI*lam*(Z*N_r/box->V+2.0*cs));
    cofka=cof/ka;
    dpcof=D+cofka;
    cofka2=cofka*cofka;
    dpcof2=dpcof*dpcof;
    shift=exp(-cof)/cofka;

    i2phi = int(2*phi);

    //Nint=int(16*L/pitch);
    NintInv=1.0/Nint;
    Zpc=Z/(Nint+1);
    tzl=0.25*Zpc*Zpc*lam;
    nt=N_r*tempr;
    
    // data arrays
    
    meanClusterDist = new int[N_r+1];
    clusterDist = new int[N_r+1];
    for (int i=0; i<=N_r; i++){
	meanClusterDist[i] = 0;
    }
    

    N_l = int(N_r*(AR+1+2.0/ka)*(AR+1+2.0/ka)*(AR+1+2.0/ka)/box->V);
    N_l = N_r;
      nlist = new int[N_r*N_l];

  for(int i = 0; i < N_r; i++){
  for(int j = 0; j < N_l; j++){
nlist[i*N_l+j]=0;
  }
  }

  nbn = new int[N_r+1];
    for (int i=0; i<N_r; i++){
	    nbn[i]=0;
		}
    
    ener = new double[N_r];
         for (int i=0; i<N_r; i++){
	    ener[i]=0;
	    }

    locOP = new int[N_r];
         for (int i=0; i<N_r; i++){
	    locOP[i]=0;
	    }

    delta = new double*[N_r];
    for (int i=0; i<N_r; i++){
      delta[i] = new double[3]();
    }

    clusters = new Cluster[N_r];
    clearClusters(); 
    
    distances = new int[nbins+1];
    meanDistances = new int[nbins+1];
    gofr = new int[nbins+1];
    g2ofr = new double[nbins+1];

    gofrpp = new int*[nbins+1];
    g2ofrpp = new double*[nbins+1];

    sHist = new long[sbins+1];
    lX = 0.5*box->x[0];
    if (0.5*box->x[1] > lX) lX = 0.5*box->x[1];
    if (0.5*box->x[2] > lX) lX = 0.5*box->x[2];
    binWidth = lX/(nbins-1);
    for (int i=0; i<nbins+1; i++){
	distances[i] = 0;
	meanDistances[i] = 0;
        gofr[i] = 0;
        g2ofr[i] = 0.0;
	gofrpp[i] = new int[nbins+1]();
	g2ofrpp[i] = new double[nbins+1]();
    }
    for (int i=0; i<=sbins; i++){
        sHist[i] = 0;
    }

    meanOP = 0; 
    meanV = 0; 
    nPerc = 0;
    clustFrac = 0.0;
    allAveCont = clustAveCont = 0.0;
    endCont = 0;

    setupList();

    for (int i=0; i<N_r; i++){
	rods[i].neighbours = (Rod **)malloc(N_r*sizeof(Rod *));
    }
						 
    out << "# Using " << box->nxCell[0] << " x "
	 << box->nxCell[1] << " x " << box->nxCell[2] << " = "
	 << box->nCells << " cells " << endl;
    
    bool overlap = false;
    double d[3], dm[3];
    double distance;
    double s0, t0; 

    double save;
    int ch;
    Rod p1, p2;
    double p1R[3], p1M[3];

    for (int i=0; i<N_r; i++){
      for (int j=0; j<i; j++){
	overlap = checkOverlapRods(rods[i], rods[j]);
	if (overlap){
	  cerr << "Overlapping start configuration!" <<my_id<<" i "<<i<<" j "<<j<< endl;
	  cerr <<rods[j].R[0]<<" "<<rods[j].R[1]<<" "<<rods[j].R[2]<<" " 
	       <<rods[j].M[0]<<" "<<rods[j].M[1]<<" "<<rods[j].M[2]
	       << endl;
	  cerr <<rods[i].R[0]<<" "<<rods[i].R[1]<<" "<<rods[i].R[2]<<" " 
	       <<rods[i].M[0]<<" "<<rods[i].M[1]<<" "<<rods[i].M[2]
	       << endl;


		for (int k=0; k<3; k++){
		p1.R[k] = rods[i].R[k];
		p2.R[k] = rods[j].R[k];
		p1.M[k] = rods[i].M[k];
		p2.M[k] = rods[j].M[k];
		}

ch = int((p1.R[2] - p2.R[2])*box->halfxInv[2]);

if (pbc==2 && ch !=0){

// phi=pi/2:
	p1R[0] = int(2*phi)*ch*p1.R[1];
	p1R[1] = -int(2*phi)*ch*p1.R[0];
	p1R[2] = p1.R[2];

	save = p1.M[0];
	p1M[0] = int(2*phi)*ch*p1.M[1];
	p1M[1] = -int(2*phi)*ch*save;
	p1M[2] = p1.M[2];

}
else{
	for (int i =0; i<3; i++){
		p1R[i] = p1.R[i];
		p1M[i] = p1.M[i];
	}
}	





	  for (int k=0; k<3; k++){
	    d[k] = rods[j].R[k] - p1R[k];
	    if (fabs(d[k]) > box->halfx[k]) 
	      d[k] -= copysign(box->x[k], d[k]);
	  } 
	  
	  if ((d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) < AR12){
	    
	    for (int k=0; k<3; k++){
	      dm[k] = 0.5*(rods[j].M[k] - p1M[k]);
	      d[k] -= dm[k];
	    }
	    
//	    distance = ld.getDistance(d, rods[j].M, p1M, &s0, &t0);
	    distance = getDistance(d, rods[j].M, p1M, &s0, &t0);
	    
	    cout << "Distance " << distance << endl;
	    cout << "Parameters " << s0 << " " << t0 << endl;
	  }
	  exit(8);
	}
      }
    }
    
}


void hYr::clearClusters(){

    for (int i=0; i<N_r; i++){
	clusters[i].npart = 0;
	clusters[i].firstParticle = 0;
	rods[i].cluster = 0;
	rods[i].nextCl = 0;
	rods[i].prevCl = 0;
	clusterDist[i] = 0;
    }
    clusterDist[N_r] = 0;
}


void hYr::setupList(){

    box->update();
      // dimensions
      for (int i=0;i<3;i++){
	  box->nxCell[i] = maxCells;
	  box->xCell[i] = box->x[i]/box->nxCell[i];
	  if (box->xCell[i] < 1) {
	      box->nxCell[i] = int(box->x[i]);
	      box->xCell[i] = box->x[i]/box->nxCell[i];
	  }
	  box->xCellInv[i] = 1/box->xCell[i];
    }
    box->nCells = box->nxCell[0]*box->nxCell[1]*box->nxCell[2];
    cells = new Cell[box->nCells];

    if (AR<1.0) maxPieces = 6;
    else maxPieces = 6*int(AR);

    Tested = new Rod *[N_r];
    PBuffer = new Piece[N_r*maxPieces];
    for (int i=0; i<N_r; i++){
	rods[i].pieces = new Piece *[maxPieces];
	for (int j=0; j<maxPieces; j++){
	    rods[i].pieces[j] = &(PBuffer[i*maxPieces+j]);
	    rods[i].pieces[j]->rod = &rods[i];
	}

    }

    // put particles into cells
    for (int i=0; i<N_r; i++){
	insertRodIntoCells(rods[i]);
    }
   
    
    // find neighbouring cells
    setUpNeighbours();

//    cout << "# Cells: " << sizeof(int)*28*box->nCells/1000000 << endl;
//    cout << "# PieceBuffer: " << maxPieces*N_r*7*sizeof(int)/1000000 << endl;
//    cout << "# Rods: " << N_r*(maxPieces+7)*sizeof(int)/1000000  << endl;
//    cout << "# Neighbour distance threshold: " << A << endl; 
	    
}



void hYr::setUpNeighbours(){

  int x,y,z,count,cl,cln;
   
  for (int k=1; k<box->nxCell[2]-1; k++){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){
	count = 0;
	for (int zs = -1; zs <= 1; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCell[1])%box->nxCell[1];
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		x = (i+xs + 100*box->nxCell[0])%box->nxCell[0];
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
      }
    }
  }

if (pbc==1){
  for (int k=0; k<box->nxCell[2]; k+=box->nxCell[2]-1){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){
	count = 0;
	for (int zs = -1; zs <= 1; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCell[1])%box->nxCell[1];
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		x = (i+xs + 100*box->nxCell[0])%box->nxCell[0];
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
      }
    }
  }
}
else{
  for (int k=0; k<1; k++){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){
	count = 0;
	for (int zs = -1; zs <0; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    if(phi>0) x = box->nxCell[0]-1-(j+ys + 100*box->nxCell[1])%box->nxCell[1];				//phi=pi/2
	    else x = (j+ys + 100*box->nxCell[1])%box->nxCell[1];										//phi=-pi/2
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		if(phi>0) y = (i+xs + 100*box->nxCell[0])%box->nxCell[0];								//phi=pi/2
		else y = box->nxCell[1]-1-(i+xs + 100*box->nxCell[0])%box->nxCell[0];				//phi=-pi/2
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
	for (int zs = 0; zs <2; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCell[1])%box->nxCell[1];
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		x = (i+xs + 100*box->nxCell[0])%box->nxCell[0];
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
      }
    }
  }
  for (int k=box->nxCell[2]-1; k<box->nxCell[2]; k++){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){
	count = 0;
	for (int zs = 1; zs <2; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    if(phi>0) x = (j+ys + 100*box->nxCell[1])%box->nxCell[1];								//phi=pi/2
	    else x = box->nxCell[0]-1-(j+ys + 100*box->nxCell[1])%box->nxCell[1];					//phi=-pi/2
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		if(phi>0) y = box->nxCell[1]-1-(i+xs + 100*box->nxCell[0])%box->nxCell[0];			//phi=pi/2
		else y = (i+xs + 100*box->nxCell[0])%box->nxCell[0];									//phi=-pi/2
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
	for (int zs = -1; zs <1; zs++){
	  // Modulokonstr. period. Randbed.
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCell[1])%box->nxCell[1];
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (!((xs==0)&&(ys==0)&&(zs==0))){
		x = (i+xs + 100*box->nxCell[0])%box->nxCell[0];
		cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;
	      }
	    }
	  }
	}
      }
    }
  }
}


} 


// reads simulation parameters in case -f<configFile> ("set up from 
// configuration file") is chosen and reads particle positions and 
// directors from configFile
  
void hYr::setUpFromFile() {

  ifstream infile, infile2;
  infile2.precision(14);
  char line[200];
  char quantity[200];
  float number, ARnew;
  int check;
  char outdirt[200];
  char temp[200];
  int doubs;

  cof=2.0;

  // read new simulation parameters from ParameterFile

  infile.open(ParameterFile);
  if (infile.bad()) {
      cerr << "Error " << ParameterFile << " not found\n";
      exit(8);
  }
  while (infile.peek() != EOF) {
    infile.getline(line,sizeof(line));
    check = sscanf(line, "%s%f", quantity, &number);
    if (check != 2) {
      cerr << "Error: wrong line format in " << ParameterFile << endl;
      exit(8);
    }

    else if (strcmp(quantity,"CutoffFactor")==0) cof=number;
    else if (strcmp(quantity,"McSteps")==0) mcSteps=int(number);
    else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
    else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
    else if (strcmp(quantity,"EqSteps")==0) eqSteps=int(number); 
    else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
    else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);
    else if (strcmp(quantity,"MaxDisplacement")==0) maxDisplace=number;
    else if (strcmp(quantity,"MaxRotation")==0) maxRot=number;
	else if (strcmp(quantity,"MaxAxRotation")==0) maxARot=number;
	else if (strcmp(quantity,"MaxVolumeCh")==0) dVmax=number;
    else if (strcmp(quantity,"Seed")==0) seed=int(number);
    else if (strcmp(quantity,"Distance")==0) A=number;
    else if (strcmp(quantity,"NInt")==0) Nint=int(number);
    	else if (strcmp(quantity,"Temperature")==0) tempr=number;
	else if (strcmp(quantity,"Charges")==0) Z=int(number); 
	else if (strcmp(quantity,"BjLength")==0) lam=number;
	else if (strcmp(quantity,"ScreenConst")==0) ka=number;  
	else if (strcmp(quantity,"Diameter")==0) D=number;
	else if (strcmp(quantity,"PitchFac")==0) pf=number;
	else if (strcmp(quantity,"kWaveCalc")==0) kw=int(number);
	else if (strcmp(quantity,"Ensemble")==0) en=int(number);  
	else if (strcmp(quantity,"Pressure")==0) press=number;
	else if (strcmp(quantity,"PBC")==0) pbc=int(number);  
 	else if (strcmp(quantity,"torsionAngle")==0) phi=number;
	else if (strcmp(quantity,"SaltConc")==0) cs=number;
	else if (strcmp(quantity,"L/D")==0) ARnew=number;
	else if (strcmp(quantity,"PBCCh")==0) pbcch=int(number);
	else if (strcmp(quantity,"DoubleSize")==0) doubs=int(number);
	else if (strcmp(quantity,"Twist")==0) twist=int(number);
	else if (strcmp(quantity,"Pi/MaxTwist")==0) maxtwi=number;
	else if (strcmp(quantity,"PressGrad")==0) pressgrad=number;
   }

  infile.close();




   if (maxDisplace >= 1){
      cerr << "Maxdisplace must be smaller than 1" << endl;
      exit(8);
    }

  // read configuration from ConfigFile

  int allSet = 0;
  int count = 0;
  double Rin[3], Min[3];
  double Pin;
  double dump[100000][7];

  infile2.open(ConfigFile);
  
  if (infile2.bad()) {
      cerr << "Error " << ConfigFile << " not found\n";
      exit(8);
  }
  box = new Box;
 
  while (infile2.peek() != EOF) {
      infile2.getline(line,sizeof(line));
      if(allSet == 6){
	  
	  check = sscanf(line, "%le %le %le %le %le %le %le", 
			 &Rin[0], &Rin[1], &Rin[2], &Min[0], &Min[1], &Min[2],&Pin);
	  
	  if (check == 7) { 
	      for (int j=0; j<3; j++){
		  dump[count][j] = Rin[j];
		  dump[count][j+3] = Min[j];
	      }
		dump[count][6] = Pin;
	      count++; 
	  }
	  
	  else { 
	      cerr << "Error: wrong line format in " << ConfigFile << endl;
	      exit(8);
	  }
      }	
      
      else {
	  check = sscanf(line, "%s%f", quantity, &number);
	  if (check < 2) {
	      cerr << "Error: wrong line format in " << ConfigFile << endl;
	      exit(8);
	  }
	  if (strcmp(quantity,"#NumberR")==0) {N_r=int(number); allSet++;}
	  else if (strcmp(quantity,"#DensityR")==0) {allSet++;} 
	  else if (strcmp(quantity,"#L/D")==0) {AR=number; allSet++;}
	  else if (strcmp(quantity,"#Boxx")==0) {box->x[0]=number; allSet++;}
	  else if (strcmp(quantity,"#Boxy")==0) {box->x[1]=number; allSet++;}
	  else if (strcmp(quantity,"#Boxz")==0) {box->x[2]=number; allSet++;}
    }
  }
  infile2.close();

    cof=pow((box->x[0]*box->x[1]*box->x[2]*(doubs+1))/(N_r+doubs*N_r)/M_PI/AR*4.0,1.0/2.0)*sqrt(4.0*M_PI*lam*(Z*(N_r+doubs*N_r)/(box->x[0]*box->x[1]*box->x[2]*(doubs+1))+2.0*cs));
//    cof=max(2,int(2.0*D*sqrt(4.0*M_PI*lam*(Z*(N_r+doubs*N_r)/(box->x[0]*box->x[1]*box->x[2]*(doubs+1))+2.0*cs)))+1);

//    if (en==1 && pbc==1) strcpy(outdirt, "NVT_pbc/");
//    else if (en==1 && pbc==2) strcpy(outdirt, "NVT_tpbc/");
//    else if (en==2 && pbc==1) strcpy(outdirt, "NpT_pbc/");
//    else if (en==2 && pbc==2) strcpy(outdirt, "NpT_tpbc/");
    strcpy(outdirt, "");
    strcat(outdirt, "N");
    sprintf(temp, "%.i",N_r+doubs*N_r);
    strcat(outdirt, temp);
    if (en==1 || en==3) {
	sprintf(temp, "%.*f",0, box->x[0]*box->x[1]*box->x[2]*(doubs+1));
	strcat(outdirt, "_V");
	strcat(outdirt, temp);
    }
    else {
	sprintf(temp, "%.*f",4, press);
	strcat(outdirt, "_p");
	strcat(outdirt, temp);
	if (pressgrad > 0.0) {
	    sprintf(temp, "%.*e",2, pressgrad);
	    strcat(outdirt, "_pg");
	    strcat(outdirt, temp);
	}

    }
    strcat(outdirt, "_AR");
    sprintf(temp, "%.*f", 0, ARnew);
    strcat(outdirt, temp);
    strcat(outdirt, "_T");
    sprintf(temp, "%.*f", 1, tempr);
    strcat(outdirt, temp);
    strcat(outdirt, "_Z");
    sprintf(temp, "%.i", Z);
    strcat(outdirt, temp);
    strcat(outdirt, "_npc");
    sprintf(temp, "%.i", Nint+1);
    strcat(outdirt, temp);
    strcat(outdirt, "_l");
    sprintf(temp, "%.*f", 2, lam);
    strcat(outdirt, temp);
    strcat(outdirt, "_cs");
    sprintf(temp, "%.*f", 3, cs);
    strcat(outdirt, temp);
    strcat(outdirt, "_pf");
    sprintf(temp, "%.*f", 1, pf);
    strcat(outdirt, temp);
    if (pbc==2){
     	strcat(outdirt, "_phi");
   	sprintf(temp, "%.*f", 1, phi);
    	strcat(outdirt, temp);
    }

	struct stat st = {0};
	int status;
strcpy(outdir,outdirt);
	if (stat(outdir, &st) == -1) {

   	 status = mkdir(outdir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	}
 

    ofstream out(strcat(outdirt, "/Outfile"));


  if (count==0) {
      cerr << "Error: parameter missing in " << ConfigFile << endl;
      exit(8);
  }
  
  else if (count!=N_r) {
      cerr << "Error: wrong number of coordinates in " << ConfigFile << endl;
      exit(8);
  }
  
  out << "# read " << ConfigFile << endl;
  
if( doubs == 0){
  rods = new Rod[N_r];
  for (int i=0; i<N_r; i++){
      for (int j=0; j<3; j++){
	  rods[i].R[j] = dump[i][j];
	  rods[i].M[j] = dump[i][j+3]*(ARnew/AR);
      }
	rods[i].P = dump[i][6];
  }
}
else{
  rods = new Rod[2*N_r];
  for (int i=0; i<N_r; i++){
      for (int j=0; j<2; j++){
	  rods[i].R[j] = dump[i][j];
	  rods[i].M[j] = dump[i][j+3]*(ARnew/AR);
	  rods[i+N_r].R[j] = dump[i][j];
	  rods[i+N_r].M[j] = dump[i][j+3]*(ARnew/AR);
      }
      for (int j=2; j<3; j++){
	  rods[i].R[j] = dump[i][j];
	  rods[i].M[j] = dump[i][j+3]*(ARnew/AR);
	  rods[i+N_r].R[j] = dump[i][j]+box->x[2];
	  rods[i+N_r].M[j] = dump[i][j+3]*(ARnew/AR);
      }
	rods[i].P = dump[i][6];
	rods[i+N_r].P = dump[i][6];
  }
  N_r *= 2;
  box->x[2] *= 2;
}

/*if (pbcch==12 || pbcch==21){
  double rad, rad1, vphi, vphi1, psi0, psi1, philoc;
  int ch;
ch = -int((float(pbcch)-16.5)/4.5);
	for (int i=0; i<N_r; i++){
		philoc = phi/box->x[2]*rods[i].R[2];
		
		rad = sqrt(rods[i].R[0]*rods[i].R[0]+rods[i].R[1]*rods[i].R[1]);
		rad1=sqrt((rods[i].R[0]-0.5*rods[i].M[0])*(rods[i].R[0]-0.5*rods[i].M[0])+(rods[i].R[1]-0.5*rods[i].M[1])*(rods[i].R[1]-0.5*rods[i].M[1]));
		vphi = atan2(rods[i].R[1],rods[i].R[0]);
		vphi1=atan2((rods[i].R[1]-0.5*rods[i].M[1]),(rods[i].R[0]-0.5*rods[i].M[0]));
		rods[i].R[0] = rad*cos(vphi+ch*philoc*M_PI);
		rods[i].R[1] = rad*sin(vphi+ch*philoc*M_PI);
		
		psi0 = acos(-rods[i].M[1]/sqrt(rods[i].M[1]*rods[i].M[1]+rods[i].M[2]*rods[i].M[2]));
		rods[i].M[0] = 2.0*(rods[i].R[0]-rad1*cos(vphi1+ch*philoc*M_PI));
		rods[i].M[1] = 2.0*(rods[i].R[1]-rad1*sin(vphi1+ch*philoc*M_PI));
	
		psi1 = acos(-rods[i].M[1]/sqrt(rods[i].M[1]*rods[i].M[1]+rods[i].M[2]*rods[i].M[2]));
		rods[i].P = rods[i].P+psi1-psi0;
	}
}
*/
 
  
  iran = &seed;
  
box->update();

if (ARnew != AR){
double oldVol, newVol, ARold;
    setupList();
    ARold = AR;
    AR = ARnew;
    oldVol = box->V;
    newVol = oldVol*((ARnew+1)/(ARold+1))*((ARnew+1)/(ARold+1));
    compressBox(oldVol,newVol,1);
    newVol = oldVol*((ARnew+1)/(ARold+1));
    compressBox(oldVol,newVol,2);
}

 initialize();

  bool overlap = false;
  
  for (int i=0; i<N_r; i++){
      for (int j=0; j<i; j++){
	  overlap = checkOverlapRods(rods[i], rods[j]);
	  if (overlap){
	      cerr << "Error: Overlapping configuration!" << endl;
	      cerr <<rods[j].R[0]<<" "<<rods[j].R[1]<<" "<<rods[j].R[2]<<" " 
		   <<rods[j].M[0]<<" "<<rods[j].M[1]<<" "<<rods[j].M[2]<< endl;
	      cerr <<rods[i].R[0]<<" "<<rods[i].R[1]<<" "<<rods[i].R[2]<<" " 
		   <<rods[i].M[0]<<" "<<rods[i].M[1]<<" "<<rods[i].M[2]<< endl;
	      exit(8);
	  }
      }
  }
 


  out << "# read parameters and configuration" << endl;
  out << "# L/D " << AR << " N rods " << N_r << endl;
  out << "# Eq " << eqSteps << " EqEval " << eqEvalSteps << " EqSnap " 
       << eqSnapSteps << endl;
  out << "# Mc " << mcSteps << " McEval " << mcEvalSteps << " McSnap " 
       << mcSnapSteps << endl;
  out << "# seed " << seed << " Displ " << maxDisplace << " Rot " 
       << maxRot << endl;
out.close();
}





void hYr::insertRodIntoCells(Rod &p){

    double s[3], c, t, save;
    double delta = 1e-2;
    int Offs, ch;
    double Minv;
    p.npieces = 0;
    
    // loop over (y,z)-planes
    if (p.M[0]*p.M[0] > 1e-8){
	Minv=1/p.M[0];
	Offs = int(AR*box->xCellInv[0]);
	for (int i=-Offs; i<box->nxCell[0]+Offs; i++){
	    c = i*box->xCell[0]-box->halfx[0];
	    t = (c-p.R[0])*Minv;
	    if (t*t < 0.25){
		s[0] = c+delta;
		s[1] = p.R[1] + t*p.M[1];
		s[2] = p.R[2] + t*p.M[2];
ch = int(s[2]*box->halfxInv[2]);
if(pbc==2 &&  ch != 0){
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
		s[0] = c-delta;
if(pbc==2 &&  ch != 0){
	s[1] = p.R[1] + t*p.M[1];
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
	    }
	}
    }
    
    // loop over (x,z)-planes
    if (p.M[1]*p.M[1] > 1e-8){
	Minv=1/p.M[1];
	Offs = int(AR*box->xCellInv[1]);
	for (int i=-Offs; i<box->nxCell[1]+Offs; i++){
	    c = i*box->xCell[1]-box->halfx[1];
	    t = (c-p.R[1])*Minv;
	    if (t*t < 0.25){
		s[0] = p.R[0] + t*p.M[0];
		s[1] = c+delta;
		s[2] = p.R[2] + t*p.M[2];
ch = int(s[2]*box->halfxInv[2]);
if(pbc==2 &&  ch != 0){
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
		s[1] = c-delta;
if(pbc==2 &&  ch != 0){
	s[0] = p.R[0] + t*p.M[0];
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
	    }
	}
    }

    // loop over (x,y)-planes
    if (p.M[2]*p.M[2] > 1e-8){
	Minv=1/p.M[2];
	Offs = int(AR*box->xCellInv[2]);
	for (int i=-Offs; i<box->nxCell[2]+Offs; i++){
	    c = i*box->xCell[2]-box->halfx[2];
	    t = (c-p.R[2])*Minv;
	    if (t*t < 0.25){
		s[0] = p.R[0] + t*p.M[0];
		s[1] = p.R[1] + t*p.M[1];
		s[2] = c+delta;
ch = int(s[2]*box->halfxInv[2]);
if(pbc==2 &&  ch != 0){
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
		s[2] = c-delta;
ch = int(s[2]*box->halfxInv[2]);
if(pbc==2 &&  ch != 0){
	s[0] = p.R[0] + t*p.M[0];
	s[1] = p.R[1] + t*p.M[1];
	save = s[0];
	s[0] = int(2*phi)*ch*s[1];
	s[1] = -int(2*phi)*ch*save;
}

		tryInsert(p,s);
	    }
	}
    }
}


void hYr::tryInsert(Rod &p, double s[3]){
    
    int k[3],icell;
double s_in[3];
    
    for (int i=0; i<3; i++){
	s[i] = -int(s[i]*box->halfxInv[i])*box->halfx[i] 
	    + fmod(s[i],box->halfx[i]);
	k[i] = int((s[i]+box->halfx[i])*box->xCellInv[i]);
    }
    icell = k[0] + k[1]*box->nxCell[0] 
	+ k[2]*box->nxCell[0]*box->nxCell[1];  

    if ((icell < box->nCells)&&(icell >= 0)) {
	if (p.npieces < maxPieces){

	    p.addPiece(cells[icell]);
	}
	else {
	    cerr << "Maximum number of pieces!" << endl;
	    exit(8);
	}
    }
    else {
	cerr << "Error: problem with rod ->cell assignment!" << endl;
	cerr << icell << " " << k[0] << " " << k[1] << " " << k[2] << endl;
	cerr<<p.R[0] << " " <<p.R[1] << " " <<p.R[2] << " " <<endl;
	cerr<<s_in[0] << " " <<s_in[1] << " " <<s_in[2] << " " <<endl;
	cerr<<box->halfx[0]<<" "<<box->xCellInv[0]<<endl;
	exit(8);
    }
    
}
