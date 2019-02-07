/*************************************************************************
 * hscnpt: NVT Monte Carlo for hard spherocylinders                      *  
 * with a cell system for large aspect ratios                            *
 * Measures distribution of clusters of rods with distance < A           *
 * and distribution of distances between particles                       *
 *                                                                       *
 * Options:                                                              *
 * -n: start a new Configuration.                                        *
 * -f<filename>: start from a configuration file (output of an earlier   *
 *               simulation).                                            *
 *               Filename contains: Nrods,           L/D, rho,           *
 *               center of mass vectors and directors of all rods        *
 *                                                                       *
 * -p<parameterfile>: in case -n is chosen, parameterfile contains:      *
 *                    N_r per dimension, box dimensions, L/D, A          *
 *                    #Blocks equilibration, #Blocks MC,                 *
 *                    #evaluations during MC,#evaluations during equil.  *
 *                    #Snapshot Configurations, maximal                  *
 *                    displacement, maximal rotation,and a seed.         *
 *                    In case -f is chosen, parameterfile contains:      *
 *                    all of the above apart from N, L/D, rho            *
 * 140107 Tanja                                                          *
 *************************************************************************/


#define maxCells 500  //maximum number of cells per dimension
#define nbins 250    //bins for distance distribution 
#define sbins 100    //bins for distribution of contact points

class hYrtls;
class Rod;		
class Cell;
class Piece;
class Box;
class Cluster;
//class lineDist;

class hYr {

 /* This module computes the distance between two linesegments. The segments
   are represented by two vectors each, R0, M0 and  R1, M1 with 
   line 0: L0=R0+s*M0 und line 1: L1=R0+t*M1, s,t from [0,1].
   getDistance takes R0-R1, M0 and M1 and returns square of distance. 
   051117
*/ 

//class lineDist {

 private:
 
  double M0[3], M1[3], DD[3]; //directions of lines and distance of startpoints
  double a,b,c,d,e,f;        //coefficients of distance function
  double smin, tmin;         //coordinates of minimum in distance function
  double ddelta;              //determinant
  double distance;

  void dinitialize();
  void parallel();
  void nonparallel();
  void calcdist();
  void quad0();
  void quad1();
  void quad2();
  void quad3();
  void quad4();
  void quad5();
  void quad6();
  void quad7();
  void quad8();

 public:
 
  int argc;
  char **argv;

  bool newSetUp;
  const char *ParameterFile, *ConfigFile;  // Input files
  char outdir[200];
  std::ofstream out;

  void usage();
  void run();

  Rod *rods;
  Cell *cells;
  Box *box;
  hYrtls *TLS;

  // Simulation Data

  int transAccept, rotAccept, axrotAccept, volAccept;
  int transMoves, rotMoves, axrotMoves, volMoves;
  double transRate, rotRate, axrotRate, volRate;
  long seed;          //seed for ran3
  long *iran;         //pointer to seed

  int eqSteps,reSteps;          
  int eqEvalSteps;      // # blockaverages taken during equilibration.
  int eqSnapSteps;      // # configuration snapshots taken during equil.
  int mcSteps;
  int mcEvalSteps;      // # measurements for averages
  int mcSnapSteps;      // # snapshots during mc

  double maxDisplace, maxRot, dVmax, maxARot;     

  // Physical data

  // input
  int Nin_r;  // suggests number of particles per dimension on input
  double AR;         // Aspect ratio (length/diameter)
  double AR2, AR12, ARInv, hAR;  // AR^2, (AR+1)^2, 1/AR, AR/2
  double A, A2, ARA2, A12; // distance, A^2, (A+AR)^2, (A+1)^2

  // simulation
  int N_r;       // particle number
  double Nrinv;       // 1/N_r
  double NN1hInv;           // 2/N_r*(N_r+1)
  double Vr;   // Volumes of rods and spheres
  double boxz0;
  double currOP, director[3]; // current S2, current director
  long allContacts;           // Number of contacts in system
  long clustContacts;         // Number of contacts in largest cluster
  long nlargest;              // Number of rods in largest cluster
//  lineDist ld;
  Rod **Tested;
  Piece *PBuffer;
  int maxPieces;      // maximum number of pieces per rod
  int Nint;
  double NintInv;
  double L,Linv;
  double D,Dhalf;
//  double Dhalf;
  double pf;
  double cof,dpcof,cofka,dpcof2,cofka2;
  int kw;
  double kwave,tkwave;
  double ka;
  double cs;
  double tempr;
  int Z;
  double Zpc;
  double lam;
  double shift,tzl,nt;
  double pitch;
  double press;
  double pressgrad;
  double phi;
  int i2phi;
  int en;
  int pbc;
  int pbcch;
  int twist;
  double maxtwi;
  int N_l;
  int *nlist;
  int *nbn;
  double *ener;
//  double *func;
//  double *func2;
  double Ecalc;
  double twAngle;
  double meantwAngle;
  double eqPitch;
  double *twistprob;
  int *locOP;
  double **delta;

  // current configuration
  int *clusterDist;          // distribution of clusters
  Cluster *clusters;         //Clusterarray
  int *distances;
 
  // averaged data
  int measureCount;
  double meanOP; 
  double meanV;
  int *meanClusterDist;      // mean distribution of clusters
  int *gofr;                 // Radial distribution function histogram
  double *g2ofr;             // Legendre P_2 as a function of separation of rods
  int **gofrpp;                 // Radial distribution function histogram
  double **g2ofrpp;             // Legendre P_2 as a function of separation of rods
  long *sHist;               // histogram of contact positions
  int *meanDistances;
  double binWidth;
  long nPerc;                // Number of percolating configurations
  long endCont;              // Number of contacts at end of rods
  double clustFrac;          // Mean fraction of rods that belong to the largest cluster
  double allAveCont;         // Contacts per rod in whole system
  double clustAveCont;       // Contacts per rod in largest cluster

  //set up
  void setUpNew();
  void setUpFromFile();
  void setUpPar();
  void initialize();

  void setupList();
  void setupListRun();
  void setUpNeighbours();
  void insertRodIntoCells(Rod &p);
  void tryInsert(Rod &p, double s[3]);

  void clearClusters();
  void findClusters();         // identifies clusters of rods closer than A
  void measureClusterDist();   // measures cluster distribution
  bool Neighbours(Rod &p1, Rod &p2);
  void measureDistances();

  void packAAAExpandAndShift();

  // simulation

  void equil();
  void MC();
//  void workerFunc(int pNr, boost::thread_group &g, boost::thread *t);
  void workerFunc(int pNr);
  void displaceParticle(Rod &p);
  bool particleMove(int pNr);
  bool orientationalMove(int pNr);
  bool axisrotMove(int pNr);

  void randomVector(double v[3]);

  void volumeMove();
  void neighborlist();
  void neighborlistTwist();
  double twistPitch();
  double twistBox(double twPhi);
  double calculateEnergyTwist(int pNr, double twphi);

//  double twistPitch();

  void nematicOPconfig();      // measures P2
  void locnematicOPconfig();      // measures P2
  int percolate();             // test whether rods percolate 

  bool checkOverlapRodAll(Rod &p);
  bool checkOverlapRods(Rod &p1, Rod &p2);
  double integrate(double func[], int Nint);
  double calculateEnergy(int pNr);
  double calculateEnergyPair(int pNr, int pNr2);
  void calculateEnergyAll();
  void calculateVirial();
  void compressBox(double oV, double nV, int dir);
  void compressBoxV(double oV, double nV, int dir);

  // io

  void writeOutConfig(const char name[16]);
  void writeOutSimData(const char name[16]);
  void writeOutLargestCluster(const char name[16]);
  void writeOutDistances(const char name[16]);
  void writeOutContacts(const char name[16]);

  //from lineDist
  double getDistance(double d[3], double m0[3], double m1[3], double *s0, double *t0);


};

class hYrtls : public hYr {

 public:
//  int N_r;


//  double *func;
//  double *func2;
  double flip;
  bool accepted;
  int move;

  double calculateEnergy(int pNr);
/*  void displaceParticle(Rod &p);
  bool particleMove(int pNr);
  bool orientationalMove(int pNr);
  bool axisrotMove(int pNr);

  void randomVector(double v[3]);

  double calculateEnergy(int pNr);
  bool checkOverlapRodAll(Rod &p);
  bool checkOverlapRods(Rod &p1, Rod &p2);
*/
};


class Cell {

 public:
 
  Piece *firstPiece;
  Cell *neighbours[26];      // neighbouring cells 
  Cell(void){
      for (int i=0; i<26; i++){
	  neighbours[i] = 0;
      }
      firstPiece = 0;
  }
  
};

class Cluster {

 public:

  int npart;
  Rod *firstParticle;
  Cluster(void){
      firstParticle = 0;
  }
  

};


class Piece {
    
 public:
    
    Piece *next, *prev;    // pieces are linked in list per cell
    Cell *cell;
    Rod *rod;

    Piece(){
	next = 0;
	prev = 0;
	cell = 0;
    }

    //void insertToCell(Cell &c, int in[3]){
    void insertToCell(Cell &c){

	prev = 0;               // initialize pointer
	next = c.firstPiece;    // -> firstPiece must be initialized to 0!
	if (next) next->prev=this;
	c.firstPiece=this;
	cell=&c;
    };

    void removeFromCell(){
	  
	if (prev) prev->next=next;
	else cell->firstPiece=next;
	if (next) next->prev=prev;
	next = 0;
	prev = 0;
	cell = 0;
    };

};




// does not create and delete pieces -- this must be done by hsc.
class Rod {

 public:

    double R[3], M[3];     // center of mass, orientation
    double P;			// azimuthal angle
    Piece **pieces;        // pieces of particle in different cells
    int npieces;           // current number of pieces

    Rod *nextCl, *prevCl; //next, previous in Cluster
    Cluster *cluster;
    Rod **neighbours; // closer than A 
    int nneb;

    // Next three elements for percolation detection  
    double copyR[3];   // Position of rod in contiguous copy of cluster
    bool done;         // true = rod has already been placed in copy of cluster
    Rod *nextPlace;    // Next rod in order of placement in copy of cluster

    Rod(){
	npieces = 0;
	nneb = 0;
    }

    void clearPieces(){
	for (int i=0; i<npieces; i++){
	    pieces[i]->removeFromCell();
	}
	npieces = 0;
    }

    void addPiece(Cell &c){
	for (int i=0; i<npieces; i++){
	    if (pieces[i]->cell == &c) return;
	}
	pieces[npieces]->insertToCell(c);
	npieces ++;
    }

    void insertToCluster(Cluster &c){
	
	prevCl = 0;
	nextCl = c.firstParticle; 
	if (nextCl) nextCl->prevCl=this;
	c.firstParticle=this;
	cluster=&c;
	c.npart++;
	
    };
    
    void moveBetweenClusters(Cluster &from, Cluster &to){
	
	if (prevCl) prevCl->nextCl=nextCl;
	else from.firstParticle=nextCl;
	if (nextCl) nextCl->prevCl=prevCl;
	
	prevCl=0;
	nextCl=to.firstParticle;
	if (nextCl) nextCl->prevCl=this;
	to.firstParticle=this;
	cluster=&to;
	
	from.npart--;
	to.npart++;
	
    };

    
};

class Box {

 public:

  double x[3], halfx[3], halfxInv[3];  // dimensions
  double V;      // Volume
  double xCell[3]; //Cell dimensions 
  double xCellInv[3];
  int nxCell[3]; // number of cells
  int nCells;

  void update(){
      for (int i=0; i<3; i++){
	  halfx[i] = 0.5*x[i];
	  halfxInv[i] = 1.0/halfx[i];
      }
      V = x[0]*x[1]*x[2];
  }
  
};

