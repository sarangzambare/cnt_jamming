/* This module computes the distance between two linesegments. The segments
   are represented by two vectors each, R0, M0 and  R1, M1 with 
   line 0: L0=R0+s*M0 und line 1: L1=R0+t*M1, s,t from [0,1].
   GetDistance takes R0-R1, M0 and M1 and returns square of distance. 
   051117
*/ 

class lineDist {

 public:

  double getDistance(double d[3], double m0[3], double m1[3], double *s0, double *t0);
 
 private:
 
  double M0[3], M1[3], D[3]; //directions of lines and distance of startpoints 
  double a,b,c,d,e,f;        //coefficients of distance function
  double smin, tmin;         //coordinates of minimum in distance function
  double delta;              //determinant
  double distance;

  void initialize();
  void parallel();
  void nonparallel();
  void calcdist ();
  void quad0();
  void quad1();
  void quad2();
  void quad3();
  void quad4();
  void quad5();
  void quad6();
  void quad7();
  void quad8();

};






