#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "circletest.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
/*
   Framework code for using brute force method to triangulate points
*/

double sumsquare(double x, double y){
  double sumsqr;
  sumsqr=x*x+y*y;
  
  return sumsqr;
}

double distance(double x1, double x2){
  double dist;
  dist=x2-x1;
  
  return dist;
}

double magnitude( double x0, double y0, double x1, double y1){
  double vectormagnitude;

  vectormagnitude=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
  
  return vectormagnitude;
}

double determinant2(double e11, double e12, double e21, double e22){
  double det;
  det=e11*e22-e12*e21;
  return det;
}

double determinant3(double e11, double e12, double e13, 
		    double e21, double e22, double e23,
		    double e31, double e32, double e33){
  double det;
  det=e11*(e22*e33-e32*e23)-e12*(e21*e33-e31*e23)+e13*(e21*e32-e31*e22);
  return det;
}

// Luke Mcculloch
// Circle test
// Nov 2011
//*test point is (x4,y4)*/
//*circle center is (xp,yp)*/
int circle_test(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
  int flag;
  double ssqr1, ssqr2, ssqr3;
  double a;
  double bx;
  double by;
  double c;
  double xo,yo;
  double radius, dist;

  ssqr1=sumsquare(x1,y1);
  ssqr2=sumsquare(x2,y2);
  ssqr3=sumsquare(x3,y3);
  
  a=determinant3(x1,y1,1.0, x2,y2,1.0, x3,y3,1.0);


  bx=determinant3(ssqr1,y1,1.0, ssqr2,y2,1.0, ssqr3,y3,1.0);

  by=determinant3(ssqr1,x1,1.0, ssqr2,x2,1.0, ssqr3,x3,1.0);

  c=determinant3(ssqr1,x1,y1, ssqr2,x2,y2, ssqr3,x3,y3);

  xo=bx/(2.0*a);
  yo=-by/(2.0*a);

  /*radius of circle*/
  radius=(sqrt((bx*bx)+(by*by)+(4.0*a*c)))/(2.0*ABS(a));

  if (ABS(a)<10e-12){
    //printf("\na= %e", a);
    xo=(x1+x2+x3)/3.0;
    yo=(y1+y2+y3)/3.0;
    radius=10e20;
  }

  /*distance from circle center to test point*/
  dist=magnitude(xo,yo,x4,y4);




  if(dist>=(radius)){
    /*Pt is not in the circle, return false*/
    flag=0;
  }
  else if(dist<(radius)){
    /*Pt is in the circle, return true*/
    flag=1;
  }
  return flag;

}
