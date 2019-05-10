#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))




double sumsquares(double x, double y){
  double sumsqr;
  sumsqr=x*x+y*y;
  
  return sumsqr;
}

double distances(double x1, double x2){
  double dist;
  dist=x2-x1;
  
  return dist;
}

double magnitudes( double x0, double y0, double x1, double y1){
  double vectormagnitude;

  vectormagnitude=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
  
  return vectormagnitude;
}

double dot_products(double xo, double yo, double xp, double yp){
  double dp;
  //dp=0.0;
  //for (i=0; i < 3; i++){
  //dp=dp+testv1[i]*testv2[i];
    dp=xo*xp+yo*yp;
    //}  

  return dp;
}


int search(int seed, double xt, double yt, double x[], double y[], int tri[][3], int nbr[][3]){
  int n0, n1, i, t;
  int xp, yp;
  int outflag, grad, inflag, count, loop;
  double mag, xnormal, ynormal, maxnew, maxold, flag;
  double midx, midy, xtv, ytv;

  /*find the triangle that contains the pt given*/
  maxold=0.0;
  maxnew=0.0;
  flag=0.0;
  loop=1;

  /*Start a loop at the beggining*/

  t=0;
  while (loop!=0){

    //printf("\n\nWe move to triangle %d.",t);
    inflag=0;
    outflag=0;
    count=0;
    //printf("\ninflag = %d",inflag);
    for (i=0; i<3; i++){
      count=count+1;
      /*Each n (n0 and n1) uniquely determines an x,y element node*/
      n0=tri[t][i];
      n1=tri[t][(i+1)%3];
      
      //printf("\nCurrent Test Points");
      //printf("\nx,y= %e,%e",x[n0],y[n0]);
      //printf("\nx,y= %e,%e",x[n1],y[n1]);
      
      /*Get vector normals on each triangle line segment*/
      mag=magnitudes(x[n1],y[n1],x[n0],y[n0]);
      xnormal = distances(x[n1],x[n0])/mag;
      ynormal = -distances(y[n1],y[n0])/mag;
      //printf("\n -ynormal = %e, xnormal = %e", ynormal, xnormal);
      /*Compute the Dot Product*/
      //1. Get the x,y Midpoint of the segment
      midx=x[n0]+0.5*distances(x[n0],x[n1]);
      midy=y[n0]+0.5*distances(y[n0],y[n1]);
      //2. Get the Vector to the test point
      xtv=xt-midx;
      ytv=yt-midy;
      
      
      flag=dot_products(ynormal,xnormal,xtv,ytv);////error compute vectors
      //printf("\ndot_products give %e", flag);
      
      if (flag<=0.0){
	/*The Pt is inside the triangle so far*/
	inflag=inflag+1;
	//printf("\ninflag = %d",inflag);
	if (inflag>2){
	  /*We have found the triangle*/
	  //printf("\nLuke Really Found the Triangle %d", t);
	  seed=t;
	  loop=0;
	}
      }
      else if (flag>0.0){
	/*The point is not in the triangle, invoke tri-gradient*/
	/*Determine which face points most directly to the point*/
	maxold=flag;
	maxnew=MAX(maxnew,flag);
	outflag=outflag+1;
	//printf("\nGot down here with no result");
	if(maxnew>=maxold){
	  t=nbr[t][i];
	  //printf("\ntest print, not the tri 3, t= %d", t);
	}
	
	// switch (){
	// case1:
	
	//   break;
	// case2:
	
	//   break;
	// case3:
	
	//   break;
	// default:
	//   //printf("\nerror, min max variable 'loc' not specified correctly ");	 
	// }
      
	if (count>2){
	  /*We are outside this triangle*/
	  //printf("\ntest print, not a tri 3");
	  
	  continue;
	  /*Find the direction of largest gradient*/
	}
      }
    
      
    }
  }
  return(t);
}
