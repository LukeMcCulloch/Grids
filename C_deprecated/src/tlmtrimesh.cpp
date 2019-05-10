#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "List.h"
#include "circletest.h"


//----------------------------------------------------------------------------------------------------
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))


//----------------------------------------------------------------------------------------------------
extern int search(int seed, double xt, double yt, double x[], double y[], int tri[][3], int nbr[][3]);
extern void make_nbrs(int nn, int nt, int tri[][3], int nbr[][3]);
//int circle_test(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
int circle_test(double (&a)[2], double (&b)[2], 
                      double (&c)[2], double (&d)[2]);



//----------------------------------------------------------------------------------------------------
double cross_product(double x1,double y1, double x2, double y2){
  double cp;
  cp=x1*y2-x2*y1;
  
  return cp;
}

double dotproduct(double xo, double yo, double xp, double yp){
  double dp;
  //dp=0.0;
  //for (i=0; i < 3; i++){
  //dp=dp+testv1[i]*testv2[i];
    dp=xo*xp+yo*yp;
    //}  

  return dp;
}

double mag( double x0, double y0, double x1, double y1){
  double vectormagnitude;

  vectormagnitude=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
  
  return vectormagnitude;
}




//----------------------------------------------------------------------------------------------------
/*Function to find the node in triangle j which is not shared with neighbor triangle m */
// Function is reverse denoted from Lawson test for confusion. swap j and m in function call.
// We want the point that is a neighbor node that is not a seed node.
// Usage: ( m=seed, j=neighbor, tri is the tri array )
int outnode(int m, int j, int tri[][3]){

  if( (tri[j][0] != tri[m][0]) && (tri[j][0] != tri[m][1])  && (tri[j][0] != tri[m][2]) ){
    return 0;
  }
  else if( (tri[j][1] != tri[m][0]) && (tri[j][1] != tri[m][1])  && (tri[j][1] != tri[m][2]) ){
    return 1;
  }
  else{
    return 2;
  }

}
//----------------------------------------------------------------------------------------------------
/*EdgeFlipper Function*/
void boundary_flip(int j,int m,int testflip[6],int tri[][3],int nbr[][3], double tempx[], double tempy[], int &areaflag){
  int aa,bb;
  int i, itestj, itestm, savei, l, flag;
  int p0,p1,p2,p3;
  int a,b,c,d;
  int cflagj, cflagm;
  double area1,area2;
  double xt, yt;
  double x1,x2,x3,y1,y2,y3;
  
  // get the node in m that's not in j :: a = either 0,1,or2
  aa=outnode(j,m,tri); /// Print this right fool
  //printf("\nThe node in m that's not in j is aa= %d, tri[m][aa] = %d ",aa, tri[m][aa]);
  
  //get the node in j that's not in m :: b = either 0,1,or2
  bb=outnode(m,j,tri);
  //printf("\nThe node in j that's not in m is bb= %d, tri[j][bb] = %d ",bb, tri[j][bb]);

  //--------------------------------------------------
  // if(tri[m][aa]==tri[j][bb]){
  //   //printf("\nTotal              ");
  //   //printf("\n      Logic        ");
  //   //printf("\n            Failure");
  //   //printf("\nOutNodes found equal to each other");
  //   exit(0);
  // }
  //--------------------------------------------------

  // The outnode a in m becomes part of j as well
  // The outnode b in j becomes part of m as well
  // I believe that either way the winding is perserved.
  // Need to find a way to test this
  //for (t=0; t<ntri; t++){

  //populate testflip with tri j and m
  for (i=0; i<3; i++){
    testflip[i]=tri[j][i];
    testflip[i+3]=tri[m][i];
  }

  // Index the seed nodes (0,1,2,3)
  p0=b;
  p1=(b+1)%4;
  p2=(b+2)%4;//a
  p3=(b+3)%4;


  switch (bb){
  case 0:
    a=0;
    b=1;
    d=2;
    
    break;
  case 1:
    a=1;
    b=2;
    d=0;
    
    break;
  case 2:
    a=2;
    b=0;
    d=1;
    
    break;
  default:
    printf("\nError, switch j not specified correctly j= %d",j);	 
  }

  
  c=aa;

  // not flipped yet
  a=tri[j][a];
  d=tri[j][d];
  b=tri[j][b];
  c=tri[m][c];

  //Area Test Combined into Edge flipper
  
  area1 = (tempx[c]-tempx[b])*(tempy[a]-tempy[b])-(tempx[a]-tempx[b])*(tempy[c]-tempy[b]);
  area2 = (tempx[a]-tempx[d])*(tempy[c]-tempy[d])-(tempx[c]-tempx[d])*(tempy[a]-tempy[d]);


  if ( (area1>0.0) && (area2>0.0)){

     areaflag=0;

    tri[m][0]=c;
    tri[m][1]=d;
    tri[m][2]=a;

    tri[j][0]=a;
    tri[j][1]=b;
    tri[j][2]=c;

    }
  //printf("\n");
  return;
}
//----------------------------------------------------------------------------------------------------
/*EdgeFlipper Function*/
void edge_flip(int j,int m,int testflip[6],int tri[][3],int nbr[][3], double tempx[], double tempy[], int &areaflag){
  int aa,bb;
  int i, itestj, itestm, savei, l, flag;
  int p0,p1,p2,p3;
  int a,b,c,d;
  int cflagj, cflagm;
  double area1,area2;
  double xt, yt;
  double x1,x2,x3,y1,y2,y3;
  
  // get the node in m that's not in j :: a = either 0,1,or2
  aa=outnode(j,m,tri); /// Print this right fool
  //printf("\nThe node in m that's not in j is aa= %d, tri[m][aa] = %d ",aa, tri[m][aa]);
  
  //get the node in j that's not in m :: b = either 0,1,or2
  bb=outnode(m,j,tri);
  //printf("\nThe node in j that's not in m is bb= %d, tri[j][bb] = %d ",bb, tri[j][bb]);

  //populate testflip with tri j and m
  for (i=0; i<3; i++){
    testflip[i]=tri[j][i];
    testflip[i+3]=tri[m][i];
  }

  // Index the seed nodes (0,1,2,3)
  p0=b;
  p1=(b+1)%4;
  p2=(b+2)%4;//a
  p3=(b+3)%4;


  switch (bb){
  case 0:
    a=0;
    b=1;
    d=2;
    
    break;
  case 1:
    a=1;
    b=2;
    d=0;
    
    break;
  case 2:
    a=2;
    b=0;
    d=1;
    
    break;
  default:
    printf("\nError, switch j not specified correctly j= %d",j);	 
  }

  
  switch (aa){
  case 0:

    c=0;
    
    break;
  case 1:

    c=1;
    
    break;
  case 2:

    c=2;
    
    break;
  default:
    printf("\nError, switch m not specified correctly m = %d", m);	 
  }
  // not flipped yet
  a=tri[j][a];
  d=tri[j][d];
  b=tri[j][b];
  c=tri[m][c];

  //Area Test Combined into Edge flipper
  
  area1 = (tempx[c]-tempx[b])*(tempy[a]-tempy[b])-(tempx[a]-tempx[b])*(tempy[c]-tempy[b]);
  area2 = (tempx[a]-tempx[d])*(tempy[c]-tempy[d])-(tempx[c]-tempx[d])*(tempy[a]-tempy[d]);

  //printf("\narea1= %e, area2= %e\n",area1, area2);


  // Now do the Circle Test
  // on the new J with outlier m
  x1=tempx[a];
  x2=tempx[b];
  x3=tempx[c];
	  
  y1=tempy[a];
  y2=tempy[b];
  y3=tempy[c];
  //then get the test point - unshared m
  xt=tempx[d];
  yt=tempy[d];
  //cflagj=circle_test_tlm(x1,y1,x2,y2,x3,y3,xt,yt);
  double tp1[2] = {x1,y1};
  double tp2[2] = {x2,y2};
  double tp3[2] = {x3,y3};
  double tpt[2] = {xt,yt};
  cflagj=circle_test( tp1,tp2,tp3,tpt );
  

  // Now do the Circle Test
  // on the new m with outlier j
  x1=tempx[c];
  x2=tempx[d];
  x3=tempx[a];
	  
  y1=tempy[c];
  y2=tempy[d];
  y3=tempy[a];
  //then get the test point - unshared j
  xt=tempx[b];
  yt=tempy[b];
  //cflagm=circle_test_tlm(x1,y1,x2,y2,x3,y3,xt,yt);
  tp1[0] = x1;
  tp1[1] = y1;
  tp2[0] = x2;
  tp2[1] = y2;
  tp3[0] = x3;
  tp3[1] = y3;
  tpt[0] = xt;
  tpt[1] = yt;
  cflagm=circle_test( tp1,tp2,tp3,tpt );


  areaflag=0;

  if ( (area1>0.0) && (area2>0.0) && (!cflagj) && (!cflagm) ){

    areaflag=1;

    tri[m][0]=c;
    tri[m][1]=d;
    tri[m][2]=a;

    tri[j][0]=a;
    tri[j][1]=b;
    tri[j][2]=c;
  }
  
  //printf("\n");
  return;
}

//----------------------------------------------------------------------------------------------------

int trimesh(int nn, int tdim, int nb, int nbs[], int ***bs, double x[], double y[], int tri[][3]){
  FILE*fp;
  int i,j,k,m,n,q,t, ntg;
  int nt, seed, nopt;
  int (*nbr)[3]; //[]x[3] array, contiguous
  int p0,p1,p2,p3;
  int nbrtemp0,nbrtemp1,nbrtemp2;
  int flag, *testflip;
  int areaflag;
  int b;
  double xo,yo, xmin, xmax, ymin, ymax;
  double *tempx, *tempy, xt, yt;
  double x1,x2,x3,y1,y2,y3;


  seed=0;
  /*Set outer limits*/
  xmin=10e20;
  ymin=10e20;
  xmax=10e-20;
  ymax=10e-20;
  for (n=0; n < nn; n++){
    xmin=MIN(xmin,x[n]);
    ymin=MIN(ymin,y[n]); 
    xmax=MAX(xmax,x[n]);
    ymax=MAX(ymax,y[n]);     
  }

  // Errored for all + pts
  xmin=1.2*xmin;
  ymin=1.2*ymin;
  xmax=1.2*xmax;
  ymax=1.2*ymax;

  ////printf("\nMin x,y = %e, %e",xmin,ymin);
  ////printf("\nMax x,y = %e, %e",xmax,ymax);


  // allocate temp coordinates
  tempx = new double[nn+4];
  tempy = new double[nn+4];
  
  for(n=0;n<nn;n++){
    tempx[n+4] = x[n];
    tempy[n+4] = y[n];
  }

  // Create the superpoints
  // Put them at the start
  tempx[0]=xmin; tempx[1]=xmax; tempx[2]=xmax; tempx[3]=xmin;
  tempy[0]=ymin; tempy[1]=ymin; tempy[2]=ymax; tempy[3]=ymax;


  // Create the supernodes
  // Put them at the start
  tri[0][0]=0;
  tri[0][1]=1;
  tri[0][2]=2;
  tri[1][0]=0;
  tri[1][1]=2;
  tri[1][2]=3;
  //ntg=2;
  nt=2;  //index of tri's shall be nt-1

  // for (t=0; t < tdim; t++){  
  //   //printf("\n %d, tri[t][0]=%d, tri[t][1]=%d, tri[t][2]=%d",t, tri[t][0], tri[t][1],tri[t][2]);
  // }


  // make neighbors using your routine
  nbr = new int[tdim][3];
  make_nbrs(nn+4,tdim,tri,nbr);
  
  /* print out neighbors for each triangle */
  for (t=0; t < nt; t++){
    //printf("\nElement %d, has neighbors %d, %d, %d",t,nbr[t][0],nbr[t][1],nbr[t][2]);
    //printf("\n");
  }






  // ===================================================
  // store the boundary triangles you find along the way.
  //int ib, jb;
  //List *boundary = new List;


  //======================
  /*Start inserting Points*/

  for (n=4; n<(nn+4); n++){
    //printf("\n|-----------------------|");
    //printf("\n| Inserting Point %d     ", n);
    //printf("\n|                       |");

    if (nt>nn*3){
      nt=-1;
      break;
      
    }

    
    // Initialize List of Triangles to Optimize
    /*Create Optimize List - of Triangles to be tested as we see fit (ergo Lawson algo)*/
    List *Optimize = new List;
    
    xt=tempx[n];
    yt=tempy[n];

    //printf("\n| x,y = %g,%g       ", xt,yt);
    //printf("\n|-----------------------|");   

 
    // Find Triangle Containing Point
    // get the seed triangle, stores last value
    seed = search(seed,xt,yt,tempx,tempy,tri,nbr);
    //printf("\nThe                     ");
    //printf("\n    Search              ");
    //printf("\n           Has          ");
    //printf("\n               Returned!");
    //printf("\nWe found a pt in triangle %d",seed);  ////////////////////////////Check this triangle for sense making


    // save nodes and neighbors of old triangle
    p0=tri[seed][0];
    p1=tri[seed][1];
    p2=tri[seed][2];  
    nbrtemp0=nbr[seed][0];
    nbrtemp1=nbr[seed][1];
    nbrtemp2=nbr[seed][2];

    // Add new point
    // move J node at P2 to n for efficiency -no leave it be for easy.
    p3=n;
    
    // Subdivide Triangle
    // Add 4 points for 3 triangles to the optimize list
    //1
    //printf("\nReassign 1st Triangle");
    tri[seed][0]=p0;
    tri[seed][1]=p1;
    tri[seed][2]=p3;
    Optimize -> Add_To_List(seed); //pop seed triangle onto optimimize list
    //2
    //printf("\nAssign 2nd Triangle");
    tri[nt][0]=p1;
    tri[nt][1]=p2;
    tri[nt][2]=p3;
    //printf("\n!-----------------------------------------------");
    //printf("\n!nt at the time when tri[2][] is created = %d   ",nt);
    //printf("\n!-----------------------------------------------");
    Optimize -> Add_To_List(nt);
    nt=nt+1;
    //3
    //printf("\nAssign 3rd Triangle");
     tri[nt][0]=p2;
     tri[nt][1]=p0;
     tri[nt][2]=p3;
     Optimize -> Add_To_List(nt);
     nt++; 
     
     // run make neighbors - this will be sllllooooowww
     //make_nbrs((n+1),nt,tri,nbr);
     make_nbrs(nn+4,tdim,tri,nbr);
     









     /**/
     // Optimize Triangles in the list
     nopt=3;  //gauranteed # of triangles at the start of each run
     
     // What do I do to nopt if I add a triangle to the optimize list later??




  
     // ==================================
     // Begin the Loop Over Optimize Points
     // Begin the Loop Over Optimize Points 
     // Begin the Loop Over Optimize Points 
     // Begin the Loop Over Optimize Points
     // Begin the Loop Over Optimize Points     
     // Begin the Loop Over Optimize Points
     //
     for(i=0; i<nopt; i++){ 
       //printf("\n opt # i = %d",i);
       j=Optimize->list[i];    //Grab the ith triangle in the optimize list and put it on j
       //printf("\n\nOptimize------------------------------------Optimize");
       //printf("\n       Starting Optimize of triangle %d ", j);
       //printf("\nOptimize------------------------------------Optimize\n");
       //use j to peal the correct tri[j][] values
       
       // If we flip an edge, come back to this point - or somewhere near here.
       // Add the new triangle to the optimize loop
       // increase the loop as needed
       
       // check seed, run circletest on the seed triangle against seed's 3 neighbors
       for (k=0;k<3; k++){
	 m=nbr[j][k];  //tri[j][k] of the k'th neighbor of j

	 // If triangle is on the edge,
	 if (m==-1){
	 }

	 else if (m!=-1){
	   
	   // If those neighbors' outnode points pass the circle test leave it alone

	   // If it fails circle test - 
	   //      put in test array and flip.
	   //      Check the area
	   //      Do circle test
	   //      If either of the two neighbors fails either test abort flip
	   //      If it passes accept the flip and pass into real container

	   // What I do:
	   // If the neighbor exists
	   // call unshared node function "outnode"
	   // Find the neighbor node that isn't a seed node.
	   // usage: (seed, neighbor, tri-array)

	   // return: the node in triangle m that isn't in seed triangle j
	   // This will be my test point
	   q=outnode(j,m,tri);
	   
	   // Get the circle test points
	   //1st get the seed points, i.e. the j's
	   x1=tempx[tri[j][0]];
	   x2=tempx[tri[j][1]];
	   x3=tempx[tri[j][2]];	  
	   y1=tempy[tri[j][0]];
	   y2=tempy[tri[j][1]];
	   y3=tempy[tri[j][2]];
	   //then get the test point - unshared m
	   xt=tempx[tri[m][q]];
	   yt=tempy[tri[m][q]];
	   //flag=circle_test(x1,y1,x2,y2,x3,y3,xt,yt);
     double tp1[2] = {x1,y1};
     double tp2[2] = {x2,y2};
     double tp3[2] = {x3,y3};
     double tpt[2] = {xt,yt};
     flag=circle_test( tp1,tp2,tp3,tpt );
	   if (!flag){
	     //nt++;
	     //printf("\nTriangle # %d, with nodes tri[j][0] = %d, tri[j][1] = %d, tri[j][2] = %d",j,tri[j][0], tri[j][1],tri[j][2] );
	     //printf("\nPassed the Circle Test with neighbor m= %d, with nodes: tri[m][0] = %d, tri[m][1] = %d, tri[m][2] = %d ", m,tri[m][0], tri[m][1],tri[m][2] );
	     //printf("\n  (%g, %g)",x1,y1);
	     //printf("\n  (%g, %g)",x2,y2);
	     //printf("\n  (%g, %g)\n",x3,y3);
	     //printf("\n test point (%g, %g)\n",xt,yt);
	     //continue; //Jump to front of the loop 
	     // GDB says this - goes to line ~ 340 - the beggining of the circletest neighbor loop - Terrible- Must go to k Loop!
	     //Checked a neighbor
	   }
	   
	   /*Edge Flip*/
	   else if(flag){
	     //printf("\nTriangle j = %d Nodes, %d,%d,%d,",j,tri[j][0], tri[j][1],tri[j][2]);
	     //printf("\nfails the circle test with neighbor m = %d, nodes %d, %d, %d\n",m,tri[m][0],tri[m][1],tri[m][2] );
	     //printf("\n\nInvoke The Edge Flipper");
	     // initialize new array 
	     testflip = new int[6];
	     // //============================================================ 
	     // fp=fopen("Insert3.dat","w");
	     // for (i=0; i < nt; i++)
	     //   {
	     // 	 int n0 = tri[i][0];
	     // 	 int n1 = tri[i][1];
	     // 	 int n2 = tri[i][2];
	     // 	 fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n0],tempy[n0]);
	     // 	 fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n1],tempy[n1]);
	     // 	 fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n2],tempy[n2]);
	     // 	 fprintf(fp,"%19.10e %19.10e 0.0\n\n",tempx[n0],tempy[n0]);
	     //   }
	     // fclose(fp);
	     // //============================================================
	     
	     edge_flip(j,m,testflip,tri,nbr,tempx,tempy, areaflag);
	     // //============================================================ 
	     // fp=fopen("Insert4.dat","w");
	     // for (i=0; i < nt; i++)
	     //   {
	     // 	 int n0 = tri[i][0];
	     // 	 int n1 = tri[i][1];
	     // 	 int n2 = tri[i][2];
	     // 	 f//printf(fp,"%19.10e %19.10e 0.0\n",  tempx[n0],tempy[n0]);
	     // 	 f//printf(fp,"%19.10e %19.10e 0.0\n",  tempx[n1],tempy[n1]);
	     // 	 f//printf(fp,"%19.10e %19.10e 0.0\n",  tempx[n2],tempy[n2]);
	     // 	 f//printf(fp,"%19.10e %19.10e 0.0\n\n",tempx[n0],tempy[n0]);
	     //   }
	     // fclose(fp);
	     // //============================================================
	     
	     if(areaflag==1){
		 //printf("\nFlip Connectivity Accepted - Yay!");
		 make_nbrs(nn+4,tdim,tri,nbr);
		 Optimize -> Add_To_List(j);  // Flipped Tri retains its j number as the same old real number.
		 Optimize -> Add_To_List(m);
		 nopt=nopt+2;
		 // Get out of the triangle
		 k=3;
              
		 //continue;
		 //and increment opt by 1??

	     }
	     else if(areaflag==0){
	       //printf("\nTri j=%d, Neighbor m= %d Cannot Satisfy Area with Flip - Aborting", j,m);
	       //m=nbr[m][k];
	      	       
	     }	     
	   }	   
	 } // if m!=-1 statement checks to see that neighbor is +




       }// k loop - circletest 3 neighbors of i loop

       // loop complete, pop j off the optimize list
       //Optimize -> Delete_From_List(nt);
       //Optimize -> Delete_From_List(j);
       // //printf("\n------------------------------------");
       // //printf("\nTriangle %d", j);
       // //printf("\nHas been Deleted from the Optimize List");
       // //printf("\n------------------------------------");

     }// i loop to optimize every j on the list

     //printf("\n|===================================");
     //printf("\n| We have gone through the entire optimize loop!");
     //printf("\n| Time To add More Triangles!");
     free(Optimize);


     printf("\n n = %d, nt = %d",n, nt);
  } // n Loop to insert all points
  

  /**/
  printf("\n======================================================");
  printf("\n==========Starting Boundary Reconstruction===========");
  printf("\n======================================================");
  //Short and slow:
  //Just loop over every triangle containing the first node
  //If you've trapped the ray,
  //flip the triangle.
  //not hard.
  int z; //=tnum # of triangles
  int b1,b2;
  int currenttriangle;
  int numtriangles;
  int boundaryflag, flipcount;
  int bnode0, bnode1, bnode2, bnode3;//actually triangle nodes, triangles which contain a boundary line.
  //x,y's for the vectors
  double xb1,xb2,xb3,yb1,yb2,yb3, xb4,yb4;
  double v1x, v1y, v2x, v2y, bxray, byray;
  flipcount=0;
  // Boundary Restoration
  // Create Hash Table.
  List **nhash;
  nhash = new List*[nn+4];
  for(n=0;n<nn+4;n++){
    nhash[n]=new List();
  }
  
  // initial loops
  // Add all points to the list
  for (t=0;t<nt;t++){//loop all triangles
    for (i=0;i<3;i++){//loop all nodes
      n=tri[t][i];
      nhash[n]->Add_To_List(t);
    }
  }
  
  // //Loop over # of boundaries
  for (i=0;i<nb;i++){   //i<nb
    printf("\n----------------------------");
    printf("\nStart of the boundaries loop");
    printf("\n");
    
  //   //Loop over segments of each boundary
    for (j=0;j<nbs[i];j++){
      printf("\n----------------------------");
      printf("\nStart of the loop over the segments of the boundary");
      printf("\n");
      b1=bs[i][j][0]+4;
      b2=bs[i][j][1]+4;
      numtriangles = nhash[b1]->max;

      printf("\nThe number of triangles to be checked = %d",numtriangles);
  //     //now loop over those triangles
  //     for (z=0;z<numtriangles;z++){  
  // 	printf("\nLoop over the Triangles from list Loop to check for exact fit");
  // 	currenttriangle=nhash[b1]->list[z];       //  get current tri, index = z
  // 	printf("\nThe Current triangle to be checked for exact fit is = %d",currenttriangle);

  // 	// loop over each node in the triangles
  //boundaryflag=0;
  // 	//loop over each node
  // 	for (m=0;m<3;m++){
  // 	  n=tri[currenttriangle][m];   //tri shows up here, getting the node.
  // 	  printf("\nLoop over each node in the triangle.  current node = %d",n);
  // 	  printf("\nStart the Node flip loop b1 = %d, b2 = %d, n= %d",b1,b2,n );


  // 	  // The equality
  // 	  if(b2==n){                    //2 found the node with the boundary.
  // 	    boundaryflag=1;              //boundary is at node n
  // 	    printf("\nBoundarydflag = %d", boundaryflag);
  // 	    //printf("\n BoundaryFlag is now 1.  b1 = %d, n= %d",b1,b2,n );
  // 	    printf("\nExit without Checking or Fliping the Edge");
  // 	    printf("\nBecause the Triangle edge is aligned with the boundary");
  // 	    m=3;//kick out of the loop for this triangle only.  It is fine.
  // 	    //z=numtriangles;
  // 	    printf("\n");
  // 	  }
	  //}	  
	  //}
      
   
      printf("\n===================== On to the test");
      printf("\n");

      
      //End of easy node finding (actually cut out)
      //while(boundaryflag==0){
      //printf("\nBoundaryFlag==0, procede to DotProducts");
      for (z=0;z<numtriangles;z++){
	printf("\n----------------------------");
	printf("\nLoop over the triangles");
	printf("\n");
	  
	currenttriangle=nhash[b1]->list[z];   //get current test tri, lindex = z
	    
	printf("\nThe current test triangle is = %d\n",currenttriangle);
	    
	for (m=0;m<3;m++){
	  printf("\n----------------------------");
	  printf("\nLoop over the nodes");
	  printf("\n");
	  printf("\nThe current test triangle =%d , tri node m = %d",currenttriangle, m);
	  n=tri[currenttriangle][m];               // n=the node
	  printf("\n Start the test loop b1 = %d, b2 = %d, n= %d",b1,b2,n );

	    // The ray origin
	    if(b1==n){
	      bnode0=m;        
	      printf("\nbnode = %d",m);

	      //get the other points on the tri
	      bnode1=(bnode0+1)%3;
	      bnode2=(bnode0+2)%3;
	    
	      //Get the x,y points on the triangle in question
	      xb1=tempx[tri[currenttriangle][bnode0] ];
	      xb2=tempx[tri[currenttriangle][(bnode0+1)%3]];
	      xb3=tempx[tri[currenttriangle][(bnode0+2)%3]];
	    

	      xb4=tempx[nbr[currenttriangle][(bnode0+1)%3]];
	      //xb4=tempx[tri[currenttriangle][(b2)%3]];
	      
	      printf("\nxb1 = %g, xb2 = %g, xb3 = %g, xb4 = %g", xb1, xb2, xb3, xb4);
	      
	      yb1=tempy[tri[currenttriangle][bnode0] ];
	      yb2=tempy[tri[currenttriangle][(bnode0+1)%3]];
	      yb3=tempy[tri[currenttriangle][(bnode0+2)%3]];
	      
	      yb4=tempy[nbr[currenttriangle][(bnode0+1)%3]];
	      //yb4=tempy[tri[currenttriangle][(b2)%3]];
	      
	      printf("\nyb1 = %g, xy2 = %g, xy3 = %g, xy4 = %g", yb1, yb2, yb3, yb4);
	      
	      v1x=xb2-xb1;
	      v2x=xb3-xb1;
	      v1y=yb2-yb1;
	      v2y=yb3-yb1;
	      
	      bxray=xb4-xb1;
	      byray=yb4-yb1;
	      
	      //get normals facing in!
	      //first -second
	      double vmag, xnormal1, ynormal1,xnormal2, ynormal2;
	      vmag=mag(xb2,yb2,xb1,yb1);
	      xnormal1 = -v1x/vmag;
	      ynormal1 = v1y/vmag;
	      vmag=mag(xb3,yb3,xb1,yb1);
	      xnormal2 = -v2x/vmag;
	      ynormal2 = v2y/vmag;
	      
	      
	      
	      
	      //Compute Dot*Products
	      double dp1, dp2;
	      dp1=dotproduct(xnormal1, ynormal1, bxray, byray);
	      dp2=dotproduct(xnormal2, ynormal2, bxray, byray);	
	      
	      printf("\nDot products 1 and 2 = %g, %g ", dp1,dp2);
	      
	      int neighbornode, rayneighbor;
	      //Compute Needed Neighbor	   
	      
	      
	      
	      
	      areaflag=0;
	      if ( (dp1>0.-10e-14) && (dp2>0.-10e-15) ){
		printf("\nSucess!");
		rayneighbor=nbr[currenttriangle][(bnode0+1)%3];
		if (rayneighbor!=-1){
		  printf("\nCalling the boundary flipper");
		  boundary_flip(currenttriangle,rayneighbor,testflip,tri,nbr,tempx,tempy, areaflag);
		  flipcount=flipcount+1;
		  printf("\nBoundaryFlipCount = %d",flipcount);
		  
		  make_nbrs(nn+4,tdim,tri,nbr);
		  //make_nbrs(nb[],[0],(q+1)%3);		  
		}
	      }
	      else {
		printf("\nDid not trap the ray.  Procede to the next triangle.");
		//boundaryflag=1;
	      }
	    }    
	}
	//boundaryflag=1;
	
	
	
	}
    }
  }
  /**/
  //End of bad boundary logic - and the end of my time to fix it.

  for(n=0;n<nn;n++){
    x[n]= tempx[n+4];
    y[n]= tempy[n+4];
  }
   /*If we run out of space return -1 (logic at top of loop)*/
   //printf("\nnt= %d",nt);


   fp=fopen("output/tlmGNUplot.dat","w");
   for (i=0; i < nt; i++)
     {
       int n0 = tri[i][0];
       int n1 = tri[i][1];
       int n2 = tri[i][2];
       fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n0],tempy[n0]);
       fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n1],tempy[n1]);
       fprintf(fp,"%19.10e %19.10e 0.0\n",  tempx[n2],tempy[n2]);
       fprintf(fp,"%19.10e %19.10e 0.0\n\n",tempx[n0],tempy[n0]);
     }
   fclose(fp);



   delete[] tempx;
   delete[] tempy;



   return (nt);
}  

