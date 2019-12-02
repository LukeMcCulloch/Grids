#include <stdio.h>
#include "List.h"

/* void make_hash(int nn, int ntri, int **tri, int **nbr){ */
/*   int t, n, i, j; */
/*   List ** nhash; */
/*   nhash = new List *[nn]; */
/*   for (n=0; n < nn; n++){ */
/*     /\*Make a list of lists*\/ */
/*     nhash[n] = new List(); */
/*   } */
/*   for (t=0; t<ntri; t++){ */
/*     for (i=0,i<3; i++){ */
/*       n=tri[t][i]; */
/*       nhash[n]->Add_To_List(t); */
/*     } */
/*   } */
/*   for (t=0; t<ntri; t++){ */
/*     for (i=0,i<3; i++){ */
/*       n0=tri[t][i]; */
/*       n1=tri[t][(i+1)%3]; //nifty trick */
/*       for (j=0,j<nhash[n0]->max; j++){ */
/* 	a=nhash[n0]->List[j]; */
/* 	if (a=t=){ */
/* 	  continue; */
/* 	} */
/* 	if (nhash[n1]->Is_In_list(a)){ */
/* 	  nbr[t][i]=a; */
/* 	break; */
/* 	} */
/*       } */
/*     } */
/*   } */
/* } */




void make_nbrs(int nn, int ntri, int tri[][3], int nbr[][3]){
  int t, n, i, j;
  int n0, n1, a;
  int check;

  printf("\nArrived in make_nbrs");
  /*Make Hash Table*/
  printf("\nMake Hash Table");
  List ** nhash;
  nhash = new List *[nn];
  for (n=0; n < nn; n++){
    printf("\n new list N %d",n);
    /*Make a list of lists*/
    nhash[n] = new List();
  }
  printf("\n ntri = %d",ntri);
  for (t=0; t<ntri; t++){
    for (i=0; i<3; i++){
      //printf("\n nhash n =Tri %d, %d",t,i);
      n=tri[t][i];
      printf("\n nhash %d =Tri %d, %d",n,t,i);
      check = nhash[n]->Add_To_List(t);
      if(check == 0){
        printf("failed to add to list");
      }
    }
  }

  /*Make Neighbors*/
  printf("\nMake Neighbors");
  for (t=0; t<ntri; t++){
    printf("\nTri %d",t);
    for (i=0; i<3; i++){
      printf("\nNode %d",i);
      nbr[t][i]=-1;
      n0=tri[t][i];
      n1=tri[t][(i+1)%3]; //nifty trick
      for (j=0; j<nhash[n0]->max; j++){
        a=nhash[n0]->list[j];
        if (a==t){
          continue;
        }
        if (nhash[n1]->Is_In_List(a) ){
          nbr[t][i]=a;
          break;
        }
      }
    }
  }
  printf("\nDone with make Neighbors");
  return;
}
