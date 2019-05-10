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

  /*Make Hash Table*/
  List ** nhash;
  nhash = new List *[nn];
  for (n=0; n < nn; n++){
    /*Make a list of lists*/
    nhash[n] = new List();
  }
  for (t=0; t<ntri; t++){
    for (i=0; i<3; i++){
      n=tri[t][i];
      nhash[n]->Add_To_List(t);
    }
  }

  /*Make Neighbors*/
  for (t=0; t<ntri; t++){
    for (i=0; i<3; i++){
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
}
