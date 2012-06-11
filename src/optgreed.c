#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "util.h"

void optgreed(double *vec, int *n, int *nrow, int *ntr, int *l2, int *l1names,  int *valid, double *validvar, double *validlb, double *validub, int *verbose, double *pairdist, int *result){
  int j, i, ii, cinf=*n, t, star[*ntr], mn[*ntr - 1], k;
  double min, r;
  if(*l2 == 1){
    for(i=0; i<*n; i++){
      j = ceil(sqrt((2*(i+1))+.25)+0.5);
      if(l1names[j-1] == l1names[((i+1)-((j-2)*(j-1)/2))-1]){
	vec[i] = HUGE_VAL;
      }
    }
  }
  if(*valid == 1){
    for(i=0; i<*n; i++){
      j = ceil(sqrt((2*(i+1))+.25)+0.5);
      if(sqrt(pow((validvar[j-1] - validvar[((i+1)-((j-2)*(j-1)/2))-1]), 2)) < *validlb || sqrt(pow((validvar[j-1] - validvar[((i+1)-((j-2)*(j-1)/2))-1]), 2)) > *validub){
	vec[i] = HUGE_VAL;
      }
    }
  }
 j=-1;
  while(cinf > 0){
  j++;
  if(*verbose == 1){
    Rprintf("Block number is %i \n", j+1);
  }
  mn[0] = 1;
  min = vec[0];
  t=0;
  for(i=0; i<*n; i++){
    if(vec[i] < min){
      t=0;
      min = vec[i];
      mn[0] = i+1;
    }
    else if(vec[i] == min){
      t++;
    }
  }
    if(t > 0){
      for(i=0; i<*n; i++){
	if(vec[i] == min){
	  GetRNGstate();
	  r = unif_rand();
	  PutRNGstate(); 
	  if(r < 1/t){
	    mn[0] = i+1;
	      }
	}
      }
    }
    star[0] = ceil(sqrt(2*mn[0]+.25)+0.5);
    star[1] =  (mn[0]-((star[0]-2)*(star[0]-1)/2));
    pairdist[j] = vec[mn[0]-1];
    vec[mn[0] - 1] = HUGE_VAL;
    for(i=0; i<*ntr; i++){
      result[j*(*ntr) + i] =  star[i];
    }
    if(*l2 == 1){
      for(i=0; i<*n; i++){
	ii = ceil(sqrt((2*(i+1))+.25)+0.5);
	if(l1names[ii-1] == l1names[star[0]-1] && (ii != star[0])){
	  vec[i] = HUGE_VAL;
	} 
	else if(l1names[ii-1] == l1names[star[1]-1] && (ii != star[1])){
	 vec[i] = HUGE_VAL;
	}
	else if(l1names[((i+1)-((ii-2)*(ii-1)/2))-1] == l1names[star[0]-1] && ((i+1)-((ii-2)*(ii-1)/2)) != star[0]){
	  vec[i] = HUGE_VAL;
	  } 
	else if(l1names[((i+1)-((ii-2)*(ii-1)/2))-1] == l1names[star[1]-1] && ((i+1)-((ii-2)*(ii-1)/2)) != star[1]){
	 vec[i] = HUGE_VAL;
	}
      }
    }
    k=2;
    while(k < *ntr){
      t = 0;
      min = vec[tri(star[0] - 2)];
      mn[k-1] = tri(star[0] - 2) + 1;
      for(i=0; i<k; i++){
	for(ii = tri(star[i] - 2); ii < tri(star[i] - 1); ii++){
	  if(vec[ii] < min){
	    t = 0;
	    min = vec[ii];
	    mn[k-1] = ii + 1;
	  }
	  else if (vec[ii] == min){
	    t++;
	  }
	}
	if(t > 0){
	 for(ii = tri(star[i] - 2); ii < tri(star[i] - 1); ii++){
	   if(vec[ii] == min){
	     	  GetRNGstate();
		  r = unif_rand();
		  PutRNGstate(); 
		  if(r < 1/t){
		    mn[k-1] = ii+1;
		  }
	   }
	 }
	}
	t=0;
	for(ii = (star[i] + 1); ii <= *nrow; ii++){
	  if(vec[(tri(ii - 2) + star[i]) - 1] < min){
	    t=0;
	    min = vec[(tri(ii - 2) + star[i]) - 1];
	    mn[k-1] = (tri(ii - 2) + star[i]);
	  }
	  else if (vec[(tri(ii - 2) + star[i]) - 1] == min){
	    t++;
	  }
	}
	if(t > 0){
	  for(ii = (star[i] + 1); ii <= *nrow; ii++){
	    if(vec[(tri(ii - 2) + star[i]) - 1] == min){
	     	  GetRNGstate();
		  r = unif_rand();
		  PutRNGstate(); 
		  if(r < 1/t){
		    mn[k-1] = (tri(ii - 2) + star[i]);
		  }
	    }
	  }
	}
      }
      t=0;
      for(i=0; i<k; i++){
	if(ceil(sqrt(2*mn[k-1]+.25)+0.5) == star[i]){
	  star[k] = (mn[k-1]-((star[i]-2)*(star[i]-1)/2));
	  t++;
	}
      }
      if(t==0){
	star[k] = ceil(sqrt(2*mn[k-1]+.25)+0.5);
      }
      for(i=0; i<k; i++){
	if(star[i] < star[k] && vec[(tri(star[k]-2) + star[i]) - 1] > pairdist[j] && vec[(tri(star[k]-2) + star[i]) - 1] < HUGE_VAL){
	  pairdist[j] = vec[(tri(star[k]-2) + star[i]) - 1];
	}
	if(star[i] > star[k] && vec[(tri(star[i]-2) + star[k]) - 1] > pairdist[j] && vec[(tri(star[i]-2) + star[k]) - 1] < HUGE_VAL){
	  pairdist[j] = vec[(tri(star[i]-2) + star[k]) - 1];
	}
      }
      if(*l2 == 1){
	for(i=0; i<*n; i++){
	  ii = ceil(sqrt(2*(i+1)+.25)+0.5);
	  if(l1names[((i+1)-((ii-2)*(ii-1)/2)) - 1] == l1names[star[k] - 1] && ((i+1)-((ii-2)*(ii-1)/2)) != star[k]){
	    vec[i] = HUGE_VAL;
	  }
	  else if(l1names[ii - 1] == l1names[star[k] - 1] && ii != star[k]){
	    vec[i] = HUGE_VAL;
	  }
	}
      }
      if(vec[mn[k-1]-1] == HUGE_VAL){
	star[k] = 0;
      }
      vec[mn[k-1] - 1] = HUGE_VAL;
      t=0;
      if(star[k] >0 ){
	for(i=0; i<k; i++){
	  if(star[i] == star[k]){
	    t++;
	  }
	}
	if(t == 0){
	  k++;
	}
      }
      else if(star[k] == 0){
	for(i = k; i<=*ntr; i++){
	  star[i] = 0;
	}
	k = *ntr;
      }
    }
    for(i=0; i<*ntr; i++){
      result[j*(*ntr) + i] =  star[i];
      if(star[i] != 0){
	for(ii = tri(star[i] - 2); ii < tri(star[i] - 1); ii++){
	  vec[ii] = HUGE_VAL;
	}
	for(ii= (star[i] + 1); ii <= *nrow; ii++){
	  vec[(tri(ii-2) + star[i]) - 1] = HUGE_VAL;
	}
      }
    }
    cinf = 0;
    for(i=0; i< *n; i++){
      if(vec[i] < HUGE_VAL){
	cinf++;
      }
    }
  }
}
