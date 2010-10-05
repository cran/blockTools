#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "util.h"

void mahalnaive2(double *data, int *nrow, int *ncol, double *vcovi,  int *ntr, int *l2, int *l1names, int *valid, double *validvar, double *validlb, double *validub, int *verbose, double *pairdist, int *result){
  int n = choose(*nrow, 2);
  int j, mn[*ntr - 1], i, k=0, ii, cinf=n, t, star[*ntr], kk;
  double min, *vec = calloc(n, sizeof(double)), X[*nrow][*ncol], tmp1[*ncol], tmp2[*ncol], r;
  for(j=0; j<*ncol; j++){
    for(i=0; i<*nrow; i++){
      X[i][j] = data[k];
      k++;
	}
  }
  k=0;
  for(i=1; i<*nrow; i++){
    for(j=0; j<i; j++){
      for(ii=0; ii<*ncol; ii++){
	tmp1[ii] = X[i][ii];
	tmp2[ii] = X[j][ii];
      }
      vec[k] = cmahal(tmp1, vcovi, ncol, tmp2);
      k++;
    }
  }
  if(*l2 == 1){
    for(i=0; i<n; i++){
      j = ceil(sqrt((2*(i+1))+.25)+0.5);
      if(l1names[j-1] == l1names[((i+1)-((j-2)*(j-1)/2))-1]){
	vec[i] = HUGE_VAL;
      }
    }
  }
  if(*valid == 1){
    for(i=0; i<n; i++){
      j = ceil(sqrt((2*(i+1))+.25)+0.5);
      if(sqrt(pow((validvar[j-1] - validvar[((i+1)-((j-2)*(j-1)/2))-1]), 2)) < *validlb && sqrt(pow((validvar[j-1] - validvar[((i+1)-((j-2)*(j-1)/2))-1]), 2)) > *validub){
	vec[i] = HUGE_VAL;
      }
    }
  }
  j = -1;
  for(i=0; i<*nrow; i++){
    j++;
    cinf = 0;
    for(ii = tri((i+1) - 2); ii < tri((i+1) - 1); ii++){
      if(vec[ii] < HUGE_VAL){
	cinf++;
      }
    }
    for(ii = ((i+1) + 1); ii <= *nrow; ii++){
      if(vec[(tri(ii - 2) + (i+1)) - 1] < HUGE_VAL){
	cinf++;
      }
    }

    if(cinf == 0){ 
      j--;
      continue;
    }
    
    if(*verbose == 1){
      Rprintf("Block number is %i \n", j+1);
    }
    mn[0] = tri((i+1) - 2) + 1;
    min = vec[tri((i+1) - 2)];
    t=0;
    for(ii = tri((i+1) - 2)+1; ii < tri((i+1) - 1); ii++){
      if(vec[ii] < min){
	t = 0;
	min = vec[ii];
	mn[0] = ii + 1;
      }
      else if (vec[ii] == min){
	t++;
      }
    }
    if(t > 0){
      for(ii = tri((i+1) - 2); ii < tri((i+1) - 1); ii++){
	if(vec[ii] == min){
	  GetRNGstate();
	  r = unif_rand();
	  PutRNGstate(); 
	  if(r < 1/t){
	    mn[0] = ii+1;
	  }
	}
      }
    }
    t=0;
    for(ii = ((i+1) + 1); ii <= *nrow; ii++){
      if(vec[(tri(ii - 2) + (i+1)) - 1] < min){
	t=0;
	min = vec[(tri(ii - 2) + (i+1)) - 1];
	mn[0] = (tri(ii - 2) + (i+1));
      }
      else if (vec[(tri(ii - 2) + (i+1)) - 1] == min){
	t++;
      }
    }
    if(t > 0){
      for(ii = ((i+1) + 1); ii <= *nrow; ii++){
	if(vec[(tri(ii - 2) + (i+1)) - 1] == min){
	  GetRNGstate();
	  r = unif_rand();
	  PutRNGstate(); 
	  if(r < 1/t){
	    mn[0] = (tri(ii - 2) + (i+1));
	  }
	}
      }
    }
    star[0] = ceil(sqrt(2*mn[0]+.25)+0.5);
    star[1] =  (mn[0]-((star[0]-2)*(star[0]-1)/2));
    pairdist[j] = vec[mn[0]-1];
    vec[mn[0] - 1] = HUGE_VAL;
    if(*l2 == 1){
      for(kk=0; kk<n; kk++){
	ii = ceil(sqrt((2*(kk+1))+.25)+0.5);
	if(l1names[ii-1] == l1names[star[0]-1] && (ii != star[0])){
	  vec[kk] = HUGE_VAL;
	} 
	else if(l1names[ii-1] == l1names[star[1]-1] && (ii != star[1])){
	 vec[kk] = HUGE_VAL;
	}
	else if(l1names[((kk+1)-((ii-2)*(ii-1)/2))-1] == l1names[star[0]-1] && ((kk+1)-((ii-2)*(ii-1)/2)) != star[0]){
	  vec[kk] = HUGE_VAL;
	  } 
	else if(l1names[((kk+1)-((ii-2)*(ii-1)/2))-1] == l1names[star[1]-1] && ((kk+1)-((ii-2)*(ii-1)/2)) != star[1]){
	 vec[kk] = HUGE_VAL;
	}
      }
    }  
    k=2;
    while(k < *ntr){
      t = 0;
      min = vec[tri(star[0] - 2)];
      mn[k-1] = tri(star[0] - 2) + 1;
      for(kk=0; kk<k; kk++){
	for(ii = tri(star[kk] - 2); ii < tri(star[kk] - 1); ii++){
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
	 for(ii = tri(star[kk] - 2); ii < tri(star[kk] - 1); ii++){
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
	for(ii = (star[kk] + 1); ii <= *nrow; ii++){
	  if(vec[(tri(ii - 2) + star[kk]) - 1] < min){
	    t=0;
	    min = vec[(tri(ii - 2) + star[kk]) - 1];
	    mn[k-1] = (tri(ii - 2) + star[kk]);
	  }
	  else if (vec[(tri(ii - 2) + star[kk]) - 1] == min){
	    t++;
	  }
	}
	if(t > 0){
	  for(ii = (star[kk] + 1); ii <= *nrow; ii++){
	    if(vec[(tri(ii - 2) + star[kk]) - 1] == min){
	     	  GetRNGstate();
		  r = unif_rand();
		  PutRNGstate(); 
		  if(r < 1/t){
		    mn[k-1] = (tri(ii - 2) + star[kk]);
		  }
	    }
	  }
	 }
      }
      t=0;
      for(kk=0; kk<k; kk++){
        if(ceil(sqrt(2*mn[k-1]+.25)+0.5) == star[kk]){
	 star[k] = (mn[k-1]-((star[kk]-2)*(star[kk]-1)/2));
	 t++;
	}
      }
      if(t==0){
	 star[k] = ceil(sqrt(2*mn[k-1]+.25)+0.5);
      }
     
      for(kk=0; kk<k; kk++){
        if(star[kk] < star[k] && vec[(tri(star[k]-2) + star[kk]) - 1] > pairdist[j] && vec[(tri(star[k]-2) + star[kk]) - 1] < HUGE_VAL){
	  pairdist[j] = vec[(tri(star[k]-2) + star[kk]) - 1];
	}
        if(star[kk] > star[k] && vec[(tri(star[kk]-2) + star[k]) - 1] > pairdist[j] && vec[(tri(star[kk]-2) + star[k]) - 1] < HUGE_VAL){
	  pairdist[j] = vec[(tri(star[kk]-2) + star[k]) - 1];
	}
      }
      if(*l2 == 1){
	for(kk=0; kk<n; kk++){
	  ii = ceil(sqrt(2*(kk+1)+.25)+0.5);
	  if(l1names[((kk+1)-((ii-2)*(ii-1)/2)) - 1] == l1names[star[k] - 1] && ((kk+1)-((ii-2)*(ii-1)/2)) != star[k]){
	    vec[kk] = HUGE_VAL;
	  }
	  else if(l1names[ii - 1] == l1names[star[k] - 1] && ii != star[k]){
	    vec[kk] = HUGE_VAL;
	  }
	}
      }
      if(vec[mn[k-1]-1] == HUGE_VAL){
	star[k] = 0;
      }
      vec[mn[k-1] - 1] = HUGE_VAL;
      k++;
    }
    for(kk=0; kk<*ntr; kk++){
       result[j*(*ntr) + kk] =  star[kk];
    }
    for(kk=0; kk<*ntr; kk++){
       if(star[kk] != 0){
	 for(ii = tri(star[kk] - 2); ii < tri(star[kk] - 1); ii++){
	   vec[ii] = HUGE_VAL;
	 }
	 for(ii= (star[kk] + 1); ii <= *nrow; ii++){
	   vec[(tri(ii-2) + star[kk]) - 1] = HUGE_VAL;
	 }
       }
    }
  }
}
    
 
