#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

int tri(int x){
  int out;
  out = ((x * x) + x)/2;
  return out;
}


double sum(double *x){
  int i, n= sizeof(*x)/sizeof(x[0]);
  double sum=0;
  for (i = 0; i < n; i++) sum += x[i];
  return sum;
}


double cmahal(double *x1, double *capX, int *n, double *x2){
  int i, j;
  double diff[*n], x1X[*n], ans=0;
  for(i=0; i<*n; i++){
    diff[i] = x1[i] - x2[i];
  }

  for(j=0;j < *n; j++){
    x1X[j]=0;
   for(i=0;i < *n;i++){
      x1X[j] += diff[i] * capX[i + j*(*n)];
    }
  }

  for(i=0; i<*n; i++){
    ans += x1X[i] * diff[i]; 
      }
  return sqrt(ans);
}
