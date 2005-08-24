/* 
   Copyright (C) 2004 Van Andel Institute
   Written by Kyle A. Furge

   stats.c: BRAINLESS library for reb calculations

   WARNING: Several things may seem strange

   0 is 999
   NA is -999 or 0
   
   This was done for efficiency reasons.   
*/

#include <stdio.h>
#include <stdlib.h>

void pretty_print(double m[],  int nr,  int nc)
{
   int r,c;

  printf("======\n");
  for(r=0;r<nr;r++) {
    for(c=0;c<nc;c++) {
      printf("% .2f ",m[r+nr*c]);
    }
    printf("\n");
  }
}

void bcount(double *data, int start, int end, int* N, int *np) {
  int i,n,p;
  double d;
  
  for(i=start,n=0,p=0;i<=end;i++) {
    d = data[i];
    if(d != -999 & d != 0) n++;
    if(d > 0) p++;
  }
  *N=n;
  *np=p;
}

double bscore(int N, int np) {
  double Z = -999;   /* this is an NA */
  
  if (N > 1) {    /* fixme make this real */
    Z = (2*(double)np - (double)N)/sqrt((double)N);
    if(Z == 0) 
      Z = 999;       /*0 is NA */
  }
  return(Z);
}

double mean(double *v, int s, int e) 
{
  double sum = 0, c = 0;
  unsigned int i = 0;
  double d;
  
  for(c=0,i=s;i<=e;i++) {
    d = v[i];
    if (d != -999 & d != 0) {   /* these are NA */
      c++;
      if(d != 999)     /* this is a 0 remember */
	sum += v[i];
    }
  }
  if (c == 0) return(-999);
  return(sum/c);
}

void mx_col_mean(double *m, int nr, int nc, double *v) {
  int start,end,c;

  start = 0;
  end = nr-1;
  for(c=0;c<nc;c++) {
    v[c] = mean(m,start,end);
    start += nr;
    end += nr;
  }
}

void mov_binom_mx(double v[], int *len, int *k, double* im) {
  double stat;
  int start, end, nr, nc, c;
  int N = 0,np = 0;

  nc = *len;
  nr = *len - *k + 1;

  start = 0;
  end = *k -1;
  bcount(v,start,end,&N,&np);
  while(end < nc) {
    /*    stat = binom(v,start,end);  */
    stat = bscore(N,np);
    for(c=start;c<=end;c++) {
      im[start + nr*c] = stat;
    }
    if(v[start] != -999) N--;
    if(v[start] > 0) np--;
    end++;
    start++;
    if(v[end] != -999) N++;
    if(v[end] > 0) np++;
  }
}

void mov_binom_test(double v[], int *len, int *k, double* iv) {

  int nr,nc;
  double *im;

  nc = *len;
  nr = *len - *k + 1;

  im = (double *) calloc(nc*nr,sizeof(double));
  mov_binom_mx(v,len,k,im);
  mx_col_mean(im,nr,nc,iv);
  free(im);
}

void mspan_mov_binom_mx(double *v, int *vlen, int *range, int *rlen, double *mx)
{
  int r,c,nc,nr;
  double *im_vec, *im_mx;

  nr = *rlen;
  nc = *vlen;

  im_vec = (double *) calloc(nc,sizeof(double));
  for(r=0;r<nr;r++) {
    mov_binom_test(v,vlen,&range[r],im_vec);
    for(c=0;c<nc;c++) {
      mx[r+nr*c] = im_vec[c];
    }
  }
  free(im_vec);
}

void mspan_mov_binom(double *v, int *vlen, int *range, int *rlen, double *sv)
{
  int r,c,nc,nr;
  double *sum_mx;

  nr = *rlen;
  nc = *vlen;

  sum_mx = (double *) calloc(nr*nc,sizeof(double));
  mspan_mov_binom_mx(v,vlen,range,rlen,sum_mx);
  mx_col_mean(sum_mx,nr,nc,sv);
  free(sum_mx);
}

int main() {
  double data[] = {0.8,1.2,2,3,5,-0.3,-1,-1.1};
  double gx[] = {0.19,0,0,-0.61,-0.53,0.18,0.67,-0.06,-0.73,0.21,-0.74,-1.35,0};
  double gx2[] = {0.19,0,0,0,0,0,0,0,0,-0.61,-0.53,0.18,0,0.67,0,0,-0.06,0,0,-0.73,0.21,0,-0.74,-1.35,0,0,0,0,-0.95,-0.19,-0.75,0,0,0,0,0,-0.28,0.14,0,0,0.82,0,0.11,0.11,0,0.26,-0.15,-0.48,0.29,0,0.35,0.59,0.16,0,0,0.51,-0.14,0.16,-0.4,0.92,1.63,0,0,0,-0.46,-1.45,-0.92,0,0,-1.47,-0.31,0,0,-1.08,0,2.19,0,-0.42,0,-0.22,-0.09,0,0,0,0,0.64,0.18,0,0.04,0.49,0.08,0.15,0,0,-0.36,0,-0.46,0.15,0,0,-0.96,0.03,0,0,-0.28,0,0,0,0,0.95,0,0,0,0,0,-0.58,0,0,0,0.34,0,-0.36,0.16,0,0,0.02,0,0,0.92,0.07,0,-0.37,0.01,0.45,1.15,-0.84,-1.04,0,0,1.17,0.91,0,0.47,0,0,0,-0.11,0,0.79,0,0,0,0,0,-0.74,0,0,0,0,-0.21,0,0,-0.28,0,0,0,0,0,0.3,0,0,-1.73,-1.61,-1.13,-1.16,0,0,0,0,0,-0.66,-1.27,0,0,-2.07,0,-0.67,-0.71,0.82,0,0,-3.28,-2.72,-0.84,0,0,-0.78,-1.05,-1.22,-0.88,0,0,0,-1.65,0,1.37,-0.61,-1.17,-0.61,-0.81,0,0,-0.71,-0.53,-0.15,0,-0.32,0,-0.56,-0.17,0,-0.31,0,0.29,0,0.62,0,0,-0.29,-0.82,0,0,0,0,0.54,0.64,0,0,0,0,0,-0.67,-0.5,-0.75,-0.59,0,0.13,0,0,0,0.06,0,0,0,-0.57,0,0,0,-0.87,0.05,0.26,-0.5,-2.78,0,0,0,0,0,0.39,0,0,-0.06,0,0.06,0,-0.4};
  double *im_vec, *im_mx;
  double *sum_vec, *sum_mx;
  int range[] = {4,5,6,7};

  /* gx values */
  int len = 12;  int k = 7; int rlen = 4;

  int nr = len - k + 1;
  int nc = len;

  
  /*  printf("t: %.3f\n",ttest(gx,0,275)); */

  pretty_print(gx,1,nc);  

  /*
  im_vec = (double *) malloc(sizeof(double)*len); 
  im_mx = (double *) malloc(sizeof(double)*nr*nc);
  mov_binom_test(gx,&len,&k,im_mx,im_vec); 
  pretty_print(im_mx,nr,nc);
  pretty_print(im_vec,1,nc);
  free(im_vec);
  free(im_mx);

  sum_mx = (double *) malloc(sizeof(double)*len*rlen);
  mspan_mov_binom_mx(gx,&len,range,&rlen,sum_mx);
  pretty_print(sum_mx,rlen,len);
  free(sum_mx);
  */
}
