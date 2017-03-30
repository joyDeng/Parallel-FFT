/*cilk program for computing the FFT.

   main function tests by making a spike vector, calling FFT, and
   checking the result.  Args to main function are lg of the problem
   size, followed by an optional spike position (default is 1).

   To compile as the C elision, use the C flag -DNOCILK

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef NOCILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#ifdef NOCILK
#define cilk
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#endif

#include "fft.h"
#include "timer.h"


void Bit_Reverse(complex_t *, int);
/* Compute the FFT of an array a of length N, given lgN. */
void FFT(complex_t *a, int lgN) {
  unsigned int n = (unsigned int)pow(2,lgN);
 // unsigned int h = (unsigned int)pow(2,lgN-1);
  

 // double *sine = (double *)calloc(h,sizeof(double));
  //double *cosi = (double *)calloc(h,sizeof(double));

  
  // cilk_for(int x = 0 ; x < h ; x++){
  //   sine[x] = sin(M_2PI * (x)/(n));
  //   cosi[x] = cos(M_2PI * (x)/(n));
  //  // printf("Index:%d: sin:%f, cos:%f \n", x, sine[x], cosi[x]);
  // }

  Bit_Reverse(a,lgN);

  
  
  for(int l = 0 ; l < lgN ; l++){
    unsigned int m = pow(2,l+1);
    //unsigned int b = pow(2, lgN - l - 1);
    cilk_for(unsigned int g = 0 ; g < n ; g += m){
	    cilk_for( unsigned int k = 0 ; k < m/2 ; k++){
      		complex_t t,u,w;
      //    int index = b*k;
      //    (w).imag = sine[index];
      //    (w).real = cosi[index];


      		COMPLEX_ROOT(w,m,k);
  
      		COMPLEX_MULT(t,w,a[g+k+m/2]);
      		u = a[g+k];
      	    COMPLEX_ADD(a[g+k],u,t);
      		COMPLEX_SUB(a[g+k+m/2],u,t);
	    }

   }
  }

}

void Bit_Reverse(complex_t *a, int lgN){
 
  unsigned int max_index = (unsigned int)pow(2,lgN);

 cilk_for(unsigned int i = 0 ; i < max_index ; i++){
        unsigned int num = i;
	unsigned int reverse = 0;
	int count = lgN;
	  while(num){
	    reverse <<= 1;
	    reverse |= (num  & 1);
 	    num >>= 1;
	    count--;
	  }
	  reverse <<= count;
	  if(reverse > i){
	    complex_t temp;
	    temp= a[reverse];
	    a[reverse] = a[i];
	    a[i] = temp;
	  }
	
  }
}

