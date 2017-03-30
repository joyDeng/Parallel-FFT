#include "fft.h"
#include "bit_matrix_fns.cpp"
#include "bit_matrix_types.h"
#include "bmmc_mpi.h"
#include "math.h"


void FFT(complex_t *a, int n, int p, int my_rank){
    //each process do bit reverse first;
    int depth = n - p;
    long N = 1<<n;
    long P = 1<<p;
    int entry_n = N/P;
    int r = n%(n-p);
    int superlevel = n/(n-p);
    int i;
    int super,n_mini;//exp of twiddle factor;
    //int t_sub = 1 << (n-1);
    //int limit = n-1;
    complex_t t,u,w;
    //printf("lllllllll");
    complex_t *temp1;
    complex_t *temp2;
    complex_t *wlist = new complex_t[1<<n];
  
    
    bit_matrix I = allocate_bit_matrix(n);
    bit_matrix r_m = allocate_bit_matrix(n);
    bit_matrix target_A = allocate_bit_matrix(n);
    bit_matrix target_B = allocate_bit_matrix(n);
    
    I = identity_matrix(I,n);//Identity matrix
    
    r_m = reverse_bit_matrix(r_m,I,n);//bit reverse matrix
    target_A = rotate_bit_matrix(target_A, I, n, n-p);//rotation matrix
    
    
    
    temp1 = new complex_t[entry_n];
    temp2 = new complex_t[entry_n];
    
    
    BMMC_MPI_proc_major(r_m,(matrix_column) 0,n,p,my_rank,MPI_COMM_WORLD,sizeof(complex_t),a,temp1);
    //    r_m = NULL;
    
    BMMC_MPI_factor_info *info;
    BMMC_MPI_factor_info *info2;
    info = new BMMC_MPI_factor_info();
    info2 = new BMMC_MPI_factor_info();
    
    factor_BMMC_MPI_proc_major(target_A,(matrix_column) 0,n,p,info);
    
    
    //if the last level is not full
    if(r != 0){
        target_B = rotate_bit_matrix(target_B, I, n, r);
        factor_BMMC_MPI_proc_major(target_B,(matrix_column) 0,n,p,info2);
        superlevel++;
    }
    
    
    
    for(i = 0 ; i < superlevel; i++){//for all the superlevel
        if(i == superlevel-1 && r != 0){//last superlevel
            //rotation for last level
            //n_mini is the number of minibutterfly pre process
            n_mini = (1 << n-p-r);
            
            for(int offset = 0 ; offset < n_mini; offset++){
                //offset: #of minibuterfly
                int begin = offset*(1 << r);
                for(int l = 0 ; l < r ; l++){
                   
                    long m = (1 << (l+1));
                     //precompute part of exponent
                    int super_left = (my_rank * n_mini + offset)/(1 << (n-i*(n-p)-r));
                    long super_right = pow(entry_n,i);
                    int sub = 1 << (i * depth + l + 1);
                    //int mod = 1 << (i*depth+l-t_sub);
                    for(int g = 0 ; g < (1 << r) ; g+=m){
                        for(int k = 0 ; k < m/2 ; k++){
                            super =  super_left+k*super_right;
                            COMPLEX_ROOT(w,sub,super);
                      
                            COMPLEX_MULT(t,w,a[begin+g+k+m/2]);
                            u = a[begin+g+k];
                            COMPLEX_ADD(a[begin+g+k],u,t);
                            COMPLEX_SUB(a[begin+g+k+m/2],u,t);
                        }
                    }
                }
            }
            perform_BMMC_MPI(info2,n,p,my_rank,MPI_COMM_WORLD,sizeof(complex_t),a,temp2);
            
            MPI_Barrier(MPI_COMM_WORLD);
            
        }else{//the level before last level
            for(int l = 0 ; l < depth ; l++){
                long m = (1 << (l+1));
                int super_left = my_rank/(1 << (n-i*(n-p)-depth));
                long super_right = pow(entry_n,i);
                int sub = 1 << (i * depth + l + 1);
                //int mod = 1 << (i*depth+l-t_sub);
                
                for(int k = 0 ; k < m/2 ; k++){
                    super = super_left + k*super_right;
                    COMPLEX_ROOT(wlist[k],sub,super);
                }
                
                for(int g = 0 ; g < entry_n ; g+=m){
                    //int begin = my_rank * entry_n;
                    for(int k = 0 ; k < m/2 ; k++){
                        
                        
                        //                        COMPLEX_ROOT(w,1 << (i * depth + l + 1)
                        //                                     ,super);
                        w = wlist[k];
                        COMPLEX_MULT(t,w,a[g+k+m/2]);
                        u = a[g+k];
                        COMPLEX_ADD(a[g+k],u,t);
                        COMPLEX_SUB(a[g+k+m/2],u,t);
                    }
                }
            }
            
            //perform right n-p bit permutation at the end of each super level
            perform_BMMC_MPI(info,n,p,my_rank,MPI_COMM_WORLD,sizeof(complex_t),a,temp2);
            MPI_Barrier(MPI_COMM_WORLD);
            
            
        }
        
    }
    
    //free allocated space
    if(temp1 != NULL)
        delete [] temp1;
    if(temp2 != NULL)
        delete [] temp2;
    free_bit_matrix(target_A);
    free_bit_matrix(target_B);
    free_bit_matrix(I);
    free_BMMC_MPI_factor_info(info);
    free_BMMC_MPI_factor_info(info2);
    free_bit_matrix(r_m);
    
}
