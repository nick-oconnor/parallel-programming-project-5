/* 
 * File:   main.c
 * Author: fjosee
 *
 * Created on March 8, 2014, 1:15 PM
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define M_SIZE 16384

#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}

#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int result=0;
    unsigned long int upper, lower,tmp;
    __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
	);
    result = upper;
    result = result<<32;
    result = result|lower;
    return(result);
}

#endif

/***********************************************************************/
/* START: MT 19937******************************************************/
/***********************************************************************/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/***********************************************************************/
/* END: MT 19937 *******************************************************/
/***********************************************************************/

/* used for reference

 void matrix_multiply( double **A, double **B, double **C, int size )
{
    int i=0, j=0, k=0;
    unsigned long long last=rdtsc();
    unsigned long long current=rdtsc();

    for( i = 0; i < size; i++ )
      {
	if( i != 0 && i % 16 == 0 )
	  {
	    last = current;
	    current = rdtsc();
	  }
        for( j = 0; j < size; j++ )
	  for( k = 0; k < size; k++ )
	    C[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
      }
}*/
struct void_star{
    double** A,**B,**C;
    int threads,myrank,slice_size;
    pthread_t* t_id ;
};

void * thread_entry(void * vars){
    struct void_star * var; var= vars;
    int i ;
    for (i = 0; i < var->threads; i++)
        if (var->t_id[i]==pthread_self())
            break;
    //printf("past break\n");
    //printf("%i: %i: A: %p\n",var->myrank,i,var->A);
    //printf("%i: %i: B: %p\n",var->myrank,i,var->B);
    //printf("%i: %i: C: %p\n",var->myrank,i,var->C);
    int comps = 0;
    int s,j,k,offset=var->myrank*var->slice_size;
    for (s = i*(var->slice_size/var->threads); s < (i+1)*(var->slice_size/var->threads); s++){\
        for (j = 0; j < var->slice_size; j++){
            for (k = 0; k < M_SIZE; k++){
                //printf("%i: %i: a\n",var->myrank,i);
                //printf("s: %i  j: %i  offset: %i  k: %i\n",s,j,offset,k);
                //printf("A: %lf\n",var->A[s][k]);
                //printf("B: %lf\n",var->B[k][j]);
                var->C[s][(j+offset)] += var->A[s][k] * var->B[k][j];
                //printf("C: %lf\n",var->C[s][(j+offset)]);
                //printf("z\n");
            }
            comps++;
        }
    }
    //printf("%i: %i: done\n",var->myrank,i);
    pthread_exit(NULL);
    return NULL;
}

int main(int argc, char** argv) {
    if (argc != 2){
        printf("Use an interger for threads!\n");
        return 1;
    }
    
    unsigned long long start = rdtsc();
    //Timer has started!
    int i, j;
    int myrank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Status status;
    int slice_size = M_SIZE/numprocs;

    // WHEN USING MPI DO: 
    //if(myrank!=0) return MPI_SUCCESS;
    unsigned long rng_init_seeds[6]={0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
    rng_init_seeds[0] = myrank;
    unsigned long rng_init_length=6;
    init_by_array(rng_init_seeds, rng_init_length);
    
    //slice_size is already divided by the number of threads
    double **A,**B,**C,**D;
    //this is done to have double** with contiguous memory
    A = (double **)calloc( slice_size, sizeof(double*));
    for( i = 0; i < slice_size; i++ ) 
        A[i] = (double *)calloc( M_SIZE, sizeof(double));

    B = (double **)calloc( M_SIZE, sizeof(double*));
    double *contiguousB=(double *)calloc(M_SIZE*slice_size,sizeof(double));
    for( i = 0; i < M_SIZE; i++ ) 
        B[i] = &(contiguousB[slice_size*i]);

    C = (double **)calloc( slice_size, sizeof(double*));
    for( i = 0; i < slice_size; i++ ) 
        C[i] = (double *)calloc( M_SIZE, sizeof(double));

    D = (double **)calloc( M_SIZE, sizeof(double*));
    double *contiguousD=(double *)calloc(M_SIZE*slice_size,sizeof(double));
    for( i = 0; i < M_SIZE; i++ ) 
        D[i] = &(contiguousD[slice_size*i]);

    if (A == NULL){
      printf("A is Null!\n");
      return 0;
    }
    if (B == NULL){
      printf("B is Null!\n");
      return 0;
    }
    if (contiguousB == NULL){
      printf("Contiguous B is Null!\n");
      return 0;
    }
    //Filling A&B matrices
    for( i = 0; i < slice_size; i++ ){
        for( j = 0; j < M_SIZE; j++ )
          {
            //printf("I: %i J: %i\n",i,j);
            A[i][j] = genrand_res53(); // fill A
            B[j][i] = genrand_res53(); // fill B
            //printf("%i,%i,%lf ",j,i,B[j][i]);
          }
        //printf("\n");
    }
    //A is slice_size * M_SIZE
    //B is M_SIZE * slice_size
    //C is slice_size * M_SIZE
    //D is M_SIZE * slice_size
    //slice_size is N/P, M_SIZE is N
    //N is size of matrix, P is number of ranks
    
    MPI_Request request_out;//for checking sends
    MPI_Request request_in; //for checking receives
    
    unsigned long long mm_tmp,isend_tmp,irecv_tmp;
    long long int mm_total=0,isend_total=0,irecv_total=0;
    
    //creation of new threads
    int threads = atoi(argv[1]);
    pthread_t t_id[threads];
    struct void_star vars;
    vars.A = (double**)A; vars.B = (double**)B; vars.C = (double**)C; 
    vars.myrank = myrank; vars.threads=threads; vars.slice_size = slice_size;
    vars.t_id = t_id;
    //printf("C: %p\n",C);
    //printf("Rank %i: B: %p\n",myrank,B);printf("Rank %i: b: %p\n",myrank,vars.B);printf("Rank %i: A: %lf\n",myrank,A[0][0]);printf("Rank %i: a: %lf\n",myrank,vars.A[0][0]);
    
    //runs multiplication for each message
    for (i = 0; i < numprocs; i++)
    {
        //post receive message
        irecv_tmp = rdtsc();
        if (i != 0){
            int flag;
            do {
                MPI_Test(&request_in,&flag,&status);
            } while (!flag);
            //debugging for testing send/receive data validity
            //printf("%i:Testing received message:\n",myrank);
            //if(myrank==0)printf("%i:D:%lf   ",myrank,D[0][0]);
            //swap contents of B and D
            int g,h;
            for (g = 0; g < M_SIZE; g++)
                for (h = 0; h < slice_size; h++){
                    double tmp = D[g][h];
                    D[g][h] = B[g][h];
                    B[g][h] = tmp;
                }
            //if(myrank==0)printf("%i:B:%lf\n",myrank,B[0][0]);
        }
        MPI_Irecv(contiguousD, M_SIZE*slice_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request_in); 
        irecv_total+=rdtsc()-irecv_tmp;
        
        
        //now multiplying the matrices
        //creating threads to do the computation
        pthread_attr_t attr;
        pthread_attr_init(&attr );
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        int t;
        for (t = 0; t < threads; t++){
            //printf("New Thread!\n");
            if (pthread_create(&t_id[t],&attr,&thread_entry,(void*)&vars)) printf("\nError during thread creation!\n");
        }
        pthread_attr_destroy(&attr );
        //printf("%i: attr destroyed\n",myrank);
        //Threads join with original again
        for (t = 0; t < threads; t++){
             //if(!pthread_join(t_id[i],NULL)) printf("%i: Thread %i exited with error.\n",myrank,i);
            int rv = pthread_join(t_id[t],NULL);
            if (rv!=0){
                printf("i: Thread %i exited with error %i .\n",myrank,i,rv);
            }
            else{
                //printf("%i: %i: thread joined successfully\n",myrank,i);
            }
        }
        //printf("%i: All threads closed.\n",myrank);
        mm_tmp = rdtsc();
        
        mm_total+=rdtsc()-mm_tmp;
        
        //posting send message
        isend_tmp = rdtsc();
        if (i != numprocs-1){
            MPI_Isend(contiguousB, M_SIZE*slice_size, MPI_DOUBLE, (myrank+1)%numprocs, 0, MPI_COMM_WORLD, &request_out);
        }
        //printf("%i:Waiting for send to %i\n",myrank,(myrank+1)%numprocs);
        //Blocks sending until previous message is sent
        int flag;
        do {
            MPI_Test(&request_out,&flag,&status);
        } while (!flag);
        
        isend_total+=rdtsc()-isend_tmp;
        //printf("%i:Sent Message!\n",myrank);
    }
    
    //Finalizing Program
    MPI_Finalize();
    if (myrank != 0) return MPI_SUCCESS;
    unsigned long long total = rdtsc()-start;
    
    
    //printf("%i\n",M_SIZE);
    //printf("%i\n",numprocs);
    //multiplication, isend, irecv
    printf("%i size, %i ranks, %i threads, %lld seconds\n", M_SIZE, numprocs, threads, total / 1600000000);
    // MPI_Finalize();
    return (EXIT_SUCCESS);
}
