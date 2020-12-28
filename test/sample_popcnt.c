
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>


#define M 30
#define N 400
#define K N-M

#define MAX_THREADS omp_get_max_threads()

typedef unsigned long long int ulonglong;

int popcnt(ulonglong x);
int parity(ulonglong x);
void sample_popcnt();

void before_f(ulonglong *h, ulonglong *b, ulonglong MbitMax);

/*function popcnt*/
int popcnt(ulonglong x){
	return (__builtin_popcountll(x));	
}

/*function parity*/
int parity(ulonglong x){
	return (__builtin_parityll(x));	
}


void sample_popcnt(){
    int m=M, l=0;
    ulonglong MbitMax = (ulonglong)(pow(2, m)) - 1;
    ulonglong b[N+1]={0}, h[N+1]={0};
    ulonglong temp=0, j=0;
    double st=0.0, ed=0.0;

    srand((unsigned)time(NULL));
    for(int i=0 ; i<N ; i++){
	if(i==0){
            h[0]=1;
	}else if(i>=1 && i<=M-1){
            h[i]=h[i-1]*2;
	}else{
            h[i]=rand()%MbitMax;

	}
	//printf("%llu, ",h[i]);
    }
        //printf("\n");
    st=omp_get_wtime();
    #ifdef PARA
    #pragma omp parallel num_threads(MAX_THREADS)
    {
        #pragma omp for private(l, temp) reduction(+:j) schedule(static)
	for(ulonglong t=1 ; t<MbitMax ; t++){
		temp=popcnt(t);
		j=0;
		for(l=M ; l<N ; l++){
			j+=parity(t & h[l]);
		}
        #pragma omp atomic
		b[temp+j]+=1;
	}
    }
    #endif
    #ifndef PARA
    for(ulonglong t=1 ; t<MbitMax ; t++){
        temp=popcnt(t);
        j=0;
	for(l=M ; l<N ; l++){
		j+=parity(t & h[l]);
	}
	b[temp+j]+=1;
    }
    #endif
    ed=omp_get_wtime();


    printf("-----result !!-----\n");

    for(int j=0 ; j<=N ; j++){
        printf("b[%d]: %lld\n",j,b[j]);
    }

    printf("time: %fs\n",ed-st);

}

/*function before_f: frequency b[j]*/
void before_f(ulonglong *h, ulonglong *b, ulonglong MbitMax){
	ulonglong j=0;
	for(ulonglong t=1 ; t<MbitMax ; t++){
		j=popcnt(t);
		//printf("j's popcnt: %lld\n",j);
		for(int l=M ; l<N ; l++){
			if(parity(t & h[l])==1){
				j++;
			}
		}
		b[j]++;
		//printf("b[j]: %lld\n",j);
	}
}

int main(){
    sample_popcnt();
    return 0;
}
