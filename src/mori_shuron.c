/*

This program is based on a reproduction experiment of the following research.

https://uec.repo.nii.ac.jp/?action=pages_view_main&active_action=repository_view_main_item_detail&item_id=8712&item_no=1&page_id=13&block_id=21

*/

/*include*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <gmp.h>
#include "../headers/mt64.h"

/*typedef*/
typedef unsigned long int ulong;
typedef unsigned long long int ulonglong;
typedef long double ldouble;

/*define*/
#define L 50

#define N 400
#define M 40
#define K N-M

#define LIMIT_M 60

/*ifdef: preset*/
#if N/K==2 && N%K==0
#define THEORY 0.0624
#elif N==400 && K==370
#define THEORY 0.0117
#elif N==60 && K==20
#define THEORY 0.0936
/*preet to add here*/
#endif


/*function defitions*/
ulonglong popcnt(ulonglong x);
ulonglong parity(ulonglong x);
void genrand(ulonglong *h, ulonglong MbitMax);
ulonglong beki();
void before_f(ulonglong *h, ulonglong *b, ulonglong MbitMax);
void after_f(ulonglong *b, ulonglong MbitMax);

/*function popcnt*/
ulonglong popcnt(ulonglong x){
	return (ulonglong)(__builtin_popcountll(x));	
}

/*function parity*/
ulonglong parity(ulonglong x){
	return (ulonglong)(__builtin_parityll(x));	
}

/*function genrand*/
void genrand(ulonglong *h, ulonglong MbitMax){
	//srand((unsigned)time(NULL));
	init_genrand64((unsigned)time(NULL));
	for(int i=0 ; i<N ; i++){
		if(i<M){
			h[i]=(ulonglong)((i==0)?1:h[i-1]*2);
		}else{
	        //printf("%lld\n",genrand64_int64()%10);
		    h[i]=genrand64_int64()%MbitMax;
		}
	}
}

/*function beki*/
ulonglong beki(){
	//printf("beki=%llu\n",(ulonglong)(1 << M));
	return (M >=1 && M<30)?((ulonglong)(1 << M)):((ulonglong)pow(2.0, (double)M));
}

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

void after_f(ulonglong *b, ulonglong MbitMax){
    double a=0.0;
	ldouble width=0.5/L;
	ldouble temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0;
	
	printf("a=%f, 2-norm=%.20Le, 1-norm=%.20Le\n"
	       ,0.0,sqrtl((ldouble)(MbitMax-1)/(ldouble)(MbitMax))
	       ,sqrtl((ldouble)(MbitMax-1)));
	
	for(a=width ; a<0.5 ; a+=width){
		temp4=0.0;
		for(int j=1 ; j<=N ; j++){
			if(j==1){
				temp3=2.0*logl(-2.0*a+1.0);
				temp4=temp3;
			}else{
				temp4+=temp3;
			}
			if(b[j]!=0){
				temp2=(ldouble)(b[j]*expl(temp4));
				temp1+=temp2;
			}
		}
		//printf("a=%lf\n",a);
		printf("a=%f, 2-norm=%.20Le, 1-norm=%.20Le\n"
	       ,a,sqrtl((ldouble)(temp1))/sqrtl((ldouble)(MbitMax))
		   ,sqrtl((ldouble)(temp1)));
	}
	
	printf("a=%f, 2-norm=%.20e, 1-norm=%.20e\n",a,0.0,0.0);

}


/*function main*/
int main(){
        /*ifndef: preset is not defined*/
        #ifndef THEORY
            printf("-----preset is not defined !!-----\n");
            printf("-----please run ./solv !!-----\n");
            return EXIT_FAILURE;
        #endif
	
	/*parameter check*/
	if(M > LIMIT_M){
		printf("-----M is too large.-----\n");
		return EXIT_FAILURE;
	}else if((N <= 0) || (M <= 0) || (K <= 0) || (N <= M) || (N <= K)){
		printf("-----setting error!!-----\n");
		return EXIT_FAILURE;
	}else{
		printf("-----setting is ok!!-----\n");
	}
	
	ulonglong h[N]={0}, b[N+1]={0};
	ulonglong MbitMax=beki();
	//printf("%llu\n",MbitMax);
	
	/*generate Parity Check H=[I, H2]*/
	printf("-----start to generate H=[I, H2]-----\n");
	genrand(h, MbitMax);
	printf("-----end to generate H=[I, H2]-----\n");
	
	/*print Parity Check H=[I, H2]*/
	printf("-----print H=[I, H2]-----\n");
	for(int i=0 ; i<N ; i++){
		if(i==0){
			printf("[");
		}
		/*
		if(i==M){
		    printf("----------------------\n");
		}
		*/
		printf("%lld",h[i]);
		if(i==N-1){
			printf("]\n");
		}else{
			printf(",");
		}
	}
	
	printf("-----start before_f-----\n");
	before_f(h, b, MbitMax);
	printf("-----end before_f-----\n");
	
	/*
	for(int j=1 ; j<=N ; j++)
	    printf("b[%d]=%lld\n",j,b[j]);
	*/
	
	printf("-----start after_f-----\n");
	after_f(b, MbitMax);
	printf("-----end after_f-----\n");
	
	
    return EXIT_SUCCESS;
}
