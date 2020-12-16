/*

This program is based on a reproduction experiment of the following research.

https://uec.repo.nii.ac.jp/?action=pages_view_main&active_action=repository_view_main_item_detail&item_id=8712&item_no=1&page_id=13&block_id=21

This program is using gcc built in function. Compile to use gcc.

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
#define L 60

#define N 60
#define M 30
#define K N-M
#define LIMIT_M 60

#define DUMP_NAME "../result/60_30_result.log"
#define TIMES 2

/*ifdef: preset*/
#if (N==80 && K==40) || (N==70 && K==35) || (N==60 && K==30) || (N==50 && K==25) || (N==40 && K==20) || (N==30 && K==15) || (N==20 && K==10) || (N==10 && K==5)
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
double bound_theory_f(double x, double y, double z);
void exchange_f(ldouble *a, ldouble *b, int c);
void after_f(ulonglong *b, ulonglong MbitMax, ldouble *Result2Norm, ldouble *Result1Norm, int *ind_theory, double *alpha);

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

/*function bound_theory_f: boundary THEORY*/
double bound_theory_f(double x, double y, double z){
    double temp1=fabs(x-z), temp2=fabs(y-z); 
    return (temp1<=temp2)?x:y;
}

/*function exchange_f: select min(a, b), return *a*/
void exchange_f(ldouble *a, ldouble *b, int c){
    ldouble temp=0.0;
    if(a[c] > b[c]){
        for(int i=0 ; i<=L ; i++){
            temp=a[i];
	    a[i]=b[i];
	    b[i]=temp;

        }
    }
}


/*function copy_f*/
void copy_f(ldouble *a, ldouble *b){
    for(int i=0 ; i<=L ; i++){
        a[i]=b[i];
    }
}

/*function after_f: calculate norm to evalulate*/
void after_f(ulonglong *b, ulonglong MbitMax, ldouble *Result2Norm, ldouble *Result1Norm, int *ind_theory, double *alpha){
        int i=0;
	double a=0.0, theory=THEORY, v=0.0, width=0.5/L;
	ldouble temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0;

        /*	
	printf("a=%f, 2-norm=%.20Le, 1-norm=%.20Le\n"
	       ,0.0,sqrtl((ldouble)(MbitMax-1)/(ldouble)(MbitMax))
	       ,sqrtl((ldouble)(MbitMax-1)));
        */

	Result2Norm[0]=sqrtl((ldouble)(MbitMax-1)/(ldouble)(MbitMax));
	Result1Norm[0]=sqrtl((ldouble)(MbitMax-1));
	alpha[0]=0.0;
	
	
	for(a=width ; a<0.5 ; a+=width){
		temp4=0.0;
		if(a<=theory && a+width>theory){
                    v=bound_theory_f(a, a+width, theory);
		    *ind_theory=(int)(v/width);
		}
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
		/*
		printf("a=%f, 2-norm=%.20Le, 1-norm=%.20Le\n"
	       ,a,sqrtl((ldouble)(temp1))/sqrtl((ldouble)(MbitMax))
	       ,sqrtl((ldouble)(temp1)));
		*/
                i++;
		Result2Norm[i]=sqrtl((ldouble)(temp1))/sqrtl((ldouble)(MbitMax));
		Result1Norm[i]=sqrtl((ldouble)(temp1));
		alpha[i]=a;
	}
	/*
	printf("a=%f, 2-norm=%.20e, 1-norm=%.20e\n",a,0.0,0.0);
	*/
	Result2Norm[L]=0.0;
	Result1Norm[L]=0.0;
	alpha[L]=0.5;

	printf("i+1=%d, L=%d\n",i+1,L);
	printf("theory=%lf, v=%lf, width=%lf\n",theory,v,width);

}


/*function main*/
int main(){
        /*ifndef: preset is not defined*/
        #ifndef THEORY
            printf("-----preset is not defined !!-----\n");
            printf("-----please run ./solv !!-----\n");
            return EXIT_FAILURE;
        #endif
	
	#ifndef DUMP_NAME
	    printf("-----DUMP_NAME is not defined !!-----\n");
	    printf("-----please set to DUMP_NAME !!-----\n");
	    return EXIT_FAILURE;
	#endif
	
	/*parameter check*/
	if(M > LIMIT_M){
		printf("-----M is too large.-----\n");
		return EXIT_FAILURE;
	}else if((N <= 0) || (M <= 0) || (K <= 0) || (N <= M) || (N <= K)){
		printf("-----setting error!!-----\n");
		return EXIT_FAILURE;
	}
	else{
		printf("-----setting is ok!!-----\n");
	}
	
	ulonglong h[N]={0}, b[N+1]={0};
	ulonglong MbitMax=beki();
	int ind_theory=0;

	ldouble Result2Norm[L+1]={0.0}, Result1Norm[L+1]={0.0};
	ldouble Result1Norm_data1[L+1]={0.0}, Result1Norm_data2[L+1]={0.0};

	double a[L+1]={0.0};
	FILE *outputfile;

	#ifdef DUMP_NAME
	    outputfile = fopen(DUMP_NAME, "w"); 
	    if (outputfile == NULL) { 
	        printf("-----cannot open-----\n");         
                return EXIT_FAILURE;                       
            }else{
                printf("-----file pointer is ok !!-----\n");
	    }
        #endif

	//printf("%llu\n",MbitMax);
	
	/*exchange for loop*/
	for(int t=0 ; t<TIMES ; t++){
	     printf("\n-----%d times: random coding, Generate H=[I, H2]-----\n\n",t+1);

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
	     after_f(b, MbitMax, Result2Norm, Result1Norm, &ind_theory, a);
	     printf("ind_theory=%d\n",ind_theory);
	     printf("-----end after_f-----\n");


	     if(t == 0){
	        copy_f(Result1Norm_data1, Result1Norm);

	     }
	     if(t >= 1){
                 exchange_f(Result1Norm_data1, Result1Norm, ind_theory);
	     }
	     printf("\n");
	}
	  
        
	     /*dump data*/
	     printf("-----start data dump-----\n");
	     for(int i=0 ; i<=L ; i++){
	         fprintf(outputfile, "%f, %f, %.23LE, %.23LE \n",a[i],THEORY,Result2Norm[i],Result1Norm[i]); 
	     }
	     printf("-----end data dump-----\n");

             /*close file*/
	     fclose(outputfile);
	
	     printf("-----success to run !!!-----\n");
             return EXIT_SUCCESS;
}
