/*

This program is based on a reproduction experiment of the following research.

https://uec.repo.nii.ac.jp/?action=pages_view_main&active_action=repository_view_main_item_detail&item_id=8712&item_no=1&page_id=13&block_id=21

This program is using gcc (clang) built in function. Compile to use gcc or clang.

*/

/*include*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <gmp.h>
#include "../headers/mt64.h"

/*typedef*/
typedef unsigned long int ulong;
typedef unsigned long long int ulonglong;
typedef long double ldouble;

/*define*/
/*check: user check data N, K, DUMP_NAME, H_DUMP_NAME, TIMES*/
//==============================
#define N 400  //check
#define K 370 //check
#define M N-K  

/*name: ./result/N_K_result.log*/
#define DUMP_NAME "./result/400_370_result.log" //check
#define H_DUMP_NAME "./result/400_370_result_h.log" //check

#define TIMES 5  //check
//==============================

#define LIMIT_M 60

#define L 200

/*if: preset*/
#if (N==80 && K==40) || (N==70 && K==35) || (N==60 && K==30) || (N==50 && K==25) || (N==40 && K==20) || (N==30 && K==15) || (N==20 && K==10) || (N==10 && K==5) 
#define THEORY 0.11

#elif N==400 && K==370
#define THEORY 0.009

#elif N==60 && K==20
#define THEORY 0.174

/*preset to add here*/
#endif

/*ifndef: random or not random*/
#ifndef ENABLE_RANDOM
#define ENABLE_RANDOM 1
#endif

/*if: ENABLE_RANDOM*/
#if defined(ENABLE_RANDOM) && ENABLE_RANDOM==0
#undef TIMES
#define TIMES 1
#endif

/*function defitions*/
int popcnt(ulonglong x);
int parity(ulonglong x);
void genrand(ulonglong *h, ulonglong MbitMax);
ulonglong beki();
void before_f(ulonglong *h, ulonglong *b, ulonglong MbitMax);
double bound_theory_f(double x, double y, double z);
void exchange_f(ldouble *a, ldouble *b, int c, ulonglong *h1, ulonglong *h2, int norm);
void copy_f(ldouble *a, ldouble *b, ulonglong *c, ulonglong *d, int mode);
void after_f(ulonglong *b, ulonglong MbitMax, ldouble *Result2Norm, ldouble *Result1Norm, int *ind_theory, double *alpha);

/*function popcnt*/
int popcnt(ulonglong x){
	return (__builtin_popcountll(x));	
}

/*function parity*/
int parity(ulonglong x){
	return (__builtin_parityll(x));	
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
		    h[i]=(ENABLE_RANDOM==1)?genrand64_int64()%MbitMax:0;
		}
	}
}

/*function beki*/
ulonglong beki(){
	//printf("beki=%llu\n",(ulonglong)(1 << M));
	return (M >=1 && M<30)?((ulonglong)(1 << (M))):((ulonglong)pow(2.0, (double)M));
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
void exchange_f(ldouble *a, ldouble *b, int c, ulonglong *h1, ulonglong *h2, int norm){
    ldouble temp=0.0;
    ulonglong h_temp=0;
    if(a[c] > b[c]){
        for(int i=0 ; i<=L ; i++){
            temp=a[i];
	    a[i]=b[i];
	    b[i]=temp;

        }
	if(norm==1){
	    /*a: h1, b: h2*/
	    for(int i=0 ; i<N ; i++){
                h_temp=h1[i];
	        h1[i]=h2[i];
	        h2[i]=h_temp;
	    }
	}
    }
}


/*function copy_f*/
void copy_f(ldouble *a, ldouble *b, ulonglong *c, ulonglong *d, int mode){
    for(int i=0 ; i<=L ; i++){
        a[i]=b[i];
    }
    if(mode==1){
        for(int i=0 ; i<N ; i++){
            c[i]=d[i];
        }
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
		temp1=0.0;
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
		Result2Norm[i]=sqrtl((ldouble)(temp1)/(ldouble)(MbitMax));
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
	ulonglong hdata1[N]={0.0};

	int ind_theory=0;

	ldouble Result2Norm[L+1]={0.0}, Result1Norm[L+1]={0.0};
	ldouble Result2Norm_data1[L+1]={0.0}, Result1Norm_data1[L+1]={0.0};


	double a[L+1]={0.0};
	FILE *outputfile, *hdata_outputfile;

	#ifdef DUMP_NAME
	    outputfile = fopen(DUMP_NAME, "w"); 
	    if (outputfile == NULL) { 
	        printf("-----DUMP_NAME: can not open-----\n");         
                return EXIT_FAILURE;                       
            }else{
                printf("-----DUMP_NAME: file pointer is ok !!-----\n");
	    }
        #endif
	#ifndef DUMP_NAME
	    printf("-----DUMP_NAME: not defined-----\n");
	    return EXIT_FAILURE;
        #endif

	#ifdef H_DUMP_NAME
	    hdata_outputfile = fopen(H_DUMP_NAME, "w"); 
	    if (hdata_outputfile == NULL) { 
	        printf("-----H_DUMP_NAME: can not open-----\n");         
                return EXIT_FAILURE;                       
            }else{
                printf("-----H_DUMP_NAME: file pointer is ok !!-----\n");
	    }
        #endif
	#ifndef H_DUMP_NAME
	    printf("-----H_DUMP_NAME: not defined ----\n");
	    return EXIT_FAILURE;
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
	        copy_f(Result1Norm_data1, Result1Norm, hdata1, h, 1);
	        copy_f(Result2Norm_data1, Result2Norm, hdata1, h, 0);
		//copy_f(hdata1, h);

	     }
	     if(t >= 1){
                 exchange_f(Result1Norm_data1, Result1Norm, ind_theory, hdata1, h, 1);
                 exchange_f(Result2Norm_data1, Result2Norm, ind_theory, hdata1, h, 2);
	     }
	     printf("\n");
	}
	  
        
	     /*dump data*/
	     printf("-----start data dump-----\n");
	     for(int i=0 ; i<=L ; i++){
	         fprintf(outputfile, "%f %f %f %f %.23LE %.23LE \n"
				 ,a[i],2.0*a[i],THEORY,(double)(MbitMax*THEORY)
				 ,Result2Norm_data1[i],Result1Norm_data1[i]); 
	         //fprintf(outputfile, "%f %f %.23LE %.23LE \n",a[i],THEORY,Result2Norm_data1[i],Result1Norm_data1[i]); 
	     }
	     printf("-----end data dump-----\n");

	     printf("-----start H=[I, H2] dump-----\n");
	     
	     fprintf(hdata_outputfile, "-----This is H=[I, H2], N=%d, M=%d, K=%d, H: M x N .-----\n",N,M,K); 
	     for(int i=0 ; i<N ; i++){
		 if(i==0){
	             fprintf(hdata_outputfile, "[%lld, ",hdata1[i]); 

		 }else if(i==N-1){
	             fprintf(hdata_outputfile, "%lld]\n",hdata1[i]); 

		 }else{
	             fprintf(hdata_outputfile, "%lld, ",hdata1[i]);

		 }
	     }
	     printf("-----end H=[I, H2] dump-----\n");

             /*close file*/
	     fclose(outputfile);
	     fclose(hdata_outputfile);
	
	     printf("-----success to run !!!-----\n");
             return EXIT_SUCCESS;
}
