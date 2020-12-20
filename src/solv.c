#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*define*/

//MODE 1: nibun, MODE 0: newton
#define MODE 1

//newton
#define EP pow(10, -5)
#define INIT 0.01

//nibun
#define XLEFT 0.01
#define XRIGHT 0.4

#define N 60
#define K 30

double f(double x);
double bibun_f(double x);
double newton_sub_f(double x);
void newton_f(double x);
void nibun_f(double x_start1, double x_start2);
double ent(double x);

double f(double x){
    //double y=pow(x, 2)-3;
	double n=(double)N, k=(double)K, y=1.0-ent(x)-(k/n);
    return y;
}

double bibun_f(double x){
    //double y=2*x;
    double y=(log(x)-log(1.0-x))/log(2.0);
    //double y=0.0;
    return y;
}

double newton_sub_f(double x){
    return x-(f(x)/bibun_f(x));
}

void newton_f(double x){
    double x1=x;
    double x2=newton_sub_f(x1);
    int cnt=1;
    
    while(fabs(x1-x2) > EP){
	x1=x2;
	x2=newton_sub_f(x1);
	cnt++;
    }
    
   printf("newton: cnt=%d, solv=%.10lf\n",cnt,x2);

/*
   for(int i=0 ; i<20 ; i++){
       x1=newton_sub_f(x);
       x=x1;
       printf("x1=%lf\n",x1);
   }
*/
}

void nibun_f(double x_start1, double x_start2){
  double a = x_start1;
  double b = x_start2;
  double fc = 0.0, c = 0.0;
  int cnt = 0;
  
  printf("f(%lf)=%lf\n",a,f(a));
  printf("f(%lf)=%lf\n",b,f(b));
	
	
  while(fabs(a-b) > EP){
  	  cnt++;
          c = (a+b)/2.0;
  	  fc = f(c);
	  if(fc == 0.0){
              break;
	  }
	  if(f(a)*fc < 0.0){
              b=c;
	  }
	  if(f(a)*fc > 0.0){
              a=c;
	  }
  	  //printf("nibun: %lf\n",c);
  }

  printf("nibun: cnt=%d, solv=%.10lf\n",cnt,c);
}

double ent(double x){
	double y=-x*log2(x)-(1.0-x)*log2(1.0-x);
	return y;
}

int main(){
        //newton_f(INIT);

	//nibun_f(XLEFT, XRIGHT);
	int n=50, flag=0;
	double l=1.0*(XRIGHT-XLEFT)/n;
	double temp1=0.0, temp2=0.0, temp3=0.0;
	
        #if MODE == 1	
	for(int i=1 ; i<=n ; i++){
		temp1=i*l;
		//printf("x=%lf\n",temp1);
		//newton_f(temp1);
		for(int j=i ; i<=n ; j++){
			temp2=j*l;
			if(f(temp1)*f(temp2) < 0){
				flag=1;
				printf("-----found temp1=%lf, temp2=%lf !!!-----\n",temp1,temp2);
				break;
			}
		}
		if(flag==1){
			break;
		}
	}
	//nibun_f(XLEFT, XRIGHT);
	
	
	if(flag==1){
	    if(temp1 > temp2){
		    temp3=temp1;
		    temp1=temp2;
		    temp2=temp3;
	    }
	    nibun_f(temp1, temp2);
	}else{
		printf("nibun is failed!!\n");
	}
        #endif

        #if MODE == 0
	newton_f(0.0001);
        #endif
	
    return 0;
}
