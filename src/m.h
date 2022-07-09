#define ln5 1.609437
#define ln2 0.693147
#define ln10 2.302585
#define e 2.71828
#define ede2 0.367879
#define pi 3.14159
#define pid2 1.57079
#define pit2 6.28309
#define gra 0.01570796327
#define deg 0.01745329252

#define W productlog
#define lb ld

static double floor(double x);
static double round(double x);
static double ceil(double x);
static double mod(double x, double MOD);
static double sgn(double x);
static double abs(double x);
static double fac(double x);
static double sqrt(double x);
static double cbrt(double x);
static double ln(double x);
static double ld(double x);
static double lg(double x);
static double log(double base, double x);
static double exp(double x);
static double pow(double x,double exponent);

static double sin(double x);
static double cos(double x);
static double tan(double x);

static double arcsin(double x);
static double arccos(double x);
static double arctan(double x);

static double sinh(double x);
static double cosh(double x);
static double tanh(double x);

static double arsinh(double x);
static double arcosh(double x);
static double artanh(double x);

static double csc(double x);
static double sec(double x);
static double cot(double x);

static double csch(double x);
static double sech(double x);
static double coth(double x);

static double arccsc(double x);
static double arcsec(double x);
static double arccot(double x);

static double arcsch(double x);
static double arsech(double x);
static double arcoth(double x);

static double sinc(double x);
static double prodcutlog(int branch,double x);

static double POWN(double x, double exponent);
static double STIRLING_GAMMA(double x);

static double precision=20;

const double NaN = 0.0/0.0;
const double Neg_Inf = 0.0/-1.0;
const double Pos_Inf = 0.0/-1.0;

float REZP[21] = {NaN,1.000000,0.500000,0.333333,0.250000,0.200000,0.166667,0.142857,0.125000,0.111111,0.100000,0.090909,0.083333,0.076923,0.071429,0.066667,0.062500,0.058824,0.055556,0.052632,0.050000};

double FAC[21] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,121645100408832000,2432902008176640000};

static double floor(double x){
	if(x>0){
	return((long long int)x);
	}
	else if(x==(long long int)x){
	return(x);
	}
	else if(x<0){
	return((long long int)(x-1));
	}
	return(NaN);
}
	
static double round(double x){
	return(floor(x+0.5));
}
	
static double ceil(double x){
	if(x==(long long int)(x)){
	return(x);
	}
	else if(x>0){
	return((long long int)(x+1));
	}
	else if(x<0){
	return((long long int)(x));
	}
	return(NaN);
}

static double mod(double x, double MOD){
	return(abs(x-MOD*floor(x/MOD)));
}

static double sgn(double x){
	if(x<0){
		return(-1);
	}
	else if(x>0){
		return(1);
	}
	else if(x==0){
		return(0);
	}
	return(NaN);
}

static double abs(double x){
	if(x<0){
		return(-x);
	}
	else if(x>=0){
		return(x);
	}
	return(NaN);
}

static double fac(double x){
	if(x==0){
		return(1);
	}
	else if(x<-1){
		return(NaN);
	}
	else if(floor(x)==x){
		double X=x;
		for(int i = 1; i < X; i++){
			x=x*i;
		}
		return(x);
	}
	else {
		return(STIRLING_GAMMA(x));
	}
}

static double STIRLING_GAMMA(double x){
	if(x<-1){
		return(NaN);
	}
	else if(x==0){
		return(1);
	}
	else {
		x=x+1;
		return((sqrt((2*pi)/x))*pow((1.0/e)*(x+(1.0/(12*x-(1/(10*x))))),x));
	}
}
	
static double sqrt(double x){
	double X=x/2;
	if(x<0){
		return(NaN);
	}
	else if(x==0){
		return(0);
	}
	else if(x>0){
		for(int i = 1; i <= precision; i++){
			X=((x/X)+X)/2;
		}
		return(X);
	}
	return(NaN);
}

static double cbrt(double x){
	double X=x/3;
	if(x==0){
		return(0);
	}
	else {
		for(int i = 1; i <= precision; i++){
			X=X-(X-(x/(X*X)))/3;	
		}
		return(X);
	}
	return(NaN);
}

static double ln(double x){
	double X=0.0;
	double j=0.0;
	double k=0.0;
	double TMP=1.0;
	if(x<=0){
	return(NaN);
	}
	else if(x<0.5){
	while(j<0.5){
		TMP=5.0*TMP;
		j=x*TMP;
		k++;
	}
	return(ln(j)-k*ln5);
	}
	else if(x>5){
	do {
		TMP=5.0*TMP;
		j=x/TMP;
		k++;
	} while(j>5);
	return(ln(j)+k*ln5);
	}
	else {
	double PRE_REZP=1;
	double PRE=1;
	double PRE_MIDDLE=((x-1)/x);
	for(int i = 1; i < precision; i++){
		PRE_REZP=REZP[i];
		PRE=PRE*PRE_MIDDLE;	
		X=X+PRE_REZP*PRE;
	}
	return(X);
	}
	return(NaN);
}

static double ld(double x){
	return(ln(x)/ln2);
}

static double lg(double x){
	return(ln(x)/ln10);
}

static double log(double base,double x){
	if(base==1 || base <= 0){
		return(NaN);
	}
	else {
		return(ln(x)/ln(base));	
	}
}
	
static double exp(double x){
	double EXP_N=floor(abs(x));
	double epown=POWN(e,EXP_N);
	double EXP_R=abs(x)-EXP_N;
	double HEADER=EXP_R;
	double epowr=1;
	for(int i = 2; i <= 10; i++){
		epowr=epowr+(HEADER)/FAC[i-1];
		HEADER=HEADER*EXP_R;
	}
	if(x==0){
		return(1);
	}
	else if(x>0){
		return(epown*epowr);
	}
	else if(x<0){
		return(1/(epown*epowr));
	}
	return(NaN);
}

static double pow(double x, double exponent){
	double EXP_E=exponent*ln(x);
	return(exp(EXP_E));
}
	
static double POWN(double x, double exponent){
	if(x == 0 && exponent <= 0){
		return(NaN);
	}
	else if(exponent == 1){
		return(x);
	}
	else if(exponent == 0){
		return(1);
	}
	else {
		double X=x;
		for(int i = 1; i < abs(exponent); i++){
			x=x*X;
		}
		if(exponent>0){
			return(x);
		}
		else if(exponent<0){
			return(1/x);
		}
	}
	return(NaN);
}

static double sin(double x){
	double CHECK=(-x+pit2*(floor(x/pit2))+pi)/2;
	if(CHECK>0){
		x=mod(x,pi);
	}
	else if(CHECK<0){
		x=-mod(x,pi);
	}	
	double HEADER=-1;
	double FOOTER=1;
	double X=0;
	double HEADER2=x;
	double XSQUARED=x*x;
	for(int i = 0; i <= 5; i++){
		FOOTER=FAC[2*i+1];
		HEADER=HEADER*-1;
		X=X+((HEADER)/(FOOTER))*HEADER2;
		HEADER2=HEADER2*XSQUARED;
	}
	return(X);
}

static double cos(double x){
	return(sin(x+pi/2));
}

static double tan(double x){
	return(sin(x)/cos(x));
}

static double arcsin(double x){
	if(abs(x)>1){
		return(NaN);
	}
	if(abs(x)==1){
		return(x*pid2);
	}
	return(arctan(x/(sqrt(1-x*x))));
}

static double arccos(double x){
	return(arcsin(-x)+pid2);
}

static double arctan(double x){
	double X=0;
	double x0=abs(x);;
	double dx=0.01;
	double ddx=0;
	if(x>100){
		x=100;
	}
	for(double i = 0; i <= x0; i+=dx){
		ddx=1/(1+i*i);
		X=X+ddx*dx;
	}
	if(x<0){
		return(-X);
	}
	else {
		return(X);
	}
	return(NaN);
}

static double sinh(double x){
	return((exp(x)-exp(-x))/2);
}

static double cosh(double x){
	return((exp(x)+exp(-x))/2);
}

static double tanh(double x){
	return(sinh(x)/cosh(x));
}

static double arsinh(double x){
	return(ln(x+sqrt(x*x+1)));
}

static double arcosh(double x){
	if(x<1){
		return(NaN);
	}
	return(ln(x+sqrt(x*x-1)));
}

static double artanh(double x){
	return(ln((1+x)/(1-x))/2);
}

static double csc(double x){
	return(1/sin(x));
}

static double sec(double x){
	return(1/cos(x));
}

static double cot(double x){
	return(cos(x)/sin(x));
}

static double csch(double x){
	return(1/sinh(x));
}

static double sech(double x){
	return(1/cosh(x));
}

static double coth(double x){
	return(1/tanh(x));
}

static double arccsc(double x){
	return(arcsin(1/x));
}

static double arcsec(double x){
	return(arccos(1/x));
}

static double arccot(double x){
	return(arctan(1/x));
}

static double arcsch(double x){
	return(arsinh(1/x));
}

static double arsech(double x){
	return(arcosh(1/x));
}

static double arcoth(double x){
	return(artanh(1/x));
}

static double sinc(double x){
	if(x==0){
		return(1);
	}
	return(sin(x)/x);
}

static double productlog(int branch,double x){
	double X=5;
	if(x< -ede2){
		return(NaN);
	}
	if(branch==-1 && x>=0){
		return(NaN);
	}
	if(branch==0){
	for(int i = 1; i<= precision; i++){
		X=X-((X*exp(X)-x)/((X+1)*exp(X)));
	}
	return(X);
	}
	else if(branch==-1){
	X=-5;
	for(int i = 1; i<= precision; i++){
		X=X-((X*exp(X)-x)/((X+1)*exp(X)));
		}
	}
	return(X);
}
