/* 2020ass06_JunshengZhu */
/* two_way analysis - RCB model - two fixed effect*/

proc iml;
reset print;
/* 					Trts 			   sex
	mu 		1   2   3   4   5   	C   F	M
*/
x ={1		1	0	0	0	0		1	0	0,
	1		1	0	0	0	0		0	1	0,
	1		1	0	0	0	0		0	0	1,
	1		0	1	0	0	0		1	0	0,
	1		0	1	0	0	0		0	1	0,
	1		0	1	0	0	0		0	0	1,
	1		0	0	1	0	0		1	0	0,
	1		0	0	1	0	0		0	0	1,
	1		0	0	0	1	0		1	0	0,
	1		0	0	0	1	0		0	1	0,
	1		0	0	0	1	0		0	0	1,
	1		0	0	0	0	1		1	0	0,
	1		0	0	0	0	1		0	1	0,
	1		0	0	0	0	1		0	0	1};
	
y = {21.4,
	21.1,
	23.0,
	23.0, 
	22.7,
	23.8,
	25.7,
	27.6,
	27.2,
	26.8,
	28.8,
	28.4,
	23.6,
	26.3};

	xtx = x`* x;
	xty = x`* y;

	invxtx = ginv(xtx);
	b = invxtx * xty;
	tss = y` * y;
	sumy = sum(y);
	ssr = b`* xty;
	nobs = nrow(x);
	ybar = sumy/ nobs;
	cf = nobs * ybar *ybar;
	ssrm = ssr - cf;
	dfd = 4;
	dfp = 2;
	rx = 1 + dfd + dfp;
	dfe = nobs - rx;
	sse = tss - ssr;
	mse = sse/dfe;

/* Type III Sums of Squares, Marginal, for Diet */	
print /;
print " Type III Sums of Squares, Marginal, for Diets ";
kp= {   0	1 -1 0 0 0	  0 0 0,
		0 	1 0 -1 0 0	  0 0 0,
		0 	1 0 0 -1 0	  0 0 0,
		0 	1 0 0 0 -1 	  0 0 0};
k = kp`;
df = nrow(kp);
kb = k` * b;
kxxk = k` * invxtx * k;
invkk = ginv(kxxk);
ssd = kb` * invkk * kb;
msd = ssd/df;
fd = msd/mse;
pr = 1 - probf(fd,df,dfe);
print "Fcal(diets) =" fd "pr =" pr;

/* Type III Sums of Squares, Marginal, for Sex */
print /;
print " Type III Sums of Squares, Marginal, for Sex ";
kp= {0	0 0 0 0 0 	 1 -1 0 ,
	 0 	0 0 0 0 0 	 1 0 -1 };
k = kp`;
df = nrow(kp);
kb = k` * b;
kxxk = k` * invxtx * k;
invkk = ginv(kxxk);
ssp = kb` * invkk * kb;
msp = ssp/df;
fp = msp/mse;
pr = 1 - probf(fp,df,dfe);
print "Fcal(sex) =" fp "pr =" pr;


/* Type III Sums of Squares for Diet and Sex */
print /;
print " Type III Sums of Squares for Diets and Sexes ";
kp = {	0	1 -1 0 0 0	  0 0 0,
		0 	1 0 -1 0 0	  0 0 0,
		0 	1 0 0 -1 0	  0 0 0,
		0 	1 0 0 0 -1 	  0 0 0,
		0	0 0 0 0 0    1 -1 0,
	 	0 	0 0 0 0 0    1 0 -1};

k = kp`;
df = nrow(kp);
kb = k` * b;
kxxk = k` * invxtx * k;
invkk = ginv(kxxk);
ssp = kb` * invkk * kb;
msp = ssp / df;
fp = msp / mse;
pr = 1 - probf(fp, df, dfe);
print "Fcal(SSRm) =" fp "pr =" pr;


/* Q8 estimate the missing values */
kp = {1		0 0 1 0 0 	0 1 0};
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print " estimate of missing value Y_3F =" kb "+ -" se;

/* Q10 estimates of fitted values, using the general k'b approach */

/* diet1, lsmean */
kp = {3		3 0 0 0 0  	1 1 1 }/3;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  diet1, lsmean =" kb "+ -" se;


/* diet2, lsmean */
kp = {3		0 3 0 0 0  	1 1 1 }/3;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  diet2, lsmean =" kb "+ -" se;

/* diet3, lsmean */
kp = {3		0 0 3 0 0  	1 1 1 }/3;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  diet3, lsmean =" kb "+ -" se;

/* diet4, lsmean */
kp = {3		0 0 0 3 0  	1 1 1 }/3;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  diet4, lsmean =" kb "+ -" se;

/* diet5, lsmean */
kp = {3		0 0 0 0 3  	1 1 1 }/3;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  diet5, lsmean =" kb "+ -" se;

/* sexC, lsmean */
kp = {5		1 1 1 1 1  	5 0 0 }/5;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  sexC, lsmean =" kb "+ -" se;

/* sexF, lsmean */
kp = {5		1 1 1 1 1  0 5 0}/5;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  sexF, lsmean =" kb "+ -" se;


/* sexM,lsmean */
kp = {5		1 1 1 1 1  	0 0 5}/5;
k = kp`;
kb = k` * b;
sv = k` * invxtx * k * mse;
se = sqrt(sv);
print "  sexM,lsmean =" kb "+ -" se;

/*compute tabulated F values */
data Fvalues;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
cards;
7 7 0.05
1 7 0.05
6 7 0.05
4 7 0.05
2 7 0.05
;
proc print data = Fvalues;
var ndf ddf prob ftab;
run;





data ass06_twoway;
input d sex$ y;
cards;
1	C	21.4
1	F	21.1
1	M	23.0 
2	C	23.0 
2	F	22.7
2	M	23.8
3	C	25.7
3	M	27.6
4	C	27.2
4	F	26.8
4	M	28.8
5	C	28.4
5	F	23.6
5	M	26.3
;
proc glm data = ass06_twoway;
class d sex;
model y  = d sex/ xpx i solution;

/* Estimate each fitted value */

estimate  	'diet1, lsmean' intercept	3	d	3	0	0	0	0	sex	1	1	1/divisor = 3;
estimate  	'diet2, lsmean' intercept	3	d	0	3	0	0	0	sex	1	1	1/divisor = 3;
estimate  	'diet3, lsmean' intercept	3	d	0	0	3	0	0	sex	1	1	1/divisor = 3;
estimate  	'diet4, lsmean' intercept	3	d	0	0	0	3	0	sex	1	1	1/divisor = 3;
estimate  	'diet5, lsmean' intercept	3	d	0	0	0	0	3	sex	1	1	1/divisor = 3;
estimate  	'sexC, lsmean' intercept	5	d	1 	1	1	1	1	sex	5	0	0/divisor = 5;
estimate  	'sexF, lsmean' intercept	5	d	1 	1	1	1	1	sex	0	5	0/divisor = 5;
estimate  	'sexM, lsmean' intercept	5	d	1 	1	1	1	1	sex	0	0	5/divisor = 5;

lsmeans d/pdiff stderr;
lsmeans sex/pdiff stderr;
run;
quit;
