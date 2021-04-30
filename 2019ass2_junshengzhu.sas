/* AEMA 610 ASS2 */

data ass02;
input plot diameter thickness y;
radius = diameter thickness y;
pi = constant("pi");
vol = 4 * Pi *(redius**3)/3.0;
vol2 = vol * vol;

drop pi radiusl;
cards;
1 5  8 11.9
2 5  9 15.1
3 5  10 14.5
4 5
5 5
6 6
7 6
8 6
9 6
10 6
11 7
12 7
13 7
14 7
15 7
16 8
17 8
18 8
19 8
20 8
21 9
22 9
23 9
24 9
25 9
26 10
27 10
28 10
29 10
30 10
;
/* check by reading data */
proc print data = ass02;
var plot diameter vol vol2 thickness y;
run;

/* calculate means, SD, range */
proc means data = ass02;
var diameter vol vol2 thickness y;
run;

/* using IML */
proc iml;
reset print;

/*columns*/
/*		1 Vol Vol2 Thickness */
X = {1 65.450 4283.68  8,
	 1 65.450 4283.68  9,
	 1 65.450 4283.68  10};

y = {11.9, 
	 15.2,
	 14.5};

print "X transpose X";
xtx = x` * x;
print " ";
print "X transpose Y";
xty = x` * y;

tss = y` *y;

invxtx = inv(xtx);

bhat = invxtx * xty;

sumy = sum(y);
nobs = nrow(x);
rx = 4;
dfe = nobs - rx;
ybar = sumy/nobs;
cf = nobs * ybar * ybar;

ssr = bhat` * xty;
ssrm = ssr - cf;

sse = tss - ssr
mse = sse/dfe;

covb = inxtx * mse;

sv0 = covb[1,1];
se0 = sqrt(sv0);
t1 = bhat[1]/se0;

sv1 = covb[2,2];
se1 = sqrt(sv1);
t2 = bhat[2]/se1;

sv2 = covb[3,3];
se2 = sqrt(sv2);
t3 = bhat[3]/se2;

sv3 = covb[4,4];
se3 = sqrt(sv3);
t4 = bhat[4]/se3;

/* generate SS by K matrix */
print " ";
print " F tests ";

/* SS b1 */
kp = {0 1 0 0};
kb = kp * bhat;
kinvk = kp*invxtx*kp`;
invkk = inv(kinvk);
ss1 = kb`*invkk*kb;

/* SS b2 */
kp = {0 0 1 0};
kb = kp * bhat;
kinvk = kp*invxtx*kp`;
invkk = inv(kinvk);
ss2 = kb`*invkk*kb;
print" F calc for b2 = ";
f2 = ss2/mse;

/* SS b3 */
kp = {0 0 0 1};
kb = kp*bhat;
kinvk = kp*inxtx*kp`;
invkk = inv(kinvk);
ss3 = kb`*invkk*kb;
print "F calc for b3 = ";
f3 = ss3/mse;

yhat = x * bhat;

ehat = y - yhat;

/* calculate V optimum, from derivatives */
vopt = -bhat[2]/(2*bhat[3]);

quit;

/*compute tabulated F values */
data F values;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
cards;
4 26 0.05
1 26 0.05
3 26 0.05
;

proc print data = Fvalues;
var ndf ddf prob ftab;
run;

/* compute chi_square for random effect */
data Chivalues;
input df prob;
chitab = cinv(1 - prob, df);
cards;
26 0.025
26 0.975
;

proc print data = Chivalues;
var df prob chitab;
run;

proc gplot data = ass02;
plot y * vol;
run;
quit;

/* Check with GLM */
ods output Estimates = est1;
proc glm data = ass02;
model y = thickness vol vol*vol/xpx inverse solution;
output out=ass02r p=yhat r=ehat;
estimate 'intercept' intercept 1;
estimate 'vol' vol 1;
estimate 'vol^2' vol*vol 1;
estimate 'thickness' thickness 1;
run;
quit;

/* to see the names and types of variables  */
proc contents data = est1;
run;

/* to print the values of observations */
proc print data = est1;
run;

/* to collect four variables in one records */
data est2;
retain  b0 b1 b2 b3;
set est1;
if (Parameter eq 'intercept') then b0 = Estimate;
if (Parameter eq 'vol') then b2 = Estimate;
if (Parameter eq 'vol^2') then b3 = Estimate;
if (Parameter eq 'thickness') then b1 = Estimate;
if (Parameter eq 'thickness') then output;
keep b0 b1 b2 b3;

run;

proc print data = est2;
run;

/* using the model to estimate */
data est3;
set est2;
Pi = constant("pi");
do d = 5 to 12 by 0.25;
vol = (4.0 * Pi * (d/2.0)**3) / 3.0;
Ypredict = b0 + b1*10.0 + b2*vol + b3*vol*vol;
output;
end;
run;

title'Ypredict, average thickness';
proc gplot data = est3;
plot Ypredict*vol;
run;
quit;

title;
proc contents data = ass02r;
run;

/*plot residuals against volume and biofilm thickness */
title 'Residuals';
proc gplot data = ass02r;
plot ehat*vol;
plot ehat*thickness;
run;
quit;





