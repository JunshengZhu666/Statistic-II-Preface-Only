/* 2020ass5_JunshengZhu */

/* Hail Dr.Roger */

proc iml;
reset print;
x = 
{1		1	0	0	0	0	0	,
1		1	0	0	0	0	0	,
1		1	0	0	0	0	0	,
1		1	0	0	0	0	0	,
1		1	0	0	0	0	0	,
1		1	0	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	1	0	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	1	0	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	1	0	0	,
1		0	0	0	0	1	0	,
1		0	0	0	0	1	0	,
1		0	0	0	0	1	0	,
1		0	0	0	0	1	0	,
1		0	0	0	0	1	0	,
1		0	0	0	0	0	1	,
1		0	0	0	0	0	1	,
1		0	0	0	0	0	1	,
1		0	0	0	0	0	1	,
1		0	0	0	0	0	1	,
1		0	0	0	0	0	1};

y = 
{15.9	,
15.6	,
15.4	,
14.0 	,
17.7	,
15.6	,
18.7	,
18.9	,
20.4	,
19.8	,
19.3	,
19.3	,
21.8	,
20.5	,
18.7	,
20.8	,
20.6	,
19.7	,
19.9	,
21.3	,
21.4	,
17.1	,
20.3	,
19.8	,
20.2	,
20.6	,
19.7	,
21.2	,
21.2	,
14.3	,
15.3	,
15.6	,
14.1	,
12.9	,
15.0 };

print" X transpose X";
xtx = x` * x;
print" X transpose Y";
xty = x` * y;

tss = y` * y;

invxtx = ginv(xtx);
btilde = invxtx * xty;

sumy = sum(y);
nobs = nrow(x);
rx = 6;
dfe = nobs - rx;
ybar = sumy/nobs;
cf = nobs * ybar *ybar;

ssr = btilde` * xty;
ssrm = ssr - cf;

msr = ssr/rx;
msrm = ssrm/(rx - 1);

sse = tss - ssr;
mse = sse/dfe;

fm = msr/mse;
fmean = cf/mse;
fmrm = msrm/mse;



/* SS for trts */
kp = { 0 1 -1 0 0 0,
	   0 1 0 -1 0 0,
	   0 1 0 0 -1 0,
	   0 1 0 0 0 -1};
kb = kp* btilde;
kinvk = kp * invxtx *kp`;
invkk = inv(kinvk);
sstrt = kb` * invkk * kb;
dftrt = nrow(kp);
mstrt = sstrt/dftrt;

/* differences between trts and control */
kp = {0 1 -1 0 0 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-2 " kb se;

kp = {0 1 0 -1 0 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-3 " kb se;

kp = {0 1 0 0 -1 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-4 " kb se;

kp = {0 1 0 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-5 " kb se;

/* differences between company A and B */
kp = {0 0 0.5 -0.5 0.5 -0.5};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2.4. vs 3.5. " kb se;


/*Q10 individual t_test */
proc iml;
reset print;
dfe = 23;
pr = 0.01;
prnominal = pr/4;
ttab = tinv(1- prnominal/2, dfe);


quit;

/*compute tabulated F values */
data Fvalues;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
cards;
1 23 0.01
4 23 0.01
5 23 0.01
;
proc print data = Fvalues;
var ndf ddf prob ftab;
run;


data ass03;
input people diet$ wt;
cards;
1	A	11.3
2	A	13.0 
3	A	12.4
4	A	14.0 
5	A	13.6
6	A	12.3
7	B	17.7
8	B	17.2
9	B	19.8
10	B	16.7
11	B	18.6
12	C	14.2
13	C	15.6
14	C	14.5
15	C	15.0 
16	C	15.7
17	C	16.8
18	D	15.4
19	D	15.1
20	D	15.6
21	D	16.5
22	D	16.0 
23	D	15.2
24	E	17.0 
25	E	19.1
26	E	15.3
27	E	17.3
28	E	18.0 
;
title '2020ass03 trt = diet, Wt';
proc glm data = ass03;
class diet;
model wt = diet/xpx inverse solution;
lsmeans diet/stderr;
estimate 'diet A-B' diet  1 -1 0 0 0;
estimate 'diet A-C' diet  1 0 -1 0 0;
estimate 'diet A-D' diet  1 0 0 -1 0;
estimate 'diet A-E' diet  1 0 0 0 -1;
estimate 'diet B.D.-C.E.' diet 0 0.5 -0.5 0.5 -0.5;
run;
quit;

proc glm;
class diet;
model wt = diet;
lsmeans diet/stderr pdiff adjust = bon;
run;
 



