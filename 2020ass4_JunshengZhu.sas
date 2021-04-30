/* 2020ass4_JunshengZhu */

proc iml;
reset print;
x = {1		1	0	0	0	0,
	 1		1	0	0	0	0,
	 1		1	0	0	0	0,
	 1		1	0	0	0	0,
	 1		1	0	0	0	0,
	 1		1	0	0	0	0,
	 1		0	1	0	0	0,
	 1		0	1	0	0	0,
	 1		0	1	0	0	0,
	 1		0	1	0	0	0,
	 1		0	1	0	0	0,
	 1		0	0	1	0	0,
	 1		0	0	1	0	0,
 	 1		0	0	1	0	0,
	 1		0	0	1	0	0,
	 1		0	0	1	0	0,
	 1		0	0	1	0	0,
	 1		0	0	0	1	0,
	 1		0	0	0	1	0,
	 1		0	0	0	1	0,
	 1		0	0	0	1	0,
	 1		0	0	0	1	0,
	 1		0	0	0	1	0,
	 1		0	0	0	0	1,
	 1		0	0	0	0	1,
	 1		0	0	0	0	1,
	 1		0	0	0	0	1,
	 1		0	0	0	0	1};

y = {11.3,
13.0,
12.4,
14.0,
13.6,
12.3,
17.7,
17.2,
19.8,
16.7,
18.6,
14.2,
15.6,
14.5,
15.0,
15.7,
16.8,
15.4,
15.1,
15.6,
16.5,
16.0, 
15.2,
17.0, 
19.1,
15.3,
17.3,
18.0};

print" X transpose X";
xtx = x` * x;
print" X transpose Y";
xty = x` * y;

tss = y` * y;

invxtx = ginv(xtx);
btilde = invxtx * xty;

sumy = sum(y);
nobs = nrow(x);
rx = 5;
dfe = nobs - rx;
ybar = sumy/nobs;
cf = nobs * ybar *ybar;

ssr = btilde` * xty;
ssrm = ssr - cf;

msr = ssr/rx;
msrm = ssrm/(rx - 1);

sse = tss - ssr;
mse = sse/dfe;

/* SS for trts */
kp = { 0 	1 -1 0 0 0,
	   0 	1 0 -1 0 0,
	   0 	1 0 0 -1 0,
	   0 	1 0 0 0 -1};
kb = kp* btilde;
kinvk = kp * invxtx *kp`;
invkk = inv(kinvk);
sstrt = kb` * invkk * kb;
dftrt = nrow(kp);
mstrt = sstrt/dftrt;

/*boniferroni tablulate_t */
dfe = 23;
pr = 0.05;
prnominal = pr/10;
ttab = tinv(1- prnominal/2, dfe);
print ttab;

/* differences between different trts */
kp = {0 1 -1 0 0 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-2 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 1-2 statistically significant##########";
else print"##########trt 1-2 not statistically significant##########";


kp = {0 1 0 -1 0 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-3 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"########## trt 1-3  statistically significant##########";
else print" ##########trt 1-3  not statistically significant##########";

kp = {0 1 0 0 -1 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 1-4 statistically significant##########";
else print"##########trt 1-4 not statistically significant##########";

kp = {0 1 0 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 1-5 statistically significant##########";
else print"##########trt 1-5 not statistically significant##########";

kp = {0 0 1 -1 0 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-3 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 2-3  statistically significant##########";
else print"##########trt 2-3  not statistically significant##########";

kp = {0 0 1 0 -1 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 2-4 statistically significant##########";
else print"##########trt 2-4 not statistically significant##########";

kp = {0 0 1 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 2-5 statistically significant##########";
else print"##########trt 2-5 not statistically significant##########";

kp = {0 0 0 1 -1 0 };
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 3-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"########## trt 3-4 statistically significant##########";
else print"##########trt 3-4 not statistically significant##########";

kp = {0 0 0 1 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 3-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 3-5 statistically significant##########";
else print"##########trt 3-5 not statistically significant##########";

kp = {0 0 0 0 1 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 4-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"##########trt 4-5 statistically significant##########";
else print"##########trt 4-5 not statistically significant##########";


/* differences between company A and B */
kp = {0 0 1 -1 1 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2.4. vs 3.5. " kb se;
cal_t = kb/se;
tab_t = 3.104;
if cal_t > tab_t 
then;
print"##########trt 2.4. vs 3.5. statistically significant##########";
else print"##########trt 2.4. vs 3.5. not statistically significant##########";
quit;



data ass04;
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
title '2020ass04 trt = diet, Wt';
proc glm data = ass04;
class diet;
model wt = diet/xpx inverse solution;
lsmeans diet/stderr;
estimate 'diet A-B' diet  1 -1 0 0 0;
estimate 'diet A-C' diet  1 0 -1 0 0;
estimate 'diet A-D' diet  1 0 0 -1 0;
estimate 'diet A-E' diet  1 0 0 0 -1;
estimate 'diet B-C' diet  0 1 -1 0 0;
estimate 'diet B-D' diet  0 1 0 -1 0;
estimate 'diet B-E' diet  0 1 0 0 -1;
estimate 'diet C-D' diet  0 0 1 -1 0;
estimate 'diet C-E' diet  0 0 1 0 -1;
estimate 'diet D-E' diet  0 0 0 1 -1;
estimate 'diet B.D.-C.E.' diet 0 1  -1  1  -1;
run;
quit;
