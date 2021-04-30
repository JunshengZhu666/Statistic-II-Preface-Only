/* 2020ass5_JunshengZhu */

/* Hail Dr.Roger */

proc iml;
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
kp = { 0 1 -1 0 0 0 0,
	   0 1 0 -1 0 0 0,
	   0 1 0 0 -1 0 0,
	   0 1 0 0 0 -1 0,
	   0 1 0 0 0 0 -1};
kb = kp* btilde;
kinvk = kp * invxtx *kp`;
invkk = inv(kinvk);
sstrt = kb` * invkk * kb;
dftrt = nrow(kp);
mstrt = sstrt/dftrt;

/*bonferroni tablulate_t */
dfe = 29;
pr = 0.01;
prnominal = pr/15;
ttab = tinv(1- prnominal/2, dfe);
print ttab;

/* differences between different trts */
kp = {0     1 -1 0 0 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-2 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 1-2  significant==========";
else print"==========trt 1-2 not  significant==========";


kp = {0     1 0 -1 0 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-3 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"========== trt 1-3   significant==========";
else print" ==========trt 1-3  not  significant==========";

kp = {0     1 0 0 -1 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 1-4  significant==========";
else print"==========trt 1-4 not  significant==========";

kp = {0     1 0 0 0 -1 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 1-5  significant==========";
else print"==========trt 1-5 not  significant==========";

kp = {0     1 0 0 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 1-6 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 1-6  significant==========";
else print"==========trt 1-6 not  significant==========";

kp = {0     0 1 -1 0 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-3 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 2-3   significant==========";
else print"==========trt 2-3  not  significant==========";

kp = {0     0 1 0 -1 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 2-4  significant==========";
else print"==========trt 2-4 not  significant==========";

kp = {0     0 1 0 0 -1 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 2-5  significant==========";
else print"==========trt 2-5 not  significant==========";

kp = {0     0 1 0 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 2-6 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 2-6  significant==========";
else print"==========trt 2-6 not  significant==========";

kp = {0     0 0 1 -1 0 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 3-4 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"========== trt 3-4  significant==========";
else print"==========trt 3-4 not  significant==========";

kp = {0     0 0 1 0 -1 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 3-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 3-5  significant==========";
else print"==========trt 3-5 not  significant==========";

kp = {0     0 0 1 0 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 3-6 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 3-6  significant==========";
else print"==========trt 3-6 not  significant==========";

kp = {0     0 0 0 1 -1 0};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 4-5 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 4-5  significant==========";
else print"==========trt 4-5 not  significant==========";

kp = {0     0 0 0 1 0 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 4-6 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 4-6  significant==========";
else print"==========trt 4-6 not  significant==========";

kp = {0     0 0 0 0 1 -1};
kb = kp * btilde;
sv = kp * invxtx * kp` * mse;
se  = sqrt(sv);
print" trt 5-6 " kb se;
cal_t = kb/se;
tab_t = ttab;
if abs(cal_t) > tab_t 
then;
print"==========trt 5-6  significant==========";
else print"==========trt 5-6 not  significant==========";


k1 = { 1,	1, 0, 0, 0, 0, 0};
kb = k1` * btilde;
kgk = k1` * invxtx * k1;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt1 =' kb '+-' se;

k2 = { 1,	0, 1, 0, 0, 0, 0};
kb = k2` * btilde;
kgk = k2` * invxtx * k2;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt2 =' kb '+-' se;

k3 = { 1,	0, 0, 1, 0, 0, 0};
kb = k3` * btilde;
kgk = k3` * invxtx * k3;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt3 =' kb '+-' se;

k4 = { 1,	0, 0, 0, 1, 0, 0};
kb = k4` * btilde;
kgk = k4` * invxtx * k4;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt4 =' kb '+-' se;

k5 = { 1,	0, 0, 0, 0, 1, 0};
kb = k5` * btilde;
kgk = k5` * invxtx * k5;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt5 =' kb '+-' se;

k6 = { 1,	0, 0, 0, 0, 0, 1};
kb = k6` * btilde;
kgk = k6` * invxtx * k6;
sv = kgk * mse;
se = sqrt(sv);
print 'estimate of mu + trt6 =' kb '+-' se;




quit;





/*compute tabulated F values */
data Fvalues;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
cards;
1 29 0.01
3 29 0.01
5 29 0.01
6 29 0.01
;
proc print data = Fvalues;
var ndf ddf prob ftab;
run;




data ass05;
input plots trt$ wt orgma;
cards;
1	A	15.9	1
2	A	15.6	1
3	A	15.4	1
4	A	14.0 	1
5	A	17.7	1
6	A	15.6	1
7	B	18.7	2
8	B	18.9	2
9	B	20.4	2
10	B	19.8	2
11	B	19.3	2
12	B	19.3	2
13	C	21.8	3
14	C	20.5	3
15	C	18.7	3
16	C	20.8	3
17	C	20.6	3
18	C	19.7	3
19	D	19.9	4
20	D	21.3	4
21	D	21.4	4
22	D	17.1	4
23	D	20.3	4
24	D	19.8	4
25	E	20.2	5
26	E	20.6	5
27	E	19.7	5
28	E	21.2	5
29	E	21.2	5
30	F	14.3	6
31	F	15.3	6
32	F	15.6	6
33	F	14.1	6
34	F	12.9	6
35	F	15.0 	6
;
title '2020ass05 trt = trts, Wt';
proc glm data = ass05;
class trt;
model wt = trt/xpx inverse solution;
lsmeans trt/stderr;
estimate 'trt A-B' trt  1 -1 0 0 0 0;
estimate 'trt A-C' trt  1 0 -1 0 0 0;
estimate 'trt A-D' trt  1 0 0 -1 0 0;
estimate 'trt A-E' trt  1 0 0 0 -1 0;
estimate 'trt A-F' trt  1 0 0 0 0 -1;
estimate 'trt B-C' trt  0 1 -1 0 0 0;
estimate 'trt B-D' trt  0 1 0 -1 0 0;
estimate 'trt B-E' trt  0 1 0 0 -1 0;
estimate 'trt B-F' trt  0 1 0 0 0 -1;
estimate 'trt C-D' trt  0 0 1 -1 0 0;
estimate 'trt C-E' trt  0 0 1 0 -1 0;
estimate 'trt C-F' trt  0 0 1 0 0 -1;
estimate 'trt D-E' trt  0 0 0 1 -1 0;
estimate 'trt D-F' trt  0 0 0 1 0 -1;
estimate 'trt E-F' trt  0 0 0 0 1 -1;
run;
quit;

proc glm data = ass05;
class trt;
model wt = orgma orgma*orgma trt;
run;




