proc iml;
reset print;

x = {1 3.05 1.45 5.67,
	 1 4.22 1.35 4.86,
	 1 3.34 0.26 4.19,
	 1 3.77 0.23 4.42,
	 1 3.52 1.10 3.17,
	 1 3.54 0.76 2.76,
	 1 3.74 1.59 3.81,
	 1 3.78 0.39 3.23,
	 1 2.92 0.39 5.44,
	 1 3.10 0.64 6.16,
	 1 2.86 0.82 5.48,
	 1 2.78 0.64 4.62,
	 1 2.22 0.85 4.49,
	 1 2.67 0.90 5.59,
	 1 3.12 0.92 5.86,
	 1 3.03 0.97 6.60,
	 1 2.45 0.18 4.51,
	 1 4.12 0.62 5.31,
	 1 4.61 0.51 5.16,
	 1 3.94 0.45 4.45,
	 1 4.12 1.79 6.17,
	 1 2.93 0.25 3.38,
	 1 2.66 0.31 3.51,
	 1 3.17 0.20 3.08,
	 1 2.79 0.24 3.98,
	 1 2.61 0.20 3.64,
	 1 3.74 2.27 6.50,
	 1 3.13 1.48 4.28,
	 1 3.49 0.25 4.71,
	 1 2.94 2.22 4.58};

y = {0.34,        
	 0.11,        
	 0.38,        
	 0.68,        
	 0.18,        
 	 0.00,       
	 0.08,        
	 0.11,        
	 1.53,        
	 0.77, 
	 1.17,        
	 1.01,       
	 0.89,        
	 1.40,        
	 1.05,        
	 1.15,        
	 1.49,        
	 0.51,        
	 0.18,       
	 0.34,        
 	 0.36,        
 	 0.89,        
	 0.91,       
	 0.92,        
	 1.35,       
	 1.33,        
	 0.23,    
	 0.26,     
	 0.73,    
 	 0.23};

xtx = x` * x;
xty = x` * y;

invxtx = inv(xtx);
bhat = invxtx * xty;


sumy = sum(y);
tss = y` * y;
ssr = bhat` * xty;
sse = tss - ssr;

nobs = nrow(x);	
rx = 4 ;
dfe = nobs - rx;

sumy = sum(y);
ybar = sumy/nobs;
cf = nobs * ybar * ybar;
ssrm = ssr - cf;

mse = sse/dfe; 
covb = invxtx * mse;

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

yhat = x * bhat;
ehat = y - yhat;

/* SS b1 */
kp = {0 1 0 0};
kb = kp * bhat;
kinvk = kp * invxtx * kp`;
invkk = inv(kinvk);
ss1 = kb` * invkk * kb;
print "F calc for b1";
f1 = ss1/mse;

/* SS b2 */
kp = {0 0 1 0};
kb = kp * bhat;
kinvk = kp * invxtx * kp`;
invkk = inv(kinvk);
ss2 = kb` * invkk * kb;
print "F calc for b2";
f2 = ss2/mse;

/* SS b3 */
kp = {0 0 0 1};
kb = kp * bhat;
kinvk = kp * invxtx * kp`;
invkk = inv(kinvk);
ss3 = kb` * invkk * kb;
print "F calc for b3";
f3 = ss3/mse;

ndf = 4;
ddf = 26;
probability = 0.01;
ftab = finv(1 - probability, ndf, ddf);

quit;

/* use glm to plot */
data reg1;
input obs x1 x2 x3 y;
cards;
1 3.05 1.45 5.67 0.34
2 4.22 1.35 4.86 0.11
3 3.34 0.26 4.19 0.38
4 3.77 0.23 4.42 0.68
5 3.52 1.10 3.17 0.18
6 3.54 0.76 2.76 0.00
7 3.74 1.59 3.81 0.08
8 3.78 0.39 3.23 0.11
9 2.92 0.39 5.44 1.53
10 3.10 0.64 6.16 0.77
11 2.86 0.82 5.48 1.17
12 2.78 0.64 4.62 1.01
13 2.22 0.85 4.49 0.89
14 2.67 0.90 5.59 1.40
15 3.12 0.92 5.86 1.05
16 3.03 0.97 6.60 1.15
17 2.45 0.18 4.51 1.49
18 4.12 0.62 5.31 0.51
19 4.61 0.51 5.16 0.18
20 3.94 0.45 4.45 0.34
21 4.12 1.79 6.17 0.36
22 2.93 0.25 3.38 0.89
23 2.66 0.31 3.51 0.91
24 3.17 0.20 3.08 0.92
25 2.79 0.24 3.98 1.35
26 2.61 0.20 3.64 1.33
27 3.74 2.27 6.50 0.23
28 3.13 1.48 4.28 0.26
29 3.49 0.25 4.71 0.73
30 2.94 2.22 4.58 0.23
;
proc glm data=reg1;
model y = x1 x2 x3/XPX inverse solution;
output out=reg1out p=yhat r=ehat stdp=se;
run;
quit;

proc gplot data= reg1out;
plot ehat*x1;
plot ehat*x2;
plot ehat*x3;
plot ehat*yhat;
run;


