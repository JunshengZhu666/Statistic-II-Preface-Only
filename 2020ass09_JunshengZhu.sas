data ass09;
input observations factorA factorB y;
cards;
1 1 1 15.5
2 1 1 16.4
3 1 1 15.7
4 1 1 16.9
5 1 2 16.6
6 1 2 17.4
7 1 2 16.8
8 1 2 18.3
9 1 3 21.6
10 1 3 17.9
11 1 3 19.8
12 1 3 19.5
13 2 1 18.3
14 2 1 17.5
15 2 1 18.9
16 2 1 20.6
17 2 2 17.9
18 2 2 18.9
19 2 2 20.1
20 2 2 19.3
21 2 3 19.2
22 2 3 19.8
23 2 3 19.7
24 2 3 19.4
25 3 1 19.9
26 3 1 19.5
27 3 1 19.6
28 3 1 18.6
29 3 2 22.3
30 3 2 22.1
31 3 2 18.1
32 3 2 18.7
33 3 3 21.1
34 3 3 19.7
35 3 3 21.5
36 3 3 22.2
37 4 1 20.1
38 4 1 18.6
39 4 1 21.5
40 4 1 20.3
41 4 2 20.5
42 4 2 20.3
43 4 2 20.4
44 4 2 20
;

/*proc print data = ass09;
var observations factorA factorB y;
run;*/

ods output estimate=est1;
proc glm data=ass09;
class factorA factorB;
model y = factorA factorB factorA*factorB/xpx i;
estimate 'A1-A2' factorA 1 -1 0 0;
estimate 'A1-A3' factorA 1 0 -1 0;
estimate 'A2-A3' factorA 0 1 -1 0;
estimate 'B1-B2' factorB 1 -1 0;
lsmeans factorA/stderr pdiff adjust = bon;
lsmeans factorB/stderr pdiff;
run;
quit;


/*compute tabulated F and Chi values */
data tabulated;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
chitab = cinv(1 - prob, ndf);
cards;
10 33 0.05
3 33 0.05
2 33 0.05
5 33 0.05
33 . 0.025
33 . 0.975
;
proc print data = tabulated;
title'Tab_F and Tab_Chi^2';
var ndf ddf prob ftab chitab;
run;
