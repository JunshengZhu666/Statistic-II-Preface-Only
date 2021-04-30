data ass10;
input Observation Trt Sow Sex$ Wt;
cards;
1  1 1  female	 21.8
2  1 1  male  	 23.4
3  1 2  female 	 24.3
4  1 2  male   	 25.1
5  1 3  female 	 19.4
6  1 3  male  	 23.1
7  1 4  female 	 21.3
8  1 4  male 	 23.1
9  1 5  female   22.5
10 1 5  male 	 24.4
11 1 6  female 	 22.2
12 1 6  male	 22.8
13 1 7  male 	 21.8
14 1 8  female   23.1
15 1 8  male 	 24.0
16 1 9  female   19.1
17 1 9  male     22.4
18 1 10 female   20.0
19 1 10 male     22.9
20 2 11 female   25.0
21 2 11 male     23.1
22 2 12 female   22.7
23 2 12 male     25.9
24 2 13 female 	 23.9
25 2 13 male 	 24.1
26 2 14 female   24.1
27 2 14 male	 24.1
28 2 15 female   22.8
29 2 15 male 	 24.7
30 2 16 female   23.7
31 2 16 male     25.2
32 2 17 female   24.2
33 2 17 male     24.9
34 2 18 female   23.4
35 2 18 male     22.3
36 2 19 female   25.2
37 2 19 male     25.6
38 2 20 female   25.0
39 2 20 male     24.3
40 3 21 female   25.0
41 3 21 male 	 23.5
42 3 23 female   23.5
43 3 23 male 	 24.7
44 3 24 female   23.6
45 3 24 male	 25.0
46 3 25 female   22.7
47 3 25 male 	 23.4
48 3 26 female   21.7
49 3 26 male     22.0
50 3 27 female   23.0
51 3 27 male     22.2
52 3 28 female   20.8
53 3 28 male     22.6
54 3 29 female   23.5
55 3 29 male 	 21.7
56 3 30 female   24.4
57 3 30 male	 23.9
;

proc glm;
title"only for d.f.";
class trt sow sex;
model wt = trt sow(trt) sex trt*sex;
random sow(trt)/test;
run;

proc mixed;
title"FROM PROC MIXED";
class trt sow sex;
model wt = trt sex trt*sex/ ddfm = kr;
random sow(trt);
/*lsmeans*/
lsmeans trt*sex/pdiff adjust = bon;

run;

proc mixed;
title"FROM PROC MIXED - WITHOUT SOW EFFECT";
class trt sow sex;
model wt = trt sex trt*sex/ ddfm = kr;
run;


/*compute tabulated F and Chi values */
data tabulated;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
chitab = cinv(1 - prob, ndf);
cards;
26 . 0.025
26 . 0.975
25 . 0.025
25 . 0.975
;
proc print data = tabulated;
title'Tab_F and Tab_Chi^2';
var ndf ddf prob ftab chitab;
run;
