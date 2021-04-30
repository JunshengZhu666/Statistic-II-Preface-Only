/* 2020ass08_JunsehngZhu */
data ass08;
input Obs	Individual	plate	tissue	Expression;
cards;
1	1	1	1	21.3
2	1	1	2	22.3
3	1	1	3	20.5
4	1	1	4	22.2
5	2	2	1	16.9
6	2	2	2	21.7
7	2	2	3	22.0
8	2	2	4	24.1
9	3	3	1	22.5
10	3	3	2	21.3
11	3	3	3	22.1
12	3	3	4	23.4
13	4	4	1	19.0
14	4	4	2	21.5
15	4	4	3	19.3
16	4	4	4	21.2
17	5	5	1	21.8
18	5	5	2	22.4
19	5	5	3	19.3
20	5	5	4	21.0
21	6	6	1	20.6
22	6	6	2	20.3
23	6	6	3	21.0
24	6	6	4	22.9
25	7	7	1	21.5
26	7	7	2	21.9
27	7	7	3	27.0
28	7	7	4	20.2
29	8	8	1	18.9
30	8	8	2	19.6
31	8	8	3	17.5
32	8	8	4	20.2
33	9	9	1	22.4
34	9	9	2	23.6
35	9	9	3	22.1
36	9	9	4	22.9
37	10	10	1	20.8
38	10	10	2	20.7
39	10	10	3	22.5
40	10	10	4	20.1
41	11	11	1	20.9
42	11	11	2	23.6
43	11	11	3	23.1
44	11	11	4	22.8
45	12	12	1	17.4
46	12	12	2	20.1
47	12	12	3	19.4
48	12	12	4	19.4
49	13	13	1	22.2
50	13	13	2	23.5
51	13	13	3	22.3
52	13	13	4	21.9
53	14	14	1	19.1
54	14	14	2	20.0
55	14	14	3	21.5
56	14	14	4	20.2
57	15	15	1	18.4
58	15	15	2	21.3
59	15	15	3	19.4
60	15	15	4	18.9
61	16	16	1	15.5
62	16	16	2	19.0
63	16	16	3	19.3
64	16	16	4	20.8
65	17	17	1	21.3
66	17	17	2	21.3
67	17	17	3	23.7
68	17	17	4	21.6
69	18	18	1	21.6
70	18	18	2	21.6
71	18	18	3	20.0
72	18	18	4	23.1
73	19	19	1	20.8
74	19	19	2	22.3
75	19	19	3	23.0
76	19	19	4	21.7
77	20	20	1	21.9
78	20	20	2	20.1
79	20	20	3	21.4
80	20	20	4	22.2
;
/* using GLM to clearify d.f. */
proc glm data = ass08;
title'GLM';
class tissue plate;
model expression = tissue plate/xpx;
random plate/test;
lsmeans tissue/stderr pdiff adjust = bon;
lsmeans plate/stderr pdiff adjust = bon;
run;
quit;

proc mixed data = ass08;
title'MIXED';
class tissue plate;
model expression = tissue/ddfm = kr;
random plate/solution;
lsmeans tissue/pdiff adjust = bon;
run;
quit;

proc mixed data = ass08;
title'MIXED without random(plate)';
class tissue plate;
model expression = tissue/ddfm = kr;
run;
quit;

/*compute tabulated F and Chi values */
data tabulated;
input ndf ddf prob;
ftab = finv(1 - prob, ndf, ddf);
chitab = cinv(1 - prob, ndf);
cards;
3 57 0.05
19 57 0.05
19 . 0.025
19 . 0.975
57 . 0.025
57 . 0.975
;
proc print data = tabulated;
var ndf ddf prob ftab chitab;
run;
