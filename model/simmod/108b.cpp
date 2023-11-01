[ SET ] delta = 0.5

[ PROB ] 
- Date: `r Sys.Date()`

This model requires mrgsolve >= 1.0.3

2 cpt model, with a few complicated covariate relationship

[ PLUGIN ] autodec nm-vars

[ PKMODEL ] cmt = "GUT CENT PERIPH", depot = TRUE

[ PARAM ]
WT = 70.0
SEX = 0.0
X1 = 0.0
X2 = 0.0
X3 = 0.0
X4 = 0.0
X5 = 0.0
X6 = 0.0
THETA1 = 0.0
THETA2 = 0.0
THETA3 = 0.0
THETA4 = 0.0
THETA5 = 0.0
THETA6 = 0.0
THETA7 = 0.0
THETA8 = 0.0
THETA9 = 0.0
THETA10 = 0.0


[ PK ]
double V2WT   = log(WT/70.0);
double CLWT   = log(WT/70.0)*0.75;
double CLEFF1 = exp(X1/10.0)*THETA6; 
double CLEFF2  = (X2 > 0.0) ?  X2*THETA7 : 0.0 ;
double CLEFF3  = X3 * THETA8 + X3 * SEX * THETA9 ;
double CLEFF4 = X4 * THETA10;
double V3WT   = log(WT/70.0);
double QWT    = log(WT/70.0)*0.75;

double KA     = exp(THETA1 + ETA(1));
double V2     = exp(THETA2 + V2WT + ETA(2));
double TVLCL  = THETA3 + CLEFF1 + CLEFF2 + CLEFF3 + CLEFF4;
double CL     = exp(TVLCL + ETA(3));
double V3     = exp(THETA4 + V3WT);
double Q      = exp(THETA5 + QWT);  
double S2     = V2/1000.0; //; dose in mcg, conc in mcg/mL

[OMEGA] 
0.1

[OMEGA] @block
0.04 0.005 0.04 

[SIGMA]
0.01

[ ERROR ] 
F = CENT/S2;
//Y = F*(1 + EPS(1)); // this gives negative values at low concentrations
Y = F* exp(EPS(1)); 
IPRED = F; 

[ CAPTURE ]
TVLCL CL Q V2 V3 KA IPRED Y //ETA(1) ETA(2) ETA(3)
