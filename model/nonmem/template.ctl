$PROBLEM Model 216: SAEM fixed covs

$INPUT C ID TIME MDV EVID AMT DV NUM DOSE SEX WT
       X1 X2 X3 X4 X5 X6 CLI V1I QI V2I SIGMAI

$DATA ../template.csv
IGNORE=(C='C')
IGNORE(ID>250) ;; only fit on 250 subjects

$SUBROUTINE ADVAN3 TRANS4

$PK

;; --------- Fixed allometric scaling ------------
CLWT = LOG(WT/70)*THETA(5)
V1WT = LOG(WT/70)*THETA(6)
QWT  = LOG(WT/70)*THETA(5)
V2WT = LOG(WT/70)*THETA(6)

;; -------- Full covariate for CL ------------
CLSEX = SEX*THETA(7)
CLX1 = X1*THETA(8)
CLX2 = X2*THETA(9)
CLX3 = X3*THETA(10)
CLX4 = X4*THETA(11)
CLX5 = X5*THETA(12)
CLX6 = X6*THETA(13)

CLCOV = CLWT + CLSEX + CLX1 + CLX2 + CLX3 + CLX4 + CLX5 + CLX6

;; -------- PK Params ---------------------
TVV1 = EXP(THETA(1)+V1WT)
V1   = EXP(THETA(1)+V1WT+ETA(1))

TVCL = EXP(THETA(2)+CLCOV)
CL   = EXP(THETA(2)+CLCOV+ETA(2))

V2   = EXP(THETA(3)+V2WT)
Q    = EXP(THETA(4)+QWT) 

S1 = V1/1000 ; dose in mg, conc in mg/mL

$ERROR
IPRED = F 
Y = F * EXP(EPS(1))

$THETA
(4)       ;  1 V1 (L) - 60
(0.75)    ;  2 CL (L/hr) - 3.5
(4)       ;  3 V2 (L) - 70
(0.75)    ;  4 Q  (L/hr) - 4

0.75 FIX  ;  5  WT~CL ()
1 FIX     ;  6  WT~V ()
(0.1)     ;  7  SEX~CL ()
(0 FIX)   ;  8  X1~CL () ;; unable to converge when included
(-0.5)    ;  9  X2~CL ()
(0.2)     ;  10 X3~CL ()
(0.1)     ;  11 X4~CL ()
(0 FIX)   ;  12 X5~CL () ;; unable to converge when included
(0 FIX)   ;  13 X6~CL () ;; unable to converge when included

$OMEGA BLOCK(2)
0.1   ;ETA(V1)
0.01  0.1   ;ETA(CL)

$SIGMA
0.05     ; 1 exp error

$EST METHOD=SAEM INTERACTION NBURN=2000 NITER=1000 
     RANMETHOD=P PRINT=10 FILE=template_saem.ext
$EST METHOD=IMP EONLY=1 ISAMPLE=1000 NITER=5 MAPITER=0 
     RANMETHOD=P PRINT=1 FILE=template.ext

$COV PRINT=E RANMETHOD=P UNCONDITIONAL MATRIX=R

$TABLE NUM IPRED NPDE CWRES NOPRINT ONEHEADER RANMETHOD=P 
       FILE=template.tab FORMAT=s1PE14.7
$TABLE NUM CL TVCL V1 TVV1 Q V2 ETAS(1:LAST) NOAPPEND NOPRINT ONEHEADER 
       FILE=templatepar.tab FORMAT=s1PE14.7
