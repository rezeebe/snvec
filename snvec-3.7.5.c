#define VER "snvec.c VERSION: 3.7.5 12/2021"
/****************************************************************
 
   file: snvec.c

   Richard E. Zeebe
   School of Ocean and Earth 
   Science and Technology
   University of Hawaii at Manoa
   1000 Pope Road, MSB 504
   Honolulu, HI 96822, USA
   email: zeebe@soest.hawaii.edu

   gcc -Wall -o snvec.x snvec-3.7.5.c -lm
   optimize: -Ofast -O2 -O3
   check values:
   $ time ./snvec.x -1.e3 1.0 0.0

   purpose:
   ODE solver to integrate vec{s}:
       calc obliquity & precession (Ward74,79,Bills90).
   this version contains one driver:
   (1)  driver() C standalone version

   updates:

   12/19/21 revised epl0, phi0, cp0
   12/01/21 added FGP factor
   11/26/21 updated test 1/(1-e^2)^3/2 term
   11/16/21 read input args. Ed, Td MACROS => globals (speed: same). 
   11/12/21 resumed. abvec => snvec
   03/04/21 added ClimPrec. vars: eei, lphi, lphu, cp.
		    check cp output. OK. see prectilt.m flag12.
   02/04/21 TD. SI units (AU,GM,R0). rm old constants
   02/03/21 V3.0 rm QUINN, mdriver, bvec > derivs (see V2.9).    
   02/01/21 qinterp(): abs -> fabs!!! (minor). kmax scaled.
   01/31/21 checked numerical noise. OK, see PTMan.
   01/12/21 euler(), unwrap(). precession OK.
   01/02/21 qinterp(). quicker! smoother?
   12/29/20 resumed, fixed #def NB, freadHNB, fvei.
   05/11/20 main() => driver(). main only calls driver().
            double **abvec calls solver and returns [t,yy].
   
   check values:    
   $ time ./snvec.x -1.e3 1.0 0.0

    Ubuntu 20.04.3 LTS Linux 5.10.0-1052, gcc 9.3.0
    Ubuntu 16.04.7 LTS Linux 4.4.0-210,   gcc 5.4.0   

    WARD79 KMNEW QINTERP fabs spin[001] SI(AU,GM,R0) kmax...*2.656 FGP NEW0
    0.404360548210294 -0.053622586775944  0.913026378223150	
    0.404360548209264 -0.053622586779551  0.913026378223374	

    @ Final values obliquity, precession (rad):
    0.413049154286417 -0.561907144352942
    0.413049154286600 -0.561907144362299		
	   
    orbital elements:
    [1] t   : time
    [2] aa  : SemiMajorAxis
    [3] ee  : Eccentricity
    [4] inc : Inclination
    [5] lph : LongPerihelion
    [6] lan : LongAscNode
    [7] arp : ArgPerihelion
    [8] mna : MeanAnomaly

    Note: [lph,lan] = [-PI PI] (jump). Hence unwrap
    lph, lan before interpolating (lphu, lphi & lanu, lani)

   counters:
   ls  : length OrbitSolution (read-in)
   lt  : length PTsolution output values (lt = kount)
   kmax: length requested output values, lt ~< kmax

****************************************************************/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "fun/utils.c"

#define NEQ 3   /* #Equations = dim(svec) */
#define WRT 1   /* Write output 0/1       */
#define WRTB 1  /* Write output 0/1 bin   */
#define QINTERP /* Quick interpolation    */

/* some function declarations */
void derivs(double t, double *y, double *yp);
void jacobn(double x, double *y, double *dfdx, double **dfdy, int n);

/* solver control */
#define TEND   (-1.e3)   /* default time final kyr   */
#define EPSLVR (1.e-7)   /* accuracy 1e-7 2.2e-7/8.5e-7 La  */
#define CSCAL  (1.0)     /* error scaling, set !?!   */
#define MAXSTP 100000000 /* max # of steps 100000000 */
int kmax,kount;          /* max # output values,        see solver */
double *tmv,**yy,dxsav;  /* output vars, save interval, see solver */
int nok,nbad;
double h1,hmin;
/* include solver routines */
#include "fun/solver.c"

/* read HNBody solution: allocate max matrix dimensions */
#define NS 8            /* 8 elements incl time         */
#define MS 7800000      /* max 8*7,800,000*16 =~ 1GB    */
/* include routine to read HNBody orbital elements      */
#include "fun/freadHNB.c" 

/* constants */
#define AU   (1.49597870700e11) /* m              */   
#define GM   (1.32712440041e20) /* m3/s2          */
#define OM   (7.292115e-5)      /* 1/s EarthRot   */
#define R0   (3.8440e8)         /* m Moon R0      */
#define GK   (0.9925194)        /* Kinoshita75,77 */
#define ED0  (0.0032738134)     /* DynEll (C-A)/C */
#define FGP  (0.99961908)
#define AU3  (AU*AU*AU)
#define R03  (R0*R0*R0)
#define PI    3.1415926535897932         /* Pi */
#define R2D  (180./PI)       /* radians to deg */

#define D2S  (3600.*24.)
#define Y2D  (365.25)
#define KY2D (1.e3*Y2D)

/* set default Ed, Td */
#define ED  (1.0000)         /* set factor 1.0 */
#define TD  (0.0000)         /* set factor 0.0 */

/* mass ratios */
#define MSEL (328900.5596)       /* MS/(ME+ML) */
#define MEL  (81.300568)         /* ME/ML      */
#define MLS  (1./(MSEL*(1+MEL))) /* ML/MS      */
/* K0, beta0 for torques */
#define K0   ((3./2.)*GM*ED0*ED/(OM*AU3))
#define K0D  (K0*D2S)            /* 1/s => 1/d */
#define BET0 (GK*MLS*AU3/R03)
#define K0B0 (K0D*(1.+BET0))
/* Moon mean motion */
#define N0   sqrt(GM/MSEL/R03)
#define NW0 (N0/OM) /* ratio (n/om)_0 */
/* Tidal dissipation Quinn91 Eqs. (3, 11) */
/* NDN = (dndt/n)_0, WDW = (domdt/om)_0 */
#define NDN (-4.6e-18*D2S*TD)   /* 1/s => 1/d */ 
#define WDW (51.*NDN*NW0)       /* Lambeck80  */
/* tidal effect on obliquity */
#define UEPSDOT (-4.17e-19)

/* SunRot Angles (Transform to HCI) */
#define OMT  (75.5940)
#define INCT (7.155)
#define EP0  (23.439291111111110)  /* Obliquity t0 */

/* a few globals */
int ls;
/* needed in derivs() */
double dts,*ts,*hh,*kk,*pp,*qq,*cc,*dd,**nn;
double ndn,wdw,k0d,k0b0,tdg;

/*===================== qinterp() ============================*/
/*
 quick interpolation. interpolates OrbSoln var y (dt=ds) at
 time step m to time t = ts[m]+dx. 
 */
double qinterp(double *y, double ds, double dx, int m)
{
double yi=y[m],dy=0.0,dsa=fabs(ds),dxa=fabs(dx);
int mm=1;
	
if(dxa > DBL_EPSILON){
 mm = m - (int)(copysign(1.0,dx));
 dy = y[mm] - y[m];
 yi += dy*dxa/dsa;
}
  
return(yi);	

}
/*===================== qinterp() END ========================*/

/*====================== unwrap() ============================*/
/*
 unwrap angle. maps jumps greater than pi to their 2pi complement.
 */
void unwrap(double *yu, double *y, int n)
{
int i;
double dy,c=0.0,*cv;

cv = dvector(1,n);	

/* unwrap */
cv[1] = 0.0;	
for(i=2;i<=n;i++){
  dy = (y[i]-y[i-1])/R2D; 
    if(dy > PI){
      c = c - 2.*PI;
    } else
	if(dy < -PI){
      c = c + 2.*PI;		
    }
  cv[i] = c;
}
for(i=1;i<=n;i++)
  yu[i] = y[i] + cv[i]*R2D;

free_dvector(cv,1,n);

}
/*====================== unwrap() END ========================*/

/*======================= euler() ============================*/
/*
 Euler transformation. s* = A * s, where spin vector s is in 
 invariable plane and s* in instant orbit plane. inv = 1 gives
 inverse transformation (A^-1 = A' = transpose(A)).
 */
void euler(double *vp, double *v, double inc, double lan, int inv)
{
int k;
double a1[4]={0,          cos(lan),           sin(lan), 0.0     };
double a2[4]={0,-cos(inc)*sin(lan),  cos(inc)*cos(lan), sin(inc)};
double a3[4]={0, sin(inc)*sin(lan), -sin(inc)*cos(lan), cos(inc)};
double  u[4]={0,v[1],v[2],v[3]};

vp[1] = vp[2] = vp[3] = 0.0;

/* inverse = transpose */
if(inv){
  ftransp(a1,a2,a3);
}

for(k=1;k<=3;k++){
  vp[1] += a1[k]*u[k];
  vp[2] += a2[k]*u[k];
  vp[3] += a3[k]*u[k];
}
}
/*======================= euler() END =========================*/

/*====================== fvei() ==============================*/
/*
 calculates global h,k,p,q etc. from ecc,inc etc.
 */
void fvei(double *ee, double *inc, double *lph, double *lan, int ls)
{
int i;

for(i=1;i<=ls;i++){
 hh[i] = ee[i]*sin(lph[i]/R2D); 
 kk[i] = ee[i]*cos(lph[i]/R2D);
 pp[i] = 2.*sin(0.5*inc[i]/R2D)*sin(lan[i]/R2D);
 qq[i] = 2.*sin(0.5*inc[i]/R2D)*cos(lan[i]/R2D);
 cc[i] = cos(inc[i]/R2D);
 dd[i] = cos(inc[i]/R2D/2.);
 /* nn = nvec(t): normal to orbit */
 nn[1][i] =  sin(inc[i]/R2D)*sin(lan[i]/R2D); 
 nn[2][i] = -sin(inc[i]/R2D)*cos(lan[i]/R2D);
 nn[3][i] =  cos(inc[i]/R2D);		
}
}
/*====================== fvei() END ==========================*/

/*====================== finargs() ===========================*/
/*
 parse input arguments. arg list:
 [1] tend
 [2] Ed
 [3] Td
 [4] dir  OrbitSoln
 [5] file OrbitSoln
 */
void finargs(int argc, char **argv, double *argd, char *dir, char *foo)
{
char mssg[BUFSIZ];

if(argc < 3+1){ 
 sprintf(mssg,"finargs(): too few input arguments.\
  provide tend[-x kyr] Ed Td.\n           using defaults:\
  %.4e %.6f %.6f",TEND,ED,TD);
 fwarn(mssg);
 argd[1] = TEND;
 argd[2] = ED;
 argd[3] = TD;
}else{
 argd[1] = strtod(argv[1],NULL);
 argd[2] = strtod(argv[2],NULL);
 argd[3] = strtod(argv[3],NULL);
}

/* dir, file HNBody Elements  */
if(argc < 5+1){
 strcpy(dir,"/dat/PrecTilt/OS/ZB18a");
 strcpy(foo,"ems-plan3.dat");
 sprintf(mssg,"finargs(): no input argument(s) for dir & file\
  of orbital solution.\n  provide tend[-x kyr] Ed Td dir foo.\
  using defaults:\n\
  dir=%s \n  foo=%s",dir,foo);
 fwarn(mssg);
}else{
 strcpy(dir,argv[4]);	
 strcpy(foo,argv[5]);	
}
}
/*====================== finargs() END =======================*/

/*====================== fedtd() =============================*/
/*
 calculates global vars ndn,wdw,k0d from Td,Ed
 */
void fedtd(double ed, double td)
{
/* K0, beta0 for torques */	
k0d  = ((3./2.)*GM*ED0*ed/(OM*AU3))*D2S; /* 1/s => 1/d */
k0b0 = k0d*(1.+BET0);	
/* Tidal dissipation Quinn91 Eqs. (3, 11) */
/* NDN = (dndt/n)_0, WDW = (domdt/om)_0   */
ndn = -4.6e-18*D2S*td;        /* 1/s => 1/d           */ 
wdw =  51.*ndn*NW0;           /* Lambeck80, see PTman */
tdg = td;                     /* global Td */
}
/*====================== fedtd() END =========================*/

/*====================== finits() ============================*/
/*
 init spin vector, transform to HCI
 s,n in HCI. s',n' in ECLIPJ2000
 */
void finits(double *s0, double **nn, double ep0, double inct, double omt)
{
int i;	
double a,b,c,cs=cos(ep0),n[3+1],np[3+1],s0p[3+1];

/* orbit normal at t=0 */	
for(i=1;i<=3;i++)
  n[i] = nn[i][1];

/* transform n => n' */
euler(np,n,inct,omt,1);

/* solve quadratic equation for s0'y */
a = np[2]*np[2] + np[3]*np[3];
b = -2.*cs*np[2];
c = cs*cs - np[3]*np[3];

s0p[2] = (-b + sqrt(b*b-4.*a*c))/(2.*a);
s0p[3] = sqrt(1. - s0p[2]*s0p[2]);
s0p[1] = 0.;

/* transform s0' => s0 */
euler(s0,s0p,inct,omt,0);

}
/*====================== finits() END ========================*/

/*====================== derivs() ============================*/
/*
 derivatives. RHS of DEQs for spin vector s = y
 */
void derivs(double t, double *y, double *yp)
{
int m;
double fac,sx=y[1],sy=y[2],sz=y[3],kb,ff;
#ifdef EPSDOT
double *s=y,dotab,tmp;
#endif	

/* K0, beta0 changing with Td, Ed */
kb = k0d*(1.+1.*wdw*t)*(1. + BET0*(1.+2.*ndn*t));

/* set time index of solution */
m = intmin((int)(round(fabs(t/dts))+1),ls);

#ifdef QINTERP
double dx,qqi,ppi,cci,ddi;	

dx = t-ts[m];
qqi = qinterp(qq,dts,dx,m);
ppi = qinterp(pp,dts,dx,m);
cci = qinterp(cc,dts,dx,m);
ddi = qinterp(dd,dts,dx,m);

/* 1/(1-e^2)^3/2 term */
ff  = (1.-hh[m]*hh[m]-kk[m]*kk[m]);
ff  = 1./sqrt(ff*ff*ff);
kb = k0d*(1.+1.*wdw*t)*(ff + BET0*(1.+2.*ndn*t));

fac = FGP*kb*(ddi*(ppi*sx - qqi*sy) + cci*sz);	

yp[1] =  fac*( cci*sy + ddi*qqi*sz);
yp[2] =  fac*(-cci*sx + ddi*ppi*sz);
yp[3] = -fac*( qqi*sx + ppi*sy)*ddi;

#else
fac = FGP*k0b0*(dd[m]*(pp[m]*sx - qq[m]*sy) + cc[m]*sz);

yp[1] =  fac*( cc[m]*sy + dd[m]*qq[m]*sz);
yp[2] =  fac*(-cc[m]*sx + dd[m]*pp[m]*sz);
yp[3] = -fac*( qq[m]*sx + pp[m]*sy)*dd[m];
#endif

#ifdef EPSDOT
/* orbit normal at ti = m is nn[j][m] */
/* dot product and dydt:  */
/* dotab = cos(epl), sqrt = sin(epl=acos(dotab)) */
dotab = s[1]*nn[1][m]+s[2]*nn[2][m]+s[3]*nn[3][m];
tmp = tdg*EPSDOT*D2S/sqrt(1.-dotab*dotab);

yp[1] += tmp*(nn[1][m] - dotab*s[1]);
yp[2] += tmp*(nn[2][m] - dotab*s[2]);
yp[3] += tmp*(nn[3][m] - dotab*s[3]);
#endif

}
/*==================== derivs() END ==========================*/

/*======================== driver() ==========================*/
/*
 driver routine solving DEQs for spin vector s = y.
 */
void driver(int argc, char **argv)
{
int i,k,m,lt;	
double t0,tfin,tfink,sckx,*nv,*y0,*u,*up,*epl,*phi,tmp,**nni,dx;
double *eei,*inci,*lphi,*lphu,*lani,*lanu,*cp;
double argd[4],tend,ed,td;
char dir[BUFSIZ],foo[BUFSIZ],mssg[BUFSIZ];
FILE *fpout;

printf("\n#===================================================# ");
printf("\nThis is %s",VER);
printf("\nRichard E. Zeebe\n");
	
/* parse input arguments */
finargs(argc,argv,argd,dir,foo);
tend = argd[1];
ed   = argd[2];
td   = argd[3];
printf("\n@ Integration parameters:");
printf("\n@ tend = %+.4e kyr",tend);
printf("\n@ Ed   = % .6f    ",ed);
printf("\n@ Td   = % .6f  \n",td);
/* calc global vars from Td, Ed */
fedtd(ed,td);

/* orbital elements (see header), ts is global */
double *aa,*ee,*inc,*lph,*lan,*arp,*mna;

/* clock vars */
clock_t tick_start, tick_end;
double wall_time;	

/* allocate, ls yet unknown */
ts  = dvector(1,MS);
aa  = dvector(1,MS);
ee  = dvector(1,MS);
inc = dvector(1,MS);
lph = dvector(1,MS);
lan = dvector(1,MS);
arp = dvector(1,MS);
mna = dvector(1,MS);

/* read HNBody Elements, get ls  */
freadHNB(dir,foo,&ls,ts,aa,ee,inc,lph,lan,arp,mna);
dts = ts[2]-ts[1];

/* allocate and calc global vars, nn(t), e-i vectors */
hh  = dvector(1,ls);
kk  = dvector(1,ls);
pp  = dvector(1,ls);
qq  = dvector(1,ls);
cc  = dvector(1,ls);
dd  = dvector(1,ls);
nn  = dmatrix(1,NEQ,1,ls);
fvei(ee,inc,lph,lan,ls);

/* init spin vector, transform to HCI */
double s0[3+1] = {0,0.,0.,1.};
finits(s0,nn,EP0/R2D,INCT/R2D,OMT/R2D);
	
/* integration parameters */
t0    = 0.0;
tfink = tend;       /* kyr */
tfin  = tfink*KY2D; /* kyr => days (negative!) */
sckx  = fabs(tfink/1.e3); /* scale factor kmax */

if(tend >= 0.0){ 
 sprintf(mssg,"driver(): expecting tend<0 but tend = %.2e kyr",tend);
 ferrx(mssg);
}
if(tend < ts[ls]/KY2D){ 
 sprintf(mssg,"driver(): tend<tend(OS). limit tend = %.2e kyr",ts[ls]/KY2D);
 ferrx(mssg);
}

/* solver control */
kmax  = (int)(floor(1000.*2.656*sckx)); /* *2.656 *10 *5 */
dxsav = (tfin-t0)/kmax;
h1    = 0.1*dxsav;
hmin  = 0.0;

/* allocate */
nv  = dvector(1,3); /* nvec, single row of nn */
y0  = dvector(1,3);
u   = dvector(1,3); 
up  = dvector(1,3); 
epl = dvector(1,kmax); 
phi = dvector(1,kmax); 
tmv = dvector(1,kmax); 
yy  = dmatrix(1,NEQ,1,kmax);
nni = dmatrix(1,NEQ,1,kmax);
eei  = dvector(1,kmax);
inci = dvector(1,kmax);
lphi = dvector(1,kmax);
lani = dvector(1,kmax);
lphu = dvector(1,ls);   /* 1,ls! lphu = unwrap(lph) */
lanu = dvector(1,ls);   /* 1,ls! lanu = unwrap(lan) */
cp   = dvector(1,kmax); /* climatic precession      */

/* start values */
for(k=1;k<=3;k++)
 y0[k] = s0[k];

/*%%%%%%%%%%%%%%%%%%%%%%%% solver %%%%%%%%%%%%%%%%%%%%%%%%%%%*/
tick_start = clock(); /* start clock */
odeint(y0,NEQ,t0,tfin,EPSLVR,h1,hmin,&nok,&nbad,derivs,stiff);
printf("\n");
tick_end  = clock(); /* end clock */
wall_time = (double)(tick_end-tick_start)/CLOCKS_PER_SEC;
printf("\n@ Wallclock Time = %f sec \n",wall_time);
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* print step info */
printf("\n@ Step Info:\n");
printf("     %6d steps OK\n",nok);
printf("     %6d steps FIXED\n",nbad);
printf("     %6d steps TOTAL\n",nok+nbad);
printf("  kmax=%d, kount=%d\n\n",kmax,kount);
if(kount < 2)
   printf("main(): kount < 2. increase tfinal?");
lt = kount;

/* final values */
printf("@ Final values s[1][2][3]; s-error = |s|-1:\n");
for(k=1;k<=NEQ;k++){
 printf("%18.15f ",yy[k][lt]);
 u[k] = yy[k][lt];
}
printf("\n %e",sqrt(vvdot(u,u,3))-1.0);
printf("\n");

/* unwrap lph, lan, interpolate various on tmv */
unwrap(lphu,lph,ls);
unwrap(lanu,lan,ls);
for(i=1;i<=kount;i++){
 m  = intmin((int)(round(fabs(tmv[i]/dts))+1),ls);
 dx = tmv[i]-ts[m];
 nni[1][i] = qinterp(nn[1],dts,dx,m);
 nni[2][i] = qinterp(nn[2],dts,dx,m);
 nni[3][i] = qinterp(nn[3],dts,dx,m);
    eei[i] = qinterp(ee,dts,dx,m);	
   inci[i] = qinterp(inc,dts,dx,m);	
   lphi[i] = qinterp(lphu,dts,dx,m);
   lani[i] = qinterp(lanu,dts,dx,m);
}

/* calc obliquity */
for(i=1;i<=kount;i++){
 for(k=1;k<=NEQ;k++){ 
    u[k] = yy[k][i];
   nv[k] = nni[k][i];
 }
 tmp = vvdot(u,nv,3);
 epl[i] = acos(tmp);
}

/* calc precession & climatic precession */
for(i=1;i<=kount;i++){
 for(k=1;k<=NEQ;k++){
   u[k] = yy[k][i];
  }
  /* coords: fixed HCI => moving orbit plane */
  euler(up,u,inci[i]/R2D,lani[i]/R2D,0);
  /* coords: relative to phi(t=0)=0 at J2000 */
  euler(up,up,0.0,-(lani[i]+OMT)/R2D+PI/2.,0);
  phi[i] = atan2(up[2],up[1]);
}
tmp = phi[1];	
for(i=1;i<=kount;i++){
  phi[i] -= tmp;
   cp[i]  = eei[i]*sin((lphi[i]+OMT)/R2D - phi[i]);
}

printf("\n@ Final values obliquity, precession (rad):\n");
printf("%18.15f %18.15f\n",epl[kount],phi[kount]);

/* write output to files */
if(WRT){
char frmt[BUFSIZ]="%22.15e ";
printf("\n@ Writing ascii output\n");	
fpout = fopen("out.dat","w");	 
for(i=1;i<=kount;i++){
 fprintf(fpout,frmt,tmv[i]/KY2D);
 if(0){
 for(k=1;k<=NEQ;k++)
  fprintf(fpout,frmt,yy[k][i]);
 }
 //fprintf(fpout,frmt,eei[i]);
 fprintf(fpout,frmt,epl[i]);
 fprintf(fpout,frmt,phi[i]);
 fprintf(fpout,frmt,cp[i]);
 //fprintf(fpout,frmt,lani[i]);
 fprintf(fpout,"\n");
}
fclose(fpout);
}

if(WRTB){
printf("\n@ Writing binary output\n");	
/* arrays: index 0 not used. pass address of 2nd element */
fpout = fopen("out.bin","wb");
fwrite(&kount,sizeof(int),1,fpout);
fwrite(tmv+1,sizeof(tmv),kount,fpout);
fwrite(epl+1,sizeof(epl),kount,fpout);
fwrite(phi+1,sizeof(phi),kount,fpout);
fwrite( cp+1,sizeof(cp) ,kount,fpout);
fclose(fpout);
}

free_dvector(ts,1,MS);
free_dvector(aa,1,MS);
free_dvector(ee,1,MS);
free_dvector(inc,1,MS);
free_dvector(lph,1,MS);
free_dvector(lan,1,MS);
free_dvector(arp,1,MS);
free_dvector(mna,1,MS);
free_dvector(hh,1,ls);
free_dvector(kk,1,ls);
free_dvector(pp,1,ls);
free_dvector(qq,1,ls);
free_dvector(cc,1,ls);
free_dvector(dd,1,ls);
free_dmatrix(nn,1,NEQ,1,ls);
free_dvector(nv,1,3);
free_dvector(y0,1,3);
free_dvector(u,1,3);
free_dvector(up,1,3);
free_dvector(epl,1,kmax);
free_dvector(phi,1,kmax);
free_dvector(tmv,1,kmax);
free_dmatrix(yy,1,NEQ,1,kmax);
free_dmatrix(nni,1,NEQ,1,kmax);
free_dvector(eei,1,kmax);
free_dvector(inci,1,kmax);
free_dvector(lphi,1,kmax);
free_dvector(lani,1,kmax);
free_dvector(lphu,1,ls);
free_dvector(lanu,1,ls);
free_dvector(cp,1,kmax);

printf("\n@ All done, exiting to system.\n");
printf("\n#===================================================# \n");

}
/*==================== driver() END ==========================*/
	
/*============================================================*/
/*==================== int main() ============================*/
/*============================================================*/
int main(int argc, char **argv)
{

driver(argc,argv);
	
return(0);	
}
/*==================== int main END ==========================*/
