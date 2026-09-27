#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "mhd.h"

/* State arrays may be stored in single precision (all arithmetic stays in
   double): the scheme is memory-bandwidth bound and float storage nearly
   halves the traffic.  Compile with -DUSE_FLOAT.                          */
#ifdef USE_FLOAT
typedef float real;
#else
typedef double real;
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ */
/* globals                                                             */
static Params P;
static real *U;    /* conserved, NV * ncell */
static real *U0;   /* stage copy */
static real *DU;   /* dU/dt accumulator */

#define IDX(i,j,k) ((((long)(i)*P.nt + (j))*P.nt) + (k))
#define UU(v,i,j,k) U[(long)(v)*P.ncell + IDX(i,j,k)]

/* ------------------------------------------------------------------ */
/* RNG: xoshiro256** + gaussian                                        */
static unsigned long rs[4];
static inline unsigned long rotl(const unsigned long x, int k){return (x<<k)|(x>>(64-k));}
static unsigned long next_rand(void){
    const unsigned long r = rotl(rs[1]*5,7)*9;
    const unsigned long t = rs[1]<<17;
    rs[2]^=rs[0]; rs[3]^=rs[1]; rs[1]^=rs[2]; rs[0]^=rs[3]; rs[2]^=t;
    rs[3]=rotl(rs[3],45);
    return r;
}
static void seed_rng(unsigned long s){
    for(int i=0;i<4;i++){ s+=0x9E3779B97f4A7C15UL; unsigned long z=s;
        z=(z^(z>>30))*0xBF58476D1CE4E5B9UL; z=(z^(z>>27))*0x94D049BB133111EBUL;
        rs[i]=z^(z>>31); }
    for(int i=0;i<64;i++) next_rand();
}
static double urand(void){ return (next_rand()>>11)*0x1.0p-53; }
static double grand(void){
    static int have=0; static double g2;
    if(have){have=0;return g2;}
    double u1=urand(), u2=urand();
    if(u1<1e-300) u1=1e-300;
    double r=sqrt(-2.0*log(u1)), th=2.0*M_PI*u2;
    g2=r*sin(th); have=1; return r*cos(th);
}

/* ------------------------------------------------------------------ */
/* isothermal HLLD Riemann solver (Mignone 2007) + Dedner GLM          */
/* w = [rho, un, ut1, ut2, Bn, Bt1, Bt2, psi]                          */
static inline void riemann(const double *wl, const double *wr, double *F)
{
    const double cs2 = P.cs*P.cs, ch = P.ch;
    /* --- GLM subsystem: exact solution of the (Bn,psi) linear system --- */
    double bn  = 0.5*(wl[4]+wr[4]) - 0.5*(wr[7]-wl[7])/ch;
    double psi = 0.5*(wl[7]+wr[7]) - 0.5*ch*(wr[4]-wl[4]);
    F[4] = psi;
    F[7] = ch*ch*bn;

    double rl=wl[0], ul=wl[1], vl=wl[2], wl2=wl[3], byl=wl[5], bzl=wl[6];
    double rr=wr[0], ur=wr[1], vr=wr[2], wr2=wr[3], byr=wr[5], bzr=wr[6];

    /* fast magnetosonic speeds */
    double b2l=(bn*bn+byl*byl+bzl*bzl)/rl, b2r=(bn*bn+byr*byr+bzr*bzr)/rr;
    double tl=cs2+b2l, tr=cs2+b2r;
    double dl=tl*tl-4.0*cs2*bn*bn/rl; if(dl<0) dl=0;
    double dr=tr*tr-4.0*cs2*bn*bn/rr; if(dr<0) dr=0;
    double cfl_=sqrt(0.5*(tl+sqrt(dl))), cfr_=sqrt(0.5*(tr+sqrt(dr)));
    double cmax = cfl_>cfr_?cfl_:cfr_;
    double SL = (ul<ur?ul:ur) - cmax;
    double SR = (ul>ur?ul:ur) + cmax;

    /* conserved and fluxes of the MHD subsystem (7 comps, Bn excluded) */
    double UL[7], UR[7], FL[7], FR[7];
    double ptl = cs2*rl + 0.5*(bn*bn+byl*byl+bzl*bzl);
    double ptr = cs2*rr + 0.5*(bn*bn+byr*byr+bzr*bzr);
    UL[0]=rl; UL[1]=rl*ul; UL[2]=rl*vl; UL[3]=rl*wl2; UL[5]=byl; UL[6]=bzl;
    UR[0]=rr; UR[1]=rr*ur; UR[2]=rr*vr; UR[3]=rr*wr2; UR[5]=byr; UR[6]=bzr;
    FL[0]=rl*ul; FL[1]=rl*ul*ul+ptl-bn*bn; FL[2]=rl*ul*vl-bn*byl; FL[3]=rl*ul*wl2-bn*bzl;
    FL[5]=ul*byl-vl*bn; FL[6]=ul*bzl-wl2*bn;
    FR[0]=rr*ur; FR[1]=rr*ur*ur+ptr-bn*bn; FR[2]=rr*ur*vr-bn*byr; FR[3]=rr*ur*wr2-bn*bzr;
    FR[5]=ur*byr-vr*bn; FR[6]=ur*bzr-wr2*bn;

    if(SL>=0.0){ F[0]=FL[0];F[1]=FL[1];F[2]=FL[2];F[3]=FL[3];F[5]=FL[5];F[6]=FL[6]; return; }
    if(SR<=0.0){ F[0]=FR[0];F[1]=FR[1];F[2]=FR[2];F[3]=FR[3];F[5]=FR[5];F[6]=FR[6]; return; }

    double idS=1.0/(SR-SL);
    double Uh[7], Fh[7];
    for(int q=0;q<7;q++){
        if(q==4) continue;
        Uh[q]=(SR*UR[q]-SL*UL[q]+FL[q]-FR[q])*idS;
        Fh[q]=(SR*FL[q]-SL*FR[q]+SR*SL*(UR[q]-UL[q]))*idS;
    }
    double rs_ = Uh[0];                 /* rho* (single value: no contact wave) */
    if(!(rs_>0.0)){ F[0]=Fh[0];F[1]=Fh[1];F[2]=Fh[2];F[3]=Fh[3];F[5]=Fh[5];F[6]=Fh[6]; return; }
    double SM = Uh[1]/rs_;

    double dnl = rl*(SL-ul)*(SL-SM) - bn*bn;
    double dnr = rr*(SR-ur)*(SR-SM) - bn*bn;
    double scl = fabs(rl*(SL-ul)*(SL-SM)) + bn*bn + 1e-300;
    double scr = fabs(rr*(SR-ur)*(SR-SM)) + bn*bn + 1e-300;
    if(fabs(dnl) < 1e-9*scl || fabs(dnr) < 1e-9*scr){
        F[0]=Fh[0];F[1]=Fh[1];F[2]=Fh[2];F[3]=Fh[3];F[5]=Fh[5];F[6]=Fh[6]; return;
    }
    double vsl = vl  - bn*byl*(SM-ul)/dnl;
    double wsl = wl2 - bn*bzl*(SM-ul)/dnl;
    double bysl= byl*(rl*(SL-ul)*(SL-ul)-bn*bn)/dnl;
    double bzsl= bzl*(rl*(SL-ul)*(SL-ul)-bn*bn)/dnl;
    double vsr = vr  - bn*byr*(SM-ur)/dnr;
    double wsr = wr2 - bn*bzr*(SM-ur)/dnr;
    double bysr= byr*(rr*(SR-ur)*(SR-ur)-bn*bn)/dnr;
    double bzsr= bzr*(rr*(SR-ur)*(SR-ur)-bn*bn)/dnr;

    double USL[7], USR[7];
    USL[0]=rs_; USL[1]=rs_*SM; USL[2]=rs_*vsl; USL[3]=rs_*wsl; USL[5]=bysl; USL[6]=bzsl;
    USR[0]=rs_; USR[1]=rs_*SM; USR[2]=rs_*vsr; USR[3]=rs_*wsr; USR[5]=bysr; USR[6]=bzsr;

    double ca = fabs(bn)/sqrt(rs_);
    double SLs = SM - ca, SRs = SM + ca;

    if(SLs>=0.0){
        for(int q=0;q<7;q++){ if(q==4)continue; F[q]=FL[q]+SL*(USL[q]-UL[q]); }
        return;
    }
    if(SRs<=0.0){
        for(int q=0;q<7;q++){ if(q==4)continue; F[q]=FR[q]+SR*(USR[q]-UR[q]); }
        return;
    }
    /* rotational (Alfven) middle state */
    double sg = (bn>0.0)?1.0:((bn<0.0)?-1.0:0.0);
    double sr_ = sqrt(rs_);
    double vss = 0.5*(vsl+vsr) + 0.5*sg*(bysr-bysl)/sr_;
    double wss = 0.5*(wsl+wsr) + 0.5*sg*(bzsr-bzsl)/sr_;
    double byss= 0.5*(bysl+bysr) + 0.5*sg*sr_*(vsr-vsl);
    double bzss= 0.5*(bzsl+bzsr) + 0.5*sg*sr_*(wsr-wsl);
    double USS[7];
    USS[0]=rs_; USS[1]=rs_*SM; USS[2]=rs_*vss; USS[3]=rs_*wss; USS[5]=byss; USS[6]=bzss;
    if(SM>=0.0){ for(int q=0;q<7;q++){ if(q==4)continue; F[q]=FL[q]+SL*(USL[q]-UL[q])+SLs*(USS[q]-USL[q]); } }
    else       { for(int q=0;q<7;q++){ if(q==4)continue; F[q]=FR[q]+SR*(USR[q]-UR[q])+SRs*(USS[q]-USR[q]); } }
}

/* ------------------------------------------------------------------ */
static inline double vanleer(double a, double b){
    double s=a*b; return (s>0.0)? 2.0*s/(a+b) : 0.0;
}

static void fill_ghosts(real *A)
{
    const int n=P.n, nt=P.nt;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int v=0;v<NV;v++)
    for(int j=0;j<nt;j++){
        real *a=A+(long)v*P.ncell;
        for(int k=0;k<nt;k++){
            for(int g=0;g<NG;g++){
                a[IDX(g,j,k)]      = a[IDX(g+n,j,k)];
                a[IDX(n+NG+g,j,k)] = a[IDX(NG+g,j,k)];
            }
        }
    }
    #pragma omp parallel for collapse(2) schedule(static)
    for(int v=0;v<NV;v++)
    for(int i=0;i<nt;i++){
        real *a=A+(long)v*P.ncell;
        for(int k=0;k<nt;k++){
            for(int g=0;g<NG;g++){
                a[IDX(i,g,k)]      = a[IDX(i,g+n,k)];
                a[IDX(i,n+NG+g,k)] = a[IDX(i,NG+g,k)];
            }
        }
    }
    #pragma omp parallel for collapse(2) schedule(static)
    for(int v=0;v<NV;v++)
    for(int i=0;i<nt;i++){
        real *a=A+(long)v*P.ncell;
        for(int j=0;j<nt;j++){
            for(int g=0;g<NG;g++){
                a[IDX(i,j,g)]      = a[IDX(i,j,g+n)];
                a[IDX(i,j,n+NG+g)] = a[IDX(i,j,NG+g)];
            }
        }
    }
}

/* sweep in direction d (0=x,1=y,2=z); accumulate -dF/dx into DU */
static void sweep_pencil(int d)
{
    const int n=P.n, nt=P.nt;
    const long ncell=P.ncell;
    const double idx=1.0/P.dx;
    /* component permutation: normal, t1, t2 */
    const int mn[3]={IMX,IMY,IMZ}, bn_[3]={IBX,IBY,IBZ};
    const int cn=d, c1=(d+1)%3, c2=(d+2)%3;
    long stride;
    if(d==0) stride=(long)nt*nt; else if(d==1) stride=nt; else stride=1;

    #pragma omp parallel
    {
        int np=n+2*NG;
        double *w  = malloc(sizeof(double)*8*np);
        double *wL = malloc(sizeof(double)*8*np);
        double *wR = malloc(sizeof(double)*8*np);
        double *Fc = malloc(sizeof(double)*8*(np+1));
        #pragma omp for collapse(2) schedule(static)
        for(int a=0;a<n;a++)
        for(int b=0;b<n;b++){
            int ia=a+NG, ib=b+NG;
            long base;
            if(d==0)      base=IDX(0,ia,ib);
            else if(d==1) base=IDX(ib,0,ia);
            else          base=IDX(ia,ib,0);
            /* gather primitives along the pencil */
            for(int i=0;i<np;i++){
                long o=base+(long)i*stride;
                double r=U[(long)IRHO*ncell+o];
                w[8*i+0]=r;
                w[8*i+1]=U[(long)mn[cn]*ncell+o]/r;
                w[8*i+2]=U[(long)mn[c1]*ncell+o]/r;
                w[8*i+3]=U[(long)mn[c2]*ncell+o]/r;
                w[8*i+4]=U[(long)bn_[cn]*ncell+o];
                w[8*i+5]=U[(long)bn_[c1]*ncell+o];
                w[8*i+6]=U[(long)bn_[c2]*ncell+o];
                w[8*i+7]=U[(long)IPSI*ncell+o];
            }
            /* PLM reconstruction (van Leer limiter on primitives) */
            for(int i=1;i<np-1;i++){
                for(int q=0;q<8;q++){
                    double dm=w[8*i+q]-w[8*(i-1)+q];
                    double dp=w[8*(i+1)+q]-w[8*i+q];
                    double s=0.5*vanleer(dm,dp);
                    wL[8*i+q]=w[8*i+q]-s;   /* state at left face of cell i  */
                    wR[8*i+q]=w[8*i+q]+s;   /* state at right face of cell i */
                }
                if(wL[8*i+0]<=P.rho_floor || wR[8*i+0]<=P.rho_floor){
                    for(int q=0;q<8;q++){ wL[8*i+q]=w[8*i+q]; wR[8*i+q]=w[8*i+q]; }
                }
            }
            /* fluxes at faces i+1/2 for i = NG-1 .. n+NG-1 */
            for(int i=NG-1;i<=n+NG-1;i++)
                riemann(&wR[8*i], &wL[8*(i+1)], &Fc[8*i]);
            /* update interior */
            for(int i=NG;i<n+NG;i++){
                long o=base+(long)i*stride;
                double *fm=&Fc[8*(i-1)], *fp=&Fc[8*i];
                DU[(long)IRHO*ncell+o]    -= idx*(fp[0]-fm[0]);
                DU[(long)mn[cn]*ncell+o]  -= idx*(fp[1]-fm[1]);
                DU[(long)mn[c1]*ncell+o]  -= idx*(fp[2]-fm[2]);
                DU[(long)mn[c2]*ncell+o]  -= idx*(fp[3]-fm[3]);
                DU[(long)bn_[cn]*ncell+o] -= idx*(fp[4]-fm[4]);
                DU[(long)bn_[c1]*ncell+o] -= idx*(fp[5]-fm[5]);
                DU[(long)bn_[c2]*ncell+o] -= idx*(fp[6]-fm[6]);
                DU[(long)IPSI*ncell+o]    -= idx*(fp[7]-fm[7]);
            }
        }
        free(w);free(wL);free(wR);free(Fc);
    }
}


/* Cache-blocked sweep for the x and y directions: VL contiguous z-columns are
   carried through reconstruction and the Riemann solve together, so every
   cache line that is touched is fully used (the naive strided pencil wastes
   7/8 of each line).                                                       */
#define VL 8
static void sweep_xy(int d)
{
    const int n=P.n, nt=P.nt; const long ncell=P.ncell; const double idx=1.0/P.dx;
    const int mn[3]={IMX,IMY,IMZ}, bn_[3]={IBX,IBY,IBZ};
    const int cn=d, c1=(d+1)%3, c2=(d+2)%3;
    const long sp=(d==0)?(long)nt*nt:(long)nt;
    const int np=n+2*NG;
    #pragma omp parallel
    {
        double *w =malloc(sizeof(double)*(size_t)np*8*VL);
        double *wl=malloc(sizeof(double)*(size_t)np*8*VL);
        double *wr=malloc(sizeof(double)*(size_t)np*8*VL);
        double *fc=malloc(sizeof(double)*(size_t)np*8*VL);
        #pragma omp for collapse(2) schedule(static)
        for(int oa=0;oa<n;oa++)
        for(int kb=0;kb<n;kb+=VL){
            int io=oa+NG;
            int vl=(n-kb<VL)?(n-kb):VL;
            long base=(d==0)?IDX(0,io,NG+kb):IDX(io,0,NG+kb);
            for(int i=0;i<np;i++){
                long off=base+(long)i*sp;
                double *W=w+(size_t)i*8*VL;
                const real *Ur=U+(long)IRHO*ncell+off;
                const real *Un=U+(long)mn[cn]*ncell+off;
                const real *U1=U+(long)mn[c1]*ncell+off;
                const real *U2=U+(long)mn[c2]*ncell+off;
                const real *Bn=U+(long)bn_[cn]*ncell+off;
                const real *B1=U+(long)bn_[c1]*ncell+off;
                const real *B2=U+(long)bn_[c2]*ncell+off;
                const real *Ps=U+(long)IPSI*ncell+off;
                for(int l=0;l<vl;l++){
                    double r=Ur[l], ir=1.0/r;
                    W[0*VL+l]=r;
                    W[1*VL+l]=Un[l]*ir;
                    W[2*VL+l]=U1[l]*ir;
                    W[3*VL+l]=U2[l]*ir;
                    W[4*VL+l]=Bn[l];
                    W[5*VL+l]=B1[l];
                    W[6*VL+l]=B2[l];
                    W[7*VL+l]=Ps[l];
                }
            }
            for(int i=1;i<np-1;i++){
                double *W=w+(size_t)i*8*VL, *Wm=w+(size_t)(i-1)*8*VL, *Wp=w+(size_t)(i+1)*8*VL;
                double *L=wl+(size_t)i*8*VL, *R=wr+(size_t)i*8*VL;
                for(int q=0;q<8;q++)
                for(int l=0;l<vl;l++){
                    double dm=W[q*VL+l]-Wm[q*VL+l];
                    double dp=Wp[q*VL+l]-W[q*VL+l];
                    double sgn=dm*dp;
                    double sl=(sgn>0.0)? sgn/(dm+dp) : 0.0;   /* 0.5*vanleer */
                    L[q*VL+l]=W[q*VL+l]-sl;
                    R[q*VL+l]=W[q*VL+l]+sl;
                }
                for(int l=0;l<vl;l++)
                    if(L[0*VL+l]<=P.rho_floor||R[0*VL+l]<=P.rho_floor)
                        for(int q=0;q<8;q++){L[q*VL+l]=W[q*VL+l];R[q*VL+l]=W[q*VL+l];}
            }
            for(int i=NG-1;i<=n+NG-1;i++){
                double *R=wr+(size_t)i*8*VL, *L=wl+(size_t)(i+1)*8*VL, *F=fc+(size_t)i*8*VL;
                double a[8],b[8],f[8];
                for(int l=0;l<vl;l++){
                    for(int q=0;q<8;q++){a[q]=R[q*VL+l];b[q]=L[q*VL+l];}
                    riemann(a,b,f);
                    for(int q=0;q<8;q++) F[q*VL+l]=f[q];
                }
            }
            for(int i=NG;i<n+NG;i++){
                long off=base+(long)i*sp;
                double *Fp=fc+(size_t)i*8*VL, *Fm=fc+(size_t)(i-1)*8*VL;
                real *Dr=DU+(long)IRHO*ncell+off;
                real *Dn=DU+(long)mn[cn]*ncell+off;
                real *D1=DU+(long)mn[c1]*ncell+off;
                real *D2=DU+(long)mn[c2]*ncell+off;
                real *Gn=DU+(long)bn_[cn]*ncell+off;
                real *G1=DU+(long)bn_[c1]*ncell+off;
                real *G2=DU+(long)bn_[c2]*ncell+off;
                real *Gp=DU+(long)IPSI*ncell+off;
                for(int l=0;l<vl;l++){
                    Dr[l]-=idx*(Fp[0*VL+l]-Fm[0*VL+l]);
                    Dn[l]-=idx*(Fp[1*VL+l]-Fm[1*VL+l]);
                    D1[l]-=idx*(Fp[2*VL+l]-Fm[2*VL+l]);
                    D2[l]-=idx*(Fp[3*VL+l]-Fm[3*VL+l]);
                    Gn[l]-=idx*(Fp[4*VL+l]-Fm[4*VL+l]);
                    G1[l]-=idx*(Fp[5*VL+l]-Fm[5*VL+l]);
                    G2[l]-=idx*(Fp[6*VL+l]-Fm[6*VL+l]);
                    Gp[l]-=idx*(Fp[7*VL+l]-Fm[7*VL+l]);
                }
            }
        }
        free(w);free(wl);free(wr);free(fc);
    }
}
static inline void sweep(int d){ if(d==2) sweep_pencil(2); else sweep_xy(d); }

static double max_signal(void)
{
    const int n=P.n; double vmax=0.0;
    #pragma omp parallel for collapse(2) reduction(max:vmax) schedule(static)
    for(int i=NG;i<n+NG;i++)
    for(int j=NG;j<n+NG;j++)
    for(int k=NG;k<n+NG;k++){
        long o=IDX(i,j,k);
        double r=U[(long)IRHO*P.ncell+o];
        double vx=U[(long)IMX*P.ncell+o]/r, vy=U[(long)IMY*P.ncell+o]/r, vz=U[(long)IMZ*P.ncell+o]/r;
        double bx=U[(long)IBX*P.ncell+o], by=U[(long)IBY*P.ncell+o], bz=U[(long)IBZ*P.ncell+o];
        double cf=sqrt(P.cs*P.cs+(bx*bx+by*by+bz*bz)/r);
        double s=fabs(vx); if(fabs(vy)>s)s=fabs(vy); if(fabs(vz)>s)s=fabs(vz);
        s+=cf; if(s>vmax) vmax=s;
    }
    return vmax;
}

/* ------------------------------------------------------------------ */
/* Ornstein-Uhlenbeck stochastic driving (Federrath et al. 2010)       */
#define MAXMODE 64
static int nmode=0;
static double kvec[MAXMODE][3];
static double amp_r[MAXMODE][3], amp_i[MAXMODE][3];
static double *Fx,*Fy,*Fz;              /* acceleration pattern */
static double g_vrms=0.0;               /* rms |v| measured during forcing */
static double *ctab,*stab;              /* per-mode cos/sin tables: [m][axis][i] */

static void driving_init(void)
{
    nmode=0;
    for(int a=-2;a<=2;a++)for(int b=-2;b<=2;b++)for(int c=-2;c<=2;c++){
        double kk=sqrt((double)(a*a+b*b+c*c));
        if(kk<1.0-1e-9||kk>2.0+1e-9) continue;
        /* keep only half-space (conjugate partner implied) */
        if(a<0) continue;
        if(a==0 && b<0) continue;
        if(a==0 && b==0 && c<0) continue;
        kvec[nmode][0]=a; kvec[nmode][1]=b; kvec[nmode][2]=c;
        nmode++;
    }
    ctab=malloc(sizeof(double)*nmode*3*P.n);
    stab=malloc(sizeof(double)*nmode*3*P.n);
    for(int m=0;m<nmode;m++)for(int ax=0;ax<3;ax++)for(int i=0;i<P.n;i++){
        double x=(i+0.5)*P.dx;
        double ph=2.0*M_PI*kvec[m][ax]*x;
        ctab[(m*3+ax)*P.n+i]=cos(ph);
        stab[(m*3+ax)*P.n+i]=sin(ph);
    }
    Fx=calloc(P.n*(size_t)P.n*P.n,sizeof(double));
    Fy=calloc(P.n*(size_t)P.n*P.n,sizeof(double));
    Fz=calloc(P.n*(size_t)P.n*P.n,sizeof(double));
    /* initialise amplitudes from the stationary distribution */
    for(int m=0;m<nmode;m++)for(int q=0;q<3;q++){amp_r[m][q]=0;amp_i[m][q]=0;}
    for(int m=0;m<nmode;m++){
        double gr[3],gi[3];
        for(int q=0;q<3;q++){gr[q]=grand();gi[q]=grand();}
        double k2=kvec[m][0]*kvec[m][0]+kvec[m][1]*kvec[m][1]+kvec[m][2]*kvec[m][2];
        double kdr=0,kdi=0;
        for(int q=0;q<3;q++){kdr+=kvec[m][q]*gr[q];kdi+=kvec[m][q]*gi[q];}
        for(int q=0;q<3;q++){
            amp_r[m][q]=P.zeta*gr[q]+(1.0-2.0*P.zeta)*kvec[m][q]*kdr/k2;
            amp_i[m][q]=P.zeta*gi[q]+(1.0-2.0*P.zeta)*kvec[m][q]*kdi/k2;
        }
    }
}

static void driving_ou(double dt)
{
    double f=exp(-dt/P.tdrive), g=sqrt(1.0-f*f);
    for(int m=0;m<nmode;m++){
        double gr[3],gi[3];
        for(int q=0;q<3;q++){gr[q]=grand();gi[q]=grand();}
        double k2=kvec[m][0]*kvec[m][0]+kvec[m][1]*kvec[m][1]+kvec[m][2]*kvec[m][2];
        double kdr=0,kdi=0;
        for(int q=0;q<3;q++){kdr+=kvec[m][q]*gr[q];kdi+=kvec[m][q]*gi[q];}
        for(int q=0;q<3;q++){
            double pr=P.zeta*gr[q]+(1.0-2.0*P.zeta)*kvec[m][q]*kdr/k2;
            double pi=P.zeta*gi[q]+(1.0-2.0*P.zeta)*kvec[m][q]*kdi/k2;
            amp_r[m][q]=f*amp_r[m][q]+g*pr;
            amp_i[m][q]=f*amp_i[m][q]+g*pi;
        }
    }
}

/* evaluate the acceleration pattern on the grid */
static void driving_pattern(void)
{
    const int n=P.n;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++){
        double cxy[MAXMODE],sxy[MAXMODE];
        for(int m=0;m<nmode;m++){
            double cx=ctab[(m*3+0)*n+i], sx=stab[(m*3+0)*n+i];
            double cy=ctab[(m*3+1)*n+j], sy=stab[(m*3+1)*n+j];
            cxy[m]=cx*cy-sx*sy; sxy[m]=sx*cy+cx*sy;
        }
        for(int k=0;k<n;k++){
            double ax=0,ay=0,az=0;
            for(int m=0;m<nmode;m++){
                double cz=ctab[(m*3+2)*n+k], sz=stab[(m*3+2)*n+k];
                double c=cxy[m]*cz-sxy[m]*sz;
                double s=sxy[m]*cz+cxy[m]*sz;
                ax+=amp_r[m][0]*c-amp_i[m][0]*s;
                ay+=amp_r[m][1]*c-amp_i[m][1]*s;
                az+=amp_r[m][2]*c-amp_i[m][2]*s;
            }
            size_t o=((size_t)i*n+j)*n+k;
            Fx[o]=2.0*ax; Fy[o]=2.0*ay; Fz[o]=2.0*az;
        }
    }
}

/* apply forcing with the amplitude that injects exactly edot*dt of kinetic
   energy (Mac Low 1999; Mac Low & Klessen 2004)                       */
static void driving_apply(double dt)
{
    const int n=P.n;
    const double dV=P.dx*P.dx*P.dx;
    double s2=0.0, sv=0.0, sf2=0.0, v2=0.0;
    #pragma omp parallel for collapse(2) reduction(+:s2,sv,sf2,v2) schedule(static)
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
    for(int k=0;k<n;k++){
        size_t o=((size_t)i*n+j)*n+k;
        long q=IDX(i+NG,j+NG,k+NG);
        double r=U[(long)IRHO*P.ncell+q];
        double f2=Fx[o]*Fx[o]+Fy[o]*Fy[o]+Fz[o]*Fz[o];
        sf2+=f2;
        s2+=r*f2;
        sv+=U[(long)IMX*P.ncell+q]*Fx[o]+U[(long)IMY*P.ncell+q]*Fy[o]+U[(long)IMZ*P.ncell+q]*Fz[o];
        double mx=U[(long)IMX*P.ncell+q],my=U[(long)IMY*P.ncell+q],mz=U[(long)IMZ*P.ncell+q];
        v2+=(mx*mx+my*my+mz*mz)/(r*r);
    }
    g_vrms=sqrt(v2/((double)n*n*n));
    double frms=sqrt(sf2/((double)n*n*n));
    if(frms<=0) return;
    /* normalise pattern to unit rms so 'A' has units of velocity */
    double a2=0.5*s2*dV/(frms*frms), b1=sv*dV/frms, c1=P.edot*dt;
    double A=(-b1+sqrt(b1*b1+4.0*a2*c1))/(2.0*a2);
    double scal=A/frms;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
    for(int k=0;k<n;k++){
        size_t o=((size_t)i*n+j)*n+k;
        long q=IDX(i+NG,j+NG,k+NG);
        double r=U[(long)IRHO*P.ncell+q];
        U[(long)IMX*P.ncell+q]+=r*scal*Fx[o];
        U[(long)IMY*P.ncell+q]+=r*scal*Fy[o];
        U[(long)IMZ*P.ncell+q]+=r*scal*Fz[o];
    }
    /* remove any net momentum */
    double px=0,py=0,pz=0,mt=0;
    #pragma omp parallel for collapse(2) reduction(+:px,py,pz,mt) schedule(static)
    for(int i=NG;i<n+NG;i++)
    for(int j=NG;j<n+NG;j++)
    for(int k=NG;k<n+NG;k++){
        long q=IDX(i,j,k);
        px+=U[(long)IMX*P.ncell+q]; py+=U[(long)IMY*P.ncell+q]; pz+=U[(long)IMZ*P.ncell+q];
        mt+=U[(long)IRHO*P.ncell+q];
    }
    double ux=px/mt, uy=py/mt, uz=pz/mt;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i=NG;i<n+NG;i++)
    for(int j=NG;j<n+NG;j++)
    for(int k=NG;k<n+NG;k++){
        long q=IDX(i,j,k);
        double r=U[(long)IRHO*P.ncell+q];
        U[(long)IMX*P.ncell+q]-=r*ux; U[(long)IMY*P.ncell+q]-=r*uy; U[(long)IMZ*P.ncell+q]-=r*uz;
    }
}

/* ------------------------------------------------------------------ */
/* Protect the state after an update.  Clamping the density while leaving the
   momentum untouched is fatal: v = m/rho then explodes by the same factor the
   density was raised by, the time step collapses and the run is destroyed.
   We clamp the *velocity* instead - i.e. rescale the momentum with the
   density - and zero the momentum outright in cells whose density came out
   non-positive, where the velocity carries no meaning.  A velocity ceiling
   well above anything physical catches the remaining pathologies.        */
static long g_nfix=0;
static void apply_floor(void)
{
    const int n=P.n;
    const double vcap=40.0*(P.ms_target>1.0?P.ms_target:1.0)*P.cs;
    const double vcap2=vcap*vcap;
    long nfix=0;
    #pragma omp parallel for collapse(2) reduction(+:nfix) schedule(static)
    for(int i=NG;i<n+NG;i++)
    for(int j=NG;j<n+NG;j++)
    for(int k=NG;k<n+NG;k++){
        long o=IDX(i,j,k);
        real *r=&U[(long)IRHO*P.ncell+o];
        real *mx=&U[(long)IMX*P.ncell+o], *my=&U[(long)IMY*P.ncell+o],
             *mz=&U[(long)IMZ*P.ncell+o];
        if(!(*r>P.rho_floor)){
            if(*r>0.0){ double f=P.rho_floor/(*r); *mx*=f; *my*=f; *mz*=f; }
            else      { *mx=0.0; *my=0.0; *mz=0.0; }
            *r=P.rho_floor; nfix++;
        }
        double ir=1.0/(*r);
        double v2=((double)(*mx)*(*mx)+(double)(*my)*(*my)+(double)(*mz)*(*mz))*ir*ir;
        if(v2>vcap2){
            double f=vcap/sqrt(v2); *mx*=f; *my*=f; *mz*=f; nfix++;
        }
    }
    g_nfix+=nfix;
}

static void rk2_step(double dt)
{
    const long tot=(long)NV*P.ncell;
    memcpy(U0,U,sizeof(real)*tot);
    /* stage 1 */
    memset(DU,0,sizeof(real)*tot);
    fill_ghosts(U);
    sweep(0); sweep(1); sweep(2);
    #pragma omp parallel for schedule(static)
    for(long q=0;q<tot;q++) U[q]+=dt*DU[q];
    apply_floor();   /* stage 1 */
    /* stage 2 */
    memset(DU,0,sizeof(real)*tot);
    fill_ghosts(U);
    sweep(0); sweep(1); sweep(2);
    #pragma omp parallel for schedule(static)
    for(long q=0;q<tot;q++) U[q]=0.5*U0[q]+0.5*(U[q]+dt*DU[q]);
    apply_floor();
    /* Dedner parabolic damping of psi */
    double dec=exp(-P.cfl/P.glm_alpha);
    #pragma omp parallel for schedule(static)
    for(long q=0;q<P.ncell;q++) U[(long)IPSI*P.ncell+q]*=dec;
}

typedef struct { double vrms,vrms_mw,ms,ma,rmin,rmax,ekin,emag,divb,sigma_lnr; } Diag;

static Diag diagnose(void)
{
    const int n=P.n; Diag D; double v2=0,v2m=0,mt=0,rmin=1e300,rmax=-1e300,ek=0,em=0,db=0,bmean=0;
    double sl=0,sl2=0;
    #pragma omp parallel for collapse(2) reduction(+:v2,v2m,mt,ek,em,db,bmean,sl,sl2) \
        reduction(min:rmin) reduction(max:rmax) schedule(static)
    for(int i=NG;i<n+NG;i++)
    for(int j=NG;j<n+NG;j++)
    for(int k=NG;k<n+NG;k++){
        long o=IDX(i,j,k);
        double r=U[(long)IRHO*P.ncell+o];
        double vx=U[(long)IMX*P.ncell+o]/r, vy=U[(long)IMY*P.ncell+o]/r, vz=U[(long)IMZ*P.ncell+o]/r;
        double bx=U[(long)IBX*P.ncell+o], by=U[(long)IBY*P.ncell+o], bz=U[(long)IBZ*P.ncell+o];
        double vv=vx*vx+vy*vy+vz*vz;
        v2+=vv; v2m+=r*vv; mt+=r; ek+=0.5*r*vv; em+=0.5*(bx*bx+by*by+bz*bz);
        if(r<rmin)rmin=r; if(r>rmax)rmax=r;
        double lr=log(r); sl+=lr; sl2+=lr*lr;
        bmean+=sqrt(bx*bx+by*by+bz*bz);
        double d=(U[(long)IBX*P.ncell+IDX(i+1,j,k)]-U[(long)IBX*P.ncell+IDX(i-1,j,k)]
                 +U[(long)IBY*P.ncell+IDX(i,j+1,k)]-U[(long)IBY*P.ncell+IDX(i,j-1,k)]
                 +U[(long)IBZ*P.ncell+IDX(i,j,k+1)]-U[(long)IBZ*P.ncell+IDX(i,j,k-1)])/(2.0*P.dx);
        db+=d*d;
    }
    double nc=(double)n*n*n;
    D.vrms=sqrt(v2/nc); D.vrms_mw=sqrt(v2m/mt); D.rmin=rmin; D.rmax=rmax;
    D.ekin=ek/nc; D.emag=em/nc;
    D.ms=D.vrms/P.cs;
    double b0=P.ms_target/P.ma_target;   /* mean field, rho_mean = 1 */
    D.ma=(b0>0)?D.vrms/b0:1e30;
    D.divb=sqrt(db/nc)*P.dx/(bmean/nc+1e-30);
    double m=sl/nc; D.sigma_lnr=sqrt(sl2/nc-m*m);
    return D;
}

static void write_snapshot(const char *dir,int num,double t)
{
    char fn[700]; snprintf(fn,sizeof(fn),"%s/snap_%04d.bin",dir,num);
    FILE *f=fopen(fn,"wb"); if(!f){fprintf(stderr,"cannot open %s\n",fn);return;}
    int n=P.n; float *buf=malloc(sizeof(float)*(size_t)n*n*n);
    int vmap[7]={IRHO,IMX,IMY,IMZ,IBX,IBY,IBZ};
    for(int v=0;v<7;v++){
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++){
            long o=IDX(i+NG,j+NG,k+NG);
            double val=U[(long)vmap[v]*P.ncell+o];
            if(v>=1&&v<=3) val/=U[(long)IRHO*P.ncell+o];   /* store velocity */
            buf[((size_t)i*n+j)*n+k]=(float)val;
        }
        fwrite(buf,sizeof(float),(size_t)n*n*n,f);
    }
    free(buf); fclose(f);
    snprintf(fn,sizeof(fn),"%s/snap_%04d.txt",dir,num);
    f=fopen(fn,"w");
    Diag D=diagnose();
    fprintf(f,"n %d\ntime %.8e\nfields rho vx vy vz bx by bz\ndtype float32\n"
              "ms %.6f\nma %.6f\nvrms %.6f\nb0 %.6f\ncs %.6f\n",
              n,t,D.ms,D.ma,D.vrms,P.ms_target/P.ma_target,P.cs);
    fclose(f);
}


/* ------------------------------------------------------------------ */
/* checkpoint / restart (needed for long cluster runs)                  */
static int write_checkpoint(const char *dir,double t,int step,int isnap)
{
    char fn[700],tmp[700];
    snprintf(fn,sizeof(fn),"%s/checkpoint.bin",dir);
    snprintf(tmp,sizeof(tmp),"%s/checkpoint.tmp",dir);
    FILE *f=fopen(tmp,"wb"); if(!f) return 0;
    int n=P.n,rsz=(int)sizeof(real);
    fwrite(&n,sizeof(int),1,f); fwrite(&rsz,sizeof(int),1,f);
    fwrite(&t,sizeof(double),1,f); fwrite(&step,sizeof(int),1,f);
    fwrite(&isnap,sizeof(int),1,f);
    fwrite(U,sizeof(real),(size_t)NV*P.ncell,f);
    fwrite(rs,sizeof(unsigned long),4,f);
    fwrite(&nmode,sizeof(int),1,f);
    fwrite(amp_r,sizeof(double),MAXMODE*3,f);
    fwrite(amp_i,sizeof(double),MAXMODE*3,f);
    fwrite(&P.edot,sizeof(double),1,f);
    fclose(f); rename(tmp,fn); return 1;
}
static int read_checkpoint(const char *dir,double *t,int *step,int *isnap)
{
    char fn[700]; snprintf(fn,sizeof(fn),"%s/checkpoint.bin",dir);
    FILE *f=fopen(fn,"rb"); if(!f) return 0;
    int n,rsz,nm;
    if(fread(&n,sizeof(int),1,f)!=1||fread(&rsz,sizeof(int),1,f)!=1){fclose(f);return 0;}
    if(n!=P.n||rsz!=(int)sizeof(real)){fclose(f);fprintf(stderr,"checkpoint mismatch\n");return 0;}
    size_t ok=1;
    ok&=fread(t,sizeof(double),1,f)==1;
    ok&=fread(step,sizeof(int),1,f)==1;
    ok&=fread(isnap,sizeof(int),1,f)==1;
    ok&=fread(U,sizeof(real),(size_t)NV*P.ncell,f)==(size_t)NV*P.ncell;
    ok&=fread(rs,sizeof(unsigned long),4,f)==4;
    ok&=fread(&nm,sizeof(int),1,f)==1;
    ok&=fread(amp_r,sizeof(double),MAXMODE*3,f)==MAXMODE*3;
    ok&=fread(amp_i,sizeof(double),MAXMODE*3,f)==MAXMODE*3;
    ok&=fread(&P.edot,sizeof(double),1,f)==1;
    fclose(f);
    if(!ok){fprintf(stderr,"checkpoint truncated\n");return 0;}
    printf("restarted from checkpoint: t=%.6f step=%d snap=%d\n",*t,*step,*isnap);
    return 1;
}

/* ------------------------------------------------------------------ */
static void alloc_state(void)
{
    P.nt=P.n+2*NG; P.ncell=(long)P.nt*P.nt*P.nt; P.dx=P.L/P.n;
    U =malloc(sizeof(real)*NV*P.ncell);
    U0=malloc(sizeof(real)*NV*P.ncell);
    DU=malloc(sizeof(real)*NV*P.ncell);
    if(!U||!U0||!DU){fprintf(stderr,"allocation failed\n");exit(1);}
    memset(U,0,sizeof(real)*NV*P.ncell);
}
static inline void setcell(int i,int j,int k,double r,double vx,double vy,double vz,
                           double bx,double by,double bz)
{
    long o=IDX(i+NG,j+NG,k+NG);
    U[(long)IRHO*P.ncell+o]=r;
    U[(long)IMX*P.ncell+o]=r*vx; U[(long)IMY*P.ncell+o]=r*vy; U[(long)IMZ*P.ncell+o]=r*vz;
    U[(long)IBX*P.ncell+o]=bx; U[(long)IBY*P.ncell+o]=by; U[(long)IBZ*P.ncell+o]=bz;
    U[(long)IPSI*P.ncell+o]=0.0;
}


/* ------------------------------------------------------------------ */
/* initialise from an upsampled lower-resolution snapshot: the large
   scales are already in statistical equilibrium, so only a short
   adjustment is needed for the new small scales to fill in.           */
static void init_upsample(const char *fn)
{
    FILE *f=fopen(fn,"rb");
    if(!f){fprintf(stderr,"cannot open %s\n",fn);exit(1);}
    fseek(f,0,SEEK_END); long sz=ftell(f); fseek(f,0,SEEK_SET);
    long ncell_src=sz/(7*4);
    int ns=(int)(round(cbrt((double)ncell_src)));
    if((long)ns*ns*ns!=ncell_src){fprintf(stderr,"bad snapshot size %ld\n",sz);exit(1);}
    printf("upsampling %s: %d^3 -> %d^3\n",fn,ns,P.n);
    float *src=malloc(sizeof(float)*7*ncell_src);
    if(fread(src,sizeof(float),7*ncell_src,f)!=(size_t)(7*ncell_src)){
        fprintf(stderr,"short read\n");exit(1);}
    fclose(f);
    #define SRC(v,i,j,k) src[(long)(v)*ncell_src + (((long)(i)*ns+(j))*ns+(k))]
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i=0;i<P.n;i++)
    for(int j=0;j<P.n;j++)
    for(int k=0;k<P.n;k++){
        double xs=((i+0.5)/(double)P.n)*ns-0.5;
        double ys=((j+0.5)/(double)P.n)*ns-0.5;
        double zs=((k+0.5)/(double)P.n)*ns-0.5;
        int i0=(int)floor(xs), j0=(int)floor(ys), k0=(int)floor(zs);
        double fx=xs-i0, fy=ys-j0, fz=zs-k0;
        int ip=((i0%ns)+ns)%ns, jp=((j0%ns)+ns)%ns, kp=((k0%ns)+ns)%ns;
        int in=(ip+1)%ns, jn=(jp+1)%ns, kn=(kp+1)%ns;
        double q[7];
        for(int v=0;v<7;v++){
            double c00=SRC(v,ip,jp,kp)*(1-fx)+SRC(v,in,jp,kp)*fx;
            double c01=SRC(v,ip,jp,kn)*(1-fx)+SRC(v,in,jp,kn)*fx;
            double c10=SRC(v,ip,jn,kp)*(1-fx)+SRC(v,in,jn,kp)*fx;
            double c11=SRC(v,ip,jn,kn)*(1-fx)+SRC(v,in,jn,kn)*fx;
            double c0=c00*(1-fy)+c10*fy, c1=c01*(1-fy)+c11*fy;
            q[v]=c0*(1-fz)+c1*fz;
        }
        if(q[0]<P.rho_floor) q[0]=P.rho_floor;
        setcell(i,j,k,q[0],q[1],q[2],q[3],q[4],q[5],q[6]);
    }
    #undef SRC
    free(src);
}

/* --- circularly polarised Alfven wave (exact nonlinear solution) --- */
static double run_cpaw(int n,double tend)
{
    P.n=n; P.L=1.0; P.cs=1.0; P.cfl=0.4; P.glm_alpha=0.18; P.rho_floor=1e-8;
    alloc_state();
    double A=0.1,B0=1.0;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++){
        double x=(i+0.5)*P.dx, ph=2.0*M_PI*x;
        double by=A*sin(ph), bz=A*cos(ph);
        setcell(i,j,k,1.0,0.0,-by,-bz,B0,by,bz);   /* travelling in +x */
    }
    real *init=malloc(sizeof(real)*NV*P.ncell);
    memcpy(init,U,sizeof(real)*NV*P.ncell);
    double t=0;
    while(t<tend-1e-12){
        fill_ghosts(U);
        double smax=max_signal();
        double dt=P.cfl*P.dx/smax; if(t+dt>tend)dt=tend-t;
        P.ch=P.cfl*P.dx/dt;
        rk2_step(dt); t+=dt;
    }
    double err=0; long cnt=0;
    for(int i=NG;i<n+NG;i++)for(int j=NG;j<n+NG;j++)for(int k=NG;k<n+NG;k++){
        long o=IDX(i,j,k);
        err+=fabs(U[(long)IBY*P.ncell+o]-init[(long)IBY*P.ncell+o])
            +fabs(U[(long)IBZ*P.ncell+o]-init[(long)IBZ*P.ncell+o]);
        cnt+=2;
    }
    free(init);free(U);free(U0);free(DU);
    return err/cnt;
}

/* --- isothermal shock tube (hydro or MHD) --- */
static void run_tube(int n,double tend,int mhd,const char *out)
{
    P.n=n; P.L=1.0; P.cs=1.0; P.cfl=0.4; P.glm_alpha=0.18; P.rho_floor=1e-8;
    alloc_state();
    int ny=4;   /* thin in y,z but code is cubic: use n^3 with uniform y,z */
    (void)ny;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++){
        double x=(i+0.5)*P.dx;
        int left = (x<0.5);
        double r  = left?1.0:0.125;
        double vx = 0.0;
        double bx = mhd?0.75:0.0;
        double by = mhd?(left?1.0:-1.0):0.0;
        setcell(i,j,k,r,vx,0.0,0.0,bx,by,0.0);
    }
    double t=0;
    while(t<tend-1e-12){
        fill_ghosts(U);
        double smax=max_signal();
        double dt=P.cfl*P.dx/smax; if(t+dt>tend)dt=tend-t;
        P.ch=P.cfl*P.dx/dt;
        rk2_step(dt); t+=dt;
    }
    FILE *f=fopen(out,"w");
    fprintf(f,"# x rho vx vy by\n");
    for(int i=0;i<n;i++){
        long o=IDX(i+NG,NG,NG);
        double r=U[(long)IRHO*P.ncell+o];
        fprintf(f,"%.8e %.8e %.8e %.8e %.8e\n",(i+0.5)*P.dx,r,
                U[(long)IMX*P.ncell+o]/r,U[(long)IMY*P.ncell+o]/r,U[(long)IBY*P.ncell+o]);
    }
    fclose(f);
    free(U);free(U0);free(DU);
}

/* --- divergence cleaning test: random B with large div B --- */
static void run_divtest(int n)
{
    P.n=n;P.L=1.0;P.cs=1.0;P.cfl=0.4;P.glm_alpha=0.18;P.rho_floor=1e-8;
    P.ms_target=1;P.ma_target=1;
    alloc_state(); seed_rng(12345);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++){
        double x=(i+0.5)*P.dx,y=(j+0.5)*P.dx,z=(k+0.5)*P.dx;
        setcell(i,j,k,1.0,0.0,0.0,0.0,
                1.0+0.3*sin(2*M_PI*x),0.3*cos(2*M_PI*y),0.3*sin(2*M_PI*z));
    }
    fill_ghosts(U);
    printf("# divB test: step  time  |divB|dx/|B|\n");
    double t=0;
    for(int s=0;s<200;s++){
        fill_ghosts(U);
        Diag D=diagnose();
        if(s%20==0) printf("%5d %10.4e %12.5e\n",s,t,D.divb);
        double smax=max_signal(); double dt=P.cfl*P.dx/smax;
        P.ch=P.cfl*P.dx/dt; rk2_step(dt); t+=dt;
    }
    fill_ghosts(U); Diag D=diagnose();
    printf("%5d %10.4e %12.5e\n",200,t,D.divb);
    free(U);free(U0);free(DU);
}

/* ------------------------------------------------------------------ */
static double getarg(int argc,char**argv,const char*key,double def)
{
    char buf[128]; snprintf(buf,sizeof(buf),"%s=",key); size_t L=strlen(buf);
    for(int i=1;i<argc;i++) if(!strncmp(argv[i],buf,L)) return atof(argv[i]+L);
    return def;
}
static const char* getargs(int argc,char**argv,const char*key,const char*def)
{
    char buf[128]; snprintf(buf,sizeof(buf),"%s=",key); size_t L=strlen(buf);
    for(int i=1;i<argc;i++) if(!strncmp(argv[i],buf,L)) return argv[i]+L;
    return def;
}

int main(int argc,char**argv)
{
    const char *test=getargs(argc,argv,"test","");
    if(!strcmp(test,"cpaw")){
        printf("# circularly polarised Alfven wave, one period\n# N  L1(B_perp)  order\n");
        double prev=0;
        for(int n=16;n<=64;n*=2){
            double e=run_cpaw(n,1.0);
            if(prev>0) printf("%4d %12.5e  %6.3f\n",n,e,log2(prev/e));
            else       printf("%4d %12.5e     -\n",n,e);
            prev=e;
        }
        return 0;
    }
    if(!strcmp(test,"tube")){
        int n=(int)getarg(argc,argv,"n",512);
        double tend=getarg(argc,argv,"tend",0.1);
        int mhd=(int)getarg(argc,argv,"mhd",0);
        run_tube(n,tend,mhd,getargs(argc,argv,"out","tube.txt"));
        return 0;
    }
    if(!strcmp(test,"divb")){ run_divtest((int)getarg(argc,argv,"n",64)); return 0; }

    /* ---------------- turbulence production run ---------------- */
    P.n=(int)getarg(argc,argv,"n",64);
    P.L=1.0; P.cs=1.0;
    P.cfl=getarg(argc,argv,"cfl",0.4);
    P.glm_alpha=getarg(argc,argv,"glm",0.18);
    P.rho_floor=getarg(argc,argv,"rhofloor",1e-5);
    P.ms_target=getarg(argc,argv,"ms",3.0);
    P.ma_target=getarg(argc,argv,"ma",1.0);
    P.zeta=getarg(argc,argv,"zeta",0.5);
    P.seed=(unsigned long)getarg(argc,argv,"seed",42);
    double edot_fac=getarg(argc,argv,"edotfac",1.0);
    double ntdyn=getarg(argc,argv,"ntdyn",10.0);
    double tsnap0=getarg(argc,argv,"tsnap0",5.0);
    int nsnap=(int)getarg(argc,argv,"nsnap",21);
    snprintf(P.outdir,sizeof(P.outdir),"%s",getargs(argc,argv,"out","run"));
    mkdir(P.outdir,0755);

    double tdyn=0.5*P.L/P.ms_target;          /* L0/v_rms, L0 = L/2 */
    P.tdrive=tdyn;
    /* Initial guess for the injection rate.  In equilibrium Edot balances
       the dissipation rate, which drops steeply as the field gets stronger
       (sub-Alfvenic turbulence cascades slowly across the mean field), hence
       the empirical M_A^2 factor.  The controller removes the residual. */
    /* Empirical calibration: the dissipation rate a driven box settles to
       depends on the field strength as well as on M_S (a strong mean field
       makes the cascade markedly less dissipative), so Edot ~ M_S^3 alone
       misses the target M_S by up to 60%.  The M_A^1.5 factor was measured
       from a first pass of the suite.  Achieved Mach numbers still differ
       somewhat from nominal, as in the paper, and every result is reported
       against the measured values.                                     */
    P.edot=edot_fac*pow(P.ms_target,3.0)*pow(P.ma_target,1.5)*P.L*P.L;
    double edot0=getarg(argc,argv,"edot0",0.0);
    if(edot0>0.0) P.edot=edot0;
    P.tmax=ntdyn*tdyn;

    alloc_state();
    seed_rng(P.seed);
    double b0=P.ms_target/P.ma_target;        /* v_A = B/sqrt(rho), rho=1 */
    const char *ups=getargs(argc,argv,"upsample","");
    for(int i=0;i<P.n;i++)for(int j=0;j<P.n;j++)for(int k=0;k<P.n;k++)
        setcell(i,j,k,1.0,0,0,0,b0,0,0);
    driving_init();

    double t=0,tstart=0; int step=0,isnap=0;
    int restart=(int)getarg(argc,argv,"restart",1);
    int have_ck=0;
    if(restart) have_ck=read_checkpoint(P.outdir,&t,&step,&isnap);
    if(!have_ck && ups[0]) init_upsample(ups);
    tstart=t;
    double chkint=getarg(argc,argv,"chkint",1000);
    int adapt=(int)getarg(argc,argv,"adapt",0);
    double tadapt=getarg(argc,argv,"tadapt",0.583);
    double tauctl=getarg(argc,argv,"tauctl",3.0);
    double eavg=0.0,ewt=0.0;

    char logn[700]; snprintf(logn,sizeof(logn),"%s/history.txt",P.outdir);
    FILE *lg=fopen(logn,have_ck?"a":"w");
    if(!have_ck) fprintf(lg,"# step time dt Ms Ma Ms_mw rho_min rho_max Ekin Emag divB sigma_lnrho edot nfix\n");
    char cfgn[700]; snprintf(cfgn,sizeof(cfgn),"%s/params.txt",P.outdir);
    FILE *cf=fopen(cfgn,"w");
    fprintf(cf,"n %d\nms_target %g\nma_target %g\nb0 %g\nedot %g\ntdyn %g\ntmax %g\n"
               "cfl %g\nzeta %g\nseed %lu\n",P.n,P.ms_target,P.ma_target,b0,P.edot,tdyn,P.tmax,
               P.cfl,P.zeta,P.seed);
    fclose(cf);

    double dsnap=(nsnap>1)?((ntdyn-tsnap0)*tdyn/(nsnap-1)):1e30;
    double tsnap=tsnap0*tdyn+isnap*dsnap;
    (void)tstart;
    while(t<P.tmax-1e-12){
        fill_ghosts(U);
        double smax=max_signal();
        double dt=P.cfl*P.dx/smax;
        if(t+dt>P.tmax) dt=P.tmax-t;
        if(isnap<nsnap && t+dt>tsnap) dt=tsnap-t;
        if(dt<=0) dt=1e-12;
        P.ch=P.cfl*P.dx/dt;
        rk2_step(dt);
        driving_ou(dt); driving_pattern(); driving_apply(dt);
        t+=dt; step++;
        /* During the first half of the run the energy injection rate is
           relaxed towards the value that yields the targeted M_S (the
           equilibrium v_rms depends on the field strength as well as on
           Edot).  It is frozen afterwards, so the analysis interval has a
           strictly constant injection rate, as in the paper.            */
        if(adapt && g_vrms>0.0){
            if(t<tadapt*P.tmax){
                double err=log(g_vrms/(P.ms_target*P.cs));
                P.edot*=exp(-3.0*err*dt/(tauctl*tdyn));
                if(t>0.6*tadapt*P.tmax){ eavg+=log(P.edot)*dt; ewt+=dt; }
            } else if(ewt>0.0){       /* freeze at the time-averaged value */
                P.edot=exp(eavg/ewt); ewt=-1.0;
                printf("  Edot frozen at %.4g (t=%.3f)\n",P.edot,t);
            }
        }
        if(step%10==0||step<5){
            Diag D=diagnose();
            /* Die loudly on a blow-up instead of quietly writing a suite of
               NaN snapshots.  (Build with -fno-finite-math-only: plain
               -ffast-math implies -ffinite-math-only, under which the
               compiler assumes NaN cannot occur and deletes checks like
               this one and the density-floor guard in apply_floor.)      */
            if(!isfinite(D.vrms)||!isfinite(D.rmin)||!(D.rmin>0.0)){
                fprintf(stderr,"FATAL: solution went non-finite at step %d, "
                        "t=%.6e (vrms=%g rho_min=%g nfix=%ld)\n",
                        step,t,D.vrms,D.rmin,g_nfix);
                fflush(stderr); fclose(lg);
                return 2;
            }
            fprintf(lg,"%d %.8e %.6e %.6f %.6f %.6f %.4e %.4e %.6e %.6e %.4e %.4f %.6e %ld\n",
                    step,t,dt,D.ms,D.ma,D.vrms_mw/P.cs,D.rmin,D.rmax,D.ekin,D.emag,D.divb,D.sigma_lnr,P.edot,g_nfix);
            fflush(lg);
            if(step%200==0)
                printf("step %6d t=%.4f (%.1f tdyn) dt=%.2e Ms=%.3f Ma=%.3f rho=[%.2e,%.2e] divB=%.2e\n",
                       step,t,t/tdyn,dt,D.ms,D.ma,D.rmin,D.rmax,D.divb);
        }
        if(step%(int)chkint==0) write_checkpoint(P.outdir,t,step,isnap);
        if(isnap<nsnap && t>=tsnap-1e-12){
            fill_ghosts(U);
            write_snapshot(P.outdir,isnap,t);
            write_checkpoint(P.outdir,t,isnap==nsnap-1?step:step,isnap+1);
            printf("  snapshot %d at t=%.4f (%.2f tdyn)\n",isnap,t,t/tdyn);
            fflush(stdout);
            isnap++; tsnap=tsnap0*tdyn+isnap*dsnap;
        }
    }
    fclose(lg);
    { char pn[700]; snprintf(pn,sizeof(pn),"%s/params.txt",P.outdir);
      FILE *pf=fopen(pn,"a"); if(pf){fprintf(pf,"edot_final %.8g\nsteps %d\n",P.edot,step);fclose(pf);} }
    printf("done: %d steps (Edot=%.4g)\n",step,P.edot);
    return 0;
}
