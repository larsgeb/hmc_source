
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aux.hpp"

const double PI = 3.14159265358979323846264338327;

/*=================================================================================*/
/* Modelling parameters class. ----------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/
parameters::parameters()
{

}

/* Copy constructor. --------------------------------------------------------------*/

parameters::parameters(const parameters &a)
{
    vp=a.vp;
    vs=a.vs;
    rho=a.rho;
    
    for (int i=0; i<20; i++)
    {
        q[i]=a.q[i];
        sigma_q[i]=a.sigma_q[i];
        mean_q[i]=a.mean_q[i];
    }
}

/* Assignment operator. -----------------------------------------------------------*/

parameters & parameters::operator=(const parameters &a)
{
    /* Check for self-assignment. */
    if (this == &a) return *this;
    
    vp=a.vp;
    vs=a.vs;
    rho=a.rho;
    
    for (int i=0; i<20; i++)
    {
        q[i]=a.q[i];
        sigma_q[i]=a.sigma_q[i];
        mean_q[i]=a.mean_q[i];
    }
    
    return *this;
}

/* Destructor. --------------------------------------------------------------------*/
parameters::~parameters()
{
    
}

/* Read input from file. ----------------------------------------------------------*/
void parameters::read_input(const char *filename)
{
    FILE *fid;
    char str[1000];
    
    fid=fopen(filename,"r");
    
    /* Medium parameters. */
    fgets(str,1000,fid);
    fscanf(fid,"%lg %lg %lg",&vp,&vs,&rho);
    fgets(str,1000,fid);
    
    /* Source coordinates (x, y, z). */
    fgets(str,1000,fid);
    for (int i=6; i<9; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=6; i<9; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);

    /* Moment tensor components. */
    fgets(str,1000,fid);
    for (int i=0; i<6; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=0; i<6; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);
    
    /* Origin time. */
    fgets(str,1000,fid);
    fscanf(fid,"%lg",&(mean_q[9]));
    fgets(str,1000,fid);
    fscanf(fid,"%lg",&(sigma_q[9]));
    fgets(str,1000,fid);

    /* Source time function coefficients. */
    fgets(str,1000,fid);
    for (int i=10; i<20; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=10; i<20; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);
    
    /* Assign values equal to prior means. */
    for (int i=0; i<20; i++) q[i]=mean_q[i];
    
    fclose(fid);
}

/* Set unit moment tensor. --------------------------------------------------------*/
void parameters::set_m(int idx)
{
    for (int i=0; i<6; i++) q[i]=0.0;
    q[idx]=1.0;
}

/* Set unit source-time function. -------------------------------------------------*/
void parameters::set_s(int idx)
{
    for (int i=10; i<20; i++) q[i]=0.0;
    q[idx]=1.0;
}

/* Addition. ----------------------------------------------------------------------*/

parameters operator+(const parameters &a, const parameters &b)
{
    parameters c;
    
    c.vp=a.vp;
    c.vs=a.vs;
    c.rho=a.rho;
    
    for (int i=0; i<20; i++)
    {
        c.q[i]=a.q[i]+b.q[i];
        c.mean_q[i]=a.mean_q[i];
        c.sigma_q[i]=a.sigma_q[i];
    }
    
    return c;
}

/* Subtraction. -------------------------------------------------------------------*/

parameters operator-(const parameters &a, const parameters &b)
{
    parameters c;
    
    c.vp=a.vp;
    c.vs=a.vs;
    c.rho=a.rho;
    
    for (int i=0; i<20; i++)
    {
        c.q[i]=a.q[i]-b.q[i];
        c.mean_q[i]=a.mean_q[i];
        c.sigma_q[i]=a.sigma_q[i];
    }
    
    return c;
}


/*=================================================================================*/
/* Seismic data class. ------------------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/
data::data()
{
    FILE *fid;
    char str[1000];
    
    /* Read time series setup. */
    fid=fopen("INPUT/setup.txt","r");
    
    fgets(str,1000,fid);
    fscanf(fid,"%d %lg %d",&nt,&dt,&nrec);
    fgets(str,1000,fid);
    
    recx=new double[nrec];
    recy=new double[nrec];
    recz=new double[nrec];
    
    fgets(str,1000,fid);
    for (int i=0; i<nrec; i++) fscanf(fid,"%lg %lg %lg",&(recx[i]),&(recy[i]),&(recz[i]));
    
    fclose(fid);
    
    /* Allocate memory for the time series and covariances. */
    ux=new double*[nrec]; uy=new double*[nrec]; uz=new double*[nrec];
    cx=new double*[nrec]; cy=new double*[nrec]; cz=new double*[nrec];

    for (int i=0; i<nrec; i++)
    {
        ux[i]=new double[nt]; uy[i]=new double[nt]; uz[i]=new double[nt];
        cx[i]=new double[nt]; cy[i]=new double[nt]; cz[i]=new double[nt];
    }
    
}

/* Destructor. --------------------------------------------------------------------*/
data::~data()
{
    if (recx) delete[] recx;
    if (recy) delete[] recy;
    if (recz) delete[] recz;
    
    if (ux)
    {
        for (int i=0; i<nrec; i++) delete[] ux[i];
        delete[] ux;
    }
    
    if (uy)
    {
        for (int i=0; i<nrec; i++) delete[] uy[i];
        delete[] uy;
    }
    
    if (uz)
    {
        for (int i=0; i<nrec; i++) delete[] uz[i];
        delete[] uz;
    }
    
    if (cx)
    {
        for (int i=0; i<nrec; i++) delete[] cx[i];
        delete[] cx;
    }
    
    if (cy)
    {
        for (int i=0; i<nrec; i++) delete[] cy[i];
        delete[] cy;
    }
    
    if (cz)
    {
        for (int i=0; i<nrec; i++) delete[] cz[i];
        delete[] cz;
    }
}

/* Fill time series arrays with artificial data for prior mean parameters. --------*/
void data::make_synthetics(parameters &q)
{
    double *outx, *outy, *outz;
    outx=new double[nt];
    outy=new double[nt];
    outz=new double[nt];
    
    for (int i=0; i<nrec; i++)
    {
        ffp(0,q,recx[i],recy[i],recz[i],dt,nt,outx);
        ffp(1,q,recx[i],recy[i],recz[i],dt,nt,outy);
        ffp(2,q,recx[i],recy[i],recz[i],dt,nt,outz);
        
        for (int j=0; j<nt; j++)
        {
            /* Synthetic time series. */
            ux[i][j]=outx[j];
            uy[i][j]=outy[j];
            uz[i][j]=outz[j];
            /* Synthetic inverse covariances zero. */
            cx[i][j]=0.0;
            cy[i][j]=0.0;
            cz[i][j]=0.0;
        }
    }
    
    delete[] outx; delete[] outy; delete[] outz;
}


/* Derivative of u. ---------------------------------------------------------------*/
void data::du(parameters &q, int component, int iq, int irec, double *out)
{
    /* Local variables. */
    double d=1000.0;    /* Finite-difference space increment [m]. */
    parameters qp=q;
    parameters qm=q;
    double *outp=new double[nt];
    double *outm=new double[nt];

    /* Derivatives with respect to moment tensor components. */
    if (iq<6)
    {
        qp.set_m(iq);
        ffp(component,qp,recx[irec],recy[irec],recz[irec],dt,nt,out);
    }
    /* Derivatives with respect to source location. */
    else if (iq<9 && iq>5)
    {
        qp.q[iq]+=d;
        qm.q[iq]-=d;
        ffp(component,qp,recx[irec],recy[irec],recz[irec],dt,nt,outp);
        ffp(component,qm,recx[irec],recy[irec],recz[irec],dt,nt,outm);
        
        for (int i=0; i<nt; i++) out[i]=(outp[i]-outm[i])/(2.0*d);
        
    }
    /* Derivative with respect to origin time. */
    else if (iq==9)
    {
        qp.q[9]+=dt;
        qm.q[9]-=dt;
        ffp(component,qp,recx[irec],recy[irec],recz[irec],dt,nt,outp);
        ffp(component,qm,recx[irec],recy[irec],recz[irec],dt,nt,outm);
        
        for (int i=0; i<nt; i++) out[i]=-(outp[i]-outm[i])/(2.0*dt);
    }
    /* Derivatives with respect to source-time function coefficient. */
    else if (iq>=10 && iq<20)
    {
        qp.set_s(iq);
        ffp(component,qp,recx[irec],recy[irec],recz[irec],dt,nt,out);
    }
    
    /* Clean up. */
    delete[] outp; delete[] outm;
}


/* Fill time series by reading data from a file. ----------------------------------*/
void data::read_data(const char *filename)
{
    /* Local variables. */
    FILE *fid;
    char str[1000];
    
    /* Open file and read header. */
    fid=fopen(filename,"r");
    
    fscanf(fid,"%d %lg %d",&nt,&dt,&nrec);
    fgets(str,1000,fid);
    
    /* Read traces. */
    for (int i=0; i<nrec; i++)
    {
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(ux[i][j]));
        fgets(str,1000,fid);
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(cx[i][j]));
        fgets(str,1000,fid);
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(uy[i][j]));
        fgets(str,1000,fid);
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(cy[i][j]));
        fgets(str,1000,fid);
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(uz[i][j]));
        fgets(str,1000,fid);
        for (int j=0; j<nt; j++) fscanf(fid,"%lg",&(cz[i][j]));
        fgets(str,1000,fid);
    }
    
    /* Clean up. */
    fclose(fid);
}

/* Write time series to a file. ---------------------------------------------------*/
void data::write(const char *filename)
{
    FILE *fid;
    
    fid=fopen(filename,"w");
    
    fprintf(fid,"%d %lg %d\n",nt,dt,nrec);
    
    for (int i=0; i<nrec; i++)
    {
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",ux[i][j]);
        fprintf(fid,"\n");
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",cx[i][j]);
        fprintf(fid,"\n");
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",uy[i][j]);
        fprintf(fid,"\n");
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",cy[i][j]);
        fprintf(fid,"\n");
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",uz[i][j]);
        fprintf(fid,"\n");
        for (int j=0; j<nt; j++) fprintf(fid,"%lg ",cz[i][j]);
        fprintf(fid,"\n");
    }
    
    fclose(fid);
}

/* Print data summary on screen. --------------------------------------------------*/
void data::print()
{
    printf("nt=%d, dt=%lg, nrec=%d\n",nt,dt,nrec);
    for (int i=0; i<nrec; i++)
    {
        printf("recx=%lg, recy=%lg, recz=%lg\n",recx[i],recy[i],recz[i]);
    }
}

/*=================================================================================*/
/* Moment rate function. ----------------------------------------------------------*/
/*=================================================================================*/

double mrf(parameters &q, double t)
{
    double s=0.0;
    
    /* March through parameters 11 to 20. */
    for (int i=10; i<20; i++) s+=q.q[i]*box(t-(double)(i-10));
    
    return s;
}

/*=================================================================================*/
/* Box function. ------------------------------------------------------------------*/
/*=================================================================================*/

double box(double t)
{
    if (t>= 0 && t<1.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

/*=================================================================================*/
/* Far-field P wave. --------------------------------------------------------------*/
/*=================================================================================*/

void ffp(int component, parameters &q, double recx, double recy, double recz, double dt, int nt, double *out)
{
    /* Local variables. */
    
    double t0, t=0.0;
    double Ap=1.0/(4.0*PI*q.rho*q.vp*q.vp*q.vp);    /* P wave amplitude. */
    double r, gx, gy, gz;                           /* Distance and direction cosines. */
    double A, Ax, Ay, Az;                           /* Radiation pattern. */
    
    double m[6], srcx, srcy, srcz;                  /* Moment tensor components and source coordinates. */
    
    /* Assign values. */
    for (int i=0; i<6; i++) m[i]=q.q[i];
    srcx=q.q[6];
    srcy=q.q[7];
    srcz=q.q[8];
    t0=q.q[9];
    
    /* Compute synthetics. */
    r=sqrt(pow(recx-srcx,2.0)+pow(recy-srcy,2.0)+pow(recz-srcz,2.0));
    gx=(recx-srcx)/r;
    gy=(recy-srcy)/r;
    gz=(recz-srcz)/r;
    
    A=1.0e15*Ap*(gx*gx*m[0]+gy*gy*m[1]+gz*gz*m[2]+2.0*gx*gy*m[3]+2.0*gx*gz*m[4]+2.0*gy*gz*m[5])/r;
    
    if (component==0) A=A*gx;
    else if (component==1) A=A*gy;
    else if (component==2) A=A*gz;
    
    for (int i=0; i<nt; i++)
    {
        out[i]=A*mrf(q,(t-t0)-r/q.vp);
        t+=dt;
    }
}

/*=================================================================================*/
/* Little helpers. ----------------------------------------------------------------*/
/*=================================================================================*/

/* Uniformly distributed, double-valued random numbers. ---------------------------*/
double randf(double min, double max)
{
    return (max-min)*(double)rand()/RAND_MAX+min;
}

/* Normally distributed, double-valued random numbers. ----------------------------*/
/* This function implements the Box-Muller transform to obtain a pair of
 normally distributed random numbers with a given mean and standard deviation. */
void randn(double mean, double stdv, double *x1, double *x2)
{
    double z1=(double)rand()/RAND_MAX;
    double z2=(double)rand()/RAND_MAX;
    
    *x1=sqrt(-2.0*log(z1))*cos(2.0*PI*z2);
    *x2=sqrt(-2.0*log(z1))*sin(2.0*PI*z2);
    
    *x1=stdv*(*x1)+mean;
    *x2=stdv*(*x2)+mean;
}

double randn(double mean, double stdv)
{
    double x;
    
    double z1=(double)rand()/RAND_MAX;
    double z2=(double)rand()/RAND_MAX;
    
    x=sqrt(-2.0*log(z1))*cos(2.0*PI*z2);
    x=stdv*x+mean;
    
    return x;
}


