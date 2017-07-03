#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mc.hpp"
#include "aux.hpp"

/* ./sampling metropolis iterations [verbose] */
/* ./sampling hamilton iterations nt dt [verbose] */

int main(int argc, char *argv[])
{
    /* Check input. -------------------------------------------------------------------*/
    
    int nt, iterations;
    double dt;
    bool verbose;
    
    if (!strcmp(argv[1],"metropolis"))
    {
        nt=0;
        dt=0.0;
        iterations=atoi(argv[2]);
        if (argv[3] && !strcmp(argv[3],"verbose")) verbose=true;
    }
    else if (!strcmp(argv[1],"hamilton"))
    {
        iterations=atoi(argv[2]);
        nt=atoi(argv[3]);
        dt=atof(argv[4]);
        if (argv[5] && !strcmp(argv[5],"verbose")) verbose=true;
    }
    
    /* Local variables. ---------------------------------------------------------------*/
    
    mc m(iterations,nt,dt,verbose);
    parameters params;
    clock_t start=clock();
    int accepted=0;
    double x, x_new;
    
    FILE *pfile;
    pfile=fopen("OUTPUT/samples.txt","w");
    
    /* Initial values. ----------------------------------------------------------------*/
    
    if (!strcmp(argv[1],"metropolis")) x=m.chi();
    else if (!strcmp(argv[1],"hamilton")) x=m.energy();
    
    m.write_sample(pfile,x,0);
    
    /* Random walk. -------------------------------------------------------------------*/
    
    for (int it=0; it<m.iterations; it++)
    {
        /* Make a model proposition and compute misfit/energy. */
        if (!strcmp(argv[1],"metropolis"))
        {
            m.propose_metropolis();
            x_new=m.chi();
        }
        else if (!strcmp(argv[1],"hamilton"))
        {
            m.propose_hamilton();
            x_new=m.energy();
        }
        
        /* Check Metropolis rule. */
        if ((x_new<x) || (exp(x-x_new)>randf(0.0,1.0)))
        {
            x=x_new;
            m.q=m.q_new;
            accepted++;
        }
        
        m.write_sample(pfile,x,it+1);
    }

    printf("accepted: %d\n",accepted);
    printf("elapsed time: %f\n",(double)(clock()-start)/CLOCKS_PER_SEC);
    
    /* Clean up. ----------------------------------------------------------------------*/
    
    fclose(pfile);
}

