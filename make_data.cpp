
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aux.hpp"

int main()
{
    /* Local variables. */
    data d;
    parameters q;
    
    /* Compute artificial data. */
    q.read_input("INPUT/parameters.txt");
    d.make_synthetics(q);
    
    /* Compute artificial inverse covariances. */
    for (int i=0; i<d.nrec; i++)
    {
        for (int j=0; j<d.nt; j++)
        {
            d.cx[i][j]=4.0e14;
            d.cy[i][j]=4.0e14;
            d.cz[i][j]=4.0e14;
        }
    }
    
    /* Write to file. */
    d.write("DATA/data.txt");
}
