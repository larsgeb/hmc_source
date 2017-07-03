
#ifndef aux_hpp
#define aux_hpp

/*=================================================================================*/
/* Modelling parameters. ----------------------------------------------------------*/
/*=================================================================================*/

class parameters
{
public:
    
    /** Modelling parameters. */
    
    double vp, vs, rho;                             /**< vp [m/s], vs [m/s], rho [kg/m**3]. */
    
    double q[20], sigma_q[20], mean_q[20];          /**< Inversion parameters (values, prior stdev, prior mean).
                                                     q[0-5]: moment tensor components, mxx, myy, mzz, mxy, mxz, mzz,
                                                     q[6-8]: source coordinates, x, y, z,
                                                     q[9]: origin time,
                                                     q[10-19]: source-time function parameters. */
    
    /** Constructor and destructor. */
    
    parameters();         /**< Constructor. */
    ~parameters();        /**< Destructor. */
    parameters(const parameters &a);                            // Copy constructor.
    parameters &operator=(const parameters &a);                 // Assignment operator.

    /** Member functions. */
    
    void read_input(const char *filename);  /**< Fill values by reading input file. */
    void set_m(int idx);                    /**< Set all moment tensor components to 0.0, except m[idx]=1.0. */
    void set_s(int idx);                    /**< Set all source-time function coefficients to 0.0, except s[idx]=1.0. */
};

/* Operators. ---------------------------------------------------------------------*/

parameters operator+(const parameters &a, const parameters &b);
parameters operator-(const parameters &a, const parameters &b);

/*=================================================================================*/
/* Seismic data class. ------------------------------------------------------------*/
/*=================================================================================*/

class data
{
public:
    
    /** Time series setup. */
    
    int nt;         /**< Number of time samples. */
    double dt;      /**< Time increment [s]. */
    int nrec;       /**< Number of receivers. */
    
    double *recx, *recy, *recz;         /**< x,y,z-coordinates of receivers. */
    double **ux, **uy, **uz;            /**< x,y,z-components displacement [m]. */
    double **cx, **cy, **cz;            /**< x,y,z-components of inverse displacement covariance. */
    
    /** Constructor and destructor. */
    
    data();         /**< Constructor. */
    ~data();        /**< Destructor. */
    
    /** Member functions. */
    
    /**< Derivative of ux evaluated at the prior mean of the parameters. */
    void du(
             parameters &q,                         /**< Parameters for synthetics. */
             int component,                         /**< Component, 0=x, 1=y, 2=z. */
             int iq,                                /**< Index of parameter. */
             int irec,                              /**< Index of receiver. */
             double *out                            /**< Output. Derivative of x-component. */
    );
    
    /**< Compute synthetics ux, uy, uz for the prior mean values of q.q. */
    void make_synthetics(parameters &q);            /**< Fill the time series by computing synthetics for parameters q. */
    void read_data(const char *filename);           /**< Fill time series by reading data from a file. */
    
    void write(const char *filename);               /**< Write data to a file. */
    void print();                                   /**< Print data summary on screen. */
    
};

/*=================================================================================*/
/* Moment rate function. ----------------------------------------------------------*/
/*=================================================================================*/

double mrf(
            parameters &q,          /**< Parameters for synthetics. */
            double t                /**< Time [s]. */
           );

/*=================================================================================*/
/* Box function. ------------------------------------------------------------------*/
/*=================================================================================*/

double box(
            double t    /**< Time [s]. */
            );


/*=================================================================================*/
/* Far-field P wave. --------------------------------------------------------------*/
/*=================================================================================*/

void ffp(
            int component,                                  /**< Component, 0=x, 1=y, 2=z. */
            parameters &q,                                  /**< Modelling parameters. */
            double recx, double recy, double recz,          /**< Receiver coordinates. */
            double dt, int nt,                              /**< Time increment and number of samples. */
            double *out                                     /**< Output vector. */
);

/*=================================================================================*/
/* Little helpers. ----------------------------------------------------------------*/
/*=================================================================================*/

/** Double-valued, uniformly distributed random numbers. */
double randf(
             double min,    /**< Minimum value. */
             double max     /**< Maximum value. */
);

/** Double-valued, normally distributed random numbers. */
void randn(
           double mean,        /**< Mean. */
           double stdv,        /**< Standard deviation. */
           double *x1,         /**< Pointer to first random number. */
           double *x2          /**< Pointer to second random number. */
);

double randn(
             double mean,       /**< Mean. */
             double stdv        /**< Standard deviation. */
);

#endif /* aux_hpp */
