//
//  jaco.h
//  cudaProj
//
//  Created by Tyler Chapman on 4/3/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef cudaProj_jaco_h
#define cudaProj_jaco_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
// Dark matter density models
#define MASS_HERNQUISTLIKE    1  // r^(-n) (1+r)^(n-4)
#define MASS_POWLAWWITHCORE   2  // (1+r)^(-n)
#define MASS_NFWWITHCORE      3  // (r0+r)^(-1) (r1+r)^(-2)
#define MASS_ISOWITHCORE      4  // (r0+r)^(-1) (r1+r)^(-1)
#define MASS_UNIVERSAL        5  // r^(-n) (1+r)^(3-n)
#define MASS_SERSIC           6  // exp(-(r/r1)^a)
#define MASS_FOLLOWSLIGHT     7  // exp(-(r/r1)^a)

#include <stdio.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <sciutils.h>
#include <hrothgar.h>
#include <params.h>

#include <sz_fitting.h>

#define DEFAULTCONFIGFILE "jaco.cfg"

#ifndef SZC
#define SZC "Compiled without SZ support.\n"
#endif

#define WELCOMESTRING "Joint Analysis of Cluster Observations v" VERSION "\nCopyright (C) 2009 Andisheh Mahdavi\n" SZC "Velocity dispersion code Copyright (C) 2009 Alison Mansheim, Andisheh Mahdavi\n" "\n"

// All defines regarding instruments should go here and ONLY here
// to facilitate adding instruments in the future.
#define NINSTR 8
#define MOS1 0
#define MOS2 1
#define PN 2
#define ACISI 3
#define ACISS 4
#define SZ 5
#define WL 6
#define VELS 7

// Projection types
#define PROJ_GAS    0
#define PROJ_DARK   1
#define PROJ_STARS  2
#define PROJ_ALL    3

//#ifdef JACO_INIT
static const char *idcode[] = {"m1","m2","pn","ai","as","wl","sz","vd"};
//#endif

// End of Instrument defines

// Physical constants
#define TMASS   4.4907052    // G in units of kev Mpc / 1e14  msun mp
#define GDYN    430241       // G in units of Mpc km^2 / (1e14 msun s^2)
#define SZNORM  0.004017065  // sigma_T Mpc keV / (cm^3 m_e c^2)
#define MPCM3   0.0040476904 // protonmass / cm^3 in units of 1e14 msun/Mpc^3
#define MPSQCM5 505551.16    // 1e-14 (1e14 msun)^2/Mpc^5 in units of mp^2/cm^5
#define RHOC    2.7749438E-7 // Critical density of a universe with
// hubble constant = 1 km / s / Mpc,
// in units of 1e14 Msun / Mpc^3
#define SIGMAC  16626.634    // c^2/(4 Pi G) in units of 1e14 msun / Mpc
#define GYR     3.1557e16    // Gigayear in seconds
#define KEV     1.60219E-9   // keV in ergs
#define MPC     3.0856776e24 // Mpc in cm

// Other constants
#define FOURPI   12.5663706
#define PI       3.14159265358979
#define ARCMINTORAD 0.0002908882086657216
#define MAXANN   1000        // Maximum number of annuli
#define INTMAX   5000        // Maximum # iterations for convergence
#define INFO_SIZE 15
#define DELTA_SIZE 8

// Number of contrasts to evaluate masses at.
#define NDELTA 4

#define DIFFER(x,y) (fabs((x)-(y)) > 1.E-6)

#ifdef BSD
#define FINITE finite
#else
#define FINITE isfinite
#endif

#define jaco_printf(format...) if (js->mpirank == 0) printf(format); 


//typedef float Real;
//char pf[4] = "%f";

struct particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float logtemp;
    float hsmooth;
    float metal;
    float logmetal;
    float tform;
    float eps;
    float phi;
    float opz[3];
    float nenp;
    float spec;
    int   type;
};

#define DARK 0
#define GAS  1
#define STAR 2


// Forward declaration; see below
struct jaco_vars;

// This structure will contain values that need to be globally
// accessible by all threads. These values should not be written
// to during parallel sections without a CRITICAL or ATOMIC block.
struct jaco_state {
    
    struct hrothgar_setup  *setup; // Link to hrothgar parameters
    
    double xrayra, xraydec;      // Cluster center in degrees
    char *jacofilename;          // List of input files
    char *psflist;               // List of PSF files
    char *profname;              // Profile file base name
    
    
    // Variables related to the cooling function
    char *spectralcode;          // File containing spectral code !!
    float ***mastercooling;      // Cooling function Lambda(E,T,Z) (CF) !!
    double *energyaxis;           // Energy axis of cooling function !!
    double *tempaxis;             // Temperature axis of cooling function !!
    double *metalaxis;            // Metallicity axis of cooling function !!
    double **bolometric;          // Bolometric flux at given temp and metal
    int egridsize;                  
    int tgridsize;               // Length of E,T,Z axes !!
    int mgridsize;
    double ***rebinnedcooling;   // Rebinned Lambda(Bin,T,Z) !!
    double ebinsize;             // Size of energy bins in keV
    long nlastbin;               // Total known number of energy bins
    double *lastbin;             // Centers of energy bins
    int   l1;                    // Beginning and ending bin for current
    int   l2;                    // spectrum
    int tempwarn;                // TRUE if we have warned that T/Z is being
    int metalwarn;               // extrapolated outside given bounds
    int entropywarning;          // Have we warned re: convective instability?
    int specwarn ;               // TRUE if spectrum was LARGE
    int binningchanged;          // Did CF binning change since last time?
    int nrebin;                  // # of times we've had to rebing the CF
    double efit1,efit2;          // Energy fitting range
    double specnorm;             // Spectrum normalization
    
    // Input miscellaneous parameters
    double syserr;               // X-ray systematic error
    double cutradius;            // Radius to cut out of analysis
    double contrast;             // Density contrast for mass evaluation
    double rshock;               // Termination shock in units of rContr
    double tshock;               // Efficiency of termination shock
    double redshift;             // Redshift of system
    double tprojpower;           // T-weighting for projected T profile. Cosmetic
    // only; the above NOT used when fitting
    double tbias;
    double hubble;               // Hubble constant in km/s/Mpc
    double omegal;               // Omega matter
    double omegam;               // Omega lambda
    int    projtype;             // Project GAS, DARK, STARS, or ALL
    
    int havecxo,havexmm;
    int calcxray;                // Whether to fit X-ray data
    int calcsz;                  // Whether to fit SZ data
    int calcwl;                  // Whether to fit weak lensing data
    int calcvel;                 // Whether to fit velocity dispersion data
    int samepars;                // Are physical params. same as last time?
    int changedpar[NTOTPARAMS];  // Which parameters changed?
    FILE  *logfile;              // Logfile for debugging
    int  mpirank;                // Rank of MPI CPU
    int  standalone;             // Whether we are running standlone or not
    int  nfitparams;             // # Parameters being fit.
    
    
    
    // Input metallicity profile parameters
    double z0,zinf;              // Metallicity at r=0 and r=infinity
    double rz;                   // Metallicity transition radius
    
    // Input gas profile parameters
    int gasmodel;                // GAS_TRIPLEBETA or GAS_BROKENPOWERLAW
    double Mg1Contr;             // 1st Gas mass at rContr (10^14 Msun)
    double Mg2Contr;             // 2nd Gas mass at rContr (10^14 Msun)
    double Mg3Contr;             // 3rd Gas mass at rContr (10^14 Msun)
    double Mg4Contr;             // 4th Gas mass at rContr (10^14 Msun)
    double alpha;                // Inner slope power law of gas model
    double b1;                   // Three slopes of the triple beta model
    double b2;                   // gas profile. Or, b1 and b2 are the broken
    double b3,b4;                // power law slopes if that is selected
    double rx1;                  // Transition radii for the triple beta model;
    double rx2;                  // for the broken power law, rx1 is used as
    double rx3;                  // the transition radius.
    double rx4;                  
    
    // Input dark profile parameters
    int darkmodel;               // Choice of a large number of dark profiles
    int totalmass;               // Dark mass refers to total or dark only mass?
    double MdContr;              // Dark mass at rContr (10^14 Msun)
    double ndark;                // Shape of profile (can have various meanings)
    double rdm0;                 // Dark matter core radius
    double rdm1;                 // Characteristic/scale radius of profile
    int fitmc;                   // Whether to interpret concentration as
    // actual concentration or normalization of
    // mass-concentration relation.
    
    // Input star profile parameters
    double MsContr;              // Stellar mass at rContr (10^14 Msun)
    double nstar;                // Slope of S'ersic profile (0.25=deVaucoul.)
    double rstar;                // Scale height: exp(-(r/rstar)^nstar)
    
    // Various useful derived and stored parameters
    double contr;                // Contrast for overdensity radius routine
    double opz;                  // One plus z
    double dist;                 // Comoving distance
    double angdist;              // Angular diameter distance
    double lumdistancesq;        // Square of luminosity distance
    double rhocrit;              // Critical density of universe at z
    double rContr;               // Radius at which density/rho_crit = contrast
    double Md;                   // Dark model normalization
    double gasdensity_norm;      // Gas density model normalization
    double Ms;                   // Star model normalization
    double nm1,nm2,nm3;          // ndark-1, 2, and 3
    double rdm3mn,rdm1sq,rdm1cu; // rdm1^(-nm3), rdm1^2, rdm1^3
    double rx1sq,rx2sq,rx3sq,rx4sq;    // rx(1,2,3)^2
    double rdm0sq,rdiff;         // rdm0^2,rdm0-rdm1
    double Mg1,Mg2,Mg3,Mg4;          // Gas mass model normalizations
    double rhog1,rhog2,rhog3,rhog4;    // Gas density model normalizations
    double pshock;               // Pressure at termination shock
    
    // Weak lensing information
    char *wldata;                // Weak lensing data file
    double wlbeta,wlbetasq;      // Weak lensing redshift-dependent correction
    double wllss;                // Lensing chi-square weight (LSS contribution)
    double wlcorrfac;            // fractors, see Hoekstra 2007
    double intR1sq,intR2sq;      // Weak lensing integration: temporary vars
    
    // Integration parameters
    double precision;            // Precision to use for integration
    int debug;                   // Whether to output debug information
    int nsplines;                // Number of splines for high-accuracy int.
    double *rspline;             // Radii of splines for high-accuracy int.
    int dofast;                  // Whether to do faster or higher-accuracy int.
    
    // Spectrum and Compton y storage variables
    double aR1[MAXANN];          // Inner and outer radii of stored
    double aR2[MAXANN];          // spectra
    double **emissivity;         // Emissivity(spectral bin, annulus)
    double **wgt_emissivity;     // Emissivity weighted by T^tprojpower
    double **unwgt_emissivity;   // Emissivity weighted by T^(tprojpower-1)
    double *pressure;            // Pressure at each annulus
    double *background;          // Unvignetted background
    char *backlist;              // List of residual background files
    char **backsyserr;
    struct pha residual_background[NINSTR];
    
    
    int nannuli;                 // Number of stored spectra
    
    // SZ information
    char *szdata;                // SZ data file name
    char *szconfig;              // List of SZ configuration files
    double szy[MAXANN];          // Tabulated Compton y
    double szrad[MAXANN];        // Radius at which it's tabulated
    int sz0;                     // First SZ l-value to fit
    int lastszannulus;
    
    // GSL storage space
    gsl_interp_accel *coolaccel; // CF interpolation accelerator
    gsl_root_fsolver *solver;    // Nonlinear solver for rContrast calculation
    gsl_interp **interp_coefs;   // Array of interpolators for emissivity  
    gsl_interp **wgt_interp_coefs; // Ditto for weighted emissivity
    gsl_interp 
    *szinterp_coefs;          // Interpolators for pressure
    
    
    // Velocity anisotropy parameters
    double aniso0,aniso1,raniso;
    char *veldata;
    
    // Simulation parameters
    void *rng;                   // Random number generator
    void *data;                  // Back pointer to data superstructure
    int dosim;                   // Whether to simulate or not. 0=no,
    // 1=add noise to JACO model; 2=read data
    // from N-body simulation; 3=also create
    // FITS image of the model
    char *simfile,*simpars;      // Simulation data & parameter file
    struct particle *ptype[3];   // Pointer to gas, dark, and star particles
    struct particle *part;       // All particles from N-body simulation
    struct particle *gas;        // Gas particles from N-body simulation
    struct particle *star;       // star particles
    struct particle *dark;       // dark matter particles
    int axis;                    // Axis along which to do the observation
    long simcount, np[3];        // # of particles of each type.
    double **simmofr;            // Simulated mass profile
    size_t *order2D,*order3D;    // Particles sorted by distance
    double *sortdist2D,*sortdist3D;
    float boxmin[3][3],
    boxmax[3][3];              // Size of each simulation box
    double **center;             // X,Y,Z centers for 3 
    double mpc,massunit,velunit; // Basic units in which simulation was done
    double **distsq;             // Matrix containing projected/3D distances 
    char *simimage;              // Input image to guide simulated outputs
    unsigned long nx,ny;         // Simulated image dimensions
    float pixscale;              // Pixel scale of simulated image
    
    // Triaxial mode parameters
    int n_ell;
    double ell_max;
    double *Tprof,*Zprof,*nenhprof;
    double a_ell,b_ell;
    int ***corners;
    int n_box_spectra;
    double theta;
    double phi;
    double epsilon;
    
    // Nongravity mode parameters
    
    int Tpar[NTOTPARAMS];
    int Zpar[NTOTPARAMS];
    int nepar[NTOTPARAMS];
    double **singleTspec;
    double **geomfactors;
    
    // Chosen gas, dark, and stellar profiles
    double (*gasdensity)(double,struct jaco_state *js);
    double (*gasmass)(double,struct jaco_state *js);
    double (*darkmass)(double,struct jaco_state *js);
    double (*darkdens)(double,struct jaco_state *js);
    double (*stellarmass)(double,struct jaco_state *js);
    double (*stellardensity)(double,struct jaco_state *js);
    double (*mass)(double,struct jaco_state *js);
    
    // Integrator (fast, slow, shocked, or non-shocked [= to infinite radius] )
    double (*ninteg)(double (*func)(double r, void *par),double R,double Rb,
                     struct jaco_vars *jv);
    
    BigSZControlParams szparams;
    
    int  fitwrite;
    char **fitsname,**fitcomp;
    int  xraycount;
    
    // Extra parameters to be tracked along fit parameters
    int ninfo;              // # of info parameters
    char **infoname;        // Names of info parameters
    double *infoarray;      // List of info values.
    
    double  *deltas;             // Contrasts at which to output.
    int     ndelta;             // # of contrasts
    int *annuluscount;
    
    double *r1list,*r2list,*ringarea,*effrad,*r1full,*r2full,*goodarea;
    int **annulusid,**psfrow,*instrument;
    double **psf;
    int *annuli, ebins;
};

// This structures contains values that change throughout the integration
// process. A separate copy needs to be allocated for each thread.
struct jaco_vars {
    
    int  calctproj;             // Whether to calculate weighted temperature
    char integstate[1000];      // Informative string for integration errors
    double **coolmatrix;        // Selects the correct CF matrix
    double *szyinteg;           // Pressure to integrate along line of sight
    double Hfrac;               // Fraction of gas that's hydrogen by mass
    double mu;                  // Mean molecular weight of gas
    double nenh;                // number of electrons per hydrogen nucleus
    double breakradius;         // Radius at which to split the integration
    double R1deg,R2deg;         // Inner and outer annulus radii in degrees
    double R1sq,R2sq;           // Square of R1 and R2 (in Mpc^2)
    long specbin;               // Spectral bin currently being integrated
    int integcounter;           // Radial integration loop variable
    int nabs;                   // Number of abscissae for Gauss-Legendre
    double *xabs;               // Abscissae for Gauss-Legendre
    double *kernel;             // Projection factor for Gauss-Legendre
    double *tint;               // Temperature profile for Gauss-Legendre
    double *logT;               // log Temperature profile for Gauss-Legendre
    double *logZ;               // log Metallicity  profile for Gauss-Legendre
    double *entropy;
    double R1,R2,Reff;          // Inner, outer, effective annulus radius
    double bolocorr;            // Bolometric flux correction.
    
    gsl_integration_workspace   // High-precision integration workspace
    *workspace[2];
    int depth;                  // High-precision integration depth
    gsl_interp_accel *accel;    // Emissivity interpolation accelerator
    gsl_interp_accel *szaccel;  // SZ Pressure interpolator
    
    struct jaco_state *state;   // Carry over the static variables as well
};

// Shortcuts for defining local variables named after structure members
// and for transferring their values back to the struct.
#define JVS(x) typeof(jv->state->x) x = jv->state->x
#define JS(x) typeof(js->x) x = js->x
#define JV(x) typeof(jv->x) x = jv->x
#define JSET(x,y) js->x = y; typeof(js->x) x; x = js->x
#define JRESET(x,y) js->x = y; x = js->x

#define GAS_TRIPLEBETA        2  // Sum of two beta models.
#define GAS_BROKENPOWERLAW    1  // Broken power law.


void  jaco_xray(int annulus, int nbin, double *spectrum, struct jaco_vars *jv);
void read_cooling_function(struct jaco_state *js);
double cooling(double logT, double logZ, struct jaco_vars *jv);
void deallocatevectors();
int rebincoolingfunction(double *bincenter, int nbin, struct jaco_state *js);
void generate_ymap(struct jaco_state *js);
int  initialize_wl(struct jaco_state *js);
void tangential_shear(double *R, double *shear, int nrad, 
                      struct jaco_state *js);
void jaco_update_models(double *params, struct jaco_state *js);
int photoabs(double nH, double
             *ener, double dener, int ne, double *abscoef);


void jaco_init(struct jaco_state *js, 
               struct hrothgar_setup *setup, int argc, char **argv);
int jaco_calc_models(double *params, int *ndata, double *outValues,
                     double *x0, struct jaco_state *js);
int jaco_read_data(struct jaco_state *js);
void jaco_thread_init(struct jaco_vars *jv);
void jaco_thread_shutdown(struct jaco_vars *jv);
void emissivity_profile(struct jaco_state *js);

void simulate_spec(int annulus, int nbin, double *spectrum, struct jaco_vars *jv);
int init_sim(struct jaco_state *js);
void make_simimages(struct jaco_state *js);
void dispersion_profile(double *Rvec, double *profile, int ndata, 
                        struct jaco_vars *jv);
double velocity_dispersion(double R, struct jaco_vars * jv);
void make_xrayimage(struct jaco_state *js);
double metallicity(double r, struct jaco_vars *jv);



#endif


