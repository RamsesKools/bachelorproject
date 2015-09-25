#include <math.h>

/*------------------INLINE SUBSTITUTIONS-------------------------------------*/

//Binary definitions
#define TRUE 1
#define FALSE 0
#define CANCELLED -99

//not so necessary defitions
#define P_INDEX(ptr)  (ptr - pts)
#define F_COMP(a,b) (fabs( a-b ) < 1e-10)
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define MESSAGE(a) printf("Message:"#a"\n")
#define SIGN(a)   ( 2*(a>0) -1)

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })
#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

//Math shorthands + special numbers
#define PI           3.141592653589793
#define HALFPI       1.5707963267948966
#define SQRT2        1.4142136
#define BIGNUM       1e99
#define EPS          1e-10
#define NVT          1
#define cubic(x) ((x)*(x)*(x))
#define square(x) ((x)*(x))

//Print shorthands
#define dprint(expr) printf(#expr " = %d\n",expr)
#define gprint(expr) printf(#expr " = %g\n",expr)
#define vprint(expr) printf(#expr " = ( %16.14g %16.14g %16.14g ) \n",expr.x, expr.y, expr.z)

//Vector definitions
#define vector_inp(a, b) (a.x * b.x + a.y * b.y + a.z * b.z )
#define vector_cross(a,b,h)   h.x = a.y * b.z - a.z * b.y; h.y = a.z * b.x - a.x * b.z; h.z = a.x * b.y - a.y * b.x;
#define vector_times(a,b,h)   h.x = a.x * b.x; h.y = a.y * b.y; h.z = a.z * b.z;
#define vector_divide(a,b,h)   h.x = a.x / b.x; h.y = a.y / b.y; h.z = a.z / b.z;
#define vector_plustimes(a,b,h)   h.x += a.x * b.x; h.y += a.y * b.y; h.z += a.z * b.z;
#define vector_mintimes(a,b,h)   h.x -= a.x * b.x; h.y -= a.y * b.y; h.z -= a.z * b.z;
#define scalar_times(a,b,h)  h.x = a.x *b; h.y =a.y * b; h.z = a.z * b;
#define scalar_divide(a,b,h)  h.x = a.x /b; h.y =a.y / b; h.z = a.z / b;
#define scalar_plustimes(a,b,h)  h.x += a.x *b; h.y +=a.y * b; h.z += a.z * b;
#define scalar_mintimes(a,b,h)  h.x -= a.x *b; h.y -=a.y * b; h.z -= a.z * b;
#define vector_add(a,b,h)  {h.x = a.x +b.x; h.y =a.y + b.y; h.z = a.z + b.z;  }
#define vector_minus(a,b,h) { h.x = a.x - b.x; h.y =a.y - b.y; h.z = a.z - b.z;  }
#define matrix_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.y = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.z = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}
#define matrixT_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.y.x*a.y + m.z.x*a.z;\
    h.y = m.x.y*a.x + m.y.y*a.y + m.z.y*a.z;\
    h.z = m.x.z*a.x + m.y.z*a.y + m.z.z*a.z;}
#define normvec(a,b) {b.x=a.x/sqrt(vector_inp(a,a)); b.y=a.y/sqrt(vector_inp(a,a)); b.z=a.z/sqrt(vector_inp(a,a));}

//Quaternion definitions
#define quat_add(a,b,h) { h.q0=a.q0+b.q0; h.q1=a.q1+b.q1; h.q2=a.q2+b.q2; h.q3=a.q3+b.q3; }
#define quat_minus(a,b,h) { h.q0=a.q0-b.q0; h.q1=a.q1-b.q1; h.q2=a.q2-b.q2; h.q3=a.q3-b.q3;}
#define quat_inp(a,b) (a.q0*b.q0 + a.q1*b.q1 + a.q2*b.q2 + a.q3*b.q3)
#define sctimes_quat(a,b,h) { h.q0=a.q0*b; h.q1=a.q1*b; h.q2=a.q2*b; h.q3=a.q3*b;}
#define scdivide_quat(a,b,h) { h.q0=a.q0/b; h.q2=a.q2/b; h.q3=a.q3/b; h.q3=a.q3/b;}
#define quat_times(a,b,h) {\
    h.q0 = a.q0*b.q0 - a.q1*b.q1 - a.q2*b.q2 - a.q3*b.q3;\
    h.q1 = a.q0*b.q1 + a.q1*b.q0 + a.q2*b.q3 - a.q3*b.q2;\
    h.q2 = a.q0*b.q2 - a.q1*b.q3 + a.q2*b.q0 + a.q3*b.q1;\
    h.q3 = a.q0*b.q3 + a.q1*b.q2 - a.q2*b.q1 + a.q3*b.q0;}
#define quat_inverse(a,h) {h.q0 = a.q0; h.q1 = -a.q1; h.q2 = -a.q2; h.q3 = -a.q3;}
#define quatmatrix_x_vec(m,a,h) {\
    h.q0 = m.w.x*a.x + m.w.y*a.y + m.w.z*a.z;\
    h.q1 = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.q2 = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.q3 = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}

//Update definitions
#define update_average(a,b) {a.now =b; a.sum += a.now; a.sumsq += a.now*a.now; a.n++;}
#define update_blockaver(a) {a.sum /= a.n; a.sumsq = sqrt(a.sumsq/a.n - a.sum*a.sum);}
#define update_finaver(a,b) {a.sum += b.sum; a.sumsq += b.sum*b.sum;  b.sum = 0; b.sumsq = 0; b.now=0; b.n=0;}
#define update_blockmc(a,b) {a.acc += b.acc; a.tries += b.tries;  b.acc = 0; b.tries = 0;}


#define SWAP(a,b,c) {c=a;a=b;b=c;}

//Periodic Boundary conditions
#define pbc(a,b) {\
    a.x -= b.x*rint(a.x/b.x);\
    a.y -= b.y*rint(a.y/b.y);\
    a.z -= b.z*rint(a.z/b.z);}

#define NPART 2
#define NSITES 20
#define MAXSLICES 500000
#define MAXSTATES 3
#define MAXREPLICA 15
#define MAXSETS 0
#define NACC 4
#define NSTAT 2
#define MAXBIN 10000

/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct vector_type {

    double        x,y,z;

} vector;


typedef struct quaternion_type {

    double        q0,q1,q2,q3;

} quaternion;


typedef struct tensor_type {

    vector        x,y,z;

} tensor;


typedef struct quattensor_type {

    vector        w,x,y,z;

} quattensor;


typedef struct particle_type {

    vector        r,
                  f,
                  t;

    quaternion    q;

} Pts;


typedef struct site_type {

    vector        r;

    double        eps,
                  delta;

} Site; // hoe moet deze struct heten? toegevoegd 10-9-15


typedef struct bonds_mat {

    double        mat[NPART][NPART][NSITES][NSITES];

} Bonds;


typedef struct hist1d_type {

    double        *bin,
                  dbin,
                  offset,
                  weight;
    int           maxbin;

    char          name[100];

} Hist1d;


typedef struct hist2d_type {

    double        **bin,
                  dbin1,
                  dbin2,
                  offset1,
                  offset2,
                  weight;
    int           maxbin1,
                  maxbin2;

    char          name[100];

} Hist2d;



typedef struct dop_type {

    int           **flag_dop,
                  *flag_norm,
                  npaths,
                  maxbin1,
                  maxbin2,
                  type;

    double        **dop,
                  ***vf,
                  **escape,
                  *norm,
                  **dopen,
                  *dopennorm,
                  weight,
                  offset1,
                  offset2,
                  dbin1,
                  dbin2;

    char          name[40];

    FILE          *pathAB;

} Dop;


typedef struct aver_type {

    char          name[100];

    double        now,
                  n,
                  sum,
                  sumsq;

} Average;


typedef struct mcacc_type {

    char          name[100];

    int           acc,
                  tries;

    double        ratio;

} Mcacc;


typedef struct stats_type {

    Average       aver[NSTAT];

    Mcacc         mcacc[NACC];

} Stats;



typedef struct replica_type {


    double        lambda,
                  weight,
                  logcrossprob,
                  lnweight,
                  dos;

    int           pathlen,
                  index,
                  swapindex,
                  ntotal,
                  navlen,
                  type;

    long int      avlen,
                  avlensq;

    Dop           dop;

    Hist2d        fe;//,commit[MAXSTATES][MAXSTATES];

    vector        string;


    //FILE          *pathfilep[2];

} Replica;



typedef struct slice_type {

   Pts           pts[NPART];

    double        order_parameter,
                  rijdist,
                  energyN,
                  energyT,
                  energy;

    int           minisite,
                  minjsite;

} Slice;



typedef struct state_type {

    Replica       srep[MAXREPLICA];

    double        weight,
                  lnweight,
                  dos,
                  crosshist[MAXREPLICA][MAXBIN],
                  pathtypenumbers[MAXREPLICA][MAXSTATES][MAXREPLICA],
                  flux0,
                  flux1,
                  rate[MAXSTATES],
                  lambda[MAXREPLICA],
                  logcrossprob,
                  min,
                  mindist,
                  scalefactor,
                  volume_op;

    int           nrep,
                  type_mat[MAXREPLICA][2],
                  mstis_mat[MAXSTATES],
                  maxpaths,
                  nflux0,
                  nflux1,
                  fcount,
                  n;

    Slice         target;

} State;



typedef struct langevin_type {

    int           ninter;

    double        timestep,
                  dtD,
                  dtBeta,
                  friction;

} Langevin;



typedef struct path_type {

    int           nslices,
                  ninter,
                  type,
                  nstates,
                  nshoot,
                  nreverse,
                  nrepswap,
                  current_initial_state,
                  current_final_state,
                  nswapstates,
                  initial_state,
                  final_state,
                  current_replica,
                  current_set,
                  stateswapbias,
                  nreplica,
                  fixedbias,
                  maxlength;

    double        scalefactor,
                  current_gsen,
                  enbond,
                  energy;

    FILE          *filepath[MAXSTATES][MAXSTATES],
                  *filepatherr[MAXSTATES],
                  *fileswap[MAXSTATES],
                  *filestateswap,
                  *fprc;

    Stats         block_stats[MAXSTATES][MAXREPLICA],
                  final_stats[MAXSTATES][MAXREPLICA];


} Path;



typedef struct system_type {

    int           ncycle1,
                  ncycle2,
                  npart,
                  nsites,
                  //nsitesN, //misschien later nog toevoegen
                  start_type,
                  sim_type,
                  freq_graphics,
                  graphics;

    double        energy,
                  beta,
                  temp,
                  delta,
                  deltaN,
                  cosdelta,
                  mobilityT,
                  sqrtmobilityT,
                  mobilityR,
                  sqrtmobilityR,
                  oneover_cosdelta,
                  lambda,
                  fSE,
                  sigma,
                  sigmaLJ,
                  sigmaLJsq,
                  epsilonC,
                  epsilonP,
                  epsilonN;

    vector        boxl;
    Site          site[NSITES];

    FILE          *filep;

} System;


/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED VARIABLES-------------------------------*/

extern Slice         *slice,*trial,*fwtrial,*bwtrial;
extern System        sys;
extern State         state[MAXSTATES];
extern Path          path;
extern vector        sites[NPART][NSITES];
extern Replica       *replica[MAXREPLICA];
extern Langevin      langevin;
extern Stats         nulstat;
extern vector        nulvec;

/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/

extern void setup_simulation( void );
extern void setup_positions( Slice * );
extern void read_input(void);
extern void print_input(void);
extern void init_model(void);
extern double potential_energy( Slice * );
extern int check_bond(Pts *, Pts *, int *);
extern void calculate_forces( Slice * );
extern void propagate_bd( Slice * );
extern tensor getrotmatrix(quaternion);
extern quattensor getquatmatrix(quaternion);
extern double langrange_multiplier_quat(quaternion, quaternion);
extern void sampleangles(int);
extern void samplediffusion(int);
extern void bmdcycle();
extern void tiscycle();
extern void printstatistis();

extern void print_tisstats();
extern void reweight_crosshist();
extern void reweight_pathtypenumbers();
extern void traj_output();
extern void traj_input();
extern void dos_input();
extern void dos_output();
extern void sites_output();
extern void sites_input();
extern void warm_start_gen();

extern void update_wanglandau();
extern void crossinghistogram(Replica *);
extern void get_flux(Replica *);
extern void get_omiflux(Replica *);
extern void write_swap();
extern void write_paths();
extern void write_errorstate(Slice *, int);
extern void write_ratematrix();

extern void conf_output(Slice *);
extern void conf_input(Slice *);

extern void readslice(Slice *);
extern int readpath();

extern void terminate_block();
extern void finalstat();
extern void print_pathacceptances();
extern void print_finpathacceptances();
extern void print_pathnumbers();

extern int clusteranalysis(Slice *);
extern void create_all_rc(Slice *);
extern double get_rc(Slice *);
extern double print_rc(Slice *,int);
extern int in_state(Slice *);
extern int trajectory(Slice *, int);
extern int in_upper_window(Slice *, Replica *, int);
extern int trajectory_state_i(Slice *, Replica *, int, int);
extern int analyse(Slice *, Replica *, int, int);
extern int analyse_state_i(Slice *, Replica *, int, int);
extern int shoot_oneway(Replica *);
extern int swap_replica(int, int);
extern int swap_replica_0(int, int);
extern int swap_states(int, int);
extern int reverse_replica(int);

extern void state_init();
extern void replica_init();
extern void read_lambda();
extern quaternion quatVecToVec(vector, vector);
extern void print_tis_setup();
extern void stats_init();

extern void printstatusbmd();
extern void printstatustis();

extern void Init_Graphics();
extern void mainloop_for_graphics();
extern void reset_center(Slice *);

extern double RandomNumber(void);
extern void InitializeRandomNumberGenerator(double);
extern double RandomGaussianNumber();
extern vector RandomBrownianVector(double);
extern vector RandomVector(double);
extern double RandomVelocity(double);
extern vector RandomUnitVector(void);
extern quaternion RandomQuaternion();
extern quaternion RandomQuaternionRange(double);

