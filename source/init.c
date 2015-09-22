#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"



void setup_simulation() {

    int istate;
    printf("Setting up the system\n");

    InitializeRandomNumberGenerator(time(0l));

    read_input();
    init_model();

    printf("Allocating memory for the pathways\n");
    //for bmdcycle only initial positions need to be set
    if(sys.sim_type==0) {
        slice = (Slice *)calloc(1,sizeof(Slice));
        setup_positions(&slice[0]);
        printf("Initial Energy system %lf\n",potential_energy(&slice[0]));
    }
    else if(sys.sim_type==1) {
        slice = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        trial = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        fwtrial = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        bwtrial = (Slice *)calloc(MAXSLICES,sizeof(Slice));

        setup_positions(&slice[0]);
        read_lambda();
        state_init();
        replica_init();
        if(sys.start_type==1) {
            traj_input();
            dos_input();
        }
        stats_init();

        print_tis_setup();
    }

    printf("Setting up the system is done\n");

    return;
}


void init_model() {

    int isite;
    double psifactor,psi;

    sys.boxl.y = sys.boxl.z = sys.boxl.x;

    sys.temp = 1.0/sys.beta;

    sys.sigmaLJ = pow(2.0,(1.0/6.0));
    printf("sigmaLJ   %lf\n",sys.sigmaLJ);
    sys.sigmaLJsq = sys.sigmaLJ*sys.sigmaLJ;
    printf("sigmaLJsq %lf\n",sys.sigmaLJsq);
    sys.delta=PI*sys.delta/180.0;
    sys.cosdelta=cos(sys.delta);
    sys.oneover_cosdelta = 1.0/(1.0-sys.cosdelta);
    langevin.dtD=sqrt(2.0*langevin.timestep);
    langevin.dtBeta=langevin.timestep*sys.beta;

    sys.sqrtmobilityT = sqrt(sys.mobilityT);
    sys.sqrtmobilityR = sqrt(sys.mobilityR);
    
    if(sys.nsites>0) {
        sys.site[0].r.x=0.;
        sys.site[0].r.y=0.;
        sys.site[0].r.z=1.;
        
        sys.site[0].eps=sys.epsilonP;
        sys.site[0].delta=sys.delta;
    }  

    //OTHER patch definitions should be done randomly
    //Possible create separate structure Site which contains: position patch, width, epsilon
    //  then every patch can be different in character


    // sites 1 to nsites will be given a random vector
    
    int overlap, jsite, failsafe_count;
    for(isite=1; isite<sys.nsites; isite++) {
        do {
            overlap = 0;
            sys.site[isite].r = RandomUnitVector();   
            for(jsite=0; jsite<isite; jsite++) {
                if(vector_inp(sys.site[isite].r ,sys.site[jsite].r) > sys.cosdelta) {
                    overlap = 1;
                    failsafe_count ++;
                    if(failsafe_count > 100000) {
                        exit(1);
                    }                    
                    break;
                }
            }

            // deze regel stond er al maar ik weet niet wat die doet
            //printf("length patch vector %d: %lf\n",isite, vector_inp(sys.site[isite],sys.site[isite]));
        } while (overlap == 1);

        sys.site[isite].eps=sys.epsilonN;
        sys.site[isite].delta=sys.deltaN;
     }


    return;
}




void setup_positions(Slice *psl) {

    double boxlength=sys.boxl.x, r2;
    int ipart, jpart, overlap;
    vector dr;
    Pts *psi, *psj;

    printf("Setting the positions for all the particles\n");
    //random configuration
    //
    if(sys.npart==1) {
        printf("Only one particle: placing it at the centre\n");
        psi=&psl->pts[0];
        psi->r=nulvec;
        psi->q=RandomQuaternion();
    }

    else if(sys.npart==2) {
        printf("Only two particles: placing them towards each other\n");
        for( ipart=0; ipart<sys.npart; ipart++ ) {
            psi=&psl->pts[ipart];
            psi->r=nulvec;
            psi->r.z -= 0.5 - sys.sigmaLJ*ipart*sys.sigma;
            //psi->q=RandomQuaternion();
            if(ipart==0) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 0.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 1.0;
            }
            if(ipart==1) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 1.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
        }
    }
    


    if(sys.sim_type==0) {
        if(sys.start_type==1) {
            printf("Reading from conf.inp\n");
            conf_input(&slice[0]);
        }
    }


    if(sys.start_type==3) {
        printf("Randomly placing particles with random orientation\n");
        for( ipart=0; ipart<sys.npart; ipart++) {
            psi=&psl->pts[ipart];
            do {
                overlap=0;
                psi->r=RandomVector(boxlength);
                psi->q=RandomQuaternion();
                for( jpart=0; jpart<ipart; jpart++) {
                    psj=&psl->pts[jpart];
                    vector_minus(psj->r,psi->r,dr);
                    pbc(dr,sys.boxl);
                    r2=vector_inp(dr,dr);
                    if(r2<=sys.sigmaLJsq) {
                        overlap=1;
                    }
                }
            } while(overlap==1);
        }
    }

    printf("Setting the positions for all the particles is done\n");
    return;
}


#define MAXLINE 100
void read_input() {

    FILE *fp;
    char *pt,line[MAXLINE];

    if((fp = fopen("path.inp","r"))==NULL) {
        printf("ERROR: could not read input file\n");
        exit(1);
    }

    printf("\nReading input from path.inp\n");

    while(fgets(line,MAXLINE, fp) != NULL) {
        pt = strtok(line," ");
        if(strcmp(pt,"ncycle1")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle1);
        } else if( strcmp(pt,"ncycle2")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle2);
        } else if( strcmp(pt,"npart")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.npart);
        } else if( strcmp(pt,"graphics")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.graphics);
        } else if( strcmp(pt,"delta")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.delta);
        } else if( strcmp(pt,"deltaN")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.deltaN);
        } else if( strcmp(pt,"epsilonN")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonN);
        } else if( strcmp(pt,"nsites")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.nsites);
        //} else if( strcmp(pt,"nsitesN")==0) { //deze regel leest hoeveel extra random patches ermoeten worden gemaakt
        //    pt = strtok(NULL," ");
        //    sscanf(pt,"%d",&sys.nsitesN);
        } else if( strcmp(pt,"beta")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.beta);
        } else if( strcmp(pt,"sigma")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.sigma);
        } else if( strcmp(pt,"mobilityT")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityT);
        } else if( strcmp(pt,"mobilityR")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityR);
        } else if( strcmp(pt,"epsilonC")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonC);
        } else if( strcmp(pt,"epsilonP")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonP);     
        } else if( strcmp(pt,"epsilonN")==0) { //  regel toegevoegd om epsilon variabelen voor extra patches in te lezen
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonN);
        } else if( strcmp(pt,"boxl")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.boxl.x);
        } else if( strcmp(pt,"timestep")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&langevin.timestep);
        } else if( strcmp(pt,"ninter")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&langevin.ninter);
        } else if( strcmp(pt,"sim_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.sim_type);
        } else if( strcmp(pt,"start_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.start_type);
        } else if( strcmp(pt,"nshoot")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nshoot);
        } else if( strcmp(pt,"nrepswap")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nrepswap);
        } else if( strcmp(pt,"nstateswap")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nswapstates);
        } else if( strcmp(pt,"nreverse")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nreverse);
        } else if( strcmp(pt,"stateswapbias")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.stateswapbias);
        } else if( strcmp(pt,"fixedbias")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.fixedbias);
        } else if( strcmp(pt,"\n")==0) {
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Simulation\n")==0) {
            printf("Reading simulation parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Potential\n")==0) {
            printf("Reading potential parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"System\n")==0) {
            printf("Reading system parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Path\n")==0) {
            printf("Reading path parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"TIS\n")==0) {
            printf("Reading TIS parameters\n");
        } else {
            printf("Keyword unknown: %s\n",pt);
        }
    }

    fclose(fp);

    printf("Done reading path.inp\n");

    print_input();

    return;
}


void print_input(void) {

    printf("\nInput read from path.inp:\n");
    printf("Simulation\n");
    printf("sim_type        %d\n", sys.sim_type);
    printf("graphics        %d\n", sys.graphics);
    printf("ncycle1         %d\n", sys.ncycle1);
    printf("ncycle2         %d\n", sys.ncycle2);

    printf("\n");
    printf("Potential\n");
    printf("sigma           %lf\n", sys.sigma);
    printf("epsilonC        %lf\n", sys.epsilonC);
    printf("epsilonP        %lf\n", sys.epsilonP);
    printf("epsilonN        %lf\n", sys.epsilonN); //toegevoegd 10-9-15

    printf("\n");
    printf("System\n");
    printf("npart           %d\n", sys.npart);
    printf("nsites          %d\n", sys.nsites);
    //printf("nsitesN         %d\n", sys.nsitesN);
    printf("delta           %lf\n", sys.delta);
    printf("deltaN          %lf\n", sys.deltaN); //toegevoegd 10-9-15
    printf("beta            %lf\n", sys.beta);
    printf("mobilityT       %lf\n", sys.mobilityT);
    printf("mobilityR       %lf\n", sys.mobilityR);
    printf("boxl            %lf\n", sys.boxl.x);

    printf("\n");
    printf("Path\n");
    printf("timestep        %lf\n", langevin.timestep);
    printf("ninter          %d\n", langevin.ninter);

    printf("\n");
    printf("TIS\n");
    printf("nshoot          %d\n", path.nshoot);
    printf("nrepswap        %d\n", path.nrepswap);
    printf("nstateswap      %d\n", path.nswapstates);
    printf("nreverse        %d\n", path.nreverse);
    printf("stateswapbias   %d\n", path.stateswapbias);
    printf("fixedbias       %d\n", path.fixedbias);

    return;
}





void state_init() {

    //define the state volume boundaries
    int islice,ipart;
    Pts *psi;

    printf("\nDefining state defintions\n");

    path.initial_state=1;
    path.nstates=2;

    printf("Only two particles: placing them towards each other\n");
    for( ipart=0; ipart<sys.npart; ipart++ ) {
        psi=&slice[0].pts[ipart];
        psi->r=nulvec;
        psi->r.z -= 0.5 - sys.sigmaLJ*ipart*sys.sigma;
        //psi->q=RandomQuaternion();
        if(ipart==0) {
            psi->q.q0 = 0.0;
            psi->q.q1 = 0.0;
            psi->q.q2 = 0.0;
            psi->q.q3 = 1.0;
        }
        if(ipart==1) {
            psi->q.q0 = 0.0;
            psi->q.q1 = 1.0;
            psi->q.q2 = 0.0;
            psi->q.q3 = 0.0;
        }
    }



    //keep the state as small as possible
    state[0].target.energy = -2.0*sys.epsilonP;
    state[0].volume_op = state[0].lambda[0];
    printf("Definition of bound state:\n");
    printf("         Ground-state energy:     %lf\n",state[0].target.energy);
    printf("         Maximum energy boundary: %lf\n",state[0].target.energy+state[0].volume_op);
    printf("         State volume set at:     %lf\n",state[0].volume_op);
    printf("\n");

    state[1].target.energy = 0;
    state[1].mindist = 4*sys.sigma*sys.sigma;
    state[1].volume_op = state[1].lambda[0];
    printf("Definition of unbound state:\n");
    printf("         Ground-state energy:     %lf\n",state[1].target.energy);
    printf("         Maximum energy boundary: %lf\n",state[1].target.energy+state[1].volume_op);
    printf("         State volume set at:     %lf\n",state[1].volume_op);
    printf("\n");
    
    if(path.nstates>2) {
        state[2].volume_op = state[2].lambda[0];
        printf("Definition of non-specific state:\n");
        printf("         Ground-state energy not set yet\n");
        printf("         Maximum energy boundary not set yet\n");
        printf("         State volume set at:     %lf\n",state[2].volume_op);
        printf("\n");
    }


    path.current_gsen = state[path.initial_state-1].target.energy;
    create_all_rc(&(slice[0]));
    for(islice=1; islice<MAXSLICES; islice++) {
        slice[islice] = slice[0];
    }
    print_rc(&(slice[0]),path.initial_state);

    printf("\nDone defining state defintions\n");

    return;
}



void replica_init() {


    Replica *prep;
    int islice,istate,irep,jrep,pathlen,maxlength,len,type;

    //define interface volumes for every state


    printf("Defining replicas, interfaces and bootstrapping initial path\n");


    printf("Replica definition per state:\n");
    for(istate=0; istate<path.nstates; istate++) {
        printf("      State %d:\n",istate);
        path.nreplica = state[istate].nrep;
        printf("            Number of replicas %d\n",path.nreplica);
        for( irep=0; irep<path.nreplica; irep++) {
            state[istate].srep[irep].index = irep;
            state[istate].srep[irep].swapindex = irep;
            state[istate].srep[irep].lambda = state[istate].lambda[irep];
            printf("            %d %lf\n",irep,state[istate].srep[irep].lambda);
        }
    }

    path.nreplica = state[path.initial_state-1].nrep;
    for(irep=0; irep<path.nreplica; irep++) {
        replica[irep] = &state[path.initial_state-1].srep[irep];
    }


    path.scalefactor = 0.01;
    for(istate=0; istate<path.nstates; istate++) {
        state[istate].scalefactor = 0.01;
    }
    printf("Bootstrapping initial path for minus interface initial state\n");
    prep = replica[0];
    maxlength= MAXSLICES/2;
    len = trajectory_state_i(&slice[0],prep,maxlength,path.initial_state);

    printf("Length initial path from target to first interface: %d slices\n",len);

    for( islice=0 ;islice<=len; islice++ ) {
        slice[islice+len+1]=slice[islice];  
    }

    for( islice=0; islice<=len; islice++) {
        slice[len -islice] = slice[islice+len+1];
    }
  
    pathlen = 2*len+2;
    printf("Total length initial path minus interface: %d slices\n",pathlen);
    for( islice=0; islice<pathlen; islice++) {
        print_rc(&slice[islice],path.initial_state);
    }

    type =analyse_state_i(slice,prep,pathlen, path.initial_state);
    printf("Type initial path in minus interface: %d\n",type);

    replica[0]->type=type;
    replica[0]->pathlen= pathlen;
    path.nslices= pathlen;
    path.energy=slice[0].energy;
    path.current_replica=0;

    printf("Finished with replica initialization\n");

    return;
}



void read_lambda() {

    FILE *fplam;
    int irep;

    printf("Reading interface values from lambda_b.inp, lambda_u.inp, lambda_n.inp\n");


    if((fplam = fopen("lambda_b.inp","r"))==NULL){
        printf("Warning: lambda_b.inp not found\n");
        exit(1);
    }
    fscanf(fplam,"%d",&state[0].nrep);
    for( irep=0; irep<state[0].nrep; irep++) {
        fscanf(fplam,"%lf",&state[0].lambda[irep]);
    }
    fclose(fplam);


    if((fplam = fopen("lambda_u.inp","r"))==NULL){
        printf("Warning: lambda_u.inp not found\n");
        exit(1);
    }
    fscanf(fplam,"%d",&state[1].nrep);
    for( irep=0; irep<state[1].nrep; irep++) {
        fscanf(fplam,"%lf",&state[1].lambda[irep]);
    }
    fclose(fplam);


    if(path.nstates==3) {
        if((fplam = fopen("lambda_n.inp","r"))==NULL){
            printf("Warning: lambda_n.inp not found\n");
            exit(1);
        }
        fscanf(fplam,"%d",&state[2].nrep);
        for( irep=0; irep<state[2].nrep; irep++) {
            fscanf(fplam,"%lf",&state[2].lambda[irep]);
        }
        fclose(fplam);
    }

    printf("Done reading interface values\n");

    return;
}







void print_tis_setup() {

    int istate, totalrep,type;

    printf("TIS setup:\n");

    totalrep = 0;
    printf("Total number of states defined: %d\n", path.nstates);
    printf("Total number of replicas defined: ");
    for(istate=0; istate<path.nstates; istate++) {
        printf("%d + ",state[istate].nrep);
        totalrep+=state[istate].nrep;
    }
    printf("= %d\n",totalrep);

    printf("Current initial state:                         %d\n",path.initial_state);
    printf("Current initial state according to in_state(): %d\n", in_state(&slice[0]));
    printf("Current middle state according to in_state():  %d\n", in_state(&slice[(int)((path.nslices-1)/2.0)]));
    printf("Current final state according to in_state():   %d\n", in_state(&slice[path.nslices-1]));
    printf("Current replica: %d\n", path.current_replica);
    if(path.current_replica==0) {
        type =analyse_state_i(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    else {
        type =analyse(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    printf("Type of current path: %d\n", type);
    printf("Length current path: %d\n", path.nslices);


    printf("\n");

    return;
}





void stats_init() {

    int istate,irep,iaver,iacc;

    printf("Initializing stats\n");
    for(istate=0; istate<MAXSTATES; istate++) {
        for(irep=0; irep<MAXREPLICA; irep++) {
            for(iaver=0; iaver<NSTAT; iaver++) {
                path.block_stats[istate][irep].aver[iaver].now=0;
                path.block_stats[istate][irep].aver[iaver].sum=0;
                path.block_stats[istate][irep].aver[iaver].sumsq=0;
                path.block_stats[istate][irep].aver[iaver].n=0;
                path.final_stats[istate][irep].aver[iaver].now=0;
                path.final_stats[istate][irep].aver[iaver].sum=0;
                path.final_stats[istate][irep].aver[iaver].sumsq=0;
                path.final_stats[istate][irep].aver[iaver].n=0;
            }
        }
    }

    for(istate=0; istate<MAXSTATES; istate++) {
        for(irep=0; irep<MAXREPLICA; irep++) {
            for(iacc=0; iacc<NACC; iacc++) {
                path.block_stats[istate][irep].mcacc[iacc].acc=0;
                path.block_stats[istate][irep].mcacc[iacc].tries=0;
                path.block_stats[istate][irep].mcacc[iacc].ratio=0;
                path.final_stats[istate][irep].mcacc[iacc].acc=0;
                path.final_stats[istate][irep].mcacc[iacc].tries=0;
                path.final_stats[istate][irep].mcacc[iacc].ratio=0;
            }
        }
    }

    for(istate=0; istate<path.nstates; istate++) {
        for(irep=0; irep<state[istate].nrep; irep++) {
            sprintf(path.block_stats[istate][irep].mcacc[0].name,"shots replica %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[1].name,"swaps replica %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[2].name,"swap  states  %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[3].name,"revs  replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[0].name,"shots replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[1].name,"swaps replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[2].name,"swap  states  %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[3].name,"revs  replica %3d ",irep);
        }
    }


    printf("Done initializing stats\n");

    return;
}







void traj_input() {

    char dum[40];
    char filename[100];
    FILE *fp;
    int i,j,type;
	Slice *psl;

	printf("\nIn traj_input()\n");
	printf("Reading in previous trajectory\n");
 
    sprintf(filename,"trajectory.inp");
    if ((fp = fopen(filename,"r"))==NULL){
        printf("input: can't open %s\n",filename);
        return;
    }
    else {
        fscanf(fp,"%d %s",&path.nslices,dum); 
        fscanf(fp,"%d %s",&path.current_replica,dum); 
        fscanf(fp,"%d %s",&path.initial_state,dum); 
        fscanf(fp,"%lf %lf %lf  %s",&sys.boxl.x,&sys.boxl.y,&sys.boxl.z,dum);
        for(j=0;j<path.nslices;j++) {
            for(i=0;i<sys.npart;i++){
                fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
                  &(slice[j].pts[i].r.x),&(slice[j].pts[i].r.y),&(slice[j].pts[i].r.z),
                  &(slice[j].pts[i].q.q0),&(slice[j].pts[i].q.q1),&(slice[j].pts[i].q.q2),&slice[j].pts[i].q.q3);
            }
        }
        fclose(fp);
    }  


    //retrieve all additional TIS parameters from read input
    path.nreplica = state[path.initial_state-1].nrep;
	for(i=0; i<path.nreplica; i++) {
		replica[i] = &state[path.initial_state-1].srep[i];
	}


    if(path.initial_state == 1 || path.initial_state==2) {
	    path.current_gsen = state[path.initial_state-1].target.energy;
    }

    if(path.initial_state == 3) {
        potential_energy(&slice[0]);
        path.current_gsen = -2. * sqrt(sys.site[slice[0].minisite].eps *sys.site[slice[0].minjsite].eps);
    }

    for(i=0;i<path.nslices;i++) {
        psl = &(slice[i]);
        create_all_rc(psl);
        printf("Slice %d  op %lf state %d energyT %lf rijdist %lf energyN %lf energy %lf inwindow %d\n", 
         i, 
         psl->order_parameter, 
         in_state(psl),
         psl->energyT,
         psl->rijdist,
         psl->energyN,
         psl->energy,
         in_upper_window(psl,replica[path.current_replica],path.initial_state));
    }
    path.final_state = in_state(&slice[path.nslices-1]);

    if(path.current_replica==0) {
        type =analyse_state_i(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    else {
        type =analyse(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    replica[path.current_replica]->pathlen=path.nslices;
    replica[path.current_replica]->type=type;


    printf("Warm start information:\n");
    printf("       number of slices:       %d\n",path.nslices);
    printf("       initial state:          %d\n",in_state(&slice[0]));
    printf("       middle state:           %d\n",in_state(&slice[(int)(path.nslices/2.0)]));
    printf("       final state:            %d\n",in_state(&slice[path.nslices-1]));
    printf("       current initial state:  %d\n",path.initial_state);
    printf("       current replica:        %d\n",path.current_replica);
    printf("       current gsen:           %lf\n",path.current_gsen);
    printf("       type of path:           %d\n",type);
        
	printf("Done with traj_input\n\n");

    return;
}



void traj_output() {

    FILE *fp;
    char filename[240];
    int i,j; 
    Replica *prep;

    prep = replica[path.current_replica];
  
    sprintf(filename,"trajectory.out");
    if ((fp = fopen(filename,"w"))==NULL) {
        printf("output: can't open %s\n",filename);
        return;
    }
    else {
        fprintf(fp,"%d slices\n",prep->pathlen); 
        fprintf(fp,"%d current_replica\n",path.current_replica); 
        fprintf(fp,"%d initial_state\n",path.initial_state); 
        fprintf(fp,"%.16lf %.16lf %.16lf boxl,dr\n",sys.boxl.x,sys.boxl.y,sys.boxl.z);
        for(j=0;j<path.nslices;j++) {
            for(i=0;i<sys.npart;i++){
                fprintf(fp,"%12.16lf %12.16lf %12.16lf %12.16lf %12.16lf %12.16lf %12.16lf\n",
                  slice[j].pts[i].r.x,slice[j].pts[i].r.y,slice[j].pts[i].r.z,
                  slice[j].pts[i].q.q0,slice[j].pts[i].q.q1,slice[j].pts[i].q.q2,slice[j].pts[i].q.q3);
            }
        }
        fclose(fp);
    }
	
    return;
}






void dos_input() {

    FILE *fp;
    char filename[30];
    int i,j,idum; 

    //Wang-Landau Dos
    sprintf(filename,"dos_all.inp");
    if ((fp = fopen(filename,"r"))==NULL){
        printf("output: can't open %s\n",filename);
        return;
    } 
    else {
        printf("reading dos data\n");
        for(i=0;i<path.nstates;i++) {
            fscanf(fp,"%lf",&state[i].scalefactor);
            for(j=0;j<state[i].nrep;j++) {
                fscanf(fp,"%d %lf\n",&idum,&state[i].srep[j].dos);
                printf("%d %lf\n",j,state[i].srep[j].dos);
            }
            fscanf(fp,"\n");
        }
        fclose(fp);
    }



	return;
}





void dos_output() {

	FILE *fp;
    int i,j;
	
    if ((fp = fopen("dos_all.dat","w"))==NULL){
        printf("output: can't open file.dat\n");
        return;
    }
    else {
        for(i=0;i<path.nstates;i++) {
            fprintf(fp,"%lf\n",state[i].scalefactor);
            for(j=0;j<state[i].nrep;j++) {
                fprintf(fp,"%d %g\n",j,state[i].srep[j].dos);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }


    if ((fp = fopen("dos_lambda.dat","w"))==NULL){
        printf("output: can't open file.dat\n");
        return;
    }
    else {
        for(i=0;i<path.nstates;i++) {
            for(j=1;j<state[i].nrep;j++) {
                fprintf(fp,"%lf %g\n",state[i].srep[j].lambda,state[i].srep[j].dos-state[i].srep[1].dos);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

	
    return;
}





void warm_start_gen() {

	traj_output();
	dos_output();

	return;
}



void conf_output(Slice *psl) {

    int ipart, isite;
    FILE *fp;


    if((fp=fopen("conf.out","w"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fprintf(fp,"%d %d %lf\n", sys.npart, sys.nsites, sys.boxl.x);
        for(ipart=0; ipart<sys.npart; ipart++) {
            fprintf(fp,"%lf %lf %lf ", psl->pts[ipart].r.x, psl->pts[ipart].r.y, psl->pts[ipart].r.z);
            fprintf(fp,"%lf %lf %lf %lf\n", psl->pts[ipart].q.q0, psl->pts[ipart].q.q1, psl->pts[ipart].q.q2, psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    return;
}




void conf_input(Slice *psl) {

    int ipart, isite,npart,nsites;
    double boxl;
    FILE *fp;


    if((fp=fopen("conf.inp","r"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fscanf(fp,"%d %d %lf\n", &npart, &nsites, &boxl);
        for(ipart=0; ipart<npart; ipart++) {
            fscanf(fp,"%lf %lf %lf ", &psl->pts[ipart].r.x, &psl->pts[ipart].r.y, &psl->pts[ipart].r.z);
            fscanf(fp,"%lf %lf %lf %lf\n", &psl->pts[ipart].q.q0, &psl->pts[ipart].q.q1, &psl->pts[ipart].q.q2, &psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    if(npart!=sys.npart) {
        printf("Warning: number of particles in system not same as in conf.inp\n");
    }
    if(nsites!=sys.nsites) {
        printf("Warning: number of sites in system not same as in conf.inp\n");
    }
    if(boxl!=sys.boxl.x) {
        printf("Warning: boxlength in system not same as in conf.inp\n");
    }


    return;
}
