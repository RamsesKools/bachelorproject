#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"


System sys;
State state[MAXSTATES];
Replica *replica[MAXREPLICA];
vector patches[NPART][NSITES];
Langevin langevin;
Slice *slice,*trial,*fwtrial,*bwtrial;
Path path;
Stats nulstat={0};
vector nulvec={0};


int main(int argc, char **argv) {

    int icycle,jcycle;
    double energy;

    //change this to setup simulation
    //  TIS simulation
    //      definition states
    //      order parameter
    //      interfaces
    setup_simulation();

    if(sys.graphics==1) {
        Init_Graphics(argc,argv);
    }

    for(icycle=0; icycle<sys.ncycle1; icycle++) {
        for(jcycle=0; jcycle<sys.ncycle2; jcycle++) {
            if(sys.sim_type==0) {
                bmdcycle();
            }
            else if(sys.sim_type==1) {
                tiscycle();
            }
        }
        printf("\nBLOCK %d\n",icycle);
        terminate_block();
    }

    finalstat();

    return 0;
}



void tiscycle() {

    int i,j,iwhich,irep,jrep,drep;
    double en;

    for(i=0;i<(path.nshoot+path.nrepswap+path.nreverse+path.nswapstates);i++) {  
        iwhich = (int)(RandomNumber()*(path.nshoot+path.nrepswap+path.nreverse+path.nswapstates));
        if (iwhich<path.nshoot) {
            j= path.current_replica;
            shoot_oneway(replica[j]);
        } 
        else if(iwhich<(path.nshoot + path.nrepswap)) {  
            irep=  path.current_replica;
            drep =2*((int)(2*RandomNumber()))-1;
            if (abs(drep)!= 1) {
                printf("error:drep = %d\n",drep);
            }
            jrep= irep+drep;
            swap_replica(irep,jrep);
        } 
        else if(iwhich<(path.nshoot+path.nrepswap+path.nswapstates)) {
            irep = path.current_replica;
            swap_states(irep,jrep);
        } 
        else if(iwhich < (path.nshoot + path.nrepswap + path.nswapstates + path.nreverse)) {
            irep=path.current_replica; 
            reverse_replica(irep);
        }
        crossinghistogram(replica[path.current_replica]);
        get_flux(replica[path.current_replica]);
    }
    printstatustis();

    return;
}


void printstatustis() {

    int istate,jstate;
  
    istate=in_state(&slice[0]);
    jstate=in_state(&slice[path.nslices-1]);

    printf("status: state= %d replica=%d path-type=%d path-length %d initial %d, final %d\n",
        path.initial_state,path.current_replica,replica[path.current_replica]->type, path.nslices,istate,jstate);

    return;
}



void bmdcycle() {

    int istep;

    propagate_bd(&slice[0]);
    if(sys.npart==1) {
        sampleangles(0); 
        samplediffusion(0);
    }

    return;
}


void mainloop_for_graphics() {

    static int icycle,jcycle,initial;

    if(initial==0) {
        initial=1;
        icycle=jcycle=1;
    }

    if(sys.sim_type==0) {
        bmdcycle();
    }

    else if(sys.sim_type==1) {
        tiscycle();
    }

    if(jcycle%sys.ncycle2==0) {
        printf("\nBLOCK %d\n",icycle);
        terminate_block();
        jcycle=0;
        if(icycle%sys.ncycle1==0) {
            printf("Finished Simulation\n");
            finalstat();
            exit(1);
        }
        icycle++;
    }

    jcycle++;

    return;
}


void printstatusbmd() {

    double energy,l;

    energy = potential_energy(&slice[0]);
    printf("Energy system %lf\n",energy);

    return;
}



void terminate_block() {

    int i=0,j=0,k=0,istate=0,irep=0,iaver=0;
    
    if(sys.sim_type==1) {
        printf("Updating block averages\n");
        for(istate=0; istate<path.nstates; istate++) {
            for(irep=0; irep<state[istate].nrep; irep++) {
                for (iaver=0; iaver<NSTAT; iaver++) {
                    update_blockaver(path.block_stats[istate][irep].aver[iaver]);
                }
            }
        }
        print_pathacceptances();
        print_pathnumbers();
    }
   
    else if(sys.sim_type==0) {
        printstatusbmd();
        conf_output(&slice[0]);
    }


    if(sys.sim_type==1) {
        printf("Updating final averages\n");
        for(istate=0; istate<path.nstates; istate++) {
            for(irep=0; irep<state[istate].nrep; irep++) {
                for(iaver=0; iaver<NSTAT; iaver++) {
                    update_finaver(path.final_stats[istate][irep].aver[iaver],path.block_stats[istate][irep].aver[iaver]);
                }
            }
        }

 
 		for(istate=0; istate<path.nstates; istate++) {
			for(irep=0; irep<state[istate].nrep; irep++) {
			    for(k=0; k<NACC; k++) {
                    update_blockmc(path.final_stats[istate][irep].mcacc[k],path.block_stats[istate][irep].mcacc[k]);
				}
			}
		}
        update_wanglandau();
        warm_start_gen();
	}

	return;
}



void finalstat() {

    int i,istate,irep;

    if(sys.sim_type==1) {
        for(istate=0; istate<path.nstates; istate++) {
            for(irep=0; irep<state[istate].nrep; irep++) {
                for (i=0; i<NSTAT; i++) {
                    update_blockaver(path.final_stats[istate][irep].aver[i]);
                }
            }
        }
    }

    printf("\nFinished Simulation\n");
    print_finpathacceptances();
    print_pathnumbers();

    if(sys.sim_type==1) {
        print_tisstats();
        reweight_crosshist();
        reweight_pathtypenumbers();
    }

    if(sys.sim_type==0) {
        if(sys.npart==1) {
            sampleangles(1); 
            samplediffusion(1);
        }
    }


    return;
}



void print_pathacceptances() {

    int i,j,k,nrep,ntot[MAXSTATES];
    int istate,jstate,irep,jrep;
    double f0,f1,flux[MAXSTATES];

    if(sys.sim_type==1) {
        for(istate=0; istate<path.nstates; istate++) {
            printf("\nPath moves Acceptances state %d\n",istate);
            for(k=0; k<NACC; k++) {
                for(irep=0; irep<state[istate].nrep; irep++) {
                    if(path.block_stats[istate][irep].mcacc[k].tries>0) {
                        path.block_stats[istate][irep].mcacc[k].ratio = 
                            (double)path.block_stats[istate][irep].mcacc[k].acc/(double)path.block_stats[istate][irep].mcacc[k].tries;
                    }
                    printf("%7d accepted %sout of %7d trial moves, ratio = %g\n",
                            path.block_stats[istate][irep].mcacc[k].acc,
                            path.block_stats[istate][irep].mcacc[k].name,
                            path.block_stats[istate][irep].mcacc[k].tries,
                            path.block_stats[istate][irep].mcacc[k].ratio);
                    path.block_stats[istate][irep].mcacc[k].ratio=0;
                }
            }
        }

  	    printf("\nSlice number averages\n");
        for (k=0;k<MAXREPLICA;k++) { 
            printf("     Path length replica %2d ",k); 
            for (i=0;i<path.nstates;i++) { 
                printf("%6.1lf ",(double) path.block_stats[i][k].aver[1].sum);
            } 
            printf("\n");
        }
    }
}


void print_finpathacceptances() {

    int i,j,k,nrep,ntot[MAXSTATES];
    int istate,jstate,irep,jrep;
    double f0,f1,flux[MAXSTATES];

    if(sys.sim_type==1) {
        for(istate=0; istate<path.nstates; istate++) {
            printf("\nPath moves Acceptances state %d\n",istate);
            for(k=0; k<NACC; k++) {
                for(irep=0; irep<state[istate].nrep; irep++) {
                    if(path.final_stats[istate][irep].mcacc[k].tries>0) {
                        path.final_stats[istate][irep].mcacc[k].ratio = 
                            (double)path.final_stats[istate][irep].mcacc[k].acc/(double)path.final_stats[istate][irep].mcacc[k].tries;
                    }
                    printf("%7d accepted %sout of %7d trial moves, ratio = %g\n",
                            path.final_stats[istate][irep].mcacc[k].acc,
                            path.final_stats[istate][irep].mcacc[k].name,
                            path.final_stats[istate][irep].mcacc[k].tries,
                            path.final_stats[istate][irep].mcacc[k].ratio);
                    path.final_stats[istate][irep].mcacc[k].ratio=0;
                }
            }
        }

  	    printf("\nSlice number averages\n");
        for (k=0;k<MAXREPLICA;k++) { 
            printf("     Path length replica %2d ",k); 
            for (i=0;i<path.nstates;i++) { 
                printf("%6.1lf ",(double) path.final_stats[i][k].aver[1].sum);
            } 
            printf("\n");
        }
    }
}


void print_pathnumbers() {

    int i,j,k,nrep,ntot[MAXSTATES];
    int istate,jstate,irep,jrep;
    double f0,f1,flux[MAXSTATES];

    if(sys.sim_type==1) {
        printf("\nPathlength averages        ");
        for (i=1;i<=path.nstates;i++) { 
            printf("%6d ",i);
        }
        printf("\n");

        for (k=0;k<MAXREPLICA;k++) { 
            printf("     Path length replica %2d ",k); 
            for (i=0;i<path.nstates;i++) { 
                printf("%6.1lf ",(double) state[i].srep[k].avlen/(double)state[i].srep[k].navlen);
            }
            printf("\n");
        }

        printf("\nPathtype averages        ");
        for (i=1;i<=path.nstates;i++) { 
            printf("  %d->%d  %d->j",i,i,i);
        }
        printf("\n");

        for (k=0;k<MAXREPLICA;k++) { 
            printf("     Path type replica %2d ",k); 
            for (i=0;i<path.nstates;i++) { 
                for (j=0;j<2;j++) { 
                    printf("%5d ", state[i].type_mat[k][j]);
                }
            }
            printf("\n");
        }


        printf("\nMSTIS matrix  ");
        for (i=0;i<path.nstates;i++) { 
            printf("    %2d",i+1);
        }
        printf("\n");
        for (i=0;i<path.nstates;i++) {
            ntot[i]=0;
        }
        for (k=0;k<path.nstates;k++) { 
            printf("     state  %2d ",k+1); 
            for (i=0;i<path.nstates;i++) { 
                printf("%5d ", state[i].mstis_mat[k]);
                ntot[i]+=state[i].mstis_mat[k];
            }
            printf("\n");   
        }
  	    printf("     total     ");
        for (i=0;i<path.nstates;i++) printf("%5d ",  ntot[i]);
        printf("\n"); 

        printf("\nNormalized matrix   ");
        for (i=0;i<path.nstates;i++) { 
            printf("%2d      ",i+1);
        }
        printf("\n");
        for (k=0;k<path.nstates;k++) { 
            printf("     state  %2d ",k+1); 
            for (i=0;i<path.nstates;i++) { 
                printf("%6.5lf ", (double)state[i].mstis_mat[k]/ ntot[i]);
            }
            printf("\n");
        }

  	    
        printf("\nFluxes     0        1     total     flux\n");
        for (i=0;i<path.nstates;i++) { 
            f0 = (double) state[i].flux0/ state[i].nflux0;
            f1 = (double) state[i].flux1/ state[i].nflux1;
            flux[i]=  1./(f0+f1);
            printf("%2d      ",i+1);
            printf("%7.2lf %7.2lf %7.2lf %7.6lf \n", f0, f1, (f0+f1), 1./(f0+f1));
        }
        printf("\n");

        printf("\nCrossing prob from dos at last replica for state 1 and state 2: %d and %d\n",state[0].nrep-1,state[1].nrep-1);
        for (i=0;i<path.nstates;i++) { 
            nrep = state[i].nrep;
            printf("%2d      ",i+1);
            printf("%7.2lf %7.2lf %12.5g \n", 
              state[i].srep[nrep-1].dos, state[i].srep[nrep-1].dos-state[i].srep[1].dos, exp(state[i].srep[nrep-1].dos-state[i].srep[1].dos));
        }
        printf("\n");

        printf("\nRate matrix   ");
        for (i=0;i<path.nstates;i++) { 
            printf("%2d      ",i+1);
        }
        printf("\n");
        for (k=0;k<path.nstates;k++) { 
            printf("     state  %2d ",k+1); 
            for (i=0;i<path.nstates;i++) { 
                nrep = state[i].nrep;
                state[i].rate[k] =  exp(state[i].srep[nrep-1].dos-state[i].srep[1].dos)*flux[i]*state[i].mstis_mat[k]/ ntot[i];
                printf("%12.5g ", (double)state[i].rate[k]);
            }  
            printf("\n");
        }

        printf("\n\n");

	}


    fflush(NULL);
 
    return;
}




