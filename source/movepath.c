#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "path.h"


void create_all_rc(Slice *psl) {


    //bound energies (bound state), non-specific energies (non-specific state), or system energy (unbound state) are calculated in potential_energy()

    psl->energy = potential_energy(psl);
    //for target state: use only energies for target patches
    if(path.initial_state==1) {
        psl->order_parameter = psl->energyT - path.current_gsen;
    }

    //for unbound state: use entire potential energy
    else if(path.initial_state==2) {
        psl->order_parameter = -psl->energy;
    }

    //for non-specific state: use only energies for non-specific patches
    else if(path.initial_state==3) {
        if(psl->rijdist >= sys.sigmaLJ) {
            psl->order_parameter = psl->energyN - path.current_gsen;
        }
    }

    return;
}
    




double get_rc(Slice *psl) {

    return psl->order_parameter;

}
  



double print_rc(Slice *psl,int state_index) {

    double op;
 
    printf("op = %g\n",(op=get_rc(psl)));

    return op;
}




int in_state(Slice *psl) {

    //define state boundaries
    //return 1 if in state 1, the bound state
    if ((psl->energyT - state[0].target.energy) < state[0].volume_op) {
        return 1;
    }

    //return 2 if in state 2, the unbound state
    if (psl->rijdist > state[1].mindist ) {
        return 2;
    }

    if(path.nstates>2) {
        //return 3 if in state 3, the non-specific state
        if ((psl->energyN + sqrt(sys.site[psl->minisite].eps*sys.site[psl->minjsite].eps)) < state[2].volume_op) {
            return 3;
        }
    }

    return 0;
}


 

int trajectory(Slice *psl, int maxlength) {

    int i,istate;

    for (i=1;i<maxlength;i++) {
        psl[i] = psl[i-1]; 
        propagate_bd(&psl[i]);
        create_all_rc(&(psl[i]));
        if (in_state(&psl[i])>0) {
            return i;
        }
    }

    printf("maxlength %d reached ",maxlength);
    printf("psl %d energy %lf order_parameter %lf mindist %lf\n",i,psl[i-1].energy,psl[i-1].order_parameter,psl[i-1].rijdist);

    return 0;
}



int in_upper_window(Slice *psl, Replica *prep, int state_index) {

    //printf("in upper window\n");
    //printf("psl op %lf\n",psl->order_parameter);
    //printf("prep lambda %lf\n",prep->lambda);
    return (psl->order_parameter > prep->lambda);
}


int trajectory_state_i(Slice *psl, Replica *prep, int maxlength, int state_index) {

    int i,istate;
    static int count,tries;

    for (i=1;i<maxlength;i++) {
        psl[i] = psl[i-1]; 
        propagate_bd(&psl[i]);
        create_all_rc(&(psl[i]));
        if (in_upper_window(&psl[i], prep,state_index)) {
            return i;
        }
    }

    return 0;
}



int analyse ( Slice *psl, Replica *prep, int len, int state_index) {

    int i,initial_state,final_state;
 
    initial_state=in_state(&psl[0]);

    final_state=in_state(&psl[len-1]);

    if (initial_state==state_index) {
        if (final_state==state_index) {
            for (i=0;i<len;i++) {
                if (in_upper_window(&psl[i], prep,state_index)){
                    return 1;
                }
            }
        //printf("path did not reach lambda %g\n",prep->lambda);
        //print_rc(&psl[0],state_index);
        //print_rc(&psl[len-1],state_index);
            return 0;
        } 
        else {
            if (final_state==0) {
                printf("final state unknown\n");
                print_rc(&psl[0],state_index);
                print_rc(&psl[len-1],state_index);
                return 0;
            }
            for (i=0;i<len;i++) {
                if (in_upper_window(&psl[i], prep,state_index)) {
                    //printf("intermediate slice %d is larger than lambda %g \n",i,prep->lambda);
                    return 2;
                }
            }
            //printf("warning: type 2 but path did not reach lambda\n");
            //dprint(initial_state);
            //dprint(final_state);
            //dprint(path.current_replica);
            //for(i=0; i<prep->pathlen; i++) {
            //  print_rc(&(trial[i]),state_index);  
            //  }
            return 0;
        }
    }

    printf("analyse: path of length %d corrupted in replica %d, initial %d  final %d, state index %d \n",len, prep->index,initial_state,final_state,state_index);

    print_rc(&psl[0],state_index);
    print_rc(&psl[len-1],state_index);

    return 0; 
}
  



int analyse_state_i(Slice *psl, Replica *prep, int len, int state_index) {

    int i,visitedA,visited2;
 
    visitedA=visited2 =0;
    if (in_upper_window(&psl[0], prep,state_index)){
        for (i=1;i<len;i++) {
            if (in_state(&psl[i])==state_index){
                visitedA =1;
            }
            if ((visitedA) && (in_upper_window(&psl[i], prep,state_index))) {
                visited2 =1;
            }
        }
        if ((visitedA) && (in_upper_window(&psl[len-1], prep,state_index))){
            //printf("path OK\n");
            return 1;
        }
        if ((visitedA) && (visited2) ) {
            printf ("path starts in 1, visits A and then 1 but doesn't end in 1\n");
            return 1;
        }
    }

    for(i=0;i<prep->pathlen;i++) {
        print_rc(&slice[i],state_index);
        printf(" lambda = %g ",prep->lambda );
        dprint(in_state(&slice[i]));
    }     
    printf("analyse i: path corrupted in replica %d   \n",prep->index);
    dprint(path.initial_state);

    return 0; 
}




int shoot_oneway(Replica *prep)
{
    int i,j,index,len,pathlen,maxlength,trial_pathlen,type;
    int initial_state,final_state;
   

    if (prep->index==0) {
        return 0;
    }
  
    //do not want to keep minus interface statistics for shooting
    path.block_stats[path.initial_state-1][prep->index].mcacc[0].tries++;

    index =0;
    for (i=0;i<path.nslices;i++) {
        if (in_upper_window(&slice[i], prep, path.initial_state)) {
            index = i;
            i = path.nslices;
        }
    }
  
    //sanity check
    if (in_upper_window(&slice[index], prep, path.initial_state)==0) {
        return 0;
    }
  
    for(i=0;i<=index;i++) {
        trial[i]=slice[i];
    }

    maxlength= MAXSLICES-index;
    len = trajectory(&trial[index],maxlength);
    if (len==0) {
        printf("Shoot oneway: trajectory too long. State %d\n",path.initial_state);
        return 0;
    }
    
    trial_pathlen = index+len+1;
    if(trial_pathlen<3) {
        printf("path is too short: %d, state %d\n",trial_pathlen,path.initial_state);
        for(i=0; i<path.nslices; i++) {
            print_rc(&slice[i],path.initial_state);
        }
        return 0;
    }

    type = analyse(trial,prep,trial_pathlen, path.initial_state);
    if ( type==0 ) {
        printf("Shoot Oneway: wrong type after shooting.\n");
        dprint(path.initial_state);
        initial_state = in_state(&(trial[0]));
        gprint(trial[0].energy);
        dprint(initial_state);
        dprint(final_state);
        dprint(prep->index);
        return 0;
    }


    for(i=index+1;i<trial_pathlen;i++) {
        slice[i] = trial[i];
    }
  
    final_state=in_state(&slice[trial_pathlen-1]);
    path.final_state = final_state;

    prep->pathlen= trial_pathlen;
    path.nslices = trial_pathlen;
    prep->type =type;

    path.block_stats[path.initial_state-1][prep->index].mcacc[0].acc++;

    return 1;
}




int swap_replica(int irep, int jrep) {

    Replica *prepi,*prepj;
    Slice *swaptrial, h;
    int i,pathlen,ABexchOK,ireptype,jreptype,type;
    //double  aux;

    prepi = replica[irep];
    prepj = replica[jrep];

    path.block_stats[path.initial_state-1][prepi->index].mcacc[1].tries++;

    if(path.fixedbias==0) {
        prepi->dos += state[path.initial_state-1].scalefactor;
    }
  
    if (jrep<0) {
        return 0;
    }

    if ( (irep==0) || (jrep==0) ) {
        return swap_replica_0(irep, jrep);
    }

    if (( irep==path.nreplica-1 ) && ( jrep==path.nreplica )) {
        return 0;
    }
    
    if ( jrep>path.nreplica-1 ) {
        return 0;
    }

    //aux = prepi->dos - prepj->dos;
 
    if ( RandomNumber() > exp(prepi->dos - prepj->dos) ) {
        return 0; 
    }
 
    type =analyse(slice, prepj, path.nslices, path.initial_state);
    if (type==0) {
        return 0;
    }

    prepj->pathlen =prepi->pathlen;
    prepj->type =  prepi->type;
    SWAP(replica[irep]->swapindex,replica[jrep]->swapindex,i);
    path.current_replica=jrep;

    path.block_stats[path.initial_state-1][prepi->index].mcacc[1].acc++;

    return 1;
}






int swap_replica_0(int irep, int jrep) {

    Replica *prep;
    Slice *trial0;
    Slice *trial1;
    int i,j,pathlen0,pathlen1,inA,in1,index,maxlength,start,len,type,type0,type1,reversal;
    int final_state;

    trial0=fwtrial;
    trial1=bwtrial;

    if (irep ==0) {
        prep =replica[0];
        type =analyse_state_i(slice,prep,path.nslices, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 0->1 error : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("error: replica 0 has wrong type\n");
            return 0;
        }
        reversal = (RandomNumber()<0.5) ? 1 : 0;
        
        if (reversal) {
            pathlen0 = path.nslices;
            for(j=0;j<pathlen0;j++) {  
                trial1[j]= slice[pathlen0 - 1 -  j];// use trial1 for temp storage, to not affect slice
            }
        } 
        else {
            for(i=0;i< path.nslices;i++) {
                trial1[i]= slice[i];// use trial1 for temp storage, to not affect slice
            }
        }

        for(i=path.nslices-1;i>=0;i--) {
            if (in_state(&trial1[i])) {
                start = i;
                break;
            }
        }

        for(i=start;i<path.nslices;i++) {
            trial0[i-start]=trial1[i];
        }

        index = path.nslices - start - 1;
        maxlength = MAXSLICES-index;
        len = trajectory(&trial0[index],maxlength);
        if (len==0) {
            printf(" Swap replica 0: trajectory too long\n");
            return 0;
        }

        pathlen0 = index +len+1;         
        if (pathlen0> MAXSLICES) {
            printf("error: pathlen = %d\n",pathlen0);
            return 0;
        }

        type= analyse(trial0,replica[1],pathlen0,path.initial_state );
        if (type==0) {
            return 0;
        }

        type0=type;
        for(i=0;i<pathlen0;i++) {
            slice[i]=trial0[i];
        }

        replica[1]->pathlen= pathlen0;
        replica[1]->type = type0;
        path.current_replica=jrep;
        path.nslices = pathlen0;
        prep = replica[1];
        type = analyse(slice,prep,prep->pathlen, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 0->1 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("swap 0->1: replica 0 has wrong type\n");
            return 0;
        }

        final_state=in_state(&slice[path.nslices-1]);
        path.final_state = final_state;
    } 


    if (jrep==0) {
        prep =replica[1];
        type =analyse(slice, prep, path.nslices, path.initial_state);  
        if (prep->type!=type) {
            printf("error 1->0 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }

        if(type==0) {
            //printf("reject: replica 1 has wrong type\n");
            return 0;
        }
    
        reversal = (RandomNumber()<0.5) ? 1 : 0;
        if ((type == 2) && (reversal ==0)) {
            reversal = 1;
            //printf("type is 2, and reversal =0, still accept\n");
            //return 0;
        }

        if (reversal) {
            pathlen1 =path.nslices;
            for(j=0;j< pathlen1;j++) { 
                trial0[j]= slice[pathlen1 - 1 -  j]; // use trial0 for temp storage
            }
        }
        else {
            for(i=0; i<path.nslices; i++) {
                trial0[i]= slice[i];// use trial0 for temp storage, to not affect slice
            }
        }
  
        for(i=path.nslices-1;i>=0;i--) {
            if (in_upper_window(&trial0[i],replica[0],path.initial_state)) {
                start = i;
                break;
            }
        }
        
        for(i=start;i<path.nslices;i++) {
            trial1[i-start]=trial0[i];
        }

        index= path.nslices - start-1;
        maxlength= MAXSLICES-index;
        len = trajectory_state_i(&trial1[index],prep,maxlength,path.initial_state);
        if (len==0) {
            printf(" Swap replica 0: trajectory too long\n");
            return 0;
        }

        pathlen1 = index +len+1;
        if (pathlen1 > MAXSLICES) {
            //printf("error: pathlen = %d\n",pathlen1);
            return 0;
        }
        
        type =analyse_state_i(trial1,replica[0],pathlen1, path.initial_state);
        if ( type==0) {
            printf(" Swap replica 0: trajectory too long\n");
            return 0;
        } 
        type1=type;

        for(i=0;i<pathlen1;i++) {
            slice[i]=trial1[i];
        }
        replica[0]->pathlen= pathlen1;
        replica[0]->type = type1;
        path.current_replica=jrep;
        path.nslices=pathlen1;
        prep = replica[0];
    
        type =analyse_state_i(slice,prep,prep->pathlen, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 1->0 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("swap 1->0: replica 0 has wrong type\n");
            return 0;
        }
    }
  
    SWAP(replica[irep]->swapindex,replica[jrep]->swapindex,i);

    path.block_stats[path.initial_state-1][prep->index].mcacc[1].acc++;

    return 1;
}






int swap_states(int irep,int jrep) {

    Replica *prepi,*prepj;
    int i,j,initial_state,final_state,pathlen,type,start,index,maxlength,len;
    double aux;

    if(path.current_replica!=(state[path.initial_state-1].nrep-1)) {
        return 0;
    }
    
    if(path.stateswapbias==2) {
        state[path.initial_state-1].dos+=path.scalefactor;
    }
    
    prepi=replica[irep];

    if (prepi->type!=2) {
        return 0;
    }

    pathlen =prepi->pathlen;
 
    initial_state=in_state(&slice[0]);
    final_state=in_state(&slice[pathlen-1]);

    if(final_state==initial_state) {
        printf("Warning swap states: final state is not initial state despite being type 2\n");
        return 0;
    }

    if (initial_state != path.initial_state ) {
        printf("error: initial state not correct\n");
        return 0;
    }

    if(final_state==0) {
        printf("error state swap: final state unknown\n");
        return 0;
    }

    if(initial_state==0) {
        printf("error state swap: initial state unknown\n");
        return 0;
    }

    path.block_stats[path.initial_state-1][prepi->index].mcacc[2].tries++;

    ////pick random replica of the final state;
    //jrep = (int)(RandomNumber()*(state[final_state-1].nrep-1)) +1;
    //if(jrep==0) {
    //    printf("Warning swap states: jrep is 0, should at least be 1\n");
    //}
    //prepj=&state[final_state-1].srep[jrep];  
    
    jrep = state[final_state-1].nrep-1;
    prepj =&state[final_state-1].srep[jrep];  
    
    double nrepratio;
    //nrepratio = (double)(state[final_state-1].nrep-1.)/(double)(state[initial_state-1].nrep-1.);
    nrepratio=1;
    if(path.stateswapbias==1) {
        aux = prepi->dos - prepj->dos;
        if (RandomNumber() > (nrepratio*exp(aux))) {
            return 0; 
        }
    }
    else if(path.stateswapbias==2) {
        aux = state[initial_state-1].dos-state[final_state-1].dos;
        if (RandomNumber() > (nrepratio*exp(aux))) {
            return 0; 
        }
    }
    else {
        if (RandomNumber() > nrepratio) {
            return 0; 
        }
    }
    //printf("Rejected state swap I %d J %d: nrepratio: %lf dosdiff: %lf\n",initial_state,final_state,nrepratio,aux);
    
    if (final_state == 1 || final_state == 2 ) {
        path.current_gsen = state[final_state-1].target.energy;
    }

    if (final_state == 3) {
        path.current_gsen =-sqrt(sys.site[slice[pathlen-1].minisite].eps*sys.site[slice[pathlen-1].minjsite].eps);
    }
    path.initial_state = final_state;
    path.final_state = initial_state;
    for(i=0;i< pathlen;i++) {
        trial[i]= slice[pathlen - 1 -  i];
        create_all_rc(&(trial[i]));
    }


    type = analyse(trial,prepj,pathlen,final_state);
    if ( type!=2 ) {
        //printf("Rejected state swap:wrong type -1 \n");
        
        //if(initial_state==6) {
        //  printf("this happens when final_state = 7\n");
        //  printf("not right type after swap\n");
        //  dprint(type);
        //  dprint(path.initial_state);
        //  dprint(initial_state);
        //  dprint(final_state);
        //  dprint(irep);
        //  dprint(jrep);
        //  gprint(path.current_gsen);
        //  gprint(state[initial_state-1].target.energy);
        //  gprint(state[final_state-1].target.energy);

        //  print_rc(&trial[0],path.initial_state);
        //  print_rc(&trial[len-1],path.initial_state);
        //  print_rc(&slice[0],path.initial_state);
        //  print_rc(&slice[len-1],path.initial_state);
        //  }
        path.current_gsen = state[initial_state-1].target.energy;
        path.initial_state = initial_state;
        path.final_state = final_state;
        return 0;
    }
  
    path.nreplica = state[final_state-1].nrep;

    for(i=0; i<path.nreplica; i++) {
        replica[i] = &state[path.initial_state-1].srep[i];
    }

    prepj->pathlen = pathlen;
    for(j=0; j<pathlen; j++) {
        slice[j] = trial[j];  
    //create_all_rc(&(slice[j]));
    }

    type = analyse(slice,prepj, prepj->pathlen, path.initial_state);
    prepj->type = type;
    path.current_replica=jrep;
 
    path.block_stats[path.initial_state-1][prepi->index].mcacc[2].acc++;

    return 1;
}



int reverse_replica(int irep) {

    Replica *prep;
    int i,j,pathlen,type;

    prep=replica[irep];
    path.block_stats[path.initial_state-1][prep->index].mcacc[3].tries++;
    if (irep==0) {
        return 0;  
    }

    type = analyse(slice,prep, path.nslices, path.initial_state);
    if (prep->type != type) {
        printf("error swap rep: type %d and replica type %d do not match\n",type,prep->type);
    }

    if (prep->type!=1) {
        return 0; //reversal can only happen for trajectories who end and start in the same state
    }

    pathlen =prep->pathlen;
    for(i=0;i< pathlen;i++) {
        trial[i]= slice[pathlen - 1 -  i];
    }

    type = analyse(trial,prep,pathlen,path.initial_state);
    if ( type==0) {
        return 0;
    }
    if ( type==2) {
        return 0;
    }

    for(j=0;j< pathlen;j++) {
        slice[j] = trial[j];
    }
  
    path.block_stats[path.initial_state-1][prep->index].mcacc[3].acc++;

    return 1;
}





