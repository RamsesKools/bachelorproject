#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "path.h"



void update_wanglandau() {
	
    int i,j,k,l;
    int max,min;

	printf("Updating WL histograms\n");
	printf("Scalefactor \n");
	for(i=0; i<path.nstates; i++) {
		printf("%g   ",state[i].scalefactor);
		printf("%g:  ",state[i].dos);
		for(j=1; j<state[i].nrep; j++) {
			printf(" %5.2lf",state[i].srep[j].dos);
		}
		printf("\n");
	}
	printf("\n");
	
    if(path.stateswapbias==1) {
        max =0; min =10000000;
        for(l=0;l<path.nstates;l++) {
            for(k=1;k<state[l].nrep;k++) { 
                if (state[l].srep[k].ntotal > max) {
                    max = state[l].srep[k].ntotal;
                }
                if (state[l].srep[k].ntotal < min) {
                    min = state[l].srep[k].ntotal;
                }
		  	}
		}
        if ( ((double)max)/((double)min)  < 1.1) {
            printf("histogram flat within 10 percent now reducing scalefactor\n");
            for(i=0; i<path.nstates; i++) {
                state[i].scalefactor/=2.0;
                gprint(state[i].scalefactor);
            }
            printf("resetting histogram\n");
            for(l=0;l<path.nstates;l++) {
                for(k=0;k<state[l].nrep;k++) {
                    state[l].srep[k].ntotal=0;
                }
            } 
        }
        else {
            printf("histogram not flat yet:  min = %d, max = %d \n",min,max);
        }
    }
    else {
        for(i=0;i<path.nstates;i++) {
            printf("State %d\n",i);
            max =0; min =10000000;
            for(k=1; k<state[i].nrep; k++) {
                if (state[i].srep[k].ntotal > max) {
                    max = state[i].srep[k].ntotal;
                }
                if (state[i].srep[k].ntotal < min) {
                    min = state[i].srep[k].ntotal;
                }
            }
            if ( ((double)max)/((double)min)  < 1.1) {
                printf("histogram flat within 10 percent now reducing scalefactor\n");
                state[i].scalefactor/=2.0;
                gprint(state[i].scalefactor);
                printf("resetting histogram\n");
                for(k=0;k<state[i].nrep;k++) {
                    state[i].srep[k].ntotal=0;
                }
            }
            else {
                printf("histogram not flat yet:  min = %d, max = %d \n",min,max);
            }
        }
    }
		
	return;

}




void crossinghistogram(Replica *prep) {

    int irep,type,i,imax,imin,ibin,k,maxrep,final_state,initial_state;
	double lambdamax,lambdamin,z,weight;
	static int count;
  
	int state_index,istate,nrep,islice;
	Slice *psl;

	prep->ntotal++;
	prep->navlen++;
	prep->avlen+= prep->pathlen;
	prep->avlensq+= (prep->pathlen*prep->pathlen);
	irep = prep->index;

	type =  prep->type;
	state_index = path.initial_state-1;

	state[state_index].type_mat[irep][type-1]++; 

	update_average(path.block_stats[path.initial_state-1][prep->index].aver[1],prep->pathlen);

	if(type ==0) {
		printf("error: incorrect path in crossinghist replica %d\n",irep);
	}
	
	maxrep = path.current_replica;
	lambdamax = -1000;
	//if(prep->index==1) {
    //    lambdamax = prep->lambda;
    //}
	nrep = state[state_index].nrep;
	for (i=0; i<path.nslices; i++) {
		z = slice[i].order_parameter;
		if (z > lambdamax) {
			lambdamax =z;
			if((maxrep+1)<nrep) {
				if(z>state[state_index].lambda[maxrep+1]) {
					maxrep++;
				}
			}
		}
	}

	if(prep->index==0) {
        maxrep=0;
    }

	if(prep->index>0) {
		final_state = path.final_state-1;
		if(final_state<0) {
            printf("Cross: final state unknown: %d\n", final_state); 
        }
		state[state_index].pathtypenumbers[irep][final_state][maxrep]++;
	}

	k=irep;
	imin = (int) (replica[k]->lambda*50);
	imax = (int) (lambdamax*50);

	if(path.initial_state==path.nstates) {
		lambdamax= (lambdamax>0) ? (log(lambdamax) + 100) : 0 ;
		imax = (int)(lambdamax*50);
		lambdamin = (replica[k]->lambda >0) ? (log(replica[k]->lambda) +100) : 0;
		imin = (int)(lambdamin*50);
	}

	for(i=imin;i<=imax;i++) {     
		if((i>=0)&&(i<MAXBIN)) {
			state[state_index].crosshist[k][i]++;
		}
	}
  
	return;				  
}


void get_flux(Replica *prep) {

    int i, len,state_index;

    state_index=path.initial_state-1;

    if (prep->index ==0) {
        for(i=path.nslices-1;i>=0;i--) {
            if (in_state(&slice[i]) ){
                len = i;
                break;
            }
        }
        state[state_index].flux0 += len;
        state[state_index].nflux0++;
  	}

    if (prep->index ==1) {
        for(i=path.nslices-1;i>=0;i--) {
            if (in_upper_window(&slice[i],replica[1],path.initial_state)) {
				len = i;
	 			break;
            }
        }
        state[state_index].flux1+=len;
		state[state_index].nflux1++;
    }

    if (prep->index==path.nreplica-1 ) {
        get_omiflux(prep);
    }

	return;
}



void get_omiflux(Replica *prep) {

    int i,initial_state,final_state,state_index;
 
    initial_state=in_state(&slice[0]);
    final_state=in_state(&slice[path.nslices-1]);
  
    if (initial_state!=path.initial_state) {
        printf("error in path, intial state not equal to stateindex\n");
        dprint(initial_state);
        dprint(path.initial_state);
    }
  
    state[path.initial_state-1].mstis_mat[final_state-1]++;

    return;
}




void print_tisstats() {

    FILE *fp;
    char filename[100];
    double x;
    int istate,jstate,irep,jrep,ibin; 

    int mkdir = system("mkdir crosshists pathtypenumbers");

    for(istate=0;istate<path.nstates;istate++) {
        sprintf(filename,"crosshists/crosshist%d.dat",istate);
        if ((fp = fopen(filename,"w"))==NULL){
            printf("output: can't open file.dat\n");
        }
        else {
            for(irep=1;irep<state[istate].nrep;irep++) {
                for(ibin=0;ibin<MAXBIN;ibin++){
                    x = (double)(ibin)/50.;
                    fprintf(fp,"%lf %g\n",x,state[istate].crosshist[irep][ibin]);
                }
                fprintf(fp,"\n");
            }
            fclose(fp);
        }
    }

    for(istate=0; istate<path.nstates; istate++) {
        sprintf(filename,"pathtypenumbers/pathtype%d.dat",istate);
        if ((fp = fopen(filename,"w"))==NULL){
            printf("output: can't open file.dat\n");
            return;
        }
        else {
            for(irep=1;irep<state[istate].nrep;irep++) {
                for(jstate=0;jstate<path.nstates;jstate++) {
                    for(jrep=1; jrep<state[istate].nrep; jrep++) {
                        fprintf(fp,"%d %.8g %g\n",jstate,state[istate].lambda[jrep],state[istate].pathtypenumbers[irep][jstate][jrep]);
                    }
                }
                fprintf(fp,"\n");
            }
            fclose(fp);
        }
    }

    
    if ((fp = fopen("dos_all.dat","w"))==NULL){
        printf("output: can't open file.dat\n");
    }
    else {
        for(istate=0; istate<path.nstates; istate++) {
            fprintf(fp,"%.16lf\n",state[istate].scalefactor);
            for(irep=0; irep<state[istate].nrep; irep++) {
                fprintf(fp,"%d %.16lf\n",irep,state[istate].srep[irep].dos);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }


    return;
}






void reweight_pathtypenumbers() {
	
	int istate,ibin,irep,jrep,j,k,nrep;
	int jstate,jbin,kstate,kbin,krep;
	double sumH[MAXSTATES][MAXREPLICA];
	double nsum[MAXREPLICA];
	double H[MAXREPLICA][MAXSTATES][MAXREPLICA];
	double difference,iteration;
	double finalH[MAXREPLICA][MAXSTATES][MAXREPLICA];
	double combinedH[MAXSTATES][MAXREPLICA];
	double crossingprob[MAXREPLICA];
	double totpathtypes=0;
	double ratematrix[MAXSTATES][MAXSTATES];
	double sumdum,weight,norm,maxweight;
	double value,maxvalue,maxlambda;
	double flux;
	int maxrep,maxstate;
	double weights[MAXREPLICA];
	double weightstate[MAXSTATES];
	FILE *fp;
	char filename[100];	
	double x,y;
	double pathnummatrix[path.nstates][path.nstates];

    int mkdir = system("mkdir rpe_pathtypes");
	//initiate values from crossing histogram
	printf("Reweighting path-type number\n");
	for(istate=0; istate<path.nstates; istate++) {
		printf("State %d\n",istate);
		nrep=1;

		//set everything to zero to be sure
		for(irep=1; irep<MAXREPLICA; irep++) {
			nsum[irep]=0;
			maxvalue=0;

			for(jstate=0; jstate<MAXSTATES; jstate++) {
				for(krep=1; krep<MAXREPLICA; krep++) {
					H[irep][jstate][krep]=0;
					sumH[jstate][krep]=0;
				}
			}
		}

		//read in pathtypenumbers in H and obtain maxvalue
		for(irep=1; irep<state[istate].nrep; irep++) {
			nsum[irep]=0;
			maxvalue=0;

			for(jstate=0; jstate<path.nstates; jstate++) {
				for(krep=1; krep<state[istate].nrep; krep++) {
					H[irep][jstate][krep] = 0;
					value = state[istate].pathtypenumbers[irep][jstate][krep];
					H[irep][jstate][krep] = value;
					if(maxvalue<value) maxvalue=value;
				}
			}

			//only allow numbers which have good statistics
			for(jstate=0; jstate<path.nstates; jstate++) {
				for(krep=1; krep<=state[istate].nrep; krep++) {
					if(H[irep][jstate][krep]<0.1*maxvalue) {
						H[irep][jstate][krep]=0;
					}
					//nsum = Mi
					nsum[irep]+=H[irep][jstate][krep];
				}
			}

			//only count this replica if it actually has a significant amount of data
			if(nsum[irep]>0) {
				nrep++;
			}

		}

		//printf("state %d nrep %d\n",istate,nrep);
		//printf("made 2D histogram from crossing hist pathtype\n");
		for(jstate=0; jstate<path.nstates; jstate++) {
			for(krep=1; krep<state[istate].nrep; krep++) {
				sumH[jstate][krep]=0;
				for(irep=1; irep<nrep; irep++) {
					sumH[jstate][krep]+=H[irep][jstate][krep];
				}
			}
		}

		//printf("done sumH\n");


		//Weights are from the crossing probability WHAM
		for(irep=1; irep<state[istate].nrep; irep++) {
			weights[irep] = state[istate].srep[irep].weight;
		}


		//instead of exp(lnZ[irep])/nsum[irep] use weight_irep
		for(irep=1; irep<nrep; irep++) {
			for(jstate=0; jstate<path.nstates; jstate++) {
				for(krep=1; krep<state[istate].nrep; krep++) {
					// final pi(Q) = Hi(Q) * Zi / Mi eq. p184 Frenkel
					finalH[irep][jstate][krep]=H[irep][jstate][krep]*weights[irep];
				}
			}
		}
		//printf("created final hist\n");
			


		//printf("Writing Path-type histograms\n");
		sprintf(filename,"rpe_pathtypes/allhists%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			for(jstate=0; jstate<path.nstates; jstate++) {
				for(krep=1; krep<state[istate].nrep; krep++) {
					fprintf(fp,"%d %d %12.9lf %12.9lf\n",
					 jstate,krep,state[istate].srep[krep].lambda,finalH[irep][jstate][krep]);
				}
				fprintf(fp,"\n");
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		//printf("print all hist\n");





		
		//Reweight pathtypenumber matrix according to equation 7: Bolhuis, Du 2014
		//w_k = 1 / ( Sum_l=1_k  1/w_l )
		//multiply every point with w_k
		//
		//again do not use nsump[irep]/exp(lnZ[irep]) but weights from crossing prob
		double sumweight;
		for(jstate=0; jstate<path.nstates; jstate++) {

			for(krep=1; krep<state[istate].nrep; krep++) {
				sumweight=0;
				//important that irep only goes until krep

				for(irep=1; irep<=krep; irep++) {
					//maybe the following is not usefull anymore
					//if(H[irep][jstate][krep]>0) {
					//	weight=1;
					//	}
					//else {
					//	weight =0;
					//	}
					weight=1;
					//sumdum = Sum ( exp(-Beta * Wi) * Mi / (Zi/Z0) ) for i 0 to n
					//sumdum+=weight*nsum[irep]/exp(lnZ[irep]);
					if(weights[irep]>0) {
						sumweight += 1./weights[irep];
					}
				}

				if(sumweight>0) {
					combinedH[jstate][krep] = weight*sumH[jstate][krep]/sumweight;
				}
				else {
					combinedH[jstate][krep]=0;
				}

			}
		}
		//printf("combined histograms\n");



		//ADD COMAXBINEDH to PATHNUMMATRIX
		//combinedH[jstate][krep] = wbark SUM_i nIJ^i(lambda_k)
		//pathnummatrix[istate][jstate] = SUM_k combinedH[jstate][krep]
		//eq. 7-8 Weina Bolhuis article
		for(jstate=0; jstate<path.nstates; jstate++) {
			pathnummatrix[istate][jstate]=0;
			for(krep=1; krep<state[istate].nrep; krep++) {
				pathnummatrix[istate][jstate]+=combinedH[jstate][krep];
			}
		}




		norm = 0;
		for(jstate=0; jstate<path.nstates; jstate++) {
			for(krep=1; krep<state[istate].nrep; krep++) {
				if(combinedH[jstate][krep]>norm) norm = combinedH[jstate][krep];
				}
			}
		if(norm<=0) printf("something wrong with combinedH, norm is %lf\n",norm);
				
		sprintf(filename,"rpe_pathtypes/hist%d.dat",istate);
		fp=fopen(filename,"w");
		for(jstate=0; jstate<path.nstates; jstate++) {
			for(krep=1; krep<state[istate].nrep; krep++) {
				//combinedH[jstate][krep]/=norm;
				y = combinedH[jstate][krep];
				fprintf(fp,"%d %d %12.9lf %12.9lf\n",jstate,krep,state[istate].srep[krep].lambda,y);
				}
			}
		fclose(fp);
		//printf("wrote combihistogram\n");




		// calculate crossing probability
		totpathtypes=0;
		for(jstate=0; jstate<path.nstates; jstate++) {
			for(krep=1; krep<state[istate].nrep; krep++) {
				totpathtypes+=combinedH[jstate][krep];
				}
			}
		for(irep=1; irep<nrep; irep++) {
			crossingprob[irep]=0;
			for(jstate=0; jstate<path.nstates; jstate++) {
				for(krep=irep; krep<state[istate].nrep; krep++) {
					crossingprob[irep]+=combinedH[jstate][krep];
					}
				}
			}
		for(irep=1; irep<nrep; irep++) {
			crossingprob[irep]/=totpathtypes;
			}




		sprintf(filename,"rpe_pathtypes/crossprob%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			fprintf(fp,"%12.9g %12.9g\n", state[istate].srep[irep].lambda, crossingprob[irep]);
			}
		fclose(fp);
		//printf("wrote crossingprob\n");

		sprintf(filename,"rpe_pathtypes/logcrossprob%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			if(crossingprob>0) x=log(crossingprob[irep]);
			else{ x=-10;}
			fprintf(fp,"%12.5g %12.5g\n", state[istate].srep[irep].lambda, x);
			}
		fclose(fp);
		//printf("wrote logcrossingprob\n");
		//all states are reweighted
		}





		//calculate pathnumber matrix
		//wham of the pathnumber matrices
		//check to see if they become symmetric
		//normalization pathnumbers and multiplication with flux gives rate
		// calculate rates


		double Hist[path.nstates][path.nstates*path.nstates];
		double pathnumsym[path.nstates][path.nstates][path.nstates];
		int ii,jj;
		double nstatesum[path.nstates];
		double sumstatesH[path.nstates*path.nstates];
		double Znew,lnZnew[path.nstates],lnZ[path.nstates];

		for(istate=0; istate<path.nstates; istate++) {
			for(jstate=0; jstate<path.nstates; jstate++) {
				for(kstate=0; kstate<path.nstates; kstate++) {
					pathnumsym[istate][jstate][kstate]=0;
					}
				}
			}

		//fill up the rows
		for(istate=0; istate<path.nstates; istate++) {
			for(jstate=0; jstate<path.nstates; jstate++) {
				pathnumsym[istate][istate][jstate]=pathnummatrix[istate][jstate];
				}
			}
		//fill up the columns
		for(istate=0; istate<path.nstates; istate++) {
			for(jstate=0; jstate<path.nstates; jstate++) {
				pathnumsym[istate][jstate][istate]=pathnumsym[istate][istate][jstate];
				}
			}

		sprintf(filename,"rpe_pathtypes/pathnumcolumnrows.dat");
		fp=fopen(filename,"w");
		for(ii=0; ii<path.nstates; ii++) {
			fprintf(fp," \n");
			for(istate=0; istate<path.nstates; istate++) fprintf(fp,"              %d",istate+1);
			fprintf(fp,"\n");
			for(istate=0; istate<path.nstates; istate++) {
				fprintf(fp,"%d ",istate+1);
				for(jstate=0; jstate<path.nstates; jstate++) {
					fprintf(fp,"& %12.5g ", pathnumsym[ii][istate][jstate]);
					}
				fprintf(fp,"\\\\ \n");
				}
			}
		fclose(fp);


		for(ii=0; ii<path.nstates; ii++) {
			nstatesum[ii]=0;
			lnZnew[ii]=0;
			lnZ[ii]=0;
			jj=0;
			for(istate=0; istate<path.nstates; istate++) {
				for(jstate=0; jstate<path.nstates; jstate++) {
					Hist[ii][jj]=pathnumsym[ii][istate][jstate];
					nstatesum[ii]+=Hist[ii][jj];
					jj++;
					}
				}
			}
		//printf("Converted matrix into one dimensional matrix\n");

		sprintf(filename,"rpe_pathtypes/histpathnum.dat");
		fp=fopen(filename,"w");
		for(ii=0; ii<path.nstates; ii++) {
			fprintf(fp," \n");
			for(istate=0; istate<path.nstates; istate++) fprintf(fp,"              %d",istate+1);
			fprintf(fp,"\n");
			jj=0;
			for(istate=0; istate<path.nstates; istate++) {
				fprintf(fp,"%d ",istate+1);
				for(jstate=0; jstate<path.nstates; jstate++) {
					fprintf(fp,"& %12.5g ", Hist[ii][jj]);
					jj++;
					}
				fprintf(fp,"\\\\ \n");
				}
			}
		fclose(fp);



		for(istate=0; istate<(path.nstates*path.nstates); istate++) {
			sumstatesH[istate]=0;
			for(jstate=0; jstate<path.nstates; jstate++) {
				sumstatesH[istate]+=Hist[jstate][istate];
				}
			}
		//printf("and created nstatesum and sumstatesH\n");



		fp=fopen("rpe_pathtypes/sumstatesH.dat","w");
		for(istate=0; istate<path.nstates; istate++) fprintf(fp,"              %d",istate+1);
		fprintf(fp,"\n");
		jj=0;
		for(istate=0; istate<path.nstates; istate++) {
			fprintf(fp,"%d ",istate+1);
			for(jstate=0; jstate<path.nstates; jstate++) {
				fprintf(fp,"& %12.5g ", sumstatesH[jj]);
				jj++;
				}
			fprintf(fp,"\\\\ \n");
			}
		fprintf(fp,"\n\n");
		for(istate=0; istate<(path.nstates*path.nstates); istate++) {
			fprintf(fp,"%lf\n",sumstatesH[istate]);
			}
		fclose(fp);

		fp=fopen("rpe_pathtypes/nstatesum.dat","w");
		for(istate=0; istate<(path.nstates); istate++) {
			fprintf(fp,"%lf\n",nstatesum[istate]);
			}
		fclose(fp);

		difference=1000;
		iteration=0;
		while((difference>0.0000001) && (iteration<10000000)) {

			for(istate=0; istate<(path.nstates); istate++) {
				Znew=0;
				//integrate over all the bins

				for(jstate=0; jstate<(path.nstates*path.nstates); jstate++) {

					if(Hist[istate][jstate]>0) {
						sumdum=0;
						
						for(kstate=0; kstate<path.nstates; kstate++) {

							if(Hist[kstate][jstate]>0) {
								//exp(-Beta Wi) only 1 if correct replica
								weight=1.;
								}
							else {
								weight =0.;
								}

							//sumdum = Sum ( exp(-Beta * Wk) * Mk / Zki ) for k 0 to n
							sumdum+=(weight*nstatesum[kstate]/exp(lnZ[kstate]));
							}
						//Znew = Sum ( Sum(Hk) for k 0 to n / sumdum )
						if(sumdum>0) {
							Znew+=(sumstatesH[jstate]/sumdum);
							}
						}
					}

				if(Znew>0) {
					lnZnew[istate]=log(Znew);
					}

				}
			difference=0;
			for(istate=0; istate<path.nstates; istate++) {
				difference+=fabs(lnZ[istate] - lnZnew[istate]);
				//always relative to Z1
				lnZ[istate] = lnZnew[istate]-lnZnew[0];
				if(lnZ[istate]!=lnZ[istate]) {
					printf("lnZ[%d] = %lf is not a number\n",istate,lnZ[istate]);
					printf("lnZnew[%d] = %lf\n",istate,lnZnew[istate]);
					}
				//printf("lnZnew[%d]\n",irep);
				}
			iteration++;
			}
		printf("Iteration %lf\n",iteration);
		printf("difference %lf\n",difference);
		printf("performed iteration\n");


		double combinedstateH[path.nstates*path.nstates];
		for(istate=0; istate<(path.nstates*path.nstates); istate++) {
			sumdum=0;

			for(jstate=0; jstate<path.nstates; jstate++) {

				if(Hist[jstate][istate]>0) {
					weight=1;
					}
				else {
					weight =0;
					}

				//sumdum = Sum ( exp(-Beta * Wi) * Mi / (Zi/Z0) ) for i 0 to n
				sumdum+=(weight*nstatesum[jstate]/exp(lnZ[jstate]));
				}

			if(sumdum>0) {
				combinedstateH[istate] = sumstatesH[istate]/sumdum;
				}
			else{
				combinedstateH[istate]=0;
				}

			}


		maxweight = -1;
		for(istate=0; istate<path.nstates; istate++) {
			//Weighti = Zi / Mi
			weightstate[istate] = exp(lnZ[istate])/nstatesum[istate];
			if(maxweight<weightstate[istate]) maxweight = weightstate[istate];
			}

		sprintf(filename,"rpe_pathtypes/stateweights.dat");
		fp=fopen(filename,"w");
		for(istate=0; istate<path.nstates; istate++) {
			weightstate[istate]/=maxweight;
			}

		for(istate=0; istate<path.nstates; istate++) {
			x = log(weightstate[istate]);
			fprintf(fp,"%d %lf %lf %lf\n",
			 istate, x, weightstate[istate], nstatesum[istate]);
			state[istate].weight = weightstate[istate];
			state[istate].lnweight = x;
			}
		fclose(fp);




		double statematrix[path.nstates][path.nstates];
		double totpaths[path.nstates];
		ii=0; 
		for(istate=0; istate<path.nstates; istate++) {
			totpaths[istate]=0;
			for(jstate=0; jstate<path.nstates; jstate++) {
				statematrix[istate][jstate]=combinedstateH[ii];
				totpaths[istate]+=statematrix[istate][jstate];
				ii++;
				}
			}

		double f0,f1;

		sprintf(filename,"rpe_pathtypes/flux.dat");
		fp=fopen(filename,"w");
		fprintf(fp," \n");
		fprintf(fp,"Fluxes    ");
		fprintf(fp,"0     1      total        flux\n");
		for(istate=0; istate<path.nstates; istate++) {
			f0 = (double)state[istate].flux0/state[istate].nflux0;		
			f1 = (double)state[istate].flux0/state[istate].nflux0;		
			flux = 1./(f0+f1);
			fprintf(fp,"%d %lf %lf %lf %lf\n",istate,f0,f1,f0+f1,1./(f0+f1));
			}


		sprintf(filename,"rpe_pathtypes/pathnummatrix.dat");
		fp=fopen(filename,"w");
		fprintf(fp," \n");
		for(istate=0; istate<path.nstates; istate++) fprintf(fp,"              %d",istate+1);
		fprintf(fp,"\n");
		for(istate=0; istate<path.nstates; istate++) {
			fprintf(fp,"%d ",istate+1);
			for(jstate=0; jstate<path.nstates; jstate++) {
				fprintf(fp," %12.5g ", statematrix[istate][jstate]);
				}
			fprintf(fp,"\n");
			}
		fclose(fp);


		sprintf(filename,"rpe_pathtypes/ratematrix.dat");
		fp=fopen(filename,"w");
		//fprintf(fp," \n");
		//for(istate=0; istate<path.nstates; istate++) fprintf(fp,"           %d",istate+1);
		//fprintf(fp,"\n");
		for(istate=0; istate<path.nstates; istate++) {
			f0 = (double)state[istate].flux0/state[istate].nflux0;		
			f1 = (double)state[istate].flux0/state[istate].nflux0;		
			flux = 1./(f0+f1);
			//fprintf(fp,"%d ",istate+1);
			for(jstate=0; jstate<path.nstates; jstate++) {
				fprintf(fp,"%12.5g ", flux*statematrix[istate][jstate]/totpaths[istate]);
				}
			fprintf(fp,"\n");
			}
		fclose(fp);
		//printf("wrote rate matrix\n");

		sprintf(filename,"rpe_pathtypes/ratematrixlatex.dat");
		fp=fopen(filename,"w");
		fprintf(fp,"\n");
		for(istate=0; istate<path.nstates; istate++) {
			f0 = (double)state[istate].flux0/state[istate].nflux0;		
			f1 = (double)state[istate].flux0/state[istate].nflux0;		
			flux = 1./(f0+f1);
			fprintf(fp,"%d ",istate+1);
			for(jstate=0; jstate<path.nstates; jstate++) {
				fprintf(fp,"& %12.5g ", flux*statematrix[istate][jstate]/totpaths[istate]);
				}
			fprintf(fp,"\\\\ \n");
			}
		fclose(fp);


	return;
}



void reweight_crosshist() {
	
	int istate,ibin,irep,jrep,j,k,nrep;
	double sumH[MAXBIN];
	double nsum[MAXREPLICA];
	double H[MAXREPLICA][MAXBIN];
	double difference,iteration;
	double Znew,lnZnew[MAXREPLICA],lnZ[MAXREPLICA];
	double finalH[MAXREPLICA][MAXBIN];
	double combinedH[MAXBIN];
	double weights[MAXREPLICA];
	double sumdum,weight,norm,maxweight;
	double value,maxvalue,maxlambda;
	FILE *fp;
	char filename[100];	
	double x,y;

//initiate values from crossing histogram
	int mkdir = system("mkdir rpe_crosshist");
	printf("Reweighting Crossing Probabilities\n");
	for(istate=0; istate<path.nstates; istate++) {
		printf("State %d\n",istate);
		nrep=1;
		for(irep=1; irep<MAXREPLICA; irep++) {
			lnZnew[irep]=0;
			lnZ[irep]=0;
			nsum[irep]=0;
			maxvalue=0;
			for(ibin=0; ibin<MAXBIN; ibin++) {
				H[irep][ibin] = 0;
				}
			}
		for(irep=1; irep<state[istate].nrep; irep++) {
			for(ibin=0; ibin<MAXBIN; ibin++) {
				value = state[istate].crosshist[irep][ibin];
				H[irep][ibin] = value;
				if(maxvalue<value) {
					maxvalue=value;
					state[istate].maxpaths=maxvalue;
					}
				}
			for(ibin=0; ibin<MAXBIN; ibin++) {
				if(H[irep][ibin]<0.1*maxvalue) {
					H[irep][ibin]=0;
					}
				//nsum = Mi
				nsum[irep]+=H[irep][ibin];
				}
			//only count this replica if it actually has a significant amount of data
			if(nsum[irep]>0) {
				nrep+=1;
				}
			}
		//printf("Considering %d replica's for state %d\n",nrep,istate);
		//printf("made histogram from crossing hist\n");

		for(ibin=0; ibin<MAXBIN; ibin++) {
			//sumH is Sum Hi for i is 0 to n
			sumH[ibin]=0;
			for(irep=1; irep<nrep; irep++) {
				sumH[ibin]+=H[irep][ibin];
				}
			}

		sprintf(filename,"rpe_crosshist/sumH%d.dat",istate);
		fp = fopen(filename,"w");
		for(ibin=0; ibin<MAXBIN; ibin++) {
			fprintf(fp,"%lf %lf\n",(double)ibin/50.,sumH[ibin]);
			}
		fclose(fp);
		//printf("done sumH\n");

		sprintf(filename,"rpe_crosshist/nsum%d.dat",istate);
		fp = fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			fprintf(fp,"%d %lf\n",irep,nsum[irep]);
			}
		fclose(fp);
		//printf("done nsum\n");

//calculate partition function
		difference=1000;
		iteration=0;
		while((difference>0.0000001) && (iteration<10000)) {
			for(irep=1; irep<nrep; irep++) {
				Znew=0;
				//integrate over all the bins
				for(ibin=0; ibin<MAXBIN; ibin++) {
					if(H[irep][ibin]>0) {
						sumdum=0;
						for(jrep=1; jrep<nrep; jrep++) {
							if(H[jrep][ibin]>0) {
								//exp(-Beta Wi) only 1 if correct replica
								weight=1.;
								}
							else {
								weight =0.;
								}
							//sumdum = Sum ( exp(-Beta * Wk) * Mk / Zki ) for k 0 to n
							sumdum+=(weight*nsum[jrep]/exp(lnZ[jrep]));
							}
						//Znew = Sum ( Sum(Hk) for k 0 to n / sumdum )
						Znew+=(sumH[ibin]/sumdum);
						}
					}
				if(Znew==0) {
					printf("Something wrong with summation\n");
					dprint(irep);
					gprint(Znew);
					}
				lnZnew[irep]=log(Znew);
				}
			difference=0;
			for(irep=1; irep<nrep; irep++) {
				difference+=fabs(lnZ[irep] - lnZnew[irep]);
				//always relative to Z1
				lnZ[irep] = lnZnew[irep]-lnZnew[1];
				if(lnZ[irep]!=lnZ[irep]) {
					printf("lnZ[%d] = %lf\n",irep,lnZ[irep]);
					printf("lnZnew[%d] = %lf\n",irep,lnZnew[irep]);
					}
				//printf("lnZnew[%d]\n",irep);
				}
			iteration++;
			}
		printf("Iteration %lf\n",iteration);
		printf("difference %lf\n",difference);
		//printf("performed iteration\n");

		for(irep=1; irep<nrep; irep++) {
			for(ibin=0; ibin<MAXBIN; ibin++) {
				// final pi(Q) = Hi(Q) * Zi / Mi eq. p184 Frenkel
				finalH[irep][ibin]=H[irep][ibin]*exp(lnZ[irep])/nsum[irep];
				}
			}
		//printf("created final hist\n");
			

		sprintf(filename,"rpe_crosshist/allhists%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			for(ibin=0; ibin<MAXBIN; ibin++) {
				x = (double)(ibin)/50.;
				fprintf(fp,"%lf %lf\n",x,finalH[irep][ibin]);
				}
			fprintf(fp,"\n");
			}
		fclose(fp);
		//printf("print all hist\n");
				
		sprintf(filename,"rpe_crosshist/allloghists%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			for(ibin=0; ibin<MAXBIN; ibin++) {
				x = (double)(ibin)/50.;
				if(finalH[irep][ibin]>0) y = log(finalH[irep][ibin]);
				else{  y=-10;}
				fprintf(fp,"%lf %lf\n",x,y);
				}
			fprintf(fp,"\n");
			}
		fclose(fp);
		//printf("print alllog hist\n");

		
		//combining the histograms according to eq 7.3.9 page 186 Frenkel
		for(ibin=0; ibin<MAXBIN; ibin++) {
			sumdum=0;
			for(irep=1; irep<nrep; irep++) {
				if(H[irep][ibin]>0) {
					weight=1;
					}
				else {
					weight =0;
					}
				//sumdum = Sum ( exp(-Beta * Wi) * Mi / (Zi/Z0) ) for i 0 to n
				sumdum+=(weight*nsum[irep]/exp(lnZ[irep]));
				}
			if(sumdum>0) {
				combinedH[ibin] = sumH[ibin]/sumdum;
				}
			else{
				combinedH[ibin]=0;
				}
			}
		//printf("combined histograms\n");
				
		sprintf(filename,"rpe_crosshist/hist%d.dat",istate);
		fp=fopen(filename,"w");
		norm=0;
		for(ibin=0; ibin<MAXBIN; ibin++) {
			if(combinedH[ibin]>norm) norm=combinedH[ibin];
			}
		for(ibin=0; ibin<MAXBIN; ibin++) {
			x = (double)(ibin)/50.;
			combinedH[ibin]/=norm;
			y=combinedH[ibin];
			fprintf(fp,"%lf %lf\n",x,y);
			}
		fclose(fp);
		//printf("wrote combihistogram\n");
				
		sprintf(filename,"rpe_crosshist/loghist%d.dat",istate);
		fp=fopen(filename,"w");
		norm=0;
		for(ibin=0; ibin<MAXBIN; ibin++) {
			if(combinedH[ibin]>0) {
				x = (double)(ibin)/50.;
				y=log(combinedH[ibin]);
				fprintf(fp,"%lf %lf\n",x,y);
				}
			}
		fclose(fp);
		//printf("wrote logcombihistogram\n");
		//

		int ilambda;
		double lambdarep;
		for(irep=1; irep<state[istate].nrep; irep++) {
			if(istate<8) {
				ilambda=(int)(state[istate].srep[irep].lambda*50);
				if(log(combinedH[ilambda])==log(combinedH[ilambda])) {
					state[istate].srep[irep].logcrossprob = log(combinedH[ilambda]);
					}
				}
			else {
				lambdarep= (state[istate].srep[irep].lambda>0) ? (log(state[istate].srep[irep].lambda) +100):0;
				ilambda=(int)(lambdarep*50);
				if(log(combinedH[ilambda])==log(combinedH[ilambda])) {
					state[istate].srep[irep].logcrossprob = log(combinedH[ilambda]);
					}
				}
			}

		sprintf(filename,"rpe_crosshist/zvalues%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			fprintf(fp,"%lf\n",lnZ[irep]);
			}
		fclose(fp);
		//printf("wrote zvalues\n");


		maxweight = -1;
		for(irep=1; irep<nrep; irep++) {
			//Weighti = Zi / Mi
			weights[irep] = exp(lnZ[irep])/nsum[irep];
			if(maxweight<weights[irep]) maxweight = weights[irep];
			}

		sprintf(filename,"rpe_crosshist/weights%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=1; irep<nrep; irep++) {
			weights[irep]/=maxweight;
			}

		for(irep=1; irep<nrep; irep++) {
			x = log(weights[irep]);
			fprintf(fp,"%d %lf %lf %lf %lf\n",
			 state[istate].srep[irep].index, state[istate].srep[irep].lambda, x, weights[irep], nsum[irep]);

			state[istate].srep[irep].weight = weights[irep];
			state[istate].srep[irep].lnweight = x; 
			}

		//add bias or not?  NOPE get the bias from the pathtype analysis.
		//replica's are chosen randomly, not only on outermost interface
		state[istate].srep[0].weight = state[istate].srep[1].weight;
		state[istate].srep[0].lnweight = state[istate].srep[1].lnweight;
		fclose(fp);
		//printf("wrote weigts\n");

		sprintf(filename,"rpe_crosshist/maxheights%d.dat",istate);
		fp=fopen(filename,"w");
		for(irep=0; irep<nrep; irep++) {
			maxvalue=0;
			maxlambda=0;
			for(ibin=0; ibin<MAXBIN; ibin++) {
				if(H[irep][ibin]>maxvalue) {
					maxvalue=H[irep][ibin];
					maxlambda = (double)(ibin)*50.;
					}
				}
			fprintf(fp,"%lf %lf\n",maxlambda,maxvalue);
			}
		fclose(fp);
		//printf("wrote maxheights\n");
		}

		sprintf(filename,"rpe_crosshist/logcrossprob.dat");
		fp=fopen(filename,"w");
		for(istate=0; istate<path.nstates; istate++) {
			for(irep=1; irep<state[istate].nrep; irep++) {
				fprintf(fp,"%lf\n",state[istate].srep[irep].logcrossprob);
				}
			fprintf(fp,"\n");
			}
		fclose(fp);



	return;
}







//void write_paths() {
//
//    int irep,j,i,k,type;
//    int initial_state,final_state,is,fs;
//    Replica *prep;
//
//    irep = path.current_replica;
//    prep=replica[irep];
//    type = prep->type;
//    initial_state = path.initial_state;
//    is = initial_state-1;
//    //final_state = in_state(&slice[prep->pathlen-1]);
//    fs = path.final_state-1;
//    if(irep==0) {
//        fs = is;
//    }
//
//    if((fs>=0)||(is>=0)) {
//        fprintf(path.filepath[is][fs],"\n");
//        fprintf(path.filepath[is][fs],"%d %d %d %d slices,currentrep,is,fs,type\n",
//			path.nslices,path.current_replica,path.initial_state,type);
//		fprintf(path.filepath[is][fs],"%lf %lf %lf boxl\n",sys.boxl.x,sys.boxl.y,sys.boxl.z);
//
//        for(j=0;j<prep->pathlen;j++) {
//            fprintf(path.filepath[is][fs],"Energy %lf\n",slice[j].energy);
//            for(i=0;i<sys.npart;i++){
//                fprintf(path.filepath[is][fs],"%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",
//                  slice[j].pts[i].r.x,slice[j].pts[i].r.y,slice[j].pts[i].r.z,
//                  slice[j].pts[i].q.r,slice[j].pts[i].q.u.x,slice[j].pts[i].q.u.y,slice[j].pts[i].q.u.z);
//			}
//		}
//		fflush(NULL);
//    }
//
//    return;
//}
//  
//
//
//void write_errorstate(Slice *psl, int index) {
//
//    int irep,j,i,k,type;
//    int initial_state,final_state,is;
//
//    fprintf(path.filepatherr[index],"\n");
//    fprintf(path.filepatherr[index],"%d %d %d %d slices,currentrep,is,fs,type\n",
//        1,path.current_replica,path.initial_state,0);
//	fprintf(path.filepatherr[index],"%lf %lf %lf boxl\n",sys.boxl.x,sys.boxl.y,sys.boxl.z);
//
//	fprintf(path.filepatherr[index],"Energy %lf\n",psl->energy);
//	for(i=0;i<sys.npart;i++) {
//		fprintf(path.filepatherr[index],"%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",
//	 	 psl->pts[i].r.x,psl->pts[i].r.y,psl->pts[i].r.z,
//	 	 psl->pts[i].q.r,psl->pts[i].q.u.x,psl->pts[i].q.u.y,psl->pts[i].q.u.z);
//    }
//    fflush(NULL);
//
//    return;
//}
//
//
//void write_swap() {
//
//	int k;
//	int istate;
//	static int count[9],countstate;
//
//	istate = path.initial_state-1;
//
//	fprintf(path.fileswap[istate],"%d ",count[istate]++);
//	for(k=0;k<path.nreplica;k++) {
//        fprintf(path.fileswap[istate],"%d ",replica[k]->swapindex);
//    }
//	fprintf(path.fileswap[istate],"\n");   
//
//	fprintf(path.filestateswap,"%d %d\n",countstate++,istate);
//
//	return;
//}
//
///
