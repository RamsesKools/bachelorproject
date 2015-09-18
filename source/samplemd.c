#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"


#define maxbinangles 200

void sampleangles(int choice) {

    static int Pcostheta[3][maxbinangles],Pphi[3][maxbinangles];
    static int costhetatotal[3],phitotal[3],costhetabin,phibin,initial;
    static double costhetadbin, phidbin;
    static vector site[3];
    double aux1,aux2,aux3,sum[3];
    FILE *fp;
    int ibin,isite;

    tensor rotmat;
    vector u[3];

    if(choice==0) {
        if(initial==0) {
            initial=1;
            for( isite=0; isite<3; isite++) {
                for( ibin=0; ibin<maxbinangles; ibin++) {
                    Pcostheta[isite][ibin]=0;
                    Pphi[isite][ibin]=0;
                }
                costhetatotal[isite]=0;
                phitotal[isite]=0;
            }
            costhetadbin = 2.0/maxbinangles;
            phidbin = (2.0*PI)/maxbinangles;

            site[0].x=1.0;
            site[0].y=0.0;
            site[0].z=0.0;
            site[1].x=0.0;
            site[1].y=1.0;
            site[1].z=0.0;
            site[2].x=0.0;
            site[2].y=0.0;
            site[2].z=1.0;
        }

        
        rotmat = getrotmatrix(slice[0].pts[0].q);
        for( isite=0; isite<3; isite++) {
            matrix_x_vector(rotmat,site[isite],u[isite]);
        }

        for( isite=0; isite<3; isite++) {
            //measure the polar angle as the z component of the patch vector
            costhetabin = (int)((u[isite].z+1.0)/costhetadbin);
            if((costhetabin<maxbinangles)&&(costhetabin>=0)) {
                Pcostheta[isite][costhetabin]++;
                costhetatotal[isite]++;
            }
            else {
                printf("Warning: costhetabin exceeds total number of bins: costhetabin %d uz %lf\n", costhetabin, u[isite].z);
            }

            //measure the azimuthal angle as the tangent of u.y over u.x
            if(!((u[isite].x==0)&&(u[isite].y==0))) {
                //phibin = (int)((atan2(1.0,u[isite].x/u[isite].y))/phidbin);
                phibin = (int)((atan2(u[isite].y,u[isite].x)+PI)/phidbin);
                if((phibin<maxbinangles)&&(phibin>=0)) {
                    Pphi[isite][phibin]++;
                    phitotal[isite]++;
                }
                else {
                    printf("Warning: phibin exceeds total number of bins: phibin %d phi uy %lf ux %lf\n",phibin, u[isite].y, u[isite].x);
                }
            }
        }
    }



    if(choice==1) {

        fp=fopen("cospolarhist.dat","w");
        for( isite=0; isite<3; isite++) {
            sum[isite]=0;
        }
        for( ibin=0; ibin<maxbinangles; ibin++) {
            aux1 = (double)Pcostheta[0][ibin]/(double)costhetatotal[0];
            sum[0]+=aux1;
            aux2 = (double)Pcostheta[1][ibin]/(double)costhetatotal[1];
            sum[1]+=aux2;
            aux3 = (double)Pcostheta[2][ibin]/(double)costhetatotal[2];
            sum[2]+=aux3;
            fprintf(fp,"%lf %lf %lf %lf %d %d %d\n", (ibin+0.5)*costhetadbin - 1.0, aux1, aux2, aux3, Pcostheta[0][ibin], Pcostheta[1][ibin], Pcostheta[2][ibin]);
        }
        fclose(fp);

        printf("Total sum for Pcostheta in all three directions: %lf %lf %lf\n",sum[0],sum[1],sum[2]);

        for( isite=0; isite<3; isite++) {
            sum[isite]=0;
        }
        fp=fopen("azimuthhist.dat","w");
        for( ibin=0; ibin<maxbinangles; ibin++) {
            aux1 = (double)Pphi[0][ibin]/(double)phitotal[0];
            sum[0]+=aux1;
            aux2 = (double)Pphi[1][ibin]/(double)phitotal[1];
            sum[1]+=aux2;
            aux3 = (double)Pphi[2][ibin]/(double)phitotal[2];
            sum[2]+=aux3;
            fprintf(fp,"%lf %lf %lf %lf %d %d %d\n", (ibin+0.5)*phidbin, aux1, aux2, aux3, Pphi[0][ibin], Pphi[1][ibin], Pphi[2][ibin]);
        }
        fclose(fp);

        printf("Total sum for Pphi in all three directions: %lf %lf %lf\n",sum[0],sum[1],sum[2]);

    }

    return;
}





#define MAXT 7500
#define MAXT0 250
#define FREQT0 50
void samplediffusion(int choice) {


    int t0index, itime,CorrelTime;
    static int initial;
    static int time,t0time[MAXT0], t0Counter, SampleCounter[MAXT];
    static double MSD[MAXT],MSAD[MAXT];
    static vector R0[MAXT0],Psi0[MAXT0];
    static vector nonpbcr,rcurrent,rprev,rstep;
    static vector Psi,pprev;
    vector pcurrent, protax, protaxnorm, Psicurrent;
    double pmagangle, l;
    int isite,ibin;
    vector site[3],dr;
    tensor rotmat;
    double aux;
    double MSDCumulant, MSADCumulant,DCT,DCR;

    FILE *fmsd, *fmsad,*fdiff;


    if(choice==0) {
        if(initial==0) {
            initial=1;

            for( ibin=0; ibin<MAXT; ibin++) {
                SampleCounter[ibin] = 0;
                MSD[ibin] = 0;
                MSAD[ibin]= 0;
            }

            for( ibin=0; ibin<MAXT0; ibin++) {
                R0[ibin] = nulvec;
                Psi0[ibin] = nulvec;
                t0time[ibin] = 0;
            }

            time=0;
            t0Counter=0;

            nonpbcr = slice[0].pts[0].r;
            rprev = nonpbcr;
            rcurrent = nonpbcr;

            rotmat = getrotmatrix(slice[0].pts[0].q);
            matrix_x_vector(rotmat,sys.site[0].r,pprev);
            Psi = nulvec;
        }


        //get position without PBC
        rcurrent=slice[0].pts[0].r;
        vector_minus(rcurrent,rprev,rstep);
        pbc(rstep,sys.boxl);
        vector_add(nonpbcr,rstep,nonpbcr);
        rprev=rcurrent;


        //get Psi: unbounded step via norm rotation axis times angle
        rotmat=getrotmatrix(slice[0].pts[0].q);
        matrix_x_vector(rotmat,sys.site[0].r,pcurrent);
        pmagangle = acos(vector_inp(pprev,pcurrent));
        if((!(pmagangle!=pmagangle)) && (pmagangle>0.0000001)) {
            vector_cross(pprev,pcurrent,protax);
            l = sqrt(vector_inp(protax,protax));
            scalar_divide(protax,l,protaxnorm);
            scalar_times(protaxnorm,pmagangle,Psicurrent);
            vector_add(Psi,Psicurrent,Psi);
        }
        pprev=pcurrent;


        //printf("pmagangle %lf\n", pmagangle);
        //printf("protax %lf %lf %lf\n", protax.x, protax.y, protax.z);
        //printf("Psi %lf %lf %lf\n", Psi.x, Psi.y, Psi.z);
        if(Psi.x!=Psi.x) {
            exit(1);
        }


        if((time%FREQT0)==0) {
            t0index = t0Counter%MAXT0;
            t0Counter++;
            t0time[t0index] = time;
            R0[t0index] = nonpbcr;
            Psi0[t0index] = Psi;
        }
        time++;

        for( ibin=0; ibin<min(t0Counter,MAXT0); ibin++) {
            CorrelTime = time-t0time[ibin];
            if(CorrelTime<MAXT) {
                SampleCounter[CorrelTime]++;
                //update MSD
                vector_minus(nonpbcr,R0[ibin],dr);
                MSD[CorrelTime] += vector_inp(dr,dr);
                //update MSAD
                vector_minus(Psi,Psi0[ibin],dr);
                MSAD[CorrelTime] += vector_inp(dr,dr);
            }
        }
    }


    if(choice==1) {
        fmsd = fopen("msd.dat","w");
        MSDCumulant = 0;
        MSADCumulant = 0;

        for( itime=1; itime<(MAXT-1); itime++ ) {
            if(SampleCounter[itime]>0) {
                aux = MSD[itime] / (double)SampleCounter[itime];
                MSDCumulant += aux/(double)itime;
            }
            else {
                aux = 0;
            }
            fprintf(fmsd, "%d %lf %lf %lf %d\n", itime, itime*langevin.timestep, aux, MSD[itime],SampleCounter[itime]);
        }
        fclose(fmsd);

        fmsad = fopen("msad.dat","w");
        for( itime=1; itime<(MAXT-1); itime++ ) {
            if(SampleCounter[itime]>0) {
                aux = MSAD[itime] / (double)SampleCounter[itime];
                MSADCumulant += aux/(double)itime;
            }
            else {
                aux = 0;
            }
            fprintf(fmsad, "%d %lf %lf %lf %d\n", itime, itime*langevin.timestep, aux, MSAD[itime],SampleCounter[itime]);
        }
        fclose(fmsad);

        DCT = MSDCumulant/((MAXT-1.0)*langevin.timestep);
        DCT/= 6.0;
        DCR = MSADCumulant/((MAXT-1.0)*langevin.timestep);
        DCR/= 4.0;

        fdiff = fopen("diffusionconstants.dat","w");
        fprintf(fdiff,"%lf %lf %lf %lf %lf\n", langevin.dtD, sys.mobilityT, sys.mobilityR, DCT, DCR);
        fclose(fdiff);
    }

       



    return;
}














