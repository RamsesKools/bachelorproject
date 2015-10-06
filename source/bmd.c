#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"




double potential_energy(Slice *psl) {

    Pts *psi,*psj;
    vector rij,p[sys.npart][sys.nsites],ui,uj,rnorm;
    double r2,r6,r6inv,r12inv,potential_energy,potential_energy_attP=0,potential_energy_attC=0,potential_energy_rep=0,potential_energy_attP_d=0,r,cositheta,cosjtheta,epsmix;
    double phi_i,phi_j;
    int ipart,jpart,isite,jsite;
    tensor rotmati,rotmatj;

    psl->energyT=0;
    psl->energyN=0;
    psl->rijdist=0;
    psl->minisite=0;
    psl->minjsite=0;

    for( ipart=0; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        rotmati = getrotmatrix(psi->q);
        for( isite=0; isite<sys.nsites; isite++) {
            matrix_x_vector(rotmati,sys.site[isite].r,p[ipart][isite]);
        }
    }
    //calculate Lennard-Jones interaction between all particles
    for( ipart=1; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        for( jpart=0; jpart<ipart; jpart++ ) {
            psj = &psl->pts[jpart];
            vector_minus(psi->r,psj->r,rij);
            pbc(rij,sys.boxl);
            r2 = vector_inp(rij,rij);
            psl->rijdist=r2;
            if(r2<4.0) {
                r6 = r2*r2*r2;
                r6inv=1.0/r6;
                r12inv = r6inv*r6inv;
                if(r2<sys.sigmaLJsq) {
                    potential_energy_rep = r12inv-r6inv+0.25;
                }
                //calculate isotropic potential
                //potential_energy_attC += r12inv-r6inv;
                r=sqrt(r2);
                psl->rijdist = r;
                scalar_divide(rij,r,rnorm);
                //calculate angular part of potential
                for( isite=0; isite<sys.nsites; isite++ ) {
                    cositheta=-vector_inp(rnorm,p[ipart][isite]);
                    if(cositheta<sys.cosdelta) {
                        continue;
                    }
                    phi_i=0.5*(1.0-cos(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta));
                    for( jsite=0; jsite<sys.nsites; jsite++ ) {
                        cosjtheta=vector_inp(rnorm,p[jpart][jsite]);
                        if(cosjtheta<sys.cosdelta) {
                            continue;
                        }

                        phi_j=0.5*(1.0-cos(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta));
                        //Berthelot mixing rule is used for determining epsilon
                        epsmix = sqrt(sys.site[isite].eps * sys.site[jsite].eps);
                        //potential_energy_attP_d = 2.0*phi_i*phi_j*(r12inv - r6inv)*sqrt(sys.site[isite].eps * sys.site[jsite].eps);
                        potential_energy_attP_d = -2.0 * phi_i * phi_j * r6inv * epsmix;

                        //bewaar door welke patches de potentiele energie tot stand komt:
                        if( (isite==0) && (jsite==0)) {
                            psl->energyT = potential_energy_attP_d;
                            psl->minisite = isite;
                            psl->minjsite = jsite;
                        }
                        else {
                            psl->energyN = potential_energy_attP_d;
                            psl->minisite = isite;
                            psl->minjsite = jsite;
                        }
                        potential_energy_attP_d = 4.0*phi_i*phi_j*(r12inv - r6inv)*epsmix;

                        potential_energy_attP += potential_energy_attP_d;

                    }
                }
            }
        }
    }

    potential_energy_rep*=4.0;

    //potential_energy_attP*=4.0;//*sys.epsilonP; // dit wordt nu in de loop al gedaan

    //potential_energy_attC*=4.0*sys.epsilonC;

    //potential_energy = potential_energy_rep + potential_energy_attC + potential_energy_attP;
    potential_energy = potential_energy_rep + potential_energy_attP;

    return potential_energy;
}



void calculate_forces(Slice *psl) {

    Pts *psi,*psj;
    vector rij,p[sys.npart][sys.nsites],ui,uj,rnorm;
    vector rcrosspi,rcrosspj,piperpr,pjperpr,fangi,fangj;
    double r2,r6,r6inv,r12inv,r2inv,rinv,potential_energy=0,potential_energy_rep=0,r,cositheta,cosjtheta;
    double phi_i,phi_j,fmag,fmagP,fmagrinv,UmagP,epsmix;
    int ipart,jpart,isite,jsite;
    tensor rotmati,rotmatj;

    //initialize all forces to zero

    for( ipart=0; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        psi->f = nulvec;
        psi->t = nulvec;
        rotmati = getrotmatrix(psi->q);
        for( isite=0; isite<sys.nsites; isite++) {
            matrix_x_vector(rotmati,sys.site[isite].r,p[ipart][isite]);
        }
    }

    //calculate forces based on Lennard-Jones interaction for all particles
    //add torques
    for( ipart=1; ipart<sys.npart; ipart++ ) {
        psi =&psl->pts[ipart];
        for( jpart=0; jpart<ipart; jpart++ ) {
            psj = &psl->pts[jpart];
            vector_minus(psi->r,psj->r,rij);
            pbc(rij,sys.boxl);
            r2 = vector_inp(rij,rij);
            if(r2<4.0) {
                r2inv=1.0/r2;
                r6inv=r2inv*r2inv*r2inv;
                r12inv = r6inv*r6inv;

                if(r2<sys.sigmaLJsq) {
                    fmag = r2inv*(48.*r12inv - 24.*r6inv);
                    scalar_plustimes(rij,fmag,psi->f);
                    scalar_mintimes(rij,fmag,psj->f);
                }

                fmag = sys.epsilonC*r2inv*(48.*r12inv - 24.*r6inv);
                scalar_plustimes(rij,fmag,psi->f);
                scalar_mintimes(rij,fmag,psj->f);

                UmagP = 4.0*(r12inv - r6inv); // hoe fixen we dit met die epsilon???
                fmag = r2inv*(48.*r12inv - 24.*r6inv); // hoe fixen???

                r=sqrt(r2);
                //printf("distance particles %lf\n",r);
                rinv=1.0/r;
                scalar_times(rij,rinv,rnorm);
                //calculate angular part of force
                for( isite=0; isite<sys.nsites; isite++ ) {
                    cositheta=-vector_inp(rnorm,p[ipart][isite]);
                    if(cositheta<sys.cosdelta) {
                        continue;
                    }
                    //calculate phi for particle i
                    phi_i=0.5*(1.0-cos(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta));
                    //printf("phi_i %d %lf\n",ipart, phi_i);

                    for( jsite=0; jsite<sys.nsites; jsite++ ) {
                        cosjtheta=vector_inp(rnorm,p[jpart][jsite]);
                        if(cosjtheta<sys.cosdelta) {
                            continue;
                        }

                        epsmix = sqrt(sys.site[isite].eps*sys.site[jsite].eps);

                        //calculate phi for particle j
                        phi_j=0.5*(1.0-cos(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta));
                        //printf("phi_j %d %lf\n",jpart, phi_j);

                        fmagP = fmag*phi_i*phi_j*epsmix;
                        //printf("fmagP attraction %lf\n", fmagP);
                        scalar_plustimes(rij,fmagP,psi->f);
                        scalar_mintimes(rij,fmagP,psj->f);

                        vector_cross(rnorm,p[ipart][isite],rcrosspi);
                        vector_cross(rnorm,p[jpart][jsite],rcrosspj);

                        vector_cross(rnorm,rcrosspi,piperpr);
                        vector_cross(rnorm,rcrosspj,pjperpr);

                        //CALCULATE DEL U / DEL COSTHETA
                        fmag = UmagP*phi_j*HALFPI*sin(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta*epsmix;
                        //printf("fmag due to patch i %lf\n", fmag);
                        fmagrinv=fmag*rinv;
                        scalar_mintimes(piperpr,fmagrinv,psi->f);
                        scalar_plustimes(piperpr,fmagrinv,psj->f);
                        //scalar_plustimes(piperpr,fmagrinv,psi->f);
                        //scalar_mintimes(piperpr,fmagrinv,psj->f);

                        scalar_mintimes(rcrosspi,fmag,psi->t);

                        fmag = UmagP*phi_i*HALFPI*sin(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta*epsmix;
                        //printf("fmag due to patch j %lf\n", fmag);
                        fmagrinv=fmag*rinv;
                        scalar_plustimes(pjperpr,fmagrinv,psi->f);
                        scalar_mintimes(pjperpr,fmagrinv,psj->f);

                        scalar_plustimes(rcrosspj,fmag,psj->t);

                        //printf("torque %d %lf %lf %lf\n",ipart,psi->t.x,psi->t.y,psi->t.z);
                        //printf("torque %d %lf %lf %lf\n",jpart,psj->t.x,psj->t.y,psj->t.z);
                    }
                }
            }
        }
    }

    //exit(1);

    return ;
}



void propagate_bd(Slice *psl)  {

    int ipart,inter;
    Pts *psi;
    vector theta,f,t,u;
    quaternion qu1,qu2,qu3,qprime;
    tensor rotmat;
    quattensor Bmat;
    double dt,lambdaq;



    for( inter=0; inter<langevin.ninter; inter++) {
        calculate_forces(psl);

        for( ipart=0; ipart<sys.npart; ipart++) {
            psi = &psl->pts[ipart];

            //translational part
            //muT*F*dt*Beta
            scalar_times(psi->f,langevin.dtBeta,f);
            scalar_times(f,sys.mobilityT,f);
            vector_add(psi->r,f,psi->r);

            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityT,theta);
            vector_add(psi->r,theta,psi->r);

            pbc(psi->r,sys.boxl);

            //rotational part
            //get matrices for q(t)
            rotmat = getrotmatrix(psi->q);
            Bmat = getquatmatrix(psi->q);

            //Baalpha * (muR) * rotmatA * torque * delt = qu1
            scalar_times(psi->t,langevin.dtBeta,t);
            //do not forget, multiply with the inverse rotation matrix to convert to body-fixed torque
            matrixT_x_vector(rotmat,t,u);
            scalar_times(u,sys.mobilityR,u);
            quatmatrix_x_vec(Bmat,u,qu1);

            //Baalpha * (muR) * theta = qu2
            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityR,theta);
            quatmatrix_x_vec(Bmat,theta,qu2);

            //qprime = qu1+qu2+q(t) for the first time
            quat_add(qu1,qu2,qprime);

            quat_add(qprime,psi->q,qprime);

            //find lambdaq, I guess it is th smallest...check which one keeps q(t+delt)^2=1
            //woohoo it works for the smaller lambdaq, maybe also simply for the bigger...
            lambdaq=langrange_multiplier_quat(qprime, psi->q);
            if(lambdaq>1e20) {
                //langrange multiplier did not work, simply renormalize quaternion
                scdivide_quat(qprime,sqrt(quat_inp(qprime,qprime)),psi->q);
            }
            else {
                //lambdaq*q(t) = qu3
                sctimes_quat(psi->q,lambdaq,qu3);
                //q(t+delt) = qprime+qu3
                quat_add(qprime,qu3,psi->q);
            }

        }
    }

    return;
}





tensor getrotmatrix(quaternion q) {

    tensor mat;
    double q0,q1,q2,q3,q0sq,q1sq,q2sq,q3sq;

    q0=q.q0;
    q1=q.q1;
    q2=q.q2;
    q3=q.q3;

    q0sq=q0*q0;
    q1sq=q1*q1;
    q2sq=q2*q2;
    q3sq=q3*q3;

    mat.x.x = q0sq + q1sq - q2sq - q3sq ;  mat.x.y = 2.0*(q1*q2 - q0*q3)       ;  mat.x.z = 2.0*(q1*q3 + q0*q2)       ;
    mat.y.x = 2.0*(q1*q2 + q0*q3)       ;  mat.y.y = q0sq - q1sq + q2sq - q3sq ;  mat.y.z = 2.0*(q2*q3 - q0*q1)       ;
    mat.z.x = 2.0*(q1*q3 - q0*q2)       ;  mat.z.y = 2.0*(q2*q3 + q0*q1)       ;  mat.z.z = q0sq - q1sq - q2sq + q3sq ;

    return mat;
}




quattensor getquatmatrix(quaternion q) {

    quattensor qmat;
    double q0,q1,q2,q3;

    q0=0.5*q.q0;
    q1=0.5*q.q1;
    q2=0.5*q.q2;
    q3=0.5*q.q3;

    qmat.w.x = -q1 ;  qmat.w.y = -q2 ;  qmat.w.z = -q3 ;
    qmat.x.x =  q0 ;  qmat.x.y = -q3 ;  qmat.x.z =  q2 ;
    qmat.y.x =  q3 ;  qmat.y.y =  q0 ;  qmat.y.z = -q1 ;
    qmat.z.x = -q2 ;  qmat.z.y =  q1 ;  qmat.z.z =  q0 ;

    return qmat;
}




double langrange_multiplier_quat(quaternion qprime, quaternion q) {

    double lambdaq,a,b,c,det,root1,root2;

    //solving for equation 15 Ilie Briels den Otter
    //lambdaq^2 + 2*lambdaq*qprime*q + qprime*qprime == 1

    a = 1.0;
    b = 2.0*quat_inp(q,qprime);
    c = quat_inp(qprime,qprime)-1.0;

    det = b*b - 4.0*a*c;
    if(det<0.0) {
        //printf("Warning: determinant lower than 0, can not find roots\nWhat to do?\n");
        //printf("b value %lf c value %lf\n",b,c);
        return BIGNUM;
    }
    else {
        //hmmm which root to choose?
        root1=(-b + sqrt(det))/(2.0*a);
        root2=(-b + sqrt(det))/(2.0*a);
    }

    //go for minimum for now...makes more sense
    //if this does not work try just renormalizing it...

    if(fabs(root1)<fabs(root2)) {
        lambdaq=root1;
    }
    else {
        lambdaq=root2;
    }

    //printf("lambdaq %lf\n",lambdaq);
    return lambdaq;
}

void confoutput(Slice *psl) {

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




void confinput(Slice *psl) {

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








