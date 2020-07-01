package com.metzner.enrico.JKriging.probability;

import com.metzner.enrico.JKriging.helper.Rotation;

public class Covariance {

	public static double[] cova3(double x1, double y1, double z1, double x2, double y2, double z2,
			int ivarg, int[] nst, int maxnst, double[] c0, int[] it, double[] cc, double[] aa, int irot, int maxrot, double[][][] rotmat) {
		double cmax, cova;
//        c-----------------------------------------------------------------------
//        c
//        c                    Covariance Between Two Points
//        c                    *****************************
//        c
//        c This subroutine calculated the covariance associated with a variogram
//        c model specified by a nugget effect and nested varigoram structures.
//        c The anisotropy definition can be different for each nested structure.
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   x1,y1,z1         coordinates of first point
//        c   x2,y2,z2         coordinates of second point
//        c   nst(ivarg)       number of nested structures (maximum of 4)
//        c   ivarg            variogram number (set to 1 unless doing cokriging
//        c                       or indicator kriging)
//        c   MAXNST           size of variogram parameter arrays
//        c   c0(ivarg)        isotropic nugget constant
//        c   it(i)            type of each nested structure:
//        c                      1. spherical model of range a;
//        c                      2. exponential model of parameter a;
//        c                           i.e. practical range is 3a
//        c                      3. gaussian model of parameter a;
//        c                           i.e. practical range is a*sqrt(3)
//        c                      4. power model of power a (a must be gt. 0  and
//        c                           lt. 2).  if linear model, a=1,c=slope.
//        c                      5. hole effect model
//        c   cc(i)            multiplicative factor of each nested structure.
//        c                      (sill-c0) for spherical, exponential,and gaussian
//        c                      slope for linear model.
//        c   aa(i)            parameter "a" of each nested structure.
//        c   irot             index of the rotation matrix for the first nested 
//        c                    structure (the second nested structure will use
//        c                    irot+1, the third irot+2, and so on)
//        c   MAXROT           size of rotation matrix arrays
//        c   rotmat           rotation matrices
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   cmax             maximum covariance
//        c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
//        c
//        c
//        c
//        c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
//        c                      rotmat    computes rotation matrix for distance
//        c-----------------------------------------------------------------------
		//    parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
		double PMX = 999.0, EPSLON = 0.00001d;
		//    integer   nst(*),it(*)
		//    real      c0(*),cc(*),aa(*)
		//    real*8    rotmat(MAXROT,3,3),hsqd,sqdist
		double hsqd;
//        c
//        c Calculate the maximum covariance value (used for zero distances and
//        c for power model covariance):
//        c
		int istart = 1 + (ivarg-1)*maxnst,
			ist;
		cmax   = c0[ivarg-1];
		for(int is=0; is<nst[ivarg-1]; is++) {
			ist = istart + is - 1;
			if(it[ist]==4) {
				cmax += PMX;
			} else {
				cmax += cc[ist];
			}
		}
//        c
//        c Check for "zero" distance, return with cmax if so:
//        c
		hsqd = Rotation.sqdist(x1,y1,z1,x2,y2,z2,irot,rotmat);
		if(hsqd < EPSLON) {
			cova = cmax;
			return new double[] {cova, cmax};
		}
//        c
//        c Loop over all the structures:
//        c
		cova = 0d;
		int ir;
		double h, hr;
		for(int is=0; is<nst[ivarg-1]; is++) {
			ist = istart + is - 1;
//        c
//        c Compute the appropriate distance:
//        c
			if(ist!=0) {
				ir   = Math.min((irot+is-1),maxrot);
				hsqd = Rotation.sqdist(x1,y1,z1,x2,y2,z2,ir,rotmat);
			}
			h = Math.sqrt(hsqd);
//        c
//        c Spherical Variogram Model?
//        c
			if(it[ist]==1) {
				hr = h / aa[ist];
				if(hr<1d) cova += cc[ist]*(1.d-hr*(1.5d-.5d*hr*hr));
//        c
//        c Exponential Variogram Model?
//        c
			} else if(it[ist]==2) {
				cova += cc[ist]*Math.exp(-3d*h/aa[ist]);
//        c
//        c Gaussian Variogram Model?
//        c
			} else if(it[ist]==3) {
				cova += cc[ist]*Math.exp(-3.d*(h/aa[ist])*(h/aa[ist]));
//        c
//        c Power Variogram Model?
//        c
			} else if(it[ist]==4) {
				cova += cmax - cc[ist]*Math.pow(h,aa[ist]);
//        c
//        c Hole Effect Model?
//        c
			} else if(it[ist]==5) {
//        c                 d = 10.0 * aa(ist)
//        c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
				cova += cc[ist]*Math.cos(h/aa[ist]*Math.PI);
			}
		}
//        c
//        c Finished:
//        c
		return new double[] {cova, cmax};
	}
	
}
