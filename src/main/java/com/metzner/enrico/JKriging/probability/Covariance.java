package com.metzner.enrico.JKriging.probability;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.helper.MathHelper;

public class Covariance {

	public static final int VARIOGRAM_SPHERICAL      = 1;
	public static final int VARIOGRAM_EXPONENTIAL    = 2;
	public static final int VARIOGRAM_GAUSSIAN       = 3;
	public static final int VARIOGRAM_POWER          = 4;
	public static final int VARIOGRAM_HOLE_EFFECT    = 5;
	public static final int VARIOGRAM_LINEAR         = 6;
	public static final int VARIOGRAM_SINGLE_AXIS_HE = 7;

	/**
	 * Covariance between two points (2-D version)
	 * <p>
	 *     This function returns the covariance associated with a variogram model
	 *     that is specified by a nugget effect and possibly four different
	 *     nested varigoram structures.  The anisotropy definition can be
	 *     different for each of the nested structures (spherical, exponential,
	 *     gaussian, or power).
	 * </p>
	 * <p>
	 *     External references:
	 *         sqdist2D     computes anisotropic squared distance in 2D space (<a href="../helper.MathHelper">MathHelper</a>)
	 * </p>
	 * 
	 * @param x1            X coordinate for first point
	 * @param y1            Y coordinate for first point
	 * @param x2            X coordinate for second point
	 * @param y2            Y coordinate for second point
	 * @param nst           Number of nested structures (max. 4)
	 * @param c0            Nugget constant (isotropic)
	 * @param it            Type of each nested structure
	 * @param cc            Multiplicative factor of each nested structure
	 * @param aa            Parameter "a" of each nested structure
	 * @param rotmat        Rotation matrices
	 * @return
	 */
	public static double cova2(double x1, double y1, double x2, double y2,
			int nst, double c0, int[] it, double[] cc, double[] aa, double[][] rotmat) {
		double cmax, cova2;
//        c-----------------------------------------------------------------------
//        c
//        c              Covariance Between Two Points (2-D Version)
//        c              *******************************************
//        c
//        c This function returns the covariance associated with a variogram model
//        c that is specified by a nugget effect and possibly four different
//        c nested varigoram structures.  The anisotropy definition can be
//        c different for each of the nested structures (spherical, exponential,
//        c gaussian, or power).
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   x1,y1            Coordinates of first point
//        c   x2,y2            Coordinates of second point
//        c   nst              Number of nested structures (max. 4).
//        c   c0               Nugget constant (isotropic).
//        c   PMX              Maximum variogram value needed for kriging when
//        c                      using power model.  A unique value of PMX is
//        c                      used for all nested structures which use the
//        c                      power model.  therefore, PMX should be chosen
//        c                      large enough to account for the largest single
//        c                      structure which uses the power model.
//        c   cc(nst)          Multiplicative factor of each nested structure.
//        c   aa(nst)          Parameter "a" of each nested structure.
//        c   it(nst)          Type of each nested structure:
//        c                      1. spherical model of range a;
//        c                      2. exponential model of parameter a;
//        c                           i.e. practical range is 3a
//        c                      3. gaussian model of parameter a;
//        c                           i.e. practical range is a*sqrt(3)
//        c                      4. power model of power a (a must be gt. 0  and
//        c                           lt. 2).  if linear model, a=1,c=slope.
//        c   ang(nst)         Azimuth angle for the principal direction of
//        c                      continuity (measured clockwise in degrees from Y)
//        c   anis(nst)        Anisotropy (radius in minor direction at 90 degrees
//        c                      from "ang" divided by the principal radius in 
//        c                      direction "ang")
//        c   first            A logical variable which is set to true if the
//        c                      direction specifications have changed - causes
//        c                      the rotation matrices to be recomputed.
//        c
//        c
//        c
//        c OUTPUT VARIABLES: returns "cova2" the covariance obtained from the
//        c                   variogram model.
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		double PMX = 999.0d, EPSLON_S = 0.0000001d;
		double hsqd;
		cmax   = c0;
		for(int is=0; is<nst; is++) {
			if(it[is]==4) {
				cmax += PMX;
			} else {
				cmax += cc[is];
			}
		}
//        c
//        c Check for very small distance:
//        c
		hsqd = MathHelper.sqdist2D(x1, y1, x2, y2, new double[] {1d,0d,0d,1d});
		if(hsqd < EPSLON_S) {
			return cmax;
		}
//        c
//        c Non-zero distance, loop over all the structures:
//        c
		cova2 = 0d;
		for(int is=0; is<nst; is++) {
//        c
//        c Compute the appropriate structural distance:
//        c
			hsqd = MathHelper.sqdist2D(x1, y1, x2, y2, rotmat[is]);
			double h = Math.sqrt(hsqd);
			cova2 += covariance(it[is], PMX, cc[is], h, aa[is]);
		}
		return cova2;
	}

	/**
	 * Covariance bewteen two points (3-D version)
	 * <p>
	 *     This subroutine calculated the covariance associated with a variogram
	 *     model specified by a nugget effect and nested varigoram structures.
	 *     The anisotropy definition can be different for each nested structure.
	 * </p>
	 * 
	 * @param x1            X coordinate for first point
	 * @param y1            Y coordinate for first point
	 * @param z1            Z coordinate for first point
	 * @param x2            X coordinate for second point
	 * @param y2            Y coordinate for second point
	 * @param z2            Z coordinate for second point
	 * @param ivarg         Variogram number (set to 1 unless doing cokriging or indicator kriging)
	 * @param nst           Number of nested structures (maximum of 4)
	 * @param maxnst        Size of variogram parameter arrays
	 * @param c0            Isotropic nugget constant
	 * @param it            Type of each nested structure
	 * @param cc            Multiplicative factor of each nested structure.
	 *                      (sill-c0) for spherical, exponential and gaussian
	 *                      slope for linear model
	 * @param aa            Parameter "a" of each nested structure
	 * @param irot          Index of the rotation matrix for the first nested structure
	 *                      (the second nested structure will use irot+1, the third irot+2, and so on)
	 * @param maxrot        Size of rotation matrix arrays
	 * @param rotmat        Rotation matricies
	 * @return
	 */
	public static double[] cova3(double x1, double y1, double z1, double x2, double y2, double z2,
			int ivarg, int[] nst, int maxnst, double[] c0, int[] it, boolean[] sahe, double[] cc, double[] aa, int irot, int maxrot, double[][][] rotmat) {
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
		double PMX = 999.0d, EPSLON = 0.00001d;
		//    integer   nst(*),it(*)
		//    real      c0(*),cc(*),aa(*)
		//    real*8    rotmat(MAXROT,3,3),hsqd,sqdist
		double hsqd;
//        c
//        c Calculate the maximum covariance value (used for zero distances and
//        c for power model covariance):
//        c
		int istart = 1 + (ivarg-1)*maxnst, ist;
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
		hsqd = MathHelper.sqdist3D(x1,y1,z1,x2,y2,z2,rotmat[irot-1]);
		if(hsqd < EPSLON) {
			cova = cmax;
			return new double[] {cova, cmax};
		}
//        c
//        c Loop over all the structures:
//        c
		cova = 0d;
		int ir = Math.min(irot-1, maxrot);
		double h;
		for(int is=0; is<nst[ivarg-1]; is++) {
			ist = istart + is - 1;
//        c
//        c Compute the appropriate distance:
//        c
			if(ist!=0) {
				ir   = Math.min((irot+is-1),maxrot);
				hsqd = MathHelper.sqdist3D(x1,y1,z1,x2,y2,z2,rotmat[ir-1]);
			}
			if(sahe[ir]) {
				double[] hh = MathHelper.matmul(rotmat[ir-1], new double[] {x2-x1,y2-y1,z2-z1});
				h = Math.sqrt(hh[0]*hh[0]+hh[1]*hh[1]);
				cova += covariance(it[ist], cmax, cc[ist], h, aa[ist]) * covariance(VARIOGRAM_SINGLE_AXIS_HE, cmax, 1d, hh[2]*Constants.SAHE_FAC, 1d);
			} else {
				h = Math.sqrt(hsqd);
				cova += covariance(it[ist], cmax, cc[ist], h, aa[ist]);
			}
		}
//        c
//        c Finished:
//        c
		return new double[] {cova, cmax};
	}
	
	public static double covariance(int cova_type, double cc_cova_maximum, double cc_cova_coeff, double h_lag, double aa_range) {
		double hr = h_lag / aa_range;
		switch(cova_type) {
			case VARIOGRAM_SPHERICAL:
				return cc_cova_coeff * (hr<1d ? 1d-hr*(1.5d-0.5d*hr*hr) : 0d);
			case VARIOGRAM_EXPONENTIAL:
				return cc_cova_coeff * Math.exp(-3d*hr);
			case VARIOGRAM_GAUSSIAN:
				return cc_cova_coeff * Math.exp(-3d*hr*hr);
			case VARIOGRAM_POWER:
				return cc_cova_maximum - cc_cova_coeff * Math.pow(h_lag, aa_range);
			case VARIOGRAM_HOLE_EFFECT:
				return cc_cova_coeff * Math.exp(-0.3d*hr) * Math.cos(hr*Math.PI);
			case VARIOGRAM_LINEAR:
				return cc_cova_coeff * (hr<1d ? 1d-hr : 0d);
			case VARIOGRAM_SINGLE_AXIS_HE:
				return cc_cova_coeff * Math.cos(hr*Math.PI);
			default: return 0d;
		}
	}

	public static double[] covarianceDerivative(int cova_type, double cc_cova_maximum, double cc_cova_coeff, double h_lag, double aa_range) {
		double hr = h_lag / aa_range;
		double hr_d = -h_lag / (aa_range * aa_range);
		switch(cova_type) {
			case VARIOGRAM_SPHERICAL:
				return new double[]{
							(hr<1d ? 1d-hr*(1.5d-0.5d*hr*hr) : 0d),
							cc_cova_coeff * (hr<1d ? 1.5d*hr*hr-1.5d : 0d) * hr_d
						};
			case VARIOGRAM_EXPONENTIAL:
				return new double[] {
							Math.exp(-3d*hr),
							cc_cova_coeff * Math.exp(-3d*hr) * (-3d) * hr_d
						};
			case VARIOGRAM_GAUSSIAN:
				return new double[] {
							Math.exp(-3d*hr*hr),
							cc_cova_coeff * Math.exp(-3d*hr*hr) * (-6d*hr) * hr_d
						};
			case VARIOGRAM_POWER:
				return new double[] {
							-Math.pow(h_lag, aa_range),
							-cc_cova_coeff * Math.pow(h_lag, aa_range) * Math.log(h_lag)
						};
			case VARIOGRAM_HOLE_EFFECT:
				return new double[] {
							Math.exp(-0.3d*hr) * Math.cos(hr*Math.PI),
							-cc_cova_coeff * Math.exp(-0.3d*hr) * (0.3d * Math.cos(hr*Math.PI) + Math.PI*Math.sin(hr*Math.PI)) * hr_d
						};
			case VARIOGRAM_LINEAR:
				return new double[] {
							(hr<1d ? 1d-hr : 0d),
							cc_cova_coeff * (hr<1d ? -hr_d : 0d)
						};
			case VARIOGRAM_SINGLE_AXIS_HE:
				return new double[] {
							Math.cos(hr*Math.PI),
							-cc_cova_coeff * Math.PI * Math.sin(hr*Math.PI) * hr_d
						};
			default: return new double[] {0d, 0d};
		}
	}
}
