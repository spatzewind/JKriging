package com.metzner.enrico.JKriging.helper;

import com.metzner.enrico.JKriging.data.Constants;

public class Rotation {

	/**
	 * Sets up an Anisotropic Rotation Matrix
	 * <p>
	 *     Sets up the matrix to transform cartesian coordinates to coordinates
	 *     accounting for angles and anisotropy (see manual for a detailed
	 *     definition):
	 * </p>
	 * <p>
	 *     NO EXTERNAL REFERENCES
	 * </p>
	 * <p>
	 *     ALgorithm:
	 *     
	 *     Converts the semi-axis lengths into anisotropy informations:
	 *              anis1   Ratio of minor axis length to major axis length
	 *              anis2   Ratio of second minor axis length to major axis length
	 *     
	 *     Converts the input angles to three angles which make more mathematical sense:
	 *              alpha   angle between the major axis of anisotropy and the
	 *                      E-W axis. Note: Counter clockwise is positive.
	 *              beta    angle between major axis and the horizontal plane.
	 *                      (The dip of the ellipsoid measured positive down)
	 *              theta   Angle of rotation of minor axis about the major axis
	 *                      of the ellipsoid with positive clockwise in the
	 *                      direction of the major axis (from inside to outside)
	 *     
	 *     After that following scheme is used:
	 *     
	 *                 -> * Scale(1,anis1,anis2)     -> * Rot(x-axis,theta)  -> * Rot(y-axis,beta)  _-> * Rot(z-axis,alpha) ->
	 *         Sphere                                                                                                          Ellipse
	 *                 <- * Scale(1,1/anis1,1/anis2) <- * Rot(x-axis,-theta) <- * Rot(y-axis,-beta) <- * Rot(z-axis,-alpha) <-
	 *     
	 *     The final transformation matrix got the following look:
	 *     
	 *             [ 1    0       0    ]   [ 1     0           0       ]   [ cos(-beta) 0 sin(-beta) ]   [ cos(-alpha) sin(-alpha) 0 ]
	 *         M = [ 0 1/anis1    0    ] * [ 0 cos(-theta) sin(-theta) ] * [     0      1     0      ] * [ sin(-alpha) cos(-alpha) 0 ]
	 *             [ 0    0    1/anis2 ]   [ 0 sin(-theta) cos(-theta) ]   [ sin(-beta) 0 cos(-beta) ]   [     0           0       1 ]
	 *             
	 *             [                            ]
	 *         M = [ ()/anis1 ()/anis1 ()/anis1 ]
	 *             [ ()/anis2 ()/anis2 ()/anis2 ]
	 * </p>
	 * @param ang1      Azimuth angle for principal direction
	 * @param ang2      Dip angle for principal direction
	 * @param ang3      Third rotation angle
	 * @param anis1     First anisotropy ratio
	 * @param anis2     Second anisotropy ratio
	 * @param ind       matrix indicator to initialize
	 * @param maxrot    maximum number of rotation matrices dimensioned
	 * @param rotmat    rotation matrices
	 * @return
	 */
	public static double[][][] setrot(double ang1, double ang2, double ang3, double anis1, double anis2, int ind, int maxrot, double[][][] rotmat) {
//        c-----------------------------------------------------------------------
//        c
//        c              Sets up an Anisotropic Rotation Matrix
//        c              **************************************
//        c
//        c Sets up the matrix to transform cartesian coordinates to coordinates
//        c accounting for angles and anisotropy (see manual for a detailed
//        c definition):
//        c
//        c
//        c INPUT PARAMETERS:
//        c
//        c   ang1             Azimuth angle for principal direction
//        c   ang2             Dip angle for principal direction
//        c   ang3             Third rotation angle
//        c   anis1            First anisotropy ratio
//        c   anis2            Second anisotropy ratio
//        c   ind              matrix indicator to initialize
//        c   MAXROT           maximum number of rotation matrices dimensioned
//        c   rotmat           rotation matrices
//        c
//        c
//        c NO EXTERNAL REFERENCES
//        c
//        c
//        c-----------------------------------------------------------------------
		//add own extra
		if(ind>maxrot) throw new ArrayIndexOutOfBoundsException(ind+" of maxrot="+maxrot);
		//end own addition
		double afac1,afac2,sina,sinb,sint,cosa,cosb,cost;
		double alpha,beta,theta;
//        c
//        c Converts the input angles to three angles which make more
//        c  mathematical sense:
//        c
//        c         alpha   angle between the major axis of anisotropy and the
//        c                 E-W axis. Note: Counter clockwise is positive.
//        c         beta    angle between major axis and the horizontal plane.
//        c                 (The dip of the ellipsoid measured positive down)
//        c         theta   Angle of rotation of minor axis about the major axis
//        c                 of the ellipsoid.
//        c
		if(ang1>=0d && ang1<270d) {
			alpha = (90d  - ang1) * Constants.DEG2RAD;
		} else {
			alpha = (450d - ang1) * Constants.DEG2RAD;
		}
		beta  = -1d * ang2 * Constants.DEG2RAD;
		theta = ang3 * Constants.DEG2RAD;
//        c
//        c Get the required sines and cosines:
//        c
		sina  = Math.sin(alpha);
		sinb  = Math.sin(beta);
		sint  = Math.sin(theta);
		cosa  = Math.cos(alpha);
		cosb  = Math.cos(beta);
		cost  = Math.cos(theta);
//        c
//        c Construct the rotation matrix in the required memory:
//        c
		afac1 = 1d / Math.max(anis1,Constants.D_EPSLON);
		afac2 = 1d / Math.max(anis2,Constants.D_EPSLON);
		rotmat[ind-1][0][0] =       (cosb * cosa);
		rotmat[ind-1][0][1] =       (cosb * sina);
		rotmat[ind-1][0][2] =       (-sinb);
		rotmat[ind-1][1][0] = afac1*(-cost*sina + sint*sinb*cosa);
		rotmat[ind-1][1][1] = afac1*(cost*cosa + sint*sinb*sina);
		rotmat[ind-1][1][2] = afac1*( sint * cosb);
		rotmat[ind-1][2][0] = afac2*(sint*sina + cost*sinb*cosa);
		rotmat[ind-1][2][1] = afac2*(-sint*cosa + cost*sinb*sina);
		rotmat[ind-1][2][2] = afac2*(cost * cosb);
//        c
//        c Return to calling program:
//        c
		return rotmat;
	}
	
	/**
	 * Squared Anisotropic Distance Calculation Given Matrix Indicator
	 * <p>
	 *     This routine calculates the anisotropic distance between two points
	 *     given the coordinates of each point and a definition of the
	 *     anisotropy.
	 * </p>
	 * <p>
	 *     NO EXTERNAL REFERENCE
	 * </p>
	 * @param x1	    x-coordinate of first point
	 * @param y1        y-coordinate of first point
	 * @param z1        z-coordinate of first point
	 * @param x2	    x-coordinate of second point
	 * @param y2	    y-coordinate of second point
	 * @param z2	    z-coordinate of second point
	 * @param rotmat    the rotation matrix
	 * @return          The squared distance accounting for the anisotropy and the rotation of coordinates (if any).
	 */
	public static double sqdist(double x1, double y1, double z1, double x2, double y2, double z2,
			double[][] rotmat) {
//        c-----------------------------------------------------------------------
//        c
//        c    Squared Anisotropic Distance Calculation Given Matrix Indicator
//        c    ***************************************************************
//        c
//        c This routine calculates the anisotropic distance between two points
//        c  given the coordinates of each point and a definition of the
//        c  anisotropy.
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   x1,y1,z1         Coordinates of first point
//        c   x2,y2,z2         Coordinates of second point
//        c   ind              The rotation matrix to use
//        c   MAXROT           The maximum number of rotation matrices dimensioned
//        c   rotmat           The rotation matrices
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   sqdist           The squared distance accounting for the anisotropy
//        c                      and the rotation of coordinates (if any).
//        c
//        c
//        c NO EXTERNAL REFERENCES
//        c
//        c
//        c-----------------------------------------------------------------------
		//    real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
		double cont,dx,dy,dz;
//        c
//        c Compute component distance vectors and the squared distance:
//        c
		dx = (x1 - x2);
		dy = (y1 - y2);
		dz = (z1 - z2);
		double sqdist = 0d;
		for(int i=0; i<3; i++) {
			cont   = rotmat[i][0] * dx
					+rotmat[i][1] * dy
					+rotmat[i][2] * dz;
			sqdist += cont * cont;
		}
		return sqdist;
	}
	
	/**
	 * Squared Anisotropic Distance Calculation Given Matrix Indicator
	 * <p>
	 *     This routine calculates the anisotropic distance between two points
	 *     given the coordinates of each point and a definition of the
	 *     anisotropy in k dimensions.
	 * </p>
	 * <p>
	 *     NO EXTERNAL REFERENCE
	 * </p>
	 * @param pos1      coordinate vector of first point
	 * @param pos2      coordinate vector of second point
	 * @param rotmat    the rotation matrix
	 * @return          The squared distance accounting for the anisotropy and the rotation of coordinates (if any).
	 */
	public static double sqdist(double[] pos1, double[] pos2, double[][] rotmat) {
		double sqdist = 0d;
		for(int j=0; j<rotmat.length; j++) {
			double cont = 0d;
			for(int i=0; i<rotmat[j].length; i++)
				cont += rotmat[j][i] * (pos1[i] - pos2[i]);
			sqdist += cont * cont;
		}
		return sqdist;
	}
}
