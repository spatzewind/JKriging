package com.metzner.enrico.JKriging.helper;

import com.metzner.enrico.JKriging.data.Constants;

public class Rotation {

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
	
	public static double sqdist(double x1, double y1, double z1, double x2, double y2, double z2,
			int ind, double[][][] rotmat) {
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
			cont   = rotmat[ind-1][i][0] * dx
					+rotmat[ind-1][i][1] * dy
					+rotmat[ind-1][i][2] * dz;
			sqdist += cont * cont;
		}
		return sqdist;
	}
}