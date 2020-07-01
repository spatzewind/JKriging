package com.metzner.enrico.JKriging.data;

public class Constants {
	
	/* Mathematical constants */
	public static final double DEG2RAD  = Math.PI/180d;
	public static final double D_EPSLON = 1.0e-20d;

	/* Indexing constants */
	public static       int    FIRST_IDX        = 0;
	public static final int    JAVA_INDEXING    = 0;
	public static final int    FORTRAN_INDEXING = 1;
	
	public static void setFirstIndex(int _f_idx) { FIRST_IDX = _f_idx; }

	/* Versions */
	public static double GAMV_VERSION = 3.000d;
	public static String NETCDF_VERSION = "5.3.3";
}
