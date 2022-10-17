package com.metzner.enrico.JKriging.data;

import java.util.Random;

public class Constants {
	
	/* Project constants */
	public static final String NAME = "JKriging";
	public static final String VERSION = "0.0.4";
	
	
	
	/* Mathematical constants */
	public static final double RAD2DEG  = 180d/Math.PI;
	public static final double DEG2RAD  = Math.PI/180d;
	public static final double D_EPSLON = 1.0e-20d;
	public static final double SAHE_FAC = 4.5d;
	public static final Random RANDOM = new Random();



	/* Indexing constants */
	public static       int    FIRST_IDX        = 0;
	public static final int    JAVA_INDEXING    = 0;
	public static final int    FORTRAN_INDEXING = 1;
	
	public static void setFirstIndex(int _f_idx) { FIRST_IDX = _f_idx; }




	/* Data */
	public static       float  FILL_VALUE_F = 1.0e+20f;
	public static       double FILL_VALUE_D = 1.0e+20d;

	public static void setFloatFillValue(float _fill_value) { FILL_VALUE_F = _fill_value; }
	public static void setDoubleFillValue(double _fill_value) { FILL_VALUE_D = _fill_value; }




	/* Options */
	public static final boolean BNO  = false;
	public static final boolean BYES = true;
	public static final int     INO  = 0;
	public static final int     IYES = 1;




	/* Versions */
	public static double GAMV_VERSION = 3.000d;
	public static String NETCDF_VERSION = "5.5.2";

}
