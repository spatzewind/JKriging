package com.metzner.enrico.JKriging.helper;

import com.metzner.enrico.JKriging.data.Constants;

public class Interpolation {

	public static double linint(double in_s, double in_e, double out_s, double out_e, double val) {
		if(Math.abs(in_e-in_s)<Constants.D_EPSLON) return 0.5d*(out_s+out_e);
		return out_s + (out_e-out_s) * (val-in_s)/(in_e-in_s);
	}
	
	public static double powint(double in_s, double in_e, double out_s, double out_e, double val, double expo) {
		if(Math.abs(in_e-in_s)<Constants.D_EPSLON) return 0.5d*(out_s+out_e);
		return out_s + (out_e-out_s) * Math.pow((val-in_s)/(in_e-in_s), expo);
	}
	
	public static double hypint(double in_s, double in_e, double out_s, double val, double expo) {
		double lambda = Math.pow(out_s,expo)*(in_e-in_s);
		return Math.pow(lambda/(in_e-val),1d/expo);
	}
}
