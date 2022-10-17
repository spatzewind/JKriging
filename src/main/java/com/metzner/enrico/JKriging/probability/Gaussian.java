package com.metzner.enrico.JKriging.probability;

public class Gaussian {
	
	private static final double NORM = Math.sqrt(0.5d/Math.PI);
	
	public static double pdf(double _x) {
		return Math.exp(-0.5d*_x*_x) * NORM;
	}

	public static double pdf(double _x, double _sigma) {
		return pdf(_x/_sigma)/Math.abs(_sigma);
	}
	
	public static double cdf(double _x) {
		return 0.5d * (1d + erf(_x));
	}
	
	public static double cdf(double _x, double _sigma) {
		return 0.5d * (1d + erf(_x/_sigma));
	}
	
	public static double cdfc(double _x) {
		return 0.5d * erfc(_x);
	}
	
	public static double cdfc(double _x, double _sigma) {
		return 0.5d * erfc(_x/_sigma);
	}
	
	public static double cdf_inv(double _x) {
		return ierf(2d*_x-1d);
	}
	
	public static double cdf_inv(double _x, double _sigma) {
		return _sigma * ierf(2d*_x-1d);
	}
	
	private static double erf(double _x) {
		if(_x<0d) return -erf(-_x);
		double x2 = -_x*_x;
		double _p = Math.sqrt(Math.PI);
		double res = 1d;
		if(_x>4d)
			return 1d - erfc(_x);
		double[] sum = new double[512];
		double s = _x;
		sum[0] = s;
		double fac = 0d;
		for(int i=1; i<sum.length; i++) {
			fac += 1d;
			s *= x2 / fac;
			sum[i] = s / (2*fac+1d);
		}
		res = 0d;
		for(int i=sum.length-1; i>=0; i--) {
			res += sum[i];
		}
		res *= 2d / _p;
		return res;
	}
	private static double erfc(double _x) {
		if(_x<0d) return -erfc(-_x);
		
		if(_x<4d) {
			return 1d - erf(_x);
		}
		
		double x22 = 2d*_x*_x;
		double erfc = 0.5d * 511d;
		for(int n=510; n>0; n--) {
			erfc = n / ( (n%2==0 ? x22 : 2d) + erfc );
		}
		erfc = Math.exp(-0.5d*x22) / ( x22 + erfc );
		return erfc / Math.sqrt(Math.PI);
	}
	
	private static double ierf(double _x) {
		if(_x<0d) return -ierf(-_x);
		double res = 0d;
		double e = erf(res);
		double tdspi = 2d/Math.sqrt(Math.PI);
		while(Math.abs(e-_x)>1.e-10d) {
			double der = tdspi*Math.exp(-res*res);
			double delta = (_x-e) / der;
			res += delta;
			e = erf(res);
		}
		return Math.sqrt(2d)*res;
	}
}
