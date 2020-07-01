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
		if(_x>4d) {
			double erfc = 0.5d * 511d / _x;
			for(int n=510; n>0; n--) {
				erfc = n / ( (n%2==1 ? 2*_x : _x) + erfc );
			}
			erfc = Math.exp(x2) / ( _x + erfc );
			return 1d - erfc / Math.sqrt(Math.PI);
		} else {
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
		}
		return res;
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
