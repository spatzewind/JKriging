package com.metzner.enrico.JKriging.helper;

public class MathHelper {
	
	public static double[][] adjunkte(double[][] mat, int i, int k) {
	    int mlen = mat.length-1;
	    double[][] res = new double[mlen][mlen];
	    for(int _k=0; _k<mlen; _k++) {
	        int kk=_k<k ? _k : _k+1;
	        for(int _i=0; _i<mlen; _i++) {
	            int ii=_i<i ? _i : _i+1;
	            res[_k][_i] = mat[kk][ii];
	        }
	    }
	    return res;
	}

	public static double determinante(double[][] mat) {
	    double res = 0d;
	    int mk = mat.length;
	    if(mk==1) return mat[0][0];
	    for(int _k=0; _k<mk; _k++) {
	        res += (_k%2==0 ? 1d : -1d) * mat[_k][0] * determinante(adjunkte(mat,0,_k));
	    }
	    return res;
	}

	public static double[][] inverse(double[][] mat) {
	    int mlen = mat.length;
	    double[][] res = new double[mlen][mlen];
	    double det = determinante(mat);
	    for(int _k=0; _k<mlen; _k++) {
	        for(int _i=0; _i<mlen; _i++) {
	            res[_k][_i] = ((_k+_i)%2==0 ? 1d : -1d) * determinante(adjunkte(mat,_k,_i)) / det;
	        }
	    }
	    return res;
	}

	public static double[][] diag(double[][] mat) {
	    return diag(mat, 1d);
	}
	public static double[][] diag(double[][] mat, double fac) {
	    double[][] res = new double[mat.length][mat[0].length];
	    for(int j=0; j<mat.length; j++) for(int i=0; i<mat[0].length; i++)
	        res[j][i] = (i==j ? mat[j][i]*fac : 0d);
	    return res;
	}

	public static double[][] transpose(double[][] mat) {
	    double[][] res = new double[mat[0].length][mat.length];
	    for(int j=0; j<mat.length; j++)
	        for(int i=0; i<mat[j].length; i++)
	            res[i][j] = mat[j][i];
	    return res;
	}

	public static double[][] matadd(double[][] matLeft, double[][] matRight) {
	    if(matLeft.length!=matRight.length || matLeft[0].length!=matRight[0].length)
	        throw new RuntimeException("matricies aren't addable, found m["+matLeft.length+"x"+matLeft[0].length+
	            "] * m["+matRight.length+"x"+matRight[0].length+"]");
	    int m = matLeft.length, n = matRight[0].length;
	    double[][] res = new double[m][n];
	    for(int j=0; j<m; j++) for(int i=0; i<n; i++)
	        res[j][i] = matLeft[j][i] + matRight[j][i];
	    return res;
	}
	public static double[][] matsub(double[][] matLeft, double[][] matRight) {
	    if(matLeft.length!=matRight.length || matLeft[0].length!=matRight[0].length)
	        throw new RuntimeException("matricies aren't subtractable, found m["+matLeft.length+"x"+matLeft[0].length+
	            "] * m["+matRight.length+"x"+matRight[0].length+"]");
	    int m = matLeft.length, n = matRight[0].length;
	    double[][] res = new double[m][n];
	    for(int j=0; j<m; j++) for(int i=0; i<n; i++)
	        res[j][i] = matLeft[j][i] - matRight[j][i];
	    return res;
	}
	public static double[] matmul(double[][] matLeft, double[] vecRight) {
	    int rlen = matLeft.length;
	    if(matLeft[0].length != vecRight.length)
	        throw new RuntimeException("matrix and vector aren't multipliable, found m["+
	                                   rlen+"x"+matLeft[0].length+"] * v["+vecRight.length+"]");
	    int klen = vecRight.length;
	    double[] res = new double[rlen];
	    for(int r=0; r<rlen; r++) {
	        res[r] = 0d;
	        for(int k=0; k<klen; k++)
	            res[r] += matLeft[r][k] * vecRight[k];
	    }
	    return res;
	}
	public static double[][] matmul(double[][] matLeft, double[][] matRight) {
	    int rlen = matLeft.length;
	    int clen = matRight[0].length;
	    if(matLeft[0].length != matRight.length)
	        throw new RuntimeException("matrices aren't multipliable, found ["+
	                                   rlen+"x"+matLeft[0].length+"] * ["+matRight.length+"x"+clen+"]");
	    int klen = matRight.length;
	    double[][] res = new double[rlen][clen];
	    for(int r=0; r<rlen; r++) for(int c=0; c<clen; c++) {
	        res[r][c] = 0d;
	        for(int k=0; k<klen; k++)
	            res[r][c] += matLeft[r][k] * matRight[k][c];
	    }
	    return res;
	}
}
