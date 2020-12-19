package com.metzner.enrico.JKriging.helper;

import java.util.List;

import com.metzner.enrico.JKriging.data.Constants;

public class MathHelper {
	
	// ************************************************************ //
	// ** COMPLEX ALGEBRA                                        ** //
	// ************************************************************ //
	public static double[] cplxsqrt(double[] cplx) {
		double rr = 0.5d * (cplx[0] + Math.sqrt(cplx[0]*cplx[0]+cplx[1]*cplx[1]));
		double ii = rr - cplx[0];
		return new double[] {Math.sqrt(rr),Math.sqrt(ii)};
	}
	
	// ************************************************************ //
	// ** MATRIX CALCULATIONS                                    ** //
	// ************************************************************ //
	public static double[][] identity(int dimension) {
		double[][] id = new double[dimension][dimension];
		for(int j=0; j<dimension; j++)
			for(int i=0; i<dimension; i++)
				id[j][i] = (j==i ? 1d : 0d);
		return id;
	}
	public static boolean isSymmetric(double[][] mat) {
		if(mat.length!=mat[0].length) return false;
		for(int j=1; j<mat.length; j++)
			for(int i=0; i<j; i++)
				if(mat[j][i]!=mat[i][j])
					return false;
		return true;
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
	public static double[][] matmul(double[][] mat, double factor) {
		double[][] res = new double[mat.length][mat[0].length];
		for(int mj=0; mj<mat.length; mj++)
			for(int mi=0; mi<mat[mj].length; mi++)
				res[mj][mi] = mat[mj][mi] * factor;
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
	public static double[][] transpose(double[] vec) {
	    double[][] res = new double[1][vec.length];
        for(int i=0; i<vec.length; i++)
            res[0][i] = vec[i];
	    return res;
	}
	public static double trace(double[][] mat) {
		if(mat.length!=mat[0].length) {
			System.err.println("No square matrix! return NaN as trace.");
			DataHelper.printStackTrace(System.err);
			return Double.NaN;
		}
		double t = 0d;
		for(int m=0; m<mat.length; m++) t += mat[m][m];
		return t;
	}
	public static double[][] transpose(double[][] mat) {
	    double[][] res = new double[mat[0].length][mat.length];
	    for(int j=0; j<mat.length; j++)
	        for(int i=0; i<mat[j].length; i++)
	            res[i][j] = mat[j][i];
	    return res;
	}
	public static double[][] eigenvalues(double[][] mat) {
		int dim = mat.length;
		if(mat[0].length!=dim) {
			System.err.println("Cannot calculate eigenvalues of a non-square matrix!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		if(dim==1) {
			//case 1x1 matrix -> constant equation -> constant
			return new double[][] { {mat[0][0],0d} };
		}
		if(dim==2) {
			//case 2x2 matrix -> quadratic equation -> p-q-formula
			double p = mat[0][0]+mat[1][1];
			double q = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
			double[][] res = { {0.5d*p,0d}, {0.5d*p,0d} };
			double r = Math.sqrt(Math.abs(0.25d*p*p-q));
			if(p*p<4*q) { res[0][1] = -r; res[1][1] = r; }
			else { res[0][0] -= r; res[1][0] += r; }
			return res;
		}
		if(dim==3) {
			//case 3x3 matrix -> cubic equation -> tranformation to use cosine
			//System.out.println("[MATRIX EIGENVALUE] case of 3x3 matrix");
			double a = -trace(mat);
			double b = 0.5d*trace(mat)*trace(mat)-0.5d*trace(matmul(mat,mat));
			double c = -determinante(mat);
			//System.out.println("[MATRIX EIGENVALUE] characteristic polynomial: 0=x³+"+a+"x²+"+b+"x+"+c);
			double p = b-a*a/3d;
			double q = a*(2*a*a/27d - b/3d) + c;
			//System.out.println("[MATRIX EIGENVALUE] reduced polynomial: 0=x³+"+p+"x+"+q);
			double D = q*q/4d + p*p*p/27d; //(27*c*c + 4*a*a*a*c-18*a*b*c+4*b*b*b) / 108d;
			//System.out.println("[MATRIX EIGENVALUE] discriminant: "+D);
			double[][] ev = new double[3][2];
			if(D>0d) {
				double u3 =  Math.sqrt(D)-0.5d*q;
				double v3 = -Math.sqrt(D)-0.5d*q;
				double u = Math.pow(Math.abs(u3), 1d/3d)*(u3<0d ? -1 : 1);
				double v = Math.pow(Math.abs(v3), 1d/3d)*(v3<0d ? -1 : 1);
				ev[0][0] =        u+v;  ev[0][1] = 0d;
				ev[1][0] = -0.5d*(u+v); ev[1][1] = 0.5d*(u-v)*Math.sqrt(3d);
				ev[2][0] = -0.5d*(u+v); ev[2][1] = 0.5d*(v-u)*Math.sqrt(3d);
				if(u+v<0d) {
					double t = ev[0][0]; ev[0][0] = ev[2][0]; ev[2][0] = t;
					t = ev[0][1]; ev[0][1] = ev[2][1]; ev[2][1] = t;
				}
			} else if(D<0d) {
				double phase = Math.PI/3d;
				double alpha = Math.acos(-0.5d*q*Math.sqrt(-27d/(p*p*p)))/3d;
				double fac = Math.sqrt(-4*p/3d);
				ev[0][0] = -fac*Math.cos(alpha-phase); ev[0][1] = 0d;
				ev[1][0] =  fac*Math.cos(alpha);       ev[1][1] = 0d;
				ev[2][0] = -fac*Math.cos(alpha+phase); ev[2][1] = 0d;
				if(ev[0][0]>ev[1][0] && ev[0][0]>ev[2][0]) {
					double t = ev[0][0]; ev[0][0] = ev[2][0]; ev[2][0] = t;
				} else if(ev[1][0]>ev[0][0] && ev[1][0]>ev[2][0]) {
					double t = ev[1][0]; ev[1][0] = ev[2][0]; ev[2][0] = t;
				}
				if(ev[0][0]>ev[1][0]) {
					double t = ev[0][0]; ev[0][0] = ev[1][0]; ev[1][0] = t;
				}
			} else {
				if(p==0d) {
					ev[0][0] = 0d; ev[0][1] = 0d;
					ev[1][0] = 0d; ev[1][1] = 0d;
					ev[2][0] = 0d; ev[2][1] = 0d;
				} else {
					ev[0][0] =    3d*q/p; ev[0][1] = 0d;
					ev[1][0] = -1.5d*q/p; ev[1][1] = 0d;
					ev[2][0] = -1.5d*q/p; ev[2][1] = 0d;
					if((q<0d && p>0d) || (p<0d && q>0d)) {
						double t = ev[0][0]; ev[0][0] = ev[2][0]; ev[2][0] = t;
					}
				}
			}
			ev[0][0] -= a/3d;
			ev[1][0] -= a/3d;
			ev[2][0] -= a/3d;
			return ev;
		}
		//QR-algorithm
		double[][] mA = matmul(mat,1d);
		//double shiftSum = 0d;
		for(int iter=0; iter<10*dim*dim; iter++) {
//			double[][] ev22 = eigenvalues(new double[][] {{mA[dim-2][dim-2],mA[dim-2][dim-1]},{mA[dim-1][dim-2],mA[dim-1][dim-1]}});
//			double lowrgt = mA[dim-1][dim-1];
//			double[] ev22dist = {(ev22[0][0]-lowrgt)*(ev22[0][0]-lowrgt)+ev22[0][1]*ev22[0][1],
//								 (ev22[1][0]-lowrgt)*(ev22[1][0]-lowrgt)+ev22[1][1]*ev22[1][1]};
//			int ev22idx = (ev22dist[0]<ev22dist[1] ? 0 : 1);
//			double shiftvalue = (ev22[ev22idx][1]!=0d ? lowrgt : ev22[ev22idx][0]);
			double shiftvalue = mA[dim-1][dim-1];
			double[][] mAk = matsub(mA,matmul(identity(dim),shiftvalue));
			//shiftSum += shiftvalue;
			//step b: decompose "mat" to unitary matrix Q and upper triangular matrix R via Householder-matrices
			double[][] mR = matmul(mAk, 1d);
			double[][] mQ = identity(dim);
			for(int k=1; k<=dim-2; k++) {
				double[][] mP = getHouseholderMatrix_column(mR, k, 0);
				mR = matmul(mP,mR);
				mQ = matmul(mQ,mP);
			}
			mA = matadd(matmul(mR,mQ),matmul(identity(dim),shiftvalue));
		}
		double[][] evs = new double[dim][2];
		double tr = Math.pow(determinante(mat), 1d/dim) * 1.0e-12d;
		for(int d=0; d<dim; d++) {
			if(dim-d==1) {
				evs[d][0] = mA[d][d];// - shiftSum;
				evs[d][1] = 0d;
			} else {
				if(Math.abs(mA[d+1][d])>tr && Math.abs(mA[d+1][d]+mA[d][d+1])<tr && Math.abs(mA[d][d]-mA[d+1][d+1])<tr) {
					evs[d][0] = mA[d][d];// - shiftSum;
					evs[d][1] = mA[d+1][d];
					evs[d+1][0] = mA[d+1][d+1];// - shiftSum;
					evs[d+1][1] = mA[d][d+1];
					d++;
				} else {
					evs[d][0] = mA[d][d];// - shiftSum;
					evs[d][1] = 0d;
				}
			}
		}
		return evs;
	}
	/**
	 * Apply QR-algorithm on matrix to extract eigenvalues and eigenvectors
	 * @param mat real matrix (no symmetry required)
	 * @return matrix (as double[][]) with eigenvectors columnwise and corresponding eigenvalues in the last row
	 */
	public static double[][] eigenQR(double[][] mat) {
		double eps = 1.0e-32d;
		int[] eigenvaluesID = new int[mat.length];
		double[][] eigenvalues = new double[mat.length][2];
		//step 1: transform matrix via Householder matrices to Hessenberg form
		double[][][] hessenberg = transformToSimilarHessenberg(mat);
		double[][] copymat = hessenberg[0];
		double[][] eigenvectors = hessenberg[1];
		double[][] matQ = new double[mat.length][mat.length];
		for(int j=0; j<mat.length; j++) {
			eigenvalues[j][0] = hessenberg[0][j][j];
			if(j==0) {
				eigenvalues[j][1] = hessenberg[0][j+1][j];
			} else if(j<mat.length-1) {
				eigenvalues[j][1] = 0.5d*(hessenberg[0][j+1][j]-hessenberg[0][j][j-1]);
			} else {
				eigenvalues[j][1] = -hessenberg[0][j][j-1];
			}
		}
		int n = mat.length-1;
		int startIter = -1, iter = 0;
		double[][] implShiftMat;
		while(n>=0) {
			double[][] subMat = new double[n+1][n+1];
			for(int j=0; j<=n; j++)
				for(int i=0; i<=n; i++)
					subMat[j][i] = copymat[j][i];
			if(n>0) {
				double u = copymat[n-1][n-1]+copymat[n][n];
				double v = copymat[n-1][n-1]*copymat[n][n]-copymat[n-1][n]*copymat[n][n-1];
				implShiftMat = matmul(subMat,subMat);
				implShiftMat = matsub(implShiftMat,matmul(subMat,u));
				implShiftMat = matadd(implShiftMat,matmul(identity(n+1),v));
			} else {
				implShiftMat = matsub(subMat,matmul(identity(n+1),copymat[n][n]));
			}
			double[][][] qr = qr_decompose(implShiftMat);
			for(int j=0; j<mat.length; j++)
				for(int i=0; i<mat.length; i++)
					if(j<=n && i<=n) {
						matQ[j][i] = qr[1][i][j]; //implicit transpose
					} else {
						matQ[j][i] = (j==i ? 1d : 0d);
					}
			copymat = matmul(transpose(matQ),matmul(copymat,matQ));
			eigenvectors = matmul(eigenvectors,matQ);
			//System.out.println(" ["+FormatHelper.nf(iter+1,3)+"]  p="
			//		+FormatHelper.nf(n+1,2)+"  H[p-1,p]="
			//		+copymat[n][n-1]+"  H[p-2,p-1]="+copymat[n-1][n-2]);
			int numEV = 0;
			if(n>0)
				numEV = Math.abs(copymat[n][n-1])<eps*(Math.abs(copymat[n-1][n-1])+Math.abs(copymat[n][n])) ? 1 : 0;
			if(n>1 && numEV==0)
				numEV = Math.abs(copymat[n-1][n-2])<eps*(Math.abs(copymat[n-2][n-2])+Math.abs(copymat[n-1][n-1])) ? 2 : 0;
			if(n<2) { numEV = n+1; if(startIter<0) startIter = iter+1; }
			iter++;
			if(startIter>=0 && iter<startIter+0) continue;
			switch(numEV) {
				case 1:
					if(n>0) copymat[n][n-1] = 0d;
					eigenvaluesID[n] = n;
					eigenvalues[n][0] = copymat[n][n]; eigenvalues[n][1] = 0d;
					break;
				case 2:
					if(n>1) copymat[n-1][n-2] = 0d;
					double r = 0.5d*(copymat[n-1][n-1]+copymat[n][n]);
					double d = copymat[n-1][n-1]*copymat[n][n] - copymat[n-1][n]*copymat[n][n-1];
					double i = Math.sqrt(Math.abs(d-r*r));
					if(d<r*r) {
						double p = copymat[n][n-1];
						double y = 0.5d*(copymat[n][n]-copymat[n-1][n-1]);
						double dd = Math.abs(y)+Math.sqrt(p*p+y*y);
						double rr = Math.sqrt(p*p+dd*dd);
						double c = dd / rr,
							   s =  p / rr;
						if(y<0d) s=-s;
						copymat[n][n-1] = 0d;
						for(int j=0; j<mat.length; j++) rotateGivens(c,s,copymat,j,n-1,j,n);
						for(int j=0; j<mat.length; j++) rotateGivens(c,s,copymat,n-1,j,n,j);
						for(int j=0; j<mat.length; j++) rotateGivens(c,s,eigenvectors,n-1,j,n,j);
						eigenvalues[n-1][0] = r-i; eigenvalues[n-1][1] = 0d;
						eigenvalues[n  ][0] = r+i; eigenvalues[n  ][1] = 0d;
					} else {
						eigenvalues[n-1][0] = r; eigenvalues[n-1][1] = -i;
						eigenvalues[n  ][0] = r; eigenvalues[n  ][1] =  i;
					}
					eigenvaluesID[n-1] = n;
					eigenvaluesID[n  ] = n;
					break;
				default:
					break;
			}
			n -= numEV;
			startIter = -1;
		}
		eigenvectors = matmul(eigenvectors, getEigenvectorsFromUpperTriangularMatrix(copymat));
		//DEBUG
		//FormatHelper.printMat(System.out, copymat);
		//check right order of eigenvalues:
		for(n=mat.length-1; n>0; n--) {
			if(eigenvaluesID[n-1]==eigenvaluesID[n]) {
				if(eigenvalues[n-1][1]==0d && eigenvalues[n][1]==0d) {
					double[] ev = new double[mat.length];
					for(int e=0; e<mat.length; e++)
						ev[e] = eigenvectors[e][n-1];
					double[] testvec = matmul(mat, ev);
					double scaleXY = 0d, scaleXX = 0d;
					for(int e=0; e<mat.length; e++)
						{ scaleXX+=ev[e]*ev[e]; scaleXY+=ev[e]*testvec[e]; }
					double testvalue=scaleXY/scaleXX, should=eigenvalues[n-1][0];
					if(Math.abs(testvalue-should)>Math.max(1.0e-32d, 1.0e-12d*Math.abs(should))) {
						double[] temp = {eigenvalues[n-1][0], eigenvalues[n-1][1]};
						eigenvalues[n-1][0] = eigenvalues[n][0];
						eigenvalues[n-1][1] = eigenvalues[n][1];
						eigenvalues[n][0] = temp[0];
						eigenvalues[n][1] = temp[1];
					}
				}
			}
		}
		
		double[][] res = new double[mat.length+2][mat.length];
		for(int j=0; j<mat.length; j++) {
			res[mat.length  ][j] = eigenvalues[j][0];
			res[mat.length+1][j] = eigenvalues[j][1];
			for(int i=0; i<mat.length; i++)
				res[j][i] = eigenvectors[j][i];
		}
		return res;
	}
	/**
	 * Apply Jacobi-Eigenvalue-algorithm on matrix to extract eigenvalues and eigenvectors
	 * @param mat real symmetric matrix
	 * @return matrix (as double[][]) with eigenvectors columnwise and corresponding eigenvalues in the last row
	 */
	public static double[][] eigenJacobi(double[][] mat) {
		double[][] copymat = new double[mat.length][mat.length];
		for(int j=0; j<mat.length; j++)
			for(int i=0; i<mat.length; i++)
				copymat[j][i] = mat[j][i];
		boolean[] change = new boolean[mat.length];
		int[] indices = new int[mat.length];
		double[] eigenvalues = new double[mat.length];
		double[][] eigenvectors = identity(mat.length);
		//initialization
		for(int e=0; e<mat.length; e++) {
			indices[e] = index_largestValue(e,copymat);
			eigenvalues[e] = copymat[e][e];
			change[e] = true;
		}
		int status = mat.length-1;
		int iter = 0, sweeps = 0, sweepSize = mat.length*(mat.length-1)/2;
		while(status!=0 && sweeps<1000) {
			//find pivot
			int m=0;
			for(int n=1; n<mat.length-1; n++)
				if(Math.abs(copymat[indices[n]][n])>Math.abs(copymat[indices[m]][m]))
					m = n;
			int k = m,
				l = indices[m];
			double p = copymat[l][k];
			//calculate c=cos(phi) and s=sin(phi)
			double y = 0.5d*(eigenvalues[l]-eigenvalues[k]);
			double d = Math.abs(y)+Math.sqrt(p*p+y*y);
			double r = Math.sqrt(p*p+d*d);
			double c = d / r,
				   s = p / r,
				   t = p*p/d;
			if(y<0d) { s=-s; t=-t; }
			copymat[l][k] = 0d;
			eigenvalues[k] -= t; eigenvalues[l] += t;
			if(Math.abs(t)>1.0e-32d) {
				if(!change[k]) { change[k]=true; status++; }
				if(!change[l]) { change[l]=true; status++; }
			} else {
				if(change[k]) { change[k]=false; status--; }
				if(change[l]) { change[l]=false; status--; }
			}
			for(int i=0; i<k; i++) rotateGivens(c,s,copymat,i,k,i,l);
			for(int i=k+1; i<l; i++) rotateGivens(c,s,copymat,k,i,i,l);
			for(int i=l+1; i<mat.length; i++) rotateGivens(c,s,copymat,k,i,l,i);
			for(int i=0; i<mat.length; i++) rotateGivens(c,s,eigenvectors,k,i,l,i);
			indices[k] = index_largestValue(k,copymat);
			indices[l] = index_largestValue(l,copymat);
			iter++;
			if(iter==sweepSize) { sweeps++; iter = 0; }
			//System.out.println("[JACOBI EIGEN#-ALGO] status="+status+", iteration="+(iter==0?sweepSize:iter)+", sweep="+sweeps);
			//System.out.println("[JACOBI EIGEN#-ALGO] DEBUG: "+(180d*Math.acos(c)/Math.PI));
			//FormatHelper.printMat(System.out, copymat);
		}
		//DEBUG:
		System.out.println("[JACOBI-ALGORITHM] status="+status+", iteration="+(iter==0?sweepSize:iter)+", sweep="+sweeps);
		//FormatHelper.printMat(System.out, copymat);
//		for(int m=mat.length-1; m>0; m--) {
//			int n = 0;
//			for(int o=1; o<m; o++)
//				if(Math.abs(eigenvalues[o])<Math.abs(eigenvalues[n])) n = o;
//			if(Math.abs(eigenvalues[n])<Math.abs(eigenvalues[m])) {
//				double t = eigenvalues[n];
//				eigenvalues[n] = eigenvalues[m];
//				eigenvalues[m] = t;
//				for(int e=0; e<mat.length; e++) {
//					t = eigenvectors[e][n];
//					eigenvectors[e][n] = eigenvectors[e][m];
//					eigenvectors[e][m] = t;
//				}
//			}
//		}
		double[][] res = new double[mat.length+2][mat.length];
		for(int j=0; j<mat.length; j++) {
			res[mat.length  ][j] = eigenvalues[j];
			res[mat.length+1][j] = 0d;
			for(int i=0; i<mat.length; i++)
				res[j][i] = eigenvectors[j][i];
		}
		return res;
	}
	private static void rotateGivens(double _c, double _s, double[][] _mat, int _k, int _l, int _i, int _j) {
		double tcs = _c*_mat[_l][_k] - _s*_mat[_j][_i];
		double tsc = _s*_mat[_l][_k] + _c*_mat[_j][_i];
		_mat[_l][_k] = tcs;
		_mat[_j][_i] = tsc;
	}
	private static int index_largestValue(int _column, double[][] _mat) {
		int m = _column+1;
		for(int n=_column+2; n<_mat.length; n++)
			if(Math.abs(_mat[n][_column])>Math.abs(_mat[m][_column])) m = n;
		return m;
	}
	
	public static double[] eigenVectorForEigenvalue(double[][] mat, double[] eigenvalue) {
		return eigenvectorForEigenvalue(mat, eigenvalue, null);
	}
	public static double[] eigenvectorForEigenvalue(double[][] mat, double[] eigenvalue, double[] startVector) {
		double maxAbs = Math.sqrt(eigenvalue[0]*eigenvalue[0]+eigenvalue[1]*eigenvalue[1]);
		for(int j=0; j<mat.length; j++)
			for(int i=0; i<mat[j].length; i++)
				if(maxAbs<Math.abs(mat[j][i])) maxAbs = Math.abs(mat[j][i]);
		//double minAbs = maxAbs * 1.0e-16d;
		double[] firstEV, secondEV;
		firstEV = new double[mat.length];
		secondEV = new double[mat.length];
		for(int e=0; e<mat.length; e++) { firstEV[e] = 0d; secondEV[e] = 0d; }
		double diff = 0.00001d * maxAbs;
		boolean testSimilar = false;
		while(!testSimilar) {
			double ev = eigenvalue[0] - diff;
			double[][] reducedMat = new double[mat.length][mat[0].length];
			for(int j=0; j<mat.length; j++)
				for(int i=0; i<mat[j].length; i++)
					reducedMat[j][i] = mat[j][i] - (j==i ? ev : 0d);
			double[][] imat = inverse(reducedMat);
			int iter=0, maxIter=1000;
			boolean converged = false;
			double r = 0d, dot = 0d;
			for(int e=0; e<mat.length; e++) {
				if(startVector==null) {
					secondEV[e] = Math.random();
				} else {
					secondEV[e] = startVector[e];
				}
				r += secondEV[e]*secondEV[e];
			}
			r = Math.sqrt(r);
			for(int e=0; e<mat.length; e++)
				secondEV[e] /= r;
			while(!converged && iter<maxIter) {
				double[] prod = matmul(imat,secondEV);
				r = 0d; dot = 0d;
				for(int e=0; e<mat.length; e++)
					{ r += prod[e]*prod[e]; dot += prod[e]*secondEV[e]; }
				r = Math.sqrt(r);
				converged = Math.abs(Math.abs(dot/r)-1d)<1.0e-12d;
				iter++;
				for(int e=0; e<mat.length; e++)
					secondEV[e] = prod[e] / r;
			}
			dot = 0d;
			for(int e=0; e<mat.length; e++)
				dot += firstEV[e]*secondEV[e];
			double preSign = 1d;
			for(int e=0; e<mat.length; e++)
				if(Math.abs(secondEV[e])>1.0e-12d) { preSign = (secondEV[e]<0d ? -1d : 1d); break; }
			for(int e=0; e<mat.length; e++)
				firstEV[e] = preSign * secondEV[e];
			if(Math.abs(Math.abs(dot)-1d)<1.0e-12d)
				break;
		}
		return null;
	}
//	public static double[][] eigenvectorsQR(double[][] mat) {
//		double[][] a = new double[mat.length][mat[0].length];
//		double[][] qp = new double[mat.length][mat[0].length];
//		double[][] rp = new double[mat.length][mat[0].length];
//		for(int j=0; j<mat.length; j++) for(int i=0; i<mat[j].length; i++) {
//			a[j][i] = mat[j][i];
//			qp[j][i] = (j==i ? 1d : 0d);
//			rp[j][i] = mat[j][i];
//		}
//		boolean shouldProceed = true;
//		int iter = 0;
//		while(shouldProceed && iter<30) {
//			//shift
//			double[][] b = a;
//			if(a.length>1) {
//				int ap = a.length-2;
//				double[][] a22 = { {a[ap][ap],a[ap][ap+1]},{a[ap+1][ap],a[ap+1][ap+1]} };
//				double asum = a22[0][0]+a22[1][1], aprod = a22[0][0]*a22[1][1]-a22[0][1]*a22[1][0];
//				b = matmul(a,a);
//				b = matsub(b,matmul(a,asum));
//				b = matadd(b,matmul(identity(a.length),aprod));
//			}
//			//QR decomposition and reverse multiplication
//			double[][] qr_dec = qr_decompose(b);
//			double[][] r = matmul(qr_dec,b);
//			qp = matmul(qp,transpose(qr_dec));
//			a = matmul(matmul(transpose(qr_dec),a),qr_dec);
//			String[][] o = new String[mat.length][3*mat[0].length];
//			for(int row=0; row<mat.length; row++) {
//				int of = 0;
//				for(int col=0; col<mat[row].length; col++)
//					o[row][col] = " "+FormatHelper.nf(qr_dec[row][col], 7,4)+" ";
//				of += mat[row+of].length;
//				for(int col=0; col<mat[row].length; col++)
//					o[row][col+of] = " "+FormatHelper.nf(r[row][col], 7,4)+" ";
//				of += mat[row].length;
//				for(int col=0; col<mat[row].length; col++)
//					o[row][col+of] = " "+FormatHelper.nf(a[row][col], 7,4)+" ";
//			}
//			int maxlen = 0;
//			for(int j=0; j<mat.length; j++) for(int i=0; i<o[j].length; i++)
//				if(maxlen<o[j][i].length()) maxlen = o[j][i].length();
//			for(int j=0; j<mat.length; j++) for(int i=0; i<o[j].length; i++)
//				while(maxlen>o[j][i].length()) o[j][i] = " "+o[j][i];
//			int ml = mat.length;
//			for(int j=0; j<mat.length; j++) {
//				String _in = (j==0 ? (j==ml-1 ? "(" : "/") : (j==ml-1 ? "\\" : "|"));
//				o[j][0]    = "    "+_in+o[j][0];
//				o[j][ml]   = "    "+_in+o[j][ml];
//				o[j][2*ml] = "    "+_in+o[j][2*ml];
//				String _out = (j==0 ? (j==ml-1 ? ")" : "\\") : (j==ml-1 ? "/" : "|"));
//				o[j][ml-1]   = o[j][ml-1]+_out;
//				o[j][2*ml-1] = o[j][2*ml-1]+_out;
//				o[j][3*ml-1] = o[j][3*ml-1]+_out;
//			}
//			FormatHelper.printTable(o);
//			iter++;
//			shouldProceed = false;
//			for(int j=0; j<mat.length; j++) for(int i=0; i<mat[j].length; i++) {
//				if(j==i) continue;
//				if(Math.abs(qr_dec[j][i])>1.0e-6d) shouldProceed = true;
//			}
//		}
//		double[][] res = new double[mat.length+1][mat[0].length];
//		for(int j=0; j<mat.length; j++)
//			for(int i=0; i<mat[j].length; i++)
//				res[j][i] = qp[j][i];
//		for(int i=0; i<mat[0].length; i++)
//			res[mat.length][i] = a[i][i];
//		return res;
//	}
	/**
	 * 
	 * @param mat
	 * @return
	 */
	public static double[][][] transformToSimilarHessenberg(double[][] mat) {
		return transformToSimilarUpperHessenberg(mat);
	}
	public static double[][][] transformToSimilarUpperHessenberg(double[][] mat) {
		double[][] hsB = new double[mat.length][mat[0].length];
		double[][] q_h = new double[mat.length][mat[0].length];
		for(int j=0; j<mat.length; j++) for(int i=0; i<mat[j].length; i++) {
			hsB[j][i] = mat[j][i];
			q_h[j][i] = (j==i ? 1d : 0d);
		}
		for(int j=0; j<mat.length-2; j++) {
			double[][] hm = getHouseholderMatrix_column(hsB, j, 1);
			hsB = matmul(matmul(hm,hsB),hm);
			q_h = matmul(hm,q_h);
		}
		for(int j=2; j<mat.length; j++)
			for(int i=0; i<j-1; i++)
				hsB[j][i] = 0d;
		if(isSymmetric(mat))
			for(int j=0; j<mat.length-2; j++)
				for(int i=j+2; i<mat.length; i++)
					hsB[j][i] = 0d;
		return new double[][][] {hsB,q_h};
	}
	public static double[][][] transformToSimilarLowerHessenberg(double[][] mat) {
		double[][] hsB = new double[mat.length][mat[0].length];
		double[][] q_h = new double[mat.length][mat[0].length];
		for(int j=0; j<mat.length; j++) for(int i=0; i<mat[j].length; i++) {
			hsB[j][i] = mat[j][i];
			q_h[j][i] = (j==i ? 1d : 0d);
		}
		for(int j=0; j<mat.length-2; j++) {
			double[][] hm = getHouseholderMatrix_row(hsB, j, 1);
			hsB = matmul(matmul(hm,hsB),hm);
			q_h = matmul(hm,q_h);
		}
		for(int j=0; j<mat.length-2; j++)
			for(int i=j+2; i<mat.length; i++)
				hsB[j][i] = 0d;
		if(isSymmetric(mat))
			for(int j=2; j<mat.length; j++)
				for(int i=0; i<j-1; i++)
					hsB[j][i] = 0d;
		return new double[][][] {hsB,q_h};
	}
	/**
	 * 
	 * @param mat
	 * @return
	 */
	public static double[][][] qr_decompose(double[][] mat) {
		double[][] tri = new double[mat.length][mat[0].length];
		double[][] q_m = new double[mat.length][mat[0].length];
		for(int j=0; j<mat.length; j++)
			for(int i=0; i<mat[j].length; i++) {
				tri[j][i] = mat[j][i];
				q_m[j][i] = (j==i ? 1d : 0d);
			}
		for(int j=0; j<mat.length-1; j++) {
			double[][] hm = getHouseholderMatrix_column(tri, j, 0);
			tri = matmul(hm,matmul(tri,hm));
			q_m = matmul(hm,q_m);
		}
		return new double[][][] {tri,q_m};
	}
	/**
	 * calculates a matrix to perform a Householder-transformation
	 * @param mat matrix to calculate the Householder transformation for
	 * @param _column column to reduce
	 * @param _nOffdia number of offdiagonals
	 * @return
	 */
	public static double[][] getHouseholderMatrix(double[][] mat, int _column, int _nOffdia) {
		return getHouseholderMatrix_column(mat, _column, _nOffdia);
	}
	public static double[][] getHouseholderMatrix_column(double[][] mat, int _column, int _nOffdia) {
		int dim = mat.length;
		if(mat[0].length!=dim) {
			System.err.println("Cannot calculate Householder matrix of a non-square matrix!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		if(_column+_nOffdia>=dim) {
			System.err.println("Cannot calculate "+(_column+1)+"th householder matrix with "+_nOffdia+" offdiagonals for a "+dim+"x"+dim+" matrix!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		double[][] hm = new double[dim][dim];
		double[] v = new double[dim];
		double asum = 0d;
		int piv = _column+_nOffdia;
		//System.out.println("column="+_column+", #off="+_nOffdia+", pivot=("+_column+","+piv+")");
		for(int i=piv; i<dim; i++)
			asum += mat[i][_column]*mat[i][_column];
		double alpha = (mat[piv][_column]<0d ? 1d : -1d)*Math.sqrt(asum);
		double r = Math.sqrt(2*alpha*(alpha-mat[piv][_column]))+1.0e-200d;
		//System.out.println("lambda = "+alpha+"\nmagnitude = "+r);
		for(int i=0; i<piv; i++)
			v[i] = 0d;
		v[piv] = (mat[piv][_column]-alpha)/r;
		for(int i=piv+1; i<dim; i++)
			v[i] = mat[i][_column]/r;
		//FormatHelper.printVec(System.out, v, false);
		for(int j=0; j<dim; j++)
			for(int i=0; i<dim; i++)
				hm[j][i] = (i==j ? 1d : 0d) - 2*v[i]*v[j];
		return hm;
	}
	public static double[][] getHouseholderMatrix_row(double[][] mat, int _column, int _nOffdia) {
		int dim = mat.length;
		if(mat[0].length!=dim) {
			System.err.println("Cannot calculate Householder matrix of a non-square matrix!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		if(_column+_nOffdia>=dim) {
			System.err.println("Cannot calculate "+(_column+1)+"th householder matrix with "+_nOffdia+" offdiagonals for a "+dim+"x"+dim+" matrix!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		double[][] hm = new double[dim][dim];
		double[] v = new double[dim];
		double asum = 0d;
		int piv = _column+_nOffdia;
		//System.out.println("column="+_column+", #off="+_nOffdia+", pivot=("+_column+","+piv+")");
		for(int i=piv; i<dim; i++)
			asum += mat[_column][i]*mat[_column][i];
		double alpha = (mat[_column][piv]<0d ? 1d : -1d)*Math.sqrt(asum);
		double r = Math.sqrt(2*alpha*(alpha-mat[_column][piv]))+1.0e-200d;
		//System.out.println("lambda = "+alpha+"\nmagnitude = "+r);
		for(int i=0; i<piv; i++)
			v[i] = 0d;
		v[piv] = (mat[_column][piv]-alpha)/r;
		for(int i=piv+1; i<dim; i++)
			v[i] = mat[_column][i]/r;
		//FormatHelper.printVec(System.out, v, false);
		for(int j=0; j<dim; j++)
			for(int i=0; i<dim; i++)
				hm[j][i] = (i==j ? 1d : 0d) - 2*v[i]*v[j];
		return hm;
	}
	/**
	 * calculates the eigenvectors of an upper triangular matrix
	 * @param mat upper triangular matrix
	 * @return matrix (double[][]) with eigenvectors on the columns
	 */
	public static double[][] getEigenvectorsFromUpperTriangularMatrix(double[][] mat) {
		double[][] evs = new double[mat.length][mat[0].length];
		for(int k=0; k<mat.length; k++) {
			evs[k][k] = 1d;
			for(int i=k-1; i>=0; i--) {
				double dotprod = 0d;
				for(int l=i+1; l<=k; l++)
					dotprod += mat[i][l]*evs[l][k];
				evs[i][k] = dotprod / (mat[k][k]-mat[i][i]);
			}
			double norm = 0d;
			for(int i=0; i<=k; i++)
				norm += evs[i][k]*evs[i][k];
			norm = 1d/Math.sqrt(norm);
			for(int i=0; i<=k; i++)
				evs[i][k] *= norm;
		}
		System.out.println("DEBUG: EVs of upper tri.:");
		FormatHelper.printMat(System.out, evs);
		return evs;
	}
	
	// ************************************************************ //
	// ** SETUP ROTATION MATRICIES                               ** //
	// ************************************************************ //
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
	 *              anis    Ratio of minor axis length to major axis length
	 *     
	 *     Converts the input angle to an angle which make more mathematical sense:
	 *              alpha   angle between the major axis of anisotropy and the
	 *                      E-W axis. Note: Counter clockwise is positive.
	 *     
	 *     After that following scheme is used:
	 *     
	 *                 -> * Scale(1,anis)   -> * Rot(alpha)  ->
	 *         Circle                                           Ellipse
	 *                 <- * Scale(1,1/anis) <- * Rot(-alpha) <-
	 *     
	 *     The final transformation matrix got the following look:
	 *     
	 *         M = [ 1    0   ] * [ cos(-alpha) -sin(-alpha) ]
	 *             [ 0 1/anis ]   [ sin(-alpha)  cos(-alpha) ]
	 *     
	 *         with abbreviations: cos($)=c$, sin($)=s$, alpha=a
	 *     
	 *         M = [   ca         sa       ]
	 *             [ (-sa)/anis  (ca)/anis ]
	 * </p>
	 * <p>
	 *     if axis-length is negative assume hyperbolic curve
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
	public static double[][] setrot2D(double ang, double anis, int ind, int maxrot, double[][] rotmat) {
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
//        c   ang              Azimuth angle for principal direction
//        c   anis             First anisotropy ratio
//        c   ind              matrix indicator to initialise
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
		double afac,sina,cosa;
		double alpha;
//        c
//        c Converts the input angles to three angles which make more
//        c  mathematical sense:
//        c
//        c         alpha   angle between the major axis of anisotropy and the
//        c                 E-W axis. Note: Counter clockwise is positive.
//        c
		if(ang>=0d && ang<270d) {
			alpha = (90d  - ang) * Constants.DEG2RAD;
		} else {
			alpha = (450d - ang) * Constants.DEG2RAD;
		}
//        c
//        c Get the required sines and cosines:
//        c
		sina  = Math.sin(alpha);
		cosa  = Math.cos(alpha);
//        c
//        c Construct the rotation matrix in the required memory:
//        c
		afac = 1d / Math.max(Math.abs(anis),Constants.D_EPSLON);
		if(anis<0d) afac *= -1d;
		rotmat[ind-1][0] =      ( cosa);
		rotmat[ind-1][1] =      ( sina);
		rotmat[ind-1][2] = afac*(-sina);
		rotmat[ind-1][3] = afac*( cosa);
//        c
//        c Return to calling program:
//        c
		return rotmat;
	}
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
	 *                      direction of the major axis (view from inside to outside)
	 *     
	 *     After that following scheme is used:
	 *     
	 *                 -> * Scale(1,anis1,anis2)     -> * Rot(x-axis,theta)  -> * Rot(y-axis,beta)  _-> * Rot(z-axis,alpha) ->
	 *         Sphere                                                                                                          Ellipsoid
	 *                 <- * Scale(1,1/anis1,1/anis2) <- * Rot(x-axis,-theta) <- * Rot(y-axis,-beta) <- * Rot(z-axis,-alpha) <-
	 *     
	 *     The final transformation matrix got the following look:
	 *     
	 *             [ 1    0       0    ]   [ 1     0            0       ]   [  cos(-beta) 0 sin(-beta) ]   [ cos(-alpha) -sin(-alpha) 0 ]
	 *         M = [ 0 1/anis1    0    ] * [ 0 cos(-theta) -sin(-theta) ] * [      0      1     0      ] * [ sin(-alpha)  cos(-alpha) 0 ]
	 *             [ 0    0    1/anis2 ]   [ 0 sin(-theta)  cos(-theta) ]   [ -sin(-beta) 0 cos(-beta) ]   [     0            0       1 ]
	 *     
	 *         with abbreviations: cos($)=c$, sin($)=s$, alpha=a, beta=b, theta=t
	 *     
	 *             [           ca*cb                   sa*cb            -sb        ]
	 *         M = [ (ca*sb*st-sa*ct)/anis1  (sa*sb*st+ca*ct)/anis1  (cb*st)/anis1 ]
	 *             [ (ca*sb*ct+sa*st)/anis2  (sa*sb*ct-ca*st)/anis2  (cb*ct)/anis2 ]
	 * </p>
	 * <p>
	 *     if axis-length is negative assume hyperboloid
	 *     one negative axis leads to a one-sheet hyperboloid
	 *     two negative axes leads to a two-sheet hyperboloid
	 * </p>
	 * @param ang1      Azimuth angle for principal direction in geographical notation (N=0°, E=90°, S=180°, W=270°)
	 * @param ang2      Dip angle for principal direction
	 * @param ang3      Third rotation angle
	 * @param anis1     First anisotropy ratio
	 * @param anis2     Second anisotropy ratio
	 * @param ind       matrix indicator to initialize
	 * @param maxrot    maximum number of rotation matrices dimensioned
	 * @param rotmat    rotation matrices
	 * @return
	 */
	public static double[][][] setrot3Dgeo(double ang1, double ang2, double ang3, double anis1, double anis2, int ind, int maxrot, double[][][] rotmat) {
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
//        c   ind              matrix indicator to initialise
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
	 *                      direction of the major axis (view from inside to outside)
	 *     
	 *     After that following scheme is used:
	 *     
	 *                 -> * Scale(1,anis1,anis2)     -> * Rot(x-axis,theta)  -> * Rot(y-axis,beta)  _-> * Rot(z-axis,alpha) ->
	 *         Sphere                                                                                                          Ellipsoid
	 *                 <- * Scale(1,1/anis1,1/anis2) <- * Rot(x-axis,-theta) <- * Rot(y-axis,-beta) <- * Rot(z-axis,-alpha) <-
	 *     
	 *     The final transformation matrix got the following look:
	 *     
	 *             [ 1    0       0    ]   [ 1     0            0       ]   [  cos(-beta) 0 sin(-beta) ]   [ cos(-alpha) -sin(-alpha) 0 ]
	 *         M = [ 0 1/anis1    0    ] * [ 0 cos(-theta) -sin(-theta) ] * [      0      1     0      ] * [ sin(-alpha)  cos(-alpha) 0 ]
	 *             [ 0    0    1/anis2 ]   [ 0 sin(-theta)  cos(-theta) ]   [ -sin(-beta) 0 cos(-beta) ]   [     0            0       1 ]
	 *     
	 *         with abbreviations: cos($)=c$, sin($)=s$, alpha=a, beta=b, theta=t
	 *     
	 *             [           ca*cb                   sa*cb            -sb        ]
	 *         M = [ (ca*sb*st-sa*ct)/anis1  (sa*sb*st+ca*ct)/anis1  (cb*st)/anis1 ]
	 *             [ (ca*sb*ct+sa*st)/anis2  (sa*sb*ct-ca*st)/anis2  (cb*ct)/anis2 ]
	 * </p>
	 * <p>
	 *     if axis-length is negative assume hyperboloid
	 *     one negative axis leads to a one-sheet hyperboloid
	 *     two negative axes leads to a two-sheet hyperboloid
	 * </p>
	 * @param ang1      Azimuth angle for principal direction in mathematical notation (N=90°, E=0°, S=270°, W=180°)
	 * @param ang2      Dip angle for principal direction
	 * @param ang3      Third rotation angle
	 * @param anis1     First anisotropy ratio
	 * @param anis2     Second anisotropy ratio
	 * @param ind       matrix indicator to initialize
	 * @param maxrot    maximum number of rotation matrices dimensioned
	 * @param rotmat    rotation matrices
	 * @return
	 */
	public static double[][][] setrot3Dmath(double ang1, double ang2, double ang3, double anis1, double anis2, int ind, int maxrot, double[][][] rotmat) {
		//add own extra
		if(ind>maxrot) throw new ArrayIndexOutOfBoundsException(ind+" of maxrot="+maxrot);
		//end own addition
		double afac1,afac2,sina,sinb,sint,cosa,cosb,cost;
		double alpha,beta,theta;
		alpha = ang1 * Constants.DEG2RAD;
		beta  = -1d * ang2 * Constants.DEG2RAD;
		theta = ang3 * Constants.DEG2RAD;
//      c
//      c Get the required sines and cosines:
//      c
		sina  = Math.sin(alpha);
		sinb  = Math.sin(beta);
		sint  = Math.sin(theta);
		cosa  = Math.cos(alpha);
		cosb  = Math.cos(beta);
		cost  = Math.cos(theta);
//      c
//      c Construct the rotation matrix in the required memory:
//      c
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
//      c
//      c Return to calling program:
//      c
		return rotmat;
	}

	public static double[] parameterFromAnistropicRotationMatrix(double[][] rotmat) {
		//double[][] eigen = eigenJacobi(matmul(transpose(rotmat),rotmat));
		double[] infoParams;
		switch(rotmat.length) {
			case 1: infoParams = new double[] {1d}; break;
			case 2: infoParams = new double[] {1d,1d,0d}; break;
			case 3: infoParams = new double[] {1d,1d,1d,0d,0d,0d}; break;
			default: infoParams = new double[0]; break;
		}
		//double[][] decomposition = eigenQR(rotmat);
		double[][] evs = new double[rotmat.length][2];
		boolean checkImaginary=false, checkNegative=false;
		for(int e=0; e<rotmat.length; e++) {
			if(evs[e][1]!=0d) checkImaginary = true;
			if(evs[e][0]<0d) checkNegative = true;
		}
		if(checkImaginary) {
			System.err.println("[ERROR] At least one axis is not real, probably a hyperboloid would be a better fit.");
			DataHelper.printStackTrace(System.err);
			return infoParams;
		}
		if(checkNegative) {
			System.out.println("[WARNING] At least one axis is negative, probably a hyperboloid would be a better fit.");
			DataHelper.printStackTrace(System.out);
			return infoParams;
		}

		double[] values = new double[rotmat.length];
		double[][] res_qr = eigenQR(rotmat);
		double[][] vectors = new double[rotmat.length][rotmat.length];
		for(int j=0; j<rotmat.length; j++) {
			values[j] = res_qr[rotmat.length][j];
			for(int i=0; i<rotmat.length; i++)
				vectors[j][i] = res_qr[j][i];
		}
		//order by absolute value of eigenvalue
		for(int e=0; e<rotmat.length-1; e++) {
			int largest = e;
			double le = Math.abs(values[e]);
			for(int f=e+1; f<rotmat.length; f++) {
				double fe = Math.abs(values[f]);
				if(fe>le) { largest = f; le = fe; }
			}
			if(e!=largest) {
				double t = Math.abs(values[e]);
				values[e] = Math.abs(values[largest]);
				values[largest] = t;
				for(int f=0; f<rotmat.length; f++) {
					t = vectors[e][f];
					vectors[e][f] = vectors[largest][f];
					vectors[largest][f] = t;
				}
			}
		}
		int paramlen = 0;
		switch(rotmat.length) {
			case 1: paramlen = 1; break;
			case 2: paramlen = 3; break;
			case 3: paramlen = 6; break;
			default: break;
		}
		double[] params = new double[paramlen];
		//extract ellipsoids semi axes
		for(int j=0; j<rotmat.length; j++) params[j] = Math.abs(values[j]);
		//TODO extract ellipsoids orientation:
		switch(rotmat.length) {
			case 1: //no further direction:
				break;
			case 2:
				params[2] = Math.atan2(vectors[0][1], vectors[0][0]);
				break;
			case 3:
				double alpha = Math.atan2(vectors[0][0], vectors[0][1]);
				double delta = Math.atan2(vectors[0][2], Math.sqrt(vectors[0][0]*vectors[0][0]+vectors[0][1]*vectors[0][1]));
				params[3] = 180d*alpha/Math.PI;
				params[4] = 180d*delta/Math.PI;
				double[][] copyvec = new double[1][3];
				for(int j=0; j<rotmat.length; j++) copyvec[0][j] = vectors[1][j];
				rotateGivens(Math.cos(alpha), -Math.sin(alpha),copyvec,0,0,1,0);
				rotateGivens(Math.cos(delta), -Math.sin(delta),copyvec,0,0,2,0);
				double roll = Math.atan2(copyvec[0][2], copyvec[0][1]);
				params[5] = 180d*roll/Math.PI;
				break;
			default:
				break;
		}
		return params;
	}
	
	// ************************************************************ //
	// ** CALCULATE WITH ROTATION MATRICIES                      ** //
	// ************************************************************ //
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
	public static double sqdist2D(double x1, double y1, double x2, double y2, double[] rotmat) {
		double dx = (x1 - x2);
		double dy = (y1 - y2);
		double cont   = rotmat[0] * dx +rotmat[1] * dy;
		double sqdist = cont * cont;
		       cont   = rotmat[2] * dx +rotmat[3] * dy;
		       sqdist += cont * cont;
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
	public static double sqdist2D(double[] pos1, double[] pos2, double[] rotmat) {
		double d0 = pos1[0] - pos2[0];
		double d1 = pos1[1] - pos2[1];
		double cont   = rotmat[0] * d0 + rotmat[1] * d1;
		double sqdist = cont * cont;
		       cont   = rotmat[2] * d0 + rotmat[3] * d1;
		       sqdist += cont * cont;
		return sqdist;
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
	public static double sqdist3D(double x1, double y1, double z1, double x2, double y2, double z2,
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
	public static double sqdistKD(double[] pos1, double[] pos2, double[][] rotmat) {
		double sqdist = 0d;
		for(int j=0; j<rotmat.length; j++) {
			double cont = 0d;
			for(int i=0; i<rotmat[j].length; i++)
				cont += rotmat[j][i] * (pos1[i] - pos2[i]);
			sqdist += cont * cont;
		}
		return sqdist;
	}

	// ************************************************************ //
	// ** ELLIPSOID FITTING                                      ** //
	// ************************************************************ //
	/**
	 * returns a matrix, which describes the ellipsoid
	 * @param pointlist
	 * @param withWeights
	 * @return
	 */
	public static double[][] fitEllipsoidToPointCloud(List<double[]> pointlist, boolean withWeights) {
		int mindim = 10;
		int maxdim = 0;
		for(double[] p: pointlist) {
			int dim = p.length;
			if(dim<mindim) mindim = dim;
			if(dim>maxdim) maxdim = dim;
		}
		int dim = mindim;
		if(mindim!=maxdim) {
			System.err.println("List contains points from different spaces (dimensions are not the same for each point!).");
			DataHelper.printStackTrace(System.err);
			dim -= (withWeights ? 1 : 0);
			double[][] mat = new double[dim][dim];
			for(int j=0; j<dim; j++)
				for(int i=0; i<dim; i++)
					mat[j][i] = (i==j ? 1d : 0d);
			return mat;
		}
		double[][] pl = new double[pointlist.size()][mindim];
		for(int p=0; p<pl.length; p++) {
			double[] point = pointlist.get(p);
			for(int d=0; d<dim; d++) pl[p][d] = point[d];
		}
		return fitEllipsoidToPointCloud(pl, withWeights);
	}
	public static double[][] fitEllipsoidToPointCloud(double[][] pointlist, boolean withWeights) {
		int dim = pointlist[0].length - (withWeights?1:0);
		int count = pointlist.length;
		double[] axes = new double[dim];
		double[] sqaxes = new double[dim];
		for(int d=0; d<dim; d++) axes[d] = 1d;
		double[][] zeroRotmat = new double[dim][dim];
		for(int j=0; j<dim; j++) for(int i=0; i<dim; i++) zeroRotmat[j][i] = 0d;
		double[][] rotmat = new double[dim][dim];
		double[][] zealp = new double[count][dim];
		double[][] refp = new double[count][dim];
		for(int p=0; p<count; p++) {
			double wgt = 1d;
			if(withWeights) wgt = pointlist[p][dim];
			for(int d=0; d<dim; d++)
				zealp[p][d] = pointlist[p][d] * wgt;
		}
		for(int k=0; k<10; k++) {
			for(int d=0; d<dim; d++) sqaxes[d] = axes[d]*axes[d];
			for(int p=0; p<count; p++) {
				double rad =  0d;
				for(int d=0; d<dim; d++)
					rad += pointlist[p][d]*pointlist[p][d] * sqaxes[d];
				rad = Math.sqrt(rad) + 1.0e-200d;
				double wgt = 1d;
				if(withWeights) wgt = pointlist[p][dim];
				for(int d=0; d<dim; d++)
					refp[p][d] = wgt * pointlist[p][d] * axes[d] / rad;
			}
			double[][] Xt = transpose(refp);
			double[][] XtX = matmul(Xt, refp);
			double[][] XtXinvXt = matmul(inverse(XtX),Xt);
			rotmat = matmul(XtXinvXt,zealp);
			//test, check via distance sum
			double diviation = 0d;
			for(int p=0; p<count; p++) {
				double[] test = matmul(rotmat, zealp[p]);
				double r = 0d;
				for(int d=0; d<rotmat.length; d++)
					r += (test[d]-zealp[p][d])*(test[d]-zealp[p][d]);
				diviation += r;
			}
			diviation /= count;
			System.out.println("[FIT ELLIPSOID] DEBUG: diviation in iter "+(k+1)+" is: "+diviation);
//			double[][] eigenOut = eigenJacobi(matmul(transpose(rotmat),rotmat));
//			System.out.println("[FIT ELLIPSOID] DEBUG: iteration "+(k+1)+" results in vectors:");
//			FormatHelper.printMat(System.out, rotmat);
//			//FormatHelper.printMat(System.out, eigenOut);
			//axes = eigenOut[rotmat.length];
			double[][] evs = eigenvalues(rotmat);
			boolean checkImaginary = false;
			boolean negativeAxis = false;
			for(int i=0; i<dim; i++) {
				if(evs[i][1]!=0d) {
					{ checkImaginary = true; break; }
				} else {
					if(evs[i][0]<0d) negativeAxis = true;
					axes[i] = 1d / Math.abs(evs[i][0]);
				}
			}
			if(checkImaginary) {
				System.err.println("At least one axis is not real, probably a hyperboloid would be a better fit.");
				DataHelper.printStackTrace(System.err);
				return zeroRotmat;
			}
			if(negativeAxis)
				System.out.println("[WARNING] At least one axis is negative, probably a hyperboloid would be a better fit.");
		}
		return rotmat; //inverse(rotmat);
	}
	
	
	// ************************************************************ //
	// ** MAYBE NOT INCLUDET...                                  ** //
	// ************************************************************ //
	/**
	 * Calculates the Minimum volume enclosing ellipsoid around listed ellipsoids in 3 dimensions
	 * 
	 * @param major            array of all ellipsoids major semi axis
	 * @param minor            array of all ellipsoids minor semi axis
	 * @param third_axis       array of all ellipsoids third semi axis
	 * @param azimuth          array of all ellipsoids major semi axis azimuth position angle
	 * @param declination      array of all ellipsoids major semi axis azimuth position angle
	 * @param roll             array of all ellipsoids major semi axis azimuth position angle
	 * @param angles_in_degree are angles given in degrees?
	 * @return
	 */
	public static double[] getMinimumVolumeEnclosingEllipsoid(double[] major, double[] minor, double[] third_axis,
			double[] azimuth, double[] declination, double[] roll, boolean angles_in_degree) {
		
		
		return null;
	}
	/**
	 * Calculates the Minimum volume enclosing ellipsoid around listed ellipsoids in arbitrary dimensions
	 * 
	 * @param ellipsoids matrices to describe each ellipsoid to enclose by the MVEE.
	 * @return
	 */
	public static double[] getMinimumVolumeEnclosingEllipsoid(double[][] ellipsoids) {
		return null;
	}
}
