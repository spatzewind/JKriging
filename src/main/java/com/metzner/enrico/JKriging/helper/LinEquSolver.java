package com.metzner.enrico.JKriging.helper;

public class LinEquSolver {

	public static double[] ksol(int nright, int neq, int nsb, double[] a, double[] r) {
		int ii, ij=0, in, ll, ll1, nm1, km1;
		double piv;
		double[] s = new double[r.length]; // <-- solution array
		//FormatHelper.printTable(neq*(neq+1)/2, a, r);
//        c-----------------------------------------------------------------------
//        c
//        c                Solution of a System of Linear Equations
//        c                ****************************************
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   nright,nsb       number of columns in right hand side matrix.
//        c                      for KB2D: nright=1, nsb=1
//        c   neq              number of equations
//        c   a()              upper triangular left hand side matrix (stored 
//        c                      columnwise)
//        c   r()              right hand side matrix (stored columnwise)
//        c                      for kb2d, one column per variable
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   s()              solution array, same dimension as  r  above.
//        c   ising            singularity indicator
//        c                      0,  no singularity problem
//        c                     -1,  neq .le. 1
//        c                      k,  a null pivot appeared at the kth iteration
//        c
//        c
//        c
//        c PROGRAM NOTES:
//        c
//        c   1. Requires the upper triangular left hand side matrix.
//        c   2. Pivots are on the diagonal.
//        c   3. Does not search for max. element for pivot.
//        c   4. Several right hand side matrices possible.
//        c   5. USE for Ordinary Kriging and Simple Kriging only, NOT for UK.
//        c
//        c
//        c-----------------------------------------------------------------------
		//    implicit real*8 (a-h,o-z)
//        c
//        c If there is only one equation then set ising and return:
//        c
		if(neq<=1) {
			return null;
		}
//        c
//        c Initialize:
//        c
		double  tol   = 0.1e-06d;
		//boolean ising = false;
		int     nn    = neq*(neq+1)/2;
		int     nm    = nsb*neq;
		int     m1    = neq-1;
		int     kk    = -1;
//        c
//        c Start triangulation:
//        c
		for(int k=0; k<m1; k++) { // for each pivot
			kk+=k+1;
			double ak=a[kk];
			if(Math.abs(ak)<tol) {
				System.err.println("Found null pivot at row k="+k);
				return null;
			}
			km1=k;
			for(int iv=0; iv<nright; iv++) { // number of solved variables (number of lambda-columns)
				nm1 = nm*iv;
				ii  = kk+nn*iv;
				piv = 1d / a[ii];
				int lp=0;
				for(int i=k; i<m1; i++) { 
					ll  = ii;
					ii += i+1;
					double ap = a[ii]*piv;
					lp++;
					ij=ii-km1;
					for(int j=i+1; j<=m1; j++) {
						ij += j;
						ll += j;
						a[ij]-=ap*a[ll];
					}
					for(int llb=k; llb<nm; llb+=neq) {
						in  = llb+lp+nm1;
						ll1 = llb+nm1;
						r[in] -= ap*r[ll1];
					}
				}
			}
		}
//        c
//        c Error checking - singular matrix:
//        c
		int ijm = ij - nn*(nright-1);
		//ijm--; // Index shift FORTRAN -> JAVA
		if(Math.abs(a[ijm]) < tol) {
			System.err.println("[KSOL] found a[ijm] be less than tol="+tol);
			//ising=neq
			return null;
		}
//        c
//        c Finished triangulation, start solving back:
//        c
		for(int iv=0; iv<nright; iv++) {
			nm1 = nm*iv;
			ij  = ijm+nn*iv;
			piv = 1d / a[ij];
			for(int llb=neq-1; llb<nm; llb+=neq) {
				ll1    = llb+nm1;
				s[ll1] = r[ll1] * piv;
			}
			int i=neq-1;
			kk=ij;
			for(ii=0; ii<m1; ii++) {
				kk  -= i+1;
				piv  = 1d / a[kk];
				i--;
				for(int llb=i; llb<nm; llb+=neq) {
					ll1 = llb+nm1;
					in  = ll1;
					double ap = r[in];
					ij=kk;
					for(int j=i; j<m1; j++) {
						ij += j+1;
						in++;
						ap -= a[ij]*s[in];
					}
					s[ll1] = ap*piv;
				}
			}
		}
//        c
//        c Finished solving back, return:
//        c
		return s;
	}
	
	public static double[] ktsol(int n, int ns, int nv, double[] a, double[] b, int maxeq) {
//        c-----------------------------------------------------------------------
//        c
//        c Solution of a system of linear equations by gaussian elimination with
//        c partial pivoting.  Several right hand side matrices and several
//        c variables are allowed.
//        c
//        c
//        c         NOTE: All input matrices must be in double precision
//        c
//        c
//        c INPUT/OUTPUT VARIABLES:
//        c
//        c   n                Number of equations
//        c   ns               Number of right hand side matrices
//        c   nv               Number of variables.
//        c   a(n*n*nv)        left hand side matrices versus columnwise.
//        c   b(n*ns*nv)       input right hand side matrices.
//        c   x(n*ns*nv)       solution matrices.
//        c   ktilt            indicator of singularity
//        c                      =  0  everything is ok.
//        c                      = -1 n.le.1
//        c                      =  k  a null pivot appeared at the kth iteration.
//        c   tol              used in test for null pivot. depends on machine
//        c                      precision and can also be set for the tolerance
//        c                      of an ill-defined kriging system.
//        c
//        c
//        c-----------------------------------------------------------------------
		//    implicit real*8 (a-h,o-z)
		//    real*8 x(maxeq),a(maxeq*maxeq),b(maxeq)
//        c
//        c Make sure there are equations to solve:
//        c
		if(n<=1) {
			System.err.println("Too few equations, there have to be at least n=2!");
			//ktilt = -1;
			return null;
		}
//        c
//        c Initialization:
//        c
		double tol   = 0.1e-10d;
		//int ktilt = 0;
		int ntn   = n*n;
		int nm1   = n-1;
		int kp1,kdiag=0,npiv,ipiv,i1,i2,j1,j2, nva,nvb,nvb2, nmk, kb=0;
		double t;
//        c
//        c Triangulation is done variable by variable:
//        c
		for(int iv=1; iv<=nv; iv++) {
//        c
//        c Indices of location in vectors a and b:
//        c
			nva = ntn*(iv-1);
			nvb = n*ns*(iv-1);
//        c
//        c Gaussian elimination with partial pivoting:
//        c
			for(int k=1; k<=nm1; k++) {
				kp1 = k+1;
//        c
//        c Indice of the diagonal element in the kth row:
//        c
				kdiag = nva+(k-1)*n+k;
//        c
//        c Find the pivot - interchange diagonal element/pivot:
//        c
				npiv = kdiag;
				ipiv = k;
				i1   = kdiag;
				for(int i=kp1; i<=n; i++) {
					i1++;
					if(Math.abs(a[i1-1])>Math.abs(a[npiv-1])) {
						npiv = i1;
						ipiv = i;
					}
				}
				t          = a[npiv-1];
				a[npiv-1]  = a[kdiag-1];
				a[kdiag-1] = t;
//        c
//        c Test for singularity:
//        c
				if(Math.abs(a[kdiag-1])<tol) {
					System.err.println("Diagonal element "+(kdiag)+" is zero -> singular matrix");
					//ktilt=k;
					return null;
				}
//        c
//        c Compute multipliers:
//        c
				i1 = kdiag;
				for(int i=kp1; i<=n; i++) {
					i1++;
					a[i1-1] /= -a[kdiag-1];
				}
//        c
//        c Interchange and eliminate column per column:
//        c
				j1 = kdiag;
				j2 = npiv;
				for(int j=kp1; j<=n; j++) {
					j1 += n;
					j2 += n;
					t       = a[j2-1];
					a[j2-1] = a[j1-1];
					a[j1-1] = t;
					i1    = j1;
					i2    = kdiag;
					for(int i=kp1; i<=n; i++) {
						i1++;
						i2++;
						a[i1-1] += a[i2-1]*a[j1-1];
					}
				}
//        c
//        c Interchange and modify the ns right hand matrices:
//        c
				i1 = nvb+ipiv;
				i2 = nvb+k;
				for(int i=1; i<=ns; i++) {
					t       = b[i1-1];
					b[i1-1] = b[i2-1];
					b[i2-1] = t;
					j1    = i2;
					j2    = kdiag;
					for(int j=kp1; j<=n; j++) {
						j1++;
						j2++;
						b[j1-1] += b[i2-1]*a[j2-1];
					}
					i1 += n;
					i2 += n;
				}
			}
//        c
//        c Test for singularity for the last pivot:
//        c
			kdiag = ntn*iv;
			if(Math.abs(a[kdiag-1])<tol) {
				System.err.println("Last pivot is zero/too small -> assume singular matrix");
				//ktilt = n;
				return null;
			}
		}
//        c
//        c End of triangulation. Now, solve back variable per variable:
//        c
		for(int iv=1; iv<=nv; iv++) {
//        c
//        c Indices of location in vectors a and b:
//        c
			nva  = ntn*iv;
			//nvb1 = n*ns*(iv-1)+1;
			nvb2 = n*ns*iv;
//        c
//        c Back substitution with the ns right hand matrices:
//        c
			for(int il=1; il<=ns; il++) {
				for(int k=1; k<=nm1; k++) {
					nmk = n-k;
//        c
//        c Indice of the diagonal element of the (n-k+1)th row and of
//        c the (n-k+1)th element of the left hand side.
//        c
					kdiag = nva-(n+1)*(k-1);
					kb    = nvb2-(il-1)*n-k+1;
					b[kb-1] /= a[kdiag-1];
					t     = -b[kb-1];
					i1    = kb;
					i2    = kdiag;
					for(int i=1; i<=nmk; i++) {
						i1--;
						i2--;
						b[i1-1] += a[i2-1]*t;
					}
				}
				kdiag -= n+1;
				kb--;
				b[kb-1] /= a[kdiag-1];
			}
//        c
//        c End of back substitution:
//        c
		}
//        c
//        c Restitution of the solution:
//        c
		int itot = n*ns*nv;
		double[] x = new double[itot];
		for(int i=0; i<itot; i++) {
			x[i] = b[i];
		}
		return x;
	}
}
