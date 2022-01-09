package com.metzner.enrico.JKriging.helper;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DataHelper {
	
	public static int o_floor(double value) {
		return (int) value - (value<0d ? 1 : 0);
	}
	public static double o_mod(double value, double modul) {
		return value - modul * o_floor(value/modul);
	}

	public static class JulianDate {
		public static double now() {
			return System.currentTimeMillis()/86400000.0d + 2440587.5d;
		}
		public static double cal2jd(int[] _greg) {
			return cal2jd(_greg[0], _greg[1], _greg[2], _greg[3], _greg[4], _greg[5]);
		}
		public static double cal2jd(int _d, int _m, int _y, int _h, int _n, int _s) {
		    int a = (_m>2 ? _m : _m+12) + 1;
		    int j = (_m>2 ? _y : _y-1);
		    long dt = (long)(365.25d * (j+4716) + 0.001d)-(j<-4716?1L:0L) + o_floor(30.6001d*a) + _d;
		    dt += 2L - o_floor(j/100d+0.000001d) + o_floor(j/400d+0.000001d);
		    dt -= 1524L;
		    double time = (_h-12) / 24d + _n/1440d + _s/86400d;
		    return (dt + time);
		}
		public static int[] jd2cal(double _t) {
		    long z = o_floor(_t+0.500005787037037d);
		    int  x = o_floor( (z-1867216.25d) / 36524.25d );
		    long a = (z<2299161L ? z : z+1+x-o_floor(x/4d));
		    long b = a + 1524L;
		    int  c = o_floor( (b-122.1d)/365.25d );
		    long d = (long)( 365.25d * c ) - (c<0 ? 1L : 0L);
		    int  e = o_floor( (b-d)/30.6001d );
		    double f = _t+0.500005787037037d - z;
		    
		    int[] res = {0, 0, 0, 0, 0, 0};
		    res[0] = (int) ( b-d-o_floor(30.6001*e) ); //day
		    res[1] = (e<13.5 ? e-1 : e-13);           //mon
		    res[2] = (res[1]>2.5 ? c-4716 : c-4715); //year
		    res[3] = o_floor(24.0d*f);
		    res[4] = o_floor(1440.0d*f-60.0*res[3]);
		    res[5] = o_floor(86400.0d*f-3600.0d*res[3]-60.0d*res[4]);
		    return res;
		}
		public static String datestring() {
			return datestring(now());
		}
		public static String datestring(double _t) {
			return datestring(jd2cal(_t));
		}
		public static String datestring(int[] gregorian_date) {
			return gregorian_date[2]+"."+gregorian_date[1]+"."+gregorian_date[0]+" "+
					FormatHelper.nf(gregorian_date[3],2)+":"+FormatHelper.nf(gregorian_date[4],2,'0')+":"+FormatHelper.nf(gregorian_date[5],2,'0');
		}
	}
	
//	private static int[] index_sortem(double[] arr) {
//		int alen = arr.length;
//		int[] idx = new int[alen];
//		for(int i=0; i<alen; i++) idx[i] = -1;
//		for(int i=0; i<alen; i++) {
//			
//		}
//	}
	
	public static int[] createIndexArrayInt(int _length) {
		int[] idx = new int[_length];
		for(int i=0; i<_length; i++) idx[i] = i;
		return idx;
	}
	public static long[] createIndexArrayLong(int _length) {
		long[] idx = new long[_length];
		for(int i=0; i<_length; i++) idx[i] = i;
		return idx;
	}
	public static float[] createIndexArrayFloat(int _length) {
		float[] idx = new float[_length];
		for(int i=0; i<_length; i++) idx[i] = i;
		return idx;
	}
	public static double[] createIndexArrayDouble(int _length) {
		double[] idx = new double[_length];
		for(int i=0; i<_length; i++) idx[i] = i;
		return idx;
	}
	
	public static boolean[] concat_bool_array(boolean[] arr1, boolean[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		boolean[] naeloob = new boolean[a1len+a2len];
		for(int a=0; a<a1len; a++) naeloob[a]       = arr1[a];
		for(int a=0; a<a2len; a++) naeloob[a+a1len] = arr2[a];
		return naeloob;
	}
	public static byte[]    concat_byte_array(byte[] arr1, byte[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		byte[] etyb = new byte[a1len+a2len];
		for(int a=0; a<a1len; a++) etyb[a]       = arr1[a];
		for(int a=0; a<a2len; a++) etyb[a+a1len] = arr2[a];
		return etyb;
	}
	public static short[]   concat_short_array(short[] arr1, short[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		short[] trohs = new short[a1len+a2len];
		for(int a=0; a<a1len; a++) trohs[a]       = arr1[a];
		for(int a=0; a<a2len; a++) trohs[a+a1len] = arr2[a];
		return trohs;
	}
	public static int[]     concat_int_array(int[] arr1, int[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		int[] tni = new int[a1len+a2len];
		for(int a=0; a<a1len; a++) tni[a]       = arr1[a];
		for(int a=0; a<a2len; a++) tni[a+a1len] = arr2[a];
		return tni;
	}
	public static long[]    concat_long_array(long[] arr1, long[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		long[] gnol = new long[a1len+a2len];
		for(int a=0; a<a1len; a++) gnol[a]       = arr1[a];
		for(int a=0; a<a2len; a++) gnol[a+a1len] = arr2[a];
		return gnol;
	}
	public static float[]   concat_float_array(float[] arr1, float[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		float[] taolf = new float[a1len+a2len];
		for(int a=0; a<a1len; a++) taolf[a]       = arr1[a];
		for(int a=0; a<a2len; a++) taolf[a+a1len] = arr2[a];
		return taolf;
	}
	public static double[]  concat_double_array(double[] arr1, double[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		double[] elbuod = new double[a1len+a2len];
		for(int a=0; a<a1len; a++) elbuod[a]       = arr1[a];
		for(int a=0; a<a2len; a++) elbuod[a+a1len] = arr2[a];
		return elbuod;
	}
	public static String[]  concat_string_array(String[] arr1, String[] arr2) {
		int a1len=arr1.length, a2len=arr2.length;
		String[] gnirts = new String[a1len+a2len];
		for(int a=0; a<a1len; a++) gnirts[a]       = arr1[a];
		for(int a=0; a<a2len; a++) gnirts[a+a1len] = arr2[a];
		return gnirts;
	}

	public static boolean[][] concat_bool_array(boolean[][] arr1, boolean[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		boolean[][] naeloob = new boolean[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				naeloob[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) naeloob[a+a1len][b] = arr2[a][b];
				if(dim==1) naeloob[a][b+a1wid] = arr2[a][b];
			}
		return naeloob;
	}
	public static byte[][]    concat_byte_array(byte[][] arr1, byte[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		byte[][] etyb = new byte[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				etyb[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) etyb[a+a1len][b] = arr2[a][b];
				if(dim==1) etyb[a][b+a1wid] = arr2[a][b];
			}
		return etyb;
	}
	public static short[][]   concat_short_array(short[][] arr1, short[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		short[][] trohs = new short[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				trohs[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) trohs[a+a1len][b] = arr2[a][b];
				if(dim==1) trohs[a][b+a1wid] = arr2[a][b];
			}
		return trohs;
	}
	public static int[][]     concat_int_array(int[][] arr1, int[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		int[][] tni = new int[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				tni[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) tni[a+a1len][b] = arr2[a][b];
				if(dim==1) tni[a][b+a1wid] = arr2[a][b];
			}
		return tni;
	}
	public static long[][]    concat_long_array(long[][] arr1, long[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		long[][] gnol = new long[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				gnol[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) gnol[a+a1len][b] = arr2[a][b];
				if(dim==1) gnol[a][b+a1wid] = arr2[a][b];
			}
		return gnol;
	}
	public static float[][]   concat_float_array(float[][] arr1, float[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		float[][] taolf = new float[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				taolf[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) taolf[a+a1len][b] = arr2[a][b];
				if(dim==1) taolf[a][b+a1wid] = arr2[a][b];
			}
		return taolf;
	}
	public static double[][]  concat_double_array(double[][] arr1, double[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		double[][] elbuod = new double[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				elbuod[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) elbuod[a+a1len][b] = arr2[a][b];
				if(dim==1) elbuod[a][b+a1wid] = arr2[a][b];
			}
		return elbuod;
	}
	public static String[][]  concat_string_array(String[][] arr1, String[][] arr2, int dim) {
		int a1len=arr1.length, a2len=arr2.length;
		int a1wid=arr1[0].length, a2wid=arr2[0].length;
		String[][] gnirts = new String[dim==0?a1len+a2len:Math.min(a1len,a2len)][dim==1?Math.min(a1wid,a2wid):a1wid+a2wid];
		for(int a=0; a<a1len; a++)
			for(int b=0; b<a1wid; b++)
				gnirts[a][b]    = arr1[a][b];
		for(int a=0; a<a2len; a++)
			for(int b=0; b<a2wid; b++) {
				if(dim==0) gnirts[a+a1len][b] = arr2[a][b];
				if(dim==1) gnirts[a][b+a1wid] = arr2[a][b];
			}
		return gnirts;
	}

	public static int strings_index(String[] sarr, String search_regex) {
		if(sarr==null || search_regex==null) return -1;
		if(sarr.length<1) return -1;
		for(int s=0; s<sarr.length; s++) { if(sarr[s]==null) continue; if(sarr[s].equals(search_regex)) return s; }
		return -1;
	}
	

	public static void sortem(double[] pivot_arr, double[]... other_arr) {
		sortem(pivot_arr, 0, pivot_arr.length-1, false, other_arr);
	}
	public static void sortem(double[] pivot_arr, int begin_, int end_, boolean backward, double[]... other_arr) {
		List<int[]> partitions = new ArrayList<int[]>();
		partitions.add(new int[] {begin_, end_});
		//System.out.println("[DEBUG : SORTEM]    ");
		//int iterations = 0;
		while(!partitions.isEmpty()) {
			//iterations++;
			int[] partition = partitions.get(0);
			int begin = partition[0],
				end   = partition[1];
			//System.out.print("   ("+iterations+") "+(begin+1)+"|"+(end+1)+"|"+pivot_arr.length+" ");
//			if(begin>end) {
//				//throw new RuntimeException();
//				partitions.remove(0);
//				continue;
//			}
			//if(iterations%10==0) System.out.print("\n");
			
			if(begin>=end) { partitions.remove(0); continue; }
			
			if(end-begin==1) {
				if(pivot_arr[begin]>pivot_arr[end]) {
					swapElm(pivot_arr, begin, end);
					for(double[] arr: other_arr) swapElm(arr, begin, end);
				}
				partitions.remove(0);
				continue;
			}
			
			int i = begin,
				j = end;
			int k = (i+j)/2;
			if(pivot_arr[k]>pivot_arr[j]) {
				swapElm(pivot_arr, k, j);
				for(double[] arr: other_arr) swapElm(arr, k, j);
			}
			if(pivot_arr[i]>pivot_arr[j]) {
				swapElm(pivot_arr, i, j);
				for(double[] arr: other_arr) swapElm(arr, i, j);
			}
			if(pivot_arr[i]>pivot_arr[k]) {
				swapElm(pivot_arr, i, k);
				for(double[] arr: other_arr) swapElm(arr, i, k);
			}
			double pivot = pivot_arr[k];
			j++; i--;
			while(i<j) {
				do { j--; } while(pivot_arr[j]>pivot);
				do { i++; } while(pivot_arr[i]<pivot);
				if(i<j) {
					swapElm(pivot_arr, i, j);
					for(double[] arr: other_arr) swapElm(arr, i, j);
				}
			}
//			if(pivot_arr[i]<pivot) {
//				swapElm(pivot_arr, i, end);
//				for(double[] arr: other_arr) swapElm(arr, i, end);
//			}
			partitions.add(new int[] {begin, j});
			partitions.add(new int[] {j+1, end});
			partitions.remove(0);
			//sortem(pivot_arr, begin, i-1, backward, other_arr);
			//sortem(pivot_arr, i+1, end, backward, other_arr);
		}
		if(backward) {
			for(int i=0; i*2<end_-begin_; i++) {
				swapElm(pivot_arr, begin_+i, end_-i);
				for(double[] arr: other_arr) swapElm(arr, begin_+i, end_-i);
			}
		}
		//System.out.println("\n");
	}
	private static void swapElm(double[] darr, int left, int right) {
		double dtmp = darr[left];
		darr[left] = darr[right];
		darr[right] = dtmp;
	}
	
	
	
	
	
	public static int[] setsupr(int nx, double xmn, double xsiz, int ny, double ymn, double ysiz, int nz, double zmn, double zsiz,
			double[] x, double[] y, double[] z, double[] vr, double[] tmp, int nsec, double[] sec1, double[] sec2, double[] sec3,
			int maxsbx, int maxsby, int maxsbz, int[] nisb, double[] supblckgrid) {
//        c-----------------------------------------------------------------------
//        c
//        c           Establish Super Block Search Limits and Sort Data
//        c           *************************************************
//        c
//        c This subroutine sets up a 3-D "super block" model and orders the data
//        c by super block number.  The limits of the super block is set to the
//        c minimum and maximum limits of the grid; data outside are assigned to
//        c the nearest edge block.
//        c
//        c The idea is to establish a 3-D block network that contains all the
//        c relevant data.  The data are then sorted by their index location in
//        c the search network, i.e., the index location is given after knowing
//        c the block index in each coordinate direction (ix,iy,iz):
//        c          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
//        c An array, the same size as the number of super blocks, is constructed
//        c that contains the cumulative number of data in the model.  With this
//        c array it is easy to quickly check what data are located near any given
//        c location.
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   nx,xmn,xsiz      Definition of the X grid being considered
//        c   ny,ymn,ysiz      Definition of the Y grid being considered
//        c   nz,zmn,zsiz      Definition of the Z grid being considered
//        c   nd               Number of data
//        c   x(nd)            X coordinates of the data
//        c   y(nd)            Y coordinates of the data
//        c   z(nd)            Z coordinates of the data
//        c   vr(nd)           Variable at each location.
//        c   tmp(nd)          Temporary storage to keep track of the super block
//        c                      index associated to each data (uses the same
//        c                      storage already allocated for the simulation)
//        c   nsec             Number of secondary variables to carry with vr
//        c   sec1(nd)         First secondary variable (if nsec >= 1)
//        c   sec2(nd)         Second secondary variable (if nsec >= 2)
//        c   sec3(nd)         Third secondary variable (if nsec = 3)
//        c   MAXSB[X,Y,Z]     Maximum size of super block network
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   nisb()                Array with cumulative number of data in each
//        c                           super block.
//        c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
//        c   nysup,ymnsup,ysizsup  Definition of the Y super block grid
//        c   nzsup,zmnsup,zsizsup  Definition of the Z super block grid
//        c
//        c
//        c
//        c EXTERNAL REFERENCES:
//        c
//        c   sortem           Sorting routine to sort the data
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		boolean inflag=false;
//        c
//        c Establish the number and size of the super blocks:
//        c
		int    nxsup   = Math.min(nx,maxsbx);
		int    nysup   = Math.min(ny,maxsby);
		int    nzsup   = Math.min(nz,maxsbz);
		double xsizsup = nx*xsiz/nxsup;
		double ysizsup = ny*ysiz/nysup;
		double zsizsup = nz*zsiz/nzsup;
		double xmnsup  = (xmn-0.5d*xsiz)+0.5d*xsizsup;
		double ymnsup  = (ymn-0.5d*ysiz)+0.5d*ysizsup;
		double zmnsup  = (zmn-0.5d*zsiz)+0.5d*zsizsup;
		supblckgrid[0] = nxsup + 0.0001d; supblckgrid[1] = xmnsup; supblckgrid[2] = xsizsup;
		supblckgrid[3] = nysup + 0.0001d; supblckgrid[4] = ymnsup; supblckgrid[5] = ysizsup;
		supblckgrid[6] = nzsup + 0.0001d; supblckgrid[7] = zmnsup; supblckgrid[8] = zsizsup;
//        c
//        c Initialize the extra super block array to zeros:
//        c
		for(int i=0; i<nxsup*nysup*nzsup; i++) {
			nisb[i] = 0;
		}
//        c
//        c Loop over all the data assigning the data to a super block and
//        c accumulating how many data are in each super block:
//        c
		int nd = vr.length;
		for(int i=0; i<nd; i++) {
			int ix = getindx(nxsup, xmnsup, xsizsup, x[i], inflag);
			int iy = getindx(nysup, ymnsup, ysizsup, y[i], inflag);
			int iz = getindx(nzsup, zmnsup, zsizsup, z[i], inflag);
			//call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
			//call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
			//call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
			int ii = ix-1 + (iy-1)*nxsup + (iz-1)*nxsup*nysup;
			tmp[i] = ii;
			nisb[ii]++;
		}
//        c
//        c Sort the data by ascending super block number:
//        c
		//int nsort = 4 + nsec;
		sortem(tmp, x,y,z,vr, sec1,sec2,sec3);
		//call sortem(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)
//        c
//        c Set up array nisb with the starting address of the block data:
//        c
		for(int i=0; i<nxsup*nysup*nzsup-1; i++) {
			nisb[i+1] += nisb[i];
		}
//        c
//        c Finished:
//        c
		    	  return nisb;
	}
	public static int  picksupr(double[] supblckgrid, int irot, double[][][] rotmat, double radsqd,
			int[] ixsbtosr, int[] iysbtosr, int[] izsbtosr) {
		//input
		int nxsup = (int) supblckgrid[0],
			nysup = (int) supblckgrid[3],
			nzsup = (int) supblckgrid[6];
		double xsizsup = supblckgrid[2],
			   ysizsup = supblckgrid[5],
			   zsizsup = supblckgrid[8];
		//System.out.println("[DEBUG] SupBlckStratg:  x["+nxsup+"|"+supblckgrid[1]+"|"+xsizsup+
        //       	                                "]  y["+nysup+"|"+supblckgrid[4]+"|"+ysizsup+
        //       	                                "]  z["+nzsup+"|"+supblckgrid[7]+"|"+zsizsup+"]");
		//System.out.println("[DEBUG] i[xyz]sbtosr.length = "+ixsbtosr.length+" "+iysbtosr.length+" "+izsbtosr.length);
//        c-----------------------------------------------------------------------
//        c
//        c             Establish Which Super Blocks to Search
//        c             **************************************
//        c
//        c This subroutine establishes which super blocks must be searched given
//        c that a point being estimated/simulated falls within a super block
//        c centered at 0,0,0.
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   nxsup,xsizsup    Definition of the X super block grid
//        c   nysup,ysizsup    Definition of the Y super block grid
//        c   nzsup,zsizsup    Definition of the Z super block grid
//        c   irot             index of the rotation matrix for searching
//        c   MAXROT           size of rotation matrix arrays
//        c   rotmat           rotation matrices
//        c   radsqd           squared search radius
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   nsbtosr          Number of super blocks to search
//        c   ixsbtosr         X offsets for super blocks to search
//        c   iysbtosr         Y offsets for super blocks to search
//        c   izsbtosr         Z offsets for super blocks to search
//        c
//        c
//        c
//        c EXTERNAL REFERENCES:
//        c
//        c   sqdist           Computes anisotropic squared distance
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    real*8  rotmat(MAXROT,3,3),hsqd,sqdist,shortest
		//    integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
		double hsqd,shortest;
//        c
//        c MAIN Loop over all possible super blocks:
//        c
		int nsbtosr = 0;
		for(int i=1-nxsup; i<=nxsup-1; i++) for(int j=1-nysup; j<=nysup-1; j++) for(int k=1-nzsup; k<=nzsup-1; k++) {
			double xo = i*xsizsup;
			double yo = j*ysizsup;
			double zo = k*zsizsup;
			//System.out.println("[TEST] picksupr:  i,j,k = "+i+" "+j+" "+k);
//        c
//        c Find the closest distance between the corners of the super blocks:
//        c
			shortest = 1.0e21d;
			for(int i1=-1; i1<=1; i1++) for(int j1=-1; j1<=1; j1++) for(int k1=-1; k1<=1; k1++) {
				//System.out.println("[TEST] picksupr:    i1,j1,k1 = "+i1+" "+j1+" "+k1);
				for(int i2=-1; i2<=1; i2++) for(int j2=-1; j2<=1; j2++) for(int k2=-1; k2<=1; k2++) {
					//System.out.println("[TEST] picksupr:            i2,j2,k2 = "+i2+" "+j2+" "+k2);
					if(i1!=0 && j1!=0 && k1!=0 && i2!=0 && j2!=0 && k2!=0) {
						double xdis = (i1-i2)*0.5d*xsizsup + xo;
						double ydis = (j1-j2)*0.5d*ysizsup + yo;
						double zdis = (k1-k2)*0.5d*zsizsup + zo;
						hsqd = MathHelper.sqdist3D(0.0,0.0,0.0,xdis,ydis,zdis,rotmat[irot-1]);
						if(hsqd<=shortest) shortest = hsqd;
					}
				}
			}
			//System.out.println("[DEBUG] picksupr: shortest = "+shortest+"    radsqd = "+radsqd);
//        c
//        c Keep this super block if it is close enough:
//        c
			if(shortest <= radsqd) {
				//nsbtosr++;
				ixsbtosr[nsbtosr] = i;
				iysbtosr[nsbtosr] = j;
				izsbtosr[nsbtosr] = k;
				nsbtosr++; //copied afterwards for JAVA indices
			}
		}
//        c
//        c Finished:
//        c
		return nsbtosr;
	}
	public static int[] srchsupr(double xloc, double yloc, double zloc, double radsqd, int irot, double[][][] rotmat,
			int nsbtosr, int[] ixsbtosr, int[] iysbtosr, int[] izsbtosr, int noct, int nd,
			double[] x, double[] y, double[] z, double[] tmp, int[] nisb, double[] supblckgrid, double[] close) {
		int nxsup = (int) supblckgrid[0], nysup = (int) supblckgrid[3], nzsup = (int) supblckgrid[6];
		double xmnsup = supblckgrid[1], xsizsup = supblckgrid[2],
			   ymnsup = supblckgrid[4], ysizsup = supblckgrid[5],
			   zmnsup = supblckgrid[7], zsizsup = supblckgrid[8];
//        c-----------------------------------------------------------------------
//        c
//        c              Search Within Super Block Search Limits
//        c              ***************************************
//        c
//        c
//        c This subroutine searches through all the data that have been tagged in
//        c the super block subroutine.  The close data are passed back in the
//        c index array "close".  An octant search is allowed.
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   xloc,yloc,zloc   location of point being estimated/simulated
//        c   radsqd           squared search radius
//        c   irot             index of the rotation matrix for searching
//        c   MAXROT           size of rotation matrix arrays
//        c   rotmat           rotation matrices
//        c   nsbtosr          Number of super blocks to search
//        c   ixsbtosr         X offsets for super blocks to search
//        c   iysbtosr         Y offsets for super blocks to search
//        c   izsbtosr         Z offsets for super blocks to search
//        c   noct             If >0 then data will be partitioned into octants
//        c   nd               Number of data
//        c   x(nd)            X coordinates of the data
//        c   y(nd)            Y coordinates of the data
//        c   z(nd)            Z coordinates of the data
//        c   tmp(nd)          Temporary storage to keep track of the squared
//        c                      distance associated with each data
//        c   nisb()                Array with cumulative number of data in each
//        c                           super block.
//        c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
//        c   nysup,ymnsup,ysizsup  Definition of the X super block grid
//        c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   nclose           Number of close data
//        c   close()          Index of close data
//        c   infoct           Number of informed octants (only computes if
//        c                      performing an octant search)
//        c
//        c
//        c
//        c EXTERNAL REFERENCES:
//        c
//        c   sqdist           Computes anisotropic squared distance
//        c   sortem           Sorts multiple arrays in ascending order
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    real    x(*),y(*),z(*),tmp(*),close(*)
		//    real*8  rotmat(MAXROT,3,3),hsqd,sqdist
		double hsqd;
		//    integer nisb(*),inoct(8)
		int[] inoct = new int[8];
		//    integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
		//    logical inflag
		boolean inflag=true;
//        c
//        c Determine the super block location of point being estimated:
//        c
		int ix = getindx(nxsup, xmnsup, xsizsup, xloc, inflag);
		int iy = getindx(nysup, ymnsup, ysizsup, yloc, inflag);
		int iz = getindx(nzsup, zmnsup, zsizsup, zloc, inflag);
		//      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
		//      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
		//      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
//        c
//        c Loop over all the possible Super Blocks:
//        c
		int nclose = 0, i;
		for(int isup=0; isup<nsbtosr; isup++) {
//        c
//        c Is this super block within the grid system:
//        c
			int ixsup = ix + ixsbtosr[isup];
			int iysup = iy + iysbtosr[isup];
			int izsup = iz + izsbtosr[isup];
			if(ixsup<=0 || ixsup>nxsup || iysup<=0 || iysup>nysup || izsup<=0 || izsup>nzsup) continue;
//        c
//        c Figure out how many samples in this super block:
//        c
			int nums;
			int ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup;
			if(ii==1) {
				nums = nisb[ii-1];
				i    = 0;
			} else {
				nums = nisb[ii-1] - nisb[ii-2];
				i    = nisb[ii-2];
			}
			i--; //reduce by 1 for JAVA indices;
//        c
//        c Loop over all the data in this super block:
//        c
			for(ii=1; ii<=nums; ii++) {
				i++;
//        c
//        c Check squared distance:
//        c
				hsqd = MathHelper.sqdist3D(xloc,yloc,zloc,x[i],y[i],z[i],rotmat[irot-1]);
				if(hsqd>radsqd) continue;
//        c
//        c Accept this sample:
//        c
				//nclose++;
				close[nclose] = i+0.1d;
				tmp[nclose]   = hsqd;
				nclose++; //copied afterwards for JAVA indices
			}
		}
//        c
//        c Sort the nearby samples by distance to point being estimated:
//        c
		sortem(tmp,0,nclose-1,false,close);
		//System.out.println("[DEBUG] srchsupr:");
		//FormatHelper.printTable(nclose, tmp, close);
		//      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
//        c
//        c If we aren't doing an octant search then just return:
//        c
		if(noct<=0) return new int[]{nclose, 0};
//        c
//        c PARTITION THE DATA INTO OCTANTS:
//        c
		for(i=0; i<8; i++) {
			inoct[i] = 0;
		}
//        c
//        c Now pick up the closest samples in each octant:
//        c
		int nt = 8*noct;
		int na = 0, iq;
		for(int j=0; j<nclose; j++) {
			i  = (int) close[j]-1;
			double h  = tmp[j];
			double dx = x[i] - xloc;
			double dy = y[i] - yloc;
			double dz = z[i] - zloc;
			if(dz>=0d) {
				iq=4;
				if(dx<=0d && dy>0d) iq=1;
				if(dx>0d && dy>=0d) iq=2;
				if(dx<0d && dy>=0d) iq=3;
			} else {
				iq=8;
				if(dx<=0d && dy>0d) iq=5;
				if(dx>0d && dy>=0d) iq=6;
				if(dx<0d && dy<=0d) iq=7;
			}
			inoct[iq-1]++;
//        c
//        c Keep this sample if the maximum has not been exceeded:
//        c
			if(inoct[iq-1]<=noct) {
				//na = na + 1
				close[na] = i;
				tmp[na]   = h;
				na++; // copied afterwards for JAVA indices
				if(na==nt) break;
			}
		}
//        c
//        c End of data selection. Compute number of informed octants and return:
//        c
		nclose = na;
		int infoct = 0;
		for(i=0; i<8; i++) {
			if(inoct[i]>0) infoct++;
		}
//        c
//        c Finished:
//        c
		return new int[] {nclose, infoct};
	}
	
	
	
	
	
	
	
	
	public static int getindx(int n, double mn, double siz, double loc, boolean inflag) {
//        c-----------------------------------------------------------------------
//        c
//        c     Gets the coordinate index location of a point within a grid
//        c     ***********************************************************
//        c
//        c
//        c n       number of "nodes" or "cells" in this coordinate direction
//        c min     origin at the center of the first cell
//        c siz     size of the cells
//        c loc     location of the point being considered
//        c index   output index within [1,n]
//        c inflag  true if the location is actually in the grid (false otherwise
//        c         e.g., if the location is outside then index will be set to
//        c         nearest boundary
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    integer   n,index
		//    real      min,siz,loc
		//    logical   inflag
//        c
//        c Compute the index of "loc":
//        c
		int index = (int)  ( (loc-mn)/siz + 1.5d );
//        c
//        c Check to see if in or out:
//        c
		if(index<1) {
			index  = 1;
			inflag = false;
		} else if(index>n) {
			index  = n;
			inflag = false;
		} else {
			inflag = true;
		}
//        c
//        c Return to calling program:
//        c
		return index;
	}

	public static int double_index(double[] darr, double search_regex) {
		if(darr==null || Double.isNaN(search_regex)) return -1;
		if(darr.length<1) return -1;
		for(int d=0; d<darr.length; d++) { if(Double.isNaN(darr[d])) continue; if(darr[d]==search_regex) return d; }
		return -1;
	}
	
	public static void printStackTrace(PrintStream stream) {
		String pst = Arrays.toString(Thread.currentThread().getStackTrace());
		pst = pst.substring(pst.indexOf(','),pst.length()).replaceAll(",", "\n");
		stream.println(pst);
	}
}
