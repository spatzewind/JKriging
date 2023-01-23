package com.metzner.enrico.JKriging.probability;

import com.metzner.enrico.JKriging.data.Constants;

public class StdAnalysis {

	public static byte[]   minmax(byte[] arr)   { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte b: arr) { if(b<in) in=b; if(b>ax) ax=b; }   return new byte[] {in,ax}; }
	public static short[]  minmax(short[] arr)  { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short s: arr) { if(s<in) in=s; if(s>ax) ax=s; }  return new short[] {in,ax}; }
	public static int[]    minmax(int[] arr)    { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int i: arr) { if(i<in) in=i; if(i>ax) ax=i; }    return new int[] {in,ax}; }
	public static long[]   minmax(long[] arr)   { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long l: arr) { if(l<in) in=l; if(l>ax) ax=l; }   return new long[] {in,ax}; }
	public static float[]  minmax(float[] arr)  { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float f: arr) { if(Float.isNaN(f)) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[] arr) { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double d: arr) { if(Double.isNaN(d)) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }
	public static byte[]   minmax(byte[][] arr)     { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte[] barr: arr) for(byte b: barr) { if(b<in) in=b; if(b>ax) ax=b; } return new byte[] {in,ax}; }
	public static short[]  minmax(short[][] arr)    { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short[] sarr: arr) for(short s: sarr) { if(s<in) in=s; if(s>ax) ax=s; } return new short[] {in,ax}; }
	public static int[]    minmax(int[][] arr)      { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int[] iarr: arr) for(int i: iarr) { if(i<in) in=i; if(i>ax) ax=i; } return new int[] {in,ax}; }
	public static long[]   minmax(long[][] arr)     { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long[] larr: arr) for(long l: larr) { if(l<in) in=l; if(l>ax) ax=l; } return new long[] {in,ax}; }
	public static float[]  minmax(float[][] arr)    { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float[] farr: arr) for(float f: farr) { if(Float.isNaN(f)) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[][] arr)   { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double[] darr: arr) for(double d: darr) { if(Double.isNaN(d)) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }
	public static byte[]   minmax(byte[][][] arr)     { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte[][] bbarr: arr) for(byte[] barr: bbarr) for(byte b: barr) { if(b<in) in=b; if(b>ax) ax=b; } return new byte[] {in,ax}; }
	public static short[]  minmax(short[][][] arr)    { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short[][] ssarr: arr) for(short[] sarr: ssarr) for(short s: sarr) { if(s<in) in=s; if(s>ax) ax=s; } return new short[] {in,ax}; }
	public static int[]    minmax(int[][][] arr)      { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int[][] iiarr: arr) for(int[] iarr: iiarr) for(int i: iarr) { if(i<in) in=i; if(i>ax) ax=i; } return new int[] {in,ax}; }
	public static long[]   minmax(long[][][] arr)     { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long[][] llarr: arr) for(long[] larr: llarr) for(long l: larr) { if(l<in) in=l; if(l>ax) ax=l; } return new long[] {in,ax}; }
	public static float[]  minmax(float[][][] arr)    { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float[][] ffarr: arr) for(float[] farr: ffarr) for(float f: farr) { if(Float.isNaN(f)) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[][][] arr)   { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double[][] ddarr: arr) for(double[] darr: ddarr) for(double d: darr) { if(Double.isNaN(d)) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }


	public static byte[]   minmax(byte[] arr, byte _fillvalue)     { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte b: arr) { if(b==_fillvalue) continue; if(b<in) in=b; if(b>ax) ax=b; } return new byte[] {in,ax}; }
	public static short[]  minmax(short[] arr, short _fillvalue)   { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short s: arr) { if(s==_fillvalue) continue; if(s<in) in=s; if(s>ax) ax=s; } return new short[] {in,ax}; }
	public static int[]    minmax(int[] arr, int _fillvalue)       { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int i: arr) { if(i==_fillvalue) continue; if(i<in) in=i; if(i>ax) ax=i; } return new int[] {in,ax}; }
	public static long[]   minmax(long[] arr, long _fillvalue)     { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long l: arr) { if(l==_fillvalue) continue; if(l<in) in=l; if(l>ax) ax=l; } return new long[] {in,ax}; }
	public static float[]  minmax(float[] arr, float _fillvalue)   { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float f: arr) { if(Float.isNaN(f) || f==_fillvalue) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[] arr, double _fillvalue) { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double d: arr) { if(Double.isNaN(d) || d==_fillvalue) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }
	public static byte[]   minmax(byte[][] arr, byte _fillvalue)     { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte[] barr: arr) for(byte b: barr)  { if(b==_fillvalue) continue; if(b<in) in=b; if(b>ax) ax=b; } return new byte[] {in,ax}; }
	public static short[]  minmax(short[][] arr, short _fillvalue)   { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short[] sarr: arr) for(short s: sarr)  { if(s==_fillvalue) continue; if(s<in) in=s; if(s>ax) ax=s; } return new short[] {in,ax}; }
	public static int[]    minmax(int[][] arr, int _fillvalue)       { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int[] iarr: arr) for(int i: iarr)  { if(i==_fillvalue) continue; if(i<in) in=i; if(i>ax) ax=i; } return new int[] {in,ax}; }
	public static long[]   minmax(long[][] arr, long _fillvalue)     { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long[] larr: arr) for(long l: larr)  { if(l==_fillvalue) continue; if(l<in) in=l; if(l>ax) ax=l; } return new long[] {in,ax}; }
	public static float[]  minmax(float[][] arr, float _fillvalue)   { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float[] farr: arr) for(float f: farr)  { if(Float.isNaN(f) || f==_fillvalue) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[][] arr, double _fillvalue) { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double[] darr: arr) for(double d: darr)  { if(Double.isNaN(d) || d==_fillvalue) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }
	public static byte[]   minmax(byte[][][] arr, byte _fillvalue)     { byte in=Byte.MAX_VALUE, ax=Byte.MIN_VALUE;
		for(byte[][] bbarr: arr) for(byte[] barr: bbarr) for(byte b: barr)  { if(b==_fillvalue) continue; if(b<in) in=b; if(b>ax) ax=b; } return new byte[] {in,ax}; }
	public static short[]  minmax(short[][][] arr, short _fillvalue)   { short in=Short.MAX_VALUE, ax=Short.MIN_VALUE;
		for(short[][] ssarr: arr) for(short[] sarr: ssarr) for(short s: sarr)  { if(s==_fillvalue) continue; if(s<in) in=s; if(s>ax) ax=s; } return new short[] {in,ax}; }
	public static int[]    minmax(int[][][] arr, int _fillvalue)       { int in=Integer.MAX_VALUE, ax=Integer.MIN_VALUE;
		for(int[][] iiarr: arr) for(int[] iarr: iiarr) for(int i: iarr)  { if(i==_fillvalue) continue; if(i<in) in=i; if(i>ax) ax=i; } return new int[] {in,ax}; }
	public static long[]   minmax(long[][][] arr, long _fillvalue)     { long in=Long.MAX_VALUE, ax=Long.MIN_VALUE;
		for(long[][] llarr: arr) for(long[] larr: llarr) for(long l: larr)  { if(l==_fillvalue) continue; if(l<in) in=l; if(l>ax) ax=l; } return new long[] {in,ax}; }
	public static float[]  minmax(float[][][] arr, float _fillvalue)   { float in=Float.POSITIVE_INFINITY, ax=Float.NEGATIVE_INFINITY;
		for(float[][] ffarr: arr) for(float[] farr: ffarr) for(float f: farr)  { if(Float.isNaN(f) || f==_fillvalue) continue; if(f<in) in=f; if(f>ax) ax=f; }
		if(in>ax) { in=Float.NaN; ax=Float.NaN; } return new float[] {in,ax}; }
	public static double[] minmax(double[][][] arr, double _fillvalue) { double in=Double.POSITIVE_INFINITY, ax=Double.NEGATIVE_INFINITY;
		for(double[][] ddarr: arr) for(double[] darr: ddarr) for(double d: darr)  { if(Double.isNaN(d) || d==_fillvalue) continue; if(d<in) in=d; if(d>ax) ax=d; }
		if(in>ax) { in=Double.NaN; ax=Double.NaN; } return new double[] {in,ax}; }


	public static byte[]   mean_var(byte[] arr) {
		if(arr.length==0) return new byte[] {Byte.MIN_VALUE, Byte.MIN_VALUE}; if(arr.length==1) return new byte[] {arr[0], Byte.MIN_VALUE};
		double m=0d; for(byte b: arr) m+=b; byte bm = (byte) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = bm;
		double v=0d; for(byte b: arr) v+=(m-b)*(m-b); byte bv = (byte) (v/(arr.length-1d)+0.5d); return new byte[] {bm, bv}; }
	public static short[]  mean_var(short[] arr) {
		if(arr.length==0) return new short[] {Short.MIN_VALUE, Short.MIN_VALUE}; if(arr.length==1) return new short[] {arr[0], Short.MIN_VALUE};
		double m=0d; for(short s: arr) m+=s; short sm = (short) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = sm;
		double v=0d; for(short s: arr) v+=(m-s)*(m-s); short sv = (short) (v/(arr.length-1d)+0.5d); return new short[] {sm, sv}; }
	public static int[]    mean_var(int[] arr) {
		if(arr.length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE}; if(arr.length==1) return new int[] {arr[0], Integer.MIN_VALUE};
		double m=0d; for(int i: arr) m+=i; int im = (int) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = im;
		double v=0d; for(int i: arr) v+=(m-i)*(m-i); int iv = (int) (v/(arr.length-1d)+0.5d); return new int[] {im, iv}; }
	public static long[]   mean_var(long[] arr) {
		if(arr.length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE}; if(arr.length==1) return new long[] {arr[0], Long.MIN_VALUE};
		double m=0d; for(long l: arr) m+=l; long lm = (long) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = lm;
		double v=0d; for(long l: arr) v+=(m-l)*(m-l); long lv = (long) (v/(arr.length-1d)+0.5d); return new long[] {lm, lv}; }
	public static float[]  mean_var(float[] arr) {
		if(arr.length==0) return new float[] {Float.NaN, Float.NaN};
		double m=0d, mc=0d; for(float f: arr) { if(Float.isNaN(f)) continue; m+=f; mc+=1d; } if(mc<0.5d) return new float[] {Float.NaN, Float.NaN};
		float fm = (float) (m/mc); m = fm; if(mc<1.5d) return new float[] {fm, Float.NaN}; double v=0d;
		for(float f: arr) { if(Float.isNaN(f)) continue; v+=(m-f)*(m-f); } float fv = (float) (v/(mc-1d)); return new float[] {fm, fv}; }
	public static double[] mean_var(double[] arr) {
		if(arr.length==0) return new double[] {Double.NaN, Double.NaN};
		double m=0d, mc=0d; for(double d: arr) { if(Double.isNaN(d)) continue; m+=d; mc+=1d; } if(mc<0.5d) return new double[] {Double.NaN, Double.NaN};
		m /= mc; if(mc<1.5d) return new double[] {m, Double.NaN}; double v=0d;
		for(double d: arr) { if(Double.isNaN(d)) continue; v+=(m-d)*(m-d); } v /= mc-1d; return new double[] {m, v}; }
	public static byte[] mean_var(byte[][] arr) {
		if(arr.length==0) return new byte[] {Byte.MIN_VALUE, Byte.MIN_VALUE}; if(arr[0].length==0) return new byte[] {Byte.MIN_VALUE, Byte.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1) return new byte[] {arr[0][0], Byte.MIN_VALUE};
		double m=0d; for(byte[] barr: arr) for(byte b: barr) m+=b; byte bm = (byte) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = bm;
		double v=0d; for(byte[] barr: arr) for(byte b: barr) v+=(m-b)*(m-b); byte bv = (byte) (v/(arr.length-1d)+0.5d); return new byte[] {bm, bv}; }
	public static short[] mean_var(short[][] arr) {
		if(arr.length==0) return new short[] {Short.MIN_VALUE, Short.MIN_VALUE}; if(arr[0].length==0) return new short[] {Short.MIN_VALUE, Short.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1) return new short[] {arr[0][0], Short.MIN_VALUE};
		double m=0d; for(short[] sarr: arr) for(short s: sarr) m+=s; short sm = (short) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = sm;
		double v=0d; for(short[] sarr: arr) for(short s: sarr) v+=(m-s)*(m-s); short sv = (short) (v/(arr.length-1d)+0.5d); return new short[] {sm, sv}; }
	public static int[] mean_var(int[][] arr) {
		if(arr.length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE}; if(arr[0].length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1) return new int[] {arr[0][0], Integer.MIN_VALUE};
		double m=0d; for(int[] iarr: arr) for(int i: iarr) m+=i; int im = (int) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = im;
		double v=0d; for(int[] iarr: arr) for(int i: iarr) v+=(m-i)*(m-i); int iv = (int) (v/(arr.length-1d)+0.5d); return new int[] {im, iv}; }
	public static long[] mean_var(long[][] arr) {
		if(arr.length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE}; if(arr[0].length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1) return new long[] {arr[0][0], Long.MIN_VALUE};
		double m=0d; for(long[] larr: arr) for(long l: larr) m+=l; long lm = (long) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = lm;
		double v=0d; for(long[] larr: arr) for(long l: larr) v+=(m-l)*(m-l); long lv = (long) (v/(arr.length-1d)+0.5d); return new long[] {lm, lv}; }
	public static float[] mean_var(float[][] arr) {
		if(arr.length==0) return new float[] {Float.NaN, Float.NaN}; if(arr[0].length==0) return new float[] {Float.NaN, Float.NaN};
		double m=0d, mc=0d; for(float[] farr: arr) for(float f: farr) { if(Float.isNaN(f)) continue; m+=f; mc+=1d; }
		if(mc<0.5d) return new float[] {Float.NaN, Float.NaN}; float fm = (float) (m/mc); m = fm; if(mc<1.5d) return new float[] {fm, Float.NaN};
		double v=0d; for(float[] farr: arr) for(float f: farr) { if(Float.isNaN(f)) continue; v+=(m-f)*(m-f); } float fv = (float) (v/(mc-1d));
		return new float[] {fm, fv}; }
	public static double[] mean_var(double[][] arr) {
		if(arr.length==0) return new double[] {Double.NaN, Double.NaN}; if(arr[0].length==0) return new double[] {Double.NaN, Double.NaN};
		double m=0d, mc=0d; for(double[] darr: arr) for(double d: darr) { if(Double.isNaN(d)) continue; m+=d; mc+=1d; }
		if(mc<0.5d) return new double[] {Double.NaN, Double.NaN}; m /= mc; if(mc<1.5d) return new double[] {m, Double.NaN};
		double v=0d; for(double[] darr: arr) for(double d: darr) { if(Double.isNaN(d)) continue; v+=(m-d)*(m-d); } v /= mc-1d;
		return new double[] {m, v}; }
	public static byte[] mean_var(byte[][][] arr) {
		if(arr.length==0) return new byte[] {Byte.MIN_VALUE, Byte.MIN_VALUE}; if(arr[0].length==0) return new byte[] {Byte.MIN_VALUE, Byte.MIN_VALUE};
		if(arr[0][0].length==0) return new byte[] {Byte.MIN_VALUE,Byte.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1 && arr[0][0].length==1) return new byte[] {arr[0][0][0], Byte.MIN_VALUE};
		double m=0d; for(byte[][] bbarr: arr) for(byte[] barr: bbarr) for(byte b: barr) m+=b;
		byte bm = (byte) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = bm; double v=0d;
		for(byte[][] bbarr: arr) for(byte[] barr: bbarr) for(byte b: barr) v+=(m-b)*(m-b); byte bv = (byte) (v/(arr.length-1d)+0.5d);
		return new byte[] {bm, bv}; }
	public static short[] mean_var(short[][][] arr) {
		if(arr.length==0) return new short[] {Short.MIN_VALUE, Short.MIN_VALUE}; if(arr[0].length==0) return new short[] {Short.MIN_VALUE, Short.MIN_VALUE};
		if(arr[0][0].length==0) return new short[] {Short.MIN_VALUE,Short.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1 && arr[0][0].length==1) return new short[] {arr[0][0][0], Short.MIN_VALUE};
		double m=0d; for(short[][] ssarr: arr) for(short[] sarr: ssarr) for(short s: sarr) m+=s;
		short sm = (short) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = sm; double v=0d;
		for(short[][] ssarr: arr) for(short[] sarr: ssarr) for(short s: sarr) v+=(m-s)*(m-s); short sv = (short) (v/(arr.length-1d)+0.5d);
		return new short[] {sm, sv}; }
	public static int[] mean_var(int[][][] arr) {
		if(arr.length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE}; if(arr[0].length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE};
		if(arr[0][0].length==0) return new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1 && arr[0][0].length==1) return new int[] {arr[0][0][0], Integer.MIN_VALUE};
		double m=0d; for(int[][] iiarr: arr) for(int[] iarr: iiarr) for(int i: iarr) m+=i;
		int im = (int) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = im; double v=0d;
		for(int[][] iiarr: arr) for(int[] iarr: iiarr) for(int i: iarr) v+=(m-i)*(m-i); int iv = (int) (v/(arr.length-1d)+0.5d);
		return new int[] {im, iv}; }
	public static long[] mean_var(long[][][] arr) {
		if(arr.length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE}; if(arr[0].length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE};
		if(arr[0][0].length==0) return new long[] {Long.MIN_VALUE, Long.MIN_VALUE};
		if(arr.length==1 && arr[0].length==1) return new long[] {arr[0][0][0], Long.MIN_VALUE};
		double m=0d; for(long[][] llarr: arr) for(long[] larr: llarr) for(long l: larr) m+=l;
		long lm = (long) (m/arr.length+(m<0d ? -0.5d : 0.5d)); m = lm; double v=0d;
		for(long[][] llarr: arr) for(long[] larr: llarr) for(long l: larr) v+=(m-l)*(m-l); long lv = (long) (v/(arr.length-1d)+0.5d);
		return new long[] {lm, lv}; }
	public static float[] mean_var(float[][][] arr) {
		if(arr.length==0) return new float[] {Float.NaN, Float.NaN}; if(arr[0].length==0) return new float[] {Float.NaN, Float.NaN};
		if(arr[0][0].length==0) return new float[] {Float.NaN, Float.NaN};
		double m=0d, mc=0d; for(float[][] ffarr: arr) for(float[] farr: ffarr) for(float f: farr) { if(Float.isNaN(f)) continue; m+=f; mc+=1d; }
		if(mc<0.5d) return new float[] {Float.NaN, Float.NaN}; float fm = (float) (m/mc); m = fm; if(mc<1.5d) return new float[] {fm, Float.NaN};
		double v=0d; for(float[][] ffarr: arr) for(float[] farr: ffarr) for(float f: farr) { if(Float.isNaN(f)) continue; v+=(m-f)*(m-f); }
		float fv = (float) (v/(mc-1d)); return new float[] {fm, fv}; }
	public static double[] mean_var(double[][][] arr) {
		if(arr.length==0) return new double[] {Double.NaN, Double.NaN}; if(arr[0].length==0) return new double[] {Double.NaN, Double.NaN};
		if(arr[0][0].length==0) return new double[] {Double.NaN, Double.NaN};
		double m=0d, mc=0d; for(double[][] ddarr: arr) for(double[] darr: ddarr) for(double d: darr) { if(Double.isNaN(d)) continue; m+=d; mc+=1d; }
		if(mc<0.5d) return new double[] {Double.NaN, Double.NaN}; m /= mc; if(mc<1.5d) return new double[] {m, Double.NaN};
		double v=0d; for(double[][] ddarr: arr) for(double[] darr: ddarr) for(double d: darr) { if(Double.isNaN(d)) continue; v+=(m-d)*(m-d); }
		v /= mc-1d; return new double[] {m, v}; }

	public static double[] mean_var(double[] arr, double[] wgt) {
		if(wgt==null || arr.length!=wgt.length) return mean_var(arr);
		if(arr.length==0) return new double[] {Double.NaN, Double.NaN};
		double m=0d, mw=0d, mc=0d; for(int i=0; i<arr.length; i++) { if(Double.isNaN(arr[i]) || Double.isNaN(wgt[i])) continue; m+=arr[i]*wgt[i]; mw+=wgt[i]; mc+=1d; }
		if(mc<0.5d || Math.abs(mw)<Constants.D_EPSLON) return new double[] {Double.NaN, Double.NaN};
		m /= mw; if(mc<1.5d) return new double[] {m, Double.NaN}; double v=0d;
		for(int i=0; i<arr.length; i++) { if(Double.isNaN(arr[i])||Double.isNaN(wgt[i])) continue; v+=(m-arr[i])*(m-arr[i])*wgt[i]*wgt[i]; mw+=wgt[i]*wgt[i]; }
		v /= mw*(mc-1d)/mc; return new double[] {m, v}; }

//	@SuppressWarnings("unchecked")
//	public static <T extends Number> T variance(T[] arr) {
//		double _variance = Double.NaN;
//		if(arr.length>1) {
//			double _mean = mean(arr).doubleValue();
//			_variance = 0d;
//			for (T a: arr) {
//				double v = a.doubleValue();
//				_variance += (v-_mean)*(v-_mean);
//			}
//			_variance /= arr.length-1;
//		}
//		T var = null;
//		try {
//			var = (T) ((Double) _variance);
//		} catch(ClassCastException cce) {
//		}
//		return var;
//	}
}
