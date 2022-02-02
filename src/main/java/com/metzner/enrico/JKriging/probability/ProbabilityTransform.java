package com.metzner.enrico.JKriging.probability;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.helper.DataHelper;
import com.metzner.enrico.JKriging.helper.Interpolation;

public class ProbabilityTransform {
	
	public static final int INTPOL_LINEAR     = 1;
	public static final int INTPOL_POWER      = 2;
	public static final int INTPOL_HYPERBOLIC = 4;

	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing dataframe
	 * @param variable_id      id of the variable in the dataframe
	 * @param weights_id       id of optional weight variable in the dataframe (if not use, set -1)
	 * @param name_transformed name of the transformed result to put it in the dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame df, int variable_id, int weights_id, double separation, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id),
				new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, separation, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by name "variable" var(i),i=1,...,n possibly weighted
	 *     by variable with name "weight_variable" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing dataframe
	 * @param variable         name of the variable in the dataframe
	 * @param weight_variable  name of optional weight variable in the dataframe (if not use, set null)
	 * @param name_transformed name of the transformed result to put it in the dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame df, String variable, String weight_variable, double separation, String name_transformed) {
		return nscore(df, variable, weight_variable, new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, separation, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing dataframe
	 * @param variable_id      id of the variable in the dataframe
	 * @param weights_id       id of optional weight variable in the dataframe (if not use, set -1)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the dataframe
	 * @return                 transformation table as array with var-transformed-pair
	 */
	public static double[][] nscore(DataFrame df, int variable_id, int weights_id, double[] trims, double separation, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id), trims, separation, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing dataframe
	 * @param variable         name of the variable in the dataframe
	 * @param weight_variable  name of optional weight variable in the dataframe (if not use, set null)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame df, String variable, String weight_variable, double[] trims, double separation, String name_transformed) {
//        c-----------------------------------------------------------------------
//        c
//        c              Transform Univariate Data to Normal Scores
//        c              ******************************************
//        c
//        c This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
//        c by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
//        c as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
//        c so that the data can be returned in the same order (just in case
//        c there are associated arrays like the coordinate location).
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   nd               Number of data (no missing values)
//        c   vr(nd)           Data values to be transformed
//        c   tmin,tmax        data trimming limits
//        c   iwt              =0, equal weighted; =1, then apply weight
//        c   wt(nd)           Weight for each data (don't have to sum to 1.0)
//        c   tmp(nd)          Temporary storage space for sorting
//        c   lout             if > 0 then transformation table will be written
//        c
//        c
//        c
//        c OUTPUT VARIABLES:
//        c
//        c   vrg(nd)          normal scores
//        c   ierror           error flag (0=error free,1=problem)
//        c
//        c
//        c
//        c EXTERNAL REFERENCES:
//        c
//        c   gauinv           Calculates the inverse of a Gaussian cdf
//        c   sortem           sorts a number of arrays according to a key array
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return null;
		}
		if(!df.hasVariable(variable)) {
			System.err.println("Variable for normal score transform does not exist!");
			return null;
		}
		double[] vr;
		switch(df.getVariableType(variable)) {
			case FLOAT: float[] fr = (float[]) df.getArray(variable);
				vr = new double[fr.length]; for(int fi=0; fi<fr.length; fi++) {
					vr[fi] = fr[fi]; if(Float.isNaN(fr[fi])) vr[fi] = Double.NaN; }
				break;
			case DOUBLE: vr = (double[]) df.getArray(variable);
				break;
			default:
				System.err.println("Normal score transformation expect datatype \"float\" or \"double\" but found "+df.getVariableType(variable).name());
				return null;
		}
		int nd = vr.length;
		double[] wt  = (double[]) df.getArray(weight_variable);
		if(weight_variable==null || wt==null) {
			wt = new double[nd];
			for(int w=0; w<nd; w++) wt[w] = 1d;
		} else if(nd!=wt.length) {
			System.err.println("No normal score transformation possible, weight has not the same length as the variable!");
			return null;
		}
		
		//remove redundancy and NaNs!
		int max_hit_count = 1;
		List<double[]> distribution = new ArrayList<double[]>();
		for(int i=0; i<nd; i++) {
			if(Double.isNaN(vr[i])) continue;
			if(vr[i]<trims[0] || vr[i]>=trims[1]) continue;
			boolean is_new = true;
			int entryID = 0;
			while(entryID<distribution.size()) {
				double[] entry = distribution.get(entryID);
				if(Math.abs(vr[i]-entry[0])<separation) {
					is_new = false;
					break;
				}
				entryID++;
			}
			if(is_new) {
				distribution.add(new double[] {vr[i], wt[i], i+0.1d});
			} else {
				double[] e_old = distribution.get(entryID);
				double[] e_new = new double[e_old.length+1];
				for(int e=0; e<e_old.length; e++)
					e_new[e] = e_old[e];
				e_new[e_old.length] = i+0.1d;
				distribution.remove(entryID);
				e_new[0] = ((e_old.length-2)*e_old[0] + vr[i]) / (e_new.length-2);
				e_new[1] += wt[i];
				distribution.add(e_new);
				max_hit_count = Math.max(max_hit_count, e_new.length-2);
			}
		}
		
		nd = distribution.size();
		System.out.println("[PROB-TRANSFORM] NSCORE[1d]:  nd="+nd+", hc="+max_hit_count);
		double[] dist = new double[nd];
		double[][] indices = new double[max_hit_count+2][nd];
		double twt = 0d;
		for(int e=0; e<nd; e++) {
			double[] entry = distribution.get(e);
			dist[e] = entry[0];
			for(int e1=0; e1<entry.length; e1++)
				indices[e1][e] = entry[e1];
			for(int e2=entry.length; e2<max_hit_count+2; e2++)
				indices[e2][e] = -1.1d;
			twt += indices[1][e];
		}
		
		if(nd<=0) {
			System.err.println("The variable only consists of NaNs or values outside the defined range, no transformation possible!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
//		double[] vr_  = new double[nd];
//		double[] wt_  = new double[nd];
//		double[] vrg_ = new double[nd],
//				 tmp  = new double[nd];
//		int    c   = -1;
//		double twt = 0d;
//		for(Double d: double_map.keySet()) {
//			c++; vr_[c] = d.doubleValue();
//			wt_[c] = double_map.get(d)[0];
//			twt += wt_[c];
//			tmp[c] = double_map.get(d)[1];
//		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(30, vr_, wt_, tmp, vrg_);

//        c Sort the data in ascending order and calculate total weight:
		if(twt<Constants.D_EPSLON) {
			System.err.println("The total weight is too small, so review and adjust weights! No transformation now.");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		System.out.println("[PROB-TRANSFORM] NSCORE[1d]:  sorting");
		DataHelper.sortem(dist,indices);

//        c Compute the cumulative probabilities:
		System.out.println("[PROB-TRANSFORM] NSCORE[1d]:  compute nscore");
		double oldcp = 0d,
			   cp    = 0d,
			   wttmp = 0d;
		for(int i=0; i<nd; i++) {
			cp += indices[1][i] / twt;
			wttmp = 0.5d * (cp + oldcp);
			oldcp = cp;
			indices[0][i] = Gaussian.cdf_inv(wttmp);
			//if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
		}

//        c create transformation table
		System.out.println("[PROB-TRANSFORM] NSCORE[1d]:  create transformation table");
		double[][] trnsf = new double[2][nd];
		for(int i=0; i<nd; i++) {
			trnsf[0][i] = dist[i];
			trnsf[1][i] = indices[0][i];
		}


//        c Get the arrays back in original order:
		

		System.out.println("[PROB-TRANSFORM] NSCORE[1d]:  insert new variable");
		double[] vrg = new double[vr.length];
		for(int i=0; i<vr.length; i++)
			vrg[i] = Double.NaN;
		for(int i=0; i<nd; i++) {
			double vr_ns = indices[0][i];
			for(int e=0; e<max_hit_count; e++) {
				int j = (int) indices[2+e][i];
				if(j>=0) vrg[j] = vr_ns;
			}
		}
//		DataHelper.sortem(tmp, vrg_,vr_,wt_);
//		double[] vrg = new double[vr.length];
//		int lastIdx=-1;
//		for(c=0; c<=nd; c++) {
//			int currentIdx = 0;
//			if(c<nd) currentIdx =(int) tmp[c];
//			if(c==nd) currentIdx = vr.length-1;
//			for(int i=lastIdx+1; i<=currentIdx; i++) {
//				if(Double.isNaN(vr[i])) { vrg[i] = Double.NaN; continue; }
//				if(vr[i]<trims[0] || vr[i]>=trims[1]) { vrg[i] = Double.NaN; continue; }
//				for(int cc=0; cc<=Math.min(c,vr_.length-1); cc++) {
//					if(vr[i]==vr_[cc]) { vrg[i] = wt_[cc]; break; }
//				}
//			}
//			lastIdx = currentIdx;
//		}
		df.addColumn(name_transformed, vrg);
		return trnsf;
	}

	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 2d-dataframe
	 * @param variable_id      id of the variable in the 2d-dataframe
	 * @param weights_id       id of optional weight variable in the 2d-dataframe (if not use, set -1)
	 * @param name_transformed name of the transformed result to put it in the 2d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame2D df, int variable_id, int weights_id, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id),
				new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by name "variable" var(i),i=1,...,n possibly weighted
	 *     by variable with name "weight_variable" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 2d-dataframe
	 * @param variable         name of the variable in the 2d-dataframe
	 * @param weight_variable  name of optional weight variable in the 2d-dataframe (if not use, set null)
	 * @param name_transformed name of the transformed result to put it in the 2d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame2D df, String variable, String weight_variable, String name_transformed) {
		return nscore(df, variable, weight_variable, new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 2d-dataframe
	 * @param variable_id      id of the variable in the 2d-dataframe
	 * @param weights_id       id of optional weight variable in the 2d-dataframe (if not use, set -1)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the 2d-dataframe
	 * @return                 transformation table as array with var-transformed-pair
	 */
	public static double[][] nscore(DataFrame2D df, int variable_id, int weights_id, double[] trims, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id), trims, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 2d-dataframe
	 * @param variable         name of the variable in the 2d-dataframe
	 * @param weight_variable  name of optional weight variable in the 2d-dataframe (if not use, set null)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the 2d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame2D df, String variable, String weight_variable, double[] trims, String name_transformed) {
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return null;
		}
		if(!df.hasVariable(variable)) {
			System.err.println("Variable for normal score transform does not exist!");
			return null;
		}
		double[][] vr;
		switch(df.getVariableType(variable)) {
			case FLOAT: float[][] fr = (float[][]) df.getArray(variable);
				vr = new double[fr.length][fr[0].length]; for(int fi1=0; fi1<fr.length; fi1++) for(int fi2=0; fi2<fr[0].length; fi2++) {
					vr[fi1][fi2] = fr[fi1][fi2]; if(Float.isNaN(fr[fi1][fi2])) vr[fi1][fi2] = Double.NaN; }
				break;
			case DOUBLE: vr = (double[][]) df.getArray(variable);
				break;
			default:
				System.err.println("Normal score transformation expect datatype \"float\" or \"double\" but found "+df.getVariableType(variable).name());
				return null;
		}
		int[] nd = {vr.length, vr[0].length};
		double[][] wt  = (double[][]) df.getArray(weight_variable);
		if(weight_variable==null) {
			wt = new double[nd[0]][nd[1]];
			for(int w1=0; w1<nd[0]; w1++) for(int w2=0; w2<nd[1]; w2++) wt[w1][w2] = 1d;
		} else if(nd[0]!=wt.length || nd[1]!=wt[0].length) {
			System.err.println("No normal score transformation possible, weight has not the same length as the variable!");
			return null;
		}
		
		//remove redundancy and NaNs!
		Map<Double,double[]> double_map = new HashMap<Double, double[]>();
		for(int j=0; j<nd[0]; j++) for(int i=0; i<nd[1]; i++) {
			if(Double.isNaN(vr[j][i])) continue;
			if(vr[j][i]<trims[0] || vr[j][i]>=trims[1]) continue;
			Double d = vr[j][i];
			if(double_map.containsKey(d)) {
				double[] nwt = double_map.get(d);
				nwt[0] += wt[j][i];
				double_map.put(d, nwt);
			} else {
				double_map.put(d, new double[] {wt[j][i], j*nd[1]+i+0.1d});
			}
		}
		int nc = double_map.keySet().size();
		if(nc<=0) {
			System.err.println("The variable only consists of NaNs or values outside the defined range, no transformation possible!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		double[] vr_  = new double[nc];
		double[] wt_  = new double[nc];
		double[] vrg_ = new double[nc],
				 tmp  = new double[nc];
		int    c   = -1;
		double twt = 0d;
		for(Double d: double_map.keySet()) {
			c++; vr_[c] = d.doubleValue();
			wt_[c] = double_map.get(d)[0];
			twt += wt_[c];
			tmp[c] = double_map.get(d)[1];
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(30, vr_, wt_, tmp, vrg_);

//        c Sort the data in ascending order and calculate total weight:
		if(twt<Constants.D_EPSLON) {
			System.err.println("The total weight is too small, so review and adjust weights! No transformation now.");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		DataHelper.sortem(vr_,wt_,tmp);

//        c Compute the cumulative probabilities:
		double oldcp = 0d,
			   cp    = 0d,
			   wttmp = 0d;
		for(int i=0; i<nc; i++) {
			cp += wt_[i] / twt;
			wttmp = 0.5d * (cp + oldcp);
			oldcp = cp;
			vrg_[i] = Gaussian.cdf_inv(wttmp);
			//if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
		}

//        c create transformation table
		double[][] trnsf = new double[2][nc];
		for(int i=0; i<nc; i++) {
			trnsf[0][i] = vr_[i];
			trnsf[1][i] = vrg_[i];
		}

//        c Get the arrays back in original order:
		double[][] vrg = new double[nd[0]][nd[1]];
		for(int j=0; j<nd[0]; j++) for(int i=0; i<nd[1]; i++) {
			if(Double.isNaN(vr[j][i])) { vrg[j][i] = Double.NaN; continue; }
			if(vr[j][i]<trims[0] || vr[j][i]>=trims[1]) { vrg[j][i] = Double.NaN; continue; }
			int indBot = 0,
				indTop = vr_.length-1;
			if(vr[j][i]==vr_[indBot]) { vrg[j][i] = vrg_[indBot]; continue; }
			if(vr[j][i]==vr_[indTop]) { vrg[j][i] = vrg_[indTop]; continue; }
			while(indBot!=indTop) {
				int indMid = (indTop+indBot)/2;
				if(vr[j][i]==vr_[indMid]) {
					indBot = indMid; indTop = indMid;
				} else if(vr[j][i]<vr_[indMid]) {
					indTop = indMid;
				} else {
					indBot = indMid;
				}
			}
			vrg[j][i] = vrg_[indBot];
		}
		df.addColumn(name_transformed, vrg);
		return trnsf;
	}

	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 3d-dataframe
	 * @param variable_id      id of the variable in the 3d-dataframe
	 * @param weights_id       id of optional weight variable in the 3d-dataframe (if not use, set -1)
	 * @param name_transformed name of the transformed result to put it in the 3d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame3D df, int variable_id, int weights_id, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id),
				new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by name "variable" var(i),i=1,...,n possibly weighted
	 *     by variable with name "weight_variable" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 3d-dataframe
	 * @param variable         name of the variable in the 3d-dataframe
	 * @param weight_variable  name of optional weight variable in the 3d-dataframe (if not use, set null)
	 * @param name_transformed name of the transformed result to put it in the 3d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame3D df, String variable, String weight_variable, String name_transformed) {
		return nscore(df, variable, weight_variable, new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 3d-dataframe
	 * @param variable_id      id of the variable in the 3d-dataframe
	 * @param weights_id       id of optional weight variable in the 3d-dataframe (if not use, set -1)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the 3d-dataframe
	 * @return                 transformation table as array with var-transformed-pair
	 */
	public static double[][] nscore(DataFrame3D df, int variable_id, int weights_id, double[] trims, String name_transformed) {
		return nscore(df, df.getVarname(variable_id), df.getVarname(weights_id), trims, name_transformed);
	}
	/**
	 * Transform Univariate Data to Normal Scores
	 * <p>
	 *     This subroutine takes data from variable by "variable_id" var(i),i=1,...,n possibly weighted
	 *     by variable with "weights_id" wgt(i),i=,...,n and returns the normal scores transform N(0,1)
	 *     as "double[][]" vrg(i),i=1,...,n. An extra storage array "tmp" is used internally
	 *     so that the data can be returned in the same order (just in case there are associated arrays
	 *     like the coordinate location).
	 * </p>
	 * 
	 * @param df               data containing 3d-dataframe
	 * @param variable         name of the variable in the 3d-dataframe
	 * @param weight_variable  name of optional weight variable in the 3d-dataframe (if not use, set null)
	 * @param trims            double array with lower and upper bound of valid range
	 * @param name_transformed name of the transformed result to put it in the 3d-dataframe
	 * @return                 transformation table as Array with var-transformed-pairs
	 */
	public static double[][] nscore(DataFrame3D df, String variable, String weight_variable, double[] trims, String name_transformed) {
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return null;
		}
		if(!df.hasVariable(variable)) {
			System.err.println("Variable for normal score transform does not exist!");
			return null;
		}
		double[][][] vr;
		switch(df.getVariableType(variable)) {
			case FLOAT: float[][][] fr = (float[][][]) df.getArray(variable);
				vr = new double[fr.length][fr[0].length][fr[0][0].length];
				for(int fi1=0; fi1<fr.length; fi1++) for(int fi2=0; fi2<fr[0].length; fi2++) for(int fi3=0; fi3<fr[0][0].length; fi3++) {
					vr[fi1][fi2][fi3] = fr[fi1][fi2][fi3]; if(Float.isNaN(fr[fi1][fi2][fi3])) vr[fi1][fi2][fi3] = Double.NaN; }
				break;
			case DOUBLE: vr = (double[][][]) df.getArray(variable);
				break;
			default:
				System.err.println("Normal score transformation expect datatype \"float\" or \"double\" but found "+df.getVariableType(variable).name());
				return null;
		}
		int[] nd = {vr.length, vr[0].length, vr[0][0].length};
		double[][][] wt  = (double[][][]) df.getArray(weight_variable);
		if(weight_variable==null) {
			wt = new double[nd[0]][nd[1]][nd[2]];
			for(int w1=0; w1<nd[0]; w1++) for(int w2=0; w2<nd[1]; w2++) for(int w3=0; w3<nd[2]; w3++) wt[w1][w2][w3] = 1d;
		} else if(nd[0]!=wt.length || nd[1]!=wt[0].length || nd[2]!=wt[0][0].length) {
			System.err.println("No normal score transformation possible, weight has not the same length as the variable!");
			return null;
		}
		
		//remove redundancy and NaNs!
		Map<Double,double[]> double_map = new HashMap<Double, double[]>();
		for(int k=0; k<nd[0]; k++) for(int j=0; j<nd[1]; j++) for(int i=0; i<nd[2]; i++) {
			if(Double.isNaN(vr[k][j][i])) continue;
			if(vr[k][j][i]<trims[0] || vr[k][j][i]>=trims[1]) continue;
			Double d = vr[k][j][i];
			if(double_map.containsKey(d)) {
				double[] nwt = double_map.get(d);
				nwt[0] += wt[k][j][i];
				double_map.put(d, nwt);
			} else {
				double_map.put(d, new double[] {wt[k][j][i], k*nd[1]*nd[2]+j*nd[2]+i+0.1d});
			}
		}
		int nc = double_map.keySet().size();
		if(nc<=0) {
			System.err.println("The variable only consists of NaNs or values outside the defined range, no transformation possible!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		double[] vr_  = new double[nc];
		double[] wt_  = new double[nc];
		double[] vrg_ = new double[nc],
				 tmp  = new double[nc];
		int    c   = -1;
		double twt = 0d;
		for(Double d: double_map.keySet()) {
			c++; vr_[c] = d.doubleValue();
			wt_[c] = double_map.get(d)[0];
			twt += wt_[c];
			tmp[c] = double_map.get(d)[1];
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(30, vr_, wt_, tmp, vrg_);

//        c Sort the data in ascending order and calculate total weight:
		if(twt<Constants.D_EPSLON) {
			System.err.println("The total weight is too small, so review and adjust weights! No transformation now.");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		DataHelper.sortem(vr_,wt_,tmp);

//        c Compute the cumulative probabilities:
		double oldcp = 0d,
			   cp    = 0d,
			   wttmp = 0d;
		for(int i=0; i<nc; i++) {
			cp += wt_[i] / twt;
			wttmp = 0.5d * (cp + oldcp);
			oldcp = cp;
			vrg_[i] = Gaussian.cdf_inv(wttmp);
			//if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
		}

//        c create transformation table
		double[][] trnsf = new double[2][nc];
		for(int i=0; i<nc; i++) {
			trnsf[0][i] = vr_[i];
			trnsf[1][i] = vrg_[i];
		}

//        c Get the arrays back in original order:
		double[][][] vrg = new double[nd[0]][nd[1]][nd[2]];
		for(int k=0; k<nd[0]; k++) for(int j=0; j<nd[1]; j++) for(int i=0; i<nd[2]; i++) {
			if(Double.isNaN(vr[k][j][i])) { vrg[k][j][i] = Double.NaN; continue; }
			if(vr[k][j][i]<trims[0] || vr[k][j][i]>=trims[1]) { vrg[k][j][i] = Double.NaN; continue; }
			int indBot = 0,
				indTop = vr_.length-1;
			if(vr[k][j][i]==vr_[indBot]) { vrg[k][j][i] = vrg_[indBot]; continue; }
			if(vr[k][j][i]==vr_[indTop]) { vrg[k][j][i] = vrg_[indTop]; continue; }
			while(indBot!=indTop) {
				int indMid = (indTop+indBot)/2;
				if(vr[k][j][i]==vr_[indMid]) {
					indBot = indMid; indTop = indMid;
				} else if(vr[k][j][i]<vr_[indMid]) {
					indTop = indMid;
				} else {
					indBot = indMid;
				}
			}
			vrg[k][j][i] = vrg_[indBot];
		}
		df.addColumn(name_transformed, vrg);
		return trnsf;
	}
	
	
	
	
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing dataframe
	 * @param variable_id            id of the variable in the dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame df, int variable_id, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing dataframe
	 * @param variable               name of the variable in the dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame df, String variable, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, variable, new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing dataframe
	 * @param variable_id            id of the variable in the dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame df, int variable_id, double[] trims,
			double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), trims, nscore_trnsf_table, name_transformed,
				    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing dataframe
	 * @param variable               name of the variable in the dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame df, String variable, double[] trims, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
//        c-----------------------------------------------------------------------
//        c
//        c           Back Transform Univariate Data from Normal Scores
//        c           *************************************************
//        c
//        c This subroutine backtransforms a standard normal deviate from a
//        c specified back transform table and option for the tails of the
//        c distribution.  Call once with "first" set to true then set to false
//        c unless one of the options for the tail changes.
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   vrgs             normal score value to be back transformed
//        c   nt               number of values in the back transform tbale
//        c   vr(nt)           original data values that were transformed
//        c   vrg(nt)          the corresponding transformed values
//        c   zmin,zmax        limits possibly used for linear or power model
//        c   ltail            option to handle values less than vrg(1):
//        c   ltpar            parameter required for option ltail
//        c   utail            option to handle values greater than vrg(nt):
//        c   utpar            parameter required for option utail
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return;
		}
		if(variable==null) {
			System.err.println("Variable for normal score transform does not exist!");
			return;
		}
		if(nscore_trnsf_table==null) {
			System.err.println("Normal score transformation table does not exist!\n"+
					           "    no back transformation done!");
			return;
		}
		if(lower_tail_intpol_type!=1 && lower_tail_intpol_type!=2) {
			System.err.println("For lower tail only interpolation types 1 and 2 possible!\n"+
					           "    Found "+lower_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		if(upper_tail_intpol_type!=1 && upper_tail_intpol_type!=2 && upper_tail_intpol_type!=4) {
			System.err.println("For upper tail only interpolation types 1,2 and 4 possible!\n"+
					           "    Found "+upper_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		
		
		double[] val = (double[]) df.getArray(variable);
		int nd = val.length;
		double[] res = new double[nd];
		double lambda;
		
		double[] vrg = nscore_trnsf_table[1];
		double[] vr  = nscore_trnsf_table[0];
		int btlen = vrg.length;
		
		if(upper_tail_intpol_type==4 && vr[btlen-1]<0d) {
			System.err.println("ERROR can not use hyperbolic tail (intpol type 4) for upper tail with negative values! - see manual\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[0]>vr[0]) {
			System.err.println("ERROR zmin should be lower or equal than the first entry in the transformation table\n" + 
					           "    zmin = "+trims[0]+" vr_first="+vr[0]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[1]<vr[btlen-1]) {
			System.err.println("ERROR zmax should at higher or equal than the last entry in the transformation table\n" + 
					           "    zmax = "+trims[1]+" vr_last="+vr[btlen-1]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		
		for(int i=0; i<nd; i++) {
//        i added by Enrico Metzner
//        i check whether data is fill_value and continue or it can be transformed
			if(val[i]==Constants.FILL_VALUE_D) {
				res[i] = Constants.FILL_VALUE_D;
				continue;
			}
//        c
//        c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
//        c
			double vrgs = val[i];
			double backtr;
			if(vrgs < vrg[0]) {
				backtr = vr[0];
				double cdflo  = Gaussian.cdf(vrg[0]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(lower_tail_intpol_type==1) {
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,1.0);
				} else if(lower_tail_intpol_type==2) {
					double cpow   = 1d / lower_tail_power;
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,cpow);
				}
//        c
//        c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
//        c
			} else if(vrgs >= vrg[btlen-1]) {
				backtr = vr[btlen-1];
				double cdfhi  = Gaussian.cdf(vrg[btlen-1]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(upper_tail_intpol_type==1) {
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,1d);
				} else if(upper_tail_intpol_type==2) {
					double cpow   = 1d / upper_tail_power;
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,cpow);
				} else if(upper_tail_intpol_type==4) {
					lambda = Math.pow(vr[btlen-1],upper_tail_power)*(1.0-Gaussian.cdf(vrg[btlen-1]));
					backtr = Math.pow(lambda/(1.0-Gaussian.cdf(vrgs)),1d/upper_tail_power);
				}
//        c
//        c Value within the transformation table:
//        c
			} else {
				int p;
				for(p=0; p<btlen-1; p++) if(vrgs<vrg[p+1]) break;
				p = Math.max(0, Math.min(btlen-2,p));
				backtr = Interpolation.powint(vrg[p],vrg[p+1],vr[p],vr[p+1],vrgs,1d);
			}
			res[i] = backtr;
		}
		df.addColumn(name_transformed, res);
	}

	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 2d-dataframe
	 * @param variable_id            id of the variable in the 2d-dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 2d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame2D df, int variable_id, int weights_id, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 2d-dataframe
	 * @param variable               name of the variable in the 2d-dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 2d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame2D df, String variable, String weight_variable, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, variable, new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 2d-dataframe
	 * @param variable_id            id of the variable in the 2d-dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 2d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame2D df, int variable_id, int weights_id, double[] trims,
			double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), trims, nscore_trnsf_table, name_transformed,
				    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 2d-dataframe
	 * @param variable               name of the variable in the 2d-dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 2d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame2D df, String variable, double[] trims, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return;
		}
		if(variable==null) {
			System.err.println("Variable for normal score transform does not exist!");
			return;
		}
		if(nscore_trnsf_table==null) {
			System.err.println("Normal score transformation table does not exist!\n"+
					           "    no back transformation done!");
			return;
		}
		if(lower_tail_intpol_type!=1 && lower_tail_intpol_type!=2) {
			System.err.println("For lower tail only interpolation types 1 and 2 possible!\n"+
					           "    Found "+lower_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		if(upper_tail_intpol_type!=1 && upper_tail_intpol_type!=2 && upper_tail_intpol_type!=4) {
			System.err.println("For upper tail only interpolation types 1,2 and 4 possible!\n"+
					           "    Found "+upper_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		
		
		double[][] val = (double[][]) df.getArray(variable);
		int[] nd = {val.length, val[0].length};
		double[][] res = new double[nd[0]][nd[1]];
		double lambda;
		
		double[] vrg = nscore_trnsf_table[1];
		double[] vr  = nscore_trnsf_table[0];
		int btlen = vrg.length;
		
		if(upper_tail_intpol_type==4 && vr[btlen-1]<0d) {
			System.err.println("ERROR can not use hyperbolic tail (intpol type 4) with negative values! - see manual\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[0]>vr[0]) {
			System.err.println("ERROR zmin should be no larger than the first entry in the transformation table\n" + 
					           "    zmin = "+trims[0]+" vr_first="+vr[0]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[1]<vr[btlen-1]) {
			System.err.println("ERROR zmax should be no less than the last entry in the transformation table\n" + 
					           "    zmax = "+trims[1]+" vr_last="+vr[btlen-1]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		
		for(int j=0; j<nd[0]; j++) for(int i=0; i<nd[1]; i++) {
//	        i added by Enrico Metzner
//	        i check whether data is fill_value and continue or it can be transformed
			if(val[j][i]==Constants.FILL_VALUE_D) {
				res[j][i] = Constants.FILL_VALUE_D;
				continue;
			}
//        c
//        c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
//        c
			double vrgs = val[j][i];
			double backtr;
			if(vrgs < vrg[0]) {
				backtr = vr[0];
				double cdflo  = Gaussian.cdf(vrg[0]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(lower_tail_intpol_type==1) {
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,1.0);
				} else if(lower_tail_intpol_type==2) {
					double cpow   = 1d / lower_tail_power;
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,cpow);
				}
//        c
//        c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
//        c
			} else if(vrgs >= vrg[btlen-1]) {
				backtr = vr[btlen-1];
				double cdfhi  = Gaussian.cdf(vrg[btlen-1]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(upper_tail_intpol_type==1) {
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,1d);
				} else if(upper_tail_intpol_type==2) {
					double cpow   = 1d / upper_tail_power;
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,cpow);
				} else if(upper_tail_intpol_type==4) {
					lambda = Math.pow(vr[btlen-1],upper_tail_power)*(1.0-Gaussian.cdf(vrg[btlen-1]));
					backtr = Math.pow(lambda/(1.0-Gaussian.cdf(vrgs)),1d/upper_tail_power);
				}
//        c
//        c Value within the transformation table:
//        c
			} else {
				int p;
				for(p=0; p<btlen-1; p++) if(vrgs<vrg[p+1]) break;
				p = Math.max(0, Math.min(btlen-2,p));
				backtr = Interpolation.powint(vrg[p],vrg[p+1],vr[p],vr[p+1],vrgs,1d);
			}
			res[j][i] = backtr;
		}
		df.addColumn(name_transformed, res);
	}

	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 3d-dataframe
	 * @param variable_id            id of the variable in the 3d-dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 3d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame3D df, int variable_id, int weights_id, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 3d-dataframe
	 * @param variable               name of the variable in the 3d-dataframe
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 3d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame3D df, String variable, String weight_variable, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, variable, new double[] {-1.0e21d, 1.0e21d}, nscore_trnsf_table, name_transformed,
			    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 3d-dataframe
	 * @param variable_id            id of the variable in the 3d-dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 3d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame3D df, int variable_id, int weights_id, double[] trims,
			double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		back_nscore(df, df.getVarname(variable_id), trims, nscore_trnsf_table, name_transformed,
				    lower_tail_intpol_type, upper_tail_intpol_type, lower_tail_power, upper_tail_power);
	}
	/**
	 * Back Transform Univariate Data from Normal Scores
	 * <p>
	 *     This subroutine backtransforms a standard normal deviate from a
	 *     specified back transform table and option for the tails of the
	 *     distribution.  Call once with "first" set to true then set to false
	 *     unless one of the options for the tail changes.
	 * </p>
	 * 
	 * @param df                     data containing 3d-dataframe
	 * @param variable               name of the variable in the 3d-dataframe
	 * @param trims                  double array with lower and upper bound of valid range
	 * @param nscore_trnsf_table     transformation table as Array with var-transformed-pairs
	 * @param name_transformed       name of the transformed result to put it in the 3d-dataframe
	 * @param lower_tail_intpol_type option to handle values less than minimum of transformation table
	 * @param upper_tail_intpol_type option to handle values greater than maximum of transformation table
	 * @param lower_tail_power       parameter required for option lower_tail_intpol_type
	 * @param upper_tail_power       parameter required for option upper_tail_intpol_type
	 */
	public static void back_nscore(DataFrame3D df, String variable, double[] trims, double[][] nscore_trnsf_table, String name_transformed,
			int lower_tail_intpol_type, int upper_tail_intpol_type, double lower_tail_power, double upper_tail_power) {
		if(df==null) {
			System.err.println("Dataframe with input variables does not exist!");
			return;
		}
		if(variable==null) {
			System.err.println("Variable for normal score transform does not exist!");
			return;
		}
		if(nscore_trnsf_table==null) {
			System.err.println("Normal score transformation table does not exist!\n"+
					           "    no back transformation done!");
			return;
		}
		if(lower_tail_intpol_type!=1 && lower_tail_intpol_type!=2) {
			System.err.println("For lower tail only interpolation types 1 and 2 possible!\n"+
					           "    Found "+lower_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		if(upper_tail_intpol_type!=1 && upper_tail_intpol_type!=2 && upper_tail_intpol_type!=4) {
			System.err.println("For upper tail only interpolation types 1,2 and 4 possible!\n"+
					           "    Found "+upper_tail_intpol_type+", so no back transformation performed!");
			return;
		}
		
		
		double[][][] val = (double[][][]) df.getArray(variable);
		int[] nd = {val.length, val[0].length, val[0][0].length};
		double[][][] res = new double[nd[0]][nd[1]][nd[2]];
		double lambda;
		
		double[] vrg = nscore_trnsf_table[1];
		double[] vr  = nscore_trnsf_table[0];
		int btlen = vrg.length;
		
		if(upper_tail_intpol_type==4 && vr[btlen-1]<0d) {
			System.err.println("ERROR can not use hyperbolic tail (intpol type 4) with negative values! - see manual\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[0]>vr[0]) {
			System.err.println("ERROR zmin should be no larger than the first entry in the transformation table\n" + 
					           "    zmin = "+trims[0]+" vr_first="+vr[0]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		if(trims[1]<vr[btlen-1]) {
			System.err.println("ERROR zmax should be no less than the last entry in the transformation table\n" + 
					           "    zmax = "+trims[1]+" vr_last="+vr[btlen-1]+"\n"+
					           "    No back transformation performed!");
			return;
		}
		
		for(int k=0; k<nd[0]; k++) for(int j=0; j<nd[1]; j++) for(int i=0; i<nd[2]; i++) {
//        i added by Enrico Metzner
//        i check whether data is fill_value and continue or it can be transformed
			if(val[k][j][i]==Constants.FILL_VALUE_D) {
				res[k][j][i] = Constants.FILL_VALUE_D;
				continue;
			}
//        c
//        c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
//        c
			double vrgs = val[k][j][i];
			double backtr;
			if(vrgs < vrg[0]) {
				backtr = vr[0];
				double cdflo  = Gaussian.cdf(vrg[0]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(lower_tail_intpol_type==1) {
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,1.0);
				} else if(lower_tail_intpol_type==2) {
					double cpow   = 1d / lower_tail_power;
					backtr = Interpolation.powint(0.0,cdflo,trims[0],vr[0],cdfbt,cpow);
				}
//        c
//        c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
//        c
			} else if(vrgs >= vrg[btlen-1]) {
				backtr = vr[btlen-1];
				double cdfhi  = Gaussian.cdf(vrg[btlen-1]);
				double cdfbt  = Gaussian.cdf(vrgs);
				if(upper_tail_intpol_type==1) {
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,1d);
				} else if(upper_tail_intpol_type==2) {
					double cpow   = 1d / upper_tail_power;
					backtr = Interpolation.powint(cdfhi,1.0,vr[btlen-1],trims[1],cdfbt,cpow);
				} else if(upper_tail_intpol_type==4) {
					lambda = Math.pow(vr[btlen-1],upper_tail_power)*(1.0-Gaussian.cdf(vrg[btlen-1]));
					backtr = Math.pow(lambda/(1.0-Gaussian.cdf(vrgs)),1d/upper_tail_power);
				}
//        c
//        c Value within the transformation table:
//        c
			} else {
				int p;
				for(p=0; p<btlen-1; p++) if(vrgs<vrg[p+1]) break;
				p = Math.max(0, Math.min(btlen-2,p));
				backtr = Interpolation.powint(vrg[p],vrg[p+1],vr[p],vr[p+1],vrgs,1d);
			}
			res[k][j][i] = backtr;
		}
		df.addColumn(name_transformed, res);
	}
}
