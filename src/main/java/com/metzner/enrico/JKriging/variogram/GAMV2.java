package com.metzner.enrico.JKriging.variogram;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.helper.DataHelper;
import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.MathHelper;

public class GAMV2 {
	public final static int    SEMI_VARIOGRAM   =  1;
	public final static int    CROSS_VARIOGRAM  =  2;
	public final static int    COVARIANCE       =  3;
	public final static int    CORRELOGRAM      =  4;
	public final static int    GENERAL_RELATIVE =  5;
	public final static int    PAIR_RELATIVE    =  6;
	public final static int    VARIO_OF_LOG     =  7;
	public final static int    SEMI_MADOGRAM    =  8;
	public final static int    CONT_INDICATOR   =  9;
	public final static int    CAT_INDICATOR    = 10;
	public final static double VERSION  = 3.000d;
	
	private final static double D2R = Math.PI/180d;
	//private final static double R2D = 180d/Math.PI;

	private List<double[]> direction;
	private List<double[]> variogram;
	String x_var, y_var, z_var;
	private int[]    ivc, indflag, dirnum, lagnum;
	private double[] dis, gam, hm, tm, hv, tv, np;
	private String[] names;

	private double xlag, xltol, tmin, tmax;
	private int    nlag, isill, ncut;
	
	private double progressNVario, progressNDir, progressNPairs;
	
	private DataFrame dataframe;

	private boolean stop;
	private boolean[] paramchecklist    = {false, false, false, false, false};
	private String[] param_descriptions = {"Lags", "Coordinates", "Variables", "Directions", "VariogramTypes"};



	public GAMV2() {
		
		direction = new ArrayList<double[]>();
		variogram = new ArrayList<double[]>();
		x_var = null; y_var = null; z_var = null;
		dataframe = null;
		
		reset();
	}


	public void setLagDefinition(int number_of_lags, double _lag_spacing, double _lag_tolerance) {
		nlag = number_of_lags;
		System.out.println(" number of lags = "+nlag);
		if(nlag<1) { System.err.println("nlag is too small: check parameters"); return; }
		xlag = _lag_spacing;
		System.out.println(" lag distance = "+xlag);
		if(xlag<0d) { System.err.println("xlag is too small: check parameter file"); return; }
		xltol = _lag_tolerance; if(xltol<0d) xltol = 0.5d * xlag;
		System.out.println(" lag tolerance = "+xltol);
		paramchecklist[0] = true;
	}
	public void setCoordinates(int _x_id, int _y_id, int _z_id) {
		if(dataframe == null) { System.err.println("No dataframe set!"); return; }
		int nvar = dataframe.getVariableCount();
		if(_x_id>nvar) { System.err.println("x variable number out of range!"); return; }
		if(_y_id>nvar) { System.err.println("y variable number out of range!"); return; }
		if(_z_id>nvar) { System.err.println("z variable number out of range!"); return; }
		boolean atLeastOneDim = (_x_id>0 || _y_id>0 || _z_id>0);
		if(!atLeastOneDim) { System.err.println("At least one dimension has to be set!"); return; }
		//System.out.println(" tail,head,type = "+dataframe.getVarname(tailID-1)+", "+dataframe.getVarname(headID-1)+", "+variogramTypeById(type));
		//variogram.add(new double[] {tailID*1.0001d,headID*1.0001d,type*1.0001d,_cut});
		//paramchecklist[4] = true;
		x_var = dataframe.getVarname(_x_id-1+Constants.FIRST_IDX);
		y_var = dataframe.getVarname(_y_id-1+Constants.FIRST_IDX);
		z_var = dataframe.getVarname(_z_id-1+Constants.FIRST_IDX);
		System.out.println(" columns for X,Y,Z = "+_x_id+" "+_y_id+" "+_z_id+" = "+
				(x_var!=null?"\""+x_var+"\"":"null")+" "+(y_var!=null?"\""+y_var+"\"":"null")+" "+(z_var!=null?"\""+z_var+"\"":"null"));
		paramchecklist[1] = true;
	}
	public void setCoordinates(String xVarName, String yVarName, String zVarName) {
		if(dataframe == null) { System.err.println("No dataframe set!"); return; }
		setCoordinates(dataframe.getVariableID(xVarName)-Constants.FIRST_IDX+1,
				dataframe.getVariableID(yVarName)-Constants.FIRST_IDX+1,
				dataframe.getVariableID(zVarName)-Constants.FIRST_IDX+1);
	}
	public void setDataframe(DataFrame _df) {
		dataframe = _df;
		paramchecklist[2] = true;
	}
	public void addMathDirection(double azm, double atol, double bandwh, double dip, double dtol, double bandwd) {
		addMathDirection(azm, atol, bandwh, dip, dtol, bandwd, 'd');
	}
	public void addMathDirection(double azm, double atol, double bandwh, double dip, double dtol, double bandwd, char angleType) {
		double a = azm;
		double at = atol;
		double d = dip;
		double dt = dtol;
		switch(angleType) {
			case 'd': //DEGREE->DEGREE
				break;
			case 'g': //GON->DEGREE
				a *= 0.9d; at *= 0.9d;
				d *= 0.9d; dt *= 0.9d;
				break;
			case 'r': //RADIANS->DEGREE
				a *= Constants.RAD2DEG; at *= Constants.RAD2DEG;
				d *= Constants.RAD2DEG; dt *= Constants.RAD2DEG;
				break;
			default: break;
		}
		direction.add(new double[] {a,at,bandwh, d,dt,bandwd});
		paramchecklist[3] = true;
	}
	public void addGeoDirection(double azm, double atol, double bandwh, double dip, double dtol, double bandwd) {
		addGeoDirection(azm, atol, bandwh, dip, dtol, bandwd, 'd');
	}
	public void addGeoDirection(double azm, double atol, double bandwh, double dip, double dtol, double bandwd, char angleType) {
		double a = azm;
		double at = atol;
		double d = dip;
		double dt = dtol;
		switch(angleType) {
			case 'd': //DEGREE->DEGREE
				break;
			case 'g': //GON->DEGREE
				a *= 0.9d; at *= 0.9d;
				d *= 0.9d; dt *= 0.9d;
				break;
			case 'r': //RADIANS->DEGREE
				a *= Constants.RAD2DEG; at *= Constants.RAD2DEG;
				d *= Constants.RAD2DEG; dt *= Constants.RAD2DEG;
				break;
			default: break;
		}
		direction.add(new double[] {90d-a,at,bandwh, d,dt,bandwd});
		paramchecklist[3] = true;
	}
	public void addVariogramDef(int type, int tailID, int headID, double _cut) {
		if(dataframe == null) { System.err.println("No dataframe set!"); return; }
		int nvar = dataframe.getVariableCount();
		if(tailID<1 || tailID>nvar) { System.err.println("tail variable number out of range!"); return; }
		if(headID<1 || headID>nvar) { System.err.println("head variable number out of range!"); return; }
		System.out.println(" tail,head,type = "+dataframe.getVarname(tailID-1+Constants.FIRST_IDX)+", "+
				dataframe.getVarname(headID-1+Constants.FIRST_IDX)+", "+variogramTypeById(type));
		double cut = 0d;
		if(type==CONT_INDICATOR || type==CAT_INDICATOR) {
			ncut++;
			if(tmin<0d) throw new RuntimeException("tmin interferes with indicators!");
			if(tmax<=1d) throw new RuntimeException("tmax interferes with indicators!");
//			int ii = Integer.parseInt(parts[0].trim());
//			int jj = Integer.parseInt(parts[1].trim());
//			int kk = Integer.parseInt(parts[2].trim());
			cut = _cut;
			if(type==CONT_INDICATOR) indflag = addIntToArray(indflag, 1);
			if(type==CAT_INDICATOR) indflag = addIntToArray(indflag, 0);
			ivc = addIntToArray(ivc, tailID);
			tailID = nvar + ncut;
			headID = nvar + ncut;
			names[nvar+ncut-1] = "Indicator "+ncut;
			System.out.println(" indicator threshold = "+cut);
		}
		variogram.add(new double[] {tailID*1.0001d,headID*1.0001d,type*1.0001d,cut});
		paramchecklist[4] = true;
	}
	public void addVariogramDef(int type, String tailVarName, String headVarName, double _cut) {
		addVariogramDef(type, dataframe.getVariableID(tailVarName)-Constants.FIRST_IDX+1,
				dataframe.getVariableID(headVarName)-Constants.FIRST_IDX+1, _cut);
	}
	public void standardizeSills(boolean _i_sill) {
		isill = (_i_sill ? 1 : 0);
	}
	public void setTminTmax(double _tmin, double _tmax) { tmin = _tmin; tmax = _tmax; }
	public void reset() {
		dataframe = null;
		x_var = null;
		y_var = null;
		z_var = null;
		
		stop = false;
		paramchecklist[0] = false; // lag decription
		paramchecklist[1] = false; // coordinates
		paramchecklist[2] = false; // variable definition
		paramchecklist[3] = false; // variogram directions
		paramchecklist[4] = false; // variogram definitions
		
		tmin = Double.NEGATIVE_INFINITY;
		tmax = Double.POSITIVE_INFINITY;
		
		isill = 0;
		ncut = 0;
		ivc = new int[0];
		indflag = new int[0];
		
		progressNVario = 0d;
		progressNDir = 0d;
		progressNPairs = 0d;
		
		direction.clear();
		variogram.clear();
	}
	public void readParameterfile(String param_path, DataFrame df) {
		direction.clear();
		variogram.clear();
		x_var = null; y_var = null; z_var = null;
		
		dataframe = df;
//        subroutine readparm
//        c-----------------------------------------------------------------------
//        c
//        c                  Initialization and Read Parameters
//        c                  **********************************
//        c
//        c The input parameters and data are read in from their files. Some quick
//        c error checking is performed and the statistics of all the variables
//        c being considered are written to standard output.
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		int MV=500; //parameter
		int[] ivar=new int[MV], ivc=new int[MV], indflag=new int[MV];

		int ncut = 0;
//        c Note VERSION number:
		System.out.println(" GAM Version: "+FormatHelper.nf(VERSION,5,3));

//        c Get the name of the parameter file - try the default name if no input:
		if(param_path==null) {
			param_path = "res/gamv.par";
		} else if(param_path.length()==0) {
			param_path = "res/gamv.par";
		}
		File f = new File(param_path);
		if(!f.exists()) {
			System.err.println("ERROR - the parameter file does not exist,\n"+
					           "        check for the file and try again\n");
			if(param_path.equals("res/gamv.par")) {
				System.out.println("        creating a blank parameter file in res/gamv.par");
				makepar();
			}
			return;
		}
		try(BufferedReader br = new BufferedReader(new FileReader(f))) {
			String line = br.readLine();

//        c Find Start of Parameters:
			while(!line.startsWith("STAR")) line = br.readLine();

//        c Read Input Parameters:
			line = br.readLine(); String[] parts = FormatHelper.splitBySpace(line);
			//datafl = parts[0]; FormatHelper.chknam(datafl); System.out.println(" data file = "+datafl);
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			int ixl=Integer.parseInt(parts[0].trim());
			if(ixl>0) x_var = dataframe.getVarname(ixl-1+Constants.FIRST_IDX);
			int iyl=Integer.parseInt(parts[1].trim());
			if(iyl>0) y_var = dataframe.getVarname(iyl-1+Constants.FIRST_IDX);
			int izl=Integer.parseInt(parts[2].trim());
			if(izl>0) z_var = dataframe.getVarname(izl-1+Constants.FIRST_IDX);
			System.out.println(" columns for X,Y,Z = \""+ixl+"\", \""+iyl+"\", \""+izl+"\"");
			paramchecklist[1] = true;
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			int nvar=Integer.parseInt(parts[0].trim());
			System.out.println(" number of variables = "+nvar);
			
			String o = ""; for(int i=0; i<nvar; i++) { ivar[i] = Integer.parseInt(parts[i+1].trim()); o+=" "+ivar[i]; }
			System.out.println(" columns ="+o);
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			tmin = Double.parseDouble(parts[0].trim()); tmax = Double.parseDouble(parts[1].trim());
			System.out.println(" trimming limits = "+tmin+" "+tmax);
			paramchecklist[2] = true;
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			//outfl = parts[0]; FormatHelper.chknam(outfl); System.out.println(" output file = "+outfl);
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line); nlag = Integer.parseInt(parts[0].trim());
			System.out.println(" number of lags = "+nlag); //TODO
			if(nlag<1) { System.err.println("nlag is too small: check parameters"); return; }
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line); xlag = Double.parseDouble(parts[0].trim());
			System.out.println(" lag distance = "+xlag);
			if(xlag<0d) { System.err.println("xlag is too small: check parameter file"); return; }
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line); xltol = Double.parseDouble(parts[0].trim());
			System.out.println(" lag tolerance = "+xltol);
			paramchecklist[0] = true;
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			int ndir = Integer.parseInt(parts[0].trim());
			System.out.println(" number of directions = "+ndir);
			if(ndir<1) { System.err.println("ndir is too small: check parameters"); return; }
			
			for(int i=0; i<ndir; i++) {
				line = br.readLine(); parts = FormatHelper.splitBySpace(line);
				double azm    = Double.parseDouble(parts[0].trim());
				double atol   = Double.parseDouble(parts[1].trim());
				double bandwh = Double.parseDouble(parts[2].trim());
				double dip    = Double.parseDouble(parts[3].trim());
				double dtol   = Double.parseDouble(parts[4].trim());
				double bandwd = Double.parseDouble(parts[5].trim());
				System.out.println(" azm, atol, bandwh = "+azm+" "+atol+" "+bandwh);
				System.out.println(" dip, dtol, bandwd = "+dip+" "+dtol+" "+bandwd);
				if(bandwh<0d) throw new RuntimeException(" Horizontal bandwidth is too small!");
				if(bandwd<0d) throw new RuntimeException(" Vertical bandwidth is too small!");
				direction.add(new double[] {azm,atol,bandwh, dip,dtol,bandwd});
			}
			paramchecklist[3] = true;
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line); isill = Integer.parseInt(parts[0].trim());
			System.out.println(" flag to standardize sills = "+isill);
			
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			int nvarg = Integer.parseInt(parts[0].trim());
			System.out.println(" number of variograms = "+nvarg);
			if(nvarg<1) throw new RuntimeException("nvarg is too small: check parameters");
			
			names = new String[nvar+nvarg];
			ncut = 0;
			for(int i=0; i<nvarg; i++) {
				line = br.readLine(); parts = FormatHelper.splitBySpace(line);
				int ivtail = Integer.parseInt(parts[0].trim());
				if(ivtail<1 || ivtail>nvar) { System.err.println("tail variable number out of range!"); return; }
				ivtail = ivar[ivtail-1]-1;
				int ivhead = Integer.parseInt(parts[1].trim());
				if(ivhead<1 || ivhead>nvar) { System.err.println("head variable number out of range!"); return; }
				ivhead = ivar[ivhead-1]-1;
				int ivtype = Integer.parseInt(parts[2].trim());
				System.out.println(" tail,head,type = "+dataframe.getVarname(ivtail-1+Constants.FIRST_IDX)+", "+
						dataframe.getVarname(ivhead-1+Constants.FIRST_IDX)+", "+variogramTypeById(ivtype));
				double cut = 0d;
				if(ivtype==9 || ivtype==10) {
					ncut++;
					if(tmin<0d) throw new RuntimeException("tmin interferes with indicators!");
					if(tmax<=1d) throw new RuntimeException("tmax interferes with indicators!");
//					int ii = Integer.parseInt(parts[0].trim());
//					int jj = Integer.parseInt(parts[1].trim());
//					int kk = Integer.parseInt(parts[2].trim());
					cut = Double.parseDouble(parts[3].trim());
					if(ivtype==9) indflag[ncut-1] = 1;
					if(ivtype==10) indflag[ncut-1] = 0;
					ivc[ncut-1] = ivtail;
					ivtail = nvar + ncut;
					ivhead = nvar + ncut;
					names[nvar+ncut-1] = "Indicator "+ncut;
					System.out.println(" indicator threshold = "+cut);
				}
				variogram.add(new double[] {ivtail*1.0001d,ivhead*1.0001d,ivtype*1.0001d,cut});
			}
			paramchecklist[4] = true;
		} catch(IOException | NullPointerException io_np_e) {
			System.err.println("ERROR in parameter file");
			io_np_e.printStackTrace();
			return;
		}

//        c Perform some quick error checking:
		if(xltol <= 0d) {
			System.out.println("WARNING: xltol is too small: resetting to xlag/2");
			xltol = 0.5d*xlag;
		}
    }

	public boolean hasStopped() {
		return !stop;
	}
	public double getAllProgress() {
		return Math.min(1d, progressNVario + 
				progressNDir / variogram.size() + 
				progressNPairs / (variogram.size() * direction.size()));
	}
	public double[] getProgress() {
		return new double[] {progressNVario, progressNDir, progressNPairs};
	}
	public int numberOfVariograms() {
		return variogram.size();
	}
	public int numberOfDirections() {
		return direction.size();
	}
	
	public DataFrame gamv() {
		for(int pcl=0; pcl<5; pcl++)
			if(!paramchecklist[pcl]) {
				System.err.println("Not all parameters are set! \""+param_descriptions[pcl]+"\" missing");
				return null;
			}
		double[] direct = direction.get(0);
		double[] variog = variogram.get(0);
		System.out.println("[DEBUG]\n"+
			" X,Y,Z => "+x_var+" "+y_var+" "+z_var+"\n"+
			" dir 1:   "+direct[0]+" "+direct[1]+" "+direct[2]+" "+direct[3]+" "+direct[4]+" "+direct[5]+"\n"+
			" vario 1: "+variogramTypeById((int)variog[2])+" "+dataframe.getVarname((int)variog[0]-1+Constants.FIRST_IDX)+" "+
				dataframe.getVarname((int)variog[1]-1+Constants.FIRST_IDX));
//        c-----------------------------------------------------------------------
//        c
//        c              Variogram of 3-D Irregularly Spaced Data
//        c              ****************************************
//        c
//        c This subroutine computes a variety of spatial continuity measures of a
//        c set for irregularly spaced data.  The user can specify any combination
//        c of direct and cross variograms using any of eight "variogram" measures
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   nd               Number of data (no missing values)
//        c   x(nd)            X coordinates of the data
//        c   y(nd)            Y coordinates of the data
//        c   z(nd)            Z coordinates of the data
//        c   nv               The number of variables
//        c   vr(nd,nv)        Data values
//        c   tmin,tmax        Trimming limits
//        c   nlag             Number of lags to calculate
//        c   xlag             Length of the unit lag
//        c   xltol            Distance tolerance (if <0 then set to xlag/2)
//        c   ndir             Number of directions to consider
//        c   azm(ndir)        Azimuth angle of direction (measured positive
//        c                      degrees clockwise from NS).
//        c   atol(ndir)       Azimuth (half window) tolerances
//        c   bandwh           Maximum Horizontal bandwidth (i.e., the deviation 
//        c                      perpendicular to the defined azimuth).
//        c   dip(ndir)        Dip angle of direction (measured in negative
//        c                      degrees down from horizontal).
//        c   dtol(ndir)       Dip (half window) tolerances
//        c   bandwd           Maximum "Vertical" bandwidth (i.e., the deviation
//        c                      perpendicular to the defined dip).
//        c   isill            1=attempt to standardize, 0=do not
//        c   sills            the sills (variances) to standardize with
//        c   nvarg            Number of variograms to compute
//        c   ivtail(nvarg)    Variable for the tail of each variogram
//        c   ivhead(nvarg)    Variable for the head of each variogram
//        c   ivtype(nvarg)    Type of variogram to compute:
//        c                      1. semivariogram
//        c                      2. cross-semivariogram
//        c                      3. covariance
//        c                      4. correlogram
//        c                      5. general relative semivariogram
//        c                      6. pairwise relative semivariogram
//        c                      7. semivariogram of logarithms
//        c                      8. rodogram
//        c                      9. indicator semivariogram (continuous)
//        c                     10. indicator semivariogram (categorical)
//        c
//        c
//        c
//        c OUTPUT VARIABLES:  The following arrays are ordered by direction,
//        c                    then variogram, and finally lag, i.e.,
//        c                      iloc = (id-1)*nvarg*MAXLG+(iv-1)*MAXLG+il
//        c
//        c   np()             Number of pairs
//        c   dis()            Distance of pairs falling into this lag
//        c   gam()            Semivariogram, covariance, correlogram,... value
//        c   hm()             Mean of the tail data
//        c   tm()             Mean of the head data
//        c   hv()             Variance of the tail data
//        c   tv()             Variance of the head data
//        c
//        c
//        c
//        c PROGRAM NOTES:
//        c
//        c   1. The file "gamv.inc" contains the dimensioning parameters.
//        c      These may have to be changed depending on the amount of data
//        c      and the requirements to compute different variograms.
//        c
//        c
//        c
//        c Original:  A.G. Journel                                           1978
//        c Revisions: K. Guertin                                             1980
//        c-----------------------------------------------------------------------
		//    use geostat
		double PI = Math.PI;
		
		//precheck and warnings:
		check(dataframe.getMins(),dataframe.getMaxs());
		
		//allocate variables:
		int ndir  = direction.size();
		int nvarg = variogram.size();
		int mxdlv = ndir*(nlag+2)*nvarg;
		
		dirnum = new int[mxdlv];
		lagnum = new int[mxdlv];
		dis = new double[mxdlv];
		gam = new double[mxdlv];
		hm  = new double[mxdlv];
		tm  = new double[mxdlv];
		hv  = new double[mxdlv];
		tv  = new double[mxdlv];
		np  = new double[mxdlv];
//		double[] uvxazm = new double[ndir];
//		double[] uvyazm = new double[ndir];
//		double[] csatol = new double[ndir];
//		double[] uvzdec = new double[ndir];
//		double[] uvhdec = new double[ndir];
//		double[] csdtol = new double[ndir];

		int datalength = 0;
		double[] x=null,y=null,z=null;
		if(dataframe.hasVariable(x_var)) {x = ((double[]) dataframe.getArray(x_var)).clone(); datalength = Math.max(datalength,x.length); }
		if(dataframe.hasVariable(y_var)) {y = ((double[]) dataframe.getArray(y_var)).clone(); datalength = Math.max(datalength,y.length); }
		if(dataframe.hasVariable(z_var)) {z = ((double[]) dataframe.getArray(z_var)).clone(); datalength = Math.max(datalength,z.length); }
		if(x_var==null) {x = new double[datalength]; for(int d=0; d<datalength; d++) x[d] = 0d; }
		if(y_var==null) {y = new double[datalength]; for(int d=0; d<datalength; d++) y[d] = 0d; }
		if(z_var==null) {z = new double[datalength]; for(int d=0; d<datalength; d++) z[d] = 0d; }

//      c Initialize the arrays for each direction, variogram, and lag:
		for(int i=0; i<mxdlv; i++) {
			np[i]  = 0d;
			dis[i] = 0d; gam[i] = 0d;
			hm[i]  = 0d; tm[i]  = 0d;
			hv[i]  = 0d; tv[i]  = 0d;
		}

//        c Define the distance tolerance if it isn't already:
		if(xltol <= 0d) xltol = 0.5d * xlag;

//        c Define the angles and tolerance for each direction:
		double dismxs = ((nlag + 0.5d - Constants.D_EPSLON) * xlag); dismxs *= dismxs;
		
//        c MAIN LOOP OVER ALL VARIOGRAMS:
		for(int iv=0; iv<nvarg; iv++) {
			double[] vario_param = variogram.get(iv);
			int ivtail = (int) vario_param[0];
			int ivhead = (int) vario_param[1];
			int ivtype = (int) vario_param[2];
			//double ivcut = vario_param[3];
			int varcount = dataframe.getVariableCount();
			if(ivhead<1 || ivhead>varcount) {
				System.err.println("[GAMV] Cannot find a variable with ID "+ivhead+" for \"head\"");
				continue;
			}
			if(ivtail<1 || ivtail>varcount) {
				System.err.println("[GAMV] Cannot find a variable with ID "+ivtail+" for \"tail\"");
				continue;
			}
			double[] tail = ((double[]) dataframe.getArray(ivtail-1+Constants.FIRST_IDX)).clone();
			double[] head = ((double[]) dataframe.getArray(ivhead-1+Constants.FIRST_IDX)).clone();
			DataHelper.sortem(head, tail);

//        c Definition of the direction corresponding to the current pair. All
//        c directions are considered (overlapping of direction tolerance cones
//        c is allowed):
			for(int id=0; id<ndir; id++) { // Loop ID: 5
				double[] dir = direction.get(id);
				double azm    = dir[0];
				double atol   = dir[1];
				double bandwh = dir[2];
				double dip    = dir[3];
				double dtol   = dir[4];
				double bandwd = dir[5];

//        c The mathematical azimuth is measured counterclockwise from EW and
//        c not clockwise from NS as the conventional azimuth is:
				double azmuth = azm * D2R; //reworked to use mathematical angles everywhere as possible
				double uvxazm = Math.cos(azmuth);
				double uvyazm = Math.sin(azmuth);
				double csatol = (atol<=0d ? Math.cos(45d*PI/180d) : Math.cos(atol*PI/180d));

//        c The declination is measured positive down from vertical (up) rather
//        c than negative down from horizontal:
				double declin = (90d-dip) * D2R;
				double uvzdec = Math.cos(declin);
				double uvhdec = Math.sin(declin);
				double csdtol = (dtol<=0d ? Math.cos(45d*D2R) : Math.cos(dtol*D2R));

//        c MAIN LOOP OVER ALL PAIRS:
				int irepo = Math.max(1, Math.min(datalength/10, 1000));
				for(int i=0; i<datalength; i++) { // Loop ID: 3
					if( (i+1)%irepo == 0 ) System.out.println("   currently on "+(i+1)+" of "+datalength);
					for(int j=i; j<datalength; j++) { // Loop ID: 4

//        c Definition of the lag corresponding to the current pair:
						double dx  = x[j] - x[i];
						double dy  = y[j] - y[i];
						double dz  = z[j] - z[i];
						double dxs = dx*dx;
						double dys = dy*dy;
						double dzs = dz*dz;
						double hs  = dxs + dys + dzs;
						if(hs > dismxs) continue;
						double h   = Math.sqrt(hs);

//        c Determine which lag this is and skip if outside the defined distance
//        c tolerance:
						int lagbeg,lagend;
						if(h <= Constants.D_EPSLON) {
							lagbeg = 1;
							lagend = 1;
						} else {
							lagbeg = -1;
							lagend = -1;
							for(int ilag=2; ilag<=nlag+2; ilag++) {
								if(h >= (xlag*(ilag-2d)-xltol) && h <= (xlag*(ilag-2d)+xltol)) {
									if(lagbeg<0) lagbeg = ilag;
									lagend = ilag;
								}
							}
							if(lagend<0) continue;
						}

//        c Check for an acceptable azimuth angle:
						double dxy = Math.sqrt(Math.max(dxs+dys,0d));
						double dcazm = (dxy<Constants.D_EPSLON ? 1d : (dx*uvxazm+dy*uvyazm)/dxy);
						if(Math.abs(dcazm) < csatol) continue; // goto 5

//        c Check the horizontal bandwidth criteria (maximum deviation 
//        c perpendicular to the specified direction azimuth):
						double band = uvxazm*dy - uvyazm*dx;
						if(Math.abs(band) > bandwh) continue; // goto 5

//        c Check for an acceptable dip angle:
						double dcdec;
						if(dcazm < 0d) dxy = -dxy;
						if(lagbeg==1) {
							dcdec = 0d;
						} else {
							dcdec = (dxy*uvhdec+dz*uvzdec)/h;
							if(Math.abs(dcdec) < csdtol) continue; // goto 5
						}

//        c Check the vertical bandwidth criteria (maximum deviation perpendicular
//        c to the specified dip direction):
						band = uvhdec*dz - uvzdec*dxy;
						if(Math.abs(band) > bandwd) continue; // goto 5

//        c Check whether or not an omni-directional variogram is being computed:
						boolean omni = false;
						if(atol >= 90d) omni = true;

//        c For this variogram, sort out which is the tail and the head value:
						double[] _vrth = sub_tail_head(i,j,tail,head, dcazm>=0d && dcdec>=0d,omni, ivtype);
						double[] pair = {_vrth[0], _vrth[1]};
						double[] pairpr = {_vrth[2], _vrth[3]};

//        c Reject this pair on the basis of missing values:
						if(Double.isNaN(pair[0]) || Double.isNaN(pair[1])) continue; // goto 6
						if(pair[0]<tmin || pair[1]<tmin || pair[0]>tmax || pair[1]>tmax) continue; // goto 6
						if(ivtype==CROSS_VARIOGRAM && (Double.isNaN(pairpr[0]) || Double.isNaN(pairpr[1]))) continue; // goto 6
						if(ivtype==CROSS_VARIOGRAM && (pairpr[0]<tmin || pairpr[1]<tmin || pairpr[0]>tmax || pairpr[1]>tmax)) continue; // goto 6

//        c             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
//        c
//        c
//        c The Semivariogram:
//        c
						int index_off = id*nvarg*(nlag+2)+iv*(nlag+2);
						if(ivtype==SEMI_VARIOGRAM || ivtype==GENERAL_RELATIVE || ivtype==CONT_INDICATOR || ivtype==CAT_INDICATOR) {
							sub_semivariogram(index_off, h, pair,pairpr, lagbeg,lagend, omni);
						} else if(ivtype==CROSS_VARIOGRAM) {
							sub_cross_semivariogram(index_off, h, pair,pairpr, lagbeg,lagend);
						} else if(Math.abs(ivtype)==COVARIANCE) {
							sub_covariance(index_off, h, pair,pairpr, lagbeg,lagend, omni);
						} else if(Math.abs(ivtype)==CORRELOGRAM) {
							sub_correlogram(index_off, h, pair,pairpr, lagbeg,lagend, omni);
						} else if(ivtype==PAIR_RELATIVE) {
							sub_pairwise_relative(index_off, h, pair,pairpr, lagbeg,lagend, omni);
						} else if(ivtype==VARIO_OF_LOG) {
							sub_log_variogram(index_off,h, pair,pairpr, lagbeg,lagend, omni);
						} else if(ivtype==SEMI_MADOGRAM) {
							sub_madogram(index_off, h, pair,pairpr, lagbeg,lagend, omni);
						}
					}
					progressNPairs = (double) (i+1) / (double) datalength;
				}
				int ii = id*nvarg*(nlag+2) + iv*(nlag+2);
				for(int il=0; il<nlag+2; il++) {
					lagnum[ii+il] = il+1;
					dirnum[ii+il] = id+1;
				}
				progressNDir = (double) (id+1) / (double) ndir;
			}
			progressNVario = (double) (iv+1) / (double) nvarg;
		}
//        c
//        c Get average values for gam, hm, tm, hv, and tv, then compute
//        c the correct "variogram" measure:
//        c
		for(int iv=0; iv<nvarg; iv++) {
			double[] vario_param = variogram.get(iv);
			int ivtail = (int) vario_param[0]-1 + Constants.FIRST_IDX;
			int ivhead = (int) vario_param[1]-1 + Constants.FIRST_IDX;
			int ivtype = (int) vario_param[2];
			//double ivcut = vario_param[3];

			for(int id=0; id<ndir; id++) { // Loop ID: 5

				for(int il=0; il<nlag+2; il++) {
					int i = id*nvarg*(nlag+2)+iv*(nlag+2)+il;
					if(np[i]<=0d) continue;
					double rnum   = np[i];
					dis[i] /= rnum;
					gam[i] /= rnum;
					hm[i]  /= rnum;
					tm[i]  /= rnum;
					hv[i]  /= rnum;
					tv[i]  /= rnum;
//        c
//        c Attempt to standardize:
//        c
					if(isill==1) {
						if(ivtail==ivhead) {
							//if(DEBUG_MODE) System.out.println("sills["+names[iii]+"] = "+sills[iii]);
							if((ivtype==SEMI_VARIOGRAM || ivtype==CONT_INDICATOR) && dataframe.getSill(ivtail)>0d)
								gam[i] /= dataframe.getSill(ivtail);
						}
					}
//        c
//        c  1. report the semivariogram rather than variogram
//        c  2. report the cross-semivariogram rather than variogram
//        c  3. the covariance requires "centering"
//        c  4. the correlogram requires centering and normalizing
//        c  5. general relative requires division by lag mean
//        c  6. report the semi(pairwise relative variogram)
//        c  7. report the semi(log variogram)
//        c  8. report the semi(madogram)
//        c
					if(ivtype==SEMI_VARIOGRAM || ivtype==CROSS_VARIOGRAM) {
						gam[i] *= 0.5d;
					} else if(Math.abs(ivtype)==COVARIANCE) {
						gam[i] -= hm[i]*tm[i];
						if(ivtype<0) {
							if(dataframe.getSill(ivtail)<0d || dataframe.getSill(ivhead)<0d) {
								gam[i] = Double.NaN;
							} else {
								double variance = Math.sqrt(dataframe.getSill(ivtail) * dataframe.getSill(ivhead));
								gam[i] = variance - gam[i];
							}
						}
					} else if(Math.abs(ivtype)==CORRELOGRAM) {
						hv[i]  -= hm[i]*hm[i];
						if(hv[i]<0d) hv[i] = 0d;
						hv[i]  = Math.sqrt(hv[i]);
						tv[i]  -= tm[i]*tm[i];
						if(tv[i]<0d) tv[i] = 0d;
						tv[i]  = Math.sqrt(tv[i]);
						if(hv[i]*tv[i]<Constants.D_EPSLON) {
							gam[i] = 0d;
						} else {
							gam[i] = (gam[i]-hm[i]*tm[i])/(hv[i]*tv[i]);
						}
						if(ivtype<0) gam[i] = 1d - gam[i];

//        c Square "hv" and "tv" so that we return the variance:
						hv[i]  *= hv[i];
						tv[i]  *= tv[i];
					} else if(ivtype==GENERAL_RELATIVE) {
						double htave  = 0.5d*(hm[i]+tm[i]);
						htave  *= htave;
						if(htave<Constants.D_EPSLON) {
							gam[i] = 0d;
						} else {
							gam[i] /= htave;
						}
					} else if(ivtype==PAIR_RELATIVE) {
						gam[i] *= 0.5d;
					}
				}
			}
		}
		DataFrame resdf = new DataFrame();
		resdf.addColumn("dir_num", dirnum);
		resdf.addColumn("lag_num", lagnum);
		resdf.addColumn("distance",   dis);
		resdf.addColumn("gamma",      gam);
		resdf.addColumn("points",      np);
		resdf.addColumn("mean_head",   hm);
		resdf.addColumn("mean_tail",   tm);
		resdf.addColumn("variance_head", hv);
		resdf.addColumn("variance_tail", tv);
		stop = false;
		return resdf;
	}

	public void writeData(String out_path) { writeData(out_path,false); }
	public void writeData(String out_path, boolean to_console) {
//        c-----------------------------------------------------------------------
//        c
//        c                  Write Out the Results of GAMV3
//        c                  ******************************
//        c
//        c An output file will be written which contains each directional
//        c variogram ordered by direction and then variogram (the directions
//        c cycle fastest then the variogram number).  For each variogram there
//        c will be a one line description and then "nlag" lines with:
//        c
//        c        a) lag number (increasing from 1 to nlag)
//        c        b) separation distance
//        c        c) the "variogram" value
//        c        d) the number of pairs for the lag
//        c        e) the mean of the data contributing to the tail
//        c        f) the mean of the data contributing to the head
//        c        g) IF the correlogram - variance of tail values
//        c        h) IF the correlogram - variance of head values
//        c
//        c
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    use geostat
		String title="";
		int nvarg = variogram.size(),
			ndir  = direction.size();
//        c
//        c Loop over all the variograms that have been computed:
//        c
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out_path)))) {
			for(int iv=0; iv<nvarg; iv++) {
				double[] vario_param = variogram.get(iv);
				int ivtail = (int) vario_param[0];
				int ivhead = (int) vario_param[1];
				int ivtype = (int) Math.abs(vario_param[2]);
				title += variogramTypeById(ivtype)+" tail:"+dataframe.getVarname(ivtail-1+Constants.FIRST_IDX)+
						"  head:"+dataframe.getVarname(ivhead-1+Constants.FIRST_IDX);
//        c
//        c Loop over all the directions (note the direction in the title):
//        c
				for(int id=0; id<ndir; id++) {
					String titled = title + "  direction "+(id+1);
					bw.append(titled+"\n");
					if(to_console) System.out.println("\n"+title);
//        c
//        c Write out all the lags:
//        c
					String[][] out = null;
					boolean[] coldivs = null, rowdivs = null;
					if(to_console) {
						out = new String[nlag+3][(ivtype==4 ? 8 : 6)];
						out[0][0] = " Lag num. "; out[0][1] = " distance ";
						out[0][2] = " gamma "; out[0][3] = " #points ";
						out[0][4] = " mean_h "; out[0][5] = " mean_t ";
						if(ivtype==4) { out[0][6]=" variance_h "; out[0][7]=" variance_t "; }
						coldivs = (ivtype==4 ? new boolean[]{true,false,true,true,false,false,false}
											 : new boolean[]{true,false,true,true,false});
						rowdivs = new boolean[nlag+2]; rowdivs[0] = true;
						for(int r=1; r<nlag+2; r++) rowdivs[r] = false;
					}
					for(int il=0; il<nlag+2; il++) {
						int i = id*nvarg*(nlag+2)+iv*(nlag+2)+il;
						int nump = (int) np[i];
						if(ivtype==4) {
							String line = " "+FormatHelper.nf(il+1,3)+
									      " "+FormatHelper.nf(dis[i],12,3)+
									      " "+FormatHelper.nf(gam[i],12,5)+
									      " "+FormatHelper.nf(nump,8)+
									      " "+FormatHelper.nf(hm[i],14,5)+
									      " "+FormatHelper.nf(tm[i],14,5)+
									      " "+FormatHelper.nf(hv[i],14,5)+
									      " "+FormatHelper.nf(tv[i],14,5);
							bw.append(line+"\n");
							if(to_console) {
								out[il+1][6] = " "+hv[i]+" "; out[il+1][7] = " "+tv[i]+" ";
							}
						} else {
							String line = " "+FormatHelper.nf(il+1,3)+
								      " "+FormatHelper.nf(dis[i],12,3)+
								      " "+FormatHelper.nf(gam[i],12,5)+
								      " "+FormatHelper.nf(nump,8)+
								      " "+FormatHelper.nf(hm[i],14,5)+
								      " "+FormatHelper.nf(tm[i],14,5);
							bw.append(line+"\n");
						}
						if(to_console) {
							out[il+1][0] = " "+(il+1)+" "; out[il+1][3] = " "+nump+" ";
							out[il+1][1] = " "+dis[i]+" "; out[il+1][2] = " "+gam[i]+" ";
							out[il+1][4] = " "+hm[i]+" "; out[il+1][5] = " "+tm[i]+" ";
						}
					}
					if(to_console) {
						FormatHelper.printTable(out, coldivs, rowdivs);
					}
				}
//        c
//        c End loop over variograms
//        c
			}
			bw.flush();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public void check(double[] vrmin, double[] vrmax) { 
//        c-----------------------------------------------------------------------
//        c
//        c                Error Check and Note Variogram Types
//        c                ************************************
//        c
//        c Go through each variogram type and note the type to the screen and
//        c report any possible errors.
//        c
//        c
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    use geostat
		String title="";
//        c
//        c Loop over all the variograms to be computed:
//        c
		System.out.println(" ");
		int nvarg = variogram.size();
		for(int iv=0; iv<nvarg; iv++) {
			double[] vario_param = variogram.get(iv);
			int ivtail = (int) vario_param[0]-1;
			int ivhead = (int) vario_param[1]-1;
			int ivtype = (int) vario_param[2];
			//double ivcut = vario_param[3];

//        c Note the variogram type and the variables being used:
			title = variogramTypeById(ivtype)+"  tail="+dataframe.getVarname(ivtail+Constants.FIRST_IDX)+
					"  head="+dataframe.getVarname(ivhead+Constants.FIRST_IDX);
			System.out.println(" Variogram "+(iv+1)+":  "+title);
//        c
//        c Check for possible errors or inconsistencies:
//        c
			if(ivtype==2) {
				if(ivtail==ivhead)
					System.out.println("  WARNING: cross variogram with the same variable!");
			} else if(ivtype==5) {
				if(ivtail!=ivhead)
					System.out.println("  WARNING: cross general relative variogram are difficult to interpret!");
				if(vrmin[ivtail]<0d && vrmax[ivtail]>0d)
					System.out.println("  WARNING: there are both positive and negative values - lag mean could be zero!");
				if(vrmin[ivhead]<0d && vrmax[ivhead]>0d)
					System.out.println("  WARNING: there are both positive and negative values - lag mean could be zero!");
			} else if(ivtype==6) {
				if(ivtail!=ivhead)
					System.out.println("  WARNING: cross pairwise relative variogram are difficult to interpret!");
				if(vrmin[ivtail]<0d && vrmax[ivtail]>0d)
					System.out.println("  WARNING: there are both positive and negative values - pair means could be zero!");
				if(vrmin[ivhead]<0d && vrmax[ivhead]>0d)
					System.out.println("  WARNING: there are both positive and negative values - pair means could be zero!");
			} else if(ivtype==7) {
				if(ivtail!=ivhead)
					System.out.println("  WARNING: cross logarithmic variograms may be difficult to interpret!");
				if(vrmin[ivtail]<0d || vrmin[ivhead]<0d)
					System.out.println("  WARNING: there are zero or negative values - logarithm undefined!");
			} else if(ivtype==8) {
				if(ivtail!=ivhead)
					System.out.println("  WARNING: cross rodograms may be difficult to interpret!");
			} else if(ivtype==9) {
				if(ivtail!=ivhead)
					System.out.println("  WARNING: cross madograms may be difficult to interpret!");
			}
//        c
//        c END Loop over all variograms:
//        c
		}
	}

	public void makepar() {
//        c-----------------------------------------------------------------------
//        c
//        c                      Write a Parameter File
//        c                      **********************
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File("res/gamv.par")))) {
			bw.append("                  Parameters for GAMV\n"
			         +"                  *******************\n"
		             +"START OF PARAMETERS:\n"
			         +"res/cluster.dat                   -file with data\n"
		             +"1   2   0                         -columns for X, Y, Z coordinates\n"
		             +"2   3   4                         -number of variables,col numbers\n"
		             +"-1.0e21     1.0e21                -trimming limits\n"
		             +"res/gamv.out                      -file for variogram output\n"
		             +"10                                -number of lags\n"
		             +"5.0                               -lag separation distance\n"
		             +"3.0                               -lag tolerance\n"
		             +"1                                 -number of directions\n"
		             +"0.0  90.0 50.0   0.0  90.0  50.0  -azm,atol,bandh,dip,dtol,bandv\n"
		           //+"0.0  22.5 25.0   0.0  22.5  25.0  -azm,atol,bandh,dip,dtol,bandv\n"
		           //+"90.  22.5 25.0   0.0  22.5  25.0  -azm,atol,bandh,dip,dtol,bandv\n"
		             +"0                                 -standardize sills? (0=no, 1=yes)\n"
		             +"1                                 -number of variograms\n"
		             +"1   1   1                         -tail var., head var., variogram type\n"
		           //+"1   2   2                         -tail var., head var., variogram type\n"
		           //+"2   2   1                         -tail var., head var., variogram type\n"
		              );
			bw.append("\n\n\n"
			         +"type 1 = traditional semivariogram\n"
			         +"     2 = traditional cross semivariogram\n"
			         +"     3 = covariance\n"
			         +"     4 = correlogram\n"
			         +"     5 = general relative semivariogram\n"
			         +"     6 = pairwise relative semivariogram\n"
			         +"     7 = semivariogram of logarithms\n"
			         +"     8 = semimadogram\n"
			         +"     9 = indicator semivariogram - continuous\n"
			         +"     10= indicator semivariogram - categorical\n"
			          );
			bw.flush();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}


	private String variogramTypeById(int _vario_id) {
		switch(Math.abs(_vario_id)) {
			case  1: return "Semivariogram                  ";
			case  2: return "Cross Semivariogram            ";
			case  3: return "Covariance                     ";
			case  4: return "Correlogram                    ";
			case  5: return "General Relative               ";
			case  6: return "Pairwise Relative              ";
			case  7: return "Variogram of Logarithms        ";
			case  8: return "Semimadogram                   ";
			case  9: return "Indicator 1/2 Variogram (cont.)";
			case 10: return "Indicator 1/2 Variogram (cat.) ";
		}
		System.err.println("No known variogram type for id "+_vario_id);
		return null;
	}

	private double[] sub_tail_head(int _i, int _j, double[] _tail, double[] _head, boolean _forward, boolean _omni, int _vario_type) {
		double[] sth = new double[4];
		if(_forward) {
			sth[0] = _head[_j];
			sth[1] = _tail[_i];
			if(_omni || _vario_type==2) {
				sth[2] = _head[_i];
				sth[3] = _tail[_j];
			}
		} else {
			sth[0] = _head[_i];
			sth[1] = _tail[_j];
			if(_omni || _vario_type==2) {
				sth[2] = _head[_j];
				sth[3] = _tail[_i];
			}
		}
		return sth;
	}

	private void sub_semivariogram(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			np[ii]  += 1d;
			dis[ii] += _h;
			tm[ii]  += _pair[0];
			hm[ii]  += _pair[1];
			gam[ii] += (_pair[1]-_pair[0])*(_pair[1]-_pair[0]);
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					np[ii]  += 1d;
					dis[ii] += _h;
					tm[ii]  += _pair_pr[0];
					hm[ii]  += _pair_pr[1];
					gam[ii] += (_pair_pr[1]-_pair_pr[0])*(_pair_pr[1]-_pair_pr[0]);
				}
			}
		}
	}
	private void sub_cross_semivariogram(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			np[ii]  += 1d;
			dis[ii] += _h;
			tm[ii]  += 0.5d*(_pair[0]+_pair_pr[0]);
			hm[ii]  += 0.5d*(_pair[1]+_pair_pr[1]);
			gam[ii] += (_pair_pr[1]-_pair[1])*(_pair[0]-_pair_pr[0]);
		}
	}
	private void sub_covariance(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			np[ii]  += 1d;
			dis[ii] += _h;
			tm[ii]  += _pair[0];
			hm[ii]  += _pair[1];
			gam[ii] += _pair[1]*_pair[0];
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					np[ii]  += 1d;
					dis[ii] += _h;
					tm[ii]  += _pair_pr[0];
					hm[ii]  += _pair_pr[1];
					gam[ii] += _pair_pr[1]*_pair_pr[0];
				}
			}
		}
	}
	private void sub_correlogram(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			np[ii]  += 1d;
			dis[ii] += _h;
			tm[ii]  += _pair[0];
			hm[ii]  += _pair[1];
			hv[ii]  += _pair[1]*_pair[1];
			tv[ii]  += _pair[0]*_pair[0];
			gam[ii] += _pair[1]*_pair[0];
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					np[ii]  += 1d;
					dis[ii] += _h;
					tm[ii]  += _pair_pr[0];
					hm[ii]  += _pair_pr[1];
					hv[ii]  += _pair_pr[1]*_pair_pr[1];
					tv[ii]  += _pair_pr[0]*_pair_pr[0];
					gam[ii] += _pair_pr[1]*_pair_pr[0];
				}
			}
		}
	}
	private void sub_pairwise_relative(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			if(Math.abs(_pair[0]+_pair[1])>=Constants.D_EPSLON) {
				np[ii]  += 1d;
				dis[ii] += _h;
				tm[ii]  += _pair[0];
				hm[ii]  += _pair[1];
				double gamma = 2d*(_pair[0]-_pair[1])/(_pair[0]+_pair[1]);
				gam[ii] += gamma*gamma;
			}
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					if(Math.abs(_pair_pr[0]+_pair_pr[1])>Constants.D_EPSLON) {
						np[ii]  += 1d;
						dis[ii] += _h;
						tm[ii]  += _pair_pr[0];
						hm[ii]  += _pair_pr[1];
						double gamma = 2d*(_pair[0]-_pair[1])/(_pair[0]+_pair[1]);
						gam[ii] += gamma*gamma;
					}
				}
			}
		}
	}
	private void sub_log_variogram(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			if(_pair[0]>Constants.D_EPSLON && _pair[1]>Constants.D_EPSLON) {
				np[ii]  += 1d;
				dis[ii] += _h;
				tm[ii]  += _pair[0];
				hm[ii]  += _pair[1];
				double gamma = Math.log(_pair[0])-Math.log(_pair[1]);
				gam[ii] += gamma*gamma;
			}
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					if(_pair_pr[0]>Constants.D_EPSLON && _pair_pr[1]>Constants.D_EPSLON) {
						np[ii]  += 1d;
						dis[ii] += _h;
						tm[ii]  += _pair_pr[0];
						hm[ii]  += _pair_pr[1];
						double gamma = Math.log(_pair[0])-Math.log(_pair[1]);
						gam[ii] += gamma*gamma;
					}
				}
			}
		}
	}
	private void sub_madogram(int _idxoff, double _h, double[] _pair, double[] _pair_pr, int _lagbeg, int _lagend, boolean _omni) {
		for(int il=_lagbeg-1; il<_lagend; il++) {
			int ii = _idxoff+il;
			np[ii]  += 1d;
			dis[ii] += _h;
			tm[ii]  += _pair[0];
			hm[ii]  += _pair[1];
			gam[ii] += Math.abs(_pair[1]-_pair[0]);
			if(_omni) {
				if(_pair_pr[0]>=tmin && _pair_pr[1]>=tmin && _pair_pr[0]<tmax && _pair_pr[1]<tmax) {
					np[ii]  += 1d;
					dis[ii] += _h;
					tm[ii]  += _pair_pr[0];
					hm[ii]  += _pair_pr[1];
					gam[ii] += Math.abs(_pair_pr[1]-_pair_pr[0]);
				}
			}
		}
	}

	private int[] addIntToArray(int[] arr, int new_entry) {
		int alen = arr.length;
		int[] cop = new int[alen];
		for(int i=0; i<alen; i++) cop[i] = arr[i];
		arr = new int[alen+1];
		for(int i=0; i<alen; i++) arr[i] = cop[i];
		arr[alen] = new_entry;
		return arr;
	}
}
