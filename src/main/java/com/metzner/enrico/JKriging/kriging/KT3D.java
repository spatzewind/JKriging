package com.metzner.enrico.JKriging.kriging;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.LinEquSolver;
import com.metzner.enrico.JKriging.helper.MathHelper;
import com.metzner.enrico.JKriging.probability.Covariance;

public class KT3D {

	public static final double VERSION = 3.000d;
	public static final double EPSLON = 0.000001d;
	private static final int MAXDT = 9; // maximal number of drift terms
	private static final int MAXNST = 4; // maximum number of nested variogram model structures
	private static final int MAXROT = MAXNST + 1; // maximum number of rotation matrices
	
	public static final int GRID      = 0;
	public static final int CROSS     = 1;
	public static final int JACKKNIFE = 2;
	public static final int SIMPLE_KRIGING                 = 0;
	public static final int ORDINARY_KRIGING               = 1;
	public static final int SIMPLE_NON_STATIONARY_KRIGING  = 2;
	public static final int EXTERNAL_DRIFT_KRIGING         = 3;
	public static final int INDEX_KRIGING_NO  = 0;
	public static final int INDEX_KRIGING_YES = 1;
	
	private int iktype, ncut, koption, idbg;
	private boolean[] withSAHE = new boolean[MAXNST];
	private int[] idrif = new int[MAXDT], it = new int[MAXNST];
	private double[] bv = new double[9], cc = new double[MAXNST], aa = new double[MAXNST],
			ang1 = new double[MAXNST], ang2 = new double[MAXNST], ang3 = new double[MAXNST], anis1 = new double[MAXNST],
			anis2 = new double[MAXNST];
	private double[][][] rotmat = new double[MAXROT][3][3];

	// module geostat
	private double[] cut, cdf;
	// end module geostat

	private String dbgfl;
	private double tmin, tmax;
	private int nx, ny, nz, nxdis, nydis, nzdis, ndmin, ndmax, noct, ktype, itrend;
	private double xmn, xsiz, ymn, ysiz, zmn, zsiz, radius, skmean;//, sanis1, sanis2, sang1, sang2, sang3;

	private int MAXCUT; // MXSXY,MXSX;

	int idhlj, ixlj, iylj, izlj, ivrlj, iextvj, iextve, nvarij;

	private DataFrame dataframe, jackdf, externdf;
	private String x_var, y_var, z_var, vr_var, dh_var, ext_var, ext_var_e;
	private String jack_x_var, jack_y_var, jack_z_var, jack_vr_var, jack_dh_var, jack_ext_var;
	private List<double[]> variograms;
	int num_krig_res; // number of kriging output variables (2 normal, 7 indicator)
	String[] estimate_titles; // titles for results from kriging
	private double[][] estimates; // result from kriging
	
	private int numOctants; //number of octants minimum has to distribute search
	
	private double progressNPoints;

	private boolean[] paramchecklist = new boolean[8];
	String[] param_descriptions = { "dataframe", "co-ordinate variables and kriging variable", "kriging field/points",
			"kriging search parameters", "debugging setting", "usage of jackknife file",
			"external drift file and parameters", "variogram model(s)" };

	private int numberOfWorker; //used for multithreading
	private AtomicInteger indexOffset;

	public KT3D() {
		variograms = new ArrayList<double[]>();
		indexOffset = new AtomicInteger(0);
		reset();
	}

	public void reset() {
		if (dataframe != null) {
			dataframe.clear();
			dataframe = null;
		}
		if (jackdf != null) {
			jackdf.clear();
			jackdf = null;
		}
		if(externdf != null) {
			externdf.clear();
			externdf = null;
		}
		variograms.clear();

		ncut = 1;
		idbg = 0;
		itrend = 0;
		num_krig_res = 0;
		estimates = null;
		estimate_titles = null;
		
		numOctants = 1;
		numberOfWorker = 1;

		paramchecklist[0] = false; // dataframe is declared
		paramchecklist[1] = false; // variable-coordinates and depending variables are defined
		paramchecklist[2] = false; // kriging field/point are set
		paramchecklist[3] = false; // kriging search parameters are set
		paramchecklist[4] = true;  // debugging level/file is choosen
		paramchecklist[5] = false; // jackknife file or unusage is defined
		paramchecklist[6] = false; // external drift is set/unset
		paramchecklist[7] = false; // variogram model ist set;
	}

	public KT3D setKrigingOptionAndDataframe(DataFrame _data, int kriging_option, int kriging_type, double simplekriging_mean, boolean use_indicator_kriging, double[] category_upper_bounds) {
		dataframe = _data;
		if(dataframe==null) {
			System.err.println("The dataframe does not exist!");
		} else {
			paramchecklist[0] = true;
			koption = Math.abs(kriging_option);
			System.out.println(" kriging option = " + koption);
			ktype = kriging_type;
			skmean = simplekriging_mean;
			System.out.println(" ktype, skmean = " + ktype + " " + skmean);
			iktype = (use_indicator_kriging ? Constants.IYES : Constants.INO);
			if(iktype==Constants.INO) {
				if(kriging_option==0) {
					paramchecklist[5] = true;
					paramchecklist[6] = true;
					num_krig_res = 3;
					estimate_titles = new String[] {"Estimate", "EstimateVariance", "NumberOfBasepoints"};
				} else {
					num_krig_res = 8;
					estimate_titles = new String[] {"X","Y","Z","True","Estimate","EstimateVariance","Error: est-true","NumberOfBasepoints"};
				}
			} else {
				ncut = category_upper_bounds.length;
				cut = new double[ncut];
				cdf = new double[ncut];
				num_krig_res = ncut;
				estimate_titles = new String[ncut];
				String str_cutoffs = "";
				for(int d=0; d<ncut; d++) {
					cut[d] = category_upper_bounds[d];
					str_cutoffs += " " + cut[d];
					estimate_titles[d] = "Category "+(d+1);
				}
				System.out.println(" number of cutoffs = " + ncut);
				System.out.println(" cutoffs =" + str_cutoffs);
			}
			
		}
		return this;
	}
	public KT3D setVariableAndCoordinates(int _x_id, int _y_id, int _z_id, int _var_id, int _dh_id) {
		if(dataframe==null) {
			System.err.println("First set a dataframe which may contain the variables you specify!");
		} else {
			x_var  = dataframe.getVarname(_x_id-1);
			y_var  = dataframe.getVarname(_y_id-1);
			z_var  = dataframe.getVarname(_z_id-1);
			vr_var = dataframe.getVarname(_var_id-1);
			dh_var = dataframe.getVarname(_dh_id-1);
			paramchecklist[1] = true;
			if(paramchecklist[0] && koption==1) {
				jackdf = dataframe;
				jack_dh_var = dh_var;
				jack_x_var  = x_var;
				jack_y_var  = y_var;
				jack_z_var  = z_var;
				jack_vr_var = vr_var;
				paramchecklist[5] = true;
			}
		}
		return this;
	}
	public KT3D setVariableAndCoordinates(String x_name, String y_name, String z_name, String var_name, String dh_name) {
		if(dataframe==null) {
			System.err.println("First set a dataframe which may contain the variables you specify!");
		} else {
			setVariableAndCoordinates(
					dataframe.getVariableID(x_name)+1,
					dataframe.getVariableID(y_name)+1,
					dataframe.getVariableID(z_name)+1,
					dataframe.getVariableID(var_name)+1,
					dataframe.getVariableID(dh_name)+1);
		}
		return this;
	}
	public KT3D setKrigingField(int num_x_points, double x_cell_width, double x_start,
			int num_y_points, double y_cell_width, double y_start,
			int num_z_points, double z_cell_width, double z_start) {
		nx = num_x_points; xmn = x_start; xsiz = x_cell_width;
		ny = num_y_points; ymn = y_start; ysiz = y_cell_width;
		nz = num_z_points; zmn = z_start; zsiz = z_cell_width;
		System.out.println(" nx, xmn, xsiz = " + nx + " " + xmn + " " + xsiz);
		System.out.println(" ny, ymn, ysiz = " + ny + " " + ymn + " " + ysiz);
		System.out.println(" nz, zmn, zsiz = " + nz + " " + zmn + " " + zsiz);
		paramchecklist[2] = true;
		return this;
	}
	public KT3D setSuperBlockParameters(int x_discretisation, int y_discretisation, int z_discretisation, int min_points, int max_points,
			int max_per_octant, double search_radius) {
		nxdis = x_discretisation; nydis = y_discretisation; nzdis = z_discretisation;
		ndmin = min_points; ndmax = max_points; noct = max_per_octant;
		if(max_per_octant>0) {
			numOctants = 8;
			System.out.println("WARNING: \"max per octant\" is set, it will override \"minimum search octant\" to 8.");
		}
		radius = search_radius; //sanis1 = second_radius / search_radius; sanis2 = third_radius / search_radius;
		//sang1 = anisotropy_angle; sang2 = second_angle; sang3 = third_angle;
		paramchecklist[3] = true;
		return this;
	}
	public KT3D setMinimumSearchOctant(int min_search_octant) {
		numOctants = min_search_octant;
		return this;
	}
	public KT3D setDebuggingLevel(int debugging_level, String debugging_file_path) {
		idbg  = debugging_level;
		dbgfl = debugging_file_path;
		paramchecklist[4] = true;
		File dbf = new File(dbgfl);
		if(idbg>0 && dbf.exists()) {
			try (BufferedWriter bw = new BufferedWriter(new FileWriter(dbf))) {
				bw.append(" ... init DEBUG File ...");
				bw.flush();
			} catch (IOException io_e) {
				io_e.printStackTrace();
			}
		}
		return this;
	}
	public KT3D setJackknife(DataFrame jack_dataframe, int _x_id, int _y_id, int _z_id, int _var_id, int _dh_id) {
		jackdf = jack_dataframe;
		if(jack_dataframe==null) {
			System.err.println("The dataframe with jackknife-data does not exist!");
		} else {
			jack_dh_var = jackdf.getVarname(_dh_id);
			jack_x_var  = jackdf.getVarname(_x_id);
			jack_y_var  = jackdf.getVarname(_y_id);
			jack_z_var  = jackdf.getVarname(_z_id);
			jack_vr_var = jackdf.getVarname(_var_id);
			paramchecklist[5] = true;
			if(paramchecklist[0] && koption==1) {
				jackdf = dataframe;
				jack_dh_var = dh_var;
				jack_x_var  = x_var;
				jack_y_var  = y_var;
				jack_z_var  = z_var;
				jack_vr_var = vr_var;
				paramchecklist[5] = (dataframe!=null);
			}
		}
		return this;
	}
	public KT3D setDriftOption(int _lin_X, int _lin_Y, int _lin_Z, int _quad_X, int _quad_Y, int _quad_Z, int _X_Y, int _X_Z, int _Y_Z) {
		idrif[0] = _lin_X; idrif[1] = _lin_Y; idrif[2] = _lin_Z;
		idrif[3] = _quad_X; idrif[4] = _quad_Y; idrif[5] = _quad_Z;
		idrif[6] = _X_Y; idrif[7] = _X_Z; idrif[8] = _Y_Z;
		paramchecklist[6] = true;
		return this;
	}
	public KT3D setExternalDrift(DataFrame extern_dataframe, int _ext_id) {
		externdf = extern_dataframe;
		if(externdf == null) {
			System.err.println("Dataframe with external drift data does not exist!");
		} else {
			ext_var_e = externdf.getVarname(_ext_id-1+Constants.FIRST_IDX);
			paramchecklist[6] = true;
		}
		return this;
	}
	public KT3D setExternalDrift(DataFrame extern_dataframe, String extern_variable_name) {
		if(extern_dataframe==null) {
			System.err.println("Dataframe with external drift data does not exist!");
		} else {
			setExternalDrift(extern_dataframe, extern_dataframe.getVariableID(extern_variable_name)+1-Constants.FIRST_IDX);
		}
		return this;
	}
	public KT3D addVariogramModelGeo(int variogram_type, double covariance,
			double azimuth, double dip, double roll, double h_max, double h_min, double h_vert) {
		return addVariogramModelMath(variogram_type, covariance, 90d-azimuth, dip, roll, h_max, h_min, h_vert);
	}
	public KT3D addVariogramModelMath(int variogram_type, double covariance,
			double azimuth, double dip, double roll, double h_max, double h_min, double h_vert) {
		if(h_max<0d) {
			System.err.println("search radius must be greater than zero!");
			return this;
		}
		if (variogram_type == Covariance.VARIOGRAM_POWER && (h_max < 0d || h_min > 2d)) {
			System.err.println("INVALID power variogram!");
			return this;
		}
		if (variogram_type == Covariance.VARIOGRAM_SINGLE_AXIS_HE) {
			if (variograms.size()==0) {
				System.err.println("There is no previous variogram to apply the single-axis hole effec!");
				return this;
			}
			double[] pre_vario = variograms.get(variograms.size()-1);
			if(Covariance.VARIOGRAM_HOLE_EFFECT == (int) pre_vario[0]) {
				System.err.println("Cannot apply single-axis hole effect to \"hole effect\"-variogram model!");
				return this;
			}
			if(pre_vario[0]-(int)pre_vario[0] > 0.5d) {
				System.err.println("Cannot apply second single-axis hole effect to variogram model!");
				return this;
			}
			variograms.set(variograms.size()-1, setSingleAxisHoleEffect(pre_vario,
					azimuth,dip,roll, h_max, h_min, h_vert));
			return this;
		}
		double anis1 = h_min / Math.max(h_max, EPSLON);
		double anis2 = h_vert / Math.max(h_max, EPSLON);
		System.out.println("search anisotropy angles = " + azimuth + " " + dip + " " + roll);
		System.out.println(" a1 a2 a3 =  " + h_max + " " + h_min + " " + h_vert);
		variograms.add(new double[] {variogram_type+0.0001d, covariance, azimuth, dip, roll, h_max, anis1, anis2});
		paramchecklist[7] = true;
		return this;
	}
	public KT3D setThreadNumber(int max_number_of_threads_to_use) {
		numberOfWorker = max_number_of_threads_to_use;
		return this;
	}
	public void loadParameterFile(String param_path, DataFrame _data, DataFrame _dfjack, DataFrame _dfextern) {
		reset();
		dataframe = _data;
		paramchecklist[0] = true;
		// c-----------------------------------------------------------------------
		// c
		// c Initialisation and Read Parameters
		// c **********************************
		// c
		// c The input parameters and data are read in from their files. Some quick
		// c error checking is performed and the statistics of all the variables
		// c being considered are written to standard output.
		// c
		// c
		// c
		// c-----------------------------------------------------------------------
		// use msflib
		// use geostat
		// include 'kt3d.inc'
		// parameter(MV=100)
		int idhl, ixl, iyl, izl, ivrl, iextv;
		double radius1, radius2, aa1, aa2;
		// c
		// c FORTRAN Units:
		// c
		// lin = 1
		// ldbg = 3
		// lout = 4
		// lext = 7
		// ljack = 8
		// c
		// c Note VERSION number:
		// c
		System.out.println(" KT3D Version: " + FormatHelper.nf(VERSION, 5, 3));
		// c
		// c Get the name of the parameter file - try the default name if no input:
		// c
		if (param_path == null) {
			param_path = "res/kt3d.par";
		} else if (param_path.length() == 0) {
			param_path = "res/kt3d.par";
		}
		File f = new File(param_path);
		if (!f.exists()) {
			System.err.println(
					"ERROR - the parameter file does not exist,\n" + "        check for the file and try again  ");
			if (param_path.equals("res/kt3d.par")) {
				System.out.println("        creating a blank parameter file");
				makepar();
			}
			return;
		}
		try (BufferedReader br = new BufferedReader(new FileReader(f))) {
			// c
			// c Find Start of Parameters:
			// c
			String line = br.readLine();
			while (!line.startsWith("STAR"))
				line = br.readLine();
			// c
			// c Read Input Parameters:
			// c
			line = br.readLine();
			String[] parts = FormatHelper.splitBySpace(line);
			// datafl = parts[0]; FormatHelper.chknam(datafl); System.out.println(" data
			// file = "+datafl);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			idhl = Integer.parseInt(parts[0].trim());
			ixl = Integer.parseInt(parts[1].trim());
			iyl = Integer.parseInt(parts[2].trim());
			izl = Integer.parseInt(parts[3].trim());
			ivrl = Integer.parseInt(parts[4].trim());
			iextv = Integer.parseInt(parts[5].trim());
			System.out.println(" columns = " + idhl + " " + ixl + " " + iyl + " " + izl + " " + ivrl + " " + iextv);
			dh_var  = dataframe.getVarname(idhl-1 );
			x_var   = dataframe.getVarname(ixl-1  );
			y_var   = dataframe.getVarname(iyl-1  );
			z_var   = dataframe.getVarname(izl-1  );
			vr_var  = dataframe.getVarname(ivrl-1 );
			ext_var = dataframe.getVarname(iextv-1);
			paramchecklist[1] = true;

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			tmin = Double.parseDouble(parts[0].trim());
			tmax = Double.parseDouble(parts[1].trim());
			System.out.println(" trimming limits = " + tmin + " " + tmax);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			koption = Integer.parseInt(parts[0].trim());
			System.out.println(" kriging option = " + koption);

			// c
			// c This is an undocumented feature to have kt3d construct an IK-type
			// c distribution:
			// c
			iktype = Constants.INO;
			if (koption < 0) {
				iktype = Constants.IYES;
				koption = -koption;
			}
			if (iktype == Constants.IYES) {
				line = br.readLine();
				parts = FormatHelper.splitBySpace(line);
				ncut = Integer.parseInt(parts[0].trim());
				System.out.println(" number of cutoffs = " + ncut);
				// c
				// c Find the needed parameter:
				// c
				// MAXCUT = ncut;
				// c
				// c Allocate the needed memory:
				// c21
				cut = new double[MAXCUT];
				// c22
				cdf = new double[MAXCUT];
				// c
				line = br.readLine();
				parts = FormatHelper.splitBySpace(line);
				String str_cutoffs = "";
				for (int i = 0; i < ncut; i++) {
					cut[i] = Double.parseDouble(parts[i].trim());
					str_cutoffs += " " + cut[i];
				}
				System.out.println(" cutoffs =" + str_cutoffs);
			}
			MAXCUT = ncut;

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			// jackfl = parts[0].trim(); FormatHelper.chknam(jackfl); System.out.println("jackknife data file = "+jackfl);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			ixlj = Integer.parseInt(parts[0].trim());
			iylj = Integer.parseInt(parts[1].trim());
			izlj = Integer.parseInt(parts[2].trim());
			ivrlj = Integer.parseInt(parts[3].trim());
			iextvj = Integer.parseInt(parts[4].trim());
			System.out.println(" columns = " + ixlj + " " + iylj + " " + izlj + " " + ivrlj + " " + iextvj);
			if (koption == 2)
				jackdf = _dfjack;
			paramchecklist[5] = (koption != 2 || jackdf != null);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			idbg = Integer.parseInt(parts[0].trim());
			System.out.println(" debugging level = " + idbg);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			dbgfl = parts[0].trim();
			FormatHelper.chknam(dbgfl);
			System.out.println(" debugging file = " + dbgfl);
			paramchecklist[4] = (idbg == 0 || dbgfl != null);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			// outfl = parts[0].trim(); FormatHelper.chknam(outfl); System.out.println("
			// output file = "+outfl);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nx = Integer.parseInt(parts[0].trim());
			xmn = Double.parseDouble(parts[1].trim());
			xsiz = Double.parseDouble(parts[2].trim());
			System.out.println(" nx, xmn, xsiz = " + nx + " " + xmn + " " + xsiz);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ny = Integer.parseInt(parts[0].trim());
			ymn = Double.parseDouble(parts[1].trim());
			ysiz = Double.parseDouble(parts[2].trim());
			System.out.println(" ny, ymn, ysiz = " + ny + " " + ymn + " " + ysiz);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nz = Integer.parseInt(parts[0].trim());
			zmn = Double.parseDouble(parts[1].trim());
			zsiz = Double.parseDouble(parts[2].trim());
			System.out.println(" nz, zmn, zsiz = " + nz + " " + zmn + " " + zsiz);
			paramchecklist[2] = true;

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nxdis = Integer.parseInt(parts[0].trim());
			nydis = Integer.parseInt(parts[1].trim());
			nzdis = Integer.parseInt(parts[2].trim());
			System.out.println(" block discretization: " + nxdis + " " + nydis + " " + nzdis);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ndmin = Integer.parseInt(parts[0].trim());
			ndmax = Integer.parseInt(parts[1].trim());
			System.out.println(" ndmin,ndmax = " + ndmin + " " + ndmax);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			noct = Integer.parseInt(parts[0].trim());
			System.out.println(" max per octant = " + noct);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			radius = Double.parseDouble(parts[0].trim());
			radius1 = Double.parseDouble(parts[1].trim());
			radius2 = Double.parseDouble(parts[2].trim());
			System.out.println(" search radii = " + radius + " " + radius1 + " " + radius2);
			if (radius < EPSLON) {
				System.err.println("radius must be greater than zero!");
				return;
			}
			// double radsqd = radius * radius;
			//sanis1 = radius1 / radius;
			//sanis2 = radius2 / radius;

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			//sang1 = Double.parseDouble(parts[0].trim());
			//sang2 = Double.parseDouble(parts[1].trim());
			//sang3 = Double.parseDouble(parts[2].trim());
			//System.out.println(" search anisotropy angles = " + sang1 + " " + sang2 + " " + sang3);
			paramchecklist[3] = true;

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ktype = Integer.parseInt(parts[0].trim());
			skmean = Double.parseDouble(parts[1].trim());
			System.out.println(" ktype, skmean = " + ktype + " " + skmean);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			String str_idrif = "";
			for (int i = 0; i < MAXDT; i++) {
				idrif[i] = Integer.parseInt(parts[i].trim());
				str_idrif += " " + idrif[i];
			} // MAXDT = 9
			System.out.println(" drift terms =" + str_idrif);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			itrend = Integer.parseInt(parts[0].trim());
			System.out.println(" itrend = " + itrend);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			//extfl = parts[0].trim(); FormatHelper.chknam(extfl); System.out.println(" external drift file = " + extfl);

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			iextve = Integer.parseInt(parts[0].trim());
			System.out.println(" variable in external drift file = " + iextve);
			paramchecklist[6] = (itrend == 0
					|| (itrend == 1 && (idrif[0] != 0 || idrif[1] != 0 || idrif[2] != 0 || idrif[3] != 0
							|| idrif[4] != 0 || idrif[5] != 0 || idrif[6] != 0 || idrif[7] != 0 || idrif[8] != 0)));
			if(ktype>=2) {
				externdf = _dfextern;
				if(externdf==null) {
					System.err.println("ERROR: dataframe with external drift data does not exist!");
					return;
				}
				ext_var_e = externdf.getVarname(iextve-1);
			}

			line = br.readLine();
			parts = FormatHelper.splitBySpace(line);
			double nst = Integer.parseInt(parts[0].trim());
			double c0 = Double.parseDouble(parts[1].trim());
			System.out.println(" nst, c0 = " + nst + " " + c0);

			if (nst <= 0) {
				System.err.println(" nst must be at least 1, it has been set to " + nst + "\n"
						+ " The c or a values can be set to zero");
				return;
			}
			variograms.clear();
			for (int i = 0; i < nst; i++) {
				line = br.readLine();
				parts = FormatHelper.splitBySpace(line);
				int it      = Integer.parseInt(parts[0].trim());
				double cc   = Double.parseDouble(parts[1].trim());
				double ang1 = Double.parseDouble(parts[2].trim());
				double ang2 = Double.parseDouble(parts[3].trim());
				double ang3 = Double.parseDouble(parts[4].trim());
				line = br.readLine();
				parts = FormatHelper.splitBySpace(line);
				double aa   = Double.parseDouble(parts[0].trim());
				aa1 = Double.parseDouble(parts[1].trim());
				aa2 = Double.parseDouble(parts[2].trim());
				double anis1 = aa1 / Math.max(aa, EPSLON);
				double anis2 = aa2 / Math.max(aa, EPSLON);
				System.out.println(
						" it,cc,ang[1,2,3] = " + it + " " + cc + " " + ang1 + " " + ang2 + " " + ang3);
				System.out.println(" a1 a2 a3 =  " + aa + " " + aa1 + " " + aa2);
				if (it == 4) {
					if (aa < 0d || aa > 2d) {
						System.err.println(" INVALID power variogram");
						return;
					}
				}
				variograms.add(new double[] {it+0.1d, cc, ang1,ang2,ang3, aa,anis1,anis2});
			}
			paramchecklist[7] = (nst > 0);
		} catch (IOException | NullPointerException io_np_e) {
			throw new RuntimeException("ERROR in parameter file!", io_np_e);
		}
		// c
		// c Perform some quick error checking:
		// c
		if (ktype == 3 && iextv <= 0) {
			System.err.println("must have external variable");
			return;
		}

		// c
		// c Open the debugging and output files:
		// c
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl)))) {
			bw.append(" ... init DEBUG File ...");
			bw.flush();
		} catch (IOException io_e) {
			io_e.printStackTrace();
		}
		num_krig_res = (iktype==0 ? (koption==0 ? 3 : 8) : ncut);
		if(iktype==0) { estimate_titles = (koption==0 ? new String[] {"Estimate", "EstimateVariance","NumberOfBasepoints"} :
			new String[] {"X","Y","Z","True","Estimate","EstimateVariance","Error: est-true","NumberOfBasepoints"});
		} else {
			estimate_titles = new String[ncut];
			for(int et=0; et<ncut; et++) estimate_titles[et] = "Category "+(et+1);
		}
//		try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfl)))) {
//			bw.append(title + "\n");
//			if (iktype == 0) {
//				if (koption == 0) {
//					bw.append(" " + FormatHelper.nf(2, 4) + " " + FormatHelper.nf(nx, 4) + " " + FormatHelper.nf(ny, 4)
//							+ " " + FormatHelper.nf(nz, 4) + "\n");
//					bw.append("Estimate\nEstimationVariance\n");
//				}
//				if (koption >= 1) {
//					bw.append("    7\n");
//					bw.append("X\nY\nZ\nTrue\nEstimate\nEstimationVariance\nError: est-true\n");
//				}
//			}
//			if (iktype == 1) {
//				if (koption == 0) {
//					bw.append(" " + FormatHelper.nf(ncut, 4) + " " + FormatHelper.nf(nx, 4) + " "
//							+ FormatHelper.nf(ny, 4) + " " + FormatHelper.nf(nz, 4) + "\n");
//				} else {
//					bw.append(" " + FormatHelper.nf(ncut + 1, 4) + "\n");
//				}
//				for (int i = 0; i < ncut; i++) {
//					bw.append(
//							"Threshold: " + FormatHelper.nf(i + 1, 2) + " = " + FormatHelper.nf(cut[i], 12, 5) + "\n");
//				}
//				if (koption == 1)
//					bw.append("true value\n");
//			}
//			bw.flush();
//		} catch (IOException io_e) {
//			io_e.printStackTrace();
//		}
		// c
		// c Open the external drift file if needed and position it at the
		// c first grid node in the file:
		// c
		if ((ktype == 2 || ktype == 3) && koption == 0) {
			/*
			 * theoretically there should be some work with an external drift file "extfl"
			 * but the original Fortran source code from
			 * http://www.statios.com/software/gslib90sc.zip and http://www.statios.com/GSLIB/index.html
			 * doesn't provide any meaningful work with this file
			 */
		}
		// c
		// c Set up for cross validation:
		// c
		if (koption == 1) {
			jackdf = dataframe;
			jack_dh_var = dh_var;
			jack_x_var = x_var;
			jack_y_var = y_var;
			jack_z_var = z_var;
			jack_vr_var = vr_var;
			jack_ext_var = ext_var;
		}
		// c
		// c Open the file with the jackknife data?
		// c
		if (koption > 0) {
			if (jackdf == null) {
				System.err.println("ERROR jackknife-dataframe does not exist!");
				paramchecklist[5] = false;
				return;
			}
		}
	}
	private double[] setSingleAxisHoleEffect(double[] pre_vario, double azimuth, double dip, double roll, double h_max, double h_min, double h_vert) {
		//calc cross section of ellipsoid (previous varigoram) and plane (perpendicular to he-axis)
		// vecX... is vector of ellipsoid, S... scale matrix, R... rotation matrix
		// vecN... is unit-vector of plane normal
		// vecU... vector of unit-sphere:
		// I:  vecU = S*R*x with |vecU| = 1  <-- sphere
		// II: vecN*(vecX) = 0        <-- plane
		// put sphere (I) into plane (II) equation:
		// vecN * (S^-1 * vecU) = (S^-T * vecN) * vecU = 0
		double[][] preRot = MathHelper.matmul(
				MathHelper.setrot3Dmath(pre_vario[2], pre_vario[3], pre_vario[4], pre_vario[6], pre_vario[7], 1, 1, new double[1][3][3])[0],
				1d/pre_vario[5]);
		double[][] Sm1t = { { pre_vario[5],              0d, 0d },
							{ 0d, pre_vario[5]*pre_vario[6], 0d },
							{ 0d, 0d, pre_vario[5]*pre_vario[7] } };
		double[] n0 = MathHelper.matmul(MathHelper.inverse(MathHelper.setrot3Dmath(azimuth, dip, roll, 1d, 1d, 1, 1, new double[1][3][3])[0]),
										new double[] {1d,0d,0d});
		// get vecM = (S^-T * vecN) / |S^-T * vecN| as hesse-form of plane equation
		double[] m = MathHelper.matmul(Sm1t, n0);
		// now get defining vectors for cross-ellipse: vecE1,vecE2
		// if vecM[z] = +/- 1 than: vecE1 = (1,0,0)T, vecE2 = (0,1,0)T
		// else:                    vecE1 = (vecM[y],-vecM[x],0)T / sqrt(vecM[x]²+vecM[y]²), vecE2 = vecM x vecE1
		double[] e1= {1d,0d,0d}, e2={0d,1d,0d};
		double mr = Math.sqrt(m[0]*m[0]+m[1]*m[1]);
		if(mr>1.0e-12d)
		{
			//set e1
			e1[0] = m[1]/mr; e1[1] = -m[0]/mr;
			//e2 = m x e1
			e2[0] = -m[2]*e1[1]/mr; e2[1] = m[2]*e1[0]/mr; e2[2] = (m[0]*e1[1]-m[1]*e1[0])/mr;
		}
		// back-transform to original space:
		// x = R^1 * S^-1 * vecU
		// --> vecF1 = (S*R)^-1 * vecE1
		//     vecF2 = (S*R)^-1 * vecE2
		double[][] infPreRot = MathHelper.inverse(preRot);
		double[] f1 = MathHelper.matmul(infPreRot, e1);
		double[] f2 = MathHelper.matmul(infPreRot, e2);
		// transform to "Scheitelform":
		// p(t) = -vecF1 * sin(t) + vecF2 * cos(t)
		// cot(2*t0) = (<vecF1,vecF1>-<vecF2,vecF2>) / (2 * <vecF1,vecF2>), with t0=0 if <vecF1,vecF2>=0
		double scp_f1f2 = MathHelper.matmul(MathHelper.transpose(f1), f2)[0];
		double magF1 = f1[0]*f1[0] + f1[1]*f1[1] + f1[2]*f1[2];
		double magF2 = f2[0]*f2[0] + f2[1]*f2[1] + f2[2]*f2[2];
		double t0 = 0d;
		if(Math.abs(scp_f1f2)>1.0e-12d)
		{
			t0 = 0.25d * Math.PI;
			if(Math.abs(magF1-magF2)>1.0e-12d)
				t0 = 0.5d * Math.atan(scp_f1f2 / (magF1-magF2));
		}
		// extrem - points are now:
		// p(t0+{-1,0,1,2}*pi/2)
		// or:
		// vecEx1 = +/- (vecF2 * cos(t0) - vecF1 * sin(t0))
		// vecEx2 = +/- (vecF2 * sin(t0) - vecF2 * cos(t0))
		double ct=Math.cos(t0), st=Math.sin(t0);
		double[] p0 = { f1[0]*ct+f2[0]*st, f1[1]*ct+f2[1]*st, f1[2]*ct+f2[2]*st };
		double[] p1 = { f2[0]*ct-f1[0]*st, f2[1]*ct-f1[1]*st, f2[2]*ct-f1[2]*st };
		//extract rotation matrix from new ellipse-extrem-points
		double p0r = Math.sqrt(p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2]);
		double p1r = Math.sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
		if(p0r<p1r) {
			double temp = p0r; p0r = p1r; p1r = temp;
			temp = p0[0]; p0[0] = p1[0]; p1[0] = temp;
			temp = p0[1]; p0[1] = p1[1]; p1[1] = temp;
			temp = p0[2]; p0[2] = p1[2]; p1[2] = temp;
		}
		double[] p2 = { p0[1]*p1[2]-p0[2]*p1[1], p0[2]*p1[0]-p0[0]*p1[2], p0[0]*p1[1]-p0[1]*p1[0] };
		double p2r = Math.sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]);
		double[][] invRM = { { p0[0]/p0r, p1[0]/p1r, p2[0]/p2r },
							 { p0[1]/p0r, p1[1]/p1r, p2[1]/p2r },
							 { p0[2]/p0r, p1[2]/p1r, p2[2]/p2r } };
		double[][] newRM = MathHelper.inverse(invRM);
		//get angles from rotation matrix
		double temp = Math.sqrt(newRM[0][0]*newRM[0][0]+newRM[0][1]*newRM[0][1]) * (newRM[0][0]<0d && newRM[0][1]<0d ? -1d : 1d);
		double alpha = Constants.RAD2DEG * Math.acos(Math.max(-1d,Math.min(1d,newRM[0][0]/temp))) * (newRM[0][1]/temp<0d ? -1d : 1d);
		if(alpha<-90d) alpha += 180d;
		if(alpha>90d) alpha -= 180d;
		double beta  = Constants.RAD2DEG * Math.asin(newRM[0][2]);
		//System.out.println("  ... \u03b1="+alpha+" and \u03b2="+beta);
		double[][] ab = MathHelper.setrot3Dmath(alpha, beta, 0d, 1d, 1d, 1, 1, new double[1][3][3])[0];
		//FormatHelper.printMat(System.out, ab);
		//System.out.println("{Det = "+MathHelper.determinante(ab)+"}");
		newRM = MathHelper.matmul(newRM, MathHelper.transpose(ab));
		if(newRM[0][0]<0d)
			newRM = MathHelper.matmul(newRM, -1d);
		//System.out.println("  ... resulting in matrix for gamma:");
		FormatHelper.printMat(System.out, newRM);
		System.out.println("{Det = "+MathHelper.determinante(newRM)+"}");
		double gamma = Constants.RAD2DEG * Math.acos(Math.min(1d,0.5d*(Math.abs(newRM[1][1])+Math.abs(newRM[2][2]))));
		gamma *= (Math.abs(newRM[1][2])<Math.abs(newRM[2][1]) ? -newRM[2][1] : newRM[1][2])<0d ? -1d : 1d;
//		mp[m][4] = alpha;
//		mp[m][5] = beta;
//		mp[m][6] = gamma;
		boolean check1 = newRM[0][0]>1d-1.e-10d;
		boolean check2 = Math.abs(newRM[0][1])<1.e-10d && Math.abs(newRM[0][2])<1.e-10d &&
						 Math.abs(newRM[1][0])<1.e-10d && Math.abs(newRM[2][0])<1.e-10d;
		boolean check3 = Math.abs(newRM[1][1]-newRM[2][2])<1.e-10d && Math.abs(newRM[1][2]+newRM[2][1])<1.e-10d
						&& newRM[1][1]>-1.e-10d;
		if(!check1 || !check2 || !check3)
		{
			//System.out.println("For check (k="+k3+"):");
			//FormatHelper.printMat(System.out, rm);
			if(!check1) System.err.println("Check for 1 top left failed.");
			if(!check2) System.err.println("Check for 0s top and left failed.");
			if(!check3) System.err.println("Check for cos and sin failed.");
		}
		pre_vario[0] += 0.5d;
		pre_vario[2] = alpha;
		pre_vario[3] = beta;
		pre_vario[4] = gamma;
		pre_vario[5] = p0r;
		pre_vario[6] = p1r / p0r;
		pre_vario[7] = h_max * Constants.SAHE_FAC / p0r;
		return pre_vario;
	}
	
	public double getProgress() {
		return progressNPoints;
	}

	public DataFrame kt3d() {
		if(!kt3d_df())
			return null;
		DataFrame result = new DataFrame();
		for(int n_est=0; n_est<num_krig_res; n_est++)
			result.addColumn(estimate_titles[n_est], estimates[n_est]);
		return result;
	}
	public DataFrame3D kt3d_asDataFrame3D() {
		if(!kt3d_df())
			return null;
		int xn=nx, yn=ny, zn=nz;
		if(xn*yn*zn!=estimates[0].length) { nx=estimates[0].length; ny=1; nz=1; }
		DataFrame3D result = new DataFrame3D();
		for(int n_est=0; n_est<num_krig_res; n_est++) {
			double[][][] estimates3D = new double[zn][yn][xn];
			for(int z_=0; z_<zn; z_++) for(int y_=0; y_<yn; y_++) for(int x_=0; x_<xn; x_++)
				estimates3D[z_][y_][x_] = estimates[n_est][xn*yn*z_ + xn*y_ + x_];
			result.addColumn(estimate_titles[n_est], estimates3D);
		}
		if(nx*ny*nz==estimates[0].length) {
			double[] dimension_x = new double[nx];
			for(int x_=0; x_<nx; x_++) dimension_x[x_] = xmn + x_*xsiz;
			double[] dimension_y = new double[ny];
			for(int y_=0; y_<ny; y_++) dimension_y[y_] = ymn + y_*xsiz;
			double[] dimension_z = new double[nz];
			for(int z_=0; z_<nz; z_++) dimension_z[z_] = zmn + z_*xsiz;
			result.setDimension(0+Constants.FIRST_IDX, dimension_z, z_var);
			result.setDimension(1+Constants.FIRST_IDX, dimension_y, y_var);
			result.setDimension(2+Constants.FIRST_IDX, dimension_x, x_var);
		}
		return result;
	}

	private boolean kt3d_df() { //TODO beginning of KT3D
		boolean someParamIsMissing = false;
		for(int pcl=0; pcl<8; pcl++) if(!paramchecklist[pcl]) {
			if(!someParamIsMissing)
				System.err.println("Not all parameters are set!");
			System.out.println("  \""+param_descriptions[pcl]+"\" missing");
			someParamIsMissing = true;
		}
		if(someParamIsMissing)
			return false;
		numOctants = Math.max(1, Math.min(8, numOctants));
		if(noct>0)
			numOctants = 8;
//        c-----------------------------------------------------------------------
//        c
//        c                Krige a 3-D Grid of Rectangular Blocks
//        c                **************************************
//        c
//        c This subroutine estimates point or block values of one variable by
//        c simple, ordinary, or kriging with a trend model.  It is also possible
//        c to estimate the trend directly.
//        c
//        c
//        c
//        c PROGRAM NOTES:
//        c
//        c   1. The data and parameters are passed in common blocks defined
//        c      in kt3d.inc.  Local storage is allocated in the subroutine
//        c      for kriging matrices, i.e.,
//        c         - xa,ya,za,vra   arrays for data within search neighborhood
//        c         - a,r,rr,s       kriging arrays
//        c         - xdb,ydb,zdb    relative position of discretization points
//        c         - cbb            block covariance
//        c   2. The kriged value and the kriging variance is written to Fortran
//        c      unit number "lout".
//        c
//        c
//        c
//        c
//        c Original:  A.G. Journel and C. Lemmer                             1981
//        c Revisions: A.G. Journel and C. Kostov                             1984
//        c-----------------------------------------------------------------------
		progressNPoints = 0d;
		double[] vr = (double[]) dataframe.getArray(vr_var);
		int datalength = vr.length;
		double[] x=null,y=null,z=null,dh=null,ve=null;
		if(x_var!=null)   x =  (double[]) dataframe.getArray(x_var);
		if(y_var!=null)   y =  (double[]) dataframe.getArray(y_var);
		if(z_var!=null)   z =  (double[]) dataframe.getArray(z_var);
		if(dh_var!=null)  dh = (double[]) dataframe.getArray(dh_var);
		if(ext_var!=null) ve = (double[]) dataframe.getArray(ext_var);
		if(x_var==null)   {x = new double[datalength];  for(int d=0; d<datalength; d++) x[d]  = xmn; }
		if(y_var==null)   {y = new double[datalength];  for(int d=0; d<datalength; d++) y[d]  = ymn; }
		if(z_var==null)   {z = new double[datalength];  for(int d=0; d<datalength; d++) z[d]  = zmn; }
		if(dh_var==null)  {dh = new double[datalength]; for(int d=0; d<datalength; d++) dh[d] = 0d;  }
		if(ext_var==null) {ve = new double[datalength]; for(int d=0; d<datalength; d++) ve[d] = 1d;  }
		if (x_var==null && nx > 1)
			System.out.println(" WARNING: no x variable and nx>1 !");
		if (y_var==null && ny > 1)
			System.out.println(" WARNING: no y variable and ny>1 !");
		if (z_var==null && nz > 1)
			System.out.println(" WARNING: no z variable and nz>1 !");
		
		int MAXDIS = nxdis*nydis*nzdis;
		ndmax = Math.max(ndmax, ndmin*numOctants);
		int maxSamples = ndmax + 1;
		int maxEquations = maxSamples + MAXDT + 2;
		//int maxSupBlckX = Math.max(1, Math.min(50, nx / 2));
		//int maxSupBlckY = Math.max(1, Math.min(50, ny / 2));
		//int maxSupBlckZ = Math.max(1, Math.min(50, nz / 2));
		//int maxSuperblocks = maxSupBlckX*maxSupBlckY*maxSupBlckZ;
		//MXSXY = 4 * MAXSBX * MAXSBY;
		//MXSX  = 2 * MAXSBX;
//		if (ndmax > maxSamples) {
//			System.err.println("ndmax is too big - modify .inc file");
//			return null;
//		}
		//int[]    nisb = new int[maxSuperblocks],
		//		 ixsbtosr = new int[8*maxSuperblocks],
		//		 iysbtosr = new int[8*maxSuperblocks],
		//		 izsbtosr = new int[8*maxSuperblocks];
		double[] xdb = new double[MAXDIS],
				 ydb = new double[MAXDIS],
				 zdb = new double[MAXDIS];
		int[]    ai = new int[maxEquations*maxEquations];
		
		int[] nst = {variograms.size(),0};
		double[] c0 = {1d,0d};
		double max_rad = radius;
		for(int i=0; i<nst[0]; i++) {
			double[] vario = variograms.get(i);
			it[i]       = (int) vario[0];
			withSAHE[i] = (vario[0]-it[i]>0.5d);
			cc[i]       = vario[1]; c0[0] -= cc[i];
			ang1[i]     = vario[2];
			ang2[i]     = vario[3];
			ang3[i]     = vario[4];
			aa[i]       = vario[5];
			anis1[i]    = vario[6];
			anis2[i]    = vario[7];
			if(it[i]==Covariance.VARIOGRAM_SPHERICAL && aa[i]<max_rad) max_rad = aa[i];
			if(withSAHE[i] && aa[i]*anis2[i]>max_rad) max_rad = aa[i]*anis2[i];
		}
		//to avoid including basepoints outside the range-ellipsoide:
		//reduces null-matrices;
		if(max_rad < radius) {
			System.out.println("[WARNING] Due to variogram settings a smaller search distance is used: "+max_rad+"!");
		}
		
		
		
		

		double cbb;
		//double[] sec3 = new double[datalength];
//        c
//        c Set up the rotation/anisotropy matrices that are needed for the
//        c variogram and search.  Also compute the maximum covariance for
//        c the rescaling factor:
//        c
		System.out.println("Setting up rotation matrices for variogram and search");
		double radsqd = radius * radius;
		double PMX    = 999.0d; //TODO why?
		double covmax = c0[0];
		double wgtsum = 0d;
		for(int is=0; is<nst[0]; is++) {
			rotmat = MathHelper.setrot3Dmath(ang1[is],ang2[is],ang3[is],anis1[is],anis2[is],is+1,MAXROT,rotmat);
			double cov_single = it[is]==Covariance.VARIOGRAM_POWER ? PMX : cc[is];
			covmax += cov_single;
			double wgt = cov_single / (aa[is]*Math.pow(anis1[is]*anis2[is], 1d/3d));
			for(int j=0; j<3; j++)
				for(int i=0; i<3; i++)
					rotmat[MAXNST][j][i] += rotmat[is][j][i]*wgt;
			wgtsum += wgt;
		}
		//int isrot = MAXNST + 1;
		//rotmat = MathHelper.setrot(sang1,sang2,sang3,sanis1,sanis1,isrot,MAXROT,rotmat);
		//searching ellipsoid is replaced by weighted variogram-range ellipsoids
		for(int j=0; j<3; j++)
			for(int i=0; i<3; i++)
				rotmat[MAXNST][j][i] /= wgtsum;
//        c
//        c Finish computing the rescaling factor and stop if unacceptable:
//        c
		double resc = 0d;
		if(radsqd<1d) {
			resc = 2d * radius / Math.max(covmax,0.0001d);
		} else {
			resc = 4d * radsqd / Math.max(covmax,0.0001d);
		}
		if(resc<=0d) {
			System.err.println("ERROR KT3D: The rescaling value is wrong "+resc+"\n"
					+          "            Maximum covariance: "+covmax+"\n"
					+          "            search radius:      "+radius);
			return false;
		}
		resc = 1d / resc;
//		System.out.println("[DEBUG] Maximum covariance: "+covmax);
//        c
//        c Set up for super block searching:
//        c
		//System.out.println("Setting up super block search strategy");
//		System.out.println("Setting up 3D-Tree instead of super block search strategy");
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(20, x,y,z,vr,tmp);
		//int nsec = 2;
		//double[] superblock_grid = new double[9]; // replace n[xyz]sup, [xyz]mnsup, [xyz]sizsup
		//nisb = DataHelper.setsupr(nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz, x, y, z, vr, tmp,
		//		nsec, ve, dh, sec3, maxSupBlckX, maxSupBlckY, maxSupBlckZ, nisb, superblock_grid);
		KdTree tree = new KdTree();
		tree.build(false, false, x,y,z);
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(20, x,y,z,vr,tmp);
		//int nsbtosr = DataHelper.picksupr(superblock_grid, isrot, rotmat, radsqd, ixsbtosr, iysbtosr, izsbtosr);
//        c
//        c Compute the number of drift terms, if an external drift is being
//        c considered then it is one more drift term, if SK is being considered
//        c then we will set all the drift terms off and mdt to 0):
//        c
		//System.out.println("[DEBUG] super block strategy set");
		int mdt = 1;
		for(int i=0; i<9; i++) {
			if(ktype==SIMPLE_KRIGING || ktype==SIMPLE_NON_STATIONARY_KRIGING) idrif[i] = 0;
			if(idrif[i]<0 || idrif[i]>1) {
				System.err.println("ERROR KT3D: invalid drift term "+idrif[i]);
				return false;
			}
			mdt += idrif[i];
		}
		if(ktype==EXTERNAL_DRIFT_KRIGING) mdt++;
		if(ktype==SIMPLE_KRIGING) mdt = 0;
		if(ktype==SIMPLE_NON_STATIONARY_KRIGING) mdt = 0;
//        c
//        c Set up the discretization points per block.  Figure out how many
//        c are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
//        c the offsets relative to the block center (this only gets done once):
//        c
//        c In all cases the offsets are relative to the lower left corner.
//        c This is done for rescaling the drift terms in the kriging matrix.
//        c
		if(nxdis<1) nxdis = 1;
		if(nydis<1) nydis = 1;
		if(nzdis<1) nzdis = 1;
		int ndb = nxdis * nydis * nzdis;
		if(ndb>MAXDIS) {
			System.err.println("ERROR KT3D: Too many discretization points "+ndb+"\n"
					+          "            Increase MAXDIS or lower n[xyz]dis");
			return false;
		}
		double xdis = xsiz  / Math.max(nxdis,1d);
		double ydis = ysiz  / Math.max(nydis,1d);
		double zdis = zsiz  / Math.max(nzdis,1d);
		int iii  = -1;
		double xloc = -0.5d*(xsiz+xdis),yloc = 0d,zloc = 0d;
		for(int ix=0; ix<nxdis; ix++) {
			xloc = xloc + xdis;
			yloc = -0.5d*(ysiz+ydis);
			for(int iy=0; iy<nydis; iy++) {
				yloc = yloc + ydis;
				zloc = -0.5d*(zsiz+zdis);
				for(int iz=0; iz<nzdis; iz++) {
					zloc = zloc + zdis;
					iii++;
					xdb[iii] = xloc + 0.5d*xsiz;
					ydb[iii] = yloc + 0.5d*ysiz;
					zdb[iii] = zloc + 0.5d*zsiz;
				}
			}
		}
//        c Initialize accumulators:
		//removed because of multithreading

//        c Calculate Block Covariance. Check for point kriging.
		double cov = Covariance.cova3(xdb[0], ydb[0], zdb[0], xdb[0], ydb[0], zdb[0],
				1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];

//        c Set the ``unbias'' variable so that the matrix solution is more stable
		double unbias = cov;
		cbb    = cov;
		if(ndb>1) {
			cbb = 0;
			for(int i=0; i<ndb; i++) {
				for(int j=0; j<ndb; j++) {
					cov = Covariance.cova3(xdb[i], ydb[i], zdb[i], xdb[j], ydb[j], zdb[j],
							1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
		//                  call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j),
		//     +               1,nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
					if(i==j) cov -= c0[0];
					cbb += cov;
				}
			}
			cbb /= ndb*ndb;
		}
		if(idbg>1) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				bw.append("\nBlock Covariance: "+cbb+"\n\n");
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}

//        c Mean values of the drift functions:
		for(int i=0; i<9; i++) bv[i] = 0d;
		for(int i=0; i<ndb; i++) {
			bv[0] += xdb[i]; bv[1] += ydb[i]; bv[2] += zdb[i];
			bv[3] += xdb[i]*xdb[i]; bv[4] += ydb[i]*ydb[i]; bv[5] += zdb[i]*zdb[i];
			bv[6] += xdb[i]*ydb[i]; bv[7] += xdb[i]*zdb[i]; bv[8] += ydb[i]*zdb[i];
		}
		for(int i=0; i<9; i++) bv[i] *= resc / ndb;

//        c Report on progress from time to time:
		int nxyz=1,nloop,irepo;
		int nd=vr.length;
		if(koption==0) {
			nxyz  = nx*ny*nz;
			nloop = nxyz;
			irepo = Math.max(1,Math.min(10000,nxyz/10));
		} else {
			nloop = 10000000;
			irepo = Math.max(1,Math.min(10000,nd/10));
		}
		estimates = new double[num_krig_res][nloop];
		System.out.println("\n\nWorking on the kriging\n");

//        c MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
		
		//main loop is split into multiple worker threads for parallelisation!!!
		indexOffset.set(0);
		numberOfWorker = Math.min(numberOfWorker, 1+(nloop-1)/irepo);
		Worker[] worker = new Worker[numberOfWorker];
		for(int t=0; t<numberOfWorker; t++) {
			worker[t] = new Worker(nloop, irepo, tree, maxEquations,maxSamples,
					unbias, covmax, cbb, ndb, c0, nst, resc, mdt,
					x,y,z,dh,ve,vr,xdb,ydb,zdb);
			worker[t].start();
		}
		//progressNPoints = (index+1d) / (double) nloop;
		//if((index+1)%irepo==0) System.out.println("   currently on estimate "+FormatHelper.nf(index+1,9));
		int currentIndex = 0;
		boolean workerUnfinished = true;
		while(workerUnfinished) {
			//wait
			workerUnfinished = false;
			for(Worker w: worker)
				if(!w.hasFinished())
					workerUnfinished = true;
			
			progressNPoints = (indexOffset.get()+1d) / (double) nloop;
			if(currentIndex!=indexOffset.get()) {
				currentIndex = indexOffset.get();
				System.out.println("   currently on estimate "+FormatHelper.nf(Math.min(currentIndex+1,nloop),9));
			}
			try {
				Thread.sleep(1000L);
			} catch(InterruptedException ie) {
				ie.printStackTrace();
				return false;
			}
		}

//		if(koption>0) {
//			try {
//				brjack.close();
//			} catch(IOException io_e) {
//				io_e.printStackTrace();
//			}
//		}
//        c
//        c Write statistics of kriged values:
//        c
		int    nk    = 0;
		double xk    = 0d;
		double vk    = 0d;
		double xkmae = 0d;
		double xkmse = 0d;
		for(Worker w: worker) {
			int wnk = w.getNumberOfKrigedPoints();
			nk += wnk;
			double[] acc = w.getAccumulators();
			xk += acc[0];
			vk += acc[1];
			xkmae += acc[2];
			xkmse += acc[3];
		}
		if(nk>0 && idbg>0) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				xk    /= nk;
				vk     = vk/nk - xk*xk;
				xkmae /= nk;
				xkmse /= nk;
				String str_stt = "Estimated   "+FormatHelper.nf(nk,8)+" blocks \n"
						+        "  average   "+FormatHelper.nf(xk,14,8)+"\n"
						+        "  variance  "+FormatHelper.nf(vk,14,8);
				bw.append(str_stt+"\n");
				System.out.println(str_stt);
				if(koption!=0) {
					bw.append("  mean error "+FormatHelper.nf(xkmae,14,8)+"  mean sqd err "+FormatHelper.nf(xkmse,14,8)+"\n");
				}
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
		
		//error-messaging
		List<String> errmsg = new ArrayList<String>();
		for(int aic=0; aic<ai.length; aic++)
			ai[aic] = 0;
		for(Worker w: worker) {
			int[] aiw = w.getAi();
			for(int aic=0; aic<ai.length; aic++)
				ai[aic] += aiw[aic];
		}
		for(int aic=0; aic<ai.length; aic++) {
			if(ai[aic]>0) {
				errmsg.add("for "+ai[aic]+" points could not do kriging correctly because got singular matrix (at entry "+aic+")");
				System.out.println("WARNING: "+errmsg.get(errmsg.size()-1));
			}
		}
		if(idbg>2) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				for(String e: errmsg)
					bw.append(e+"\n");
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
		
		//finish kriging
		return true;
	}

	private class Worker extends Thread {
		
		private boolean isFinished;
		private int nloop,delta, nk;
		private int maxEqu, maxSmp, ndb, mdt;
		private int[] nst;
		private double unbias, covmax, cbb, resc;
		private double xloc,yloc,zloc;
		private double xk, vk, xkmae, xkmse;
		private double[] x,y,z,dh,ve,vr;
		private double[] xdb,ydb,zdb, c0;
		private KdTree tree;
		
		//for error logging
		int[] ai;
		
		public Worker(int numberOfLoops, int numberOfPointsPerTask, KdTree searchTree, int maxNumberOfEquations, int maxNumberOfSamples,
				double unbias_in, double covmax_in, double cbb_in, int ndb_in, double[] c0_in, int[] nst_in, double resc_in, int mdt_in,
				double[] x_var, double[] y_var, double[] z_var, double[] dh_var, double[] ext_var, double[] vr_var,
				double[] xdb_in, double[] ydb_in, double[] zdb_in) {
			// TODO Auto-generated constructor stub
			nloop = numberOfLoops;
			delta = numberOfPointsPerTask;
			tree = searchTree;
			
			maxEqu = maxNumberOfEquations;
			maxSmp = maxNumberOfSamples;
			unbias = unbias_in;
			covmax = covmax_in;
			cbb    = cbb_in;
			ndb    = ndb_in;
			c0     = c0_in;
			nst    = nst_in;
			resc   = resc_in;
			mdt    = mdt_in;
			x  = x_var;
			y  = y_var;
			z  = z_var;
			dh = dh_var;
			ve = ext_var;
			vr = vr_var;
			xdb = xdb_in;
			ydb = ydb_in;
			zdb = zdb_in;
			
			ai = new int[maxEqu*maxEqu];
			for(int aic=0; aic<ai.length; aic++) ai[aic] = 0;
			
			isFinished = false;
		}
		
		public int getNumberOfKrigedPoints() {
			return nk;
		}
		public double[] getAccumulators() {
			return new double[] {xk, vk, xkmae, xkmse};
		}
		public int[] getAi() {
			return ai;
		}
		public boolean hasFinished() {
			return isFinished;
		}
		
		@Override
		public void run() { //TODO worker
			isFinished = false;
			
			//init accumulators:
			nk    = 0;
			xk    = 0d;
			vk    = 0d;
			xkmae = 0d;
			xkmse = 0d;
			
			//init some local variables
			int nxy = nx*ny;
			double[] xa = new double[maxSmp],
					 ya = new double[maxSmp],
					 za = new double[maxSmp];
			double[] vra = new double[maxSmp],
					 vea = new double[maxSmp];
			double[] r = new double[maxEqu],
					 rr = new double[maxEqu],
					 s= new double[maxEqu],
					 a = new double[maxEqu*maxEqu];
			boolean fircon = true;
			
			//search for work and do it
			while(true) {
				int index_position = indexOffset.addAndGet(delta)-delta;
				if(index_position>=nloop) break;

				double _true_=0d,//secj,
						extest=1d;
				double est=0d,estv=0d,resce=0d, ddh=0d, cov,cb,cb1, wt; // cmax;
				int n_accepted, ind, neq;
				int nclose = 0;
				
				int localDelta = Math.min(delta, nloop-index_position);
				for(int index=index_position; index<index_position+localDelta; index++) { //TODO begin of index-loop

//		        c Where are we making an estimate?
					int ix=0,iy=0,iz=0;
					if(koption==0) {
						iz   = (int) (index/nxy);
						iy   = (int) ((index-iz*nxy)/nx);
						ix   = index - iz*nxy - iy*nx;
						xloc = xmn + ix*xsiz;
						yloc = ymn + iy*ysiz;
						zloc = zmn + iz*zsiz;
					} else {
			            ddh    = 0d;
			            xloc   = xmn;
			            yloc   = ymn;
			            zloc   = zmn;
			            _true_ = Constants.FILL_VALUE_D;
//			            secj   = Constants.FILL_VALUE_D;
			            if(jackdf!=null) {
			            	if(jack_dh_var!=null)  ddh    = ((double[])jackdf.getArray(jack_dh_var))[index];
			            	if(jack_x_var!=null)   xloc   = ((double[])jackdf.getArray(jack_x_var))[index];
			            	if(jack_y_var!=null)   yloc   = ((double[])jackdf.getArray(jack_y_var))[index];
			            	if(jack_z_var!=null)   zloc   = ((double[])jackdf.getArray(jack_z_var))[index];
			            	if(jack_vr_var!=null)  _true_ = ((double[])jackdf.getArray(jack_vr_var))[index];
			            	if(jack_ext_var!=null) extest = ((double[])jackdf.getArray(jack_ext_var))[index];
			            }
			            if(_true_<tmin || _true_>=tmax) _true_ = Constants.FILL_VALUE_D;
					}

//		        i added by E. Metzner:
//		        i Initialise output fields, if something goes wrong on some kriging points (so est|estv = Constants.FILL_VALUE_D)
					if(iktype==Constants.INO) {
						if(koption==0) {
							writeEstimatedData(index, Constants.FILL_VALUE_D, Constants.FILL_VALUE_D, 0);
						} else {
							writeEstimatedData(index, Constants.FILL_VALUE_D, Constants.FILL_VALUE_D, 0, xloc, yloc, zloc, Constants.FILL_VALUE_D, Constants.FILL_VALUE_D);
						}
					} else {
						writeEstimatedData(index, s, vra, 0, Constants.FILL_VALUE_D);
					}

//		        c Read in the external drift variable for this grid node if needed:
					if(ktype==SIMPLE_NON_STATIONARY_KRIGING || ktype==EXTERNAL_DRIFT_KRIGING) {
						if(koption==0) {
							extest = ((double[])externdf.getArray(ext_var_e))[index];
//				                  read(lext,*) (var(i),i=1,iextve)
//				                  extest = var(iextve)
						}
						if(extest<tmin || extest>=tmax) {
							est  = Constants.FILL_VALUE_D;
							estv = Constants.FILL_VALUE_D;
							continue;
						}
						resce  = covmax / Math.max(extest,0.0001d);
					}

//		        c Find the nearest samples:
					//int[] supres = DataHelper.srchsupr(xloc, yloc, zloc, radsqd, isrot, rotmat,
					//		nsbtosr, ixsbtosr, iysbtosr, izsbtosr, noct,
					//		nd, x, y, z, tmp, nisb, superblock_grid, closest);
					n_accepted = 0;
					int actualUsedOctants = 0;
					int tryNumOctants = numOctants;
					while(tryNumOctants>0) {
						actualUsedOctants = tryNumOctants;
						int[] supres = tree.getIndexOfNClosestPointsTo(ndmin, ndmax, numOctants, rotmat[MAXNST], false, false, xloc, yloc, zloc); //use redefined maximum search distance
						nclose = supres[0]; //TODO
						//int infoct = supres[1];
	
	//		        c Load the nearest data in xa,ya,za,vra,vea:
						n_accepted = 0;
	//					System.out.println("[DEBUG]  pivot:                 pos "+
	//							FormatHelper.nf(xloc,12,8)+" "+FormatHelper.nf(yloc,12,8)+" "+FormatHelper.nf(zloc,12,8));
						double pre_x=0d, pre_y=0d, pre_z=0d, pre_r=0d, pre_e=0d, pre_cov;
						for(int i=0; i<nclose; i++) {
							ind    = supres[i+1];
	//						System.out.println("[DEBUG]  closest point: ind "+FormatHelper.nf(ind,3)+" pos "+
	//								FormatHelper.nf(x[ind],12,8)+" "+FormatHelper.nf(y[ind],12,8)+" "+FormatHelper.nf(z[ind],12,8));
							boolean accept = true;
							if(koption!=0 && (Math.abs(x[ind]-xloc)+Math.abs(y[ind]-yloc)+Math.abs(z[ind]-zloc))<EPSLON)
								accept = false;
							if(koption!=0 && (Math.abs(dh[ind]-ddh))<EPSLON)
								accept = false;
							if(accept) { //pre-acception
								pre_x = x[ind] - xloc + 0.5d*xsiz;
								pre_y = y[ind] - yloc + 0.5d*ysiz;
								pre_z = z[ind] - zloc + 0.5d*zsiz;
								pre_r = vr[ind];
								pre_e = ve[ind];
								if(ndb<=1) {
									pre_cov = Covariance.cova3(pre_x, pre_y, pre_z, xdb[0], ydb[0], zdb[0],
											1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
								} else {
									pre_cov  = 0d;
									for(int j=0; j<ndb; j++) {
										cov = Covariance.cova3(pre_x, pre_y, pre_z, xdb[j], ydb[j], zdb[j],
												1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
										pre_cov += cov;
										double dx = pre_x - xdb[j];
										double dy = pre_y - ydb[j];
										double dz = pre_z - zdb[j];
										if(dx*dx+dy*dy+dz*dz < EPSLON) pre_cov -= c0[0];
									}
									pre_cov /= ndb;
								}
								accept = Math.abs(pre_cov)>EPSLON;
							}
							if(accept) { //pre-acception
								if(n_accepted<ndmax) {
									xa[n_accepted]  = pre_x;
									ya[n_accepted]  = pre_y;
									za[n_accepted]  = pre_z;
									vra[n_accepted] = pre_r;
									vea[n_accepted] = pre_e;
									n_accepted++; //copy afterward for JAVA indices
								}
							}
						}
						if(n_accepted>=ndmin)
							break;
						tryNumOctants = Math.min(Math.max(1+(n_accepted-1)/ndmin, 0), tryNumOctants-1);
					}

//		        c Test number of samples found:
					int ndmino = ndmin*actualUsedOctants;
					if(n_accepted<ndmino) {
						est  = Constants.FILL_VALUE_D;
						estv = Constants.FILL_VALUE_D;
						if(idbg>=2) {
							try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
								bw.append(" Encountered a location where there were too few data\n"
										+ " for Ord.Kriging or Simple Kriging. KT3D currently\n"
										+ " leaves these locations Constants.FILL_VALUEimated.\n");
								bw.flush();
							}catch(IOException io_e) {
								io_e.printStackTrace();
							}
							System.out.println("   Too few data: No Ordinary or Simple Kriging!");
						}
						continue;
					}

//		        c Test if there are enough samples to estimate all drift terms:
					if(n_accepted>=1 && n_accepted<=mdt) {
						if(fircon) {
							try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
								bw.append(" Encountered a location where there were too few data\n"
										+ " to estimate all of the drift terms but there would be\n"
										+ " enough data for Ord.Kriging or Simple Kriging. KT3D\n"
										+ " currently leaves these locations Constants.FILL_VALUEimated.\n"
										+ " This message is only written once - the first time.\n");
								bw.flush();
							}catch(IOException io_e) {
								io_e.printStackTrace();
							}
							fircon = false;
						}
						est  = Constants.FILL_VALUE_D;
						estv = Constants.FILL_VALUE_D;
						continue;
					}

//		        c There are enough samples - proceed with estimation.
					if(n_accepted<=1) {

//		        c Handle the situation of only one sample:
						cb1 = Covariance.cova3(xa[0], ya[0], za[0], xa[0], ya[0], za[0],
								1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];

//		        c Establish Right Hand Side Covariance:
						if(ndb<=1) {
							cb = Covariance.cova3(xa[0], ya[0], za[0], xdb[0], ydb[0], zdb[0],
									1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
						} else {
							cb  = 0d;
							for(int i=0; i<ndb; i++) {
								cov = Covariance.cova3(xa[0], ya[0], za[0], xdb[i], ydb[i], zdb[i],
										1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
								cb += cov;
								double dx = xa[0] - xdb[i];
								double dy = ya[0] - ydb[i];
								double dz = za[0] - zdb[i];
								if(dx*dx+dy*dy+dz*dz < EPSLON) cb -= c0[0];
							}
							cb /= ndb;
						}

//		        c Early bug - always did OK in presence of one data.
						if(ktype==SIMPLE_NON_STATIONARY_KRIGING) skmean = extest;
						if(ktype==SIMPLE_KRIGING || ktype==SIMPLE_NON_STATIONARY_KRIGING) {
							wt   = cb / cb1;
							est  = wt * vra[0] + (1d-wt) * skmean;
							estv = cbb - wt*cb;
						} else {
							est  = vra[0];
							estv = cbb - 2d*cb + cb1;
						}
						nk++;
						xk += est;
						vk += est*est;
						continue;
					}

//		        c Go ahead and set up the OK portion of the kriging matrix:
					neq = mdt+n_accepted;

//		        c Initialize the main kriging matrix:
					for(int i=0; i<neq*neq; i++) {
						a[i] = 0d;
					}

//		        c Fill in the kriging matrix:
					for(int i=0; i<n_accepted; i++) for(int j=i; j<n_accepted; j++) {
						cov = Covariance.cova3(xa[i], ya[i], za[i], xa[j], ya[j], za[j],
								1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
						a[neq*i+j] = cov;
						a[neq*j+i] = cov;
					}

//		        c Fill in the OK unbiasedness portion of the matrix (if not doing SK):
					if(neq>n_accepted) {
						for(int i=0; i<n_accepted; i++) {
							a[neq*i+n_accepted] = unbias;
							a[neq*n_accepted+i] = unbias;
						}
					}

//		        c Set up the right hand side:
					//System.out.println("[KT3D] point to estimate:"); //TODO DEBUG, remove later
					//System.out.println("         "+FormatHelper.nf(xdb[0],10,6)+"  "+FormatHelper.nf(ydb[0],10,6)+"  "+FormatHelper.nf(zdb[0],10,6));
					//System.out.println("       points surround with value:");
					for(int i=0; i<n_accepted; i++) {
						if(ndb<=1) {
							cb = Covariance.cova3(xa[i],ya[i],za[i], xdb[0],ydb[0],zdb[0],
									1,nst,MAXNST, c0,it, withSAHE,cc,aa, 1,MAXROT,rotmat)[0];
							//System.out.println("         xyz={"+FormatHelper.nf(xa[i],10,6)+" "+FormatHelper.nf(ya[i],10,6)+" "+FormatHelper.nf(za[i],10,6)+"}  "+
							//		"abc={"+FormatHelper.nf(xdb[0],10,6)+" "+FormatHelper.nf(ydb[0],10,6)+" "+FormatHelper.nf(zdb[0],10,6)+"}  "+
							//		FormatHelper.nf(cb,10,6));
						} else {
							cb  = 0d;
							//String tempstr = "";
							for(int j=0; j<ndb; j++) {
								cov = Covariance.cova3(xa[i], ya[i], za[i], xdb[j], ydb[j], zdb[j],
										1, nst, MAXNST, c0, it, withSAHE, cc, aa, 1, MAXROT, rotmat)[0];
								cb += cov;
								//tempstr += " "+FormatHelper.nf(cov,8,6);
								double dx = xa[i] - xdb[j];
								double dy = ya[i] - ydb[j];
								double dz = za[i] - zdb[j];
								if(dx*dx+dy*dy+dz*dz < EPSLON) cb -= c0[0];
							}
							cb /= ndb;
							//System.out.println("           ,--<  "+tempstr);
							//System.out.println("         "+FormatHelper.nf(xa[i],10,6)+"  "+FormatHelper.nf(ya[i],10,6)+"  "+FormatHelper.nf(za[i],10,6)+"  "+
							//		FormatHelper.nf(cb,10,6));
						}
						r[i] = cb;
					}
					if(neq>n_accepted) r[n_accepted] = unbias;

//		        c Add the additional unbiasedness constraints:
					int im = n_accepted + 1;

//		        c (dt) drift term (linear|quadratic|cubic in "x|y"):
					for(int dt=0; dt<MAXDT; dt++) {
						if(idrif[dt]==Constants.IYES) {
							for(int k=0; k<n_accepted; k++) {
								a[neq*im+k] = drift_term(dt+1, xa[k], ya[k], za[k]) * resc;
								a[neq*k+im] = a[neq*im+k];
							}
							r[im] = bv[dt];
							im++;
						}
					}
//		        c
//		        c External drift term (specified by external variable):
//		        c
					if(ktype==EXTERNAL_DRIFT_KRIGING) {
						//im=im+1
						for(int k=0; k<n_accepted; k++) {
							a[neq*im+k] = vea[k]*resce;
							a[neq*k+im] = vea[k]*resce;
						}
						r[im] = extest*resce;
						im++; // copied afterwards for JAVA indices
					}
//		        c
//		        c Copy the right hand side to compute the kriging variance later:
//		        c
					for(int k=0; k<neq; k++) {
						rr[k] = r[k];
					}
					//int kadim = neq * neq;
					//int ksdim = neq;
					int nrhs  = 1;
					int nv    = 1;
//		        c
//		        c If estimating the trend then reset all the right hand side terms=0.0:
//		        c
					if(itrend==Constants.IYES) {
						for(int i=0; i<n_accepted; i++) {
							r[i]  = 0d;
							rr[i] = 0d;
						}
					}
//		        c
//		        c Write out the kriging Matrix if Seriously Debugging:
//		        c
					if(idbg==3) {
						try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
							bw.append("\nEstimating node index : "+(ix+1)+" "+(iy+1)+" "+(iz+1)+"\n");
							System.out.println("Estimating node index : "+(ix+1)+" "+(iy+1)+" "+(iz+1)+"\n");
				            int is = 1 - neq,ie;
				            if(neq==0) {
				            	bw.append("    number of equation is zero, no solving possible!");
				            	System.out.println("    number of equation is zero, no solving possible!");
				            }
				            for(int i=0; i<neq; i++) {
				            	is = i*neq;
				            	ie = is + neq;
				            	String str_a = "";
				            	for(int j=is; j<ie; j++) str_a += " "+FormatHelper.nf(a[j],7,4);
				            	bw.append("    r("+FormatHelper.nf(i+1,2)+")= "+FormatHelper.nf(r[i],7,4)+"  a="+str_a+"\n");
				            	System.out.println("    r("+FormatHelper.nf(i+1,2)+")= "+FormatHelper.nf(r[i],7,4)+"  a="+str_a);
				            }
				            bw.flush();
						}catch(IOException io_e) {
							io_e.printStackTrace();
						}
					}
//		        c
//		        c Solve the kriging system:
//		        c
				//      call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
					s = LinEquSolver.ktsol(neq, nrhs, nv, a, r, maxEqu, ai); //TODO ktsol
//		        c
//		        c Compute the solution:
//		        c
					if(s==null) {
						if(idbg>=3) {
							try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
								bw.append(" Singular Matrix "+(ix+1)+" "+(iy+1)+" "+(iz+1)+"\n");
								bw.flush();
							}catch(IOException io_e) {
								io_e.printStackTrace();
							}
						}
						est  = Constants.FILL_VALUE_D;
						estv = Constants.FILL_VALUE_D;
					} else {
						est  = 0d;
						estv = cbb;
						if(ktype==SIMPLE_NON_STATIONARY_KRIGING) skmean = extest;
						for(int j=0; j<neq; j++) {
							estv -= s[j]*rr[j];
							if(j<n_accepted) {
								if(ktype==SIMPLE_KRIGING) {
									est += s[j]*(vra[j]-skmean);
								} else if(ktype==SIMPLE_NON_STATIONARY_KRIGING) {
									est += s[j]*(vra[j]-vea[j]);
								} else {
									est += s[j]*vra[j];
								}
							}
						}
						//TODO is following line a right estimate for variance if not desired octant number is reached?
						estv *= numOctants / actualUsedOctants;
						if(ktype==SIMPLE_KRIGING || ktype==SIMPLE_NON_STATIONARY_KRIGING) est += skmean;
						nk++;
						xk += est;
						vk += est*est;
//		        c
//		        c Write the kriging weights and data if debugging level is above 2:
//		        c
						if(idbg>=2) {
							try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
								bw.append("\nBLOCK: "+(ix+1)+" "+(iy+1)+" "+(iz+1)+" at "+xloc+" "+yloc+" "+zloc+"\n\n");
								if(ktype!=SIMPLE_KRIGING)
									bw.append("  Lagrange : "+(s[n_accepted]*unbias)+"\n");
								bw.append("  BLOCK EST: x,y,z,vr,wt \n");
								for(int i=0; i<n_accepted; i++) {
									xa[i] += xloc - 0.5d*xsiz;
									ya[i] += yloc - 0.5d*ysiz;
									za[i] += zloc - 0.5d*zsiz;
									bw.append(" "+FormatHelper.nf(xa[i],12,3)+" "+FormatHelper.nf(ya[i],12,3)+" "+FormatHelper.nf(za[i],12,3)
											+ " "+FormatHelper.nf(vra[i],12,3)+" "+FormatHelper.nf(s[i],12,3)+"\n");
								}
								bw.append("  estimate, variance  "+est+" "+estv+"\n");
								bw.flush();
							} catch(IOException io_e) {
								io_e.printStackTrace();
							}
						}
					}
//		        c
//		        c END OF MAIN KRIGING LOOP:
//		        c
				//} // Loop Mark 1: continue;
					if(iktype==Constants.INO) { //index-kriging?
						if(koption==0) {
							writeEstimatedData(index, est, estv, n_accepted);
						} else {
							double err = Constants.FILL_VALUE_D;
							if(_true_!=Constants.FILL_VALUE_D && est!=Constants.FILL_VALUE_D) {
								err=est-_true_;
								xkmae += Math.abs(err);
								xkmse += err*err;
							}
							writeEstimatedData(index, est, estv, n_accepted, xloc, yloc, zloc, _true_, err);
						}
					} else {
						writeEstimatedData(index, s, vra, n_accepted, _true_);
					}
//				      end do
				}
			}
			
			isFinished = true;
		}
	}






	private double drift_term(int term_id, double _x, double _y, double _z) {
		switch(term_id) {
			case 1: return _x; // linear in x
			case 2: return _y; // linear in y
			case 3: return _z; // linear in z
			case 4: return _x*_x; // quadratic in x
			case 5: return _y*_y; // quadratic in x
			case 6: return _z*_z; // quadratic in x
			case 7: return _x*_y; // cross square in x and y
			case 8: return _x*_z; // cross square in x and z
			case 9: return _y*_z; // cross square in y and z
			default: return 0d;
		}
	}
	private void writeEstimatedData(int _index, double _est, double _estv, int _n_accepted) {
		estimates[0][_index] = _est;
		estimates[1][_index] = _estv;
		estimates[2][_index] = _n_accepted;
	}
	private void writeEstimatedData(int _index, double _est, double _estv, int _n_accepted,
			double _xloc, double _yloc, double _zloc, double __true_, double _err) {
		estimates[0][_index] = _xloc;
		estimates[1][_index] = _yloc;
		estimates[2][_index] = _zloc;
		estimates[3][_index] = __true_;
		estimates[4][_index] = _est;
		estimates[5][_index] = _estv;
		estimates[6][_index] = _err;
		estimates[7][_index] = _n_accepted;
	}
	private void writeEstimatedData(int _index, double[] _s, double[] _vra, int _n_accepted, double __true_) {
//      c
//      c Work out the IK-type distribution implicit to this data configuration
//      c and kriging weights:
//      c
		for(int icut=0; icut<ncut; icut++) {
			cdf[icut] = -1d;
		}
		double wtmin = 1d;
		for(int i=0; i<_n_accepted; i++) {
			if(_s[i]<wtmin) wtmin = _s[i];
		}
		double sumwt = 0d;
		for(int i=0; i<_n_accepted; i++) {
			_s[i]  -= wtmin;
			sumwt += _s[i];
		}
		for(int i=0; i<_n_accepted; i++) {
			_s[i] /= Math.max(0.00001d,sumwt);
		}
		if(_n_accepted>1 && sumwt<0.00001d) {
			for(int icut=0; icut<ncut; icut++) {
				cdf[icut] = 0d;
				for(int i=0; i<_n_accepted; i++) {
					if(_vra[i]<=cut[icut]) cdf[icut] += _s[i];
				}
			}
		}
		for(int i=0; i<ncut; i++) estimates[i][_index] = cdf[i];
		estimates[ncut][_index] = __true_;
	}

	public void writeGslibOutputfile(String output_file_path) {
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_file_path),true))) {
			int nloop = (koption==0 ? nx*ny*nz : 10000000);
			for(int index=0; index<nloop; index++) {
				if(iktype==0) {
						if(koption==0) {
							bw.append(" "+FormatHelper.nf(estimates[0][index],14,8)+" "+FormatHelper.nf(estimates[1][index],14,8)+"\n");
						} else {
							bw.append(" "+FormatHelper.nf(estimates[0][index],14,8)+" "+FormatHelper.nf(estimates[1][index],14,8)+" "+FormatHelper.nf(estimates[2][index],14,8)
									+ " "+FormatHelper.nf(estimates[3][index],14,8)+" "+FormatHelper.nf(estimates[4][index],14,8)+" "+FormatHelper.nf(estimates[5][index],14,8)
									+ " "+FormatHelper.nf(estimates[6][index],14,8)+"\n");
						}
				} else {
						String str_cdf = "";
						for(int i=0; i<ncut; i++) str_cdf+=" "+FormatHelper.nf(estimates[i][index],8,4);
						if(koption==0) {
							bw.append(str_cdf+"\n");
						} else {
							bw.append(str_cdf+" "+FormatHelper.nf(estimates[ncut][index],8,4)+"\n");
						}
				}
			}
			bw.flush();
		}catch(IOException io_e) {
			io_e.printStackTrace();
		}
	}
	
	public void makepar() {
		// c-----------------------------------------------------------------------
		// c
		// c Write a Parameter File
		// c **********************
		// c
		// c
		// c
		// c-----------------------------------------------------------------------
		// lun = 99
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File("res/kt3d.par")))) {
			bw.append("                  Parameters for KT3D\n" + "                  *******************\n"
					+ "START OF PARAMETERS:\n" + "res/cluster.dat                  -file with data\n"
					+ "0  1  2  0  3  0                 -columns for DH,X,Y,Z,var,sec var\n"
					+ "-1.0e21   1.0e21                 -trimming limits\n"
					+ "0                                -option: 0=grid, 1=cross, 2=jackknife\n"
					+ "res/xvk.dat                      -file with jackknife data\n"
					+ "1   2   0    3    0              -columns for X,Y,Z,vr and sec var\n"
					+ "3                                -debugging level: 0,1,2,3\n"
					+ "res/kt3d.dbg                     -file for debugging output\n"
					+ "res/kt3d.out                     -file for kriged output\n"
					+ "50   0.5    1.0                  -nx,xmn,xsiz\n"
					+ "50   0.5    1.0                  -ny,ymn,ysiz\n"
					+ "1    0.5    1.0                  -nz,zmn,zsiz\n"
					+ "5    5      5                    -x,y and z block discretization\n"
					+ "4    8                           -min, max data for kriging\n"
					+ "0                                -max per octant (0->not used)\n"
					+ "20.0  20.0  20.0                 -maximum search radii\n"
					+ " 0.0   0.0   0.0                 -angles for search ellipsoid\n"
					+ "1     2.302                      -0=SK,1=OK,2=non-st SK,3=exdrift\n"
					+ "0 0 0 0 0 0 0 0 0                -drift: x,y,z,xx,yy,zz,xy,xz,zy\n"
					+ "0                                -0, variable; 1, estimate trend\n"
					+ "res/extdrift.dat                 -gridded file with drift/mean\n"
					+ "4                                -column number in gridded file\n"
					+ "1    0.2                         -nst, nugget effect\n"
					+ "1    0.8  0.0   0.0   0.0        -it,cc,ang1,ang2,ang3\n"
					+ "10.0  10.0  10.0                 -a_hmax, a_hmin, a_vert\n");
			bw.flush();
		} catch (IOException io_e) {
			io_e.printStackTrace();
		}
	}
}
