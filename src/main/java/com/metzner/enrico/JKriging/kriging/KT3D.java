package com.metzner.enrico.JKriging.kriging;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.helper.DataHelper;
import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.LinEquSolver;
import com.metzner.enrico.JKriging.helper.Rotation;
import com.metzner.enrico.JKriging.probability.Covariance;

public class KT3D {

	public static final double VERSION = 3.000d;
	public static final double EPSLON = 0.000001d;
	private double UNEST = 1.0e20d;
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

	private int iktype, ncut, koption, idbg;
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
	private double xmn, xsiz, ymn, ysiz, zmn, zsiz, radius, sanis1, sanis2, sang1, sang2, sang3, skmean;

	private int MAXCUT; // MXSXY,MXSX;

	int idhlj, ixlj, iylj, izlj, ivrlj, iextvj, iextve, nvarij;

	private DataFrame dataframe, jackdf, externdf;
	private String x_var, y_var, z_var, vr_var, dh_var, ext_var, ext_var_e;
	private String jack_x_var, jack_y_var, jack_z_var, jack_vr_var, jack_dh_var, jack_ext_var;
	private List<double[]> variograms;
	int num_krig_res; // number of kriging output variables (2 normal, 7 indicator)
	String[] estimate_titles; // titles for results from kriging
	private double[][] estimates; // result from kriging

	private boolean[] paramchecklist = new boolean[8];
	String[] param_descriptions = { "dataframe", "variable-cooredinates and kriging variable", "kriging field/points",
			"kriging search parameters", "debugging setting", "usage of jackknife file",
			"external drift file and parameters", "variogram model(s)" };

	public KT3D() {
		variograms = new ArrayList<double[]>();
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
			iktype = (use_indicator_kriging ? 1 : 0);
			if(iktype==0) {
				if(kriging_option==0) {
					paramchecklist[5] = true;
					paramchecklist[6] = true;
					num_krig_res = 2;
					estimate_titles = new String[] {"Estimate", "EstimateVariance"};
				} else {
					num_krig_res = 7;
					estimate_titles = new String[] {"X","Y","Z","True","Estimate","EstimateVariance","Error: est-true"};
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
			int max_per_octant, double search_radius, double second_radius, double third_radius,
			double anisotropy_angle, double second_angle, double third_angle) {
		nxdis = x_discretisation; nydis = y_discretisation; nzdis = z_discretisation;
		ndmin = min_points; ndmax = max_points; noct = max_per_octant;
		radius = search_radius; sanis1 = second_radius / search_radius; sanis2 = third_radius / search_radius;
		sang1 = anisotropy_angle; sang2 = second_angle; sang3 = third_angle;
		paramchecklist[3] = true;
		return this;
	}
	public KT3D setDebuggingLevel(int debugging_level, String debugging_file_path) {
		idbg  = debugging_level;
		dbgfl = debugging_file_path;
		paramchecklist[4] = true;
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
			ext_var_e = externdf.getVarname(_ext_id-1);
			paramchecklist[6] = true;
		}
		return this;
	}
	public KT3D setExternalDrift(DataFrame extern_dataframe, String extern_variable_name) {
		if(extern_dataframe==null) {
			System.err.println("Dataframe with external drift data does not exist!");
		} else {
			setExternalDrift(extern_dataframe, extern_dataframe.getVariableID(extern_variable_name)+1);
		}
		return this;
	}
	public KT3D addVariogramModel(int variogram_type, double covariance,
			double azimuth, double dip, double roll, double h_max, double h_min, double h_vert) {
		if(h_max<0d) {
			System.err.println("search radius must be greater than zero!");
		} else if (variogram_type == 4 && (h_max < 0d || h_min > 2d)) {
			System.err.println("INVALID power variogram");
		} else {
			double anis1 = h_min / Math.max(h_max, EPSLON);
			double anis2 = h_vert / Math.max(h_max, EPSLON);
			System.out.println("search anisotropy angles = " + azimuth + " " + dip + " " + roll);
			System.out.println(" a1 a2 a3 =  " + h_max + " " + h_min + " " + h_vert);
			variograms.add(new double[] {variogram_type+0.1d, covariance, azimuth, dip, roll, h_max, anis1, anis2});
			paramchecklist[7] = true;
		}
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
			iktype = 0;
			if (koption < 0) {
				iktype = 1;
				koption = -koption;
			}
			if (iktype == 1) {
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
			sanis1 = radius1 / radius;
			sanis2 = radius2 / radius;

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			sang1 = Double.parseDouble(parts[0].trim());
			sang2 = Double.parseDouble(parts[1].trim());
			sang3 = Double.parseDouble(parts[2].trim());
			System.out.println(" search anisotropy angles = " + sang1 + " " + sang2 + " " + sang3);
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
		num_krig_res = (iktype==0 ? (koption==0 ? 2 : 7) : ncut);
		if(iktype==0) { estimate_titles = (koption==0 ? new String[] {"Estimate", "EstimateVariance"} :
			new String[] {"X","Y","Z","True","Estimate","EstimateVariance","Error: est-true"});
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

	public DataFrame kt3d() {
		return (DataFrame) kt3d_df(1);
	}
	public DataFrame3D kt3d_asDataFrame3D() {
		return (DataFrame3D) kt3d_df(3);
	}
	
	private Object kt3d_df(int _dataframe_dims) { //TODO beginning of KT3D
		for(int pcl=0; pcl<8; pcl++)
			if(!paramchecklist[pcl]) {
				System.err.println("Not all parameters are set! \""+param_descriptions[pcl]+"\" missing");
				return null;
			}
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
		double[] tmp = new double[datalength];
		double[] closest = new double[datalength];
		if (x_var==null && nx > 1)
			System.out.println(" WARNING: no x variable and nx>1 !");
		if (y_var==null && ny > 1)
			System.out.println(" WARNING: no y variable and ny>1 !");
		if (z_var==null && nz > 1)
			System.out.println(" WARNING: no z variable and nz>1 !");
		
		int MAXDIS = nxdis*nydis*nzdis;
		int maxSamples = ndmax + 1;
		int maxEquations = maxSamples + MAXDT + 2;
		int maxSupBlckX = Math.max(1, Math.min(50, nx / 2));
		int maxSupBlckY = Math.max(1, Math.min(50, ny / 2));
		int maxSupBlckZ = Math.max(1, Math.min(50, nz / 2));
		int maxSuperblocks = maxSupBlckX*maxSupBlckY*maxSupBlckZ;
		//MXSXY = 4 * MAXSBX * MAXSBY;
		//MXSX  = 2 * MAXSBX;
//		if (ndmax > maxSamples) {
//			System.err.println("ndmax is too big - modify .inc file");
//			return null;
//		}
		int[]    nisb = new int[maxSuperblocks],
				 ixsbtosr = new int[8*maxSuperblocks],
				 iysbtosr = new int[8*maxSuperblocks],
				 izsbtosr = new int[8*maxSuperblocks];
		double[] xa = new double[maxSamples],
				 ya = new double[maxSamples],
				 za = new double[maxSamples];
		double[] vra = new double[maxSamples],
				 vea = new double[maxSamples];
		double[] xdb = new double[MAXDIS],
				 ydb = new double[MAXDIS],
				 zdb = new double[MAXDIS];
		double[] r = new double[maxEquations],
				 rr = new double[maxEquations],
				 s= new double[maxEquations],
				 a = new double[maxEquations*maxEquations];
		
		int[] nst = {variograms.size(),0};
		double[] c0 = {1d,0d};
		for(int i=0; i<nst[0]; i++) {
			double[] vario = variograms.get(i);
			it[i]    = (int) vario[0];
			cc[i]    = vario[1]; c0[0] -= cc[i];
			ang1[i]  = vario[2];
			ang2[i]  = vario[3];
			ang3[i]  = vario[4];
			aa[i]    = vario[5];
			anis1[i] = vario[6];
			anis2[i] = vario[7];
		}
		
		
		
		

		double cbb;
		boolean fircon=true, accept;
		double[] sec3 = new double[datalength];
//        c
//        c Set up the rotation/anisotropy matrices that are needed for the
//        c variogram and search.  Also compute the maximum covariance for
//        c the rescaling factor:
//        c
		System.out.println("Setting up rotation matrices for variogram and search");
		double radsqd = radius * radius;
		double PMX    = 999.0;
		double covmax = c0[0];
		for(int is=0; is<nst[0]; is++) {
			rotmat = Rotation.setrot(ang1[is],ang2[is],ang3[is],anis1[is],anis2[is],is+1,MAXROT,rotmat);
			if(it[is]==4) {
				covmax += PMX; 
			} else {
				covmax += cc[is];
			}
		}
		int isrot = MAXNST + 1;
		rotmat = Rotation.setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat);
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
			return null;
		}
		resc = 1d / resc;
		System.out.println("[DEBUG] Maximum covariance: "+covmax);
//        c
//        c Set up for super block searching:
//        c
		System.out.println("Setting up super block search strategy");
		System.out.println("[DEBUG]");
		FormatHelper.printTable(20, x,y,z,vr,tmp);
		int nsec = 2;
		double[] superblock_grid = new double[9]; // replace n[xyz]sup, [xyz]mnsup, [xyz]sizsup
		nisb = DataHelper.setsupr(nx, xmn, xsiz, ny, ymn, ysiz, nz, zmn, zsiz, x, y, z, vr, tmp,
				nsec, ve, dh, sec3, maxSupBlckX, maxSupBlckY, maxSupBlckZ, nisb, superblock_grid);
		//      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
		//     +             vr,tmp,nsec,ve,dh,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb,
		//     +             nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
		//     +             zmnsup,zsizsup)
		System.out.println("[DEBUG]");
		FormatHelper.printTable(20, x,y,z,vr,tmp);
		int nsbtosr = DataHelper.picksupr(superblock_grid, isrot, rotmat, radsqd, ixsbtosr, iysbtosr, izsbtosr);
		//      call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
		//     +             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
		//     +             iysbtosr,izsbtosr)
//        c
//        c Compute the number of drift terms, if an external drift is being
//        c considered then it is one more drift term, if SK is being considered
//        c then we will set all the drift terms off and mdt to 0):
//        c
		System.out.println("[DEBUG] super block strategy set");
		int mdt = 1;
		for(int i=0; i<9; i++) {
			if(ktype==0 || ktype==2) idrif[i] = 0;
			if(idrif[i]<0 || idrif[i]>1) {
				System.err.println("ERROR KT3D: invalid drift term "+idrif[i]);
				return null;
			}
			mdt += idrif[i];
		}
		if(ktype==3) mdt++;
		if(ktype==0) mdt = 0;
		if(ktype==2) mdt = 0;
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
			return null;
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
		int    nk    = 0;
		double xk    = 0d;
		double vk    = 0d;
		double xkmae = 0d;
		double xkmse = 0d;

//        c Calculate Block Covariance. Check for point kriging.
		double cov = Covariance.cova3(xdb[0], ydb[0], zdb[0], xdb[0], ydb[0], zdb[0],
				1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];

//        c Set the ``unbias'' variable so that the matrix solution is more stable
		double unbias = cov;
		cbb    = cov;
		if(ndb>1) {
			cbb = 0;
			for(int i=0; i<ndb; i++) {
				for(int j=0; j<ndb; j++) {
					cov = Covariance.cova3(xdb[i], ydb[i], zdb[i], xdb[j], ydb[j], zdb[j],
							1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];
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
		int nxy=1,nxyz=1,nloop,irepo;
		int nd=vr.length;
		if(koption==0) {
			nxy   = nx*ny;
			nxyz  = nx*ny*nz;
			nloop = nxyz;
			irepo = Math.max(1,Math.min(10000,nxyz/10));
		} else {
			nloop = 10000000;
			irepo = Math.max(1,Math.min(10000,nd/10));
		}
		estimates = new double[num_krig_res][nloop];
		double ddh = 0d;
		System.out.println("\n\nWorking on the kriging\n");

//        c MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
		double _true_=0d,//secj,
				extest=1d;
		double est=0d,estv=0d,resce=0d, cb,cb1, wt; // cmax;
		int na, ind, neq;
		int nclose = closest.length;
		for(int index=0; index<nloop; index++) {
			if((index+1)%irepo==0) System.out.println("   currently on estimate "+FormatHelper.nf(index+1,9));

//        c Where are we making an estimate?
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
	            _true_ = UNEST;
//	            secj   = UNEST;
	            if(jackdf!=null) {
	            	if(jack_dh_var!=null)  ddh    = ((double[])jackdf.getArray(jack_dh_var))[index];
	            	if(jack_x_var!=null)   xloc   = ((double[])jackdf.getArray(jack_x_var))[index];
	            	if(jack_y_var!=null)   yloc   = ((double[])jackdf.getArray(jack_y_var))[index];
	            	if(jack_z_var!=null)   zloc   = ((double[])jackdf.getArray(jack_z_var))[index];
	            	if(jack_vr_var!=null)  _true_ = ((double[])jackdf.getArray(jack_vr_var))[index];
	            	if(jack_ext_var!=null) extest = ((double[])jackdf.getArray(jack_ext_var))[index];
	            }
	            if(_true_<tmin || _true_>=tmax) _true_ = UNEST;
			}

//        c Read in the external drift variable for this grid node if needed:
			if(ktype==2 || ktype==3) {
				if(koption==0) {
					extest = ((double[])externdf.getArray(ext_var_e))[index];
//		                  read(lext,*) (var(i),i=1,iextve)
//		                  extest = var(iextve)
				}
				if(extest<tmin || extest>=tmax) {
					est  = UNEST;
					estv = UNEST;
					continue;
				}
				resce  = covmax / Math.max(extest,0.0001d);
			}

//        c Find the nearest samples:
			int[] supres = DataHelper.srchsupr(xloc, yloc, zloc, radsqd, isrot, rotmat,
					nsbtosr, ixsbtosr, iysbtosr, izsbtosr, noct,
					nd, x, y, z, tmp, nisb, superblock_grid, closest);
			nclose = supres[0];
			//int infoct = supres[1];

//        c Load the nearest data in xa,ya,za,vra,vea:
			na = 0;
			for(int i=0; i<nclose; i++) {
				ind    = (int) (closest[i]+0.5d);
				accept = true;
				if(koption!=0 && (Math.abs(x[ind]-xloc)+Math.abs(y[ind]-yloc)+ Math.abs(z[ind]-zloc))<EPSLON)
					accept = false;
				if(koption!=0 && (Math.abs(dh[ind]-ddh))<EPSLON)
					accept = false;
				if(accept) {
					if(na<ndmax) {
						//na = na + 1
						xa[na]  = x[ind] - xloc + 0.5d*xsiz;
						ya[na]  = y[ind] - yloc + 0.5d*ysiz;
						za[na]  = z[ind] - zloc + 0.5d*zsiz;
						vra[na] = vr[ind];
						vea[na] = ve[ind];
						na++; //copy afterward for JAVA indices
					}
				}
			}

//        c Test number of samples found:
			if(na<ndmin) {
				est  = UNEST;
				estv = UNEST;
				continue;
			}

//        c Test if there are enough samples to estimate all drift terms:
			if(na>=1 && na<=mdt) {
				if(fircon) {
					try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
						bw.append(" Encountered a location where there were too few data\n"
								+ " to estimate all of the drift terms but there would be\n"
								+ " enough data for Ord.Kriging or Simple Kriging. KT3D\n"
								+ " currently leaves these locations unestimated.\n"
								+ " This message is only written once - the first time.\n");
						bw.flush();
					}catch(IOException io_e) {
						io_e.printStackTrace();
					}
					fircon = false;
				}
				est  = UNEST;
				estv = UNEST;
				continue;
			}

//        c There are enough samples - proceed with estimation.
			if(na<=1) {

//        c Handle the situation of only one sample:
				cb1 = Covariance.cova3(xa[0], ya[0], za[0], xa[0], ya[0], za[0],
						1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];

//        c Establish Right Hand Side Covariance:
				if(ndb<=1) {
					cb = Covariance.cova3(xa[0], ya[0], za[0], xdb[0], ydb[0], zdb[0],
							1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];
				} else {
					cb  = 0d;
					for(int i=0; i<ndb; i++) {
						cov = Covariance.cova3(xa[0], ya[0], za[0], xdb[i], ydb[i], zdb[0],
								1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];
						cb += cov;
						double dx = xa[0] - xdb[i];
						double dy = ya[0] - ydb[i];
						double dz = za[0] - zdb[i];
						if(dx*dx+dy*dy+dz*dz < EPSLON) cb -= c0[0];
					}
					cb /= ndb;
				}

//        c Early bug - always did OK in presence of one data.
				if(ktype==2) skmean = extest;
				if(ktype==0 || ktype==2) {
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

//        c Go ahead and set up the OK portion of the kriging matrix:
			neq = mdt+na;

//        c Initialize the main kriging matrix:
			for(int i=0; i<neq*neq; i++) {
				a[i] = 0d;
			}

//        c Fill in the kriging matrix:
			for(int i=0; i<na; i++) for(int j=i; j<na; j++) {
				cov = Covariance.cova3(xa[i], ya[i], za[i], xa[j], ya[j], za[j],
						1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];
				a[neq*i+j] = cov;
				a[neq*j+i] = cov;
			}

//        c Fill in the OK unbiasedness portion of the matrix (if not doing SK):
			if(neq>na) {
				for(int i=0; i<na; i++) {
					a[neq*i+na] = unbias;
					a[neq*na+i] = unbias;
				}
			}

//        c Set up the right hand side:
			for(int i=0; i<na; i++) {
				if(ndb<=1) {
					cb = Covariance.cova3(xa[i],ya[i],za[i], xdb[0],ydb[0],zdb[0],
							1,nst,MAXNST, c0,it,cc,aa, 1,MAXROT,rotmat)[0];
				} else {
					cb  = 0d;
					for(int j=0; j<ndb; j++) {
						cb = Covariance.cova3(xa[i], ya[i], za[i], xdb[j], ydb[j], zdb[j],
								1, nst, MAXNST, c0, it, cc, aa, 1, MAXROT, rotmat)[0];
						cb += cov;
						double dx = xa[i] - xdb[j];
						double dy = ya[i] - ydb[j];
						double dz = za[i] - zdb[j];
						if(dx*dx+dy*dy+dz*dz < EPSLON) cb -= c0[0];
					}
					cb /= ndb;
				}
				r[i] = cb;
			}
			if(neq>na) r[na] = unbias;

//        c Add the additional unbiasedness constraints:
			int im = na + 1;

//        c First drift term (linear in "x"):
			for(int dt=0; dt<MAXDT; dt++) {
				if(idrif[dt]==1) {
					for(int k=0; k<na; k++) {
						a[neq*im+k] = drift_term(dt+1, xa[k], ya[k], za[k]) * resc;
						a[neq*k+im] = a[neq*im+k];
					}
					r[im] = bv[dt];
					im++;
				}
			}
//        c
//        c External drift term (specified by external variable):
//        c
			if(ktype==3) {
				//im=im+1
				for(int k=0; k<na; k++) {
					a[neq*im+k] = vea[k]*resce;
					a[neq*k+im] = vea[k]*resce;
				}
				r[im] = extest*resce;
				im++; // copied afterwards for JAVA indices
			}
//        c
//        c Copy the right hand side to compute the kriging variance later:
//        c
			for(int k=0; k<neq; k++) {
				rr[k] = r[k];
			}
			//int kadim = neq * neq;
			//int ksdim = neq;
			int nrhs  = 1;
			int nv    = 1;
//        c
//        c If estimating the trend then reset all the right hand side terms=0.0:
//        c
			if(itrend==1) {
				for(int i=0; i<na; i++) {
					r[i]  = 0d;
					rr[i] = 0d;
				}
			}
//        c
//        c Write out the kriging Matrix if Seriously Debugging:
//        c
			if(idbg==3) {
				try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
					bw.append("\nEstimating node index : "+ix+" "+iy+" "+iz+"\n");
		            int is = 1 - neq,ie;
		            for(int i=0; i<neq; i++) {
		            	is = 1 + i*neq;
		            	ie = is + neq - 1;
		            	String str_a = "";
		            	for(int j=is-1; j<ie; j++) str_a += " "+FormatHelper.nf(a[j],7,4);
		            	bw.append("    r("+FormatHelper.nf(i+1,2)+")= "+FormatHelper.nf(r[i],7,4)+"  a="+str_a+"\n");
		            }
		            bw.flush();
				}catch(IOException io_e) {
					io_e.printStackTrace();
				}
			}
//        c
//        c Solve the kriging system:
//        c
		//      call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
			s = LinEquSolver.ktsol(neq, nrhs, nv, a, r, maxEquations);
//        c
//        c Compute the solution:
//        c
			if(s==null) {
				if(idbg>=3) {
					try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
						bw.append(" Singular Matrix "+(ix+1)+" "+(iy+1)+" "+(iz+1)+"\n");
						bw.flush();
					}catch(IOException io_e) {
						io_e.printStackTrace();
					}
				}
				est  = UNEST;
				estv = UNEST;
			} else {
				est  = 0d;
				estv = cbb;
				if(ktype==2) skmean = extest;
				for(int j=0; j<neq; j++) {
					estv -= s[j]*rr[j];
					if(j<na) {
						if(ktype==0) {
							est += s[j]*(vra[j]-skmean);
						} else if(ktype==2) {
							est += s[j]*(vra[j]-vea[j]);
						} else {
							est += s[j]*vra[j];
						}
					}
				}
				if(ktype==0 || ktype==2) est += skmean;
				nk++;
				xk += est;
				vk += est*est;
//        c
//        c Write the kriging weights and data if debugging level is above 2:
//        c
				if(idbg>=2) {
					try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
						bw.append("\nBLOCK: "+(ix+1)+" "+(iy+1)+" "+(iz+1)+" at "+xloc+" "+yloc+" "+zloc+"\n\n");
						if(ktype!=0)
							bw.append("  Lagrange : "+(s[na]*unbias)+"\n");
						bw.append("  BLOCK EST: x,y,z,vr,wt \n");
						for(int i=0; i<na; i++) {
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
//        c
//        c END OF MAIN KRIGING LOOP:
//        c
		//} // Loop Mark 1: continue;
			double err,wtmin,sumwt;
			if(iktype==0) {
//				try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfl),true))) {
					if(koption==0) {
//						bw.append(" "+FormatHelper.nf(est,14,8)+" "+FormatHelper.nf(estv,14,8)+"\n");
						estimates[0][index] = est;
						estimates[1][index] = estv;
					} else {
						err = UNEST;
						if(_true_!=UNEST && est!=UNEST) {
							err=est-_true_;
							xkmae += Math.abs(err);
							xkmse += err*err;
						}
//						bw.append(" "+FormatHelper.nf(xloc,14,8)+" "+FormatHelper.nf(yloc,14,8)+" "+FormatHelper.nf(zloc,14,8)
//								+ " "+FormatHelper.nf(_true_,14,8)+" "+FormatHelper.nf(est,14,8)+" "+FormatHelper.nf(estv,14,8)
//								+ " "+FormatHelper.nf(err,14,8)+"\n");
						estimates[0][index] = xloc;
						estimates[1][index] = yloc;
						estimates[2][index] = zloc;
						estimates[3][index] = _true_;
						estimates[4][index] = est;
						estimates[5][index] = estv;
						estimates[6][index] = err;
					}
//					bw.flush();
//				}catch(IOException io_e) {
//					io_e.printStackTrace();
//				}
//        c
//        c Work out the IK-type distribution implicit to this data configuration
//        c and kriging weights:
//        c
			} else {
				for(int icut=0; icut<ncut; icut++) {
					cdf[icut] = -1d;
				}
				wtmin = 1d;
				for(int i=0; i<na; i++) {
					if(s[i]<wtmin) wtmin = s[i];
				}
				sumwt = 0d;
				for(int i=0; i<na; i++) {
					s[i]  -= wtmin;
					sumwt += s[i];
				}
				for(int i=0; i<na; i++) {
					s[i] /= Math.max(0.00001d,sumwt);
				}
				if(na>1 && sumwt<0.00001d) {
					for(int icut=0; icut<ncut; icut++) {
						cdf[icut] = 0d;
						for(int i=0; i<na; i++) {
							if(vra[i]<=cut[icut]) cdf[icut] += s[i];
						}
					}
				}
				for(int i=0; i<ncut; i++) estimates[i][index] = cdf[i];
//				try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfl),true))) {
//					String str_cdf = "";
//					for(int i=0; i<ncut; i++) str_cdf+=" "+FormatHelper.nf(cdf[i],8,4);
//					if(koption==0) {
//						bw.append(str_cdf+"\n");
//					} else {
//						bw.append(str_cdf+" "+FormatHelper.nf(_true_,8,4)+"\n");
//					}
//					bw.flush();
//				}catch(IOException io_e) {
//					io_e.printStackTrace();
//				}
			}
//		      end do
		} // 2    continue
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
					bw.append("  mean error "+FormatHelper.nf(xkmae,14,8)+"  mean sqd e "+FormatHelper.nf(xkmse,14,8)+"\n");
				}
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
//        c
//        c All finished the kriging:
//        c
		if(_dataframe_dims==1) {
			DataFrame result = new DataFrame();
			for(int n_est=0; n_est<num_krig_res; n_est++)
				result.addColumn(estimate_titles[n_est], estimates[n_est]);
			return result;
		}
		if(_dataframe_dims==3) {
			int xn=nx, yn=ny, zn=nz;
			if(xn*yn*zn!=estimates[0].length) { nx=estimates[0].length; ny=1; nz=1; }
			DataFrame3D result = new DataFrame3D();
			for(int n_est=0; n_est<num_krig_res; n_est++) {
				double[][][] estimates3D = new double[xn][yn][zn];
				for(int x_=0; x_<xn; x_++) for(int y_=0; y_<yn; y_++) for(int z_=0; z_<zn; z_++)
					estimates3D[x_][y_][z_] = estimates[n_est][x_*yn*zn+y_*zn+z_];
				result.addColumn(estimate_titles[n_est], estimates3D);
			}
			if(nx*ny*nz==estimates[0].length) {
				double[] dimension_x = new double[nx];
				for(int x_=0; x_<nx; x_++) dimension_x[x_] = xmn + x_*xsiz;
				double[] dimension_y = new double[ny];
				for(int y_=0; y_<ny; y_++) dimension_y[y_] = ymn + y_*xsiz;
				double[] dimension_z = new double[nz];
				for(int z_=0; z_<nz; z_++) dimension_z[z_] = zmn + z_*xsiz;
				result.setDimension(1, dimension_x);
				result.setDimension(2, dimension_y);
				result.setDimension(3, dimension_z);
			}
			return result;
		}
		return null;
//		 96   stop 'ERROR in jackknife file!'
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