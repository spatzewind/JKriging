package com.metzner.enrico.JKriging.Tests;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.kriging.KT3D;
import com.metzner.enrico.JKriging.variogram.GAMV;

import ucar.ma2.InvalidRangeException;

public class Test_3DOrdinaryKriging_v2 {

	public static void main(String[] args) {
		
		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *
		 * +++                                                     +++ *
		 * +++  FIRST CREATE A SEMIVARIOGRAM                       +++ *
		 * +++                                                     +++ *
		 * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
		
		File f = new File("res/alex3D.par");

		try(BufferedWriter bw = new BufferedWriter(new FileWriter(f))) {
			bw.append("                  Parameters for GAMV\r\n" + 
					  "                  **************************\r\n" + 
					  "\r\n" + 
					  "START OF PARAMETERS:\r\n" + 
					  "res/alex3D.dat\r\n" + 
					  "1 2 3                                   -   columns for X, Y and Z coordinates\r\n" + 
					  "1 4 \r\n" + 
					  "-1E+21 1E+21                            -trimming limits\r\n" + 
					  "res/GV_alex3D.out\r\n" + 
					  "10                                      -number of lags\r\n" + 
					  "5.0                                     -lag separation distance\r\n" + 
					  "3.0                                     -lag tolerance\r\n" + 
					  "3                                       -number of Directions\r\n" + 
					  "0 45 10000 0 45 10000                   -azm, atol, bandh, dip, dtol, bandv\r\n" + 
					  "90 45 10000 0 45 10000                  -azm, atol, bandh, dip, dtol, bandv\r\n" + 
					  "0 90 10000 90 45 10000                  -azm, atol, bandh, dip, dtol, bandv\r\n" + 
					  "1                                       -standardize sill? (0=no, 1=yes)\r\n" + 
					  "1                                       -number of variograms\r\n" + 
					  "1 1 1                                   -Tail var, head var, vario type\r\n");
			bw.flush();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
		
		DataFrame df_alex = new DataFrame();
		df_alex.read_gam_dat("res/alex3D.dat");
		
		GAMV gamv = new GAMV();
		gamv.readParameterfile(f.getAbsolutePath(), df_alex);
		DataFrame df_gm = gamv.gamv();
		df_gm.printfull();
		
		System.out.println("\nDeleted temporary parameter file : "+f.delete());
		
		
		
		
		

		
		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *
		 * +++                                                     +++ *
		 * +++  SECOND SETUP AND CREATE 3D-KRIGING                 +++ *
		 * +++                                                     +++ *
		 * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

//		File k = new File("res/K3_alex3D.par");
//
//		try(BufferedWriter bw = new BufferedWriter(new FileWriter(k))) {
//			bw.append("                  Parameters for KT3D\r\n" + 
//					  "                  **************************\r\n" + 
//					  "\r\n" + 
//					  "START OF PARAMETERS:\r\n" + 
//					  "res/alex3D.dat\r\n" + 
//					  "0 1 2 3 4 0                             -   Columns for DH, X, Y, Z, var, sec var\r\n" + 
//					  "-1E+21 1E+21                            -trimming limits\r\n" + 
//					  "0                                       -option: 0=grid, 1=cross, 2=jackknife\r\n" + 
//					  "d:\\GEOGRAPHIE\\data\\xvk.dat\r\n" + 
//					  "1 2 0 3 0                               -   Columns for X, Y, Z, vr and sec var\r\n" + 
//					  "1                                       -debug level (0-3)\r\n" + 
//					  "res/K3_alex3D.dbg\r\n" + 
//					  "res/K3_alex3D.out\r\n" + 
//					  "10 2.5 5                                -nx, xmin, xsize\r\n" + 
//					  "10 2.5 5                                -ny, ymin, ysize\r\n" + 
//					  "10 2.5 5                                -nz, zmin, zsize\r\n" + 
//					  "1 1 1                                   -x, y and z block discretization\r\n" + 
//					  "4 8                                     -Min and max data for kriging\r\n" + 
//					  "0                                       -Max per octant (0-> not used)\r\n" + 
//					  "20 20 20                                -Maximum search radii\r\n" + 
//					  "0 0 0                                   -Angles for search ellipsoid\r\n" + 
//					  "1 0                                     -0=SK, 1=OK, 2=non-st SK,3=exdrift\r\n" + 
//					  "0 0 0 0 0 0 0 0 0                       -drift: x,y,z,xx,yy,zz,xy,xz,zy\r\n" + 
//					  "0                                       -Variable; 1, estimate; trend\r\n" + 
//					  "d:\\GEOGRAPHIE\\data\\extdrift.dat -gridded file with drift/mean\r\n" + 
//					  "4                                       -   Column number in gridded file\r\n" + 
//					  "1 0                                     -nst, nugget NOFILE\r\n" + 
//					  "3 1 0 0 0                               -it, cc, ang1, ang2, ang3\r\n" + 
//					  "15 11.5 15                              -a_hmax, a_hmin, a_vert");
//			bw.flush();
//		} catch(IOException ioe) {
//			ioe.printStackTrace();
//		}
		
		KT3D kt3d = new KT3D();
//		kt3d.readData(k.getAbsolutePath());
		kt3d.setKrigingOptionAndDataframe(df_alex, KT3D.GRID, KT3D.ORDINARY_KRIGING, 0d, false, null);
		kt3d.setVariableAndCoordinates("Xlocation", "Ylocation", "Zlocation", "Primary", null);
		kt3d.setKrigingField(10, 5.0d, 2.5d, 10, 5.0d, 2.5d, 10, 5.0d, 2.5d);
		kt3d.setSuperBlockParameters(1, 1, 1, 4, 8, 0, 20d);
		kt3d.addVariogramModelMath(3, 1.0d, 0d, 0d, 0d, 15d, 11.5d, 15d);
//		kt3d.kt3d();
		
		DataFrame3D df_K_alex = kt3d.kt3d_asDataFrame3D();
		//df_alex.read_gam_dat("res/K3_alex3D.out");
		//df_K_alex.printfull();
		df_K_alex.head();
		try {
			df_K_alex.writeToNetcdf("res/K3_alex_out.nc");
		} catch (IOException | InvalidRangeException e) {
			e.printStackTrace();
		}
		
//		System.out.println("\nDeleted temporary parameter file : "+k.delete());
		
		
		
	}

}
