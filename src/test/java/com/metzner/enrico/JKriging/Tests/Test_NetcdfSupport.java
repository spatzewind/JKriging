package com.metzner.enrico.JKriging.Tests;

import java.io.IOException;

import com.metzner.enrico.JKriging.data.DataFrame;

import ucar.nc2.NetcdfFile;

public class Test_NetcdfSupport {

	public static void main(String[] args) {
		
		/* +--------------------------------------------+ *
		 * |                                            | *
		 * |   Test Netcdf-Support                      | *
		 * |                                            | *
		 * +--------------------------------------------+ */
		
		NetcdfFile ncdf = null;
		
		try {
			ncdf = NetcdfFile.open("res/test2d.nc");
			System.out.println(ncdf);
	
			DataFrame df_nc = new DataFrame();
			df_nc.readFromNetcdf(ncdf, true, "test2d");
			
			df_nc.head(0);
			
			df_nc.writeToNetcdf("res/test2d_dataframe.nc");
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
}
