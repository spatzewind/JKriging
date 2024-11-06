package com.metzner.enrico.JKriging.Tests;

import java.io.IOException;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.reader.DataReader;
import com.metzner.enrico.JKriging.error.UnknownFileFormatException;

public class Test_NetcdfSupport {
	
	public static void main(String[] args) {
		
		/* +--------------------------------------------+ *
		 * |                                            | *
		 * |   Test Netcdf-Support                      | *
		 * |                                            | *
		 * +--------------------------------------------+ */
		
		DataFrame df = null;
		
		try(DataReader dr = DataReader.openFile("res/test2d.nc")) {
			dr.describeContent();
			
			df = dr.getVars1D("test2d");
			if(df!=null) {
				System.out.println("Description:");
				df.describe();
				System.out.println("First content:");
				df.head(0);
			}
			
//			df_nc.writeToNetcdf("res/test2d_dataframe.nc");
		} catch (IOException ioe) {
			ioe.printStackTrace();
		} catch (UnknownFileFormatException uffe) {
			uffe.printStackTrace();
		} catch (IllegalAccessException iae) {
			iae.printStackTrace();
		}
	}
}
