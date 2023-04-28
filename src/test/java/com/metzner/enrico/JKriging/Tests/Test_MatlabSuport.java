package com.metzner.enrico.JKriging.Tests;

import java.io.IOException;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.reader.DataReader;
import com.metzner.enrico.JKriging.error.UnknownFileFormatException;

public class Test_MatlabSuport {
	
	public static void main(String[] args) {
		
		/* +--------------------------------------------+ *
		 * |                                            | *
		 * |   Test Matlab-Support                      | *
		 * |                                            | *
		 * +--------------------------------------------+ */
		
		DataFrame df = null;
		try(DataReader dr = DataReader.openFile("res/data.mat")) {
			System.out.println(dr);
			
			dr.decomposeStructures(false);
			df = dr.getVars1D("var1", "var2", "var3", "var4");
			
			if(df!=null) {
				System.out.println("Description:");
				df.describe();
				System.out.println("First content:");
				df.head();
			}
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
		} catch (UnknownFileFormatException uffe) {
			uffe.printStackTrace();
		} catch (IllegalAccessException iae) {
			iae.printStackTrace();
		}
	}

}
