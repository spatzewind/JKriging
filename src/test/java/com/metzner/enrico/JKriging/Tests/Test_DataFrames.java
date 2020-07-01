package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.helper.DataHelper;

public class Test_DataFrames {

	public static void main(String[] args) {
		System.out.println(
				"    Test: Dataframe functionality\n"+
				"----=============================----");
		DataFrame df = new DataFrame();
		df.setDefaultDatatype("double");
		System.out.println("Default datatype: "+df.getDefaultDatatype());
		df.read_csv("res/sample_data_MV_biased.csv", 1, null);
		System.out.println("\n");
		df.head();
		System.out.println("\n");
		df.describe();
		
		System.out.println("\n\n\n"+
				"    Test: Sorting\n"+
				"----=============----");
		
		double[] test_arr = new double[20],
				 test_rra = new double[20],
				 temparray = new double[20];
		for(int i=0; i<20; i++) {
			test_arr[i] = Math.random();
			test_rra[i] = test_arr[i];
		}
		DataHelper.sortem(test_rra, temparray);
		DataFrame df_sort = new DataFrame();
		df_sort.addColumn("unsorted", test_arr);
		df_sort.addColumn("sorted", test_rra);
		df_sort.printfull();
	}

}
