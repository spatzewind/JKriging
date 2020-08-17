package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;

public class Test_ConditionalSubDataframe {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Constants.setFirstIndex(Constants.JAVA_INDEXING);
		
		DataFrame df_test = new DataFrame();
		double[] x = new double[100],
				 y = new double[100];
		boolean[] b = new boolean[100];
		for(int i=0; i<100; i++) {
			x[i] = 0.25d * i * Math.PI;
			y[i] = Math.cos(x[i]);
			b[i] = ((i%5)==0);
		}
		df_test.addColumn("x", x);
		df_test.addColumn("y", y);
		df_test.addColumn("fifth", b);

		df_test.head();
		df_test.describe();
		
		DataFrame subdf = df_test.filterSubDataFrame("$0>10 || @y<0.0 && ~@fifth");
		subdf.head();
		subdf.describe();
	}

}
