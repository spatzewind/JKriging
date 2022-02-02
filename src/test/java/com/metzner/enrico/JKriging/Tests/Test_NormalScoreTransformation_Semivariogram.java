package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.probability.ProbabilityTransform;
import com.metzner.enrico.JKriging.variogram.GAMV;

public class Test_NormalScoreTransformation_Semivariogram {

	public static void main(String[] args) {
		DataFrame df = new DataFrame();
		df.read_gam_dat("res/cluster.dat");
		double[][] nstrnsf = ProbabilityTransform.nscore(df, "Primary", "Declustering_Weight", 1.0e-20d, "NS_Primary");
		System.out.println("Test:");
		df.head();
		df.describe();
		System.out.println("\n\n");
		
		GAMV gamv = new GAMV();
//		gamv.readData("res/gamv.par", df);
		gamv.setDataframe(df);
		gamv.setLagDefinition(10, 5.0d, 3.0d);
		gamv.setCoordinates("Xlocation", "Ylocation", null);
		gamv.addGeoDirection(0d, 90d, 50d, 0d, 90d, 50d);
		gamv.addVariogramDef(1, "NS_Primary", "NS_Primary", 0d);
		gamv.standardizeSills(true);
		DataFrame gdf = gamv.gamv();
		ProbabilityTransform.back_nscore(gdf, "mean_head", nstrnsf, "BNS_head_mean",
				ProbabilityTransform.INTPOL_POWER, ProbabilityTransform.INTPOL_HYPERBOLIC, 1.5d, 2d);
		ProbabilityTransform.back_nscore(gdf, "mean_tail", nstrnsf, "BNS_tail_mean",
				ProbabilityTransform.INTPOL_POWER, ProbabilityTransform.INTPOL_HYPERBOLIC, 1.5d, 2d);
		//gamv.writeData("res/gamv.out", true);
		gdf.printfull();
	}

}
