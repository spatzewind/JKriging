package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.kriging.KT3D;

public class Test_3DOrdinaryKriging {

	public static void main(String[] args) {
		DataFrame df_cluster = new DataFrame();
		df_cluster.read_gam_dat("res/cluster.dat");
		
		KT3D kt3d = new KT3D();
		//kt3d.loadParameterFile("res/test_kt3d.par", df_cluster, null, null);
		kt3d.setKrigingOptionAndDataframe(df_cluster, KT3D.GRID, KT3D.ORDINARY_KRIGING, 0d, false, null);
		kt3d.setVariableAndCoordinates("Xlocation", "Ylocation", null, "Primary", null);
		kt3d.setKrigingField(2, 25.0d, 12.5d, 2, 25.0d, 12.5d, 1, 1.0d, 0.0d);
		kt3d.setSuperBlockParameters(1, 1, 1, 4, 8, 0, 20d);
		kt3d.addVariogramModelMath(1, 0.8f, 0d, 0d, 0d, 20d, 20d, 20d);
		kt3d.setDebuggingLevel(3, "./res/kt3d_debug_cluster.dbg");
		DataFrame res = kt3d.kt3d();
		
		res.head();
		res.describe();
	}

}
