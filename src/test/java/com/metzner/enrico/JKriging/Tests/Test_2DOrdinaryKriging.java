package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.kriging.KB2D;

public class Test_2DOrdinaryKriging {

	public static void main(String[] args) {
		
		KB2D kb2d = new KB2D();
		kb2d.makepar();
		kb2d.readParam("res/kb2d.par");
		kb2d.kb2d();
	}

}
