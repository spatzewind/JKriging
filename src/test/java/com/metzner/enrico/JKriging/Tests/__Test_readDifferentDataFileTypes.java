package com.metzner.enrico.JKriging.Tests;

import com.metzner.enrico.JKriging.data.reader.DataReader;

public class __Test_readDifferentDataFileTypes {

	public static void main(String[] args) {
		try(DataReader dr = DataReader.openFile("/home/emetzner/Uni/Ming/analysis/Kopie von Auswertung.xlsx")) {
			dr.describeContent();
		} catch(Exception e) {
			e.printStackTrace();
		}
	}

}
