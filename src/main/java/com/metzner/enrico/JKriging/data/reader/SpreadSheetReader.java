package com.metzner.enrico.JKriging.data.reader;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.github.miachm.sods.Range;
import com.github.miachm.sods.Sheet;
import com.github.miachm.sods.SpreadSheet;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.error.DimensionMismatchException;

public class SpreadSheetReader extends DataReader {
	
	private SpreadSheet odsfile = null;
	
	public SpreadSheetReader(File file) throws IOException {
		super(file);
		odsfile = new SpreadSheet(file);
		readStructure();
	}
	
	@Override
	public void describeContent() {
		System.out.println("Spreadsheet file");
		System.out.println(odsfile);
//		System.out.println(wT+wL+"Path: "+odsfile.);
//		System.out.println(wT+wL+"Name: "+wbfile.getName());
	}
	
	@Override
	public Map<String, Integer> getDimensionsOfVariable(String var) throws IllegalAccessException {
		return null;
	}

	@Override
	public DataFrame getVars1D(String... vars) throws IllegalAccessException {
		return null;
	}
	@Override
	public DataFrame2D getVars2D(String... vars) throws IllegalAccessException {
		return null;
	}
	@Override
	public DataFrame3D getVars3D(String... vars) throws IllegalAccessException {
		throw new DimensionMismatchException("Spreadsheets cannot contain 3D data");
	}
	@Override
	public DataFrame get1Dslice(String filtervar, int index, String... vars) throws IllegalAccessException {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public void close() throws IOException {
		odsfile = null;
	}
	
	protected static boolean hasValidExtension(String filename) {
//		if(filename.endsWith(".xls")) return true;
//		if(filename.endsWith(".xlsx")) return true;
		if(filename.endsWith(".ods")) return true;
		return false;
	}
	
	
	
	private void readStructure() {
		List<Sheet> sheets = odsfile.getSheets();
		for(Sheet sheet: sheets) {
			 System.out.println("In sheet " + sheet.getName());

             Range range = sheet.getDataRange();
             System.out.println(range.toString());
		}
	}
}
