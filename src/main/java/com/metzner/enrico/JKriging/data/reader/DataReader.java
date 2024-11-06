package com.metzner.enrico.JKriging.data.reader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.error.UnknownFileFormatException;

public abstract class DataReader implements Closeable {
	
	private static List<DataReader> reader;
	static {
		reader = new ArrayList<>();
	}
	
	public static final int MAX_DATAFRAME_LENGTH = Integer.MAX_VALUE>>3;
	protected static final String wT = "\u251C"; // |-
	protected static final String wP = "\u2502"; // |
	protected static final String wL = "\u2500"; // --
	protected static final String wE = "\u2514"; // `-
	
	protected boolean decompose;
	protected File associated_file;
	
	/**
	 protected constructor with parameter "file"
	 <p>
	   each class extending from this has to match this constructor for being initializable
	   by static methods {@link #openFile(File) openFile(File file)} and {@link #openFile(String) openFile(String filepath)}.
	 </p>
	 
	 @param file file, which to read by DataReader
	 */
	protected DataReader(File file) {
		decompose = false;
	}
	/**
	 * each class extending from this has to override this method
	 * @param name filename with extension, which is checked by this function
	 * @return true, if filename-extension is valid and the DataReader should be able to open/read the file
	 */
	protected static boolean hasValidExtension(String name) {
		return false;
	}
	
	public static DataReader openFile(String filepath) throws IOException,UnknownFileFormatException {
		if(filepath==null) throw new NullPointerException("filepath must not be null!");
		return openFile(new File(filepath));
	}
	public static DataReader openFile(File file) throws IOException,UnknownFileFormatException {
		if(file==null) throw new NullPointerException("file must not be null!");
		if(file.isDirectory()) throw new IOException("Cannot read a directory.");
		String name = file.getName();
		
		if(SpreadSheetReader.hasValidExtension(name))
			return new SpreadSheetReader(file);
		if(ExcelReader.hasValidExtension(name))
			return new ExcelReader(file);
		if(Matlab5Reader.hasValidExtension(name)) {
			int cr = Matlab5Reader.canRead(file);
			if(cr==1) return new Matlab5Reader(file);
			if(cr==-1) return new Matlab7Reader(file);
		}
		if(NetcdfReader.hasValidExtension(name)) {
			if(NetcdfReader.canRead(file))
				return new NetcdfReader(file);
		}
		throw new UnknownFileFormatException("Cannot find reader for file "+file.getAbsolutePath());
	};
//	public static DataReader openFileAsNetcdf(File file) throws IOException {
//		if(file==null) throw new NullPointerException("file must not be null!");
//		if(file.isDirectory()) throw new IOException("Cannot read a directory.");
//	}
	
	public abstract void describeContent();
	public abstract Map<String, Integer> getDimensionsOfVariable(String var) throws IllegalAccessException;
	public abstract DataFrame   getVars1D(String... vars) throws IllegalAccessException;
	public abstract DataFrame2D getVars2D(String... vars) throws IllegalAccessException;
	public abstract DataFrame3D getVars3D(String... vars) throws IllegalAccessException;
	public abstract DataFrame   get1Dslice(String filtervar, int index, String... vars) throws IllegalAccessException;
	@Override
	public abstract void close() throws IOException;
	
	public void decomposeStructures(boolean b) {
		decompose = b;
	}
}
