package com.metzner.enrico.JKriging.data.reader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.error.UnknownFileFormatException;

public abstract class DataReader implements Closeable {
	
	public static final int MAX_DATAFRAME_LENGTH = Integer.MAX_VALUE>>3;
	
	protected boolean decompose;
	
	protected DataReader() {
		decompose = false;
	}
	
	public static DataReader openFile(String filepath) throws IOException,UnknownFileFormatException {
		if(filepath==null) throw new NullPointerException("filepath must not be null!");
		return openFile(new File(filepath));
	}
	public static DataReader openFile(File file) throws IOException,UnknownFileFormatException {
		if(file==null) throw new NullPointerException("file must not be null!");
		if(file.isDirectory()) throw new IOException("Cannot read a directory.");
		String name = file.getName();
		if(MatlabReader.hasValidExtension(name)) {
			int cr = MatlabReader.canRead(file);
			if(cr==1) return new MatlabReader(file);
			if(cr==-1) return new NetcdfReader(file);
		}
		if(NetcdfReader.hasValidExtension(name)) {
			if(NetcdfReader.canRead(file))
				return new NetcdfReader(file);
		}
		throw new UnknownFileFormatException("Cannot find reader for file "+file.getAbsolutePath());
	};
	
	public abstract void describeContent();
	public abstract DataFrame   getVars1D(String... vars) throws IllegalAccessException;
	public abstract DataFrame2D getVars2D(String... vars) throws IllegalAccessException;
	public abstract DataFrame3D getVars3D(String... vars) throws IllegalAccessException;
	@Override
	public abstract void close() throws IOException;
	
	public void decomposeStructures(boolean b) {
		decompose = b;
	}
}
