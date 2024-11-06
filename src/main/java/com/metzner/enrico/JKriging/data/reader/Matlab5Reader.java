package com.metzner.enrico.JKriging.data.reader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;

import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.format.Mat5File;
import us.hebi.matlab.mat.types.AbstractCharBase;
import us.hebi.matlab.mat.types.AbstractMatrixBase;
import us.hebi.matlab.mat.types.Array;
import us.hebi.matlab.mat.types.MatlabType;
import us.hebi.matlab.mat.types.Struct;

class Matlab5Reader extends DataReader {
	
	private static final int MAX_READ_DIM_COUNT = 20;
	
	private static final String MAT5_HEADER_START = "MATLAB 5.0 MAT-file";
	private static final String MAT7_HEADER_START = "MATLAB 7.";
	
	private Mat5File mafile = null;
	
	Matlab5Reader(File file) throws IOException {
		super(file);
		mafile = Mat5.readFromFile(file);
	}
	
	@Override
	public void describeContent() {
		System.out.println(mafile);
	}
	
	@Override
	public Map<String, Integer> getDimensionsOfVariable(String var) throws IllegalAccessException {
		if(mafile==null) throw new IllegalAccessException("The Matlab file does not exist.");
		if(var==null) throw new IllegalAccessException("Variable cannot be null");
		Array a = mafile.getArray(var);
//		if(a==null) throw new IllegalAccessException("Cannot find variable <"+var+"> in the Matlab file.");
		if(a==null) return null;
		Map<String, Integer> dimmap = new HashMap<>();
		if(decompose && a instanceof Struct) {
			dimmap.put("STRUCT", 0);
		} else {
			for(int d=0; d<a.getNumDimensions(); d++)
				dimmap.put("dim"+d, a.getDimensions()[d]);
		}
		return dimmap;
	}
	
	@Override
	public DataFrame   getVars1D(String... vars) throws IllegalAccessException {
		if(mafile==null) throw new IllegalAccessException("The Matlab file does not exist.");
		List<MatArray> mavars = new ArrayList<MatArray>();
		int[] dimensions = new int[MAX_READ_DIM_COUNT];
		collectVars(vars, dimensions, mavars);
		if(mavars.size()==0)
			throw new IllegalArgumentException("Matlab file does not contain any of given variables.");
		long datasize = 1L;
		int dlen = 0;
		for(int d: dimensions) {
			datasize *= Math.max(d,1);
			if(datasize>MAX_DATAFRAME_LENGTH)
				throw new ArrayStoreException("Cannot expand combination of dimensions, maximum length of dataframe reached!");
			if(d>0) dlen++;
		}
		DataFrame df = new DataFrame();
		int datalength = 1;
		int[] dimlen = new int[dlen];
		for(int d=0; d<dlen; d++) { dimlen[d] = dimensions[d]; datalength *= dimlen[d]; }
		int[][] indices = new int[datalength][dlen];
		for(int dl=0; dl<datalength; dl++) {
			int sublen = 1;
			for(int d=dlen-1; d>=0; d--) {
				indices[dl][d] = (dl / sublen) % dimlen[d];
				sublen *= dimlen[d];
			}
		}
		for(MatArray s: mavars) {
			switch(s.getType()) {
				case Int8:
				case UInt8: byte[] arrB = new byte[datalength]; for(int i=0; i<datalength; i++)
					arrB[i] = s.getByte(s.getPos(indices[i])); df.addColumn(s.varname, arrB); break;
				case Int16:
				case UInt16: short[] arrR = new short[datalength]; for(int i=0; i<datalength; i++)
					arrR[i] = s.getShort(s.getPos(indices[i])); df.addColumn(s.varname, arrR); break;
				case Int32:
				case UInt32: int[] arrI = new int[datalength]; for(int i=0; i<datalength; i++)
					arrI[i] = s.getInt(s.getPos(indices[i])); df.addColumn(s.varname, arrI); break;
				case Int64:
				case UInt64: long[] arrL = new long[datalength]; for(int i=0; i<datalength; i++)
					arrL[i] = s.getLong(s.getPos(indices[i])); df.addColumn(s.varname, arrL); break;
				case Single: float[] arrF = new float[datalength]; for(int i=0; i<datalength; i++)
					arrF[i] = s.getFloat(s.getPos(indices[i])); df.addColumn(s.varname, arrF); break;
				case Double: double[] arrD = new double[datalength]; for(int i=0; i<datalength; i++)
					arrD[i] = s.getDouble(s.getPos(indices[i])); df.addColumn(s.varname, arrD); break;
				case Character: //convert to Strings
					String[] arrS = new String[datalength]; for(int i=0; i<datalength; i++)
					arrS[i] = s.getString(s.getPos(indices[i])); df.addColumn(s.varname, arrS); break;
//				case Cell:
//					break;
//				case Function:
//					break;
//				case Object:
//					break;
//				case Opaque:
//					break;
//				case Sparse:
//					break;
				case Structure:
					System.err.println("WARNING: could not add variable \""+s.varname+
							"\": structs have to be decomposed befor usage.");
					break;
				default:
					System.err.println("WARNING: could not add variable \""+s.varname+
							"\": unsupported data type.");
					break;
			}
		}
		return df;
	}
	@Override
	public DataFrame2D getVars2D(String... vars) throws IllegalAccessException {
		if(mafile==null) throw new IllegalAccessException("The Matlab file does not exist.");
		List<MatArray> mavars = new ArrayList<MatArray>();
		int[] dimensions = new int[MAX_READ_DIM_COUNT];
		collectVars(vars, dimensions, mavars);
		if(mavars.size()==0)
			throw new IllegalArgumentException("Matlab file does not contain any structs of given variables.");
		
		//TODO
		return null;
	}
	@Override
	public DataFrame3D getVars3D(String... vars) throws IllegalAccessException {
		if(mafile==null) throw new IllegalAccessException("The Matlab file does not exist.");
		List<MatArray> mavars = new ArrayList<MatArray>();
		int[] dimensions = new int[MAX_READ_DIM_COUNT];
		collectVars(vars, dimensions, mavars);
		if(mavars.size()==0)
			throw new IllegalArgumentException("Matlab file does not contain any structs of given variables.");
		
		//TODO
		return null;
	}
	private void collectVars(String[] names, int[] dimlengths, List<MatArray> vars) {
		for(int i=0; i<dimlengths.length; i++) dimlengths[i] = 0;
		for(String name: names) {
			Array var = mafile.getArray(name);
//			System.out.println("  Array-var: "+var);
			if(var==null) continue;
			Array[] arrays;
			String[] vnames;
			if(decompose && var instanceof Struct) { //decompose structures, only one level
				Struct s = (Struct) var;
				String[] fields = s.getFieldNames().toArray(new String[0]);
				arrays = new Array[fields.length];
				vnames = new String[fields.length];
				for(int f=0; f<fields.length; f++) {
					Array field = s.get(fields[f]);
					if(field instanceof AbstractCharBase)
						arrays[f] = new MatString((AbstractCharBase)field);
					else
						arrays[f] = field;
					vnames[f] = name+"_"+fields[f];
				}
			} else
			if(var instanceof AbstractCharBase) {
				arrays = new Array[] {new MatString((AbstractCharBase) var)};
				vnames = new String[] {name};
			} else {
				arrays = new Array[] {var};
				vnames = new String[] {name};
			}
//			System.out.println("  decomposed to: "+arrays+" ("+arrays.length+")");
			for(int i=0; i<arrays.length; i++) {
				Array a = arrays[i];
				int[] dims = new int[a.getNumDimensions()];
				int[] lens = a.getDimensions();
				boolean fit = true;
				for(int d=0; d<dims.length; d++) {
					int pos = -1;
					for(int l=0; l<dimlengths.length && pos<0; l++) {
						if(dimlengths[l]>0 && dimlengths[l]!=lens[d])
							continue;
						boolean possible = true;
						for(int k=0; k<d; k++)
							if(dims[k]==l)
								possible = false;
						if(!possible)
							continue;
						pos = l;
					}
					if(pos<0) {
						fit = false;
						break;
					}
					dims[d] = pos;
					if(dimlengths[pos]==0)
						dimlengths[pos] = lens[d];
				}
				if(!fit) continue;
				vars.add(new MatArray(vnames[i], a, dims));
			}
		}
	}
	
	public DataFrame get1Dslice(String filtervar, int index, String... vars) throws IllegalAccessException {
		System.err.println("WARNING: get1Dslice for Matlab files not implemented yet!");
		DataFrame res = new DataFrame();
		//TODO get1Dslice matlab
		return res;
	}
	
	protected static int canRead(File file) {
		try(BufferedReader br = new BufferedReader(new FileReader(file))) {
			char[] buffer = new char[116];
			int readlength = br.read(buffer);
			if(readlength==0) return 0;
			String header = new String(buffer);
			if(header.startsWith(MAT7_HEADER_START)) {
//				System.out.println("INFO: Matlab version -v7 and higher are probably true HDF-files and are read like Netcdf-files.");
				return -1;
			}
			if(header.startsWith(MAT5_HEADER_START)) {
				if(!Constants.supportMatlab) {
					System.err.println("WARNING: Matlab version -v5 and -v6 can be read, if mfl-core.jar from "+
							"<a href=\"https://github.com/HebiRobotics/MFL\">HebiRobotics</a> is in class path.");
					return 0;
				}
				return 1;
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return 0;
	}
	protected static boolean hasValidExtension(String filename) {
		if(filename.endsWith(".m")) return true;
		if(filename.endsWith(".mat")) return true;
		return false;
	}
	
	
	@Override
	public void close() throws IOException {
		if(mafile!=null) mafile.close();
		mafile = null;
	}
	
	private class MatArray {
		public String varname;
		public int[] dimids;
		public Array array;
		public MatArray(String n, Array a, int[] d) {
			varname = n;
			array = a;
			dimids = d;
		}
		
		public MatlabType getType() {
			return array.getType();
		}
		
		public int[] getPos(int[] indices) {
			if(dimids.length==1) return new int[] {indices[dimids[0]]};
			if(dimids.length==2) return new int[] {indices[dimids[0]], indices[dimids[1]]};
			int[] pos = new int[dimids.length];
			for(int i=0; i<dimids.length; i++)
				pos[i] = indices[dimids[i]];
			return pos;
		}
		public boolean getBoolean(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getBoolean(pos);
			return false;
		}
		public byte getByte(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getByte(pos);
			return Byte.MAX_VALUE;
		}
		public short getShort(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getShort(pos);
			return Short.MAX_VALUE;
		}
		public int getInt(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getInt(pos);
			return Integer.MAX_VALUE;
		}
		public long getLong(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getLong(pos);
			return Long.MAX_VALUE;
		}
		public float getFloat(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getFloat(pos);
			return Float.NaN;
		}
		public double getDouble(int... pos) {
			if(array instanceof AbstractMatrixBase) return ((AbstractMatrixBase)array).getDouble(pos);
			return Double.NaN;
		}
		public String getString(int... pos) {
			if(!(array instanceof MatString)) return "";
			return ((MatString) array).getString(pos[0]);
		}
		
		@Override
		public String toString() {
			String dimstring = "";
			for(int i=0; i<dimids.length; i++) dimstring += (i==0?"":",")+dimids[i];
			return "["+varname+", "+array.getClass().getSimpleName()+", ("+dimstring+")]";
		}
	}
	private class MatString implements Array {
		
		AbstractCharBase srcChar;
		public MatString(AbstractCharBase source) {
			srcChar = source;
		}
		
		public String getString(int index) {
			return srcChar.getRow(index);
		}
		
		@Override
		public void close() throws IOException {
			srcChar.close();
		}
		@Override
		public MatlabType getType() {
			return MatlabType.Character;
		}
		@Override
		public int[] getDimensions() {
			return new int[] {srcChar.getNumRows(),1};
		}
		@Override
		public int getNumDimensions() {
			return 2;
		}
		@Override
		public int getNumRows() {
			return srcChar.getNumRows();
		}
		@Override
		public int getNumCols() {
			return 1;
		}
		@Override
		public int getNumElements() {
			return srcChar.getNumRows();
		}
	}
}
