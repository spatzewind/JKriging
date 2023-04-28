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
import com.metzner.enrico.JKriging.error.DimensionMismatchException;

import ucar.ma2.Array;
import ucar.ma2.Index;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;

class NetcdfReader extends DataReader {
	
	private NetcdfFile ncfile = null;
	
	NetcdfReader(File file) throws IOException {
		super();
		ncfile = NetcdfFiles.open(file.getAbsolutePath());
	}
	
	@Override
	public void describeContent() {
		System.out.println(ncfile);
	}
	
	public DataFrame   getVars1D(String... vars) throws IllegalAccessException {
		if(ncfile==null) throw new IllegalAccessException("The Netcdf file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		List<Dimension> dims = new ArrayList<Dimension>();
		List<Variable> ncvars = new ArrayList<Variable>();
		collectVarsAndDims(vars, dims, ncvars);
		if(ncvars.size()==0)
			throw new IllegalArgumentException("Netcdf file does not contain any of given variables.");
		//check size!
		long datasize = 1L;
		for(Dimension d: dims) {
			datasize *= d.getLength();
			if(datasize>MAX_DATAFRAME_LENGTH)
				throw new ArrayStoreException("Cannot expand combination of dimensions, maximum length of dataframe reached!");
		}
//		for(Dimension d: dims) {
//			Variable dimvar = ncfile.findVariable(d.getName());
//			if(dimvar==null) continue;
//			if(!ncvars.contains(dimvar)) ncvars.add(dimvar);
//			dims.remove(d);
//		}
		int datalength = 1;
		int dlen = dims.size();
		int[] dimlen = new int[dlen];
		for(int d=0; d<dlen; d++) { dimlen[d] = dims.get(d).getLength(); datalength *= dimlen[d]; }
		int[][] indices = new int[datalength][dlen];
		for(int dl=0; dl<datalength; dl++) {
			int sublen = 1;
			for(int d=dlen-1; d>=0; d--) {
				indices[dl][d] = (dl / sublen) % dimlen[d];
				sublen *= dimlen[d];
			}
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(2, indices);
		DataFrame df = new DataFrame();
		int dimid = 0;
		for(Dimension dim: dims) {
			Variable var = ncfile.findVariable(dim.getName());
			if(var==null) {
				int[] arr = new int[datalength];
				for(int dl=0; dl<datalength; dl++) arr[dl] = indices[dl][dimid];
				df.addColumn(dim.getName(), arr);
			} else {
				Array a = null;
				try {
					a = var.read();
				} catch (IOException e) {
					System.out.println("WARNING: could not read dimension \""+dim.getName()+"\": does not add it to dataframe");
					continue;
				}
				Index ind = a.getIndex();
				switch(var.getDataType()) {
					case BOOLEAN:
						boolean[] bl = new boolean[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); bl[dl] = a.getBoolean(ind); }
						df.addColumn(dim.getName(), bl);
						break;
					case BYTE:
						byte[] etyb = new byte[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); etyb[dl] = a.getByte(ind); }
						df.addColumn(dim.getName(), etyb);
						break;
					case INT:
						int[] tni = new int[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); tni[dl] = a.getInt(ind); }
						df.addColumn(dim.getName(), tni);
						break;
					case SHORT:
						short[] trohs = new short[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); trohs[dl] = a.getShort(ind); }
						df.addColumn(dim.getName(), trohs);
						break;
					case LONG:
						long[] gnol = new long[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnol[dl] = a.getLong(ind); }
						df.addColumn(dim.getName(), gnol);
						break;
					case FLOAT:
						boolean hasFVf = (var.findAttribute("_FillValue")!=null);
						double ffill = Double.NaN; if(hasFVf) ffill = (float) var.findAttribute("_FillValue").getNumericValue();
						if(hasFVf) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+ffill);
						float[] taolf = new float[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); taolf[dl] = a.getFloat(ind);
							if(hasFVf && taolf[dl]==ffill) taolf[dl]=Float.NaN; }
						df.addColumn(dim.getName(), taolf);
						break;
					case DOUBLE:
						boolean hasFVd = (var.findAttribute("_FillValue")!=null);
						double dfill = Double.NaN; if(hasFVd) dfill = (double) var.findAttribute("_FillValue").getNumericValue();
						if(hasFVd) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+dfill);
						double[] elbuod = new double[datalength];
						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); elbuod[dl] = a.getDouble(ind);
							if(hasFVd && elbuod[dl]==dfill) elbuod[dl] = Double.NaN; }
						df.addColumn(dim.getName(), elbuod);
						break;
//					case CHAR:
//						break;
//					case STRING:
//						String[] gnirts = new String[datalength];
//						for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnirts[dl] = ""+a.getChar(ind); }
//						addColumn(dim.getFullName(), gnirts);
//						break;
					default:
						System.err.println("WARNING: could not add variable \""+dim.getName()+"\": unsupported data type!");
						break;
				}
			}
			dimid++;
		}
		for(Variable var: ncvars) {
			String vname = var.getFullName();
			if(df.hasVariable(vname))
				continue;
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
				//ioe.printStackTrace();
				System.err.println("WARNING: could not read variable \""+vname+"\": does not at to the dataframe!");
				continue;
			}
			Index ind = a.getIndex();
			List<Dimension> vdims = var.getDimensions();
			int[] dimids = new int[vdims.size()];
			for(int d=0; d<dimids.length; d++) {
				dimids[d] = dims.indexOf(vdims.get(d));
//				dimids[d] = 0;
//				int dd = -1;
//				for(dd=0; dd<dlen; dd++) if(dims.get(dd).equals(vdims.get(d))) {
//					dimids[d] = dd;
//					break;
//				}
//				if(dd==dlen) System.err.println("WARNING: could not find dimension \""+vdims.get(d).getName()+"\": maybe read data in wrong order!");
			}
			int[] indexindex = new int[vdims.size()];
			switch(var.getDataType()) {
				case BOOLEAN: boolean[] o_arr = new boolean[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); o_arr[dl] = a.getBoolean(ind);
					}
					df.addColumn(vname, o_arr);
					break;
				case BYTE: byte[] b_arr = new byte[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); b_arr[dl] = a.getByte(ind);
					}
					df.addColumn(vname, b_arr);
					break;
				case INT: int[] i_arr = new int[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); i_arr[dl] = a.getInt(ind);
					}
					df.addColumn(vname, i_arr);
					break;
				case SHORT: short[] r_arr = new short[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); r_arr[dl] = a.getShort(ind);
					}
					df.addColumn(vname, r_arr);
					break;
				case LONG: long[] l_arr = new long[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); l_arr[dl] = a.getLong(ind);
					}
					df.addColumn(vname, l_arr);
					break;
				case FLOAT: float[] f_arr = new float[datalength];
					boolean hasFVf = (var.findAttribute("_FillValue")!=null); double ffill = Float.NaN;
					if(hasFVf) { ffill = (float) var.findAttribute("_FillValue").getNumericValue();
						//System.out.println("    [DEBUG] \""+var.getFullName()+"\":_FillValue = "+ffill);
					}
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); f_arr[dl] = a.getFloat(ind);
						if(hasFVf && f_arr[dl]==ffill) f_arr[dl] = Float.NaN;
					}
					df.addColumn(vname, f_arr);
					break;
				case DOUBLE: double[] d_arr = new double[datalength];
					boolean hasFVd = (var.findAttribute("_FillValue")!=null); double dfill = Double.NaN;
					if(hasFVd) { dfill = (double) var.findAttribute("_FillValue").getNumericValue();
						//System.out.println("    [DEBUG] \""+var.getFullName()+"\":_FillValue = "+dfill);
					}
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); d_arr[dl] = a.getDouble(ind);
						if(hasFVd && d_arr[dl]==dfill) d_arr[dl] = Double.NaN;
					}
					df.addColumn(vname, d_arr);
					break;
//				case CHAR:
//					break;
//				case STRING:
//					break;
//				case STRUCTURE:
//					readStructs(var, vname, a, datalength, ind, indexindex, indices, dimids);
//					break;
				default:
					System.err.println("WARNING: could not add variable \""+vname+"\": unsupported data type!");
					break;
			}
		}
		return df;
	}
	public DataFrame2D getVars2D(String... vars) throws IllegalAccessException {
		if(ncfile==null) throw new IllegalAccessException("The Netcdf file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		List<Dimension> dims = new ArrayList<Dimension>();
		List<Variable> ncvars = new ArrayList<Variable>();
		collectVarsAndDims(vars, dims, ncvars);
		if(ncvars.size()==0)
			throw new IllegalArgumentException("Netcdf file does not contain any of given variables.");
		if(dims.size()!=2)
			throw new DimensionMismatchException("The 2D-dataframe can only handle exact 2 dimensions, but found "+dims.size()+".");
		long datasize = 1L;
		for(Dimension d: dims) {
			datasize *= d.getLength();
			if(datasize>MAX_DATAFRAME_LENGTH)
				throw new ArrayStoreException("Cannot expand combination of dimensions, maximum length of dataframe reached!");
		}
		int datalength = 1;
		int dlen = dims.size();
		int[] dimlen = new int[dlen];
		for(int d=0; d<dlen; d++) { dimlen[d] = dims.get(d).getLength(); datalength *= dimlen[d]; }
		int[][] indices = new int[datalength][dlen];
		for(int dl=0; dl<datalength; dl++) {
			int sublen = 1;
			for(int d=dlen-1; d>=0; d--) {
				indices[dl][d] = (dl / sublen) % dimlen[d];
				sublen *= dimlen[d];
			}
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(2, indices);
		DataFrame2D df = new DataFrame2D();
		for(Variable var: ncvars) {
			String vname = var.getFullName();
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
				//ioe.printStackTrace();
				System.err.println("WARNING: could not read variable \""+vname+"\": does not add to the dataframe!");
				continue;
			}
			Dimension[] vardims = var.getDimensions().toArray(new Dimension[0]);
			if(vardims.length==0) {
				switch(var.getDataType()) {
					case BOOLEAN: boolean bool = a.getBoolean(0); boolean[][] bool_arr = new boolean[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) bool_arr[v][u] = bool; 
						df.addColumn(vname, bool_arr); break;
					case BYTE: byte etyb = a.getByte(0); byte[][] byte_arr = new byte[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) byte_arr[v][u] = etyb; 
						df.addColumn(vname, byte_arr); break;
					case SHORT: short trohs = a.getShort(0); short[][] short_arr = new short[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) short_arr[v][u] = trohs; 
						df.addColumn(vname, short_arr); break;
					case INT: int tni = a.getInt(0); int[][] int_arr = new int[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) int_arr[v][u] = tni; 
						df.addColumn(vname, int_arr); break;
					case LONG: long gnol = a.getLong(0); long[][] long_arr = new long[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) long_arr[v][u] = gnol; 
						df.addColumn(vname, long_arr); break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float taolf = a.getFloat(0)==f_fill ? Float.NaN : a.getFloat(0); float[][] float_arr = new float[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) float_arr[v][u] = taolf; 
						df.addColumn(vname, float_arr); break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double elbuod = a.getDouble(0)==d_fill ? Double.NaN : a.getDouble(0); double[][] double_arr = new double[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) double_arr[v][u] = elbuod; 
						df.addColumn(vname, double_arr); break;
//					case CHAR:
//						break;
//					case STRING:
//						break;
//					case STRUCTURE:
//						readStructs(var, var.getFullName(), a, dimlen);
//						break;
					default:
						System.err.println("WARNING: datatype "+var.getDataType().name()+" not supported");
						break;
				}
			} else {
				int uf=0, vf=1; Dimension d = var.getDimension(0);
				if(vardims.length==1) {
					if(d.equals(dims.get(1))) { uf=1; vf=0; }
				} else {
					uf = 1; vf = dimlen[1];
					if(d.equals(dims.get(1))) { uf=dimlen[1]; vf = 1; }
				}
				switch(var.getDataType()) {
					case BOOLEAN: boolean[] arr_bool = (boolean[]) a.get1DJavaArray(ucar.ma2.DataType.BOOLEAN);
						boolean[][] bool_arr = new boolean[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) bool_arr[v][u] = arr_bool[v*vf+u*uf];
						df.addColumn(vname, bool_arr); arr_bool = null; bool_arr = null; break;
					case BYTE: byte[] arr_byte = (byte[]) a.get1DJavaArray(ucar.ma2.DataType.BYTE);
						byte[][] byte_arr = new byte[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) byte_arr[v][u] = arr_byte[v*vf+u*uf];
						df.addColumn(vname, byte_arr); arr_byte = null; byte_arr = null; break;
					case SHORT: short[] arr_short = (short[]) a.get1DJavaArray(ucar.ma2.DataType.SHORT);
						short[][] short_arr = new short[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) short_arr[v][u] = arr_short[v*vf+u*uf];
						df.addColumn(vname, short_arr); arr_short = null; short_arr = null; break;
					case INT: int[] arr_int = (int[]) a.get1DJavaArray(ucar.ma2.DataType.INT);
						int[][] int_arr = new int[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) int_arr[v][u] = arr_int[v*vf+u*uf];
						df.addColumn(vname, int_arr); arr_int = null; int_arr = null; break;
					case LONG: long[] arr_long = (long[]) a.get1DJavaArray(ucar.ma2.DataType.LONG);
						long[][] long_arr = new long[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) long_arr[v][u] = arr_long[v*vf+u*uf];
						df.addColumn(vname, long_arr); arr_long = null; long_arr = null; break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float[] arr_float = (float[]) a.get1DJavaArray(ucar.ma2.DataType.FLOAT);
						for(int ff=0; ff<arr_float.length; ff++) if(arr_float[ff]==f_fill) arr_float[ff] = Float.NaN;
						float[][] float_arr = new float[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) float_arr[v][u] = arr_float[v*vf+u*uf];
						df.addColumn(vname, float_arr); arr_float = null; float_arr = null; break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double[] arr_double = (double[]) a.get1DJavaArray(ucar.ma2.DataType.DOUBLE);
						for(int dd=0; dd<arr_double.length; dd++) if(arr_double[dd]==d_fill) arr_double[dd] = Double.NaN;
						double[][] double_arr = new double[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) double_arr[v][u] = arr_double[v*vf+u*uf];
						df.addColumn(vname, double_arr); arr_double = null; double_arr = null; break;
//					case CHAR:
//						break;
//					case STRING:
//						break;
//					case STRUCTURE:
//						readStructs(var, var.getFullName(), a, dimlen, uf, vf);
//						break;
					default:
						System.err.println("WARNING: datatype "+var.getDataType().name()+" not supported");
						break;
				}
			}
			a = null;
			var = null;
		}
		for(int dimid=0; dimid<2; dimid++) {
			String dname = dims.get(dimid).getName();
			Variable var = ncfile.findVariable(dname);
			if(var==null) {
				System.err.println("Cannot find dimension-variable with name \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
//				ioe.printStackTrace();
				System.err.println("Cannot read variable for dimension \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			double[] dimvalues = new double[dims.get(dimid).getLength()];
			for(int i=0; i<dimvalues.length; i++)
				dimvalues[i] = a.getDouble(i);
			df.setDimension(dimid+Constants.FIRST_IDX, dimvalues, dname);
			Map<String, String> attributes = new HashMap<String,String>();
			for(Attribute att: var.attributes())
				attributes.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue()));
			df.setDimensionsAttributes(dimid+Constants.FIRST_IDX, attributes);
		}
		return df;
	}
	public DataFrame3D getVars3D(String... vars) throws IllegalAccessException {
		if(ncfile==null) throw new IllegalAccessException("The Netcdf file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		List<Dimension> dims = new ArrayList<Dimension>();
		List<Variable> ncvars = new ArrayList<Variable>();
		collectVarsAndDims(vars, dims, ncvars);
		if(ncvars.size()==0)
			throw new IllegalArgumentException("Netcdf file does not contain any of given variables.");
		if(dims.size()!=3)
			throw new DimensionMismatchException("The 2D-dataframe can only handle exact 2 dimensions, but found "+dims.size()+".");
		long datasize = 1L;
		for(Dimension d: dims) {
			datasize *= d.getLength();
			if(datasize>MAX_DATAFRAME_LENGTH)
				throw new ArrayStoreException("Cannot expand combination of dimensions, maximum length of dataframe reached!");
		}
		int datalength = 1;
		int dlen = dims.size();
		int[] dimlen = new int[dlen];
		for(int d=0; d<dlen; d++) { dimlen[d] = dims.get(d).getLength(); datalength *= dimlen[d]; }
		int[][] indices = new int[datalength][dlen];
		for(int dl=0; dl<datalength; dl++) {
			int sublen = 1;
			for(int d=dlen-1; d>=0; d--) {
				indices[dl][d] = (dl / sublen) % dimlen[d];
				sublen *= dimlen[d];
			}
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(2, indices);
		DataFrame3D df = new DataFrame3D();
		for(Variable var: ncvars) {
			String vname = var.getFullName();
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
				//ioe.printStackTrace();
				System.err.println("WARNING: could not read variable \""+vname+"\": does not at to the dataframe!");
				continue;
			}
			Dimension[] vardims = var.getDimensions().toArray(new Dimension[0]);
			if(vardims.length==0) {
				switch(var.getDataType()) {
					case BOOLEAN: boolean bool = a.getBoolean(0); boolean[][][] bool_arr = new boolean[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) bool_arr[w][v][u] = bool;
						df.addColumn(vname, bool_arr); break;
					case BYTE: byte etyb = a.getByte(0); byte[][][] byte_arr = new byte[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) byte_arr[w][v][u] = etyb;
						df.addColumn(vname, byte_arr); break;
					case SHORT: short trohs = a.getShort(0); short[][][] short_arr = new short[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) short_arr[w][v][u] = trohs;
						df.addColumn(vname, short_arr); break;
					case INT: int tni = a.getInt(0); int[][][] int_arr = new int[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) int_arr[w][v][u] = tni;
						df.addColumn(vname, int_arr); break;
					case LONG: long gnol = a.getLong(0); long[][][] long_arr = new long[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) long_arr[w][v][u] = gnol;
						df.addColumn(vname, long_arr); break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float taolf = a.getFloat(0)==f_fill ? Float.NaN : a.getFloat(0); float[][][] float_arr = new float[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) float_arr[w][v][u] = taolf;
						df.addColumn(vname, float_arr); break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double elbuod = a.getDouble(0)==d_fill ? Double.NaN : a.getDouble(0); double[][][] double_arr = new double[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) double_arr[w][v][u] = elbuod;
						df.addColumn(vname, double_arr); break;
//					case CHAR:
//						break;
//					case STRING:
//						break;
//					case STRUCTURE:
//						readStructs(var, var.getFullName(), a, dimlen);
//						break;
					default:
						System.out.println("WARNING: datatype "+var.getDataType().name()+" not supported");
						break;
				}
			} else {
				int uf=0, vf=0, wf=1; Dimension d1 = var.getDimension(0);
				if(vardims.length==1) {
					if(d1.equals(dims.get(1))) { uf=0; vf=1; wf=0; }
					if(d1.equals(dims.get(2))) { uf=1; vf=0; wf=0; }
				}
				if(vardims.length==2) {
					Dimension d2 = var.getDimension(1);
					if(d1.equals(dims.get(0))) {      // w,?
						if(d2.equals(dims.get(1))) { wf=dimlen[1]; vf=1; uf=0; }   // w,v
						if(d2.equals(dims.get(2))) { wf=dimlen[2]; vf=0; uf=1; } } // w,u
					else if(d1.equals(dims.get(1))) { // v,?
						if(d2.equals(dims.get(0))) { wf=1; vf=dimlen[0]; uf=0; }   // v,w
						if(d2.equals(dims.get(2))) { wf=0; vf=dimlen[2]; uf=1; } } // v,u
					else {                           // u,?
						if(d2.equals(dims.get(0))) { wf=1; vf=0; uf=dimlen[0]; }   // u,w
						if(d2.equals(dims.get(1))) { wf=0; vf=1; uf=dimlen[1]; } } // u,v
				}
				if(vardims.length==3) {
					Dimension d2 = var.getDimension(1);
					if(d1.equals(dims.get(0))) {      // w,?,?
						if(d2.equals(dims.get(1))) { wf=dimlen[1]*dimlen[2]; vf=dimlen[2]; uf=1; }   // w,v,u
						if(d2.equals(dims.get(2))) { wf=dimlen[2]*dimlen[1]; vf=1; uf=dimlen[1]; } } // w,u,v
					else if(d1.equals(dims.get(1))) { // v,?,?
						if(d2.equals(dims.get(0))) { wf=dimlen[2]; vf=dimlen[0]*dimlen[2]; uf=1; }   // v,w,u
						if(d2.equals(dims.get(2))) { wf=1; vf=dimlen[2]*dimlen[0]; uf=dimlen[0]; } } // v,u,w
					else {                           // u,?,?
						if(d2.equals(dims.get(0))) { wf=dimlen[1]; vf=1; uf=dimlen[0]*dimlen[1]; }   // u,w,v
						if(d2.equals(dims.get(1))) { wf=1; vf=dimlen[0]; uf=dimlen[1]*dimlen[0]; } } // u,v,w
				}
				switch(var.getDataType()) {
					case BOOLEAN: boolean[] arr_bool = (boolean[]) a.get1DJavaArray(ucar.ma2.DataType.BOOLEAN);
						boolean[][][] bool_arr = new boolean[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							bool_arr[w][v][u] = arr_bool[w*wf+v*vf+u*uf];
						df.addColumn(vname, bool_arr);
						arr_bool = null; bool_arr = null; break;
					case BYTE: byte[] arr_byte = (byte[]) a.get1DJavaArray(ucar.ma2.DataType.BYTE);
						byte[][][] byte_arr = new byte[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							byte_arr[w][v][u] = arr_byte[w*wf+v*vf+u*uf];
						df.addColumn(vname, byte_arr);
						arr_byte = null; byte_arr = null; break;
					case SHORT: short[] arr_short = (short[]) a.get1DJavaArray(ucar.ma2.DataType.SHORT);
						short[][][] short_arr = new short[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							short_arr[w][v][u] = arr_short[w*wf+v*vf+u*uf];
						df.addColumn(vname, short_arr);
						arr_short = null; short_arr = null; break;
					case INT: int[] arr_int = (int[]) a.get1DJavaArray(ucar.ma2.DataType.INT);
						int[][][] int_arr = new int[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							int_arr[w][v][u] = arr_int[w*wf+v*vf+u*uf];
						df.addColumn(vname, int_arr);
						arr_int = null; int_arr = null; break;
					case LONG: long[] arr_long = (long[]) a.get1DJavaArray(ucar.ma2.DataType.LONG);
						long[][][] long_arr = new long[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							long_arr[w][v][u] = arr_long[w*wf+v*vf+u*uf];
						df.addColumn(vname, long_arr);
						arr_long = null; long_arr = null; break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float[] arr_float = (float[]) a.get1DJavaArray(ucar.ma2.DataType.FLOAT);
						for(int ff=0; ff<arr_float.length; ff++) if(arr_float[ff]==f_fill) arr_float[ff] = Float.NaN;
						float[][][] float_arr = new float[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							float_arr[w][v][u] = arr_float[w*wf+v*vf+u*uf];
						df.addColumn(vname, float_arr);
						arr_float = null; float_arr = null; break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double[] arr_double = (double[]) a.get1DJavaArray(ucar.ma2.DataType.DOUBLE);
						for(int dd=0; dd<arr_double.length; dd++) if(arr_double[dd]==d_fill) arr_double[dd] = Double.NaN;
						double[][][] double_arr = new double[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							double_arr[w][v][u] = arr_double[w*wf+v*vf+u*uf];
						df.addColumn(vname, double_arr);
						arr_double = null; double_arr = null; break;
//					case CHAR:
//						break;
//					case STRING:
//						break;
//					case STRUCTURE:
//						readStructs(var, var.getFullName(), a, dimlen, uf, vf, wf);
//						break;
					default:
						System.out.println("WARNING: datatype "+var.getDataType().name()+" not supported");
						break;
				}
			}
		}
		for(int dimid=0; dimid<3; dimid++) {
			String dname = dims.get(dimid).getName();
			Variable var = ncfile.findVariable(dname);
			if(var==null) {
				System.err.println("Cannot find dimension-variable with name \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
//				ioe.printStackTrace();
				System.err.println("Cannot read variable for dimension \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			double[] dimvalues = new double[dimlen[dimid]];
			for(int i=0; i<dimvalues.length; i++)
				dimvalues[i] = a.getDouble(i);
			df.setDimension(dimid+Constants.FIRST_IDX, dimvalues, dname);
			Map<String,String> attributes = new HashMap<String,String>();
			for(Attribute att: var.attributes())
				attributes.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue()));
			df.setDimensionsAttributes(dimid+Constants.FIRST_IDX, attributes);
		}
		return df;
	}
	private void collectVarsAndDims(String[] names, List<Dimension> dims, List<Variable> vars) {
		for(int b=0; b<names.length; b++) {
			Variable var = ncfile.findVariable(names[b]);
			if(var==null) continue;
			if(!vars.contains(var)) vars.add(var);
			for(Dimension d: var.getDimensions()) {
				if(!dims.contains(d)) dims.add(d);
			}
		}
	}
	
	
	protected static boolean canRead(File file) {
		try(BufferedReader br = new BufferedReader(new FileReader(file))) {
			char[] buffer = new char[116];
			int readlength = br.read(buffer);
			if(readlength==0) return false;
			String header = new String(buffer);
			if(header.contains("CDF")) return true;
			if(header.contains("HDF")) return true;
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return false;
	}
	protected static boolean hasValidExtension(String filename) {
		if(filename.endsWith(".nc")) return true;
		return false;
	}
	
	
	@Override
	public void close() throws IOException {
		if(ncfile!=null) ncfile.close();
		ncfile = null;
	}

}
