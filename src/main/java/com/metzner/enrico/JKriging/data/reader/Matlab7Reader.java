package com.metzner.enrico.JKriging.data.reader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableList;
import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;
import com.metzner.enrico.JKriging.error.DimensionMismatchException;

import ucar.ma2.Array;
import ucar.ma2.ArrayString;
import ucar.ma2.Index;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.Group;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;

class Matlab7Reader extends DataReader {
	private final static ucar.ma2.DataType NCSTRING = ucar.ma2.DataType.STRING;
	private final static ucar.ma2.DataType NCSTRUCT = ucar.ma2.DataType.STRUCTURE;
	
	private NetcdfFile ncfile = null;
	
	Matlab7Reader(File file) throws IOException {
		super(file);
		ncfile = NetcdfFiles.open(file.getAbsolutePath());
	}
	
	@Override
	public void describeContent() {
		System.out.println(ncfile.getFileTypeDescription());
		System.out.println(wT+wL+"Path: "+ncfile.getLocation());
		System.out.println(wT+wL+"Name: "+ncfile.getTitle());
		describeGroup(ncfile.getRootGroup(),wT,wP);
		describeAttributes(ncfile.getGlobalAttributes(),wE," ");
	}
	
	@Override
	public Map<String, Integer> getDimensionsOfVariable(String var) throws IllegalAccessException {
		if(ncfile==null) throw new IllegalAccessException("The Netcdf file does not exist.");
		if(var==null) throw new IllegalAccessException("Variable cannot be null");
		Variable v = Variable.of(ncfile.findVariable(var));
//		if(v==null) throw new IllegalAccessException("Cannot find variable <"+var+"> in the Netcdf file.");
		if(v==null) return null;
		Map<String, Integer> dimmap = new HashMap<>();
		for(Dimension d: v.getDimensions())
			dimmap.put(d.getName(), d.getLength());
		return dimmap;
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
			Variable var = Variable.of(ncfile.findVariable(dim.getName()));
			if(var==null) {
				int[] arr = new int[datalength];
				for(int dl=0; dl<datalength; dl++) arr[dl] = indices[dl][dimid];
				df.addColumn(dim.getName(), arr);
			} else {
				Array a = null;
				try {
					a = var.read();
				} catch (IOException|InvalidRangeException e) {
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
			} catch(IOException|InvalidRangeException e) {
				//ioe.printStackTrace();
				System.err.println("WARNING: could not read variable \""+vname+"\": does not at to the dataframe!");
				continue;
			}
			Index ind = a.getIndex();
			List<Dimension> vdims = var.getDimensions();
			int[] dimids = new int[vdims.size()];
			for(int d=0; d<dimids.length; d++) {
				dimids[d] = dims.indexOf(vdims.get(d));
				dimids[d] = 0;
				int dd = -1;
				for(dd=0; dd<dlen; dd++) if(dims.get(dd).equals(vdims.get(d))) {
					dimids[d] = dd;
					break;
				}
				if(dd==dlen) System.err.println("WARNING: could not find dimension \""+vdims.get(d).getName()+"\": maybe read data in wrong order!");
			}
			int[] indexindex = new int[vdims.size()];
			switch(var.getDataType()) {
				case BOOLEAN: boolean[] o_arr = new boolean[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); o_arr[dl] = a.getBoolean(ind);
						} else { o_arr[dl] = a.getBoolean(0); }
					}
					df.addColumn(vname, o_arr);
					break;
				case BYTE: byte[] b_arr = new byte[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); b_arr[dl] = a.getByte(ind);
						} else { b_arr[dl] = a.getByte(0); }
					}
					df.addColumn(vname, b_arr);
					break;
				case INT: int[] i_arr = new int[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); i_arr[dl] = a.getInt(ind);
						} else { i_arr[dl] = a.getInt(0); }
					}
					df.addColumn(vname, i_arr);
					break;
				case SHORT: short[] r_arr = new short[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); r_arr[dl] = a.getShort(ind);
						} else { r_arr[dl] = a.getShort(0); }
					}
					df.addColumn(vname, r_arr);
					break;
				case LONG: long[] l_arr = new long[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); l_arr[dl] = a.getLong(ind);
						} else { l_arr[dl] = a.getLong(0); }
					}
					df.addColumn(vname, l_arr);
					break;
				case FLOAT: float[] f_arr = new float[datalength];
					boolean hasFVf = (var.findAttribute("_FillValue")!=null); double ffill = Float.NaN;
					if(hasFVf) { ffill = (float) var.findAttribute("_FillValue").getNumericValue();
						//System.out.println("    [DEBUG] \""+var.getFullName()+"\":_FillValue = "+ffill);
					}
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); f_arr[dl] = a.getFloat(ind);
						} else { f_arr[dl] = a.getFloat(0); }
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
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); d_arr[dl] = a.getDouble(ind);
						} else { d_arr[dl] = a.getDouble(0); }
						if(hasFVd && d_arr[dl]==dfill) d_arr[dl] = Double.NaN;
					}
					df.addColumn(vname, d_arr);
					break;
//				case CHAR:
//					break;
				case STRING: String[] g_arr = new String[datalength];
					for(int dl=0; dl<datalength; dl++) {
						if(!vdims.get(0).equals(Variable.SCALER_DIM)) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); g_arr[dl] = (String) a.getObject(ind);
						} else { g_arr[dl] = (String) a.getObject(0); }
					}
					df.addColumn(vname, g_arr);
					break;
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
			} catch(IOException|InvalidRangeException e) {
//				e.printStackTrace();
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
			Variable var = Variable.of(ncfile.findVariable(dname));
			if(var==null) {
				System.err.println("Cannot find dimension-variable with name \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			Array a = null;
			try {
				a = var.read();
			} catch(IOException|InvalidRangeException e) {
//				e.printStackTrace();
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
			} catch(IOException|InvalidRangeException e) {
				//e.printStackTrace();
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
			Variable var = Variable.of(ncfile.findVariable(dname));
			if(var==null) {
				System.err.println("Cannot find dimension-variable with name \""+dname+"\", use index-array instead.");
				df.setDimension(dimid+Constants.FIRST_IDX, dname);
				continue;
			}
			Array a = null;
			try {
				a = var.read();
			} catch(IOException|InvalidRangeException e) {
//				e.printStackTrace();
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
			Variable var = Variable.of(ncfile.findVariable(names[b]));
			if(var==null) continue;
			if(!vars.contains(var)) vars.add(var);
			for(Dimension d: var.getDimensions()) {
				if(!dims.contains(d)) dims.add(d);
			}
		}
	}
	
	public DataFrame get1Dslice(String filtervar, int index, String... vars) throws IllegalAccessException {
		if(ncfile==null) throw new IllegalAccessException("The Netcdf file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		List<Dimension> dims = new ArrayList<Dimension>();
		List<Variable> ncvars = new ArrayList<Variable>();
		collectVarsAndDims(vars, dims, ncvars);
		if(ncvars.size()==0)
			throw new IllegalArgumentException("Netcdf file does not contain any of given variables.");
		if(dims.size()!=2)
			throw new DimensionMismatchException("Can only take 1d slice from collection with exactly 2 dimensions.");
		int seldim = -1;
		if(dims.get(0).getName().equals(filtervar)) seldim = 0;
		if(dims.get(1).getName().equals(filtervar)) seldim = 1;
		if(seldim<0) throw new IllegalAccessException("Cannot find dimension "+filtervar+" for choosen variables.");
		String dname   = dims.get(1-seldim).getName();
		int datalength = dims.get(1-seldim).getLength();
		int[][] indices = new int[datalength][2];
		for(int dl=0; dl<datalength; dl++) {
			indices[dl][0] = dl;
			indices[dl][1] = dl;
			indices[dl][seldim] = index;
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(2, indices);
		DataFrame df = new DataFrame();
		for(Variable var: ncvars) {
			String vname = var.getFullName();
			Array a = null;
			Dimension[] vardims = var.getDimensions().toArray(new Dimension[0]);
			int effdimcnt = -1;
			if(vardims.length==0) effdimcnt = 0;
			if(vardims.length==1) effdimcnt = vardims[0].getName().equals(filtervar) ? 0 : 1;
			if(vardims.length==2) effdimcnt = 2;
			if(effdimcnt<0) {
				System.err.println("WARNING: something went wrong while reading variable "+vname+".");
				continue;
			}
			if(effdimcnt==0) {
				try {
					a = var.read();
				} catch(IOException|InvalidRangeException e) {
//					e.printStackTrace();
					System.err.println("WARNING: could not read variable \""+vname+"\": does not add to the dataframe!");
					continue;
				}
				int selidx = vardims.length==0 ? 0 : index;
				switch(var.getDataType()) {
					case BOOLEAN: df.addColumn(vname, (boolean[])constArray('o',datalength, a.getBoolean(selidx))); break;
					case BYTE: df.addColumn(vname, (byte[])constArray('b',datalength, a.getByte(selidx))); break;
					case SHORT: df.addColumn(vname, (short[])constArray('r',datalength, a.getShort(selidx))); break;
					case INT: df.addColumn(vname, (int[])constArray('i',datalength, a.getInt(selidx))); break;
					case LONG: df.addColumn(vname, (boolean[])constArray('l',datalength, a.getLong(selidx))); break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						df.addColumn(vname, (float[])constArray('f',datalength,f_fill==a.getFloat(selidx)?Float.NaN:a.getFloat(0))); break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						df.addColumn(vname, (double[])constArray('d',datalength,d_fill==a.getDouble(selidx)?Double.NaN:a.getDouble(0))); break;
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
				try {
					if(effdimcnt==1) a = var.read();
					if(effdimcnt==2) {
						int[] st={0,0},ln={vardims[0].getLength(),vardims[1].getLength()};
						if(dims.get(0).getName().equals(filtervar)) { st[0]=index; ln[0]=1; }
						if(dims.get(1).getName().equals(filtervar)) { st[1]=index; ln[1]=1; }
						a = var.read(st,ln).reduce();
					}
				} catch(IOException|InvalidRangeException e) {
					//e.printStackTrace();
					System.err.println("WARNING: could not read variable \""+vname+"\": does not add to the dataframe!");
					continue;
				}
				switch(var.getDataType()) {
					case BOOLEAN: df.addColumn(vname, (boolean[]) a.get1DJavaArray(ucar.ma2.DataType.BOOLEAN)); break;
					case BYTE:df.addColumn(vname, (byte[]) a.get1DJavaArray(ucar.ma2.DataType.BYTE)); break;
					case SHORT: df.addColumn(vname, (short[]) a.get1DJavaArray(ucar.ma2.DataType.SHORT)); break;
					case INT: df.addColumn(vname, (int[]) a.get1DJavaArray(ucar.ma2.DataType.INT)); break;
					case LONG: df.addColumn(vname, (long[]) a.get1DJavaArray(ucar.ma2.DataType.LONG)); break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float[] arr_float = (float[]) a.get1DJavaArray(ucar.ma2.DataType.FLOAT);
						for(int ff=0; ff<arr_float.length; ff++) if(arr_float[ff]==f_fill) arr_float[ff] = Float.NaN;
						df.addColumn(vname, arr_float); break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double[] arr_double = (double[]) a.get1DJavaArray(ucar.ma2.DataType.DOUBLE);
						for(int dd=0; dd<arr_double.length; dd++) if(arr_double[dd]==d_fill) arr_double[dd] = Double.NaN;
						df.addColumn(vname, arr_double); break;
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
		
		Variable var = Variable.of(ncfile.findVariable(dname));
		while(true) {
			if(var==null) {
				System.err.println("Cannot find dimension-variable with name \""+dname+"\", use index-array instead.");
				df.setDimension(dname);
				break;
			}
			if(var.getDimensions().size()!=1) {
				System.err.println("Dimension-variable is multidimensional");
				df.setDimension(dname);
				break;
			}
			Array a = null;
			try {
				a = var.read();
			} catch(IOException|InvalidRangeException e) {
//				e.printStackTrace();
				System.err.println("Cannot read variable for dimension \""+dname+"\", use index-array instead.");
				df.setDimension(dname);
				break;
			}
			System.out.println(a);
			double[] dimvalues = new double[datalength];
			for(int i=0; i<dimvalues.length; i++)
				dimvalues[i] = a.getDouble(i);
			df.setDimension(dimvalues, dname);
			Map<String, String> attributes = new HashMap<String,String>();
			for(Attribute att: var.attributes())
				attributes.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue()));
			df.setDimensionsAttributes(attributes);
			break;
		}
		return df;
	}
	private Object constArray(char type, int len, Object value) {
		switch(type) {
			case 'o': boolean[] ores = new boolean[len]; for(int o=0; o<len; o++) ores[o] = (boolean)value; return ores;
			case 'b': byte[]    bres = new byte[len];    for(int o=0; o<len; o++) bres[o] = (byte)value;    return bres;
			case 'r': short[]   rres = new short[len];   for(int o=0; o<len; o++) rres[o] = (short)value;   return rres;
			case 'i': int[]     ires = new int[len];     for(int o=0; o<len; o++) ires[o] = (int)value;     return ires;
			case 'l': long[]    lres = new long[len];    for(int o=0; o<len; o++) lres[o] = (long)value;    return lres;
			case 'f': float[]   fres = new float[len];   for(int o=0; o<len; o++) fres[o] = (float)value;   return fres;
			case 'd': double[]  dres = new double[len];  for(int o=0; o<len; o++) dres[o] = (double)value;  return dres;
			case 's': String[]  sres = new String[len];  for(int o=0; o<len; o++) sres[o] = (String)value;  return sres;
			default: return null;
		}
	}
	
	protected static boolean canRead(File file) {
		try(BufferedReader br = new BufferedReader(new FileReader(file))) {
			char[] buffer = new char[116];
			int readlength = br.read(buffer);
			if(readlength==0) return false;
			String header = new String(buffer);
//			System.out.println("[NetcdfReader] fileheader: "+header);
			if(header.contains("CDF")) return true;
			if(header.contains("HDF")) return true;
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		System.err.println("WARNING: cannot determin format from header of file <"+file.getAbsolutePath()+">");
		try(NetcdfFile nc = NetcdfFiles.open(file.getAbsolutePath())) {
			return true;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
	}
	protected static boolean hasValidExtension(String filename) {
		if(filename.endsWith(".nc")) return true;
		if(filename.endsWith(".hdf")) return true;
		if(filename.endsWith(".grib")) return true;
		return false;
	}
	
	
	@Override
	public void close() throws IOException {
		if(ncfile!=null) ncfile.close();
		ncfile = null;
	}
	
	@SuppressWarnings("deprecation")
	private void describeGroup(Group group,String pad1,String pad2) {
		if(pad1.equals(wT))
			System.out.println(pad1+wL+"Rootgroup ("+group.getFullName()+"):");
		else
			System.out.println(pad1+wL+group.getFullName()+":");
//		System.out.println(pad2+"  "+wT+wL+"Dimensions:");
		describeDimensions(group.getDimensions(), pad2+"  ");
		describeVariables(group.getVariables(), pad2+"  ");
		List<Group> groups = group.getGroups();
		describeAttributes(group.getAttributes(), pad2+"  "+(groups.isEmpty()?wE:wT), pad2+"  "+(groups.isEmpty()?" ":wP));
		for(int g=0; g<groups.size(); g++)
			describeGroup(groups.get(g), pad2+"  "+(g+1==groups.size()?wE:wT), pad2+"  "+(g+1==groups.size()?" ":wP));
	}
	private void describeDimensions(List<Dimension> dims, String pad) {
		System.out.println(pad+wT+wL+"Dimensions:");
		for(int i=0; i<dims.size(); i++) {
			String p = pad+wP+"  "+(i+1==dims.size()?wE:wT)+wL;
			System.out.println(p+dims.get(i).getName()+" ("+dims.get(i).getLength()+")");
		}
	}
	private void describeVariables(List<ucar.nc2.Variable> vars, String pad) {
		System.out.println(pad+wT+wL+"Variables:");
		for(int i=0; i<vars.size(); i++) {
			Variable v = Variable.of(vars.get(i));
			String p = pad+wP+"  "+(i+1==vars.size()?wE:wT)+wL;
			String t = String.format("%-6s", v.getDataType().name().toLowerCase());
			String n = v.getNameAndDimensions();
			System.out.println(p+t+" "+n);
		}
	}
	private void describeAttributes(List<Attribute> attribs, String pad1, String pad2) {
		boolean global = pad1.equals(wE);
		System.out.println(pad1+wL+(global?"Gloabl a":"A")+"ttributes");
		for(int a=0; a<attribs.size(); a++) {
			Attribute att = attribs.get(a);
			System.out.println(pad2+" "+(a+1==attribs.size()?wE:wT)+wL+att.getName()+": "+att.getStringValue());
		}
	}
	
	@SuppressWarnings("unused")
	private static class Variable {
		public final static Dimension SCALER_DIM = new Dimension("ScalarDimension",1);
		private static Map<String,ucar.ma2.DataType> types;
		static {
			types = new HashMap<>();
			types.put("byte", ucar.ma2.DataType.BYTE); types.put("ubyte", ucar.ma2.DataType.UBYTE);
			types.put("int", ucar.ma2.DataType.INT); types.put("uint", ucar.ma2.DataType.UINT);
			types.put("short", ucar.ma2.DataType.SHORT); types.put("ushort", ucar.ma2.DataType.USHORT);
			types.put("long", ucar.ma2.DataType.LONG); types.put("ulong", ucar.ma2.DataType.ULONG);
			types.put("boolean", ucar.ma2.DataType.BOOLEAN);
			types.put("float", ucar.ma2.DataType.FLOAT); types.put("double", ucar.ma2.DataType.DOUBLE);
			types.put("string", ucar.ma2.DataType.STRING);
		};
		public static Variable of(ucar.nc2.Variable v) {
			if(v==null) return null;
			return new Variable(v);
		}
		
		private ucar.nc2.Variable original;
		private ucar.ma2.DataType original_datatype, matlab_datatype;
		private int[] file_dims, output_dims;
		private int[][] dim_mapping;
		
		private Variable(ucar.nc2.Variable v) {
			original = v;
			original_datatype = original.getDataType();
			matlab_datatype = original.getDataType();
			file_dims = original.getShape();
			String m_dt = original.findAttributeString("MATLAB_class", "null");
			if(m_dt.equals("char")) matlab_datatype = ucar.ma2.DataType.STRING;
			else for(String type: types.keySet())
				if(m_dt.equalsIgnoreCase(type))
					matlab_datatype = types.get(type);
//			String fd = ""; for(int i=0; i<file_dims.length; i++) fd+=(i==0?"":",")+file_dims[i];
			if(matlab_datatype==ucar.ma2.DataType.STRING) {
				int c = -1, ld = -1;
				for(int i=0; i<file_dims.length; i++)
					if(file_dims[i]>1) { c++; ld=i; }
					else if(ld<0) ld = i;
//				System.err.println("[DEBUG] "+v.getFullName()+": types="+original_datatype.name()+"/"+matlab_datatype.name()+
//						" fd=["+fd+"] c="+c+" ld="+ld);
				if(c<=0) {
					output_dims = new int[0];
					dim_mapping = new int[1][file_dims.length];
					for(int i=0; i<file_dims.length; i++) dim_mapping[0][i] = i==ld ? 1 : 0;
				} else {
					output_dims = new int[c];
					dim_mapping = new int[c][file_dims.length];
					c = 0;
					for(int j=0; j<file_dims.length; j++) {
						if(file_dims[j]==1) continue;
						output_dims[c] = file_dims[j];
						for(int i=0; i<file_dims.length; i++)
							dim_mapping[c][i] = j==i ? 1 : 0;
						c++;
					}
					dim_mapping[c-1][ld] = 1;
				}
//				output_dims = new int[file_dims.length-1];
//				dim_mapping = new int[output_dims.length][file_dims.length];
//				for(int j=0; j<output_dims.length; j++) {
//					output_dims[j] = file_dims[j];
//					for(int i=0; i<file_dims.length; i++)
//						dim_mapping[j][i] = j==i ? 1 : i==file_dims.length-1 ? 1 : 0;
//				}
			} else {
				int c = 0;
				for(int i: file_dims) if(i>1) c++;
				output_dims = new int[c];
				if(c==0) {
					dim_mapping = new int[1][file_dims.length];
					for(int i=0; i<file_dims.length; i++)
						dim_mapping[0][i] = i==0 ? 1 : 0;
				} else {
					dim_mapping = new int[c][file_dims.length];
					c = 0;
					for(int j=0; j<file_dims.length; j++) {
						if(file_dims[j]==1) continue;
						output_dims[c] = file_dims[j];
						for(int i=0; i<file_dims.length; i++)
							dim_mapping[c][i] = j==i ? 1 : 0;
						c++;
					}
				}
			}
		}
		
		public ucar.ma2.DataType getNetcdfDataType() {
			return original_datatype;
		}
		public ucar.ma2.DataType getDataType() {
			return matlab_datatype;
		}
		public String getFullName() {
			return original.getFullName();
		}
		public String getNameAndDimensions() {
			String res = original.getFullName()+"(";
			for(int d=0; d<output_dims.length; d++) res+=(d==0?"":",")+output_dims[d];
			if(output_dims.length==0) res+="S";
			return res+")";
		}
		public int getShape(int idx) {
			return output_dims[idx];
		}
		public int[] getShape() {
			return output_dims;
		}
		public Dimension getDimension(int idx) {
			if(0<=idx && idx<output_dims.length)
				return new Dimension("generic"+output_dims[idx], output_dims[idx]);
			return null;
		}
		public ImmutableList<Dimension> getDimensions() {
			Dimension[] dims = new Dimension[Math.max(1,output_dims.length)];
			if(output_dims.length==0) dims[0] = SCALER_DIM;
			else for(int i=0; i<dims.length; i++) dims[i] = getDimension(i);
			return ImmutableList.copyOf(dims);
		}
		public boolean hasAttribute(String attName) {
			return original.hasAttribute(attName);
		}
		public Attribute findAttribute(String name) {
			return original.findAttribute(name);
		}
		@SuppressWarnings("deprecation")
		public Array read() throws IOException,InvalidRangeException {
			if(matlab_datatype!=NCSTRING && original_datatype!=NCSTRUCT) {
				//normal data handling
				//implicit: matlab_datatype==original_datatype
				Array array = original.read();
				if(output_dims.length<file_dims.length)
					array = array.reduce();
				return array;
			} else
			if(matlab_datatype!=NCSTRING) {
				//handling structures
			} else
			{
				//String arrays are allways reduced, with last dimension larger than one as the string length
				Array raw = original.read().reduce();
				ArrayString strArr = new ArrayString(output_dims);
				Index index = strArr.getIndex();
				int[] start = raw.getShape().clone();
				int[] range = raw.getShape().clone();
				start[start.length-1] = 0;
				range[range.length-1] = raw.getShape()[range.length-1];
				for(int i=0; i<index.getSize(); i++) {
					int[] pos = index.getCurrentCounter();
					for(int p=0; p<pos.length; p++) {
						start[p] = pos[p];
						range[p] = 1;
					}
					short[] datas = (short[]) raw.section(start, range).get1DJavaArray(short.class);
					byte[] datab = new byte[datas.length*2];
					for(int b=0,s=0; s<datas.length; b=++s) {
						datab[b] = (byte)(datas[s]&255); datab[b+1] = (byte)((datas[s]>>8)&255);
						if(datas[s]==0) { datab[b]=0; datab[b+1] = ' '; } }
					Charset charset = Charset.defaultCharset();
					Attribute att_decode = original.findAttribute("MATLAB_int_decode");
					if(att_decode!=null) switch(att_decode.getNumericValue().intValue()) {
						case 2: charset = StandardCharsets.US_ASCII; break;
						case 3: charset = StandardCharsets.UTF_8; break; //TODO assumption MATLAB_int_decode=3 -> UTF-8 ???
						default: break;
					}
					strArr.set(index, new String(datab, charset).trim());
					index.incr();
				}
				raw = null;
				return strArr;
			}
			//TODO implement read()
			return null;
		}
		public Array read(int[] origin, int[] shape) throws IOException, InvalidRangeException {
			int datastate = 0;
			if(original_datatype==matlab_datatype) datastate |= 4;
			if(output_dims.length==file_dims.length) {
				boolean e = true;
				for(int i=0; i<file_dims.length && e; i++)
					e = output_dims[i]==file_dims[i];
				datastate |= e ? 3 : 2;
			}
			if(matlab_datatype!=NCSTRING) {
				switch(datastate) {
					case 1: return null;
					case 2: return null;
					case 3: return null;
					case 4: return null;
					case 5: return null;
					case 6: return null;
					case 7: return original.read();
				}
			}
			//TODO implement read(origin,shape)
			return null;
		}
		public ucar.nc2.AttributeContainer attributes() {
			return original.attributes();
		}
	}
}
