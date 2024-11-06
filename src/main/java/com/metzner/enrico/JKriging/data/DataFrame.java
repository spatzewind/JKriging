package com.metzner.enrico.JKriging.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.metzner.enrico.JKriging.helper.DataHelper;
import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.LogicHelper;
import com.metzner.enrico.JKriging.probability.StdAnalysis;

import ucar.ma2.Array;
import ucar.ma2.ArrayBoolean;
import ucar.ma2.ArrayByte;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.ma2.ArrayLong;
import ucar.ma2.ArrayShort;
import ucar.ma2.ArrayString;
import ucar.ma2.Index;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;
import ucar.nc2.write.NetcdfFileFormat;
import ucar.nc2.write.NetcdfFormatWriter;

//@SuppressWarnings("deprecation")
public class DataFrame {

	protected int datalength;
	protected DataType default_data_type;
	protected String[] titles;
	protected DataType[] types;
	protected double[] dimension;
	private double[][] minmax_mean_sill;
	private String dimension_name;
	private Map<String, boolean[]> bool_column;
	private Map<String, byte[]>    byte_column;
	private Map<String, int[]>     int_column;
	private Map<String, short[]>   short_column;
	private Map<String, long[]>    long_column;
	private Map<String, float[]>   float_column;
	private Map<String, double[]>  double_column;
	private Map<String, String[]>  string_column;
	protected Map<String, Object>  struct_column;
	private Map<String, String>    attribs_dim;



	public DataFrame() {
		datalength = 0;
		default_data_type = DataType.STRING;
		titles = new String[0];
		types = new DataType[0];
		minmax_mean_sill = new double[4][0];
		bool_column =   new HashMap<>();
		byte_column =   new HashMap<>();
		int_column =    new HashMap<>();
		short_column =  new HashMap<>();
		long_column =   new HashMap<>();
		float_column =  new HashMap<>();
		double_column = new HashMap<>();
		string_column = new HashMap<>();
		struct_column = new HashMap<>();
		dimension = new double[0];
		dimension_name = "";
		attribs_dim = new HashMap<>();
	}
	



	public void addColumn(String _column_name, boolean[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, boolean[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.BOOL);
		add_stdana(Double.NaN,Double.NaN,Double.NaN,Double.NaN);
		boolean[] new_data = new boolean[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = false; }
		bool_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, byte[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, byte[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.BYTE);
		byte[] inax = StdAnalysis.minmax(_column_data); byte[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		byte[] new_data = new byte[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Byte.MIN_VALUE; }
		byte_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, int[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, int[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.INT);
		int[] inax = StdAnalysis.minmax(_column_data); int[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		int[] new_data = new int[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Integer.MIN_VALUE; }
		int_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, short[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, short[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.SHORT);
		short[] inax = StdAnalysis.minmax(_column_data); short[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		short[] new_data = new short[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Short.MIN_VALUE; }
		short_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, long[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, long[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.LONG);
		long[] inax = StdAnalysis.minmax(_column_data); long[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		long[] new_data = new long[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Long.MIN_VALUE; }
		long_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, float[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, float[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.FLOAT);
		float[] inax = StdAnalysis.minmax(_column_data); float[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		float[] new_data = new float[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Float.NaN; }
		float_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, double[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, double[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.DOUBLE);
		double[] inax = StdAnalysis.minmax(_column_data); double[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		double[] new_data = new double[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Double.NaN; }
		double_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void addColumn(String _column_name, String[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, String[] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(hasVariable(_column_name_c)) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name_c, DataType.STRING);
		add_stdana(Double.NaN,-Double.NaN,Double.NaN,Double.NaN);
		String[] new_data = new String[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = ""; }
		string_column.put(_column_name_c, new_data);
		if(dimension.length==0) createDimension();
	}
	public void setDimension(String dim_name) {
		setDimension(null, dim_name);
	}
	public void setDimension(double[] dim_values) {
		setDimension(dim_values, null);
	}
	public void setDimension(double[] dim_values, String dim_name) {
		String _dim_name_c = FormatHelper.name2CFConvention(dim_name);
		boolean change_name = (_dim_name_c!=null);
		boolean change_values = (dim_values!=null);
		if(change_name)
			if(hasVariable(_dim_name_c)) {
				System.err.println("Cannot rename dimension, an other variable has the same name!");
				DataHelper.printStackTrace(System.err);
				return;
			}
		if(change_values) {
			if(dim_values.length==0) {
				System.err.println("Length of dimension has to be greater than 0! Not a valid DataFrame.");
				DataHelper.printStackTrace(System.err);
				return;
			}
			if(dim_values.length!=datalength && datalength>0) {
				System.err.println("Number of values must match corresponding dimension length:\n    expect "+datalength+" but got "+dim_values.length);
				DataHelper.printStackTrace(System.err);
				return;
			}
		}
		if(change_name)
			dimension_name = _dim_name_c;
		if(change_values) {
			if(datalength==0) {
				datalength = dim_values.length;
				dimension = new double[dim_values.length];
			} else
			if(dimension.length==0) {
				dimension = new double[datalength];
			}
			for(int dn=0; dn<datalength; dn++)
				dimension[dn] = dim_values[dn];
		}
	}
	public void addAttributeToDimension(String name, String value) {
		if(attribs_dim.containsKey(name))
			System.out.println("[WARNING] override key \""+name+"\" for dimension");
		attribs_dim.put(name, value);
	}
	public void setDimensionsAttributes(Map<String, String> attributes) {
		attribs_dim.clear(); attribs_dim.putAll(attributes);
	}
	public void copyDimensionsFrom(DataFrame otherDF) {
		this.setDimension(otherDF.getDimensionValues(), otherDF.getDimensionName());
		this.setDimensionsAttributes(otherDF.getDimensionAttributes());
	}
	
	public void renameVariable(int _old_var_id, String _new_name) {
		String _new_name_c = FormatHelper.name2CFConvention(_new_name);
		int ovi = _old_var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not rename variable with ID "+_old_var_id+", could not find variable!");
			return;
		}
		String _old_var_name = getVarname(_old_var_id);
		if(hasVariable(_new_name_c)) {
			System.err.println("Could not rename variable \""+_old_var_name+"\", a variable by name \""+_new_name_c+"\" already exist!");
			//DataHelper.printStackTrace(System.err);
			return;
		}
		switch(types[ovi]) {
			case BOOL:   addColumn(_new_name_c, bool_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case BYTE:   addColumn(_new_name_c, byte_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case SHORT:  addColumn(_new_name_c, short_column.get(_old_var_name));  removeVariable(_old_var_name); break;
			case INT:    addColumn(_new_name_c, int_column.get(_old_var_name));    removeVariable(_old_var_name); break;
			case LONG:   addColumn(_new_name_c, long_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case FLOAT:  addColumn(_new_name_c, float_column.get(_old_var_name));  removeVariable(_old_var_name); break;
			case DOUBLE: addColumn(_new_name_c, double_column.get(_old_var_name)); removeVariable(_old_var_name); break;
			case STRING: addColumn(_new_name_c, string_column.get(_old_var_name)); removeVariable(_old_var_name); break;
			default: System.err.println("An unexpected error occured: Unknown Variable type!"); DataHelper.printStackTrace(System.err); break;
		}
	}
	public void renameVariable(String _old_name, String _new_name) {
		renameVariable(getVariableID(_old_name), _new_name);
	}
	public void setVariableContent(int _var_id, boolean[] new_bools) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_bools.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case BOOL: boolean[] arro = bool_column.get(vn);
				for(int x=0; x<datalength; x++) arro[x] = new_bools[x];
				break;
			case BYTE: byte[] arrb = byte_column.get(vn);
				for(int x=0; x<datalength; x++) arrb[x] = new_bools[x]?(byte)1:(byte)0;
				break;
			case SHORT: short[] arrs = short_column.get(vn);
				for(int x=0; x<datalength; x++) arrs[x] = new_bools[x]?(short)1:(short)0;
				break;
			case INT: int[] arri = int_column.get(vn);
				for(int x=0; x<datalength; x++) arri[x] = new_bools[x]?1:0;
				break;
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = new_bools[x]?1L:0L;
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_bools[x]?1f:0f;
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_bools[x]?1d:0d;
				break;
			case STRING: String[] arrt = string_column.get(vn);
				for(int x=0; x<datalength; x++) arrt[x] = new_bools[x]?"true":"false";
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from bool-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, boolean[] new_bools) {
		setVariableContent(getVariableID(_var_name), new_bools);
	}
	public void setVariableContent(int _var_id, byte[] new_bytes) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_bytes.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case BYTE: byte[] arrb = byte_column.get(vn);
				for(int x=0; x<datalength; x++) arrb[x] = new_bytes[x];
				break;
			case SHORT: short[] arrs = short_column.get(vn);
				for(int x=0; x<datalength; x++) arrs[x] = new_bytes[x];
				break;
			case INT: int[] arri = int_column.get(vn);
				for(int x=0; x<datalength; x++) arri[x] = new_bytes[x];
				break;
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = new_bytes[x];
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_bytes[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_bytes[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from byte-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, byte[] new_bytes) {
		setVariableContent(getVariableID(_var_name), new_bytes);
	}
	public void setVariableContent(int _var_id, short[] new_shorts) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_shorts.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case SHORT: short[] arrs = short_column.get(vn);
				for(int x=0; x<datalength; x++) arrs[x] = new_shorts[x];
				break;
			case INT: int[] arri = int_column.get(vn);
				for(int x=0; x<datalength; x++) arri[x] = new_shorts[x];
				break;
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = new_shorts[x];
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_shorts[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_shorts[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from short-array.");
				break;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, short[] new_shorts) {
		setVariableContent(getVariableID(_var_name), new_shorts);
	}
	public void setVariableContent(int _var_id, int[] new_ints) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_ints.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case INT: int[] arri = int_column.get(vn);
				for(int x=0; x<datalength; x++) arri[x] = new_ints[x];
				break;
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = new_ints[x];
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_ints[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_ints[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from int-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, int[] new_ints) {
		setVariableContent(getVariableID(_var_name), new_ints);
	}
	public void setVariableContent(int _var_id, long[] new_longs) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_longs.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = new_longs[x];
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_longs[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_longs[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from long-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, long[] new_longs) {
		setVariableContent(getVariableID(_var_name), new_longs);
	}
	public void setVariableContent(int _var_id, float[] new_floats) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_floats.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = new_floats[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_floats[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from float-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, float[] new_floats) {
		setVariableContent(getVariableID(_var_name), new_floats);
	}
	public void setVariableContent(int _var_id, double[] new_doubles) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_doubles.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = (float) new_doubles[x];
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = new_doubles[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from double-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, double[] new_doubles) {
		setVariableContent(getVariableID(_var_name), new_doubles);
	}
	public void setVariableContent(int _var_id, String[] new_strings) {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		if(new_strings.length!=datalength)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+datalength+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case BOOL: boolean[] arro = bool_column.get(vn);
				for(int x=0; x<datalength; x++) arro[x] = set_boolean(new_strings[x]);
				break;
			case BYTE: byte[] arrb = byte_column.get(vn);
				for(int x=0; x<datalength; x++) arrb[x] = set_byte(new_strings[x]);
				break;
			case SHORT: short[] arrs = short_column.get(vn);
				for(int x=0; x<datalength; x++) arrs[x] = set_short(new_strings[x]);
				break;
			case INT: int[] arri = int_column.get(vn);
				for(int x=0; x<datalength; x++) arri[x] = set_int(new_strings[x]);
				break;
			case LONG: long[] arrl = long_column.get(vn);
				for(int x=0; x<datalength; x++) arrl[x] = set_long(new_strings[x]);
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				for(int x=0; x<datalength; x++) arrf[x] = set_float(new_strings[x]);
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				for(int x=0; x<datalength; x++) arrd[x] = set_double(new_strings[x]);
				break;
			case STRING: String[] arrt = string_column.get(vn);
				for(int x=0; x<datalength; x++) arrt[x] = new_strings[x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from string-array.");
				return;
		}
		recalcStats(_var_id);
	}
	public void setVariableContent(String _var_name, String[] new_strings) {
		setVariableContent(getVariableID(_var_name), new_strings);
	}
	
	public void removeVariable(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) {
			System.err.println("There exist no variable with ID "+_var_id);
			return;
		}
		removeVariable(getVarname(_var_id));
	}
	public void removeVariable(String _var_name) {
		if(!hasVariable(_var_name)) {
			System.err.println("Variable \""+_var_name+"\" does not exist in this dataframe!");
			return;
		}
		int var_id = getVariableID(_var_name)-Constants.FIRST_IDX;
		if(titles.length==1) { titles = new String[0]; types = new DataType[0]; datalength = 0; }
		else {
			String[] stemp = new String[titles.length-1];
			for(int s=0; s<stemp.length; s++) stemp[s] = titles[s+(s>=var_id?1:0)];
			titles = new String[stemp.length];
			for(int s=0; s<stemp.length; s++) titles[s] = stemp[s];
			DataType[] ttemp = new DataType[stemp.length];
			for(int t=0; t<ttemp.length; t++) ttemp[t] = types[t+(t>=var_id?1:0)];
			types = new DataType[stemp.length];
			for(int t=0; t<ttemp.length; t++) types[t] = ttemp[t];
		}
	}
	public void changeVariableType(String _var_name, String _new_type) {
		changeVariableType(getVariableID(_var_name), DataType.getDataType(_new_type, default_data_type));
	}
	public void changeVariableType(String _var_name, DataType _new_type) {
		changeVariableType(getVariableID(_var_name), _new_type);
	}
	public void changeVariableType(int _var_id, String _new_type) {
		changeVariableType(_var_id, DataType.getDataType(_new_type, default_data_type));
	}
	public void changeVariableType(int _var_id, DataType _new_type) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) {
			System.err.println("Cannot find variable, so no datatype change can be made!");
			DataHelper.printStackTrace(System.err); return;
		}
		if(types[vi]==_new_type) {
			System.out.println("NOTE: variable is already of type "+_new_type.name()+", so no change is made!");
			return;
		}
		if(_new_type==DataType.STRUCT) {
			System.out.println("Cannot convert to \"struct\".");
			return;
		}
		String variable_name = titles[vi];
		switch(types[vi]) {
			case BOOL: boolean[] bool_arr = bool_column.get(variable_name);
				switch(_new_type) {
					case BYTE: byte[] barr = new byte[bool_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte)(bool_arr[i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[bool_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short)(bool_arr[i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[bool_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (bool_arr[i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[bool_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = (bool_arr[i]?1L:0L);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[bool_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = (bool_arr[i]?1f:0f);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[bool_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = (bool_arr[i]?1d:0d);
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[bool_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = (bool_arr[i]?"true":"false");
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case BYTE: byte[] byte_arr = byte_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[] oarr = new boolean[byte_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = (byte_arr[i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case SHORT: short[] rarr = new short[byte_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = byte_arr[i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[byte_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = byte_arr[i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[byte_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = byte_arr[i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[byte_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = byte_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[byte_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = byte_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[byte_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+byte_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case SHORT: short[] short_arr = short_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[] oarr = new boolean[short_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = (short_arr[i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[short_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) short_arr[i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case INT: int[] iarr = new int[short_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = short_arr[i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[short_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = short_arr[i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[short_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = short_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[short_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = short_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[short_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+short_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case INT: int[] int_arr = int_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[] oarr = new boolean[int_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = (int_arr[i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[int_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) int_arr[i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[int_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) int_arr[i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case LONG: long[] larr = new long[int_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = int_arr[i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[int_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = int_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[int_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = int_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[int_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+int_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case LONG: long[] long_arr = long_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[] oarr = new boolean[long_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = (long_arr[i]!=0L ? true : false);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[long_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[long_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[long_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (int) long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[long_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[long_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[long_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = long_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[long_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+long_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case FLOAT: float[] float_arr = float_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert float to bool result to true except NaN!");
						boolean[] oarr = new boolean[float_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = !Float.isNaN(float_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[float_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) ((int)(float_arr[i]+0.5f)-(float_arr[i]<0f?1:0));
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[float_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) ((short)(float_arr[i]+0.5f) - (float_arr[i]<0f?1:0));
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[float_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (int) (float_arr[i]+0.5f) - (float_arr[i]<0f ? 1 : 0);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[float_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = (long) (float_arr[i]+0.5f) - (float_arr[i]<0f ? 1 : 0);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[float_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = float_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[float_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = float_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[float_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+float_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case DOUBLE: double[] double_arr = double_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert double to bool result to true except NaN!");
						boolean[] oarr = new boolean[double_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = !Double.isNaN(double_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[double_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) ((int)(double_arr[i]+0.5d) - (double_arr[i]<0d?1:0));
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[double_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) ((int)(double_arr[i]+0.5d) - (double_arr[i]<0d?1:0));
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[double_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (int) (double_arr[i]+0.5d) - (double_arr[i]<0d ? 1 : 0);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[double_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = (long) (double_arr[i]+0.5d) - (double_arr[i]<0d?1L:0L);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[double_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = (float) double_arr[i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[double_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = double_arr[i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[] sarr = new String[double_arr.length]; for(int i=0; i<sarr.length; i++) sarr[i] = ""+double_arr[i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case STRING: String[] string_arr = string_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[] oarr = new boolean[string_arr.length]; for(int i=0; i<oarr.length; i++) oarr[i] = set_boolean(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[] barr = new byte[string_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = set_byte(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[string_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = set_short(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[string_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = set_int(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[string_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = set_long(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[] farr = new float[string_arr.length]; for(int i=0; i<farr.length; i++) farr[i] = set_float(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[] darr = new double[string_arr.length]; for(int i=0; i<darr.length; i++) darr[i] = set_double(string_arr[i]);
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case STRUCT:
				System.err.println("\"struct\" cannot be converted, subclass DataFrame for own usage!");
				DataHelper.printStackTrace(System.err); break;
			default:
				System.err.println("An unexpected error occured, can not determin the datatype of the variable!");
				DataHelper.printStackTrace(System.err); break;
		}
	}
	
	public DataFrame clone() {
		DataFrame dfCopy = new DataFrame();
		for(String var: titles) {
			switch(getVariableType(var)) {
				case BOOL:   dfCopy.addColumn(var, bool_column.get(var).clone()); break;
				case BYTE:   dfCopy.addColumn(var, byte_column.get(var).clone()); break;
				case SHORT:  dfCopy.addColumn(var, short_column.get(var).clone()); break;
				case INT:    dfCopy.addColumn(var, int_column.get(var).clone()); break;
				case LONG:   dfCopy.addColumn(var, long_column.get(var).clone()); break;
				case FLOAT:  dfCopy.addColumn(var, float_column.get(var).clone()); break;
				case DOUBLE: dfCopy.addColumn(var, double_column.get(var).clone()); break;
				case STRING: dfCopy.addColumn(var, string_column.get(var).clone()); break;
//				case STRUCT: dfCopy.addColumn(var, struct_column.get(var).clone()); break;
				default: break;
			}
		}
		dfCopy.copyDimensionsFrom(this);
		return dfCopy;
	}
	public DataFrame concat(DataFrame df2) {
		if(getNumberOfDatapoints()==0) {
			for(String var: df2.allVariableNames()) {
				switch(df2.getVariableType(var)) {
					case BOOL:   addColumn(var, (boolean[])df2.getArray(var)); break;
					case BYTE:   addColumn(var, (byte[])df2.getArray(var)); break;
					case SHORT:  addColumn(var, (short[])df2.getArray(var)); break;
					case INT:    addColumn(var, (int[])df2.getArray(var)); break;
					case LONG:   addColumn(var, (long[])df2.getArray(var)); break;
					case FLOAT:  addColumn(var, (float[])df2.getArray(var)); break;
					case DOUBLE: addColumn(var, (double[])df2.getArray(var)); break;
					case STRING: addColumn(var, (String[])df2.getArray(var)); break;
					case STRUCT: System.out.println("Cannot concat struct Variable \""+var+"\" will not be included."); break;
				}
			}
			setDimension(df2.getDimensionValues(), df2.getDimensionName());
			setDimensionsAttributes(df2.getDimensionAttributes());
			return this;
		}
		
		
		boolean haveDifferentKeys = false;
		boolean haveDifferentTypes = false;
		for(String tit: df2.allVariableNames()) {
			if(this.hasVariable(tit)) {
				if(this.getVariableType(tit)!=df2.getVariableType(tit)) {
					System.err.println("Variable \""+tit+"\" is of has different types in dataframe: "+
							this.getVariableType(tit).name()+" <> "+df2.getVariableType(tit).name());
					DataHelper.printStackTrace(System.err);
					haveDifferentTypes = true; break;
				}
			} else {
				System.err.println("Found var \""+tit+"\" in second DataFrame, but not in first one.");
				DataHelper.printStackTrace(System.err);
				haveDifferentKeys = true; break;
			}
		}
		if(haveDifferentTypes)
			return this;
		
		if(!haveDifferentKeys) {
			for(String tit: this.titles) {
				if(!df2.hasVariable(tit)) {
					System.err.println("Found var \""+tit+"\" in first DataFrame, but not in second one.");
					DataHelper.printStackTrace(System.err);
					haveDifferentKeys = true; break;
				}
			}
		}
		if(haveDifferentKeys)
			return this;
		
		for(int iv=0; iv<titles.length; iv++) {
			int ivf = iv+Constants.FIRST_IDX;
			switch(types[iv]) {
				case BOOL:   bool_column.put(  titles[iv], DataHelper.concat_bool_array(  (boolean[])this.getArray(ivf), (boolean[])df2.getArray(titles[iv]))); break;
				case BYTE:   byte_column.put(  titles[iv], DataHelper.concat_byte_array(  (byte[])this.getArray(ivf),    (byte[])df2.getArray(titles[iv])));    break;
				case SHORT:  short_column.put( titles[iv], DataHelper.concat_short_array( (short[])this.getArray(ivf),   (short[])df2.getArray(titles[iv])));   break;
				case INT:    int_column.put(   titles[iv], DataHelper.concat_int_array(   (int[])this.getArray(ivf),     (int[])df2.getArray(titles[iv])));     break;
				case LONG:   long_column.put(  titles[iv], DataHelper.concat_long_array(  (long[])this.getArray(ivf),    (long[])df2.getArray(titles[iv])));    break;
				case FLOAT:  float_column.put( titles[iv], DataHelper.concat_float_array( (float[])this.getArray(ivf),   (float[])df2.getArray(titles[iv])));   break;
				case DOUBLE: double_column.put(titles[iv], DataHelper.concat_double_array((double[])this.getArray(ivf),  (double[])df2.getArray(titles[iv])));  break;
				case STRING: string_column.put(titles[iv], DataHelper.concat_string_array((String[])this.getArray(ivf),  (String[])df2.getArray(titles[iv])));  break;
				case STRUCT: System.out.println("Cannot concat struct Variable \""+titles[iv]+"\" will not be included."); break;
				default: break;
			}
		}
		dimension = DataHelper.concat_double_array(dimension, df2.getDimensionValues());
//		if(!dimension_name.equals(df2.getDimensionName()))
//			System.err.println("[WARNING] dimension names of DataFrames are not the same! Got \""+dimension_name+"\" and \""+df2.getDimensionName()+"\"");
		attribs_dim.putAll(df2.getDimensionAttributes());
		for(int iv=titles.length-1; iv>=0; iv--) {
			if(types[iv]==DataType.STRUCT)
				removeVariable(iv);
		}
		datalength += df2.datalength;
		return this;
	}
	public DataFrame2D toDataFrame2D(String dim_one_name, String dim_two_name) {
		DataFrame2D res2D = new DataFrame2D();
		//check if the DataFrame contain these dimensions!
		int dim1ID = DataHelper.strings_index(titles, dim_one_name);
		if(dim1ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_one_name+"\"!");
			DataHelper.printStackTrace(System.err); return res2D; }
		int dim2ID = DataHelper.strings_index(titles, dim_two_name);
		if(dim2ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_two_name+"\"!");
			DataHelper.printStackTrace(System.err); return res2D; }
		//collect entries for each dimension
		Map<Double,List<Integer>> values_dim_one = new HashMap<>();
		Map<Double,List<Integer>> values_dim_two = new HashMap<>();
		for(int dim=0; dim<2; dim++) {
			int dimID = (dim==0 ? dim1ID : dim2ID);
			DataType dt = (dim==0 ? types[dim1ID] : types[dim2ID]);
			Map<Double,List<Integer>> to_fill = (dim==0 ? values_dim_one : values_dim_two);
			switch(dt) {
				case BOOL: boolean[] bool_arr = bool_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = bool_arr[dl] ? 1d : 0d;
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case BYTE: byte[] byte_arr = byte_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) byte_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case SHORT: short[] short_arr = short_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) short_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case INT: int[] int_arr = int_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) int_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case LONG: long[] long_arr = long_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) long_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case FLOAT: float[] float_arr = float_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Float.isNaN(float_arr[dl])) { System.out.println("WARNING: dimension variable contains missing value"); continue;}
						Double e = (double) float_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case DOUBLE: double[] double_arr = double_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Double.isNaN(double_arr[dl])) { System.out.println("WARNING: dimension variable contains missing values"); continue;}
						Double e = double_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case STRING:
				case STRUCT:
				default:
					String dim_name = (dim==0 ? dim_one_name : dim_two_name);
					System.err.println("DataType \""+dt.name()+"\" of variable \""+dim_name+"\" is not supported as dimension!");
					return res2D;
			}
		}
		//create array for each dimension
		int[] dim_lengths = new int[2];
		int[] dim_counter_length = { 0, 0 };
		dim_lengths[0] = values_dim_one.keySet().size();
		for(Double e: values_dim_one.keySet()) dim_counter_length[0] = Math.max(dim_counter_length[0], values_dim_one.get(e).size());
		dim_lengths[1] = values_dim_two.keySet().size();
		for(Double e: values_dim_two.keySet()) dim_counter_length[1] = Math.max(dim_counter_length[1], values_dim_two.get(e).size());
		if(dim_lengths[0]==0 || dim_lengths[1]==0) {
			System.err.println("At least one dimension has length 0!\nCannot fill the dataframe!");
			return res2D;
		}
		if(dim_counter_length[0]>dim_lengths[1] || dim_counter_length[1]>dim_lengths[0]) {
			System.err.println("At least one combination of dimension-variables-values has more than one entry!\nCannot fill the dataframe!");
			return res2D;
		}
		Double[] temp = values_dim_one.keySet().toArray(new Double[0]);
		double[] dimension_one = new double[dim_lengths[0]];
		for(int dl=0; dl<dim_lengths[0]; dl++) dimension_one[dl] = temp[dl].doubleValue();
		DataHelper.sortem(dimension_one);
		temp = values_dim_two.keySet().toArray(new Double[0]);
		double[] dimension_two = new double[dim_lengths[1]];
		for(int dl=0; dl<dim_lengths[1]; dl++) dimension_two[dl] = temp[dl].doubleValue();
		DataHelper.sortem(dimension_two);
		//create index array for mapping from 1D to 3D
		int[][] indices = new int[dim_lengths[0]][dim_lengths[1]];
		for(int v=0; v<dim_lengths[0]; v++) { List<Integer> lv = values_dim_one.get(dimension_one[v]);
			for(int u=0; u<dim_lengths[1]; u++) { List<Integer> lu = values_dim_two.get(dimension_two[u]);
				indices[v][u] = -1;
				for(Integer i: lv)
					if(lu.contains(i)) { indices[v][u] = i.intValue(); break; }
		}	}
		//for each other variable than the dimensions create a 3D-array
		for(int iv=0; iv<titles.length; iv++) {
			if(iv==dim1ID) continue; if(iv==dim2ID) continue;
			switch(types[iv]) {
				case BOOL: boolean[] oarr = bool_column.get(titles[iv]);
					boolean[][] bool_arr = new boolean[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { bool_arr[v][u] = false; } else { bool_arr[v][u] = oarr[indices[v][u]]; }
					res2D.addColumn(titles[iv], bool_arr);
					break;
				case BYTE: byte[] barr = byte_column.get(titles[iv]);
					byte[][] byte_arr = new byte[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { byte_arr[v][u] = Byte.MIN_VALUE; } else { byte_arr[v][u] = barr[indices[v][u]]; }
					res2D.addColumn(titles[iv], byte_arr);
					break;
				case SHORT: short[] rarr = short_column.get(titles[iv]);
					short[][] short_arr = new short[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { short_arr[v][u] = Short.MIN_VALUE; } else { short_arr[v][u] = rarr[indices[v][u]]; }
					res2D.addColumn(titles[iv], short_arr);
					break;
				case INT: int[] iarr = int_column.get(titles[iv]);
					int[][] int_arr = new int[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { int_arr[v][u] = Integer.MIN_VALUE; } else { int_arr[v][u] = iarr[indices[v][u]]; }
					res2D.addColumn(titles[iv], int_arr);
					break;
				case LONG: long[] larr = long_column.get(titles[iv]);
					long[][] long_arr = new long[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { long_arr[v][u] = Long.MIN_VALUE; } else { long_arr[v][u] = larr[indices[v][u]]; }
					res2D.addColumn(titles[iv], long_arr);
					break;
				case FLOAT: float[] farr = float_column.get(titles[iv]);
					float[][] float_arr = new float[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { float_arr[v][u] = Float.NaN; } else { float_arr[v][u] = farr[indices[v][u]]; }
					res2D.addColumn(titles[iv], float_arr);
					break;
				case DOUBLE: double[] darr = double_column.get(titles[iv]);
					double[][] double_arr = new double[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { double_arr[v][u] = Double.NaN; } else { double_arr[v][u] = darr[indices[v][u]]; }
					res2D.addColumn(titles[iv], double_arr);
					break;
				case STRING: String[] sarr = string_column.get(titles[iv]);
					String[][] string_arr = new String[dim_lengths[0]][dim_lengths[1]];
					for(int v=0; v<dim_lengths[0]; v++) for(int u=0; u<dim_lengths[1]; u++)
						if(indices[v][u]==-1) { string_arr[v][u] = ""; } else { string_arr[v][u] = sarr[indices[v][u]]; }
					res2D.addColumn(titles[iv], string_arr);
					break;
				case STRUCT:
					System.out.println("Cannot reshape struct variable \""+titles[iv]+"\".");
					break;
				default:
					System.err.println("An unexpected datatype occure!");
					return new DataFrame2D();
			}
		}
		res2D.setDimension(0+Constants.FIRST_IDX, dimension_one, dim_one_name);
		res2D.setDimension(1+Constants.FIRST_IDX, dimension_two, dim_two_name);
		return res2D;
	}
	public DataFrame3D toDataFrame3D(String dim_one_name, String dim_two_name, String dim_three_name) {
		DataFrame3D res3D = new DataFrame3D();
		//check if the DataFrame contain these dimensions!
		int dim1ID = DataHelper.strings_index(titles, dim_one_name);
		if(dim1ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_one_name+"\"!");
			DataHelper.printStackTrace(System.err); return res3D; }
		int dim2ID = DataHelper.strings_index(titles, dim_two_name);
		if(dim2ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_two_name+"\"!");
			DataHelper.printStackTrace(System.err); return res3D; }
		int dim3ID = DataHelper.strings_index(titles, dim_three_name);
		if(dim3ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_three_name+"\"!");
			DataHelper.printStackTrace(System.err); return res3D; }
		//collect entries for each dimension
		Map<Double,List<Integer>> values_dim_one = new HashMap<>();
		Map<Double,List<Integer>> values_dim_two = new HashMap<>();
		Map<Double,List<Integer>> values_dim_three = new HashMap<>();
		for(int dim=0; dim<3; dim++) {
			int dimID = (dim==0 ? dim1ID : dim==1 ? dim2ID : dim3ID);
			DataType dt = (dim==0 ? types[dim1ID] : dim==1 ? types[dim2ID] : types[dim3ID]);
			Map<Double,List<Integer>> to_fill = (dim==0 ? values_dim_one : dim==1 ? values_dim_two : values_dim_three);
			switch(dt) {
				case BOOL: boolean[] bool_arr = bool_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = bool_arr[dl] ? 1d : 0d;
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case BYTE: byte[] byte_arr = byte_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) byte_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case SHORT: short[] short_arr = short_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) short_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case INT: int[] int_arr = int_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) int_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case LONG: long[] long_arr = long_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = (double) long_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case FLOAT: float[] float_arr = float_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Float.isNaN(float_arr[dl])) { System.out.println("WARNING: dimension variable contains missing value"); continue;}
						Double e = (double) float_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case DOUBLE: double[] double_arr = double_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Double.isNaN(double_arr[dl])) { System.out.println("WARNING: dimension variable contains missing values"); continue;}
						Double e = double_arr[dl];
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(dl); }
						else { List<Integer> l = new ArrayList<>(); l.add(dl); to_fill.put(e, l); }
					} break;
				case STRING:
				case STRUCT:
				default:
					String dim_name = (dim==0 ? dim_one_name : dim==1 ? dim_two_name : dim_three_name);
					System.err.println("DataType \""+dt.name()+"\" of variable \""+dim_name+"\" is not supported as dimension!");
					return res3D;
			}
		}
		//create array for each dimension
		int[] dim_lengths = new int[3];
		int[] dim_counter_length = { 0, 0, 0};
		dim_lengths[0] = values_dim_one.keySet().size();
		for(Double e: values_dim_one.keySet()) dim_counter_length[0] = Math.max(dim_counter_length[0], values_dim_one.get(e).size());
		dim_lengths[1] = values_dim_two.keySet().size();
		for(Double e: values_dim_two.keySet()) dim_counter_length[1] = Math.max(dim_counter_length[1], values_dim_two.get(e).size());
		dim_lengths[2] = values_dim_three.keySet().size();
		for(Double e: values_dim_three.keySet()) dim_counter_length[2] = Math.max(dim_counter_length[2], values_dim_three.get(e).size());
		if(dim_counter_length[0]>dim_lengths[1]*dim_lengths[2] || dim_counter_length[1]>dim_lengths[0]*dim_lengths[2] ||
				dim_counter_length[2]>dim_lengths[0]*dim_lengths[1]) {
			System.err.println("At least one combination of dimension-variables-values has more than one entry!\nCannot fill the dataframe!");
			return res3D;
		}
		Double[] temp = values_dim_one.keySet().toArray(new Double[0]);
		double[] dimension_one = new double[dim_lengths[0]];
		for(int dl=0; dl<dim_lengths[0]; dl++) dimension_one[dl] = temp[dl].doubleValue();
		DataHelper.sortem(dimension_one);
		temp = values_dim_two.keySet().toArray(new Double[0]);
		double[] dimension_two = new double[dim_lengths[1]];
		for(int dl=0; dl<dim_lengths[1]; dl++) dimension_two[dl] = temp[dl].doubleValue();
		DataHelper.sortem(dimension_two);
		temp = values_dim_three.keySet().toArray(new Double[0]);
		double[] dimension_thr = new double[dim_lengths[2]];
		for(int dl=0; dl<dim_lengths[2]; dl++) dimension_thr[dl] = temp[dl].doubleValue();
		DataHelper.sortem(dimension_thr);
		//create index array for mapping from 1D to 3D
		int[][][] indices = new int[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
		for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
			indices[w][v][u] = -1;
		for(int w=0; w<dim_lengths[0]; w++) { List<Integer> lw = values_dim_one.get(dimension_one[w]);
			for(int v=0; v<dim_lengths[1]; v++) { List<Integer> lv = values_dim_two.get(dimension_two[v]);
				for(int u=0; u<dim_lengths[2]; u++) { List<Integer> lu = values_dim_three.get(dimension_thr[u]);
					indices[w][v][u] = -1;
					for(Integer i: lw)
						if(lv.contains(i) && lu.contains(i)) { indices[w][v][u] = i.intValue(); break; }
		}	}	}
		//for each other variable than the dimensions create a 3D-array
		for(int iv=0; iv<titles.length; iv++) {
			if(iv==dim1ID) continue; if(iv==dim2ID) continue; if(iv==dim3ID) continue;
			switch(types[iv]) {
				case BOOL: boolean[] oarr = bool_column.get(titles[iv]);
					boolean[][][] bool_arr = new boolean[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { bool_arr[w][v][u] = false; } else { bool_arr[w][v][u] = oarr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], bool_arr);
					break;
				case BYTE: byte[] barr = byte_column.get(titles[iv]);
					byte[][][] byte_arr = new byte[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { byte_arr[w][v][u] = Byte.MIN_VALUE; } else { byte_arr[w][v][u] = barr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], byte_arr);
					break;
				case SHORT: short[] rarr = short_column.get(titles[iv]);
					short[][][] short_arr = new short[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { short_arr[w][v][u] = Short.MIN_VALUE; } else { short_arr[w][v][u] = rarr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], short_arr);
					break;
				case INT: int[] iarr = int_column.get(titles[iv]);
					int[][][] int_arr = new int[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { int_arr[w][v][u] = Integer.MIN_VALUE; } else { int_arr[w][v][u] = iarr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], int_arr);
					break;
				case LONG: long[] larr = long_column.get(titles[iv]);
					long[][][] long_arr = new long[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { long_arr[w][v][u] = Long.MIN_VALUE; } else { long_arr[w][v][u] = larr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], long_arr);
					break;
				case FLOAT: float[] farr = float_column.get(titles[iv]);
					float[][][] float_arr = new float[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { float_arr[w][v][u] = Float.NaN; } else { float_arr[w][v][u] = farr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], float_arr);
					break;
				case DOUBLE: double[] darr = double_column.get(titles[iv]);
					double[][][] double_arr = new double[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { double_arr[w][v][u] = Double.NaN; } else { double_arr[w][v][u] = darr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], double_arr);
					break;
				case STRING: String[] sarr = string_column.get(titles[iv]);
					String[][][] string_arr = new String[dim_lengths[0]][dim_lengths[1]][dim_lengths[2]];
					for(int w=0; w<dim_lengths[0]; w++) for(int v=0; v<dim_lengths[1]; v++) for(int u=0; u<dim_lengths[2]; u++)
						if(indices[w][v][u]==-1) { string_arr[w][v][u] = ""; } else { string_arr[w][v][u] = sarr[indices[w][v][u]]; }
					res3D.addColumn(titles[iv], string_arr);
					break;
				case STRUCT:
					System.out.println("Cannot reshape struct variable \""+titles[iv]+"\".");
					break;
				default:
					System.err.println("An unexpected datatype occure!");
					return new DataFrame3D();
			}
		}
		res3D.setDimension(0+Constants.FIRST_IDX, dimension_one, dim_one_name);
		res3D.setDimension(1+Constants.FIRST_IDX, dimension_two, dim_two_name);
		res3D.setDimension(2+Constants.FIRST_IDX, dimension_thr, dim_three_name);
		return res3D;
	}
	public DataFrame filterSubDataFrame(String condition) {
		return filterSubDataFrame(condition, false);
	}
	public DataFrame filterSubDataFrame(String condition, boolean inline) {
		String[][] conditions = LogicHelper.readConditions(this, condition);
		if(conditions==null || conditions.length==0) {
			System.err.println("Could not convert condition \""+condition+"\" into mask...");
			DataHelper.printStackTrace(System.err);
			return this;
		}
//		boolean[] rowdivs = new boolean[conditions.length-1]; for(int r=0; r<rowdivs.length; r++) rowdivs[r] = false;
//		boolean[] coldivs = new boolean[conditions[0].length-1]; for(int c=0; c<coldivs.length; c++) coldivs[c] = true;
//		FormatHelper.printTable(conditions,coldivs,rowdivs);
		boolean[] mask = new boolean[datalength];
		for(int e=0; e<datalength; e++)
			mask[e] = LogicHelper.evaluateConditionLines(this, e, conditions);
		return filterByMask(mask, inline);
	}
	public DataFrame filterByMask(boolean[] mask) {
		return filterByMask(mask, false);
	}
	public DataFrame filterByMask(boolean[] mask, boolean inline) {
		DataFrame res = new DataFrame();
		int count = 0;
		for (int i = 0; i < mask.length; i++) {
			if (mask[i])
				count++;
		}
		int[] ids = new int[count];
		int ct = 0;
		for (int i = 0; i < mask.length; i++) {
			if (mask[i]) {
				ids[ct] = i;
				ct++;
			}
		}
		for(String var: allVariableNames()) {
			switch (getVariableType(var)) {
				case BOOL: boolean[] newBool = new boolean[count], oldBool = bool_column.get(var);
					for (int c = 0; c < count; ) { newBool[c] = oldBool[ids[c]]; c++; }
					if(inline) bool_column.put(var, newBool); else res.addColumn(var, newBool); break;
				case BYTE: byte[] newByte = new byte[count], oldByte = byte_column.get(var);
					for (int c = 0; c < count; ) { newByte[c] = oldByte[ids[c]]; c++; }
					if(inline) byte_column.put(var, newByte); else res.addColumn(var, newByte); break;
				case SHORT: short[] newShort = new short[count], oldShort = short_column.get(var);
					for (int c = 0; c < count; ) { newShort[c] = oldShort[ids[c]]; c++; }
					if(inline) short_column.put(var, newShort); else res.addColumn(var, newShort); break;
				case INT: int[] newInt = new int[count], oldInt = int_column.get(var);
					for (int c = 0; c < count; ) { newInt[c] = oldInt[ids[c]]; c++; }
					if(inline) int_column.put(var, newInt); else res.addColumn(var, newInt); break;
				case LONG: long[] newLong = new long[count], oldLong = long_column.get(var);
					for (int c = 0; c < count; ) { newLong[c] = oldLong[ids[c]]; c++; }
					if(inline) long_column.put(var, newLong); else res.addColumn(var, newLong); break;
				case FLOAT: float[] newFloat = new float[count], oldFloat = float_column.get(var);
					for (int c = 0; c < count; ) { newFloat[c] = oldFloat[ids[c]]; c++; }
					if(inline) float_column.put(var, newFloat); else res.addColumn(var, newFloat); break;
				case DOUBLE: double[] newDouble = new double[count], oldDouble = double_column.get(var);
					for (int c = 0; c < count; ) { newDouble[c] = oldDouble[ids[c]]; c++; }
					if(inline) double_column.put(var, newDouble); else res.addColumn(var, newDouble); break;
				case STRING: String[] newString = new String[count], oldString = string_column.get(var);
					for (int c = 0; c < count; ) { newString[c] = oldString[ids[c]]; c++; }
					if(inline) string_column.put(var, newString); else res.addColumn(var, newString); break;
				default:
					System.err.println("Unknown datatype, cannot extract and mask variable <" + var + "> from DataFrame."); break; }
		}
		double[] newDim = new double[count];
		for(int c=0; c<count; c++) newDim[c] = dimension[ids[c]];
		if(inline) {
			datalength = count;
			dimension = newDim.clone();
			return this;
		} else {
			res.setDimension(newDim, dimension_name);
			res.setDimensionsAttributes(attribs_dim);
			return res;
		}
	}
	public DataFrame reorderByIndexList(int[] ids) {
		return this.reorderByIndexList(ids, false);
	}
	public DataFrame reorderByIndexList(int[] ids, boolean inline) {
		int newlength = ids.length;
		DataFrame res = new DataFrame();
		for(String var: allVariableNames()) {
			switch (getVariableType(var)) {
				case BOOL: boolean[] newBool = new boolean[newlength], oldBool = bool_column.get(var);
					for (int c = 0; c < newlength; ) { newBool[c] = oldBool[ids[c]]; c++; }
					if(inline) bool_column.put(var, newBool); else res.addColumn(var, newBool); break;
				case BYTE: byte[] newByte = new byte[newlength], oldByte = byte_column.get(var);
					for (int c = 0; c < newlength; ) { newByte[c] = oldByte[ids[c]]; c++; }
					if(inline) byte_column.put(var, newByte); else res.addColumn(var, newByte); break;
				case SHORT: short[] newShort = new short[newlength], oldShort = short_column.get(var);
					for (int c = 0; c < newlength; ) { newShort[c] = oldShort[ids[c]]; c++; }
					if(inline) short_column.put(var, newShort); else res.addColumn(var, newShort); break;
				case INT: int[] newInt = new int[newlength], oldInt = int_column.get(var);
					for (int c = 0; c < newlength; ) { newInt[c] = oldInt[ids[c]]; c++; }
					if(inline) int_column.put(var, newInt); else res.addColumn(var, newInt); break;
				case LONG: long[] newLong = new long[newlength], oldLong = long_column.get(var);
					for (int c = 0; c < newlength; ) { newLong[c] = oldLong[ids[c]]; c++; }
					if(inline) long_column.put(var, newLong); else res.addColumn(var, newLong); break;
				case FLOAT: float[] newFloat = new float[newlength], oldFloat = float_column.get(var);
					for (int c = 0; c < newlength; ) { newFloat[c] = oldFloat[ids[c]]; c++; }
					if(inline) float_column.put(var, newFloat); else res.addColumn(var, newFloat); break;
				case DOUBLE: double[] newDouble = new double[newlength], oldDouble = double_column.get(var);
					for (int c = 0; c < newlength; ) { newDouble[c] = oldDouble[ids[c]]; c++; }
					if(inline) double_column.put(var, newDouble); else res.addColumn(var, newDouble); break;
				case STRING: String[] newString = new String[newlength], oldString = string_column.get(var);
					for (int c = 0; c < newlength; ) { newString[c] = oldString[ids[c]]; c++; }
					if(inline) string_column.put(var, newString); else res.addColumn(var, newString); break;
				default:
					System.err.println("Unknown datatype, cannot extract and mask variable <" + var + "> from DataFrame.");
					break;
			}
			if(inline) recalcStats(getVariableID(var));
		}
		double[] newDim = new double[newlength];
		for(int c=0; c<newlength; c++) newDim[c] = dimension[ids[c]];
		if(inline) {
			datalength = newlength;
			dimension = newDim.clone();
			return this;
		} else {
			res.setDimension(newDim, dimension_name);
			return res;
		}
	}
	
	public void setDefaultDatatype(String _type) { default_data_type = DataType.getDataType(_type, default_data_type); }
	public String getDefaultDatatype() { return default_data_type.toString(); }
	
	/**
	 * 
	 * @param file_path
	 * @param delimiter
	 * @param header_length
	 * @param data_types
	 */
	@Deprecated
	public void read_csv(String file_path, String[] data_types) {
		read_csv(file_path, ",", 0, data_types);
	}
	@Deprecated
	public void read_csv(String file_path, String delimiter, String[] data_types) {
		read_csv(file_path, delimiter, 0, data_types);
	}
	@Deprecated
	public void read_csv(String file_path, int header_length, String[] data_types) {
		read_csv(file_path, ",", header_length, data_types);
	}
	@Deprecated
	public void read_csv(String file_path, String delimiter, int header_length, String[] data_types) {
		try(BufferedReader br = new BufferedReader(new FileReader(new File(file_path)))) {
			boolean useGivenTypes = true;
			if(data_types==null) {
				useGivenTypes = false;
			} else {
				useGivenTypes = (data_types.length==0);
			}
			String line = br.readLine();
			int linenumber = 0;
			List<String[]> data = new ArrayList<String[]>();
			String[] pre_types = new String[0];
			while(line!=null) {
				//System.out.println(line);
				linenumber++;
				if(linenumber<=header_length) {
					if(header_length>=2) {
						if(linenumber+1==header_length)
							titles = break_line(line,delimiter);
						if(linenumber+0==header_length)
							pre_types = break_line(line,delimiter);
					}
					if(header_length==1) {
						titles = break_line(line,delimiter);
						pre_types = (useGivenTypes ? data_types : null_line(line,delimiter));
					}
					line = br.readLine();
					continue;
				}
				if(linenumber==1) { titles = null_line(line,delimiter); pre_types = (useGivenTypes ? data_types : null_line(line,delimiter)); }
				data.add(break_line(line,delimiter));
				//System.out.println(line);
				line = br.readLine();
			}
			//process header/titles
			int clen = titles.length; // number of columns
			datalength = data.size();
			if(clen!=pre_types.length) throw new IllegalArgumentException("Number of columns isn't constant through the dataset!");
			types = new DataType[titles.length];
			for(int t=0; t<clen; t++) {
				if(titles[t].length()==0) titles[t] = createTitle(titles);
				if(pre_types[t].length()==0) types[t] = default_data_type;
				types[t] = DataType.getDataType(pre_types[t],default_data_type);
				switch(types[t]) {
					case BOOL:   bool_column.put(titles[t], new boolean[datalength]); break;
					case BYTE:   byte_column.put(titles[t], new byte[datalength]); break;
					case INT:    int_column.put(titles[t], new int[datalength]); break;
					case SHORT:  short_column.put(titles[t], new short[datalength]); break;
					case LONG:   long_column.put(titles[t], new long[datalength]); break;
					case FLOAT:  float_column.put(titles[t], new float[datalength]); break;
					case DOUBLE: double_column.put(titles[t], new double[datalength]); break;
					default:
					case STRING: string_column.put(titles[t], new String[datalength]); break;
				}
			}
			for(int c=0; c<clen; c++) {
				switch(types[c]) {
					case BOOL:
						boolean[] o_arr = bool_column.get(titles[c]);
						for(int de=0; de<datalength; de++) o_arr[de] = set_boolean(data.get(de)[c]);
						bool_column.put(titles[c], o_arr);
						break;
					case BYTE:
						byte[] b_arr = byte_column.get(titles[c]);
						for(int de=0; de<datalength; de++) b_arr[de] = set_byte(data.get(de)[c]);
						byte_column.put(titles[c], b_arr);
						break;
					case INT:
						int[] i_arr = int_column.get(titles[c]);
						for(int de=0; de<datalength; de++) i_arr[de] = set_int(data.get(de)[c]);
						int_column.put(titles[c], i_arr);
						break;
					case SHORT:
						short[] r_arr = short_column.get(titles[c]);
						for(int de=0; de<datalength; de++) r_arr[de] = set_short(data.get(de)[c]);
						short_column.put(titles[c], r_arr);
						break;
					case LONG:
						long[] l_arr = long_column.get(titles[c]);
						for(int de=0; de<datalength; de++) l_arr[de] = set_long(data.get(de)[c]);
						long_column.put(titles[c], l_arr);
						break;
					case FLOAT:
						float[] f_arr = float_column.get(titles[c]);
						for(int de=0; de<datalength; de++) f_arr[de] = set_float(data.get(de)[c]);
						float_column.put(titles[c], f_arr);
						break;
					case DOUBLE:
						double[] d_arr = double_column.get(titles[c]);
						for(int de=0; de<datalength; de++) d_arr[de] = set_double(data.get(de)[c]);
						double_column.put(titles[c], d_arr);
						break;
					default:
					case STRING:
						String[] s_arr = string_column.get(titles[c]);
						for(int de=0; de<datalength; de++) s_arr[de] = data.get(de)[c];
						string_column.put(titles[c], s_arr);
						break;
				}
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
			clear();
		} catch(IllegalArgumentException iae) {
			iae.printStackTrace();
			clear();
		}
	}
	
	/**
	 * Read data from a text file in simplified GeoEAS format
	 * @param _file_path path to the text file
	 */
	@Deprecated
	public void read_gam_dat(String _file_path) {
		if(!_file_path.endsWith(".dat") && !_file_path.endsWith(".out"))
			System.out.println("WARNING: Standard textfiles for GSLIB ends with '.dat' or '.out'");
		int nvari,maxdat;
		String[] names = null;
//      c Note VERSION number:
		//System.out.println(" Based on GAM Version: "+FormatHelper.nf(GAMV.VERSION,5,3));

//      c Check to make sure the data file exists, then either read in the
//      c data or write an error message and stop:
		File f = new File(_file_path);
		if(!f.exists()) throw new RuntimeException("ERROR data file "+_file_path+" does not exist!");

//      c The data file exists so open the file and read in the header
//      c information. Initialize the storage that will be used to collect
//      c the data found in the file:
		try(BufferedReader br = new BufferedReader(new FileReader(f))) {
			List<String> lines = new ArrayList<String>();
			br.mark(1024*1024);
			String line = br.readLine();
			lines.add(""+line);
			line = br.readLine(); String[] parts = FormatHelper.splitBySpace(line); nvari = Integer.parseInt(parts[0].trim());
			lines.add(""+line);
			for(int i=0; i<nvari; i++) line = br.readLine();
			maxdat = -1;
			while(line!=null) {
				maxdat++;
				line = br.readLine();
			}
			System.out.println("maxdat = "+maxdat);
			
			names = new String[nvari];
			double[][] vr = new double[nvari][maxdat];
			
			br.reset();
			line = br.readLine();
			line = br.readLine(); parts = FormatHelper.splitBySpace(line); nvari = Integer.parseInt(parts[0].trim());
			for(int i=0; i<nvari; i++) {
				line = br.readLine(); //parts = FormatHelper.splitBySpace(line);
				names[i] = line; //parts[0];
				System.out.println(" found variable name \""+names[i]+"\"");
			}
			System.out.println(" nvari="+nvari);

//      c Read all the data until the end of the file:
			int nd = -1;
			while(true) {
				line = br.readLine();
				if(line==null || line.length()==0) break;
				parts = FormatHelper.splitBySpace(line);
				nd++;
				for(int j=0; j<nvari; j++) {
					//vr[j][nd] = Double.parseDouble(parts[j].trim());
					vr[j][nd] = set_double(parts[j].trim()); //with error check for reading errors and NaNs
				}
			}
			nd++;
			System.out.println("nd     = "+nd);
			//System.out.println("maxdat = "+maxdat);
			
			this.clear();
			//titles = new String[names.length];
			//types  = new DataType[names.length];
			datalength = nd;
			//default_data_type = DataType.DOUBLE;
			for(int iv=0; iv<nvari; iv++) {
//				String ttt = names[iv];
//				if(ttt==null) ttt = createTitle(names);
//				if(ttt.length()==0) ttt = createTitle(names);
//				names[iv] = ttt;
//				add_variable(names[iv], DataType.DOUBLE);
//				//titles[iv] = names[iv];
//				//types[iv] = DataType.DOUBLE;
//				double_column.put(names[iv], vr[iv]);
//				double[] inax = StdAnalysis.minmax(vr[iv]);
//				add_stdana(inax[0], inax[1], StdAnalysis.mean(vr[iv]), StdAnalysis.variance(vr[iv]));
				addColumn(names[iv], vr[iv]);
			}
			br.close();
		} catch(IOException | NullPointerException io_e) {
			throw new RuntimeException("ERROR in data file!", io_e);
		}
		System.out.println("[DEBUG] all data read");
		System.out.print("\n"); head();
		System.out.print("\n"); describe();
		System.out.print("\n");
	}
	
	/**
	 * Read selected variable data from a Netcdf file
	 * if on variable isn't in the file, it would not be read to this dataframe
	 * and a warning would be printed to the command line
	 * @param filepath           path to the NetCDF file
	 * @param include_dimensions wether dimensions of variables should be includes as extra data columns
	 * @param variable           name(s) of variable(s)
	 * @return                   return this DataFrame
	 * @throws IOException 
	 */
	@Deprecated
	public DataFrame readFromNetcdf(String filepath, boolean include_dimensions, String... variable) throws IOException {
		NetcdfFile nc = NetcdfFiles.open(filepath);
		this.readFromNetcdf(nc, include_dimensions, variable);
		return this;
	}
	/**
	 * Read selected variable data from a Netcdf file
	 * if on variable isn't in the file, it would not be read to this dataframe
	 * and a warning would be printed to the command line
	 * @param netcdf_file        NetCDF object with the data to read
	 * @param include_dimensions wether dimensions of variables should be includes as extra data columns
	 * @param variable           name(s) of variable(s)
	 * @return                   return this DataFrame
	 */
	@Deprecated
	public DataFrame readFromNetcdf(NetcdfFile netcdf_file, boolean include_dimensions, String... variable) {
		if(netcdf_file==null) {
			System.err.println("The Netcdf file does not exist!");
			return this;
		}
		if(variable==null || variable.length<=0) {
			System.err.println("You have to specify at least one variable you want to read from the Netcdf file!");
			return this;
		}
		clear();
		boolean[] var_exist = new boolean[variable.length];
		List<Dimension> dims = new ArrayList<Dimension>();
		for(int b=0; b<var_exist.length; b++) {
			Variable var = netcdf_file.findVariable(variable[b]);
			var_exist[b] = (var!=null);
			if(var_exist[b]) {
				for(Dimension d: var.getDimensions()) {
					if(!dims.contains(d))
						dims.add(d);
				}
			}
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
		if(include_dimensions) {
			int dimid = 0;
			for(Dimension dim: dims) {
				Variable var = netcdf_file.findVariable(dim.getName());
				if(var==null) {
					int[] arr = new int[datalength];
					for(int dl=0; dl<datalength; dl++) arr[dl] = indices[dl][dimid];
					addColumn(dim.getName(), arr);
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
							addColumn(dim.getName(), bl);
							break;
						case BYTE:
							byte[] etyb = new byte[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); etyb[dl] = a.getByte(ind); }
							addColumn(dim.getName(), etyb);
							break;
						case INT:
							int[] tni = new int[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); tni[dl] = a.getInt(ind); }
							addColumn(dim.getName(), tni);
							break;
						case SHORT:
							short[] trohs = new short[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); trohs[dl] = a.getShort(ind); }
							addColumn(dim.getName(), trohs);
							break;
						case LONG:
							long[] gnol = new long[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnol[dl] = a.getLong(ind); }
							addColumn(dim.getName(), gnol);
							break;
						case FLOAT:
							boolean hasFVf = (var.findAttribute("_FillValue")!=null);
							double ffill = Double.NaN; if(hasFVf) ffill = (float) var.findAttribute("_FillValue").getNumericValue();
							if(hasFVf) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+ffill);
							float[] taolf = new float[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); taolf[dl] = a.getFloat(ind);
								if(hasFVf && taolf[dl]==ffill) taolf[dl]=Float.NaN; }
							addColumn(dim.getName(), taolf);
							break;
						case DOUBLE:
							boolean hasFVd = (var.findAttribute("_FillValue")!=null);
							double dfill = Double.NaN; if(hasFVd) dfill = (double) var.findAttribute("_FillValue").getNumericValue();
							if(hasFVd) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+dfill);
							double[] elbuod = new double[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); elbuod[dl] = a.getDouble(ind);
								if(hasFVd && elbuod[dl]==dfill) elbuod[dl] = Double.NaN; }
							addColumn(dim.getName(), elbuod);
							break;
//						case CHAR:
//						case STRING:
//							String[] gnirts = new String[datalength];
//							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnirts[dl] = ""+a.getChar(ind); }
//							addColumn(dim.getFullName(), gnirts);
//							break;
						default:
							System.out.println("WARNING: could not add variable \""+dim.getName()+"\": unsupported data type!");
							break;
					}
				}
				dimid++;
			}
		}
		for(int vi=0; vi<var_exist.length; vi++) {
			if(!var_exist[vi]) continue;
			Variable var = netcdf_file.findVariable(variable[vi]);
			Array a = null;
			try {
				a = var.read();
			} catch(IOException ioe) {
				//ioe.printStackTrace();
				System.out.println("WARNING: could not read variable \""+variable[vi]+"\": does not at to the dataframe!");
				continue;
			}
			Index ind = a.getIndex();
			List<Dimension> vdims = var.getDimensions();
			int[] dimids = new int[vdims.size()];
			for(int d=0; d<dimids.length; d++) {
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
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); o_arr[dl] = a.getBoolean(ind);
					}
					addColumn(variable[vi], o_arr);
					break;
				case BYTE: byte[] b_arr = new byte[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); b_arr[dl] = a.getByte(ind);
					}
					addColumn(variable[vi], b_arr);
					break;
				case INT: int[] i_arr = new int[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); i_arr[dl] = a.getInt(ind);
					}
					addColumn(variable[vi], i_arr);
					break;
				case SHORT: short[] r_arr = new short[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); r_arr[dl] = a.getShort(ind);
					}
					addColumn(variable[vi], r_arr);
					break;
				case LONG: long[] l_arr = new long[datalength];
					for(int dl=0; dl<datalength; dl++) {
						for(int d=0; d<dimids.length; d++) { indexindex[d] = indices[dl][dimids[d]]; }
						ind.set(indexindex); l_arr[dl] = a.getLong(ind);
					}
					addColumn(variable[vi], l_arr);
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
					addColumn(variable[vi], f_arr);
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
					addColumn(variable[vi], d_arr);
					break;
//					case CHAR:
//					case STRING:
				case STRUCTURE:
					readStructs(var, variable[vi], a, datalength, ind, indexindex, indices, dimids);
					break;
				default:
					System.out.println("WARNING: could not add variable \""+variable[vi]+"\": unsupported data type!");
					break;
			}
		}
		//createDimension();
		return this;
	}
	@Deprecated
	public void readStructs(Variable var, String varname, Array a, int datalength, Index ind, int[] indexindex, int[][] indices, int[] dimids) {
		System.out.println("WARNING: could not add variable \""+varname+"\": unsupported data type!");
	}
	/**
	 * Write all content of this dataframe to a netcdf file
	 * @param netcdf_file_path path to the netcdf file
	 * @throws IOException 
	 * @throws InvalidRangeException 
	 */
	@Deprecated
	@SuppressWarnings("rawtypes")
	public void writeToNetcdf(String netcdf_file_path) throws IOException, InvalidRangeException {
		NetcdfFormatWriter.Builder ncdfWriter = NetcdfFormatWriter.createNewNetcdf4(NetcdfFileFormat.NETCDF4_CLASSIC, netcdf_file_path, null);
		//NetcdfFileWriter ncdfWriter = NetcdfFileWriter.createNew(Version.netcdf4, netcdf_file_path);
		String dim_name = ""; boolean is_used = true; int dim_test_num = -1;
		while(is_used) {
			dim_test_num++; dim_name = "dim"+dim_test_num;
			is_used = false;
			for(String tit: titles) if(tit.equals(dim_name)) { is_used = true; break; }
		}
		Dimension dim = ncdfWriter.addDimension(dim_name, datalength);
		List<Dimension> dims = new ArrayList<>(); dims.add(dim);
		Variable.Builder[] vars = new Variable.Builder[titles.length+1];
		vars[0] = ncdfWriter.addVariable(dim_name, ucar.ma2.DataType.INT, dims);
		for(int iv=0; iv<titles.length; iv++) {
			ucar.ma2.DataType dataType = null;
			switch(types[iv]) {
				case BOOL:   dataType = ucar.ma2.DataType.BOOLEAN; break;
				case BYTE:   dataType = ucar.ma2.DataType.BYTE;    break;
				case SHORT:  dataType = ucar.ma2.DataType.SHORT;   break;
				case INT:    dataType = ucar.ma2.DataType.INT;     break;
				case LONG:   dataType = ucar.ma2.DataType.LONG;    break;
				case FLOAT:  dataType = ucar.ma2.DataType.FLOAT;   break;
				case DOUBLE: dataType = ucar.ma2.DataType.DOUBLE;  break;
				default:
				case STRING: dataType = ucar.ma2.DataType.STRING;  break;
			}
			vars[iv+1] = ncdfWriter.addVariable(FormatHelper.name2CFConvention(titles[iv]), dataType, dims);
			switch(types[iv]) {
				case SHORT:  vars[iv+1].addAttribute(new Attribute("_FillValue", Short.MIN_VALUE));   break;
				case INT:    vars[iv+1].addAttribute(new Attribute("_FillValue", Integer.MIN_VALUE)); break;
				case LONG:   vars[iv+1].addAttribute(new Attribute("_FillValue", Long.MIN_VALUE));    break;
				case FLOAT:  vars[iv+1].addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_F));  break;
				case DOUBLE: vars[iv+1].addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_D)); break;
				default: break;
			}
		}
		ncdfWriter.getRootGroup().addAttribute(new Attribute("history", "created with \""+Constants.NAME+" v"+Constants.VERSION+"\" from Dataframe - JAVA Netcdf "+Constants.NETCDF_VERSION));
		NetcdfFormatWriter ncdf = ncdfWriter.build();
		Variable dimVar = ncdf.findVariable(vars[0].getFullName());
		ArrayInt.D1 dim_arr = new ArrayInt.D1(datalength, false);
		for(int dl=0; dl<datalength; dl++) dim_arr.set(dl,1+dl);
		ncdf.write(dimVar, dim_arr);
		for(int iv=0; iv<titles.length; iv++) {
			switch(types[iv]) {
				case BOOL:
					ArrayBoolean.D1 bool_arr = new ArrayBoolean.D1(datalength);
					boolean[] bool_source = bool_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) bool_arr.set(dl, bool_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), bool_arr);
					break;
				case BYTE:
					ArrayByte.D1 byte_arr = new ArrayByte.D1(datalength, false);
					byte[] byte_source = byte_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) byte_arr.set(dl, byte_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), byte_arr);
					break;
				case SHORT:
					ArrayShort.D1 short_arr = new ArrayShort.D1(datalength, false);
					short[] short_source = short_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) short_arr.set(dl, short_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), short_arr);
					break;
				case INT:
					ArrayInt.D1 int_arr = new ArrayInt.D1(datalength, false);
					int[] int_source = int_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) int_arr.set(dl, int_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), int_arr);
					break;
				case LONG:
					ArrayLong.D1 long_arr = new ArrayLong.D1(datalength, false);
					long[] long_source = long_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) long_arr.set(dl, long_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), long_arr);
					break;
				case FLOAT:
					ArrayFloat.D1 float_arr = new ArrayFloat.D1(datalength);
					float[] float_source = float_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) float_arr.set(dl, float_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), float_arr);
					break;
				case DOUBLE:
					ArrayDouble.D1 double_arr = new ArrayDouble.D1(datalength);
					double[] double_source = double_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) double_arr.set(dl, double_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), double_arr);
					break;
				case STRING:
					ArrayString.D1 string_arr = new ArrayString.D1(datalength);
					String[] string_source = string_column.get(titles[iv]);
					for(int dl=0; dl<datalength; dl++) string_arr.set(dl, string_source[dl]);
					ncdf.write(ncdf.findVariable(vars[iv+1].getFullName()), string_arr);
					break;
				default:
					break;
			}
		}
		ncdf.flush();
		ncdf.close();
	}
	/**
	 * Write all content of these dataframe(2d/3d) objects to one netcdf file.
	 * @param path path to the netcdf file
	 * @param df DataFrame(2D/3D)s
	 * @throws IOException 
	 * @throws InvalidRangeException 
	 */
	@Deprecated
	public static void writeToNetcdf(String path, Object... df) throws IOException, InvalidRangeException {
		if(df==null) return;
		if(df.length==0) return;
		//check, if all objects in df are of type DataFrame/DataFrame2D/DataFrame3D
		int maxDim = 0;
		for(Object obj: df) {
			if(obj instanceof DataFrame) {   maxDim=Math.max(maxDim,1); continue; }
			if(obj instanceof DataFrame2D) { maxDim=Math.max(maxDim,2); continue; }
			if(obj instanceof DataFrame3D) { maxDim=Math.max(maxDim,3); continue; }
			System.err.println("Try to write non DataFrame(2D/3D) object to Netcdf");
			return;
		}
		//System.out.println("\n\n\nBefor Netcdf-creation:");
		//System.out.println("\ncreate Netcdf-file ...");
		NetcdfFormatWriter.Builder builder = NetcdfFormatWriter.createNewNetcdf4(NetcdfFileFormat.NETCDF4_CLASSIC, path, null);
		Map<String, Dimension> dimensions = new HashMap<String, Dimension>();
		Map<String, String> variables = new HashMap<String, String>();
		for(Object obj: df) {
			if(obj instanceof DataFrame)
				setupDimensions(builder, dimensions, variables, (DataFrame)obj, null, null);
			if(obj instanceof DataFrame2D)
				setupDimensions(builder, dimensions, variables, null, (DataFrame2D)obj, null);
			if(obj instanceof DataFrame3D)
				setupDimensions(builder, dimensions, variables, null, null, (DataFrame3D)obj);
		}
		for(Object obj: df) {
			if(obj instanceof DataFrame)
				setupVariables(builder, dimensions, variables, (DataFrame)obj, null, null);
			if(obj instanceof DataFrame2D)
				setupVariables(builder, dimensions, variables, null, (DataFrame2D)obj, null);
			if(obj instanceof DataFrame3D)
				setupVariables(builder, dimensions, variables, null, null, (DataFrame3D)obj);
		}
		builder.getRootGroup().addAttribute(new Attribute("history", "created with \"JKriging "+Constants.VERSION+"\" - JAVA Netcdf " + Constants.NETCDF_VERSION));
		
		NetcdfFormatWriter writer = builder.build();
		//write dimension
		for(Object obj: df) {
			if(obj instanceof DataFrame)   writeDim(writer, variables, (DataFrame)obj, null, null);
			if(obj instanceof DataFrame2D) writeDim(writer, variables, null, (DataFrame2D)obj, null);
			if(obj instanceof DataFrame3D) writeDim(writer, variables, null, null, (DataFrame3D)obj);
		}
		//write variables
		for(Object obj: df) {
			if(obj instanceof DataFrame)   write1D(writer, variables, (DataFrame)obj);
			if(obj instanceof DataFrame2D) write2D(writer, variables, (DataFrame2D)obj);
			if(obj instanceof DataFrame3D) write3D(writer, variables, (DataFrame3D)obj);
		}
		
		writer.flush();
		writer.close();
	}
	
	
	public double[] getDimensionValues() {
		if(dimension.length==0) createDimension();
		return dimension;
	}
	public String getDimensionName() {
		if(dimension.length==0) createDimension();
		return dimension_name;
	}
	public Map<String,String> getDimensionAttributes() {
		if(dimension.length==0) createDimension();
		return attribs_dim;
	}
	public int getNumberOfDatapoints() { return datalength; }
	public int getVariableCount() { return titles.length; }
	public String getVarname(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) return null;
		return titles[vi];
	}
	public int getVariableID(String _var_name) {
		return DataHelper.strings_index(titles, _var_name)+Constants.FIRST_IDX;
	}
	public boolean hasVariable(String _var_name) {
		if(_var_name==null) return false;
		return DataHelper.strings_index(titles, _var_name)>=0;
	}
	public String[] allVariableNames() { return titles; }
	public DataType getVariableType(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) return null;
		return types[vi];
	}
	public DataType getVariableType(String _var_name) {
		return getVariableType(getVariableID(_var_name));
	}
	public DataType[] allVariableTypes() { return types; }
	public double[] getMins() { return minmax_mean_sill[0]; }
	public double getMin(String var_name) { return getMin(getVariableID(var_name)); }
	public double getMin(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=minmax_mean_sill[0].length) return Double.NaN;
		return minmax_mean_sill[0][vi];
	}
	public double[] getMaxs() { return minmax_mean_sill[1]; }
	public double getMax(String var_name) { return getMax(getVariableID(var_name)); }
	public double getMax(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=minmax_mean_sill[1].length) return Double.NaN;
		return minmax_mean_sill[1][vi];
	}
	public double getMean(String var_name) { return getMean(getVariableID(var_name)); }
	public double getMean(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=minmax_mean_sill[2].length) return Double.NaN;
		return minmax_mean_sill[2][vi];
	}
	public double getSill(String var_name) { return getSill(getVariableID(var_name)); }
	public double getSill(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=minmax_mean_sill[3].length) return Double.NaN;
		return minmax_mean_sill[3][vi];
	}
	public Object getArray(int _var_id) {
		return getArray(getVarname(_var_id));
	}
	public Object getArray(String _col_name) {
		if(_col_name==null)
			return null;
		if(DataHelper.strings_index(titles, _col_name)<0)
			return null;
		switch(types[DataHelper.strings_index(titles,_col_name)]) {
			case BOOL:   return bool_column.get(_col_name);
			case BYTE:   return byte_column.get(_col_name);
			case INT:    return int_column.get(_col_name);
			case SHORT:  return short_column.get(_col_name);
			case LONG:   return long_column.get(_col_name);
			case FLOAT:  return float_column.get(_col_name);
			case DOUBLE: return double_column.get(_col_name);
			case STRING: return string_column.get(_col_name);
			case STRUCT: return struct_column.get(_col_name);
			default: return null;
		}
	}
	public void head() {
		head(10);
	}
	public void head(int number_of_lines) {
		int clen = titles.length;
		//System.out.println("number of titles,types,minmaxs = "+clen+","+types.length+","+minmax_mean_sill[0].length);
		if(clen==0) {
			System.err.println("No Table data!");
			return;
		}
		int dlen = (number_of_lines==0 ? datalength+2 : Math.min(number_of_lines+2, datalength+2));
		String[][] out = new String[clen][dlen];
		int[] col_width = new int[clen];
		for(int c=0; c<clen; c++) {
			//System.out.println("Title "+c+" = \""+titles[c]+"\"");
			out[c][0] = " "+titles[c]+" ";
			col_width[c] = titles[c].length()+2;
			switch(types[c]) {
				case BOOL:
					out[c][1] = " bool ";
					boolean[] ocol = bool_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+(ocol[d] ? "true" : "false")+" ";
					break;
				case BYTE:
					out[c][1] = " byte ";
					byte[] bcol = byte_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+bcol[d]+" ";
					break;
				case INT:
					out[c][1] = " int ";
					int[] icol = int_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+icol[d]+" ";
					break;
				case SHORT:
					out[c][1] = " short ";
					short[] rcol = short_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+rcol[d]+" ";
					break;
				case LONG:
					out[c][1] = " long ";
					long[] lcol = long_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+lcol[d]+" ";
					break;
				case FLOAT:
					out[c][1] = " float ";
					float[] fcol = float_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+fcol[d]+" ";
					break;
				case DOUBLE:
					out[c][1] = " double ";
					double[] dcol = double_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+dcol[d]+" ";
					break;
				case STRING:
					out[c][1] = " string ";
					String[] scol = string_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " "+scol[d]+" ";
					break;
				case STRUCT:
					out[c][1] = " struct ";
					for(int d=0; d+2<dlen; d++) out[c][d+2] = " ??? ";
					break;
			}
		}
		boolean[] coldivs = new boolean[clen-1];
		for(int c=0; c<clen-1; c++) coldivs[c] = true;
		boolean[] rowdivs = new boolean[dlen-1];
		for(int r=0; r<dlen-1; r++) rowdivs[r] = false;
		rowdivs[0] = true; if(dlen>2) rowdivs[1] = true;
		FormatHelper.printTable(out, coldivs, rowdivs, true);
	}
	public void printfull() { head(0); }
	public void describe() {
		int clen = titles.length;
		String[][] out = new String[10][clen+1];
		out[0][0] = " "; out[1][0] = " Type "; out[2][0] = " Count "; out[3][0] = " mean ";
		out[4][0] = " std "; out[5][0] = " min "; out[6][0] = " 25% ";
		out[7][0] = " 50% "; out[8][0] = " 75% "; out[9][0] = " max ";
		for(int c=0; c<clen; c++) {
			int cc = c+1;
			out[0][cc] = titles[c]+" ";
			out[1][cc] = types[c].name();
			switch(types[c]) {
				case BOOL:
					out[3][cc] = " "; out[4][cc] = " "; //mean and std
					out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " "; out[9][cc] = " ";
					int olen = 0;
					for(boolean o: bool_column.get(titles[c])) if(o) { olen++; }
					out[2][cc] = " "+olen+" ";
					break;
				case BYTE: double bmean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " --- "; out[9][cc] = " NaN ";
					List<Byte> bytes = new ArrayList<Byte>();
					for(byte b: byte_column.get(titles[c])) if(Byte.MIN_VALUE!=b) { bmean += b; bytes.add(b); }
					int blen = bytes.size();
					out[2][cc] = " "+blen+" ";
					if(blen==0) break;
					bmean /= blen;
					byte bmean_b = (byte) ( (int)(bmean+0.5d) - (bmean<-0.5d ? 1 : 0) );
					out[3][cc] = " "+bmean_b+" ";
					double bstd = 0d;
					for(byte b: bytes) bstd += (b-bmean)*(b-bmean);
					if(blen==1) { bstd = Byte.MIN_VALUE; } else { bstd /= (blen-1); }
					byte bstd_b = (byte) (bstd+0.99d);
					out[4][cc] = " "+bstd_b+" ";
					Byte[] a_b = bytes.toArray(new Byte[0]);
					Arrays.sort(a_b);
					out[5][cc] = " "+a_b[0]+" "; out[9][cc] = " "+a_b[blen-1]+" ";
					if(blen>3) out[6][cc] = " "+a_b[blen/4]+" "; if(blen>1) out[7][cc] = " "+a_b[(blen+1)/2]+" "; if(blen>3) out[8][cc] = " "+a_b[(3*blen)/4]+" ";
					break;
				case INT: double imean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " --- "; out[9][cc] = " NaN ";
					List<Integer> integers = new ArrayList<Integer>();
					for(int i: int_column.get(titles[c])) if(Integer.MIN_VALUE!=i) { imean += i; integers.add(i); }
					int ilen = integers.size();
					out[2][cc] = " "+ilen+" ";
					if(ilen==0) break;
					imean /= ilen;
					int imean_i = (int) (imean+0.5d) - (imean<-0.5d ? 1 : 0);
					out[3][cc] = " "+imean_i+" ";
					double istd = 0d;
					for(int i: integers) istd += (i-imean)*(i-imean);
					if(ilen==1) { istd = Integer.MIN_VALUE;; } else { istd /= (ilen-1); }
					int istd_i = (int) (istd+0.99d);
					out[4][cc] = " "+istd_i+" ";
					Integer[] a_i = integers.toArray(new Integer[0]);
					Arrays.sort(a_i);
					out[5][cc] = " "+a_i[0]+" "; out[9][cc] = " "+a_i[ilen-1]+" ";
					if(ilen>3) out[6][cc] = " "+a_i[ilen/4]+" "; if(ilen>1) out[7][cc] = " "+a_i[(ilen+1)/2]+" "; if(ilen>3) out[8][cc] = " "+a_i[(3*ilen)/4]+" ";
					break;
				case SHORT: double rmean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " --- "; out[9][cc] = " NaN ";
					List<Short> shorts = new ArrayList<Short>();
					for(short r: short_column.get(titles[c])) if(Short.MIN_VALUE!=r) { rmean += r; shorts.add(r); }
					int rlen = shorts.size();
					out[2][cc] = " "+rlen+" ";
					if(rlen==0) break;
					rmean /= rlen;
					short rmean_s = (short) ( (int)(rmean+0.5d) - (rmean<-0.5d ? 1 : 0) );
					out[3][cc] = " "+rmean_s+" ";
					double rstd = 0d;
					for(short r: shorts) rstd += (r-rmean)*(r-rmean);
					if(rlen==1) { rstd = Integer.MIN_VALUE;; } else { rstd /= (rlen-1); }
					short rstd_s = (short) (rstd+0.99d);
					out[4][cc] = " "+rstd_s+" ";
					Short[] a_r = shorts.toArray(new Short[0]);
					Arrays.sort(a_r);
					out[5][cc] = " "+a_r[0]+" "; out[9][cc] = " "+a_r[rlen-1]+" ";
					if(rlen>3) out[6][cc] = " "+a_r[rlen/4]+" "; if(rlen>1) out[7][cc] = " "+a_r[(rlen+1)/2]+" "; if(rlen>3) out[8][cc] = " "+a_r[(3*rlen)/4]+" ";
					break;
				case LONG: double lmean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " --- "; out[9][cc] = " NaN ";
					List<Long> longs = new ArrayList<Long>();
					for(long l: long_column.get(titles[c])) if(Long.MIN_VALUE!=l) { lmean += l; longs.add(l); }
					int llen = longs.size();
					out[2][cc] = " "+llen+" ";
					if(llen==0) break;
					lmean /= llen;
					long imean_l = (long) (lmean+0.5d) - (lmean<-0.5d ? 1L : 0L);
					out[3][cc] = " "+imean_l+" ";
					double lstd = 0d;
					for(long l: longs) lstd += (l-lmean)*(l-lmean);
					if(llen==1) { lstd = Long.MIN_VALUE; } else { lstd /= (llen-1); }
					long istd_l = (long) (lstd+0.99d);
					out[4][cc] = " "+istd_l+" ";
					Long[] a_l = longs.toArray(new Long[0]);
					Arrays.sort(a_l);
					out[5][cc] = " "+a_l[0]+" "; out[9][cc] = " "+a_l[llen-1]+" ";
					if(llen>3) out[6][cc] = " "+a_l[llen/4]+" "; if(llen>1) out[7][cc] = " "+a_l[(llen+1)/2]+" "; if(llen>3) out[8][cc] = " "+a_l[(3*llen)/4]+" ";
					break;
				case FLOAT: float fmean = 0f;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " NaN "; out[7][cc] = " NaN "; out[8][cc] = " NaN "; out[9][cc] = " NaN ";
					List<Float> floats = new ArrayList<Float>();
					for(float f: float_column.get(titles[c])) if(!Float.isNaN(f)) { fmean += f; floats.add(f); }
					int flen = floats.size();
					out[2][cc] = " "+flen+" ";
					if(flen==0) break;
					fmean /= flen;
					out[3][cc] = " "+fmean+" ";
					float fstd = 0f;
					for(float f: floats) fstd += (f-fmean)*(f-fmean);
					if(flen==1) { fstd = Float.POSITIVE_INFINITY; } else { fstd /= (flen-1f); }
					out[4][cc] = " "+fstd+" ";
					Float[] a_f = floats.toArray(new Float[0]);
					Arrays.sort(a_f);
					out[5][cc] = " "+a_f[0]+" "; out[9][cc] = " "+a_f[flen-1]+" ";
					if(flen>3) out[6][cc] = " "+a_f[flen/4]+" "; if(flen>1) out[7][cc] = " "+a_f[(flen+1)/2]+" "; if(flen>3) out[8][cc] = " "+a_f[(3*flen)/4]+" ";
					break;
				case DOUBLE: double dmean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " NaN "; out[7][cc] = " NaN "; out[8][cc] = " NaN "; out[9][cc] = " NaN ";
					List<Double> doubles = new ArrayList<Double>();
					for(double d: double_column.get(titles[c])) if(!Double.isNaN(d)) { dmean += d; doubles.add(d); }
					int dlen = doubles.size();
					out[2][cc] = " "+dlen+" ";
					if(dlen==0) break;
					dmean /= dlen;
					out[3][cc] = " "+dmean+" ";
					double dstd = 0d;
					for(double d: doubles) dstd += (d-dmean)*(d-dmean);
					if(dlen==1) { dstd = Double.POSITIVE_INFINITY; } else { dstd /= (dlen-1d); }
					out[4][cc] = " "+dstd+" ";
					Double[] a_d = doubles.toArray(new Double[0]);
					Arrays.sort(a_d);
					out[5][cc] = " "+a_d[0]+" "; out[9][cc] = " "+a_d[dlen-1]+" ";
					if(dlen>3) out[6][cc] = " "+a_d[dlen/4]+" "; if(dlen>1) out[7][cc] = " "+a_d[(dlen+1)/2]+" "; if(dlen>3) out[8][cc] = " "+a_d[(3*dlen)/4]+" ";
					break;
				case STRING:
					out[3][cc] = " "; out[4][cc] = " "; //mean and std
					out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " "; out[9][cc] = " ";
					int slen = 0;
					for(String s: string_column.get(titles[c])) if(!(s.equals(" ") || s.length()<1)) { slen++; }
					out[2][cc] = " "+slen+" ";
					break;
				default:
					out[3][cc] = " "; out[4][cc] = " "; //mean and std
					out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " "; out[9][cc] = " ";
			}
		}
		boolean[] coldivs = {true,false,false,false,false,false,false,false,false};
		boolean[] rowdivs = new boolean[clen]; rowdivs[0] = true;
		for(int r=1; r<clen; r++) rowdivs[r] = false;
		FormatHelper.printTable(out, coldivs, rowdivs, true);
	}
	public void describeDimension() {
		String dimstr = "(1/1): "+dimension_name+" ["+dimension.length+"] ";
		double[] rng = StdAnalysis.minmax(dimension);
		double[] std = StdAnalysis.mean_var(dimension);
		dimstr += "{"+rng[0]+" ... "+rng[1]+", mean="+std[0]+", var="+std[1]+"}";
		System.out.println(dimstr);
	}
	public void clear() {
		titles = null;
		types = null;
		titles = new String[0];
		types = new DataType[0];
		minmax_mean_sill = new double[4][0];
		datalength = 0;
		bool_column.clear();
		byte_column.clear();
		int_column.clear();
		short_column.clear();
		long_column.clear();
		float_column.clear();
		double_column.clear();
		string_column.clear();
		dimension = new double[0];
		dimension_name = "";
		attribs_dim.clear();
	}
	
	private String[] null_line(String _in, String _del) {
		int count = 1;
		int occ = _in.indexOf(_del);
		String temp = ""+_in;
		while(occ>=0) { count++; temp = temp.substring(occ+_del.length()); occ = temp.indexOf(_del); }
		String[] out = new String[count];
		for(int c=0; c<count; c++) out[c] = "";
//		System.out.println("[DEBUG] found "+(count-1)+" commas -> new String["+count+"]");
		return out;
	}
	private String[] break_line(String _in, String _del) {
		List<String> parts = new ArrayList<String>();
		String temp = ""+_in;
		int occ = temp.indexOf(_del);
		while(occ>=0) { parts.add(temp.substring(0, occ).trim()); temp = temp.substring(occ+_del.length()); occ = temp.indexOf(_del); }
		parts.add(temp.trim());
//		String debug = ""+parts.get(0);
//		for(int s=1; s<parts.size(); s++) debug +=" | "+parts.get(s);
//		System.out.println("[DEBUG] break_line: "+debug);
		return parts.toArray(new String[0]);
	}
	private boolean set_boolean(String _s) {
		boolean b = Boolean.FALSE;
		try { b = Boolean.parseBoolean(_s); } catch(NullPointerException | NumberFormatException np_nf_e) { b = Boolean.FALSE; }
		return b;
	}
	private byte set_byte(String _s) {
		byte b = Byte.MIN_VALUE;
		try { b = Byte.parseByte(_s); } catch(NumberFormatException nfe) { b = Byte.MIN_VALUE; }
		return b;
	}
	private int set_int(String _s) {
		int i = Integer.MIN_VALUE;
		try { i = Integer.parseInt(_s); } catch(NumberFormatException nfe) { i = Integer.MIN_VALUE; }
		return i;
	}
	private short set_short(String _s) {
		short s = Short.MIN_VALUE;
		try { s = Short.parseShort(_s); } catch(NumberFormatException nfe) { s = Short.MIN_VALUE; }
		return s;
	}
	private long set_long(String _s) {
		long l = Integer.MIN_VALUE;
		try { l = Long.parseLong(_s); } catch(NumberFormatException nfe) { l = Long.MIN_VALUE; }
		return l;
	}
	private float set_float(String _s) {
		float f = Float.NaN;
		try { f = Float.parseFloat(_s); } catch(NumberFormatException nfe) { f = Float.NaN; }
		return f;
	}
	private double set_double(String _s) {
		double d = Double.NaN;
		try { d = Double.parseDouble(_s); } catch(NumberFormatException nfe) { d = Double.NaN; }
		return d;
	}
	protected void add_variable(String _name, DataType _type) {
		String[] old_tit = new String[titles.length];
		DataType[] old_tp = new DataType[types.length];
		for(int t=0; t<titles.length; t++) { old_tit[t] = titles[t]; old_tp[t] = types[t]; }
		titles = null; titles = new String[old_tit.length+1];
		types = null; types = new DataType[old_tp.length+1];
		for(int t=0; t<old_tit.length; t++) { titles[t] = old_tit[t]; types[t] = old_tp[t]; }
		titles[old_tit.length] = _name;
		types[old_tp.length] = _type;
	}
	protected void add_stdana(double _min, double _max, double _mean, double _sill) {
		int dlen = minmax_mean_sill[0].length;
		double[][] old_mas = new double[dlen][4];
		for(int t=0; t<dlen; t++) {
			old_mas[t][0] = minmax_mean_sill[0][t]; old_mas[t][1] = minmax_mean_sill[1][t];
			old_mas[t][2] = minmax_mean_sill[2][t]; old_mas[t][3] = minmax_mean_sill[3][t];
		}
		minmax_mean_sill = null; minmax_mean_sill = new double[4][dlen+1];
		for(int t=0; t<dlen; t++) {
			minmax_mean_sill[0][t] = old_mas[t][0]; minmax_mean_sill[1][t] = old_mas[t][1];
			minmax_mean_sill[2][t] = old_mas[t][2]; minmax_mean_sill[3][t] = old_mas[t][3];
		}
		minmax_mean_sill[0][dlen] = _min;
		minmax_mean_sill[1][dlen] = _max;
		minmax_mean_sill[2][dlen] = _mean;
		minmax_mean_sill[3][dlen] = _sill;
	}
	protected void set_stdana(int _var_id, double _min, double _max, double _mean, double _sill) {
		minmax_mean_sill[0][_var_id] = _min;
		minmax_mean_sill[1][_var_id] = _max;
		minmax_mean_sill[2][_var_id] = _mean;
		minmax_mean_sill[3][_var_id] = _sill;
	}
	protected void recalcStats(int _var_id) {
		String vn = getVarname(_var_id);
		switch(getVariableType(_var_id)) {
			case BYTE: byte[] arrb = byte_column.get(vn);
				byte[] inaxb = StdAnalysis.minmax(arrb); byte[] mvb = StdAnalysis.mean_var(arrb);
				set_stdana(_var_id, inaxb[0], inaxb[1], mvb[0], mvb[1]);
				break;
			case SHORT: short[] arrs = short_column.get(vn);
				short[] inaxs = StdAnalysis.minmax(arrs); short[] mvs = StdAnalysis.mean_var(arrs);
				set_stdana(_var_id, inaxs[0], inaxs[1], mvs[0], mvs[1]);
				break;
			case INT: int[] arri = int_column.get(vn);
				int[] inaxi = StdAnalysis.minmax(arri); int[] mvi = StdAnalysis.mean_var(arri);
				set_stdana(_var_id, inaxi[0], inaxi[1], mvi[0], mvi[1]);
				break;
			case LONG: long[] arrl = long_column.get(vn);
				long[] inaxl = StdAnalysis.minmax(arrl); long[] mvl = StdAnalysis.mean_var(arrl);
				set_stdana(_var_id, inaxl[0], inaxl[1], mvl[0], mvl[1]);
				break;
			case FLOAT: float[] arrf = float_column.get(vn);
				float[] inaxf = StdAnalysis.minmax(arrf); float[] mvf = StdAnalysis.mean_var(arrf);
				set_stdana(_var_id, inaxf[0], inaxf[1], mvf[0], mvf[1]);
				break;
			case DOUBLE: double[] arrd = double_column.get(vn);
				double[] inaxd = StdAnalysis.minmax(arrd); double[] mvd = StdAnalysis.mean_var(arrd);
				set_stdana(_var_id, inaxd[0], inaxd[1], mvd[0], mvd[1]);
				break;
			default:
				break;
		}
	}
	private String createTitle(String[] _all_titles) {
		int tnum = 1; int clen = _all_titles.length;
		String test_title = "Column_"+Integer.toHexString(tnum);
		boolean exists = true;
		while(exists) {
			exists = false;
			for(int tt=0; tt<clen; tt++) {
				if(_all_titles[tt]==null) continue;
				if(_all_titles[tt].equals(test_title)) { exists = true; break; }
			}
			tnum++;
		}
		return test_title;
	}
	protected void createDimension() {
		dimension = new double[datalength];
		for(int n=0; n<datalength; n++)
			dimension[n] = n+1d;
		boolean isContained = true;
		while(isContained) {
			dimension_name = "dim"+Math.abs(Constants.RANDOM.nextInt(Integer.MAX_VALUE>>2));
			isContained = hasVariable(dimension_name);
		}
		attribs_dim.clear();
	}
	
	@SuppressWarnings("rawtypes")
	private static void setupDimensions(NetcdfFormatWriter.Builder builder, Map<String, Dimension> dimensions, Map<String, String> variables,
			DataFrame df1, DataFrame2D df2, DataFrame3D df3) {
		int dimCount = df1!=null?1 : df2!=null?2 : 3;
		for(int dim = 0; dim < dimCount; dim++) {
			int getter = dim+Constants.FIRST_IDX;
			String nam = df1!=null ? df1.getDimensionName() :
						 df2!=null ? df2.getDimensionName(getter) :
							 df3.getDimensionName(getter);
			if (!dimensions.containsKey(nam)) {
				System.out.println("  ... create dimension {" + nam + "}");
				int len = 	df1!=null ? df1.getNumberOfDatapoints() :
							df2!=null ? df2.getDimensionValues(getter).length :
								df3.getDimensionValues(getter).length;
				dimensions.put(nam, builder.addDimension(nam, len));
				Variable.Builder vv = builder.addVariable(nam, ucar.ma2.DataType.DOUBLE, nam);
				variables.put(nam, "empty");
				Map<String, String> attribs = 
						df1!=null ? df1.getDimensionAttributes() :
						df2!=null ? df2.getAttributesFromDimension(getter) :
							df3.getAttributesFromDimension(getter);
				for (String a : attribs.keySet()) {
					if (a.equals("_FillValue")) {
						vv.addAttribute(new Attribute(a, Double.parseDouble(attribs.get(a)))); continue;
					}
					vv.addAttribute(new Attribute(a, attribs.get(a)));
				}
			} else {
				System.err.println("[WARNING] dimension \""+nam+"\" allready exist in other DataFrame(2D/3D) object.");
			}
		}
	}
	@SuppressWarnings("rawtypes")
	private static void setupVariables(NetcdfFormatWriter.Builder builder, Map<String, Dimension> dimensions, Map<String, String> variables,
			DataFrame df, DataFrame2D df2, DataFrame3D df3) {
		List<Dimension> dims = new ArrayList<Dimension>();
		int varCount = 0;
		if(df!=null) {
			dims.add((Dimension)dimensions.get(df.getDimensionName()));
			varCount = df.getVariableCount();
		}
		if(df2!=null) {
			dims.add((Dimension)dimensions.get(df2.getDimensionName(0+Constants.FIRST_IDX)));
			dims.add((Dimension)dimensions.get(df2.getDimensionName(1+Constants.FIRST_IDX)));
			varCount = df2.getVariableCount();
		}
		if(df3!=null) {
			dims.add((Dimension)dimensions.get(df3.getDimensionName(0+Constants.FIRST_IDX)));
			dims.add((Dimension)dimensions.get(df3.getDimensionName(1+Constants.FIRST_IDX)));
			dims.add((Dimension)dimensions.get(df3.getDimensionName(2+Constants.FIRST_IDX)));
			varCount = df3.getVariableCount();
		}
		for (int iv = 0; iv < varCount; iv++) {
			String nam = 
					df!=null ? df.getVarname(iv+Constants.FIRST_IDX) :
					df2!=null ? df2.getVarname(iv+Constants.FIRST_IDX) :
						df3.getVarname(iv+Constants.FIRST_IDX);
			if (!variables.containsKey(nam)) {
				System.out.println("  ... create variable <" + nam + ">");
				ucar.ma2.DataType dataType = null;
				DataFrame.DataType vartype = 
						df!=null ? df.getVariableType(iv+Constants.FIRST_IDX) :
						df2!=null ? df2.getVariableType(iv+Constants.FIRST_IDX) :
							df3.getVariableType(iv+Constants.FIRST_IDX);
				if(vartype==DataType.STRUCT) {
					if(df!=null)  df.defineStructure(builder, nam, dims, variables, dimensions);
					if(df2!=null) df2.defineStructure(builder, nam, dims, variables, dimensions);
					if(df3!=null) df3.defineStructure(builder, nam, dims, variables, dimensions);
					continue;
				}
				switch (vartype) {
					case BOOL:   dataType = ucar.ma2.DataType.BOOLEAN; break;
					case BYTE:   dataType = ucar.ma2.DataType.BYTE; break;
					case SHORT:  dataType = ucar.ma2.DataType.SHORT; break;
					case INT:    dataType = ucar.ma2.DataType.INT; break;
					case LONG:   dataType = ucar.ma2.DataType.LONG; break;
					case FLOAT:  dataType = ucar.ma2.DataType.FLOAT; break;
					case DOUBLE: dataType = ucar.ma2.DataType.DOUBLE; break;
					default: dataType = ucar.ma2.DataType.STRING; break;
				}
				Variable.Builder vars2 = builder.addVariable(nam, dataType, dims);
				variables.put(nam, "empty");
				switch(vartype) {
					case SHORT:  vars2.addAttribute(new Attribute("_FillValue", Short.MIN_VALUE)); break;
					case INT:    vars2.addAttribute(new Attribute("_FillValue", Integer.MIN_VALUE)); break;
					case LONG:   vars2.addAttribute(new Attribute("_FillValue", Long.MIN_VALUE)); break;
					case FLOAT:  vars2.addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_F)); break;
					case DOUBLE: vars2.addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_D)); break;
					default: break;
				}
			}
		}
	}
	private static void writeDim(NetcdfFormatWriter writer, Map<String, String> variables,
			DataFrame df1, DataFrame2D df2, DataFrame3D df3) throws IOException, InvalidRangeException {
		int dimCount = df1!=null?1 : df2!=null?2 : 3;
		for(int dim = 0; dim < dimCount; dim++) {
			int getter = dim + Constants.FIRST_IDX;
			String dimname = df1!=null ? df1.getDimensionName() :
							 df2!=null ? df2.getDimensionName(getter) :
								 df3.getDimensionName(getter);
			if (!variables.containsKey(dimname)) {
				System.err.println("[ERROR] found unprepared variable! It will not be written to the Netcdf-file!");
			}
			else if (!((String)variables.get(dimname)).equals("full")) {
				double[] double_source = df1!=null ? df1.getDimensionValues() :
										 df2!=null ? df2.getDimensionValues(getter) :
											 df3.getDimensionValues(getter);
				ArrayDouble.D1 double_arr = new ArrayDouble.D1(double_source.length);
				for(int dl = 0; dl < double_source.length; dl++ )
					double_arr.set(dl, double_source[dl]);
				writer.write(writer.findVariable(dimname), double_arr);
				variables.put(dimname, "full");
			} else {
//				System.err.println("[WARNING] Dimension \""+dimname+"\" has already been written to Netcdf file from other source DataFrame(2D/3D).");
			}
		}
	}
	private static void write1D(NetcdfFormatWriter writer, Map<String, String> variables, DataFrame df) throws IOException, InvalidRangeException {
		for (int iv = 0; iv < df.getVariableCount(); iv++) {
			DataType vartype = df.getVariableType(iv+Constants.FIRST_IDX);
			String varname = df.getVarname(iv+Constants.FIRST_IDX);
			if (!variables.containsKey(varname)) {
				System.err.println("[ERROR] found unprepared variable! It will not be written to the Netcdf-file!");
			}
			else if (!((String)variables.get(varname)).equals("full")) {
				int datalength = df.getNumberOfDatapoints();
				switch(vartype) {
					case BOOL: ArrayBoolean.D1 bool_arr = new ArrayBoolean.D1(datalength); boolean[] bool_source = (boolean[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { bool_arr.set(dl, bool_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), bool_arr); break;
					case BYTE: ArrayByte.D1 byte_arr = new ArrayByte.D1(datalength, false); byte[] byte_source = (byte[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { byte_arr.set(dl, byte_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), byte_arr); break;
					case SHORT: ArrayShort.D1 short_arr = new ArrayShort.D1(datalength, false); short[] short_source = (short[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { short_arr.set(dl, short_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), short_arr); break;
					case INT: ArrayInt.D1 int_arr = new ArrayInt.D1(datalength, false); int[] int_source = (int[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { int_arr.set(dl, int_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), int_arr); break;
					case LONG: ArrayLong.D1 long_arr = new ArrayLong.D1(datalength, false); long[] long_source = (long[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { long_arr.set(dl, long_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), long_arr); break;
					case FLOAT: ArrayFloat.D1 float_arr = new ArrayFloat.D1(datalength); float[] float_source = (float[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { float_arr.set(dl, float_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), float_arr); break;
					case DOUBLE: ArrayDouble.D1 double_arr = new ArrayDouble.D1(datalength); double[] double_source = (double[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { double_arr.set(dl, double_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), double_arr); break;
					case STRING: ArrayString.D1 string_arr = new ArrayString.D1(datalength); String[] string_source = (String[])df.getArray(varname);
						for(int dl = 0; dl < datalength; ) { string_arr.set(dl, string_source[dl]); dl++; }
						writer.write(writer.findVariable(varname), string_arr); break;
					case STRUCT: df.writeStruct(writer, varname); break;
					default: break;
				}
				variables.put(varname, "full");
			}
		}
	}
	private static void write2D(NetcdfFormatWriter writer, Map<String, String> variables, DataFrame2D df2) throws IOException, InvalidRangeException {
		for (int iv = 0; iv < df2.getVariableCount(); iv++) {
			DataType vartype = df2.getVariableType(iv+Constants.FIRST_IDX);
			String varname = df2.getVarname(iv+Constants.FIRST_IDX);
			if (!variables.containsKey(varname)) {
				System.err.println("[ERROR] found unprepared variable! It will not be written to the Netcdf-file!");
			} else
			if (!((String)variables.get(varname)).equals("full")) {
				int[] datalength = { df2.getDimensionValues(0+Constants.FIRST_IDX).length, df2.getDimensionValues(1+Constants.FIRST_IDX).length };
				switch(vartype) {
					case BOOL: ArrayBoolean.D2 bool_arr = new ArrayBoolean.D2(datalength[0], datalength[1]);
						boolean[][] bool_source = (boolean[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							bool_arr.set(dl0, dl1, bool_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), bool_arr); break;
					case BYTE: ArrayByte.D2 byte_arr = new ArrayByte.D2(datalength[0], datalength[1], false);
						byte[][] byte_source = (byte[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							byte_arr.set(dl0, dl1, byte_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), byte_arr); break;
					case SHORT: ArrayShort.D2 short_arr = new ArrayShort.D2(datalength[0], datalength[1], false);
						short[][] short_source = (short[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							short_arr.set(dl0, dl1, short_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), short_arr); break;
					case INT: ArrayInt.D2 int_arr = new ArrayInt.D2(datalength[0], datalength[1], false);
						int[][] int_source = (int[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							int_arr.set(dl0, dl1, int_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), int_arr); break;
					case LONG: ArrayLong.D2 long_arr = new ArrayLong.D2(datalength[0], datalength[1], false);
						long[][] long_source = (long[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							long_arr.set(dl0, dl1, long_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), long_arr); break;
					case FLOAT: ArrayFloat.D2 float_arr = new ArrayFloat.D2(datalength[0], datalength[1]);
						float[][] float_source = (float[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							float_arr.set(dl0, dl1, Float.isNaN(float_source[dl0][dl1]) ? Constants.FILL_VALUE_F : float_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), float_arr); break;
					case DOUBLE: ArrayDouble.D2 double_arr = new ArrayDouble.D2(datalength[0], datalength[1]);
						double[][] double_source = (double[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							double_arr.set(dl0, dl1, Double.isNaN(double_source[dl0][dl1]) ? Constants.FILL_VALUE_D : double_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), double_arr); break;
					case STRING: ArrayString.D2 string_arr = new ArrayString.D2(datalength[0], datalength[1]);
						String[][] string_source = (String[][])df2.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++)
							string_arr.set(dl0, dl1, string_source[dl0][dl1]);
						writer.write(writer.findVariable(varname), string_arr); break;
					case STRUCT:
						df2.writeStruct(writer, varname);
						break;
					default: break;
				}
				variables.put(varname, "full");
			}
		}
	}
	private static void write3D(NetcdfFormatWriter writer, Map<String, String> variables, DataFrame3D df3) throws IOException, InvalidRangeException {
		for (int iv = 0; iv < df3.getVariableCount(); iv++) {
			DataType vartype = df3.getVariableType(iv+Constants.FIRST_IDX);
			String varname = df3.getVarname(iv+Constants.FIRST_IDX);
			if (!variables.containsKey(varname)) {
				System.err.println("[ERROR] found unprepared variable! It will not be written to the Netcdf-file!");
			} else
			if (!((String)variables.get(varname)).equals("full")) {
				int[] datalength = { df3.getDimensionValues(0+Constants.FIRST_IDX).length,
						df3.getDimensionValues(1+Constants.FIRST_IDX).length,
						df3.getDimensionValues(2+Constants.FIRST_IDX).length};
				switch(vartype) {
					case BOOL: ArrayBoolean.D3 bool_arr = new ArrayBoolean.D3(datalength[0], datalength[1], datalength[2]);
						boolean[][][] bool_source = (boolean[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							bool_arr.set(dl0, dl1, dl2, bool_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), bool_arr); break;
					case BYTE: ArrayByte.D3 byte_arr = new ArrayByte.D3(datalength[0], datalength[1], datalength[2], false);
						byte[][][] byte_source = (byte[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							byte_arr.set(dl0, dl1, dl2, byte_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), byte_arr); break;
					case SHORT: ArrayShort.D3 short_arr = new ArrayShort.D3(datalength[0], datalength[1], datalength[2], false);
						short[][][] short_source = (short[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							short_arr.set(dl0, dl1, dl2, short_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), short_arr); break;
					case INT: ArrayInt.D3 int_arr = new ArrayInt.D3(datalength[0], datalength[1], datalength[2], false);
						int[][][] int_source = (int[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							int_arr.set(dl0, dl1, dl2, int_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), int_arr); break;
					case LONG: ArrayLong.D3 long_arr = new ArrayLong.D3(datalength[0], datalength[1], datalength[2], false);
						long[][][] long_source = (long[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							long_arr.set(dl0, dl1, dl2, long_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), long_arr); break;
					case FLOAT: ArrayFloat.D3 float_arr = new ArrayFloat.D3(datalength[0], datalength[1], datalength[2]);
						float[][][] float_source = (float[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							float_arr.set(dl0, dl1, dl2, Float.isNaN(float_source[dl0][dl1][dl2]) ? Constants.FILL_VALUE_F : float_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), float_arr); break;
					case DOUBLE: ArrayDouble.D3 double_arr = new ArrayDouble.D3(datalength[0], datalength[1], datalength[2]);
						double[][][] double_source = (double[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							double_arr.set(dl0, dl1, dl2, Double.isNaN(double_source[dl0][dl1][dl2]) ? Constants.FILL_VALUE_D : double_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), double_arr); break;
					case STRING: ArrayString.D3 string_arr = new ArrayString.D3(datalength[0], datalength[1], datalength[2]);
						String[][][] string_source = (String[][][])df3.getArray(varname);
						for(int dl0 = 0; dl0 < datalength[0]; dl0++) for(int dl1 = 0; dl1 < datalength[1]; dl1++) for(int dl2=0; dl2<datalength[2]; dl2++)
							string_arr.set(dl0, dl1, dl2, string_source[dl0][dl1][dl2]);
						writer.write(writer.findVariable(varname), string_arr); break;
					case STRUCT:
						df3.writeStruct(writer, varname);
						break;
					default: break;
				}
				variables.put(varname, "full");
			}
		}
	}
	
	/**
	 * struct variables not supported natively by JKriging<br>
	 * this method has to be override, when use<br><br>
	 * implement what size specific struct variables use.
	 * 
	 * @param builder    used builder
	 * @param varname    name of the struct variable
	 * @param dims       map of collected dimensions
	 * @param variables  map of collected variables 
	 * @param dimensions used dimension(s) for this struct variable
	 */
	protected void defineStructure(NetcdfFormatWriter.Builder builder, String varname, List<Dimension> dims, Map<String,String> variables, Map<String,Dimension> dimensions) {
		System.out.println("struct variables aren't supported.");
		//do nothing
	}
	/**
	 * struct variables not supported natively by JKriging<br>
	 * this method has to be overridden, when use<br><br>
	 * implement how struct variables are written to netcdf files
	 * 
	 * @param writer  the netcdf writer
	 * @param varname name of the struct variable
	 */
	protected void writeStruct(NetcdfFormatWriter writer, String varname) {
		System.out.println("Writing struct variable \""+varname+"\" is unsupported.");
	}
	
	public enum DataType {
		BOOL("bool"), BYTE("byte"), INT("int"), SHORT("short"), LONG("long"),
		FLOAT("float"), DOUBLE("double"),
		STRING("string"), STRUCT("struct");
		
		private String name;
		
		private DataType(String _name) {
			name = _name;
		}
		
		@Override
		public String toString() {
			return name;
		}
		public static DataType getDataType(String _in, DataType _default) {
			String _t = _in.toLowerCase();
			if(_t.equals("boolean") || _t.equals("logical") || _t.equals("bool")) return BOOL;
			if(_t.equals("byte")    || _t.equals("bt")  || _t.equals("b")) return BYTE;
			if(_t.equals("integer") || _t.equals("int") || _t.equals("i")) return INT;
			if(_t.equals("short")   || _t.equals("sht")) return SHORT;
			if(_t.equals("long")    || _t.equals("lng") || _t.equals("l")) return LONG;
			if(_t.equals("float")   || _t.equals("flt") || _t.equals("f")) return FLOAT;
			if(_t.equals("double")  || _t.equals("dbl") || _t.equals("d")) return DOUBLE;
			if(_t.equals("string")  || _t.equals("str") || _t.equals("s")) return STRING;
			if(_t.equals("struct")  || _t.equals("stc")) return STRUCT;
			return _default;
		}
	}
}
