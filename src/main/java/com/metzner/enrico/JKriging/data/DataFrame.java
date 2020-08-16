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
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.NetcdfFileWriter.Version;
import ucar.nc2.Variable;

public class DataFrame {

	private int datalength;
	private DataType default_data_type;
	private String[] titles;
	private DataType[] types;
	private double[][] minmax_mean_sill;
	private Map<String, boolean[]> bool_column;
	private Map<String, byte[]>    byte_column;
	private Map<String, int[]>     int_column;
	private Map<String, short[]>   short_column;
	private Map<String, long[]>    long_column;
	private Map<String, float[]>   float_column;
	private Map<String, double[]>  double_column;
	private Map<String, String[]>  string_column;




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
	}




	public void addColumn(String _column_name, boolean[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, boolean[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.BOOL);
		add_stdana(Double.NaN,Double.NaN,Double.NaN,Double.NaN);
		boolean[] new_data = new boolean[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = false; }
		bool_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, byte[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, byte[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.BYTE);
		byte[] inax = StdAnalysis.minmax(_column_data); byte[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		byte[] new_data = new byte[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Byte.MIN_VALUE; }
		byte_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, int[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, int[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.INT);
		int[] inax = StdAnalysis.minmax(_column_data); int[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		int[] new_data = new int[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Integer.MIN_VALUE; }
		int_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, short[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, short[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.SHORT);
		short[] inax = StdAnalysis.minmax(_column_data); short[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		short[] new_data = new short[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Short.MIN_VALUE; }
		short_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, long[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, long[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.LONG);
		long[] inax = StdAnalysis.minmax(_column_data); long[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		long[] new_data = new long[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Long.MIN_VALUE; }
		long_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, float[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, float[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.FLOAT);
		float[] inax = StdAnalysis.minmax(_column_data); float[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		float[] new_data = new float[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Float.NaN; }
		float_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, double[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, double[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.DOUBLE);
		double[] inax = StdAnalysis.minmax(_column_data); double[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		double[] new_data = new double[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = Double.NaN; }
		double_column.put(_column_name, new_data);
	}
	public void addColumn(String _column_name, String[] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, String[] _column_data, boolean chopORfill) {
		if(hasVariable(_column_name)) {
			System.err.println("A column with title <"+_column_name+"> allready exists in this dataframe!");
			return;
		}
		if(_column_data.length!=datalength && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) datalength = _column_data.length;
		add_variable(_column_name, DataType.STRING);
		add_stdana(Double.NaN,-Double.NaN,Double.NaN,Double.NaN);
		String[] new_data = new String[datalength];
		for(int t=0; t<datalength; t++) if(t<_column_data.length) { new_data[t] = _column_data[t]; } else { new_data[t] = ""; }
		string_column.put(_column_name, new_data);
	}

	public void renameVariable(int _old_var_id, String _new_name) {
		int ovi = _old_var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not rename variable with ID "+_old_var_id+", could not find variable!");
			return;
		}
		String _old_var_name = getVarname(_old_var_id);
		if(hasVariable(_new_name)) {
			System.err.println("Could not rename variable \""+_old_var_name+"\", a variable by name \""+_new_name+"\" already exist!");
			//DataHelper.printStackTrace(System.err);
			return;
		}
		switch(types[ovi]) {
			case BOOL:   addColumn(_new_name, bool_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case BYTE:   addColumn(_new_name, byte_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case SHORT:  addColumn(_new_name, short_column.get(_old_var_name));  removeVariable(_old_var_name); break;
			case INT:    addColumn(_new_name, int_column.get(_old_var_name));    removeVariable(_old_var_name); break;
			case LONG:   addColumn(_new_name, long_column.get(_old_var_name));   removeVariable(_old_var_name); break;
			case FLOAT:  addColumn(_new_name, float_column.get(_old_var_name));  removeVariable(_old_var_name); break;
			case DOUBLE: addColumn(_new_name, double_column.get(_old_var_name)); removeVariable(_old_var_name); break;
			case STRING: addColumn(_new_name, string_column.get(_old_var_name)); removeVariable(_old_var_name); break;
			default: System.err.println("An unexpected error occured: Unknown Variable type!"); DataHelper.printStackTrace(System.err); break;
		}
	}
	public void renameVariable(String _old_name, String _new_name) {
		renameVariable(getVariableID(_old_name), _new_name);
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
					case BYTE: byte[] barr = new byte[float_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) float_arr[i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[float_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) float_arr[i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[float_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (int) float_arr[i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[float_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = (long) float_arr[i];
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
					case BYTE: byte[] barr = new byte[double_arr.length]; for(int i=0; i<barr.length; i++) barr[i] = (byte) double_arr[i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[] rarr = new short[double_arr.length]; for(int i=0; i<rarr.length; i++) rarr[i] = (short) double_arr[i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[] iarr = new int[double_arr.length]; for(int i=0; i<iarr.length; i++) iarr[i] = (int) double_arr[i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[] larr = new long[double_arr.length]; for(int i=0; i<larr.length; i++) larr[i] = (long) double_arr[i];
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
			default:
				System.err.println("An unexpected error occured, can not determin the datatype of the variable!");
				DataHelper.printStackTrace(System.err); break;
		}
	}
	
	public DataFrame concat(DataFrame df2) {
		boolean haveDifferentKeys = false;
		for(String tit: df2.allVariableNames())
			if(this.hasVariable(tit)) {
				if(this.getVariableType(tit)!=df2.getVariableType(tit)) {
					haveDifferentKeys = true; break;
				}
			} else {
				haveDifferentKeys = true; break;
			}
		if(!haveDifferentKeys)
			for(String tit: this.titles)
				if(!df2.hasVariable(tit)) {
					haveDifferentKeys = true; break;
				}
		if(haveDifferentKeys) {
			System.err.println("DataFrames have not the same variables, can't concat them!");
			return this;
		}
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
				default: break;
			}
		}
		datalength += df2.datalength;
		return this;
	}
	public DataFrame2D toDataFrame2D(String dim_one_name, String dim_two_name) {
		DataFrame2D res2D = new DataFrame2D();
		//check if the DataFrame contain these dimensions!
		int dim1ID = strings_index(titles, dim_one_name);
		if(dim1ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_one_name+"\"!");
			DataHelper.printStackTrace(System.err); return res2D; }
		int dim2ID = strings_index(titles, dim_two_name);
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
						Double e = new Double(bool_arr[dl] ? 1d : 0d);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case BYTE: byte[] byte_arr = byte_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(byte_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case SHORT: short[] short_arr = short_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(short_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case INT: int[] int_arr = int_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(int_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case LONG: long[] long_arr = long_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(long_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case FLOAT: float[] float_arr = float_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Float.isNaN(float_arr[dl])) { System.out.println("WARNING: dimension variable contains missing value"); continue;}
						Double e = new Double(float_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case DOUBLE: double[] double_arr = double_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Double.isNaN(double_arr[dl])) { System.out.println("WARNING: dimension variable contains missing values"); continue;}
						Double e = new Double(double_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case STRING:
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
		for(int v=0; v<dim_lengths[0]; v++) { List<Integer> lv = values_dim_one.get(new Double(dimension_one[v]));
			for(int u=0; u<dim_lengths[1]; u++) { List<Integer> lu = values_dim_two.get(new Double(dimension_two[u]));
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
		int dim1ID = strings_index(titles, dim_one_name);
		if(dim1ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_one_name+"\"!");
			DataHelper.printStackTrace(System.err); return res3D; }
		int dim2ID = strings_index(titles, dim_two_name);
		if(dim2ID<0) { System.err.println("DataFrame does not contain any variable \""+dim_two_name+"\"!");
			DataHelper.printStackTrace(System.err); return res3D; }
		int dim3ID = strings_index(titles, dim_three_name);
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
						Double e = new Double(bool_arr[dl] ? 1d : 0d);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case BYTE: byte[] byte_arr = byte_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(byte_arr[dl]);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case SHORT: short[] short_arr = short_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(short_arr[dl]);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case INT: int[] int_arr = int_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(int_arr[dl]);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case LONG: long[] long_arr = long_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						Double e = new Double(long_arr[dl]);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case FLOAT: float[] float_arr = float_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Float.isNaN(float_arr[dl])) { System.out.println("WARNING: dimension variable contains missing value"); continue;}
						Double e = new Double(float_arr[dl]);
						if ( to_fill.containsKey(e) ) { ((List<Integer>)to_fill.get(e)).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				case DOUBLE: double[] double_arr = double_column.get(titles[dimID]);
					for(int dl=0; dl<datalength; dl++) {
						if(Double.isNaN(double_arr[dl])) { System.out.println("WARNING: dimension variable contains missing values"); continue;}
						Double e = new Double(double_arr[dl]);
						if ( to_fill.containsKey(e) ) { to_fill.get(e).add(new Integer(dl)); }
						else { List<Integer> l = new ArrayList<>(); l.add(new Integer(dl)); to_fill.put(e, l); }
					} break;
				default:
				case STRING:
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
		for(int w=0; w<dim_lengths[0]; w++) { List<Integer> lw = values_dim_one.get(new Double(dimension_one[w]));
			for(int v=0; v<dim_lengths[1]; v++) { List<Integer> lv = values_dim_two.get(new Double(dimension_two[v]));
				for(int u=0; u<dim_lengths[2]; u++) { List<Integer> lu = values_dim_three.get(new Double(dimension_thr[u]));
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
	
	public void setDefaultDatatype(String _type) { default_data_type = DataType.getDataType(_type, default_data_type); }
	public String getDefaultDatatype() { return default_data_type.toString(); }


	/**
	 * 
	 * @param file_path
	 * @param delimiter
	 * @param header_length
	 * @param data_types
	 */
	public void read_csv(String file_path, String[] data_types) {
		read_csv(file_path, ",", 0, data_types);
	}
	public void read_csv(String file_path, String delimiter, String[] data_types) {
		read_csv(file_path, delimiter, 0, data_types);
	}
	public void read_csv(String file_path, int header_length, String[] data_types) {
		read_csv(file_path, ",", header_length, data_types);
	}
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
	 * if on variable isn't in the file, it woould not be read to this dataframe
	 * and a warning would be printed to the command line
	 * @param netcdf_file
	 * @param variables
	 */
	public void readFromNetcdf(NetcdfFile netcdf_file, boolean include_dimensions, String... variable) {
		clear();
		if(netcdf_file==null) { System.err.println("The Netcdf file does not exist!"); return; }
		if(variable==null || variable.length<=0) {
			System.err.println("You have to specify at least one variable you want to read from the Netcdf file!"); return; }
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
				Variable var = netcdf_file.findVariable(dim.getFullNameEscaped());
				if(var==null) {
					int[] arr = new int[datalength];
					for(int dl=0; dl<datalength; dl++) arr[dl] = indices[dl][dimid];
					addColumn(dim.getFullName(), arr);
				} else {
					Array a = null;
					try {
						a = var.read();
					} catch (IOException e) {
						System.out.println("WARNING: could not read dimension \""+dim.getFullName()+"\": does not add it to dataframe");
						continue;
					}
					Index ind = a.getIndex();
					switch(var.getDataType()) {
						case BOOLEAN:
							boolean[] bl = new boolean[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); bl[dl] = a.getBoolean(ind); }
							addColumn(dim.getFullName(), bl);
							break;
						case BYTE:
							byte[] etyb = new byte[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); etyb[dl] = a.getByte(ind); }
							addColumn(dim.getFullName(), etyb);
							break;
						case INT:
							int[] tni = new int[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); tni[dl] = a.getInt(ind); }
							addColumn(dim.getFullName(), tni);
							break;
						case SHORT:
							short[] trohs = new short[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); trohs[dl] = a.getShort(ind); }
							addColumn(dim.getFullName(), trohs);
							break;
						case LONG:
							long[] gnol = new long[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnol[dl] = a.getLong(ind); }
							addColumn(dim.getFullName(), gnol);
							break;
						case FLOAT:
							boolean hasFVf = (var.findAttribute("_FillValue")!=null);
							double ffill = Double.NaN; if(hasFVf) ffill = (float) var.findAttribute("_FillValue").getNumericValue();
							if(hasFVf) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+ffill);
							float[] taolf = new float[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); taolf[dl] = a.getFloat(ind);
								if(hasFVf && taolf[dl]==ffill) taolf[dl]=Float.NaN; }
							addColumn(dim.getFullName(), taolf);
							break;
						case DOUBLE:
							boolean hasFVd = (var.findAttribute("_FillValue")!=null);
							double dfill = Double.NaN; if(hasFVd) dfill = (double) var.findAttribute("_FillValue").getNumericValue();
							if(hasFVd) System.out.println("  [DEBUG] \""+var.getFullName()+"\":_FillValue = "+dfill);
							double[] elbuod = new double[datalength];
							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); elbuod[dl] = a.getDouble(ind);
								if(hasFVd && elbuod[dl]==dfill) elbuod[dl] = Double.NaN; }
							addColumn(dim.getFullName(), elbuod);
							break;
//						case CHAR:
//						case STRING:
//							String[] gnirts = new String[datalength];
//							for(int dl=0; dl<datalength; dl++) { ind.set(indices[dl][dimid]); gnirts[dl] = ""+a.getChar(ind); }
//							addColumn(dim.getFullName(), gnirts);
//							break;
						default:
							System.out.println("WARNING: could not add variable \""+dim.getFullName()+"\": unsupported data type!");
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
				if(dd==dlen) System.err.println("WARNING: could not find dimension \""+vdims.get(d).getFullName()+"\": maybe read data in wrong order!");
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
				default:
					System.out.println("WARNING: could not add variable \""+variable[vi]+"\": unsupported data type!");
					break;
			}
		}
	}
	/**
	 * Write all content of this dataframe to a netcdf file
	 * @param netcdf_file_path path to the netcdf file
	 */
	public void writeToNetcdf(String netcdf_file_path) {
		try {
			NetcdfFileWriter ncdfWriter = NetcdfFileWriter.createNew(Version.netcdf4, netcdf_file_path);
			String dim_name = ""; boolean is_used = true; int dim_test_num = -1;
			while(is_used) {
				dim_test_num++; dim_name = "dim"+dim_test_num;
				is_used = false;
				for(String tit: titles) if(tit.equals(dim_name)) { is_used = true; break; }
			}
			Dimension dim = ncdfWriter.addDimension(null, dim_name, datalength);
			List<Dimension> dims = new ArrayList<>(); dims.add(dim);
			Variable[] vars = new Variable[titles.length+1];
			vars[0] = ncdfWriter.addVariable(null, dim_name, ucar.ma2.DataType.INT, dims);
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
				vars[iv+1] = ncdfWriter.addVariable(null, FormatHelper.underscore_spaces(titles[iv]), dataType, dims);
				switch(types[iv]) {
					case SHORT:  ncdfWriter.addVariableAttribute(vars[iv+1], new Attribute("_FillValue", Short.MIN_VALUE));   break;
					case INT:    ncdfWriter.addVariableAttribute(vars[iv+1], new Attribute("_FillValue", Integer.MIN_VALUE)); break;
					case LONG:   ncdfWriter.addVariableAttribute(vars[iv+1], new Attribute("_FillValue", Long.MIN_VALUE));    break;
					case FLOAT:  ncdfWriter.addVariableAttribute(vars[iv+1], new Attribute("_FillValue", Constants.FILL_VALUE_F));  break;
					case DOUBLE: ncdfWriter.addVariableAttribute(vars[iv+1], new Attribute("_FillValue", Constants.FILL_VALUE_D)); break;
					default: break;
				}
			}
			ncdfWriter.addGroupAttribute(null, new Attribute("history", "created with \"EnMe-Kriging\" from Dataframe - JAVA Netcdf "+Constants.NETCDF_VERSION));
			ncdfWriter.create();
			Variable dimVar = vars[0];
			ArrayInt.D1 dim_arr = new ArrayInt.D1(datalength, false);
			for(int dl=0; dl<datalength; dl++) dim_arr.set(dl,1+dl);
			ncdfWriter.write(dimVar, dim_arr);
			for(int iv=0; iv<titles.length; iv++) {
				switch(types[iv]) {
					case BOOL:
						ArrayBoolean.D1 bool_arr = new ArrayBoolean.D1(datalength);
						boolean[] bool_source = bool_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) bool_arr.set(dl, bool_source[dl]);
						ncdfWriter.write(vars[iv+1], bool_arr);
						break;
					case BYTE:
						ArrayByte.D1 byte_arr = new ArrayByte.D1(datalength, false);
						byte[] byte_source = byte_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) byte_arr.set(dl, byte_source[dl]);
						ncdfWriter.write(vars[iv+1], byte_arr);
						break;
					case SHORT:
						ArrayShort.D1 short_arr = new ArrayShort.D1(datalength, false);
						short[] short_source = short_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) short_arr.set(dl, short_source[dl]);
						ncdfWriter.write(vars[iv+1], short_arr);
						break;
					case INT:
						ArrayInt.D1 int_arr = new ArrayInt.D1(datalength, false);
						int[] int_source = int_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) int_arr.set(dl, int_source[dl]);
						ncdfWriter.write(vars[iv+1], int_arr);
						break;
					case LONG:
						ArrayLong.D1 long_arr = new ArrayLong.D1(datalength, false);
						long[] long_source = long_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) long_arr.set(dl, long_source[dl]);
						ncdfWriter.write(vars[iv+1], long_arr);
						break;
					case FLOAT:
						ArrayFloat.D1 float_arr = new ArrayFloat.D1(datalength);
						float[] float_source = float_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) float_arr.set(dl, float_source[dl]);
						ncdfWriter.write(vars[iv+1], float_arr);
						break;
					case DOUBLE:
						ArrayDouble.D1 double_arr = new ArrayDouble.D1(datalength);
						double[] double_source = double_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) double_arr.set(dl, double_source[dl]);
						ncdfWriter.write(vars[iv+1], double_arr);
						break;
					case STRING:
						ArrayString.D1 string_arr = new ArrayString.D1(datalength);
						String[] string_source = string_column.get(titles[iv]);
						for(int dl=0; dl<datalength; dl++) string_arr.set(dl, string_source[dl]);
						ncdfWriter.write(vars[iv+1], string_arr);
						break;
					default:
						break;
				}
			}
			ncdfWriter.flush();
			ncdfWriter.close();
		} catch (IOException io_e) {
			io_e.printStackTrace();
		} catch (InvalidRangeException ir_e) {
			ir_e.printStackTrace();
		}
	}



	public int getVariableCount() { return titles.length; }
	public String getVarname(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) return null;
		return titles[vi];
	}
	public int getVariableID(String _var_name) {
		return strings_index(titles, _var_name)+Constants.FIRST_IDX;
	}
	public boolean hasVariable(String _var_name) {
		return strings_index(titles, _var_name)>=0;
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
		if(strings_index(titles, _col_name)<0)
			return null;
		switch(types[strings_index(titles,_col_name)]) {
			case BOOL:   return bool_column.get(_col_name);
			case BYTE:   return byte_column.get(_col_name);
			case INT:    return int_column.get(_col_name);
			case SHORT:  return short_column.get(_col_name);
			case LONG:   return long_column.get(_col_name);
			case FLOAT:  return float_column.get(_col_name);
			case DOUBLE: return double_column.get(_col_name);
			case STRING: return string_column.get(_col_name);
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
			}
		}
		boolean[] coldivs = new boolean[clen-1];
		for(int c=0; c<clen-1; c++) coldivs[c] = true;
		boolean[] rowdivs = new boolean[dlen-1];
		for(int r=0; r<dlen-1; r++) rowdivs[r] = false;
		rowdivs[0] = true; rowdivs[1] = true;
		FormatHelper.printTable(out, coldivs, rowdivs, true);
	}
	public void printfull() { head(0); }
	public void describe() {
		int clen = titles.length;
		String[][] out = new String[9][clen+1];
		out[0][0] = " "; out[1][0] = " Count "; out[2][0] = " mean ";
		out[3][0] = " std "; out[4][0] = " min "; out[5][0] = " 25% ";
		out[6][0] = " 50% "; out[7][0] = " 75% "; out[8][0] = " max ";
		for(int c=0; c<clen; c++) {
			int cc = c+1;
			out[0][cc] = titles[c]+" ";
			switch(types[c]) {
				case BOOL:
					out[2][cc] = " "; out[3][cc] = " "; //mean and std
					out[4][cc] = " "; out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " ";
					int olen = 0;
					for(boolean o: bool_column.get(titles[c])) if(o) { olen++; }
					out[1][cc] = " "+olen+" ";
					break;
				case BYTE: double bmean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " --- "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " NaN ";
					List<Byte> bytes = new ArrayList<Byte>();
					for(byte b: byte_column.get(titles[c])) if(Byte.MIN_VALUE!=b) { bmean += b; bytes.add(b); }
					int blen = bytes.size();
					out[1][cc] = " "+blen+" ";
					if(blen==0) break;
					bmean /= blen;
					byte bmean_b = (byte) ( (int)(bmean+0.5d) - (bmean<-0.5d ? 1 : 0) );
					out[2][cc] = " "+bmean_b+" ";
					double bstd = 0d;
					for(byte b: bytes) bstd += (b-bmean)*(b-bmean);
					if(blen==1) { bstd = Byte.MIN_VALUE; } else { bstd /= (blen-1); }
					byte bstd_b = (byte) (bstd+0.99d);
					out[3][cc] = " "+bstd_b+" ";
					Byte[] a_b = bytes.toArray(new Byte[0]);
					Arrays.sort(a_b);
					out[4][cc] = " "+a_b[0]+" "; out[8][cc] = " "+a_b[blen-1]+" ";
					if(blen>3) out[5][cc] = " "+a_b[blen/4]+" "; if(blen>1) out[6][cc] = " "+a_b[(blen+1)/2]+" "; if(blen>3) out[7][cc] = " "+a_b[(3*blen)/4]+" ";
					break;
				case INT: double imean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " --- "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " NaN ";
					List<Integer> integers = new ArrayList<Integer>();
					for(int i: int_column.get(titles[c])) if(Integer.MIN_VALUE!=i) { imean += i; integers.add(i); }
					int ilen = integers.size();
					out[1][cc] = " "+ilen+" ";
					if(ilen==0) break;
					imean /= ilen;
					int imean_i = (int) (imean+0.5d) - (imean<-0.5d ? 1 : 0);
					out[2][cc] = " "+imean_i+" ";
					double istd = 0d;
					for(int i: integers) istd += (i-imean)*(i-imean);
					if(ilen==1) { istd = Integer.MIN_VALUE;; } else { istd /= (ilen-1); }
					int istd_i = (int) (istd+0.99d);
					out[3][cc] = " "+istd_i+" ";
					Integer[] a_i = integers.toArray(new Integer[0]);
					Arrays.sort(a_i);
					out[4][cc] = " "+a_i[0]+" "; out[8][cc] = " "+a_i[ilen-1]+" ";
					if(ilen>3) out[5][cc] = " "+a_i[ilen/4]+" "; if(ilen>1) out[6][cc] = " "+a_i[(ilen+1)/2]+" "; if(ilen>3) out[7][cc] = " "+a_i[(3*ilen)/4]+" ";
					break;
				case SHORT: double rmean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " --- "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " NaN ";
					List<Short> shorts = new ArrayList<Short>();
					for(short r: short_column.get(titles[c])) if(Short.MIN_VALUE!=r) { rmean += r; shorts.add(r); }
					int rlen = shorts.size();
					out[1][cc] = " "+rlen+" ";
					if(rlen==0) break;
					rmean /= rlen;
					short rmean_s = (short) ( (int)(rmean+0.5d) - (rmean<-0.5d ? 1 : 0) );
					out[2][cc] = " "+rmean_s+" ";
					double rstd = 0d;
					for(short r: shorts) rstd += (r-rmean)*(r-rmean);
					if(rlen==1) { rstd = Integer.MIN_VALUE;; } else { rstd /= (rlen-1); }
					short rstd_s = (short) (rstd+0.99d);
					out[3][cc] = " "+rstd_s+" ";
					Short[] a_r = shorts.toArray(new Short[0]);
					Arrays.sort(a_r);
					out[4][cc] = " "+a_r[0]+" "; out[8][cc] = " "+a_r[rlen-1]+" ";
					if(rlen>3) out[5][cc] = " "+a_r[rlen/4]+" "; if(rlen>1) out[6][cc] = " "+a_r[(rlen+1)/2]+" "; if(rlen>3) out[7][cc] = " "+a_r[(3*rlen)/4]+" ";
					break;
				case LONG: double lmean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " --- "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " NaN ";
					List<Long> longs = new ArrayList<Long>();
					for(long l: long_column.get(titles[c])) if(Long.MIN_VALUE!=l) { lmean += l; longs.add(l); }
					int llen = longs.size();
					out[1][cc] = " "+llen+" ";
					if(llen==0) break;
					lmean /= llen;
					long imean_l = (long) (lmean+0.5d) - (lmean<-0.5d ? 1L : 0L);
					out[2][cc] = " "+imean_l+" ";
					double lstd = 0d;
					for(long l: longs) lstd += (l-lmean)*(l-lmean);
					if(llen==1) { lstd = Long.MIN_VALUE; } else { lstd /= (llen-1); }
					long istd_l = (long) (lstd+0.99d);
					out[3][cc] = " "+istd_l+" ";
					Long[] a_l = longs.toArray(new Long[0]);
					Arrays.sort(a_l);
					out[4][cc] = " "+a_l[0]+" "; out[8][cc] = " "+a_l[llen-1]+" ";
					if(llen>3) out[5][cc] = " "+a_l[llen/4]+" "; if(llen>1) out[6][cc] = " "+a_l[(llen+1)/2]+" "; if(llen>3) out[7][cc] = " "+a_l[(3*llen)/4]+" ";
					break;
				case FLOAT: float fmean = 0f;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " NaN "; out[6][cc] = " NaN "; out[7][cc] = " NaN "; out[8][cc] = " NaN ";
					List<Float> floats = new ArrayList<Float>();
					for(float f: float_column.get(titles[c])) if(!Float.isNaN(f)) { fmean += f; floats.add(f); }
					int flen = floats.size();
					out[1][cc] = " "+flen+" ";
					if(flen==0) break;
					fmean /= flen;
					out[2][cc] = " "+fmean+" ";
					float fstd = 0f;
					for(float f: floats) fstd += (f-fmean)*(f-fmean);
					if(flen==1) { fstd = Float.POSITIVE_INFINITY; } else { fstd /= (flen-1f); }
					out[3][cc] = " "+fstd+" ";
					Float[] a_f = floats.toArray(new Float[0]);
					Arrays.sort(a_f);
					out[4][cc] = " "+a_f[0]+" "; out[8][cc] = " "+a_f[flen-1]+" ";
					if(flen>3) out[5][cc] = " "+a_f[flen/4]+" "; if(flen>1) out[6][cc] = " "+a_f[(flen+1)/2]+" "; if(flen>3) out[7][cc] = " "+a_f[(3*flen)/4]+" ";
					break;
				case DOUBLE: double dmean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " NaN "; out[6][cc] = " NaN "; out[7][cc] = " NaN "; out[8][cc] = " NaN ";
					List<Double> doubles = new ArrayList<Double>();
					for(double d: double_column.get(titles[c])) if(!Double.isNaN(d)) { dmean += d; doubles.add(d); }
					int dlen = doubles.size();
					out[1][cc] = " "+dlen+" ";
					if(dlen==0) break;
					dmean /= dlen;
					out[2][cc] = " "+dmean+" ";
					double dstd = 0d;
					for(double d: doubles) dstd += (d-dmean)*(d-dmean);
					if(dlen==1) { dstd = Double.POSITIVE_INFINITY; } else { dstd /= (dlen-1d); }
					out[3][cc] = " "+dstd+" ";
					Double[] a_d = doubles.toArray(new Double[0]);
					Arrays.sort(a_d);
					out[4][cc] = " "+a_d[0]+" "; out[8][cc] = " "+a_d[dlen-1]+" ";
					if(dlen>3) out[5][cc] = " "+a_d[dlen/4]+" "; if(dlen>1) out[6][cc] = " "+a_d[(dlen+1)/2]+" "; if(dlen>3) out[7][cc] = " "+a_d[(3*dlen)/4]+" ";
					break;
				case STRING:
					out[2][cc] = " "; out[3][cc] = " "; //mean and std
					out[4][cc] = " "; out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " ";
					int slen = 0;
					for(String s: string_column.get(titles[c])) if(!(s.equals(" ") || s.length()<1)) { slen++; }
					out[1][cc] = " "+slen+" ";
					break;
				default:
					out[2][cc] = " "; out[3][cc] = " "; //mean and std
					out[4][cc] = " "; out[5][cc] = " "; out[6][cc] = " "; out[7][cc] = " "; out[8][cc] = " ";
			}
		}
		boolean[] coldivs = {true,false,false,false,false,false,false,false};
		boolean[] rowdivs = new boolean[clen]; rowdivs[0] = true;
		for(int r=1; r<clen; r++) rowdivs[r] = false;
		FormatHelper.printTable(out, coldivs, rowdivs, true);
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
	}

	private int strings_index(String[] sarr, String search_regex) {
		if(sarr==null || search_regex==null) return -1;
		if(sarr.length<1) return -1;
		for(int s=0; s<sarr.length; s++) { if(sarr[s]==null) continue; if(sarr[s].equals(search_regex)) return s; }
		return -1;
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
	private void add_variable(String _name, DataType _type) {
		String[] old_tit = new String[titles.length];
		DataType[] old_tp = new DataType[types.length];
		for(int t=0; t<titles.length; t++) { old_tit[t] = titles[t]; old_tp[t] = types[t]; }
		titles = null; titles = new String[old_tit.length+1];
		types = null; types = new DataType[old_tp.length+1];
		for(int t=0; t<old_tit.length; t++) { titles[t] = old_tit[t]; types[t] = old_tp[t]; }
		titles[old_tit.length] = _name;
		types[old_tp.length] = _type;
	}
	private void add_stdana(double _min, double _max, double _mean, double _sill) {
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
	
	public enum DataType {
		BOOL("bool"), BYTE("byte"), INT("int"), SHORT("short"), LONG("long"),
		FLOAT("float"), DOUBLE("double"),
		STRING("string");
		
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
			return _default;
		}
	}
}
