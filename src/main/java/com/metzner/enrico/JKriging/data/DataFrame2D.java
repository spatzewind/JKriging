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

import com.metzner.enrico.JKriging.data.DataFrame.DataType;
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

public class DataFrame2D {

	private int[] datalength;
	private DataType default_data_type;
	private String[] titles;
	private DataType[] types;
	private String[] dimension_names;
	private double[] dimension_one, dimension_two;
	private double[][] minmax_mean_sill;
	private Map<String, boolean[][]> bool_column;
	private Map<String, byte[][]>    byte_column;
	private Map<String, int[][]>     int_column;
	private Map<String, short[][]>   short_column;
	private Map<String, long[][]>    long_column;
	private Map<String, float[][]>   float_column;
	private Map<String, double[][]>  double_column;
	private Map<String, String[][]>  string_column;




	public DataFrame2D() {
		datalength = new int[] { 0, 0 };
		default_data_type = DataType.STRING;
		titles = new String[0];
		types = new DataType[0];
		int[] currentDate = DataHelper.JulianDate.jd2cal(DataHelper.JulianDate.now());
		String currentDateString = currentDate[0]+FormatHelper.nf(currentDate[1], 2,'0')+FormatHelper.nf(currentDate[2],2,'0');
		currentDateString += FormatHelper.nf(currentDate[3],2,'0')+FormatHelper.nf(currentDate[4],2,'0')+FormatHelper.nf(currentDate[5],2,'0');
		dimension_names = new String[] {"Dim"+currentDateString+"A0","Dim"+currentDateString+"B1"};
		dimension_one = new double[0];
		dimension_two = new double[0];
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




	public void addColumn(String _column_name, boolean[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, boolean[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.BOOL);
		add_stdana(Double.NaN,Double.NaN,Double.NaN,Double.NaN);
		boolean[][] new_data = new boolean[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++) for(int t1=0; t1<datalength[1]; t1++)
			if(t0<_column_data.length && t1<_column_data[0].length) {
				new_data[t0][t1] = _column_data[t0][t1]; } else {
				new_data[t0][t1] = false; }
		bool_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, byte[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, byte[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.BYTE);
		byte[] inax = StdAnalysis.minmax(_column_data); byte[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		byte[][] new_data = new byte[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Byte.MIN_VALUE; }
		byte_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, short[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, short[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.SHORT);
		short[] inax = StdAnalysis.minmax(_column_data); short[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		short[][] new_data = new short[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Short.MIN_VALUE; }
		short_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, int[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, int[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.INT);
		int[] inax = StdAnalysis.minmax(_column_data); int[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		int[][] new_data = new int[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Integer.MIN_VALUE; }
		int_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, long[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, long[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.LONG);
		long[] inax = StdAnalysis.minmax(_column_data); long[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		long[][] new_data = new long[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Long.MIN_VALUE; }
		long_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, float[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, float[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.FLOAT);
		float[] inax = StdAnalysis.minmax(_column_data); float[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		float[][] new_data = new float[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Float.NaN; }
		float_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, double[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, double[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.DOUBLE);
		double[] inax = StdAnalysis.minmax(_column_data); double[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		double[][] new_data = new double[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = Double.NaN; }
		double_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, String[][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, String[][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1]) && !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]); }
		add_variable(_column_name_c, DataType.STRING); add_stdana(Double.NaN,-Double.NaN,Double.NaN,Double.NaN);
		String[][] new_data = new String[datalength[0]][datalength[1]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				if(t0<_column_data.length && t1<_column_data[0].length) {
					new_data[t0][t1] = _column_data[t0][t1]; } else {
					new_data[t0][t1] = ""; }
		string_column.put(_column_name_c, new_data);
	}
	public void setDimension(int dim_num, String dim_name) {
		setDimension(dim_num, null, dim_name);
	}
	public void setDimension(int dim_num, double[] dim_values) {
		setDimension(dim_num, dim_values, null);
	}
	public void setDimension(int dim_num, double[] dim_values, String dim_name) {
		String _dim_name_c = FormatHelper.name2CFConvention(dim_name);
		int di = dim_num - Constants.FIRST_IDX;
		if(di<0 || di>1) {
			System.err.println("Number of dimension must be "+Constants.FIRST_IDX+" or "+(1+Constants.FIRST_IDX)+" for DataFrame2D!");
			DataHelper.printStackTrace(System.err);
			return;
		}
		boolean change_name = (_dim_name_c!=null);
		boolean change_values = (dim_values!=null);
		if(change_name)
			if(_dim_name_c.equals(dimension_names[(di+1)%2]) || hasVariable(_dim_name_c)) {
				System.err.println("Cannot rename dimension "+(di+Constants.FIRST_IDX)+", an other dimension or variable has the same name!");
				DataHelper.printStackTrace(System.err);
				return;
			}
		if(change_values) {
			if(dim_values.length==0) {
				System.err.println("Length of dimension "+(di+Constants.FIRST_IDX)+" has to be greater than 0! Not a valid DataFrame2D.");
				DataHelper.printStackTrace(System.err);
				return;
			}
			if(dim_values.length!=datalength[di] && datalength[di]>0) {
				System.err.println("Number of values must match corresponding dimension length:\n    expect "+datalength[di]+" but got "+dim_values.length);
				DataHelper.printStackTrace(System.err);
				return;
			}
		}
		if(change_name)
			dimension_names[di] = _dim_name_c;
		if(change_values) {
			if(datalength[di]==0) {
				datalength[di] = dim_values.length;
				if(dim_num==1) dimension_one = new double[dim_values.length];
				if(dim_num==2) dimension_two = new double[dim_values.length];
			}
			if(dim_num==1)
				for(int dn=0; dn<datalength[0]; dn++)
					dimension_one[dn] =
						dim_values[dn];
			if(dim_num==2)
				for(int dn=0; dn<datalength[1]; dn++)
					dimension_two[dn] =
						dim_values[dn];
		}
	}
	public void copyDimensionsFrom(DataFrame2D otherDF) {
		this.setDimension(0+Constants.FIRST_IDX, otherDF.getDimensionValues(0+Constants.FIRST_IDX), otherDF.getDimensionName(0+Constants.FIRST_IDX));
		this.setDimension(1+Constants.FIRST_IDX, otherDF.getDimensionValues(1+Constants.FIRST_IDX), otherDF.getDimensionName(1+Constants.FIRST_IDX));
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
		if(titles.length==1) { titles = new String[0]; types = new DataType[0]; datalength = new int[] {0, 0}; }
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
			case BOOL: boolean[][] bool_arr = bool_column.get(variable_name);
				switch(_new_type) {
					case BYTE: byte[][] barr = new byte[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte)(bool_arr[j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = (short)(bool_arr[j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = (bool_arr[j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = (bool_arr[j][i]?1L:0L);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = (bool_arr[j][i]?1f:0f);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = (bool_arr[j][i]?1d:0d);
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[bool_arr.length][bool_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = (bool_arr[j][i]?"true":"false");
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case BYTE: byte[][] byte_arr = byte_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][] oarr = new boolean[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = (byte_arr[j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case SHORT: short[][] rarr = new short[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = byte_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = byte_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = byte_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = byte_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = byte_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[byte_arr.length][byte_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+byte_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case SHORT: short[][] short_arr = short_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][] oarr = new boolean[short_arr.length][short_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = (short_arr[j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[short_arr.length][short_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte) short_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case INT: int[][] iarr = new int[short_arr.length][short_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = short_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[short_arr.length][short_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = short_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[short_arr.length][short_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = short_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[short_arr.length][short_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = short_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[short_arr.length][short_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+short_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case INT: int[][] int_arr = int_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][] oarr = new boolean[int_arr.length][int_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = (int_arr[j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[int_arr.length][int_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte) int_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[int_arr.length][int_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = (short) int_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case LONG: long[][] larr = new long[int_arr.length][int_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = int_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[int_arr.length][int_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = int_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[int_arr.length][int_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = int_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[int_arr.length][int_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+int_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case LONG: long[][] long_arr = long_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][] oarr = new boolean[long_arr.length][long_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = (long_arr[j][i]!=0L ? true : false);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[long_arr.length][long_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte) long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[long_arr.length][long_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = (short) long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[long_arr.length][long_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = (int) long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[long_arr.length][long_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[long_arr.length][long_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[long_arr.length][long_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = long_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[long_arr.length][long_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+long_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case FLOAT: float[][] float_arr = float_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert float to bool result to true except NaN!");
						boolean[][] oarr = new boolean[float_arr.length][float_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = !Float.isNaN(float_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[float_arr.length][float_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte) float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[float_arr.length][float_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = (short) float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[float_arr.length][float_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = (int) float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[float_arr.length][float_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = (long) float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[float_arr.length][float_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[float_arr.length][float_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = float_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[float_arr.length][float_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+float_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case DOUBLE: double[][] double_arr = double_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert double to bool result to true except NaN!");
						boolean[][] oarr = new boolean[double_arr.length][double_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = !Double.isNaN(double_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[double_arr.length][double_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = (byte) double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[double_arr.length][double_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = (short) double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[double_arr.length][double_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = (int) double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[double_arr.length][double_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = (long) double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[double_arr.length][double_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = (float) double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[double_arr.length][double_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = double_arr[j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][] sarr = new String[double_arr.length][double_arr[0].length];
						for(int j=0; j<sarr.length; j++) for(int i=0; i<sarr[0].length; i++) sarr[j][i] = ""+double_arr[j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case STRING: String[][] string_arr = string_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][] oarr = new boolean[string_arr.length][string_arr[0].length];
						for(int j=0; j<oarr.length; j++) for(int i=0; i<oarr[0].length; i++) oarr[j][i] = set_boolean(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][] barr = new byte[string_arr.length][string_arr[0].length];
						for(int j=0; j<barr.length; j++) for(int i=0; i<barr[0].length; i++) barr[j][i] = set_byte(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][] rarr = new short[string_arr.length][string_arr[0].length];
						for(int j=0; j<rarr.length; j++) for(int i=0; i<rarr[0].length; i++) rarr[j][i] = set_short(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][] iarr = new int[string_arr.length][string_arr[0].length];
						for(int j=0; j<iarr.length; j++) for(int i=0; i<iarr[0].length; i++) iarr[j][i] = set_int(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][] larr = new long[string_arr.length][string_arr[0].length];
						for(int j=0; j<larr.length; j++) for(int i=0; i<larr[0].length; i++) larr[j][i] = set_long(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][] farr = new float[string_arr.length][string_arr[0].length];
						for(int j=0; j<farr.length; j++) for(int i=0; i<farr[0].length; i++) farr[j][i] = set_float(string_arr[j][i]);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][] darr = new double[string_arr.length][string_arr[0].length];
						for(int j=0; j<darr.length; j++) for(int i=0; i<darr[0].length; i++) darr[j][i] = set_double(string_arr[j][i]);
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
	
/*	public DataFrame2D concat(DataFrame2D df2, int dimension_to_append) {
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
			switch(types[iv]) {
				case BOOL:   bool_column.put(  titles[iv], DataHelper.concat_bool_array(  (boolean[])this.getArray(iv), (boolean[])df2.getArray(titles[iv]))); break;
				case BYTE:   byte_column.put(  titles[iv], DataHelper.concat_byte_array(  (byte[])this.getArray(iv),    (byte[])df2.getArray(titles[iv])));    break;
				case SHORT:  short_column.put( titles[iv], DataHelper.concat_short_array( (short[])this.getArray(iv),   (short[])df2.getArray(titles[iv])));   break;
				case INT:    int_column.put(   titles[iv], DataHelper.concat_int_array(   (int[])this.getArray(iv),     (int[])df2.getArray(titles[iv])));     break;
				case LONG:   long_column.put(  titles[iv], DataHelper.concat_long_array(  (long[])this.getArray(iv),    (long[])df2.getArray(titles[iv])));    break;
				case FLOAT:  float_column.put( titles[iv], DataHelper.concat_float_array( (float[])this.getArray(iv),   (float[])df2.getArray(titles[iv])));   break;
				case DOUBLE: double_column.put(titles[iv], DataHelper.concat_double_array((double[])this.getArray(iv),  (double[])df2.getArray(titles[iv])));  break;
				case STRING: string_column.put(titles[iv], DataHelper.concat_string_array((String[])this.getArray(iv),  (String[])df2.getArray(titles[iv])));  break;
				default: break;
			}
		}
		datalength += df2.datalength;
		return this;
	} */

	public void setDefaultDatatype(String _type) { default_data_type = DataType.getDataType(_type, default_data_type); }
	public String getDefaultDatatype() { return default_data_type.toString(); }


/*	/**
	 * 
	 * @param file_path
	 * @param delimiter
	 * @param header_length
	 * @param data_types
	 * /
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
	} */

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
			int[] datarrlen = { Integer.parseInt(parts[1].trim()), Integer.parseInt(parts[2].trim()) };
			lines.add(""+line);
			for(int i=0; i<nvari; i++) line = br.readLine();
			maxdat = -1;
			while(line!=null) {
				maxdat++;
				line = br.readLine();
			}
			System.out.println("maxdat = "+maxdat+"    with shape "+datarrlen[0]+" x "+datarrlen[1]);
			if(datarrlen[0]*datarrlen[1]!=maxdat) {
				System.err.println("Shape is not compatible with the number of data entries! ("+datarrlen[0]+"*"+datarrlen[1]+"!="+maxdat+")");
				return;
			}
			
			names = new String[nvari];
			double[][][] vr = new double[nvari][datarrlen[0]][datarrlen[1]];
			
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
				nd++; int u=nd%datarrlen[1], v=nd/datarrlen[1];
				for(int j=0; j<nvari; j++) {
					//vr[j][nd] = Double.parseDouble(parts[j].trim());
					vr[j][v][u] = set_double(parts[j].trim()); //with error check for reading errors and NaNs
				}
			}
			nd++;
			System.out.println("nd     = "+nd);
			//System.out.println("maxdat = "+maxdat);
			
			this.clear();
			//titles = new String[names.length];
			//types  = new DataType[names.length];
			datalength = new int[] {datarrlen[0],datarrlen[1]};
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
	public void readFromNetcdf(NetcdfFile netcdf_file, String... variable) {
		clear();
		if(netcdf_file==null) { System.err.println("The Netcdf file does not exist!"); return; }
		if(variable==null || variable.length<=0) {
			System.err.println("You have to specify at least one variable you want to read from the Netcdf file!"); return; }
		boolean[] var_exist = new boolean[variable.length];
		List<Dimension> dims = new ArrayList<Dimension>();
		boolean at_least_one_exist = false;
		for(int b=0; b<var_exist.length; b++) {
			Variable var = netcdf_file.findVariable(variable[b]);
			var_exist[b] = (var!=null);
			if(var_exist[b]) {
				if(var.getDimensions().size()>2) {
					var_exist[b] = false;
					continue;
				}
				for(Dimension d: var.getDimensions()) {
					if(!dims.contains(d))
						dims.add(d);
				}
				at_least_one_exist = true;
			}
		}
		if(!at_least_one_exist) {
			System.err.println("At least one variable has to have two dimensions or two variables have to have two dimensions together!");
			DataHelper.printStackTrace(System.err); return;
		}
		if(dims.size()!=2) {
			System.err.println("The 2D-dataframe can only handle exact 2 dimensions, but found "+dims.size()+"!");
			DataHelper.printStackTrace(System.err); return;
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
		dimension_one = new double[dimlen[0]];
		dimension_two = new double[dimlen[1]];
		for(int dimid=0; dimid<2; dimid++) {
			Variable var = netcdf_file.findVariable(dims.get(dimid).getFullNameEscaped());
			boolean was_succesful = false;
			if(var!=null) {
				Array a = null;
				try {
					a = var.read();
					if(dimid==0)
						for(int dl=0; dl<dimlen[0]; dl++)
							dimension_one[dl] = a.getDouble(dl);
					if(dimid==1)
						for(int dl=0; dl<dimlen[1]; dl++)
							dimension_two[dl] = a.getDouble(dl);
					dimension_names[dimid] = var.getFullName();
					was_succesful = true;
				} catch (IOException e) {
					System.out.println("WARNING: could not read dimension \""+dims.get(dimid).getFullName()+"\": add as index-array to dataframe");
				} catch (NumberFormatException nfe) {
					System.out.println("WARNING: could not read dimension \""+dims.get(dimid).getFullName()+"\": add as index-array to dataframe");
				}
			}
			if(!was_succesful) {
				if(dimid==0) for(int dl=0; dl<dimension_one.length; dl++) dimension_one[dl] = dl+1d;
				if(dimid==1) for(int dl=0; dl<dimension_two.length; dl++) dimension_two[dl] = dl+1d;
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
			Dimension[] vardims = var.getDimensions().toArray(new Dimension[0]);
			if(vardims.length==0) {
				ind.set(new int[0]);
				switch(var.getDataType()) {
					case BOOLEAN: boolean bool = a.getBoolean(ind); boolean[][] bool_arr = new boolean[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) bool_arr[v][u] = bool; 
						addColumn(var.getFullName(), bool_arr); break;
					case BYTE: byte etyb = a.getByte(ind); byte[][] byte_arr = new byte[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) byte_arr[v][u] = etyb; 
						addColumn(var.getFullName(), byte_arr); break;
					case SHORT: short trohs = a.getShort(ind); short[][] short_arr = new short[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) short_arr[v][u] = trohs; 
						addColumn(var.getFullName(), short_arr); break;
					case INT: int tni = a.getInt(ind); int[][] int_arr = new int[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) int_arr[v][u] = tni; 
						addColumn(var.getFullName(), int_arr); break;
					case LONG: long gnol = a.getLong(ind); long[][] long_arr = new long[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) long_arr[v][u] = gnol; 
						addColumn(var.getFullName(), long_arr); break;
					case FLOAT: float taolf = a.getFloat(ind); float[][] float_arr = new float[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) float_arr[v][u] = taolf; 
						addColumn(var.getFullName(), float_arr); break;
					case DOUBLE: double elbuod = a.getDouble(ind); double[][] double_arr = new double[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) double_arr[v][u] = elbuod; 
						addColumn(var.getFullName(), double_arr); break;
//						case CHAR:
//						case STRING:
//							break;
					default:
						System.out.println("WARNING: does not this datatype: "+var.getDataType().name()+
								", so the variable is not added to the dataframe");
						break;
				}
			} else if(vardims.length==1) {
				int uf=0, vf=1; Dimension d = var.getDimension(0);
				if(d.equals(dims.get(1))) { uf=1; vf=0; }
				switch(var.getDataType()) {
					case BOOLEAN: boolean[][] bool_arr = new boolean[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); bool_arr[v][u] = a.getBoolean(ind); }
						addColumn(var.getFullName(), bool_arr); break;
					case BYTE: byte[][] byte_arr = new byte[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); byte_arr[v][u] = a.getByte(ind); }
						addColumn(var.getFullName(), byte_arr); break;
					case SHORT: short[][] short_arr = new short[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); short_arr[v][u] = a.getShort(ind); }
						addColumn(var.getFullName(), short_arr); break;
					case INT: int[][] int_arr = new int[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); int_arr[v][u] = a.getInt(ind); }
						addColumn(var.getFullName(), int_arr); break;
					case LONG: long[][] long_arr = new long[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); long_arr[v][u] = a.getLong(ind); }
						addColumn(var.getFullName(), long_arr); break;
					case FLOAT: float[][] float_arr = new float[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); float_arr[v][u] = a.getFloat(ind); }
						addColumn(var.getFullName(), float_arr); break;
					case DOUBLE: double[][] double_arr = new double[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) { ind.set(v*vf+u*uf); double_arr[v][u] = a.getDouble(ind); }
						addColumn(var.getFullName(), double_arr); break;
	//				case CHAR:
	//				case STRING:
	//					break;
					default:
						System.out.println("WARNING: does not this datatype: "+var.getDataType().name()+
								", so the variable is not added to the dataframe");
						break;
				}
			} else {
				int au=0,av=0, bu=0,bv=0; Dimension d0=var.getDimension(0), d1=var.getDimension(1);
				if(d0.equals(dims.get(0))) av=1; if(d0.equals(dims.get(1))) au=1;
				if(d1.equals(dims.get(0))) bv=1; if(d1.equals(dims.get(1))) bu=1;
				switch(var.getDataType()) {
					case BOOLEAN: boolean[][] bool_arr = new boolean[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); bool_arr[v][u] = a.getBoolean(ind); }
						addColumn(var.getFullName(), bool_arr); break;
					case BYTE: byte[][] byte_arr = new byte[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); byte_arr[v][u] = a.getByte(ind); }
						addColumn(var.getFullName(), byte_arr); break;
					case SHORT: short[][] short_arr = new short[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); short_arr[v][u] = a.getShort(ind); }
						addColumn(var.getFullName(), short_arr); break;
					case INT: int[][] int_arr = new int[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); int_arr[v][u] = a.getInt(ind); }
						addColumn(var.getFullName(), int_arr); break;
					case LONG: long[][] long_arr = new long[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); long_arr[v][u] = a.getLong(ind); }
						addColumn(var.getFullName(), long_arr); ;
					case FLOAT: float[][] float_arr = new float[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); float_arr[v][u] = a.getFloat(ind); }
						addColumn(var.getFullName(), float_arr); break;
					case DOUBLE: double[][] double_arr = new double[dimlen[0]][dimlen[1]];
						for(int v=0; v<dimlen[0]; v++) for(int u=0; u<dimlen[1]; u++) {
							ind.set(v*av+u*au,v*bv+u*bu); double_arr[v][u] = a.getDouble(ind); }
						addColumn(var.getFullName(), double_arr); break;
	//				case CHAR:
	//				case STRING:
	//					break;
					default:
						System.out.println("WARNING: does not this datatype: "+var.getDataType().name()+
								", so the variable is not added to the dataframe");
						break;
				}
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
//			String dim0_name = ""; boolean is_used = true; int dim_test_num = -1;
//			while(is_used) {
//				dim_test_num++; dim0_name = "dim"+dim_test_num;
//				is_used = false;
//				for(String tit: titles) if(tit.equals(dim0_name)) { is_used = true; break; }
//			}
//			String dim1_name = "";
//			while(is_used) {
//				dim_test_num++; dim1_name = "dim"+dim_test_num;
//				is_used = false;
//				for(String tit: titles) if(tit.equals(dim1_name)) { is_used = true; break; }
//			}
			Dimension dim0 = ncdfWriter.addDimension(null, dimension_names[0], datalength[0]);
			Dimension dim1 = ncdfWriter.addDimension(null, dimension_names[1], datalength[1]);
			List<Dimension> dims = new ArrayList<>(); dims.add(dim0);
			Variable[] vars = new Variable[titles.length+2];
			vars[0] = ncdfWriter.addVariable(null, dimension_names[0], ucar.ma2.DataType.DOUBLE, dims);
			dims.remove(0); dims.add(dim1);
			vars[1] = ncdfWriter.addVariable(null, dimension_names[1], ucar.ma2.DataType.DOUBLE, dims);
			dims.remove(0); dims.add(dim0); dims.add(dim1);
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
				vars[iv+2] = ncdfWriter.addVariable(null, FormatHelper.name2CFConvention(titles[iv]), dataType, dims);
				switch(types[iv]) {
					case SHORT:  ncdfWriter.addVariableAttribute(vars[iv+2], new Attribute("_FillValue", Short.MIN_VALUE));   break;
					case INT:    ncdfWriter.addVariableAttribute(vars[iv+2], new Attribute("_FillValue", Integer.MIN_VALUE)); break;
					case LONG:   ncdfWriter.addVariableAttribute(vars[iv+2], new Attribute("_FillValue", Long.MIN_VALUE));    break;
					case FLOAT:  ncdfWriter.addVariableAttribute(vars[iv+2], new Attribute("_FillValue", Constants.FILL_VALUE_F));  break;
					case DOUBLE: ncdfWriter.addVariableAttribute(vars[iv+2], new Attribute("_FillValue", Constants.FILL_VALUE_D)); break;
					default: break;
				}
			}
			ncdfWriter.addGroupAttribute(null, new Attribute("history", "created with \"EnMe-Kriging\" from Dataframe2D - JAVA Netcdf "+Constants.NETCDF_VERSION));
			ncdfWriter.create();
			Variable dim0Var = vars[0];
			ArrayDouble.D1 dim0_arr = new ArrayDouble.D1(datalength[0]);
			for(int dl=0; dl<datalength[0]; dl++) dim0_arr.set(dl,dimension_one[dl]);
			ncdfWriter.write(dim0Var, dim0_arr);
			Variable dim1Var = vars[1];
			ArrayDouble.D1 dim1_arr = new ArrayDouble.D1(datalength[1]);
			for(int dl=0; dl<datalength[1]; dl++) dim1_arr.set(dl,dimension_two[dl]);
			ncdfWriter.write(dim1Var, dim1_arr);
			for(int iv=0; iv<titles.length; iv++) {
				switch(types[iv]) {
					case BOOL:
						ArrayBoolean.D2 bool_arr = new ArrayBoolean.D2(datalength[0],datalength[1]);
						boolean[][] bool_source = bool_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								bool_arr.set(dl0,dl1, bool_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], bool_arr);
						break;
					case BYTE:
						ArrayByte.D2 byte_arr = new ArrayByte.D2(datalength[0],datalength[1],false);
						byte[][] byte_source = byte_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								byte_arr.set(dl0,dl1, byte_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], byte_arr);
						break;
					case SHORT:
						ArrayShort.D2 short_arr = new ArrayShort.D2(datalength[0],datalength[1],false);
						short[][] short_source = short_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								short_arr.set(dl0,dl1, short_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], short_arr);
						break;
					case INT:
						ArrayInt.D2 int_arr = new ArrayInt.D2(datalength[0],datalength[1],false);
						int[][] int_source = int_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								int_arr.set(dl0,dl1, int_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], int_arr);
						break;
					case LONG:
						ArrayLong.D2 long_arr = new ArrayLong.D2(datalength[0],datalength[1],false);
						long[][] long_source = long_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								long_arr.set(dl0,dl1, long_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], long_arr);
						break;
					case FLOAT:
						ArrayFloat.D2 float_arr = new ArrayFloat.D2(datalength[0],datalength[1]);
						float[][] float_source = float_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								float_arr.set(dl0,dl1, float_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], float_arr);
						break;
					case DOUBLE:
						ArrayDouble.D2 double_arr = new ArrayDouble.D2(datalength[0],datalength[1]);
						double[][] double_source = double_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								double_arr.set(dl0,dl1, double_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], double_arr);
						break;
					case STRING:
						ArrayString.D2 string_arr = new ArrayString.D2(datalength[0],datalength[1]);
						String[][] string_source = string_column.get(titles[iv]);
						for(int dl0=0; dl0<datalength[0]; dl0++)
							for(int dl1=0; dl1<datalength[1]; dl1++)
								string_arr.set(dl0,dl1, string_source[dl0][dl1]);
						ncdfWriter.write(vars[iv+2], string_arr);
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



	public double[] getDimensionValues(int dim_id) {
		int di = dim_id - Constants.FIRST_IDX;
		if(di<0 || di>1) return null;
		return (di==0 ? dimension_one : dimension_two);
	}
	public String getDimensionName(int dim_id) {
		int di = dim_id - Constants.FIRST_IDX;
		if(di<0 || di>1) return null;
		return dimension_names[di];
	}
	public int getVariableCount() { return titles.length; }
	public String getVarname(int _var_id) {
		int vi = _var_id - Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) return null;
		return titles[vi];
	}
	public int getVariableID(String _var_name) {
		return DataHelper.strings_index(titles, _var_name) + Constants.FIRST_IDX;
	}
	public boolean hasVariable(String _var_name) {
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
	public double getMin(String var_name) { return getMin(DataHelper.strings_index(titles, var_name)); }
	public double getMin(int _var_id) {
		if(_var_id<0 || _var_id>=minmax_mean_sill[0].length) return Double.NaN;
		return minmax_mean_sill[0][_var_id];
	}
	public double[] getMaxs() { return minmax_mean_sill[1]; }
	public double getMax(String var_name) { return getMax(DataHelper.strings_index(titles, var_name)); }
	public double getMax(int _var_id) {
		if(_var_id<0 || _var_id>=minmax_mean_sill[1].length) return Double.NaN;
		return minmax_mean_sill[1][_var_id];
	}
	public double getMean(String var_name) { return getMean(DataHelper.strings_index(titles, var_name)); }
	public double getMean(int _var_id) {
		if(_var_id<0 || _var_id>=minmax_mean_sill[2].length) return Double.NaN;
		return minmax_mean_sill[2][_var_id];
	}
	public double getSill(String var_name) { return getSill(DataHelper.strings_index(titles, var_name)); }
	public double getSill(int _var_id) {
		if(_var_id<0 || _var_id>=minmax_mean_sill[3].length) return Double.NaN;
		return minmax_mean_sill[3][_var_id];
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
		int dlen = (number_of_lines==0 ? datalength[0]+2 : Math.min(number_of_lines+2, datalength[0]+2));
		String[][] out = new String[clen][dlen];
		int[] col_width = new int[clen];
		for(int c=0; c<clen; c++) {
			//System.out.println("Title "+c+" = \""+titles[c]+"\"");
			out[c][0] = " "+titles[c]+" ";
			col_width[c] = titles[c].length()+2;
			switch(types[c]) {
				case BOOL:
					out[c][1] = " bool ";
					boolean[][] ocol = bool_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+(ocol[d][0] ? "true" : "false");
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+(ocol[d][datalength[1]-1] ? "true" : "false")+" ";
					}
					break;
				case BYTE:
					out[c][1] = " byte ";
					byte[][] bcol = byte_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+bcol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+bcol[d][datalength[1]-1]+" ";
					}
					break;
				case INT:
					out[c][1] = " int ";
					int[][] icol = int_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+icol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+icol[d][datalength[1]-1]+" ";
					}
					break;
				case SHORT:
					out[c][1] = " short ";
					short[][] rcol = short_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+rcol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+rcol[d][datalength[1]-1]+" ";
					}
					break;
				case LONG:
					out[c][1] = " long ";
					long[][] lcol = long_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+lcol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+lcol[d][datalength[1]-1]+" ";
					}
					break;
				case FLOAT:
					out[c][1] = " float ";
					float[][] fcol = float_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+fcol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+fcol[d][datalength[1]-1]+" ";
					}
					break;
				case DOUBLE:
					out[c][1] = " double ";
					double[][] dcol = double_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+dcol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+dcol[d][datalength[1]-1]+" ";
					}
					break;
				case STRING:
					out[c][1] = " string ";
					String[][] scol = string_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+scol[d][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+scol[d][datalength[1]-1]+" ";
					}
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
					for(boolean[] oo: bool_column.get(titles[c])) for(boolean o: oo) if(o) { olen++; }
					out[1][cc] = " "+olen+" ";
					break;
				case BYTE: double bmean = 0d;
					out[2][cc] = " NaN "; out[3][cc] = " NaN "; //mean and std
					out[4][cc] = " NaN "; out[5][cc] = " --- "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " NaN ";
					List<Byte> bytes = new ArrayList<Byte>();
					for(byte[] bb: byte_column.get(titles[c])) for(byte b: bb) if(Byte.MIN_VALUE!=b) { bmean += b; bytes.add(b); }
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
					for(int[] ii: int_column.get(titles[c])) for(int i: ii) if(Integer.MIN_VALUE!=i) { imean += i; integers.add(i); }
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
					for(short[] rr: short_column.get(titles[c])) for(short r: rr) if(Short.MIN_VALUE!=r) { rmean += r; shorts.add(r); }
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
					for(long[] ll: long_column.get(titles[c])) for(long l: ll) if(Long.MIN_VALUE!=l) { lmean += l; longs.add(l); }
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
					for(float[] ff: float_column.get(titles[c])) for(float f: ff) if(!Float.isNaN(f)) { fmean += f; floats.add(f); }
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
					for(double[] dd: double_column.get(titles[c])) for(double d: dd) if(!Double.isNaN(d)) { dmean += d; doubles.add(d); }
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
					for(String[] ss: string_column.get(titles[c])) for(String s: ss) if(!(s.equals(" ") || s.length()<1)) { slen++; }
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
		dimension_names = new String[] {"", ""};
		dimension_one = new double[0];
		dimension_two = new double[0];
		minmax_mean_sill = new double[4][0];
		datalength = new int[] {0, 0};
		bool_column.clear();
		byte_column.clear();
		int_column.clear();
		short_column.clear();
		long_column.clear();
		float_column.clear();
		double_column.clear();
		string_column.clear();
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
}
