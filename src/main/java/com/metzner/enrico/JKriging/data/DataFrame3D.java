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
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;
import ucar.nc2.write.NetcdfFileFormat;
import ucar.nc2.write.NetcdfFormatWriter;

public class DataFrame3D {

	private int[] datalength;
	private DataType default_data_type;
	private String[] titles;
	private DataType[] types;
	private String[] dimension_names;
	private double[] dimension_one, dimension_two, dimension_thr;
	private Map<String, String> attribs_dim_one, attribs_dim_two, attribs_dim_thr;
	private double[][] minmax_mean_sill;
	private Map<String, boolean[][][]> bool_column;
	private Map<String, byte[][][]>    byte_column;
	private Map<String, int[][][]>     int_column;
	private Map<String, short[][][]>   short_column;
	private Map<String, long[][][]>    long_column;
	private Map<String, float[][][]>   float_column;
	private Map<String, double[][][]>  double_column;
	private Map<String, String[][][]>  string_column;
	private Map<String, Object>        struct_column;




	public DataFrame3D() {
		datalength = new int[] {0,0,0};
		default_data_type = DataType.STRING;
		titles = new String[0];
		types = new DataType[0];
		int[] currentDate = DataHelper.JulianDate.jd2cal(DataHelper.JulianDate.now());
		String currentDateString = currentDate[0]+FormatHelper.nf(currentDate[1], 2,'0')+FormatHelper.nf(currentDate[2],2,'0');
		currentDateString += FormatHelper.nf(currentDate[3],2,'0')+FormatHelper.nf(currentDate[4],2,'0')+FormatHelper.nf(currentDate[5],2,'0');
		dimension_names = new String[] {"Dim"+currentDateString+"A0","Dim"+currentDateString+"B1","Dim"+currentDateString+"C2"};
		dimension_one = new double[0]; attribs_dim_one = new HashMap<>();
		dimension_two = new double[0]; attribs_dim_two = new HashMap<>();
		dimension_thr = new double[0]; attribs_dim_thr = new HashMap<>();
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
	}




	public void addColumn(String _column_name, boolean[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, boolean[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.BOOL);
		add_stdana(Double.NaN,Double.NaN,Double.NaN,Double.NaN);
		boolean[][][] new_data = new boolean[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++) for(int t1=0; t1<datalength[1]; t1++) for(int t2=0; t2<datalength[2]; t2++)
			if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
				new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
				new_data[t0][t1][t2] = false; }
		bool_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, byte[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, byte[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.BYTE);
		byte[] inax = StdAnalysis.minmax(_column_data); byte[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		byte[][][] new_data = new byte[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Byte.MIN_VALUE; }
		byte_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, short[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, short[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.SHORT);
		short[] inax = StdAnalysis.minmax(_column_data); short[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		short[][][] new_data = new short[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Short.MIN_VALUE; }
		short_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, int[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, int[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.INT);
		int[] inax = StdAnalysis.minmax(_column_data); int[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		int[][][] new_data = new int[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Integer.MIN_VALUE; }
		int_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, long[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, long[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.LONG);
		long[] inax = StdAnalysis.minmax(_column_data); long[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		long[][][] new_data = new long[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Long.MIN_VALUE; }
		long_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, float[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, float[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.FLOAT);
		float[] inax = StdAnalysis.minmax(_column_data); float[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		float[][][] new_data = new float[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Float.NaN; }
		float_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, double[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, double[][][] _column_data, boolean chopORfill) {
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
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.DOUBLE);
		double[] inax = StdAnalysis.minmax(_column_data); double[] mv = StdAnalysis.mean_var(_column_data);
		add_stdana(inax[0],inax[1],mv[0],mv[1]);
		double[][][] new_data = new double[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = Double.NaN; }
		double_column.put(_column_name_c, new_data);
	}
	public void addColumn(String _column_name, String[][][] _column_data) { addColumn(_column_name, _column_data, false); }
	public void addColumn(String _column_name, String[][][] _column_data, boolean chopORfill) {
		String _column_name_c = FormatHelper.name2CFConvention(_column_name);
		if(DataHelper.strings_index(titles, _column_name_c)>=0) {
			System.err.println("A column with title <"+_column_name_c+"> allready exists in this dataframe!");
			return;
		}
		if((_column_data.length!=datalength[0] || _column_data[0].length!=datalength[1] || _column_data[0][0].length!=datalength[2])
				&& !chopORfill && titles.length>0) {
			System.err.println("The new data-column is not compatible with existing data (length is "+
					datalength+" but found "+_column_data.length);
			return;
		}
		if(titles.length==0) { datalength = new int[] {_column_data.length, _column_data[0].length, _column_data[0][0].length};
			dimension_one = DataHelper.createIndexArrayDouble(datalength[0]);
			dimension_two = DataHelper.createIndexArrayDouble(datalength[1]);
			dimension_thr = DataHelper.createIndexArrayDouble(datalength[2]); }
		add_variable(_column_name_c, DataType.STRING); add_stdana(Double.NaN,-Double.NaN,Double.NaN,Double.NaN);
		String[][][] new_data = new String[datalength[0]][datalength[1]][datalength[2]];
		for(int t0=0; t0<datalength[0]; t0++)
			for(int t1=0; t1<datalength[1]; t1++)
				for(int t2=0; t2<datalength[2]; t2++)
					if(t0<_column_data.length && t1<_column_data[0].length && t2<_column_data[0][0].length) {
						new_data[t0][t1][t2] = _column_data[t0][t1][t2]; } else {
						new_data[t0][t1][t2] = ""; }
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
		int di = dim_num-Constants.FIRST_IDX;
		if(di<0 || di>2) {
			System.err.println("Number of dimension must be between 1 and 3 for DataFrame3D!");
			DataHelper.printStackTrace(System.err);
			return;
		}
		boolean change_name   = (_dim_name_c!=null);
		boolean change_values = (dim_values!=null);
		if(change_name)
			if(_dim_name_c.equals(dimension_names[(di+1)%3]) || _dim_name_c.equals(dimension_names[(di+2)%3]) || hasVariable(_dim_name_c)) {
				System.err.println("Cannot rename dimension "+(di+Constants.FIRST_IDX)+", an other dimension or variable has the same name!");
				DataHelper.printStackTrace(System.err);
				return;
			}
		if(change_values)
			if(dim_values.length!=datalength[di] && datalength[di]>0) {
				System.err.println("Number of values must match corresponding dimension length:\n    expect "+datalength[di]+" but got "+dim_values.length);
				DataHelper.printStackTrace(System.err);
				return;
			}
		if(change_name) dimension_names[di] = _dim_name_c;
		if(change_values) {
			if(datalength[di]==0) {
				datalength[di] = dim_values.length;
				if(di==0) dimension_one = new double[dim_values.length];
				if(di==1) dimension_two = new double[dim_values.length];
				if(di==2) dimension_thr = new double[dim_values.length];
			}
			if(di==0)
				for(int dn=0; dn<datalength[0]; dn++) dimension_one[dn] = dim_values[dn];
			if(di==1)
				for(int dn=0; dn<datalength[1]; dn++) dimension_two[dn] = dim_values[dn];
			if(di==2)
				for(int dn=0; dn<datalength[2]; dn++) dimension_thr[dn] = dim_values[dn];
		}
	}
	public void addAttributeToDimension(int dim_num, String name, String value) {
		int di = dim_num - Constants.FIRST_IDX;
		if(di<0 || di>2) {
			System.err.println("Number of dimension must be between "+(Constants.FIRST_IDX)+" and "+
							   (2+Constants.FIRST_IDX)+" for DataFrame3D!");
			DataHelper.printStackTrace(System.err);
			return;
		}
		Map<String,String> att = (di<2 ? (di<1 ? attribs_dim_one : attribs_dim_two) : attribs_dim_thr);
		if(att.containsKey(name))
			System.out.println("[WARNING] override key \""+name+"\" for dimension "+(di+Constants.FIRST_IDX));
		att.put(name, value);
	}
	public void setDimensionsAttributes(int dim_num, Map<String, String> attributes) {
		int di = dim_num - Constants.FIRST_IDX;
		if(di<0 || di>2) {
			System.err.println("Number of dimension must be between "+(Constants.FIRST_IDX)+" and "+
							   (2+Constants.FIRST_IDX)+" for DataFrame3D!");
			DataHelper.printStackTrace(System.err);
			return;
		}
		switch(di) {
			case 0: attribs_dim_one.clear(); attribs_dim_one.putAll(attributes); break;
			case 1: attribs_dim_one.clear(); attribs_dim_one.putAll(attributes); break;
			case 2: attribs_dim_one.clear(); attribs_dim_one.putAll(attributes); break;
		}
	}
	public void copyDimensionsFrom(DataFrame3D otherDF) {
		this.setDimension(0+Constants.FIRST_IDX, otherDF.getDimensionValues(0+Constants.FIRST_IDX), otherDF.getDimensionName(0+Constants.FIRST_IDX));
		this.setDimensionsAttributes(0+Constants.FIRST_IDX, otherDF.getAttributesFromDimension(0+Constants.FIRST_IDX));
		this.setDimension(1+Constants.FIRST_IDX, otherDF.getDimensionValues(1+Constants.FIRST_IDX), otherDF.getDimensionName(1+Constants.FIRST_IDX));
		this.setDimensionsAttributes(1+Constants.FIRST_IDX, otherDF.getAttributesFromDimension(1+Constants.FIRST_IDX));
		this.setDimension(2+Constants.FIRST_IDX, otherDF.getDimensionValues(2+Constants.FIRST_IDX), otherDF.getDimensionName(2+Constants.FIRST_IDX));
		this.setDimensionsAttributes(2+Constants.FIRST_IDX, otherDF.getAttributesFromDimension(2+Constants.FIRST_IDX));
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
	public void setVariableContent(int _var_id, byte[][][] new_bytes) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_bytes.length!=l1 || new_bytes[0].length!=l2 || new_bytes[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case BYTE: byte[][][] arrb = byte_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrb[z][y][x] = new_bytes[z][y][x];
				break;
			case SHORT: short[][][] arrs = short_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrs[z][y][x] = new_bytes[z][y][x];
				break;
			case INT: int[][][] arri = int_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arri[z][y][x] = new_bytes[z][y][x];
				break;
			case LONG: long[][][] arrl = long_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrl[z][y][x] = new_bytes[z][y][x];
				break;
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = new_bytes[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_bytes[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from byte-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, byte[][][] new_bytes) {
		setVariableContent(getVariableID(_var_name), new_bytes);
	}
	public void setVariableContent(int _var_id, short[][][] new_shorts) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_shorts.length!=l1 || new_shorts[0].length!=l2 || new_shorts[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case SHORT: short[][][] arrs = short_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrs[z][y][x] = new_shorts[z][y][x];
				break;
			case INT: int[][][] arri = int_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arri[z][y][x] = new_shorts[z][y][x];
				break;
			case LONG: long[][][] arrl = long_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrl[z][y][x] = new_shorts[z][y][x];
				break;
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = new_shorts[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_shorts[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from short-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, short[][][] new_shorts) {
		setVariableContent(getVariableID(_var_name), new_shorts);
	}
	public void setVariableContent(int _var_id, int[][][] new_ints) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_ints.length!=l1 || new_ints[0].length!=l2 || new_ints[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case INT: int[][][] arri = int_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arri[z][y][x] = new_ints[z][y][x];
				break;
			case LONG: long[][][] arrl = long_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrl[z][y][x] = new_ints[z][y][x];
				break;
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = new_ints[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_ints[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from int-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, int[][][] new_ints) {
		setVariableContent(getVariableID(_var_name), new_ints);
	}
	public void setVariableContent(int _var_id, long[][][] new_longs) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_longs.length!=l1 || new_longs[0].length!=l2 || new_longs[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case LONG: long[][][] arrl = long_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrl[z][y][x] = new_longs[z][y][x];
				break;
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = new_longs[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_longs[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from long-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, long[][][] new_longs) {
		setVariableContent(getVariableID(_var_name), new_longs);
	}
	public void setVariableContent(int _var_id, float[][][] new_floats) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_floats.length!=l1 || new_floats[0].length!=l2 || new_floats[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = new_floats[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_floats[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from float-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, float[][][] new_floats) {
		setVariableContent(getVariableID(_var_name), new_floats);
	}
	public void setVariableContent(int _var_id, double[][][] new_doubles) throws ArrayIndexOutOfBoundsException {
		int ovi = _var_id - Constants.FIRST_IDX;
		if(ovi<0 || ovi>=titles.length) {
			System.err.println("Can not find and set variable with ID "+_var_id+"!");
			return;
		}
		int l1=dimension_one.length,
			l2=dimension_two.length,
			l3=dimension_thr.length;
		if(new_doubles.length!=l1 || new_doubles[0].length!=l2 || new_doubles[0][0].length!=l3)
			throw new ArrayIndexOutOfBoundsException("Size of array does not match size of dataframe-content. Expected shape ("+l1+","+l2+","+l3+").");
		
		String vn = getVarname(ovi);
		DataType ovt = types[ovi];
		switch(types[ovi]) {
			case FLOAT: float[][][] arrf = float_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrf[z][y][x] = (float) new_doubles[z][y][x];
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				for(int z=0; z<l1; z++) for(int y=0; y<l2; y++) for(int x=0; x<l3; x++) arrd[z][y][x] = new_doubles[z][y][x];
				break;
			default: System.err.println("Cannot set content of "+ovt.name()+"-variable from double-array.");
				break;
		}
	}
	public void setVariableContent(String _var_name, double[][][] new_doubles) {
		setVariableContent(getVariableID(_var_name), new_doubles);
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
		if(titles.length==1) { titles = new String[0]; types = new DataType[0]; datalength = new int[] {0, 0, 0}; }
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
			case BOOL: boolean[][][] bool_arr = bool_column.get(variable_name);
				switch(_new_type) {
					case BYTE: byte[][][] barr = new byte[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte)(bool_arr[k][j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = (short)(bool_arr[k][j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = (bool_arr[k][j][i]?1:0);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = (bool_arr[k][j][i]?1L:0L);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = (bool_arr[k][j][i]?1f:0f);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][][] darr = new double[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = (bool_arr[k][j][i]?1d:0d);
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[bool_arr.length][bool_arr[0].length][bool_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = (bool_arr[k][j][i]?"true":"false");
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case BYTE: byte[][][] byte_arr = byte_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][][] oarr = new boolean[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = (byte_arr[k][j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case SHORT: short[][][] rarr = new short[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = byte_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = byte_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = byte_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = byte_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][][] darr = new double[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = byte_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[byte_arr.length][byte_arr[0].length][byte_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = ""+byte_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case SHORT: short[][][] short_arr = short_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][][] oarr = new boolean[short_arr.length][short_arr[0].length][short_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = (short_arr[k][j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte) short_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case INT: int[][][] iarr = new int[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = short_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = short_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = short_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][][] darr = new double[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = short_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[short_arr.length][short_arr[0].length][short_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = ""+short_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case INT: int[][][] int_arr = int_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][][] oarr = new boolean[int_arr.length][int_arr[0].length][int_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = (int_arr[k][j][i]!=0);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte) int_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = (short) int_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case LONG: long[][][] larr = new long[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = int_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = int_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][][] darr = new double[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = int_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[int_arr.length][int_arr[0].length][int_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = ""+int_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case LONG: long[][][] long_arr = long_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][][] oarr = new boolean[long_arr.length][long_arr[0].length][long_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = (long_arr[k][j][i]!=0L ? true : false);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[long_arr.length][long_arr[0].length][long_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte) long_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[long_arr.length][long_arr[0].length][long_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = (short) long_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[long_arr.length][long_arr[0].length][long_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = (int) long_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case DOUBLE: double[][][] darr = new double[long_arr.length][long_arr[0].length][long_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = long_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[long_arr.length][long_arr[0].length][long_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = ""+long_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case FLOAT: float[][][] float_arr = float_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert float to bool result to true except NaN!");
						boolean[][][] oarr = new boolean[float_arr.length][float_arr[0].length][float_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = !Float.isNaN(float_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte) float_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = (short) float_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = (int) float_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = (long) float_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case DOUBLE: double[][][] darr = new double[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = float_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, darr); break;
					case STRING: String[][][] sarr = new String[float_arr.length][float_arr[0].length][float_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++)  sarr[k][j][i] = ""+float_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case DOUBLE: double[][][] double_arr = double_column.get(variable_name);
				switch(_new_type) {
					case BOOL: System.out.println("WARNING: convert double to bool result to true except NaN!");
						boolean[][][] oarr = new boolean[double_arr.length][double_arr[0].length][double_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = !Double.isNaN(double_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = (byte) double_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = (short) double_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = (int) double_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = (long) double_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = (float) double_arr[k][j][i];
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case STRING: String[][][] sarr = new String[double_arr.length][double_arr[0].length][double_arr[0][0].length]; for(int k=0; k<sarr.length; k++)
						for(int j=0; j<sarr[0].length; j++) for(int i=0; i<sarr[0][0].length; i++) sarr[k][j][i] = ""+double_arr[k][j][i]+"";
						removeVariable(variable_name); addColumn(variable_name, sarr); break;
					default:
						System.err.println("An unexpected error occured, cannot change type of variable!");
						DataHelper.printStackTrace(System.err); break;
				} break;
			case STRING: String[][][] string_arr = string_column.get(variable_name);
				switch(_new_type) {
					case BOOL: boolean[][][] oarr = new boolean[string_arr.length][string_arr[0].length][string_arr[0][0].length];for(int k=0; k<oarr.length; k++)
						for(int j=0; j<oarr[0].length; j++) for(int i=0; i<oarr[0][0].length; i++) oarr[k][j][i] = set_boolean(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, oarr); break;
					case BYTE: byte[][][] barr = new byte[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<barr.length; k++)
						for(int j=0; j<barr[0].length; j++) for(int i=0; i<barr[0][0].length; i++) barr[k][j][i] = set_byte(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, barr); break;
					case SHORT: short[][][] rarr = new short[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<rarr.length; k++)
						for(int j=0; j<rarr[0].length; j++) for(int i=0; i<rarr[0][0].length; i++) rarr[k][j][i] = set_short(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, rarr); break;
					case INT: int[][][] iarr = new int[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<iarr.length; k++)
						for(int j=0; j<iarr[0].length; j++) for(int i=0; i<iarr[0][0].length; i++) iarr[k][j][i] = set_int(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, iarr); break;
					case LONG: long[][][] larr = new long[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<larr.length; k++)
						for(int j=0; j<larr[0].length; j++) for(int i=0; i<larr[0][0].length; i++) larr[k][j][i] = set_long(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, larr); break;
					case FLOAT: float[][][] farr = new float[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<farr.length; k++)
						for(int j=0; j<farr[0].length; j++) for(int i=0; i<farr[0][0].length; i++) farr[k][j][i] = set_float(string_arr[k][j][i]);
						removeVariable(variable_name); addColumn(variable_name, farr); break;
					case DOUBLE: double[][][] darr = new double[string_arr.length][string_arr[0].length][string_arr[0][0].length]; for(int k=0; k<darr.length; k++)
						for(int j=0; j<darr[0].length; j++) for(int i=0; i<darr[0][0].length; i++) darr[k][j][i] = set_double(string_arr[k][j][i]);
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
	public void switchDimensions(String dim_name_one, String dim_name_two, String dim_name_three) {
		int d1=-1, d2=-1, d3=-1;
		for(int i=0; i<3; i++) if(dim_name_one.equals(dimension_names[i])) d1=i+Constants.FIRST_IDX;
		if(d1<0) { System.err.println("Unknown dimension name \""+dim_name_one+"\"."); DataHelper.printStackTrace(System.err); }
		for(int i=0; i<3; i++) if(dim_name_two.equals(dimension_names[i])) d2=i+Constants.FIRST_IDX;
		if(d2<0) { System.err.println("Unknown dimension name \""+dim_name_two+"\"."); DataHelper.printStackTrace(System.err); }
		for(int i=0; i<3; i++) if(dim_name_three.equals(dimension_names[i])) d3=i+Constants.FIRST_IDX;
		if(d3<0) { System.err.println("Unknown dimension name \""+dim_name_three+"\"."); DataHelper.printStackTrace(System.err); }
		switchDimensions(d1, d2, d3);
	}
	public void switchDimensions(int dim_id_one, int dim_id_two, int dim_id_three) {
		if(titles.length==0) return;
		int fi = Constants.FIRST_IDX;
		if(dim_id_one<fi || dim_id_one>2+fi ||
			dim_id_two<fi || dim_id_two>2+fi ||
			dim_id_three<fi || dim_id_three>2+fi) {
			System.err.println("Dimension-IDs have to be between "+fi+" and "+(2+fi)+".");
			DataHelper.printStackTrace(System.err);
			return;
		}
		int dimid = (dim_id_one-fi)*100 + (dim_id_two-fi)*10 + (dim_id_three-fi) + 111;
		int uk=1,uj=0,ui=0,  vk=0,vj=1,vi=0,  wk=0,wj=0,wi=1;
		switch(dimid) {
			case 123: return;
			case 132: uk=1; uj=0; ui=0; vk=0; vj=0; vi=1; wk=0; wj=1; wi=0; break;
			case 213: uk=0; uj=1; ui=0; vk=1; vj=0; vi=0; wk=0; wj=0; wi=1; break;
			case 231: uk=0; uj=1; ui=0; vk=0; vj=0; vi=1; wk=1; wj=0; wi=0; break;
			case 312: uk=0; uj=0; ui=1; vk=1; vj=0; vi=0; wk=0; wj=1; wi=0; break;
			case 321: uk=0; uj=0; ui=1; vk=0; vj=1; vi=0; wk=1; wj=0; wi=0; break;
			default:
				System.err.println("Each dimension has to occurre exactly ones!");
				DataHelper.printStackTrace(System.err);
				return;
		}
		int klen = datalength[0];
		int jlen = datalength[1];
		int ilen = datalength[2];
		int ulen = uk*klen + uj*jlen + ui*ilen;
		int vlen = vk*klen + vj*jlen + vi*ilen;
		int wlen = wk*klen + wj*jlen + wi*ilen;
		for(String v: allVariableNames()) {
			switch(getVariableType(v)) {
				case BOOL: boolean[][][] oldo = bool_column.get(v); boolean[][][] newo = new boolean[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newo[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldo[k][j][i];
					bool_column.put(v, newo); break;
				case BYTE: byte[][][] oldb = byte_column.get(v); byte[][][] newb = new byte[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newb[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldb[k][j][i];
					byte_column.put(v, newb); break;
				case SHORT: short[][][] oldr = short_column.get(v); short[][][] newr = new short[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newr[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldr[k][j][i];
					short_column.put(v, newr); break;
				case INT: int[][][] oldi = int_column.get(v); int[][][] newi = new int[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newi[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldi[k][j][i];
					int_column.put(v, newi); break;
				case LONG: long[][][] oldl = long_column.get(v); long[][][] newl = new long[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newl[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldl[k][j][i];
					long_column.put(v, newl); break;
				case FLOAT: float[][][] oldf = float_column.get(v); float[][][] newf = new float[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newf[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldf[k][j][i];
					float_column.put(v, newf); break;
				case DOUBLE: double[][][] oldd = double_column.get(v); double[][][] newd = new double[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) newd[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = oldd[k][j][i];
					double_column.put(v, newd); break;
				case STRING: String[][][] olds = string_column.get(v); String[][][] news = new String[ulen][vlen][wlen];
					for(int k=0; k<klen; k++) for(int j=0; j<jlen; j++) for(int i=0; i<ilen; i++) news[uk*k+uj*j+ui*i][vk*k+vj*j+vi*i][wk*k+wj*j+wi*i] = olds[k][j][i];
					string_column.put(v, news); break;
				default:
					System.err.println("Cannot switch dimensions for variable of type "+getVariableType(v).name()+"!");
					DataHelper.printStackTrace(System.err); break;
			}
		}
		double[] d1 = dimension_one.clone(); String s1 = dimension_names[0];
		double[] d2 = dimension_two.clone(); String s2 = dimension_names[1];
		double[] d3 = dimension_thr.clone(); String s3 = dimension_names[2];
		Map<String,String> att1 = new HashMap<>(); att1.putAll(attribs_dim_one); attribs_dim_one.clear();
		Map<String,String> att2 = new HashMap<>(); att2.putAll(attribs_dim_two); attribs_dim_two.clear();
		Map<String,String> att3 = new HashMap<>(); att3.putAll(attribs_dim_thr); attribs_dim_thr.clear();
		switch(1+dim_id_one-fi) {
			case 1: dimension_one = d1.clone(); dimension_names[0] = s1; attribs_dim_one.putAll(att1); break;
			case 2: dimension_one = d2.clone(); dimension_names[0] = s2; attribs_dim_one.putAll(att2); break;
			case 3: dimension_one = d3.clone(); dimension_names[0] = s3; attribs_dim_one.putAll(att3); break;
			default: break;
		}
		switch(1+dim_id_two-fi) {
			case 1: dimension_two = d1.clone(); dimension_names[1] = s1; attribs_dim_two.putAll(att1); break;
			case 2: dimension_two = d2.clone(); dimension_names[1] = s2; attribs_dim_two.putAll(att2); break;
			case 3: dimension_two = d3.clone(); dimension_names[1] = s3; attribs_dim_two.putAll(att3); break;
			default: break;
		}
		switch(1+dim_id_three-fi) {
			case 1: dimension_thr = d1.clone(); dimension_names[2] = s1; attribs_dim_thr.putAll(att1); break;
			case 2: dimension_thr = d2.clone(); dimension_names[2] = s2; attribs_dim_thr.putAll(att2); break;
			case 3: dimension_thr = d3.clone(); dimension_names[2] = s3; attribs_dim_thr.putAll(att3); break;
			default: break;
		}
		datalength[0] = dimension_one.length;
		datalength[1] = dimension_two.length;
		datalength[2] = dimension_thr.length;
		att1.clear(); att2.clear(); att3.clear();
	}
	
	public DataFrame3D concat(DataFrame3D df2, int dimension_to_append) {
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
		int dimID = dimension_to_append - Constants.FIRST_IDX;
		switch(dimID) {
			case 0: if(dimension_two.length!=df2.getDimensionValues(1+Constants.FIRST_IDX).length) {
				System.err.println("DataFrames' second dimension have not the same length, can't concat them!");
				return this; }
				break;
			case 1: if(dimension_two.length!=df2.getDimensionValues(1+Constants.FIRST_IDX).length) {
				System.err.println("DataFrames' second dimension have not the same length, can't concat them!");
				return this; }
				break;
			case 2: if(dimension_thr.length!=df2.getDimensionValues(2+Constants.FIRST_IDX).length) {
				System.err.println("DataFrames' second dimension have not the same length, can't concat them!");
				return this; }
				break;
			default:
				break;
		}
		for(int iv=0; iv<titles.length; iv++) {
			switch(types[iv]) {
				case BOOL:   bool_column.put(  titles[iv], DataHelper.concat_bool_array(  (boolean[][][])this.getArray(iv), (boolean[][][])df2.getArray(titles[iv]), dimID)); break;
				case BYTE:   byte_column.put(  titles[iv], DataHelper.concat_byte_array(  (byte[][][])this.getArray(iv),    (byte[][][])df2.getArray(titles[iv]),    dimID)); break;
				case SHORT:  short_column.put( titles[iv], DataHelper.concat_short_array( (short[][][])this.getArray(iv),   (short[][][])df2.getArray(titles[iv]),   dimID)); break;
				case INT:    int_column.put(   titles[iv], DataHelper.concat_int_array(   (int[][][])this.getArray(iv),     (int[][][])df2.getArray(titles[iv]),     dimID)); break;
				case LONG:   long_column.put(  titles[iv], DataHelper.concat_long_array(  (long[][][])this.getArray(iv),    (long[][][])df2.getArray(titles[iv]),    dimID)); break;
				case FLOAT:  float_column.put( titles[iv], DataHelper.concat_float_array( (float[][][])this.getArray(iv),   (float[][][])df2.getArray(titles[iv]),   dimID)); break;
				case DOUBLE: double_column.put(titles[iv], DataHelper.concat_double_array((double[][][])this.getArray(iv),  (double[][][])df2.getArray(titles[iv]),  dimID)); break;
				case STRING: string_column.put(titles[iv], DataHelper.concat_string_array((String[][][])this.getArray(iv),  (String[][][])df2.getArray(titles[iv]),  dimID)); break;
				default: break;
			}
		}
		switch(dimID) {
			case 0: datalength[0] += df2.getDimensionValues(dimension_to_append).length;
				dimension_one = DataHelper.concat_double_array(dimension_one, df2.getDimensionValues(dimension_to_append));
				break;
			case 1: datalength[1] += df2.getDimensionValues(dimension_to_append).length;
				dimension_two = DataHelper.concat_double_array(dimension_two, df2.getDimensionValues(dimension_to_append));
				break;
			case 2: datalength[2] += df2.getDimensionValues(dimension_to_append).length;
				dimension_thr = DataHelper.concat_double_array(dimension_thr, df2.getDimensionValues(dimension_to_append));
				break;
			default:
				break;
		}
		return this;
	}
	public DataFrame3D filterByMask(int dimension_to_mask, boolean[] mask) {
		return filterByMask(dimension_to_mask, mask, false);
	}
	public DataFrame3D filterByMask(int dimension_to_mask, boolean[] mask, boolean inline) {
		int dim = dimension_to_mask - Constants.FIRST_IDX;
		DataFrame3D res = new DataFrame3D();
		int count = 0, count2 = 0, count3 = 0;
		for(int i = 0; i < mask.length; i++) {
			if (mask[i])
				count++; 
		}
		int[] ids = new int[count];
		int ct = 0;
		for(int i = 0; i < mask.length; i++) {
			if (mask[i]) {
				ids[ct] = i;
				ct++;
			}
		}
		double[] newDim = new double[count];
		double[] oldDim = getDimensionValues(dim + Constants.FIRST_IDX);
		for (int c = 0; c < count; c++)
			newDim[c] = oldDim[ids[c]]; 
		if(dim==0) {
			count2 = getDimensionValues(1 + Constants.FIRST_IDX).length; 
			count3 = getDimensionValues(2 + Constants.FIRST_IDX).length;
		}
		if(dim==1) {
			count2 = count;
			count = getDimensionValues(0 + Constants.FIRST_IDX).length;
			count3 = getDimensionValues(2 + Constants.FIRST_IDX).length;
		}
		if(dim==2) {
			count3 = count;
			count = getDimensionValues(0 + Constants.FIRST_IDX).length;
			count2 = getDimensionValues(1 + Constants.FIRST_IDX).length;
		}
		for (String var: allVariableNames()) {
			switch (getVariableType(var)) {
				case BOOL: boolean[][][] newBool = new boolean[count][count2][count3], oldBool = bool_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newBool[c][d][e] = oldBool[ids[c]][d][e];
						if(dim==1) newBool[c][d][e] = oldBool[c][ids[d]][e];
						if(dim==2) newBool[c][d][e] = oldBool[c][d][ids[e]];
					} } } if(inline) bool_column.put(var, newBool); else res.addColumn(var, newBool); break;
				case BYTE: byte[][][] newByte = new byte[count][count2][count3], oldByte = byte_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newByte[c][d][e] = oldByte[ids[c]][d][e];
						if(dim==1) newByte[c][d][e] = oldByte[c][ids[d]][e];
						if(dim==2) newByte[c][d][e] = oldByte[c][d][ids[e]];
					} } } if(inline) byte_column.put(var, newByte); else res.addColumn(var, newByte); break;
				case SHORT: short[][][] newShort = new short[count][count2][count3], oldShort = short_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newShort[c][d][e] = oldShort[ids[c]][d][e];
						if(dim==1) newShort[c][d][e] = oldShort[c][ids[d]][e];
						if(dim==2) newShort[c][d][e] = oldShort[c][d][ids[e]];
					} } } if(inline) short_column.put(var, newShort); else res.addColumn(var, newShort); break;
				case INT: int[][][] newInt = new int[count][count2][count3], oldInt = int_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newInt[c][d][e] = oldInt[ids[c]][d][e];
						if(dim==1) newInt[c][d][e] = oldInt[c][ids[d]][e];
						if(dim==2) newInt[c][d][e] = oldInt[c][d][ids[e]];
					} } } if(inline) int_column.put(var, newInt); else res.addColumn(var, newInt); break;
				case LONG: long[][][] newLong = new long[count][count2][count3], oldLong = long_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newLong[c][d][e] = oldLong[ids[c]][d][e];
						if(dim==1) newLong[c][d][e] = oldLong[c][ids[d]][e];
						if(dim==2) newLong[c][d][e] = oldLong[c][d][ids[e]];
					} } } if(inline) long_column.put(var, newLong); else res.addColumn(var, newLong); break;
				case FLOAT: float[][][] newFloat = new float[count][count2][count3], oldFloat = float_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newFloat[c][d][e] = oldFloat[ids[c]][d][e];
						if(dim==1) newFloat[c][d][e] = oldFloat[c][ids[d]][e];
						if(dim==2) newFloat[c][d][e] = oldFloat[c][d][ids[e]];
					} } }  if(inline) float_column.put(var, newFloat); else res.addColumn(var, newFloat); break;
				case DOUBLE: double[][][] newDouble = new double[count][count2][count3], oldDouble = double_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newDouble[c][d][e] = oldDouble[ids[c]][d][e];
						if(dim==1) newDouble[c][d][e] = oldDouble[c][ids[d]][e];
						if(dim==2) newDouble[c][d][e] = oldDouble[c][d][ids[e]];
					} } } if(inline) double_column.put(var, newDouble); else res.addColumn(var, newDouble); break;
				case STRING: String[][][] newString = new String[count][count2][count3], oldString = string_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newString[c][d][e] = oldString[ids[c]][d][e];
						if(dim==1) newString[c][d][e] = oldString[c][ids[d]][e];
						if(dim==2) newString[c][d][e] = oldString[c][d][ids[e]];
					} } } if(inline) string_column.put(var, newString); else res.addColumn(var, newString); break;
				default:
					System.err.println("Unknown datatype, cannot extract and mask variable <" + var + "> from DataFrame2D.");
					break;
			}
		}
		if(inline) {
			datalength[0] = dim==0 ? count  : dimension_one.length;
			datalength[1] = dim==1 ? count2 : dimension_two.length;
			datalength[2] = dim==2 ? count3 : dimension_thr.length;
			if(dim==0) dimension_one = newDim.clone();
			if(dim==1) dimension_two = newDim.clone();
			if(dim==2) dimension_thr = newDim.clone();
			return this;
		} else {
			for(int d=0; d<3; d++) {
				res.setDimension(d+Constants.FIRST_IDX,
						(dim==d) ? newDim : getDimensionValues(d+Constants.FIRST_IDX),
						getDimensionName(d+Constants.FIRST_IDX));
				res.setDimensionsAttributes(d+Constants.FIRST_IDX, getAttributesFromDimension(d+Constants.FIRST_IDX));
			}
			return res;
		}
	}
	public DataFrame3D reorderByIndexList(int dimension_to_sort, int[] ids) {
		return this.reorderByIndexList(dimension_to_sort, ids, false);
	}
	public DataFrame3D reorderByIndexList(int dimension_to_sort, int[] ids, boolean inline) {
		int dim = dimension_to_sort - Constants.FIRST_IDX;
		double[] newDim = new double[ids.length];
		double[] oldDim = getDimensionValues(dimension_to_sort);
		for (int c = 0; c < ids.length; c++)
			newDim[c] = oldDim[ids[c]];
		int count  = dim==0 ? newDim.length : getDimensionValues(0 + Constants.FIRST_IDX).length;
		int count2 = dim==1 ? newDim.length : getDimensionValues(1 + Constants.FIRST_IDX).length;
		int count3 = dim==2 ? newDim.length : getDimensionValues(2 + Constants.FIRST_IDX).length;
		DataFrame3D res = new DataFrame3D();
		for (String var: allVariableNames()) {
			switch (getVariableType(var)) {
				case BOOL: boolean[][][] newBool = new boolean[count][count2][count3], oldBool = bool_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newBool[c][d][e] = oldBool[ids[c]][d][e];
						if(dim==1) newBool[c][d][e] = oldBool[c][ids[d]][e];
						if(dim==2) newBool[c][d][e] = oldBool[c][d][ids[e]];
					} } } if(inline) bool_column.put(var, newBool); else res.addColumn(var, newBool); break;
				case BYTE: byte[][][] newByte = new byte[count][count2][count3], oldByte = byte_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newByte[c][d][e] = oldByte[ids[c]][d][e];
						if(dim==1) newByte[c][d][e] = oldByte[c][ids[d]][e];
						if(dim==2) newByte[c][d][e] = oldByte[c][d][ids[e]];
					} } } if(inline) byte_column.put(var, newByte); else res.addColumn(var, newByte); break;
				case SHORT: short[][][] newShort = new short[count][count2][count3], oldShort = short_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newShort[c][d][e] = oldShort[ids[c]][d][e];
						if(dim==1) newShort[c][d][e] = oldShort[c][ids[d]][e];
						if(dim==2) newShort[c][d][e] = oldShort[c][d][ids[e]];
					} } } if(inline) short_column.put(var, newShort); else res.addColumn(var, newShort); break;
				case INT: int[][][] newInt = new int[count][count2][count3], oldInt = int_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newInt[c][d][e] = oldInt[ids[c]][d][e];
						if(dim==1) newInt[c][d][e] = oldInt[c][ids[d]][e];
						if(dim==2) newInt[c][d][e] = oldInt[c][d][ids[e]];
					} } } if(inline) int_column.put(var, newInt); else res.addColumn(var, newInt); break;
				case LONG: long[][][] newLong = new long[count][count2][count3], oldLong = long_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newLong[c][d][e] = oldLong[ids[c]][d][e];
						if(dim==1) newLong[c][d][e] = oldLong[c][ids[d]][e];
						if(dim==2) newLong[c][d][e] = oldLong[c][d][ids[e]];
					} } } if(inline) long_column.put(var, newLong); else res.addColumn(var, newLong); break;
				case FLOAT: float[][][] newFloat = new float[count][count2][count3], oldFloat = float_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newFloat[c][d][e] = oldFloat[ids[c]][d][e];
						if(dim==1) newFloat[c][d][e] = oldFloat[c][ids[d]][e];
						if(dim==2) newFloat[c][d][e] = oldFloat[c][d][ids[e]];
					} } } if(inline) float_column.put(var, newFloat); else res.addColumn(var, newFloat); break;
				case DOUBLE: double[][][] newDouble = new double[count][count2][count3], oldDouble = double_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newDouble[c][d][e] = oldDouble[ids[c]][d][e];
						if(dim==1) newDouble[c][d][e] = oldDouble[c][ids[d]][e];
						if(dim==2) newDouble[c][d][e] = oldDouble[c][d][ids[e]];
					} } } if(inline) double_column.put(var, newDouble); else res.addColumn(var, newDouble); break;
				case STRING: String[][][] newString = new String[count][count2][count3], oldString = string_column.get(var);
					for(int c=0; c<count; c++) { for(int d=0; d<count2; d++) { for(int e=0; e<count3; e++) {
						if(dim==0) newString[c][d][e] = oldString[ids[c]][d][e];
						if(dim==1) newString[c][d][e] = oldString[c][ids[d]][e];
						if(dim==2) newString[c][d][e] = oldString[c][d][ids[e]];
					} } } if(inline) string_column.put(var, newString); else res.addColumn(var, newString); break;
				default:
					System.err.println("Unknown datatype, cannot extract and mask variable <" + var + "> from DataFrame2D.");
					break;
			}
			if(inline) recalcStats(getVariableID(var));
		}
		if(inline) {
			datalength[dim] = newDim.length;
			setDimension(dimension_to_sort, newDim);
			return this;
		} else {
			res.setDimension(0+Constants.FIRST_IDX, dim==0?newDim:dimension_one, dimension_names[0]);
			res.setDimensionsAttributes(0+Constants.FIRST_IDX, attribs_dim_one);
			res.setDimension(1+Constants.FIRST_IDX, dim==1?newDim:dimension_two, dimension_names[1]);
			res.setDimensionsAttributes(1+Constants.FIRST_IDX, attribs_dim_two);
			res.setDimension(2+Constants.FIRST_IDX, dim==2?newDim:dimension_thr, dimension_names[2]);
			res.setDimensionsAttributes(2+Constants.FIRST_IDX, attribs_dim_thr);
			return res;
		}
	}
	public DataFrame2D toDataFrame2D(String dim_name, int selectedIndex) {
		int dim = -1;
		if(dimension_names[0].equals(dim_name)) dim = 0;
		if(dimension_names[1].equals(dim_name)) dim = 1;
		if(dimension_names[2].equals(dim_name)) dim = 2;
		if(dim<0) {
			System.err.println("DataFrame3D has no dimension \""+dim_name+"\".");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		dim += Constants.FIRST_IDX;
		return toDataFrame2D(dim, selectedIndex);
	}
	public DataFrame2D toDataFrame2D(int dim_id, int selectedIndex) {
		DataFrame2D res2D = new DataFrame2D();
		int dim = dim_id - Constants.FIRST_IDX;
		if(dim<0 || dim>2) {
			System.err.println("DataFrame3D has no dimension "+dim_id+". Index has to be from "+
					(Constants.FIRST_IDX)+" to "+(2+Constants.FIRST_IDX)+".");
			DataHelper.printStackTrace(System.err);
			return res2D;
		}
		int ilen = getDimensionValues(dim_id).length;
		if(selectedIndex-Constants.FIRST_IDX<0 || selectedIndex-Constants.FIRST_IDX>=ilen) {
			System.err.println("selected index out of range!");
			DataHelper.printStackTrace(System.err);
			return res2D;
		}
		//first create new arrays as input res2D from this DataFram3D
		String dim0name = getDimensionName((dim==0?1:0)+Constants.FIRST_IDX);
		double[] dim0 = getDimensionValues((dim==0?1:0)+Constants.FIRST_IDX);
		int dim0len = dim0.length;
		String dim1name = getDimensionName((dim==2?1:2)+Constants.FIRST_IDX);
		double[] dim1 = getDimensionValues((dim==2?1:2)+Constants.FIRST_IDX);
		int dim1len = dim1.length;
		int[] wf = {dim==0 ? 0 : 1,          0,     dim==0 ? selectedIndex-Constants.FIRST_IDX : 0};
		int[] vf = {dim==0 ? 1 : 0, dim==2 ? 1 : 0, dim==1 ? selectedIndex-Constants.FIRST_IDX : 0};
		int[] uf = {         0,     dim==2 ? 0 : 1, dim==2 ? selectedIndex-Constants.FIRST_IDX : 0};
		//check if the DataFrame contain these dimensions!
		//for each other variable than the dimensions create a 3D-array
		for(int iv=0; iv<titles.length; iv++) {
			switch(types[iv]) {
				case BOOL: boolean[][][] oarr = bool_column.get(titles[iv]);
					boolean[][] bool_arr = new boolean[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						bool_arr[a][b] = oarr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], bool_arr);
					break;
				case BYTE: byte[][][] barr = byte_column.get(titles[iv]);
					byte[][] byte_arr = new byte[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						byte_arr[a][b] = barr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], byte_arr);
					break;
				case SHORT: short[][][] rarr = short_column.get(titles[iv]);
					short[][] short_arr = new short[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						short_arr[a][b] = rarr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], short_arr);
					break;
				case INT: int[][][] iarr = int_column.get(titles[iv]);
					int[][] int_arr = new int[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						int_arr[a][b] = iarr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], int_arr);
					break;
				case LONG: long[][][] larr = long_column.get(titles[iv]);
					long[][] long_arr = new long[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						long_arr[a][b] = larr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], long_arr);
					break;
				case FLOAT: float[][][] farr = float_column.get(titles[iv]);
					float[][] float_arr = new float[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						float_arr[a][b] = farr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], float_arr);
					break;
				case DOUBLE: double[][][] darr = double_column.get(titles[iv]);
					double[][] double_arr = new double[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						double_arr[a][b] = darr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], double_arr);
					break;
				case STRING: String[][][] sarr = string_column.get(titles[iv]);
					String[][] string_arr = new String[dim0len][dim1len];
					for(int a=0; a<dim0len; a++) for(int b=0; b<dim1len; b++)
						string_arr[a][b] = sarr[a*wf[0]+b*wf[1]+wf[2]][a*vf[0]+b*vf[1]+vf[2]][a*uf[0]+b*uf[1]+uf[2]];
					res2D.addColumn(titles[iv], string_arr);
					break;
				default:
					System.err.println("An unexpected datatype occure!");
					return new DataFrame2D();
			}
		}
		res2D.setDimension(0+Constants.FIRST_IDX, dim0, dim0name);
		res2D.setDimensionsAttributes(0+Constants.FIRST_IDX, getAttributesFromDimension((dim==0?1:0)+Constants.FIRST_IDX));
		res2D.setDimension(1+Constants.FIRST_IDX, dim1, dim1name);
		res2D.setDimensionsAttributes(1+Constants.FIRST_IDX, getAttributesFromDimension((dim==2?1:2)+Constants.FIRST_IDX));
		return res2D;
	}

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
			int[] datarrlen = { Integer.parseInt(parts[1].trim()), Integer.parseInt(parts[2].trim()), Integer.parseInt(parts[3].trim()) };
			lines.add(""+line);
			for(int i=0; i<nvari; i++) line = br.readLine();
			maxdat = -1;
			while(line!=null) {
				maxdat++;
				line = br.readLine();
			}
			System.out.println("maxdat = "+maxdat+"    with shape "+datarrlen[0]+" x "+datarrlen[1]+" x "+datarrlen[2]);
			if(datarrlen[0]*datarrlen[1]*datarrlen[2]!=maxdat) {
				System.err.println("Shape is not compatible with the number of data entries! ("+datarrlen[0]+"*"+datarrlen[1]+"*"+datarrlen[2]+"!="+maxdat+")");
				return;
			}
			
			names = new String[nvari];
			double[][][][] vr = new double[nvari][datarrlen[0]][datarrlen[1]][datarrlen[2]];
			
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
				nd++; int u=nd%datarrlen[2], v=(nd/datarrlen[2])%datarrlen[1], w=nd/(datarrlen[1]*datarrlen[2]);
				for(int j=0; j<nvari; j++) {
					//vr[j][nd] = Double.parseDouble(parts[j].trim());
					vr[j][w][v][u] = set_double(parts[j].trim()); //with error check for reading errors and NaNs
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
	 * if on variable isn't in the file, it would not be read to this dataframe
	 * and a warning would be printed to the command line
	 * @param filepath path to the NetCDF file
	 * @param variable name(s) of variable(s)
	 * @return         return this DataFrame3D
	 * @throws IOException 
	 */
	public DataFrame3D readFromNetcdf(String filepath, String... variable) throws IOException {
		NetcdfFile nc = NetcdfFiles.open(filepath);
		this.readFromNetcdf(nc, variable);
		nc.close();
		return this;
	}
	/**
	 * Read selected variable data from a Netcdf file
	 * if on variable isn't in the file, it woould not be read to this dataframe
	 * and a warning would be printed to the command line
	 * @param netcdf_file NetCDF object with the data to read
	 * @param variable    name(s) of variable(s)
	 * @return            return this DataFrame3D
	 */
	public DataFrame3D readFromNetcdf(NetcdfFile netcdf_file, String... variable) {
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
		boolean at_least_one_exist = false;
		for(int b=0; b<var_exist.length; b++) {
			Variable var = netcdf_file.findVariable(variable[b]);
			var_exist[b] = (var!=null);
			if(var_exist[b]) {
				if(var.getDimensions().size()>3) {
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
			System.err.println("At least one variable has to have three dimensions or mulitple variables have to share three dimensions together!");
			DataHelper.printStackTrace(System.err);
			return this;
		}
		if(dims.size()!=3) {
			System.err.println("The 3D-dataframe can only handle exact 3 dimensions, but found "+dims.size()+"!");
			DataHelper.printStackTrace(System.err);
			return this;
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
		dimension_thr = new double[dimlen[2]];
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
			Dimension[] vardims = var.getDimensions().toArray(new Dimension[0]);
			if(vardims.length==0) {
				switch(var.getDataType()) {
					case BOOLEAN: boolean bool = a.getBoolean(0); boolean[][][] bool_arr = new boolean[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) bool_arr[w][v][u] = bool;
						addColumn(var.getFullName(), bool_arr); break;
					case BYTE: byte etyb = a.getByte(0); byte[][][] byte_arr = new byte[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) byte_arr[w][v][u] = etyb;
						addColumn(var.getFullName(), byte_arr); break;
					case SHORT: short trohs = a.getShort(0); short[][][] short_arr = new short[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) short_arr[w][v][u] = trohs;
						addColumn(var.getFullName(), short_arr); break;
					case INT: int tni = a.getInt(0); int[][][] int_arr = new int[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) int_arr[w][v][u] = tni;
						addColumn(var.getFullName(), int_arr); break;
					case LONG: long gnol = a.getLong(0); long[][][] long_arr = new long[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) long_arr[w][v][u] = gnol;
						addColumn(var.getFullName(), long_arr); break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float taolf = a.getFloat(0)==f_fill ? Float.NaN : a.getFloat(0); float[][][] float_arr = new float[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) float_arr[w][v][u] = taolf;
						addColumn(var.getFullName(), float_arr); break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double elbuod = a.getDouble(0)==d_fill ? Double.NaN : a.getDouble(0); double[][][] double_arr = new double[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++) double_arr[w][v][u] = elbuod;
						addColumn(var.getFullName(), double_arr); break;
//					case CHAR:
//					case STRING:
//						break;
					case STRUCTURE:
						readStructs(var, var.getFullName(), a, dimlen);
						break;
					default:
						System.out.println("WARNING: does not this datatype: "+var.getDataType().name()+
								", so the variable is not added to the dataframe");
						break;
				}
			} else {
				int uf=0, vf=0, wf=1; Dimension d = var.getDimension(0);
				if(vardims.length==1) {
					if(d.equals(dims.get(1))) { uf=0; vf=1; wf=0; }
					if(d.equals(dims.get(2))) { uf=1; vf=0; wf=0; }
				}
				if(vardims.length==2) {
					Dimension d2 = var.getDimension(1);
					if(d.equals(dims.get(0))) {      // w,?
						if(d2.equals(dims.get(1))) { wf=dimlen[1]; vf=1; uf=0; }   // w,v
						if(d2.equals(dims.get(2))) { wf=dimlen[2]; vf=0; uf=1; } } // w,u
					else if(d.equals(dims.get(1))) { // v,?
						if(d2.equals(dims.get(0))) { wf=1; vf=dimlen[0]; uf=0; }   // v,w
						if(d2.equals(dims.get(2))) { wf=0; vf=dimlen[2]; uf=1; } } // v,u
					else {                           // u,?
						if(d2.equals(dims.get(0))) { wf=1; vf=0; uf=dimlen[0]; }   // u,w
						if(d2.equals(dims.get(1))) { wf=0; vf=1; uf=dimlen[1]; } } // u,v
				}
				if(vardims.length==3) {
					Dimension d2 = var.getDimension(1);
					if(d.equals(dims.get(0))) {      // w,?,?
						if(d2.equals(dims.get(1))) { wf=dimlen[1]*dimlen[2]; vf=dimlen[2]; uf=1; }   // w,v,u
						if(d2.equals(dims.get(2))) { wf=dimlen[2]*dimlen[1]; vf=1; uf=dimlen[1]; } } // w,u,v
					else if(d.equals(dims.get(1))) { // v,?,?
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
						addColumn(var.getFullName(), bool_arr);
						arr_bool = null; bool_arr = null; break;
					case BYTE: byte[] arr_byte = (byte[]) a.get1DJavaArray(ucar.ma2.DataType.BYTE);
						byte[][][] byte_arr = new byte[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							byte_arr[w][v][u] = arr_byte[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), byte_arr);
						arr_byte = null; byte_arr = null; break;
					case SHORT: short[] arr_short = (short[]) a.get1DJavaArray(ucar.ma2.DataType.SHORT);
						short[][][] short_arr = new short[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							short_arr[w][v][u] = arr_short[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), short_arr);
						arr_short = null; short_arr = null; break;
					case INT: int[] arr_int = (int[]) a.get1DJavaArray(ucar.ma2.DataType.INT);
						int[][][] int_arr = new int[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							int_arr[w][v][u] = arr_int[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), int_arr);
						arr_int = null; int_arr = null; break;
					case LONG: long[] arr_long = (long[]) a.get1DJavaArray(ucar.ma2.DataType.LONG);
						long[][][] long_arr = new long[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							long_arr[w][v][u] = arr_long[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), long_arr);
						arr_long = null; long_arr = null; break;
					case FLOAT: float f_fill = Float.NaN; if(var.hasAttribute("_FillValue")) f_fill = (float) var.findAttribute("_FillValue").getNumericValue();
						float[] arr_float = (float[]) a.get1DJavaArray(ucar.ma2.DataType.FLOAT);
						for(int ff=0; ff<arr_float.length; ff++) if(arr_float[ff]==f_fill) arr_float[ff] = Float.NaN;
						float[][][] float_arr = new float[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							float_arr[w][v][u] = arr_float[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), float_arr);
						arr_float = null; float_arr = null; break;
					case DOUBLE: double d_fill = Double.NaN; if(var.hasAttribute("_FillValue")) d_fill = (double) var.findAttribute("_FillValue").getNumericValue();
						double[] arr_double = (double[]) a.get1DJavaArray(ucar.ma2.DataType.DOUBLE);
						for(int dd=0; dd<arr_double.length; dd++) if(arr_double[dd]==d_fill) arr_double[dd] = Double.NaN;
						double[][][] double_arr = new double[dimlen[0]][dimlen[1]][dimlen[2]];
						for(int w=0; w<dimlen[0]; w++) for(int v=0; v<dimlen[1]; v++) for(int u=0; u<dimlen[2]; u++)
							double_arr[w][v][u] = arr_double[w*wf+v*vf+u*uf];
						addColumn(var.getFullName(), double_arr);
						arr_double = null; double_arr = null; break;
//					case CHAR:
//					case STRING:
//						break;
					case STRUCTURE:
						readStructs(var, var.getFullName(), a, dimlen, uf, vf, wf);
						break;
					default:
						System.out.println("WARNING: does not this datatype: "+var.getDataType().name()+
								", so the variable is not added to the dataframe");
						break;
				}
			}
		}
		for(int dimid=0; dimid<3; dimid++) {
			Variable var = netcdf_file.findVariable(dims.get(dimid).getName());
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
					if(dimid==2)
						for(int dl=0; dl<dimlen[2]; dl++)
							dimension_thr[dl] = a.getDouble(dl);
					dimension_names[dimid] = var.getFullName();
					//double[] dimarray = ( dimid==0 ? dimension_one : (dimid==1 ? dimension_two : dimension_thr));
					//String dim_vals = ""+dimarray[0];
					//for(int dl=1; dl<dimarray.length; dl++)
					//	dim_vals += ", "+dimarray[dl];
					//System.out.println("[DF3D - NDread] read dimension \""+dimension_names[dimid]+"\": "+dim_vals);
					switch(dimid) {
						case 0: attribs_dim_one.clear(); break;
						case 1: attribs_dim_two.clear(); break;
						case 2: attribs_dim_thr.clear(); break;
						default: break;
					}
					for(Attribute att: var.attributes()) {
						switch(dimid) {
							case 0: attribs_dim_one.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue())); break;
							case 1: attribs_dim_two.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue())); break;
							case 2: attribs_dim_thr.put(att.getName(), ""+(att.isString()?att.getStringValue():att.getNumericValue())); break;
							default: break;
						}
					}
					was_succesful = true;
				} catch (IOException e) {
					System.out.println("WARNING: could not read dimension \""+dims.get(dimid).getName()+"\": add as index-array to dataframe");
				} catch (NumberFormatException nfe) {
					System.out.println("WARNING: could not read dimension \""+dims.get(dimid).getName()+"\": add as index-array to dataframe");
				}
			}
			if(!was_succesful) {
				if(dimid==0) dimension_one = DataHelper.createIndexArrayDouble(dimension_one.length);
				if(dimid==1) dimension_two = DataHelper.createIndexArrayDouble(dimension_two.length);
				if(dimid==2) dimension_thr = DataHelper.createIndexArrayDouble(dimension_thr.length);
			}
		}
		return this;
	}
	public void readStructs(Variable var, String varname, Array a, int[] dimlen, int... coeff) {
		System.out.println("WARNING: could not add variable \""+varname+"\": unsupported data type!");
	}
	/**
	 * Write all content of this dataframe to a netcdf file
	 * @param netcdf_file_path path to the netcdf file
	 * @throws IOException 
	 * @throws InvalidRangeException 
	 */
	@SuppressWarnings("rawtypes")
	public void writeToNetcdf(String netcdf_file_path) throws IOException, InvalidRangeException {
		NetcdfFormatWriter.Builder ncdfWriter = NetcdfFormatWriter.createNewNetcdf4(NetcdfFileFormat.NETCDF4_CLASSIC, netcdf_file_path, null);
		Dimension dim0 = ncdfWriter.addDimension(dimension_names[0], datalength[0]);
		Dimension dim1 = ncdfWriter.addDimension(dimension_names[1], datalength[1]);
		Dimension dim2 = ncdfWriter.addDimension(dimension_names[2], datalength[2]);
		List<Dimension> dims = new ArrayList<>(); dims.add(dim0);
		Variable.Builder[] vars = new Variable.Builder[titles.length+3];
		vars[0] = ncdfWriter.addVariable(dimension_names[0], ucar.ma2.DataType.DOUBLE, dims);
		for(String k: attribs_dim_one.keySet())
			if(k.equals("_FillValue"))
				vars[0].addAttribute(new Attribute(k, Double.parseDouble(attribs_dim_one.get(k))));
			else
				vars[0].addAttribute(new Attribute(k, attribs_dim_one.get(k)));
		dims.remove(0); dims.add(dim1);
		vars[1] = ncdfWriter.addVariable(dimension_names[1], ucar.ma2.DataType.DOUBLE, dims);
		for(String k: attribs_dim_two.keySet())
			if(k.equals("_FillValue"))
				vars[1].addAttribute(new Attribute(k, Double.parseDouble(attribs_dim_two.get(k))));
			else
				vars[1].addAttribute(new Attribute(k, attribs_dim_two.get(k)));
		dims.remove(0); dims.add(dim2);
		vars[2] = ncdfWriter.addVariable(dimension_names[2], ucar.ma2.DataType.DOUBLE, dims);
		for(String k: attribs_dim_thr.keySet())
			if(k.equals("_FillValue"))
				vars[2].addAttribute(new Attribute(k, Double.parseDouble(attribs_dim_thr.get(k))));
			else
				vars[2].addAttribute(new Attribute(k, attribs_dim_thr.get(k)));
		dims.remove(0); dims.add(dim0); dims.add(dim1); dims.add(dim2);
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
			vars[iv+3] = ncdfWriter.addVariable(FormatHelper.name2CFConvention(titles[iv]), dataType, dims);
			switch(types[iv]) {
				case SHORT:  vars[iv+3].addAttribute(new Attribute("_FillValue", Short.MIN_VALUE));   break;
				case INT:    vars[iv+3].addAttribute(new Attribute("_FillValue", Integer.MIN_VALUE)); break;
				case LONG:   vars[iv+3].addAttribute(new Attribute("_FillValue", Long.MIN_VALUE));    break;
				case FLOAT:  vars[iv+3].addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_F));  break;
				case DOUBLE: vars[iv+3].addAttribute(new Attribute("_FillValue", Constants.FILL_VALUE_D)); break;
				default: break;
			}
		}
		ncdfWriter.getRootGroup().addAttribute(new Attribute("history", "created with \"EnMe-Kriging\" from Dataframe3D - JAVA Netcdf "+Constants.NETCDF_VERSION));
		NetcdfFormatWriter ncdf = ncdfWriter.build();
		Variable dim0Var = ncdf.findVariable(vars[0].getFullName());
		ArrayDouble.D1 dim0_arr = new ArrayDouble.D1(datalength[0]);
		for(int dl=0; dl<datalength[0]; dl++) dim0_arr.set(dl,dimension_one[dl]);
		ncdf.write(dim0Var, dim0_arr);
		Variable dim1Var = ncdf.findVariable(vars[1].getFullName());
		ArrayDouble.D1 dim1_arr = new ArrayDouble.D1(datalength[1]);
		for(int dl=0; dl<datalength[1]; dl++) dim1_arr.set(dl,dimension_two[dl]);
		ncdf.write(dim1Var, dim1_arr);
		Variable dim2Var = ncdf.findVariable(vars[2].getFullName());
		ArrayDouble.D1 dim2_arr = new ArrayDouble.D1(datalength[2]);
		for(int dl=0; dl<datalength[2]; dl++) dim2_arr.set(dl,dimension_thr[dl]);
		ncdf.write(dim2Var, dim2_arr);
		for(int iv=0; iv<titles.length; iv++) {
			switch(types[iv]) {
				case BOOL:
					ArrayBoolean.D3 bool_arr = new ArrayBoolean.D3(datalength[0],datalength[1],datalength[2]);
					boolean[][][] bool_source = bool_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								bool_arr.set(dl0,dl1,dl2, bool_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), bool_arr);
					break;
				case BYTE:
					ArrayByte.D3 byte_arr = new ArrayByte.D3(datalength[0],datalength[1],datalength[2],false);
					byte[][][] byte_source = byte_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								byte_arr.set(dl0,dl1,dl2, byte_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), byte_arr);
					break;
				case SHORT:
					ArrayShort.D3 short_arr = new ArrayShort.D3(datalength[0],datalength[1],datalength[2],false);
					short[][][] short_source = short_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								short_arr.set(dl0,dl1,dl2, short_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), short_arr);
					break;
				case INT:
					ArrayInt.D3 int_arr = new ArrayInt.D3(datalength[0],datalength[1],datalength[2],false);
					int[][][] int_source = int_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								int_arr.set(dl0,dl1,dl2, int_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), int_arr);
					break;
				case LONG:
					ArrayLong.D3 long_arr = new ArrayLong.D3(datalength[0],datalength[1],datalength[2],false);
					long[][][] long_source = long_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								long_arr.set(dl0,dl1,dl2, long_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), long_arr);
					break;
				case FLOAT:
					ArrayFloat.D3 float_arr = new ArrayFloat.D3(datalength[0],datalength[1],datalength[2]);
					float[][][] float_source = float_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								float_arr.set(dl0,dl1,dl2, Float.isNaN(float_source[dl0][dl1][dl2])?Constants.FILL_VALUE_F:float_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), float_arr);
					break;
				case DOUBLE:
					ArrayDouble.D3 double_arr = new ArrayDouble.D3(datalength[0],datalength[1],datalength[2]);
					double[][][] double_source = double_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								double_arr.set(dl0,dl1,dl2, Double.isNaN(double_source[dl0][dl1][dl2])?Constants.FILL_VALUE_D:double_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), double_arr);
					break;
				case STRING:
					ArrayString.D3 string_arr = new ArrayString.D3(datalength[0],datalength[1],datalength[2]);
					String[][][] string_source = string_column.get(titles[iv]);
					for(int dl0=0; dl0<datalength[0]; dl0++)
						for(int dl1=0; dl1<datalength[1]; dl1++)
							for(int dl2=0; dl2<datalength[2]; dl2++)
								string_arr.set(dl0,dl1,dl2, string_source[dl0][dl1][dl2]);
					ncdf.write(ncdf.findVariable(vars[iv+3].getFullName()), string_arr);
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
	 * @throws InvalidRangeException 
	 * @throws IOException 
	 */
	public static void writeToNetcdf(String path, Object... df) throws IOException, InvalidRangeException {
		DataFrame.writeToNetcdf(path, df);
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


	public double[] getDimensionValues(int dim_id) {
		int di = dim_id - Constants.FIRST_IDX;
		if(di<0 || di>2) return null;
		return (di==0 ? dimension_one : di==1 ? dimension_two : dimension_thr);
	}
	public String getDimensionName(int dim_id) {
		int di = dim_id - Constants.FIRST_IDX;
		if(di<0 || di>2) return null;
		return dimension_names[di];
	}
	public Map<String, String> getAttributesFromDimension(int dim_num) {
		int di = dim_num - Constants.FIRST_IDX;
		if(di<0 || di>2) {
			System.err.println("Number of dimension must be between "+(Constants.FIRST_IDX)+" and "+
							   (2+Constants.FIRST_IDX)+" for DataFrame3D!");
			DataHelper.printStackTrace(System.err);
			return null;
		}
		switch(di) {
			case 0: return attribs_dim_one;
			case 1: return attribs_dim_two;
			case 2: return attribs_dim_thr;
			default: return null;
		}
	}
	public int getVariableCount() { return titles.length; }
	public String getVarname(int _var_id) {
		int vi = _var_id-Constants.FIRST_IDX;
		if(vi<0 || vi>=titles.length) return null;
		return titles[vi];
	}
	public int getVariableID(String _var_name) {
		return DataHelper.strings_index(titles, _var_name);
	}
	public boolean hasVariable(String _var_name) {
		return getVariableID(_var_name)>=0;
	}
	public String[] allVariableNames() { return titles; }
	public DataType getVariableType(int _var_id) {
		if(_var_id<0 || _var_id>=titles.length) return null;
		return types[_var_id];
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
					boolean[][][] ocol = bool_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+(ocol[d][0][0] ? "true" : "false");
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+(ocol[d][datalength[1]-1][0] ? "true" : "false")+" ";
					}
					break;
				case BYTE:
					out[c][1] = " byte ";
					byte[][][] bcol = byte_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+bcol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+bcol[d][datalength[1]-1][0]+" ";
					}
					break;
				case INT:
					out[c][1] = " int ";
					int[][][] icol = int_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+icol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+icol[d][datalength[1]-1][0]+" ";
					}
					break;
				case SHORT:
					out[c][1] = " short ";
					short[][][] rcol = short_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+rcol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+rcol[d][datalength[1]-1][0]+" ";
					}
					break;
				case LONG:
					out[c][1] = " long ";
					long[][][] lcol = long_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+lcol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+lcol[d][datalength[1]-1][0]+" ";
					}
					break;
				case FLOAT:
					out[c][1] = " float ";
					float[][][] fcol = float_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+fcol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+fcol[d][datalength[1]-1][0]+" ";
					}
					break;
				case DOUBLE:
					out[c][1] = " double ";
					double[][][] dcol = double_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+dcol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+dcol[d][datalength[1]-1][0]+" ";
					}
					break;
				case STRING:
					out[c][1] = " string ";
					String[][][] scol = string_column.get(titles[c]);
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " "+scol[d][0][0];
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+scol[d][datalength[1]-1][0]+" ";
					}
					break;
				case STRUCT:
					out[c][1] = " struct ";
					for(int d=0; d+2<dlen; d++) {
						out[c][d+2] = " ???";
						if(datalength[1]>1) out[c][d+2] += (datalength[1]==2 ? ", " : " ... ")+"??? ";
					}
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
		String[][] out = new String[10][clen+1];
		out[0][0] = " "; out[1][0] = " TYPE "; out[2][0] = " Count "; out[3][0] = " mean ";
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
					for(boolean[][] ooo: bool_column.get(titles[c])) for(boolean[] oo: ooo) for(boolean o: oo) if(o) { olen++; }
					out[2][cc] = " "+olen+" ";
					break;
				case BYTE: double bmean = 0d;
					out[3][cc] = " NaN "; out[4][cc] = " NaN "; //mean and std
					out[5][cc] = " NaN "; out[6][cc] = " --- "; out[7][cc] = " --- "; out[8][cc] = " --- "; out[9][cc] = " NaN ";
					List<Byte> bytes = new ArrayList<Byte>();
					for(byte[][] bbb: byte_column.get(titles[c])) for(byte[] bb: bbb) for(byte b: bb) if(Byte.MIN_VALUE!=b) { bmean += b; bytes.add(b); }
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
					for(int[][] iii: int_column.get(titles[c])) for(int[] ii: iii) for(int i: ii) if(Integer.MIN_VALUE!=i) { imean += i; integers.add(i); }
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
					for(short[][] rrr: short_column.get(titles[c])) for(short[] rr: rrr) for(short r: rr) if(Short.MIN_VALUE!=r) { rmean += r; shorts.add(r); }
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
					for(long[][] lll: long_column.get(titles[c])) for(long[] ll: lll) for(long l: ll) if(Long.MIN_VALUE!=l) { lmean += l; longs.add(l); }
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
					for(float[][] fff: float_column.get(titles[c])) for(float[] ff: fff) for(float f: ff) if(!Float.isNaN(f)) { fmean += f; floats.add(f); }
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
					for(double[][] ddd: double_column.get(titles[c])) for(double[] dd: ddd) for(double d: dd) if(!Double.isNaN(d)) { dmean += d; doubles.add(d); }
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
					for(String[][] sss: string_column.get(titles[c])) for(String[] ss: sss) for(String s: ss) if(!(s.equals(" ") || s.length()<1)) { slen++; }
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
	public void describeDimensions() {
		String str1 = "(1/3): "+dimension_names[0]+" ["+dimension_one.length+"] ";
		double[] rng1 = StdAnalysis.minmax(dimension_one);
		double[] std1 = StdAnalysis.mean_var(dimension_one);
		str1 += "{"+rng1[0]+" ... "+rng1[1]+", mean="+std1[0]+", var="+std1[1]+"}";
		System.out.println(str1);
		String str2 = "(2/3): "+dimension_names[1]+" ["+dimension_two.length+"] ";
		double[] rng2 = StdAnalysis.minmax(dimension_two);
		double[] std2 = StdAnalysis.mean_var(dimension_two);
		str2 += "{"+rng2[0]+" ... "+rng2[1]+", mean="+std2[0]+", var="+std2[1]+"}";
		System.out.println(str2);
		String str3 = "(3/3): "+dimension_names[2]+" ["+dimension_thr.length+"] ";
		double[] rng3 = StdAnalysis.minmax(dimension_thr);
		double[] std3 = StdAnalysis.mean_var(dimension_thr);
		str3 += "{"+rng3[0]+" ... "+rng3[1]+", mean="+std3[0]+", var="+std3[1]+"}";
		System.out.println(str3);
	}
	public void clear() {
		titles = null;
		types = null;
		titles = new String[0];
		types = new DataType[0];
		dimension_names = new String[] {"", "", ""};
		dimension_one = new double[0];
		dimension_two = new double[0];
		dimension_thr = new double[0];
		attribs_dim_one.clear();
		attribs_dim_two.clear();
		attribs_dim_thr.clear();
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

	@SuppressWarnings("unused")
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
			case BYTE: byte[][][] arrb = byte_column.get(vn);
				byte[] inaxb = StdAnalysis.minmax(arrb); byte[] mvb = StdAnalysis.mean_var(arrb);
				set_stdana(_var_id, inaxb[0], inaxb[1], mvb[0], mvb[1]);
				break;
			case SHORT: short[][][] arrs = short_column.get(vn);
				short[] inaxs = StdAnalysis.minmax(arrs); short[] mvs = StdAnalysis.mean_var(arrs);
				set_stdana(_var_id, inaxs[0], inaxs[1], mvs[0], mvs[1]);
				break;
			case INT: int[][][] arri = int_column.get(vn);
				int[] inaxi = StdAnalysis.minmax(arri); int[] mvi = StdAnalysis.mean_var(arri);
				set_stdana(_var_id, inaxi[0], inaxi[1], mvi[0], mvi[1]);
				break;
			case LONG: long[][][] arrl = long_column.get(vn);
				long[] inaxl = StdAnalysis.minmax(arrl); long[] mvl = StdAnalysis.mean_var(arrl);
				set_stdana(_var_id, inaxl[0], inaxl[1], mvl[0], mvl[1]);
				break;
			case FLOAT: float[][][] arrf = float_column.get(vn);
				float[] inaxf = StdAnalysis.minmax(arrf); float[] mvf = StdAnalysis.mean_var(arrf);
				set_stdana(_var_id, inaxf[0], inaxf[1], mvf[0], mvf[1]);
				break;
			case DOUBLE: double[][][] arrd = double_column.get(vn);
				double[] inaxd = StdAnalysis.minmax(arrd); double[] mvd = StdAnalysis.mean_var(arrd);
				set_stdana(_var_id, inaxd[0], inaxd[1], mvd[0], mvd[1]);
				break;
			default:
				break;
		}
	}
}
