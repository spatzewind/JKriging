package com.metzner.enrico.JKriging.data.reader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame.DataType;

import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;

public class ExcelReader extends DataReader {
	
	private File wbfile = null;
	private ReadableWorkbook wb = null;
	private List<ExcelVariable> all_vars;
	private List<Sheet> orderedSheets;
	private Map<Sheet,ExcelVariable[]> var_table;
	
	public ExcelReader(File file) throws IOException {
		super(file);
		orderedSheets = new ArrayList<Sheet>();
		all_vars = new ArrayList<>();
		var_table = new HashMap<>();
		wb = new ReadableWorkbook(file);
		wbfile = file;
		retrieveVarTable();
	}
	
	@Override
	public void describeContent() {
		System.out.println("MS Excel file");
		System.out.println(wT+wL+"Path: "+wbfile.getParent());
		System.out.println(wT+wL+"Name: "+wbfile.getName());
		boolean isLast = false;
		for(int s=0; s<100 && !isLast; s++) {
			Optional<Sheet> next = wb.getSheet(s+1);
			isLast = !next.isPresent();
			describeSheet(wb.getSheet(s).get(), isLast?wE:wT, isLast?" ":wP);
		}
//		describeGroup(ncfile.getRootGroup(),wT,wP);
//		describeAttributes(ncfile.getGlobalAttributes(),wE," ");
	}
	
	@Override
	public Map<String, Integer> getDimensionsOfVariable(String var) throws IllegalAccessException {
		return null;
	}
	
	@Override
	public DataFrame getVars1D(String... vars) throws IllegalAccessException {
		if(wbfile==null) throw new IllegalAccessException("The Excel file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		ExcelVariable[] varlist = new ExcelVariable[vars.length];
		for(int i=0; i<vars.length; i++) {
			for(Sheet s: orderedSheets)
				for(ExcelVariable v: var_table.get(s)) {
					if(v==null || v.name==null)
						continue;
					if(!vars[i].equals(v.name))
						continue;
					if(varlist[i]==null) {
						varlist[i] = v;
					} else {
						System.out.println("[Warning] The variable "+vars[i]+" will only read from sheet "+s.getName()+".");
					}
			}
			if(varlist[i]==null)
				throw new IllegalStateException("The Excel sheet does not have a field called "+vars[i]);
		}
		//collect dimensions
		List<JKDim> dimlist = new ArrayList<JKDim>();
		for(ExcelVariable v: varlist)
			for(JKDim d: v.dims)
				if(!dimlist.contains(d))
					dimlist.add(d);
		JKDim[] dims = dimlist.toArray(new JKDim[0]);
		System.out.println("[DEBUG] collected dims: "+Arrays.toString(dims));
		if(dims.length!=1) {
			double logDimSize = 0d;
			long dimSize = 1;
			for(JKDim d: dims) {
				logDimSize += Math.log(d.len);
				dimSize *= d.len;
			}
			if(logDimSize>=Math.log(DataReader.MAX_DATAFRAME_LENGTH))
				throw new ArrayStoreException("Cannot expand combination of dimensions, maximum length of dataframe reached!");
			else
				System.out.println("[Warning] combine multiple dimension into one long dimension of size "+dimSize+".");
		}
		
		int datalength = 1;
		int[] dimlen = new int[dims.length];
		for(int d=0; d<dims.length; d++) { dimlen[d] = dims[d].len; datalength *= dimlen[d]; }
		int[][] indices = new int[datalength][dims.length];
		for(int dl=0; dl<datalength; dl++) {
			int sublen = 1;
			for(int d=dims.length-1; d>=0; d--) {
				indices[dl][d] = (dl / sublen) % dimlen[d];
				sublen *= dimlen[d];
			}
		}
//		System.out.println("[DEBUG]");
//		FormatHelper.printTable(2, indices);
		DataFrame df = new DataFrame();
		if(dims.length>1) {
			for(int di=0; di<dims.length; di++) {
				int[] coordinates = new int[datalength];
				for(int dl=0; dl<datalength; dl++)
					coordinates[dl] = indices[dl][di];
				df.addColumn(dims[di].name, coordinates);
			}
		}
		@SuppressWarnings("unchecked")
		Optional<Cell>[] cells = new Optional[datalength];
		for(ExcelVariable var: varlist) {
			if(df.hasVariable(var.name))
				continue;
			List<Row> var_data_in_rows;
			try {
				var_data_in_rows = var.sheet.read();
			} catch (IOException ioe) {
				throw (IllegalAccessException) new IllegalAccessException("Cannot read from the Excel file.").initCause(ioe);
			}
			for(int dl=0; dl<datalength; dl++) cells[dl] = null;
			int[] dimids = new int[var.dims.length];
			for(int d=0; d<dimids.length; d++) {
				for(int dd=0; dd<dims.length; dd++)
					if (dims[dd].equals(var.dims[d]))
						dimids[d] = dd;
			}
			for(int dl=0; dl<datalength; dl++) {
				int colIdx = var.columnID;
				int rowIdx = var.title_row;
				if(var.isWholeSheet) {
					colIdx += indices[dl][dimids[0]];
					rowIdx += indices[dl][dimids[1]];
				} else {
					rowIdx += indices[dl][dimids[0]];
				}
				cells[dl] = var_data_in_rows.get(rowIdx).getOptionalCell(colIdx);
			}
			switch(var.type) {
				case BOOL: boolean[] read_bools = new boolean[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_bools[dl] = cells[dl].get().asBoolean().booleanValue(); else
					read_bools[dl] = false; df.addColumn(var.name, read_bools); break;
				case BYTE:
					byte[] read_bytes = new byte[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_bytes[dl] = cells[dl].get().asNumber().byteValue(); else
					read_bytes[dl] = 0b0; df.addColumn(var.name, read_bytes); break;
				case SHORT:
					short[] read_shorts = new short[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_shorts[dl] = cells[dl].get().asNumber().shortValue(); else
					read_shorts[dl] = 0; df.addColumn(var.name, read_shorts); break;
				case INT:
					int[] read_ints = new int[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_ints[dl] = cells[dl].get().asNumber().intValue(); else
					read_ints[dl] = 0; df.addColumn(var.name, read_ints); break;
				case LONG:
					long[] read_longs = new long[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_longs[dl] = cells[dl].get().asNumber().longValue(); else
					read_longs[dl] = 0L; df.addColumn(var.name, read_longs); break;
				case FLOAT:
					float[] read_floats = new float[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_floats[dl] = cells[dl].get().asNumber().floatValue(); else
					read_floats[dl] = Float.NaN; df.addColumn(var.name, read_floats); break;
				case DOUBLE:
					double[] read_doubles = new double[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_doubles[dl] = cells[dl].get().asNumber().doubleValue(); else
					read_doubles[dl] = Double.NaN; df.addColumn(var.name, read_doubles); break;
				case STRING:
					String[] read_strings = new String[datalength]; for(int dl=0; dl<datalength; dl++)
					if(cells[dl].isPresent()) read_strings[dl] = cells[dl].get().getText(); else
					read_strings[dl] = ""; df.addColumn(var.name, read_strings); break;
				default:
					System.err.println("WARNING: could not add variable \""+var.name+"\": unsupported data type!");
					break;
			}
		}
		if(dims.length==1) {
			df.setDimension(dims[0].name);
		}
		return df;
	}
	public DataFrame2D getVars2D(String... vars) throws IllegalAccessException {
		if(wbfile==null) throw new IllegalAccessException("The Excel file does not exist.");
		if(vars==null || vars.length==0) throw new IllegalArgumentException("At least one variable has to be specified.");
		ExcelVariable[] varlist = new ExcelVariable[vars.length];
		for(int i=0; i<vars.length; i++) {
			sheetloop: for(Sheet s: var_table.keySet())
				for(ExcelVariable v: var_table.get(s)) {
					if(v!=null && v.name!=null)
						if(vars[i].equals(v.name)) {
							varlist[i] = v;
							break sheetloop;
						}
			}
			if(varlist[i]==null)
				throw new IllegalStateException("The Excel sheet does not have a field called "+vars[i]);
		}
		return null;
	}
	public DataFrame3D getVars3D(String... vars) throws IllegalAccessException {
		throw new IllegalAccessException("Cannot read 3D variables from an Excel Sheet.");
	}
	public DataFrame get1Dslice(String filtervar, int index, String... vars) throws IllegalAccessException {
		return null;
	}
	
	@Override
	public void close() throws IOException {
		all_vars.clear();
		var_table.clear();
		if(wb!=null) wb.close();
		wb = null;
	}
	
	protected static boolean hasValidExtension(String filename) {
		if(filename.endsWith(".xls")) return true;
		if(filename.endsWith(".xlsx")) return true;
		return false;
	}
	
	private Optional<Cell> getLastNonEmptyCell(Row row) {
		return row.stream()
				.filter(Objects::nonNull)
				.filter(cell -> !cell.getText().isEmpty())
				.reduce((f,s)->s);
	}
	private DataType getNumberType(List<Optional<Cell>> cells) {
		String[] s = new String[cells.size()];
		for(int i=0; i<s.length; i++)
			s[i] = !cells.get(i).isPresent()?"":cells.get(i).get().getText().trim();
		//check long with max and min value:
		try {
			long minV = Long.MAX_VALUE, maxV = Long.MIN_VALUE;
			for(int i=0; i<s.length; i++)
				if(s[i].length()>0) {
					long l = Long.parseLong(s[i]);
					if(l<minV) minV = l;
					if(l>maxV) maxV = l;
				}
			if(maxV<minV) //all cells are empty/NaN
				return DataType.FLOAT;
			if(minV>=Byte.MIN_VALUE && maxV<=Byte.MAX_VALUE)
				return DataType.BYTE;
			if(minV>=Short.MIN_VALUE && maxV<=Short.MAX_VALUE)
				return DataType.SHORT;
			if(minV>=Integer.MIN_VALUE && maxV<=Integer.MAX_VALUE)
				return DataType.INT;
			return DataType.LONG;
		} catch(NumberFormatException nfe) {
			System.err.println("  got NumberFormatException: "+nfe.getLocalizedMessage());
		};
		return DataType.DOUBLE;
	}
	private void retrieveVarTable() {
		wb.getSheets().forEach( (Sheet sheet) -> {
			orderedSheets.add(sheet);
			System.out.println(String.format("[DEBUG] Sheet %d: %s (%s)",
					sheet.getIndex(), sheet.getName(), sheet.getVisibility().toString()));
			try {
				int colmin=Integer.MAX_VALUE, colmax=Integer.MIN_VALUE, nrows=-1;
				int rowmin=Integer.MAX_VALUE, rowmax=Integer.MAX_VALUE, ri=-1;
				List<Row> rows = sheet.read();
				for(Row r: rows) {
					ri++;
					Optional<Cell> f = r.getFirstNonEmptyCell();
					if(!f.isPresent()) continue;
					if(ri<rowmin) rowmin = ri;
					rowmax = ri;
					Optional<Cell> s = getLastNonEmptyCell(r);
					int f_cidx = f.get().getColumnIndex();
					int s_cidx = s.get().getColumnIndex();
					if(f_cidx<colmin) colmin = f_cidx;
					if(s_cidx>colmax) colmax = s_cidx;
					nrows++;
				}
				ExcelVariable varSheet = new ExcelVariable(sheet.getName(), sheet, true, 0, DataType.STRING, rowmin,
						"dimXlen"+(colmax+1-colmin),"dimYlen"+(rowmax+1-rowmin));
				var_table.put(sheet, new ExcelVariable[colmax+2-colmin]);
				for(int vi=1; vi<colmax+2-colmin; vi++)
					var_table.get(sheet)[vi] = new ExcelVariable();
				all_vars.add(varSheet);
				var_table.get(sheet)[0] = varSheet;
				System.out.println(String.format("[DEBUG]     %d rows read...", nrows));
				System.out.println(String.format("[DEBUG]     %d columns read...", colmax+1-colmin));
				List<Optional<Cell>> column = new ArrayList<>();
				for(int vi=colmin; vi<=colmax; vi++) {
					column.clear();
					Optional<Cell> head = Optional.ofNullable(rows.get(rowmin).getCell(vi));
					String varname = head.isPresent() ? head.get().getText() : String.format("ExcelVar$%03d",vi);
					boolean var_exists = false;
					for(Sheet s: var_table.keySet())
						for(ExcelVariable v: var_table.get(s))
							if (v!=null && v.name != null)
								var_exists |= v.name.equals(varname);
					if(var_exists) {
						//TODO: same variable in other sheets?
						//sheet structure not applicable for reading to a DataFrame?
						System.err.println("Warning: dublicate of variable <"+varname+">\n"+
						                   "         old entries will be overwritten");
					}
					for(ri=rowmin+1; ri<=rowmax; ri++)
						column.add(Optional.ofNullable(rows.get(ri).getCell(vi)));
					boolean allNumber=true,allBool=true;
					for(Optional<Cell> oc: column)
						if(oc.isPresent())
							switch(oc.get().getType()) {
								case EMPTY: break;
								case BOOLEAN: allNumber = false; break;
								case NUMBER: allBool = false; break;
								default: allNumber = false; allBool = false; break;
							}
					DataType dt = allNumber?(allBool?DataType.FLOAT:getNumberType(column)):(allBool?DataType.BOOL:DataType.STRING);
					ExcelVariable eVar =new ExcelVariable(varname, sheet, false, vi, dt, rowmin+1, "dimYlen"+(rowmax-rowmin));
					all_vars.add(eVar);
					var_table.get(sheet)[vi+1-colmin] = eVar;
				}
			} catch(Exception e) {
				e.printStackTrace();
			}
		});
	}
	
	private void describeSheet(Sheet sheet,String pad1,String pad2) { //NetcdfReader::describeGroup
		System.out.println(pad1+wL+sheet.getName()+":");
//		System.out.println(pad2+"  "+wT+wL+"Dimensions:");
		List<ExcelVariable> vars = new ArrayList<>();
		for(ExcelVariable v: var_table.get(sheet))
			if(v.name.length()>0)
				vars.add(v);
		describeVariables(vars, pad2+"  ");
	}
//	private void describeDimensions(List<Dimension> dims, String pad) {
//		System.out.println(pad+wT+wL+"Dimensions:");
//		for(int i=0; i<dims.size(); i++) {
//			String p = pad+wP+"  "+(i+1==dims.size()?wE:wT)+wL;
//			System.out.println(p+dims.get(i).getName()+" ("+dims.get(i).getLength()+")");
//		}
//	}
	private void describeVariables(List<ExcelVariable> vars, String pad) {
//		System.out.println(pad+wT+wL+"Variables:");
		for(int i=0; i<vars.size(); i++) {
			ExcelVariable var = vars.get(i);
			String p = pad+(i+1==vars.size()?wE:wT)+wL;
			String t = String.format("%-6s", var.type.name().toLowerCase());
			String n = var.name+"(";
			for(int d=0; d<var.dims.length; d++)
				n += (d>0?",":"")+var.dims[d].name;
			n += ")";
			System.out.println(p+t+" "+n);
		}
	}
//	private void describeAttributes(List<Attribute> attribs, String pad1, String pad2) {
//		boolean global = pad1.equals(wE);
//		System.out.println(pad1+wL+(global?"Gloabl a":"A")+"ttributes");
//		for(int a=0; a<attribs.size(); a++) {
//			Attribute att = attribs.get(a);
//			System.out.println(pad2+" "+(a+1==attribs.size()?wE:wT)+wL+att.getName()+": "+att.getStringValue());
//		}
//	}

	
	private class ExcelVariable {
		public String name;
		public JKDim[] dims;
		public Sheet sheet;
		public boolean isWholeSheet;
		public int columnID;
		public DataType type;
		public int title_row;
		
		public ExcelVariable() {
			name = ""; sheet = null; isWholeSheet = false;
			columnID = -1; type = DataType.STRING; title_row = -1;
			dims = null;
		}
		public ExcelVariable(String name, Sheet sheet, boolean isWholeSheet, int columnID, DataType type, int title_row, String... dims) {
			this.name = name; this.sheet = sheet; this.isWholeSheet = isWholeSheet;
			this.columnID = columnID; this.type = type; this.title_row = title_row;
			this.dims = new JKDim[dims.length];
			for(int i=0; i<dims.length; i++)
				this.dims[i] = new JKDim(dims[i], Integer.parseInt(dims[i].substring(7)));
		}
		public String toString() {
			return "ExcelVariable@"+Integer.toHexString(this.hashCode())+"[name=\""+name+"\",sheet=\""+(sheet==null?"null":sheet.getName())+"\",type="+type.name().toLowerCase()+"]";
		}
	}
}
