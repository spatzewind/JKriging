package com.metzner.enrico.JKriging.data.reader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
	private Map<Sheet,ExcelVariable[]> var_table;
	
	public ExcelReader(File file) throws IOException {
		super(file);
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
		return null;
	}
	public DataFrame2D getVars2D(String... vars) throws IllegalAccessException {
		return null;
	}
	public DataFrame3D getVars3D(String... vars) throws IllegalAccessException {
		return null;
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
			System.out.println(String.format("Sheet %d: %s (%s)",
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
				ExcelVariable varSheet = new ExcelVariable(sheet.getName(), sheet, true, -1, DataType.STRING, rows.size(),
						"dimXlen"+(colmax+1-colmin));
				var_table.put(sheet, new ExcelVariable[colmax+2-colmin]);
				all_vars.add(varSheet);
				var_table.get(sheet)[0] = varSheet;
				System.out.println(String.format("    %d rows read...", nrows));
				System.out.println(String.format("    %d columns read...", colmax+1-colmin));
				List<Optional<Cell>> column = new ArrayList<>();
				for(int vi=colmin; vi<=colmax; vi++) {
					column.clear();
					Optional<Cell> head = Optional.ofNullable(rows.get(rowmin).getCell(vi));
					String varname = head.isPresent() ? head.get().getText() : String.format("ExcelVar$%03d",vi);
					boolean var_exists = false;
					for(ExcelVariable v: var_table)
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
					
					var_table.add(new ExcelVariable(varname, sheet, false, rowmin, dt, rowmax-rowmin));
				}
			} catch(Exception e) {
				e.printStackTrace();
			}
		});
	}
	
	private void describeSheet(Sheet sheet,String pad1,String pad2) {
		System.out.println(pad1+wL+sheet.getName()+":");
//		System.out.println(pad2+"  "+wT+wL+"Dimensions:");
		List<String> vars = new ArrayList<>();
		for(String v: var_table.keySet())
			if(var_table.get(v)[0].equals(sheet.getName()))
				vars.add(v);
		for(int i=0; i<vars.size(); i++) {
			String var = vars.get(i);
			String[] info = var_table.get(var);
			String p = pad2+"   "+(i+1==vars.size()?wE:wT)+wL;
			String t = String.format("%-6s", info[4].toLowerCase());
			String n = var+" ("+info[3].substring(4)+")";
			System.out.println(p+t+" "+n);
		}
	}
//	private void describeDimensions(List<Dimension> dims, String pad) {
//		System.out.println(pad+wT+wL+"Dimensions:");
//		for(int i=0; i<dims.size(); i++) {
//			String p = pad+wP+"  "+(i+1==dims.size()?wE:wT)+wL;
//			System.out.println(p+dims.get(i).getName()+" ("+dims.get(i).getLength()+")");
//		}
//	}
//	private void describeVariables(List<String> vars, String pad) {
////		System.out.println(pad+wT+wL+"Variables:");
//		for(int i=0; i<vars.size(); i++) {
//			String var = vars.get(i);
//			String[] info = var_table.get(var);
//			String p = pad+(i+1==vars.size()?wE:wT)+wL;
//			String t = String.format("%-6s", info[4].toLowerCase());
//			String n = var+" ("+info[3].substring(4)+")";
//			System.out.println(p+t+" "+n);
//		}
//	}
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
		public String[] dims;
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
			this.dims = dims;
		}
	}
}
