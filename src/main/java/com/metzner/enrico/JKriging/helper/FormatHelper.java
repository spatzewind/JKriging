package com.metzner.enrico.JKriging.helper;

import java.io.PrintStream;

public class FormatHelper {
	
	public static void chknam(String _possible_file_name) {
	}
	
	public static String nf(int _num, int all) {
		return nf(_num,all,' ');
	}
	public static String nf(int _num, int all, char fillchar) {
		String nf = ""+_num;
		while(nf.length()<all) nf = fillchar+nf;
		return nf;
	}
	public static String nf(double _num, int all, int flt) {
		return nf(_num,all,flt,' ');
	}
	public static String nf(double _num, int all, int flt, char fillchar) {
		double _pos = Math.abs(_num);
		long li = (long) _pos;
		long lf = (long) (Math.pow(10, flt)*(_pos-li)+0.5d);
		String nf = ""+lf;
		if(nf.length()>flt) { li++; nf = nf.substring(nf.length()-flt); }
		while(nf.length()<flt) nf = "0"+nf;
		nf = (_num<0d ? "-" : "")+li+"."+nf;
		while(nf.length()<all) nf = fillchar+nf;
		return nf;
	}

	public static String underscore_spaces(String _s) {
		String temp = _s.replaceAll(" ", "_");
		temp = temp.replaceAll("\\\\", "");
		String patterns = "+-*/:.,;#~!?§$%&(){}[]^°";
		String out = "";
		for(char r: temp.toCharArray()) {
			if(!patterns.contains(""+r)) out += r;
		}
		if(out.length()>0)
			if(out.charAt(0)>='0' && out.charAt(0)<='9')
				temp = "v"+out;
		return temp;
	}

	public static String[] splitBySpace(String _txt) {
		String string_wo_doublespace = _txt.trim();
		while(string_wo_doublespace.indexOf("  ")>=0) string_wo_doublespace = string_wo_doublespace.replace("  ", " ");
		return string_wo_doublespace.split(" ");
	}


	public static void printTable(int _trim, int[]... int_columns) {
		int clen = int_columns.length;
		if(clen<=0) return;
		int max_length=0;
		for(int c=0; c<clen; c++) max_length = Math.max(max_length, int_columns[c].length);
		if(_trim>0 && max_length>_trim) max_length = _trim;
		String[][] out = new String[max_length][clen];
		for(int c=0; c<clen; c++) {
			int[] column = int_columns[c];
			for(int l=0; l<max_length; l++) {
				if(l<column.length) {
					out[l][c] = " "+column[l]+" ";
				} else {
					out[l][c] = " ";
				}
			}
		}
		printTable(out);
	}
	public static void printTable(int _trim, double[]... double_columns) {
		int clen = double_columns.length;
		if(clen<=0) return;
		int max_length=0;
		for(int c=0; c<clen; c++) max_length = Math.max(max_length, double_columns[c].length);
		if(_trim>0 && max_length>_trim) max_length = _trim;
		String[][] out = new String[max_length][clen];
		for(int c=0; c<clen; c++) {
			double[] column = double_columns[c];
			for(int l=0; l<max_length; l++) {
				if(l<column.length) {
					out[l][c] = " "+column[l]+" ";
				} else {
					out[l][c] = " ";
				}
			}
		}
		printTable(out);
	}
	
	
	public static void printTable(String[][] strings) {
		boolean[] rowdivs = new boolean[strings.length-1];
		for(int rd=0; rd<rowdivs.length; rd++) rowdivs[rd] = false;
		boolean[] coldivs = new boolean[strings[0].length-1];
		for(int cd=0; cd<coldivs.length; cd++) coldivs[cd] = false;
		printTable(strings, coldivs, rowdivs, false);
	}
	public static void printTable(String[][] strings, boolean[] col_divider, boolean[] row_divider) {
		printTable(strings, col_divider, row_divider, false);
	}
	public static void printTable(String[][] strings, boolean[] col_divider, boolean[] row_divider, boolean transpose) {
		int nrow = strings.length,
			ncol = strings[0].length;
		String[][] txt = (transpose ? new String[ncol][nrow] : strings);
		if(transpose) {
			for(int j=0; j<nrow; j++) for(int i=0; i<ncol; i++) txt[i][j] = strings[j][i];
			nrow = txt.length;
			ncol = txt[0].length;
		}
		String[] hdivs = new String[ncol];
		for(int c=0; c<ncol; c++) {
			int cw = txt[0][c].length();
			for(int r=1; r<nrow; r++) cw = Math.max(cw, txt[r][c].length());
			hdivs[c] = "";
			while(hdivs[c].length()<cw) hdivs[c]+="\u2500";
			while(txt[0][c].length()<cw) txt[0][c]+=" ";
			for(int r=1; r<nrow; r++) while(txt[r][c].length()<cw) txt[r][c] = " "+txt[r][c];
		}
		for(int r=0; r<nrow; r++) {
			String o=txt[r][0];
			for(int c=1; c<ncol; c++) o+=(col_divider[c-1] ? "\u2502" : "")+txt[r][c];
			System.out.println(o);
			if(r+1<nrow) {
				if(row_divider[r]) {
					o = hdivs[0];
					for(int c=1; c<ncol; c++) o+=(col_divider[c-1] ? "\u253c" : "")+hdivs[c];
					System.out.println(o);
				}
			}
		}
	}

	public static void printMat(PrintStream stream, double[][] _mat) {
	    int la = _mat.length, lb = _mat[0].length;
	    String[][] mout = new String[la][lb];
	    //stream.println("    [PRINTMAT] "+la+"x"+lb+" matrix");
	    int ml = 0;
	    for(int a=0; a<la; a++) for(int b=0; b<lb; b++) {
	        mout[a][b] = " "+_mat[a][b]+" ";
	        ml = Math.max(ml,mout[a][b].length());
	    }
	    //stream.println("    [PRINTMAT] "+ml+"chars per entry");
	    for(int a=0; a<la; a++) for(int b=0; b<lb; b++) {
	        while(mout[a][b].length()<ml) { mout[a][b] = " "+mout[a][b]; }
	    }
	    for(int b=0; b<lb; b++) {
	        String l = "";
	        for(int a=0; a<la; a++) l += mout[a][b];
	        stream.println(l);
	    }
	}
}
