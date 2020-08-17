package com.metzner.enrico.JKriging.helper;

import java.util.ArrayList;
import java.util.List;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame.DataType;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;

public class LogicHelper {
	
	public static String[][] readConditions(DataFrame _df, String _condition) {
		return conditions2lowLisp(_df.allVariableNames(), _df.allVariableTypes(), _condition);
	}
	public static String[][] readConditions(DataFrame2D _df, String _condition) {
		return conditions2lowLisp(_df.allVariableNames(), _df.allVariableTypes(), _condition);
	}
	public static String[][] readConditions(DataFrame3D _df, String _condition) {
		return conditions2lowLisp(_df.allVariableNames(), _df.allVariableTypes(), _condition);
	}
	private static String[][] conditions2lowLisp(String[] var_names, DataType[] var_types, String condition) {
		List<String> sub_conds = new ArrayList<String>();
		String error_message = "";
		System.out.println("\n\n[COND-LISP] build condition tree..."); //TODO remove
		int err = buildConditionTree(var_names, var_types, condition, sub_conds, error_message);
		if(err!=-1) {
			System.err.println("An error occured while parsing the condition string:");
			System.err.println("    "+error_message);
			System.err.println(condition);
			String marker = "";
			for(int m=0; m<err; m++) marker += " ";
			System.err.println(marker+"^");
			return new String[][] { {"s", "B", "false", "false"} };
		}
		System.out.println("[COND-LISP] convert conditions to lisp code"); //TODO remove
		String[][] conds = new String[sub_conds.size()][1];
		for(int s=0; s<conds.length; s++)
			conds[s] = sub_conds.get(s).split(" ");
		return conds;
	}
	private static int buildConditionTree(String[] names, DataType[] types, String _cond, List<String> cond_list, String err_msg) {
		System.out.println("[COND-TREE] analyse expression: \""+_cond+"\"..."); //TODO remove
		System.out.println("[COND-TREE] search for brackets..."); //TODO remove
		String temp = _cond.trim();
		int lastClose = -1, dependendOpen = -1, bracketCount=0;
		char[] _cond_c = temp.toCharArray();
		for(int c=_cond_c.length-1; c>=0; c--) {
			if(_cond_c[c]==')') {
				bracketCount++;
				if(lastClose==-1) lastClose = c;
			}
			if(_cond_c[c]=='(') {
				bracketCount--;
				if(bracketCount==0) {
					dependendOpen = c;
					break;
				}
			}
		}
		int err_offset = _cond.indexOf(temp);
		if(bracketCount!=0) {
			err_msg = "malformed expression, no opening/closing bracket found";
			return lastClose+err_offset;
		}
		boolean hasBrackets = false;
		if(dependendOpen>=0) hasBrackets = true;
		int replPosB = -1;
		String replStr = "";
		if(hasBrackets) {
			if(lastClose==_cond_c.length-1 && temp.substring(0, dependendOpen).trim().length()==0) {
				System.out.println("[COND-TREE] brackets are meaning less, parse condition without brackets..."); //TODO remove
				int err = buildConditionTree(names, types, temp.substring(dependendOpen+1, lastClose), cond_list, err_msg);
				if(err>=0) err += err_offset;
				return err;
			}
			System.out.println("[COND-TREE] temporaly replace brackets with ###"); //TODO remove
			replPosB = dependendOpen;
			replStr = temp.substring(dependendOpen, lastClose+1);
			if(lastClose==_cond_c.length-1) {
				temp = temp.substring(0, dependendOpen)+"###";
			} else {
				String temp2 = temp.substring(0, dependendOpen)+"###"+temp.substring(lastClose+1);
				temp = ""+temp2;
			}
			System.out.println("[COND-TREE]     therefor analyse \""+temp+"\""); //TODO remove
		}
		//evaluate type of condition and go to sub conditions if necessary
		System.out.println("[COND-TREE] check: temp=\""+temp+"\""); //TODO remove
		//conditional expressions
		int le = temp.lastIndexOf("<=");
		int ge = temp.lastIndexOf(">=");
		int eq = temp.lastIndexOf("==");
		int lt = temp.lastIndexOf("<");
		int gt = temp.lastIndexOf(">");
		int maxp_c = -1; char ct_c = ' ';
		if(eq>maxp_c) { maxp_c=eq; ct_c = '='; }
		if(lt>maxp_c) { maxp_c=lt; ct_c = '<'; } if(gt>maxp_c) { maxp_c=gt; ct_c = '>'; }
		if(le+1>maxp_c) { maxp_c=le; ct_c = 'l'; } if(ge+1>maxp_c) { maxp_c=ge; ct_c = 'g'; }
		//logical expressions
		int an = temp.lastIndexOf("&&");
		int or = temp.lastIndexOf("||");
		int no = temp.lastIndexOf("~");
		int maxp_l = -1; char ct_l = ' ';
		if(an>maxp_l) { maxp_l=an; ct_l = '&'; } if(or>maxp_l) { maxp_l=or; ct_l = '|'; }
		int maxp = -1; char ct = ' ';
		//combined by priority
		if(maxp_c>=(ct_c=='<'||ct_c=='>'?0:1)) { maxp = maxp_c; ct = ct_c; }
		if(no>=0) { maxp=no; ct = 'n'; } // only use negation, if a variable or an expression enclosed within brackets
		if(maxp_l>=(ct_l=='~'?0:1)) { maxp = maxp_l; ct = ct_l; } // logical expression (&&,||) has higher priority than conditional expression (<,<=,==,>=,>)
		System.out.println("[COND-TREE] found "+(ct=='&'||ct=='|'||ct=='n' ? "logical" : "conditional")+" expression: \""+ct+"\""); //TODO remove
		//analyse expression by type
		if(maxp<0) { // no extra logic needed
			if(temp.equals("true")) { cond_list.add("s B true true"); return -1; }
			if(temp.equals("false")) { cond_list.add("s B false false"); return -1; }
			int v_id = -2;
			if(temp.charAt(0)=='@') {
				v_id = DataHelper.strings_index(names, temp.substring(1));
			}
			if(temp.charAt(0)=='$') try{
				v_id = Integer.parseInt(temp.substring(1))-Constants.FIRST_IDX;
			}catch(NumberFormatException nfe) {
				err_msg = "could not find variable"; System.err.println(err_msg);
				return 1+err_offset;
			}
			if(v_id==-1) {
				err_msg = "variable not found"; System.err.println(err_msg);
			}
			if(v_id>=0) {
				if(types[v_id]==DataType.BOOL) {
					cond_list.add("s B $"+v_id+" t");
				} else if(types[v_id]==DataType.BYTE || types[v_id]==DataType.SHORT || types[v_id]==DataType.INT || types[v_id]==DataType.LONG) {
					cond_list.add("= L $"+v_id+" 0");
					cond_list.add("n B r"+(cond_list.size()-1)+" false");
				}
			}
			return -1;
		} else if(ct=='n') {
			if(maxp!=0) {
				err_msg = "bad formated condition string, expect ~ as first char in sub condition"; System.err.println(err_msg);
				return 0+err_offset;
			}
			//temp = temp.substring(1).trim();
			String negSide = temp.substring(1).trim();
			System.out.println("[COND-TREE]     analyse \""+negSide+"\" for negated expression..."); //TODO remove
			int idx_logic_neg = Math.max(negSide.indexOf("<"), negSide.indexOf(">"));
			idx_logic_neg = Math.max(idx_logic_neg, negSide.indexOf("="));
			if(idx_logic_neg>=0) {
				if(ct=='&' || ct=='|') {
					err_msg = "malformed expression, only conditional expressions (<,<=,==,>=,>) or expressions within brackets can be negated"; System.err.println(err_msg);
					return idx_logic_neg+temp.indexOf(negSide)+err_offset;
				}
				int negErr = buildConditionTree(names, types, negSide, cond_list, err_msg);
				if(negErr!=-1) return 1+negErr+temp.indexOf(negSide)+err_offset;
				cond_list.add("n B r"+(cond_list.size()-1)+" false");
			} else {
				if(negSide.equals("true")) cond_list.add("s B false false");
				if(negSide.equals("false")) cond_list.add("s B true true");
				int v_id = -2;
				if(negSide.charAt(0)=='@') {
					v_id = DataHelper.strings_index(names, negSide.substring(1));
				}
				if(negSide.charAt(0)=='$') try{
					v_id = Integer.parseInt(negSide.substring(1))-Constants.FIRST_IDX;
				}catch(NumberFormatException nfe) {
					err_msg = "Could not identify variable"; System.err.println(err_msg);
					return 1+err_offset;
				}
				if(v_id==-1) {
					err_msg = "variable not found"; System.err.println(err_msg);
					return 1+err_offset;
				}
				if(v_id>=0) {
					if(types[v_id]==DataType.BOOL) {
						cond_list.add("n B $"+v_id+" false");
					} else if(types[v_id]==DataType.BYTE || types[v_id]==DataType.SHORT || types[v_id]==DataType.INT || types[v_id]==DataType.LONG) {
						cond_list.add("= L $"+v_id+" 0");
					}
				}
				if(hasBrackets && negSide.charAt(0)=='#') {
					int replErr = buildConditionTree(names, types, replStr, cond_list, err_msg);
					if(replErr!=-1) return replErr+temp.indexOf("###")+err_offset;
					cond_list.add("n B r"+(cond_list.size()-1)+" false");
				}
			}
		} else {
			int leftEnd = maxp, rightBegin = maxp+1;
			if(!(ct=='<' || ct=='>' || ct=='n')) rightBegin++;
			String leftSide = temp.substring(0, leftEnd).trim();
			if(hasBrackets && replPosB+3<=leftEnd) leftSide = leftSide.replace("###", replStr);
			String rightSide = temp.substring(rightBegin).trim();
			if(hasBrackets && replPosB>=rightBegin) rightSide = rightSide.replace("###", replStr);
			System.out.println("[COND-TREE] splitted expression in left=\""+leftSide+"\" and right=\""+rightSide+"\""); //TODO remove
			//analyse left side
			int idx_logic_left = Math.max(leftSide.indexOf("<"), leftSide.indexOf(">"));
			idx_logic_left = Math.max(idx_logic_left, Math.max(leftSide.indexOf("="), leftSide.indexOf("###")));
			idx_logic_left = Math.max(idx_logic_left, Math.max(leftSide.indexOf("&"), leftSide.indexOf("|")));
			idx_logic_left = Math.max(idx_logic_left, leftSide.indexOf("~"));
			String leftCond = ""; DataType leftType = DataType.BOOL;
			if(idx_logic_left>=0) {
				if(!(ct=='&' || ct=='|')) {
					err_msg = "malformed expression, expect <=> only in leafs of condition tree"; System.err.println(err_msg);
					return leftEnd+err_offset;
				}
				int leftErr = buildConditionTree(names, types, leftSide, cond_list, err_msg);
				if(leftErr!=-1) return leftErr+temp.indexOf(leftSide)+err_offset;
				leftCond = "r"+(cond_list.size()-1);
			} else {
				if(leftSide.equals("true")) { leftCond = "true"; leftType = DataType.BOOL; }
				if(leftSide.equals("false")) { leftCond = "false"; leftType = DataType.BOOL; }
				int v_id = -2;
				if(leftSide.charAt(0)=='@') {
					v_id = DataHelper.strings_index(names, leftSide.substring(1));
				}
				if(leftSide.charAt(0)=='$') try{
					v_id = Integer.parseInt(leftSide.substring(1))-Constants.FIRST_IDX;
				}catch(NumberFormatException nfe) {
					err_msg = "could not find variable"; System.err.println(err_msg);
					return temp.indexOf(leftSide)+err_offset;
				}
				if(v_id==-1) {
					err_msg = "variable not found"; System.err.println(err_msg);
					return temp.indexOf(leftSide)+err_offset;
				}
				if(v_id>=0) {
					leftType = types[v_id];
					leftCond = "$"+v_id;
				} else if(leftCond.length()<1) {
					int idx_num_real = Math.max(leftSide.indexOf('.'), Math.max(leftSide.indexOf('e'), leftSide.indexOf('E')));
					boolean read_successful = false;
					if(idx_num_real<0) {
						try {
							leftCond = "c"+Long.parseLong(leftSide);
							leftType = DataType.LONG;
							read_successful = true;
						}catch(NumberFormatException nfe) {
							//do nothing
						}
					}
					if(!read_successful) {
						try {
							leftCond = "c"+Double.parseDouble(leftSide);
							leftType = DataType.DOUBLE;
						}catch(NumberFormatException nfe) {
							err_msg = "could not find variable"; System.err.println(err_msg);
							return temp.indexOf(leftSide)+err_offset;
						}
					}
				}
			}
			//analyse right side
			int idx_logic_right = Math.max(rightSide.indexOf("<"), rightSide.indexOf(">"));
			idx_logic_right = Math.max(idx_logic_right, Math.max(rightSide.indexOf("="), rightSide.indexOf("###")));
			idx_logic_right = Math.max(idx_logic_right, Math.max(rightSide.indexOf("&"), rightSide.indexOf("|")));
			idx_logic_right = Math.max(idx_logic_right, rightSide.indexOf("~"));
			String rightCond = ""; DataType rightType = DataType.BOOL;
			if(idx_logic_right>=0) {
				if(!(ct=='&' || ct=='|')) {
					err_msg = "malformed expression, expect <=> only in leafs of condition tree"; System.err.println(err_msg);
					return leftEnd+err_offset;
				}
				int rightErr = buildConditionTree(names, types, rightSide, cond_list, err_msg);
				if(rightErr!=-1) return rightErr+temp.indexOf(rightSide)+err_offset;
				rightCond = "r"+(cond_list.size()-1);
			} else {
				if(rightSide.equals("true")) { rightCond = "true"; rightType = DataType.BOOL; }
				if(rightSide.equals("false")) { rightCond = "false"; rightType = DataType.BOOL; }
				int v_id = -2;
				if(rightSide.charAt(0)=='@') {
					v_id = DataHelper.strings_index(names, rightSide.substring(1));
				}
				if(rightSide.charAt(0)=='$') try{
					v_id = Integer.parseInt(rightSide.substring(1))+Constants.FIRST_IDX;
				}catch(NumberFormatException nfe) {
					return temp.indexOf(rightSide)+err_offset;
				}
				if(v_id==-1) {
					err_msg = "(numerical) variable not found"; System.err.println(err_msg);
					return temp.indexOf(rightSide)+err_offset;
				}
				if(v_id>=0) {
					rightType = types[v_id];
					rightCond = "$"+v_id;
				} else if(rightCond.length()<1) {
					int idx_num_real = Math.max(rightSide.indexOf('.'), Math.max(rightSide.indexOf('e'), rightSide.indexOf('E')));
					boolean read_successful = false;
					if(idx_num_real<0) {
						try {
							rightCond = "c"+Long.parseLong(rightSide);
							rightType = DataType.LONG;
							read_successful = true;
						}catch(NumberFormatException nfe) {
							//do nothing
						}
					}
					if(!read_successful) {
						try {
							rightCond = "c"+Double.parseDouble(rightSide);
							rightType = DataType.DOUBLE;
						}catch(NumberFormatException nfe) {
							err_msg = "could not find (numerical) variable"; System.err.println(err_msg);
							return temp.indexOf(rightSide)+err_offset;
						}
					}
				}
			}
			String condVarType = "";
			if(ct=='&' || ct=='|') {
				if(leftType!=DataType.BOOL || rightType!=DataType.BOOL)
					return leftEnd+err_offset;
				condVarType = "B";
			} else {
				condVarType = "L";
				if(leftType==DataType.FLOAT || leftType==DataType.DOUBLE) condVarType = "D";
				if(rightType==DataType.FLOAT || rightType==DataType.DOUBLE) condVarType = "D";
			}
			cond_list.add(ct+" "+condVarType+" "+leftCond+" "+rightCond);
		}
		return -1;
	}
	

	public static boolean evaluateConditionLines(DataFrame _df, int _i, String[][] conds_lowLua) {
		boolean[] res = new boolean[conds_lowLua.length];
		boolean leftB=false, rightB=false;
		long leftL=0L, rightL=0L;
		double leftD=0d, rightD=0d;
		for(int s=0; s<res.length; s++) {
			if("&|ns".contains(""+conds_lowLua[s][0].charAt(0))) {
				if(conds_lowLua[s][2].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						leftB = ((boolean[])_df.getArray(v_id))[_i];
					} else {
						leftB = false;
					}
				}
				if(conds_lowLua[s][2].charAt(0)=='r') leftB = res[Integer.parseInt(conds_lowLua[s][2].substring(1))];
				if(conds_lowLua[s][2].charAt(0)=='f') leftB = false;
				if(conds_lowLua[s][2].charAt(0)=='t') leftB = true;
				if(conds_lowLua[s][3].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						rightB = ((boolean[])_df.getArray(v_id))[_i];
					} else {
						rightB = false;
					}
				}
				if(conds_lowLua[s][3].charAt(0)=='r') rightB = res[Integer.parseInt(conds_lowLua[s][3].substring(1))];
				if(conds_lowLua[s][3].charAt(0)=='f') rightB = false;
				if(conds_lowLua[s][3].charAt(0)=='t') rightB = true;
			}
			if("<>lg=".contains(""+conds_lowLua[s][0].charAt(0))) {
				switch(conds_lowLua[s][1].charAt(0)) {
					case 'D':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftD = ((byte[])_df.getArray(v_id))[_i]; break;
								case SHORT:  leftD = ((short[])_df.getArray(v_id))[_i]; break;
								case INT:    leftD = ((int[])_df.getArray(v_id))[_i]; break;
								case LONG:   leftD = ((long[])_df.getArray(v_id))[_i]; break;
								case FLOAT:  leftD = ((float[])_df.getArray(v_id))[_i]; break;
								case DOUBLE: leftD = ((double[])_df.getArray(v_id))[_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							if(conds_lowLua[s][2].indexOf('.')<0 && conds_lowLua[s][2].indexOf('e')<0 && conds_lowLua[s][2].indexOf('E')<0) {
								leftD = (double) Long.parseLong(conds_lowLua[s][2].substring(1));
							} else {
								leftD = Double.parseDouble(conds_lowLua[s][2].substring(1));
							}
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightD = ((byte[])_df.getArray(v_id))[_i]; break;
								case SHORT:  rightD = ((short[])_df.getArray(v_id))[_i]; break;
								case INT:    rightD = ((int[])_df.getArray(v_id))[_i]; break;
								case LONG:   rightD = ((long[])_df.getArray(v_id))[_i]; break;
								case FLOAT:  rightD = ((float[])_df.getArray(v_id))[_i]; break;
								case DOUBLE: rightD = ((double[])_df.getArray(v_id))[_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							if(conds_lowLua[s][3].indexOf('.')<0 && conds_lowLua[s][3].indexOf('e')<0 && conds_lowLua[s][3].indexOf('E')<0) {
								rightD = (double) Long.parseLong(conds_lowLua[s][3].substring(1));
							} else {
								rightD = Double.parseDouble(conds_lowLua[s][3].substring(1));
							}
						}
						break;
					case 'L':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftL = ((byte[])_df.getArray(v_id))[_i]; break;
								case SHORT:  leftL = ((short[])_df.getArray(v_id))[_i]; break;
								case INT:    leftL = ((int[])_df.getArray(v_id))[_i]; break;
								case LONG:   leftL = ((long[])_df.getArray(v_id))[_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							leftL = Long.parseLong(conds_lowLua[s][2].substring(1));
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightL = ((byte[])_df.getArray(v_id))[_i]; break;
								case SHORT:  rightL = ((short[])_df.getArray(v_id))[_i]; break;
								case INT:    rightL = ((int[])_df.getArray(v_id))[_i]; break;
								case LONG:   rightL = ((long[])_df.getArray(v_id))[_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							rightL = Long.parseLong(conds_lowLua[s][3].substring(1));
						}
						break;
					default:
						break;
				}
			}
			switch(conds_lowLua[s][0].charAt(0)) {
				case '<': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<rightD); break;
						case 'L': res[s] = (leftL<rightL); break;
						default: res[s] = false;
					} break;
				case '>': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>rightD); break;
						case 'L': res[s] = (leftL>rightL); break;
						default: res[s] = false;
					} break;
				case 'l': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<=rightD); break;
						case 'L': res[s] = (leftL<=rightL); break;
						default: res[s] = false;
					} break;
				case 'g': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>=rightD); break;
						case 'L': res[s] = (leftL>=rightL); break;
						default: res[s] = false;
					} break;
				case '=': switch(conds_lowLua[s][1].charAt(0)) {
						case 'B': res[s] = (leftB==rightB); break;
						case 'D': res[s] = (leftD==rightD); break;
						case 'L': res[s] = (leftL==rightL); break;
						default: res[s] = false;
					} break;
				case '&': res[s] = (leftB && rightB); break;
				case '|': res[s] = (leftB || rightB); break;
				case 'n': res[s] = !leftB; break;
				case 's': res[s] = (conds_lowLua[s][2].charAt(0)=='t'); break;
				default: res[s] = false; break;
			}
		}
		return res[res.length-1];
	}
	public static boolean evaluateConditionLines(DataFrame2D _df, int _i, int _j, String[][] conds_lowLua) {
		boolean[] res = new boolean[conds_lowLua.length];
		boolean leftB=false, rightB=false;
		long leftL=0L, rightL=0L;
		double leftD=0d, rightD=0d;
		for(int s=0; s<res.length; s++) {
			if("&|ns".contains(""+conds_lowLua[s][0].charAt(0))) {
				if(conds_lowLua[s][2].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						leftB = ((boolean[][])_df.getArray(v_id))[_j][_i];
					} else {
						leftB = false;
					}
				} else { // conds_lowLua[s][2].charAt(0)=='r' with r='res'
					leftB = res[Integer.parseInt(conds_lowLua[s][2].substring(1))];
				}
				if(conds_lowLua[s][2].charAt(0)=='r') leftB = res[Integer.parseInt(conds_lowLua[s][2].substring(1))];
				if(conds_lowLua[s][2].charAt(0)=='f') leftB = false;
				if(conds_lowLua[s][2].charAt(0)=='t') leftB = true;
				if(conds_lowLua[s][3].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						rightB = ((boolean[][])_df.getArray(v_id))[_j][_i];
					} else {
						rightB = false;
					}
				} else { // conds_lowLua[s][3].charAt(0)=='r' with r='res'
					rightB = res[Integer.parseInt(conds_lowLua[s][3].substring(1))];
				}
				if(conds_lowLua[s][3].charAt(0)=='r') rightB = res[Integer.parseInt(conds_lowLua[s][3].substring(1))];
				if(conds_lowLua[s][3].charAt(0)=='f') rightB = false;
				if(conds_lowLua[s][3].charAt(0)=='t') rightB = true;
			}
			if("<>lg=".contains(""+conds_lowLua[s][0].charAt(0))) {
				switch(conds_lowLua[s][1].charAt(0)) {
					case 'D':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftD = ((byte[][])_df.getArray(v_id))[_j][_i]; break;
								case SHORT:  leftD = ((short[][])_df.getArray(v_id))[_j][_i]; break;
								case INT:    leftD = ((int[][])_df.getArray(v_id))[_j][_i]; break;
								case LONG:   leftD = ((long[][])_df.getArray(v_id))[_j][_i]; break;
								case FLOAT:  leftD = ((float[][])_df.getArray(v_id))[_j][_i]; break;
								case DOUBLE: leftD = ((double[][])_df.getArray(v_id))[_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							leftD = Double.parseDouble(conds_lowLua[s][2].substring(1));
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightD = ((byte[][])_df.getArray(v_id))[_j][_i]; break;
								case SHORT:  rightD = ((short[][])_df.getArray(v_id))[_j][_i]; break;
								case INT:    rightD = ((int[][])_df.getArray(v_id))[_j][_i]; break;
								case LONG:   rightD = ((long[][])_df.getArray(v_id))[_j][_i]; break;
								case FLOAT:  rightD = ((float[][])_df.getArray(v_id))[_j][_i]; break;
								case DOUBLE: rightD = ((double[][])_df.getArray(v_id))[_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							rightD = Double.parseDouble(conds_lowLua[s][3].substring(1));
						}
						break;
					case 'L':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftL = ((byte[][])_df.getArray(v_id))[_j][_i]; break;
								case SHORT:  leftL = ((short[][])_df.getArray(v_id))[_j][_i]; break;
								case INT:    leftL = ((int[][])_df.getArray(v_id))[_j][_i]; break;
								case LONG:   leftL = ((long[][])_df.getArray(v_id))[_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							leftL = Long.parseLong(conds_lowLua[s][2].substring(1));
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightL = ((byte[][])_df.getArray(v_id))[_j][_i]; break;
								case SHORT:  rightL = ((short[][])_df.getArray(v_id))[_j][_i]; break;
								case INT:    rightL = ((int[][])_df.getArray(v_id))[_j][_i]; break;
								case LONG:   rightL = ((long[][])_df.getArray(v_id))[_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							rightL = Long.parseLong(conds_lowLua[s][3].substring(1));
						}
						break;
					default:
						break;
				}
			}
			switch(conds_lowLua[s][0].charAt(0)) {
				case '<': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<rightD); break;
						case 'L': res[s] = (leftL<rightL); break;
						default: res[s] = false;
					} break;
				case '>': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>rightD); break;
						case 'L': res[s] = (leftL>rightL); break;
						default: res[s] = false;
					} break;
				case 'l': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<=rightD); break;
						case 'L': res[s] = (leftL<=rightL); break;
						default: res[s] = false;
					} break;
				case 'g': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>=rightD); break;
						case 'L': res[s] = (leftL>=rightL); break;
						default: res[s] = false;
					} break;
				case '=': switch(conds_lowLua[s][1].charAt(0)) {
						case 'B': res[s] = (leftB==rightB); break;
						case 'D': res[s] = (leftD==rightD); break;
						case 'L': res[s] = (leftL==rightL); break;
						default: res[s] = false;
					} break;
				case '&': res[s] = (leftB && rightB); break;
				case '|': res[s] = (leftB || rightB); break;
				case 'n': res[s] = !leftB; break;
				case 's': res[s] = (conds_lowLua[s][2].charAt(0)=='t'); break;
				default: res[s] = false; break;
			}
		}
		return res[res.length-1];
	}
	public static boolean evaluateConditionLines(DataFrame3D _df, int _i, int _j, int _k, String[][] conds_lowLua) {
		boolean[] res = new boolean[conds_lowLua.length];
		boolean leftB=false, rightB=false;
		long leftL=0L, rightL=0L;
		double leftD=0d, rightD=0d;
		for(int s=0; s<res.length; s++) {
			if("&|ns".contains(""+conds_lowLua[s][0].charAt(0))) {
				if(conds_lowLua[s][2].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						leftB = ((boolean[][][])_df.getArray(v_id))[_k][_j][_i];
					} else {
						leftB = false;
					}
				} else { // conds_lowLua[s][2].charAt(0)=='r' with r='res'
					leftB = res[Integer.parseInt(conds_lowLua[s][2].substring(1))];
				}
				if(conds_lowLua[s][2].charAt(0)=='r') leftB = res[Integer.parseInt(conds_lowLua[s][2].substring(1))];
				if(conds_lowLua[s][2].charAt(0)=='f') leftB = false;
				if(conds_lowLua[s][2].charAt(0)=='t') leftB = true;
				if(conds_lowLua[s][3].charAt(0)=='$') {
					int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
					if(_df.getVariableType(v_id)==DataType.BOOL) {
						rightB = ((boolean[][][])_df.getArray(v_id))[_k][_j][_i];
					} else {
						rightB = false;
					}
				} else { // conds_lowLua[s][3].charAt(0)=='r' with r='res'
					rightB = res[Integer.parseInt(conds_lowLua[s][3].substring(1))];
				}
				if(conds_lowLua[s][3].charAt(0)=='r') rightB = res[Integer.parseInt(conds_lowLua[s][3].substring(1))];
				if(conds_lowLua[s][3].charAt(0)=='f') rightB = false;
				if(conds_lowLua[s][3].charAt(0)=='t') rightB = true;
			}
			if("<>lg=".contains(""+conds_lowLua[s][0].charAt(0))) {
				switch(conds_lowLua[s][1].charAt(0)) {
					case 'D':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftD = ((byte[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case SHORT:  leftD = ((short[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case INT:    leftD = ((int[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case LONG:   leftD = ((long[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case FLOAT:  leftD = ((float[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case DOUBLE: leftD = ((double[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							leftD = Double.parseDouble(conds_lowLua[s][2].substring(1));
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightD = ((byte[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case SHORT:  rightD = ((short[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case INT:    rightD = ((int[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case LONG:   rightD = ((long[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case FLOAT:  rightD = ((float[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case DOUBLE: rightD = ((double[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							rightD = Double.parseDouble(conds_lowLua[s][3].substring(1));
						}
						break;
					case 'L':
						if(conds_lowLua[s][2].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   leftL = ((byte[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case SHORT:  leftL = ((short[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case INT:    leftL = ((int[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case LONG:   leftL = ((long[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][2].charAt(0)=='c'
							leftL = Long.parseLong(conds_lowLua[s][2].substring(1));
						}
						if(conds_lowLua[s][3].charAt(0)=='$') {
							int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
							switch(_df.getVariableType(v_id)) {
								case BYTE:   rightL = ((byte[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case SHORT:  rightL = ((short[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case INT:    rightL = ((int[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								case LONG:   rightL = ((long[][][])_df.getArray(v_id))[_k][_j][_i]; break;
								default: break;
							}
						} else { // conds_lowLua[s][3].charAt(0)=='c'
							rightL = Long.parseLong(conds_lowLua[s][3].substring(1));
						}
						break;
					default:
						break;
				}
			}
			switch(conds_lowLua[s][0].charAt(0)) {
				case '<': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<rightD); break;
						case 'L': res[s] = (leftL<rightL); break;
						default: res[s] = false;
					} break;
				case '>': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>rightD); break;
						case 'L': res[s] = (leftL>rightL); break;
						default: res[s] = false;
					} break;
				case 'l': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD<=rightD); break;
						case 'L': res[s] = (leftL<=rightL); break;
						default: res[s] = false;
					} break;
				case 'g': switch(conds_lowLua[s][1].charAt(0)) {
						case 'D': res[s] = (leftD>=rightD); break;
						case 'L': res[s] = (leftL>=rightL); break;
						default: res[s] = false;
					} break;
				case '=': switch(conds_lowLua[s][1].charAt(0)) {
						case 'B': res[s] = (leftB==rightB); break;
						case 'D': res[s] = (leftD==rightD); break;
						case 'L': res[s] = (leftL==rightL); break;
						default: res[s] = false;
					} break;
				case '&': res[s] = (leftB && rightB); break;
				case '|': res[s] = (leftB || rightB); break;
				case 'n': res[s] = !leftB; break;
				case 's': res[s] = (conds_lowLua[s][2].charAt(0)=='t'); break;
				default: res[s] = false; break;
			}
		}
		return res[res.length-1];
	}
}
