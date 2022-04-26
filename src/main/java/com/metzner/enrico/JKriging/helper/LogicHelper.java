package com.metzner.enrico.JKriging.helper;

import java.util.ArrayList;
import java.util.List;

import com.metzner.enrico.JKriging.data.Constants;
import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.DataFrame.DataType;
import com.metzner.enrico.JKriging.data.DataFrame2D;
import com.metzner.enrico.JKriging.data.DataFrame3D;

public class LogicHelper {
	
	private final static String[] mathfncts = {"abs","acos","asin","atan","atan2","cbrt","cos","exp","log","nan","pow","sin","sqrt","tan"};
	//abbreviation to mathfunctions:           A     C      S      T      U       Q      c     e     L     N     p     s     r      t
	
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
			return new String[][] { {"k", "B", "false", "false"} };
		} else {
			optimizeConditions(sub_conds);
		}
		System.out.println("[COND-LISP] convert conditions to lisp code"); //TODO remove
		String[][] conds = new String[sub_conds.size()][1];
		for(int s=0; s<conds.length; s++)
			conds[s] = sub_conds.get(s).split(" ");
		return conds;
	}
	private static void optimizeConditions(List<String> cond_list) {
		boolean foundOptimization = true;
		while(foundOptimization) {
			foundOptimization = false;
			for(int r=cond_list.size()-1; r>=0; r--) {
				String[] tokens = cond_list.get(r).split(" ");
				if(tokens[2].charAt(0)=='r') {
					int r2 = Integer.parseInt(tokens[2].substring(1));
					String[] parts = cond_list.get(r2).split(" ");
					if(parts[0].charAt(0)=='k') {
						replaceConstantInHigherConditions(r2, parts[2], cond_list);
						foundOptimization = true;
						break;
					}
				}
				if(tokens[3].charAt(0)=='r') {
					int r3 = Integer.parseInt(tokens[3].substring(1));
					String[] parts = cond_list.get(r3).split(" ");
					if(parts[0].charAt(0)=='k') {
						replaceConstantInHigherConditions(r3, parts[2], cond_list);
						foundOptimization = true;
						break;
					}
				}
			}
		}
	}
	private static void replaceConstantInHigherConditions(int r, String newContent, List<String> lowlua) {
		for(int i=lowlua.size()-1; i>=0; i--) {
			String[] tokens = lowlua.get(i).split(" ");
			boolean replacedSomething = false;
			if(tokens[2].charAt(0)=='r') {
				int r2 = Integer.parseInt(tokens[2].substring(1));
				if(r2==r) { tokens[2] = newContent; replacedSomething = true; } else
				if(r2>r) { tokens[2] = "r"+(r2-1); replacedSomething = true; }
			}
			if(tokens[3].charAt(0)=='r') {
				int r3 = Integer.parseInt(tokens[3].substring(1));
				if(r3==r) { tokens[3] = newContent; replacedSomething = true; } else
				if(r3>r) { tokens[3] = "r"+(r3-1); replacedSomething = true; }
			}
			if(replacedSomething)
				lowlua.set(i, tokens[0]+" "+tokens[1]+" "+tokens[2]+" "+tokens[3]);
		}
		lowlua.remove(r);
	}
	private static int buildConditionTree(String[] names, DataType[] types, String _cond, List<String> cond_list, String err_msg) {
		System.out.println("[COND-TREE] analyse expression: \""+_cond+"\"..."); //TODO remove
		System.out.println("[COND-TREE] search for brackets..."); //TODO remove
		String temp = _cond.trim();
		int err_offset = _cond.indexOf(temp);
		
		//*********************************
		//resolve brackets -> subconditions
		boolean hasBrackets = true;
		while(hasBrackets) {
			hasBrackets = false;
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
			if(bracketCount!=0) {
				err_msg = "malformed expression, no opening/closing bracket found";
				return lastClose+err_offset;
			}
			if(dependendOpen>=0) hasBrackets = true;
			if(hasBrackets) {
				if(lastClose==_cond_c.length-1 && dependendOpen==0) {
					System.out.println("[COND-TREE] brackets are meaning less, parse condition without brackets..."); //TODO remove
				}
				System.out.println("[COND-TREE] found bracket in: \""+temp+"\"");
				String cursor = ""; for(int i=temp.length()-1; i>=0; i--) cursor += " ";
				cursor = FormatHelper.setCharAt('1', dependendOpen, cursor);
				cursor = FormatHelper.setCharAt('2', lastClose,     cursor);
				System.out.println("[COND-TREE]                   \""+cursor+"\"");
				String sub_condition = temp.substring(dependendOpen+1, lastClose);
				//System.out.println("[COND-TREE]     analyse condition within brackets: \""+sub_condition+"\""); //TODO remove
				int err = buildConditionTree(names, types, sub_condition, cond_list, err_msg);
				if(err>=0)
					return err + temp.indexOf(sub_condition) + err_offset;
				//System.out.println("[COND-TREE] temporaly replace brackets with ###"); //TODO remove
				String temp2 = "[r"+(cond_list.size()-1)+"]";
				if(dependendOpen>0) temp2 = temp.substring(0, dependendOpen) + temp2;
				if(lastClose<_cond_c.length-1) temp2 += temp.substring(lastClose+1);
				temp = temp2;
				//System.out.println("[COND-TREE]     therefor analyse \""+temp+"\""); //TODO remove
			}
		}
		
		//evaluate type of condition and go to sub conditions if necessary
		System.out.println("[COND-TREE] check: temp=\""+temp+"\""); //TODO remove
		
		//top priority are logical expression, so check for them first
		int an = temp.lastIndexOf("&&");
		int or = temp.lastIndexOf("||");
		int maxp_l = -1; char ct_l = ' ';
		if(an>maxp_l) { maxp_l=an; ct_l = '&'; } if(or>maxp_l) { maxp_l=or; ct_l = '|'; }
		if(maxp_l>=0) {
			int err = buildConditionTree(names, types, temp.substring(0,maxp_l), cond_list, err_msg);
			if(err>=0) return err+err_offset;
			String leftR = "r"+(cond_list.size()-1);
			err = buildConditionTree(names, types, temp.substring(maxp_l+2), cond_list, err_msg);
			if(err>=0) return err+err_offset;
			String rightR = "r"+(cond_list.size()-1);
			cond_list.add((maxp_l==an?"&":"|")+" B "+leftR+" "+rightR);
			return -1;
		}
		
		//second priority are conditional expressions, so check them next
		int le = temp.lastIndexOf("<=");
		int ge = temp.lastIndexOf(">=");
		int eq = temp.lastIndexOf("==");
		int lt = temp.lastIndexOf("<");
		int gt = temp.lastIndexOf(">");
		int maxp_c = -1; char ct_c = ' ';
		if(eq>maxp_c) { maxp_c=eq; ct_c = '='; }
		if(lt>maxp_c) { maxp_c=lt; ct_c = '<';
			if(le>=maxp_c) { maxp_c=le; ct_c = 'l'; } }
		if(gt>maxp_c) { maxp_c=gt; ct_c = '>';
			if(ge>=maxp_c) { maxp_c=ge; ct_c = 'g'; } }
		if(maxp_c>=0) {
			int err = buildConditionTree(names, types, temp.substring(0, maxp_c), cond_list, err_msg);
			if(err>=0) return err+err_offset;
			String leftR = "r"+(cond_list.size()-1);
			int maxp_c2 = maxp_c + (ct_c=='<'||ct_c=='>'?1:2);
			err = buildConditionTree(names, types, temp.substring(maxp_c2), cond_list, err_msg);
			if(err>=0) return err+maxp_c2+err_offset;
			String rightR = "r"+(cond_list.size()-1);
			cond_list.add(ct_c+" D "+leftR+" "+rightR);
			return -1;
		}
		
		//next check for boolean inversion operator "~"
		int no = temp.lastIndexOf("~");
		if(no>0) {
			err_msg = "\"~\" can only be befor brackets of boolean variables!";
			return no+err_offset;
		} else
		if(no==0) {
			int err = buildConditionTree(names, types, temp.substring(1), cond_list, err_msg);
			if(err>=0) return err+1+err_offset;
			cond_list.add("n B r"+(cond_list.size()-1)+" false");
			return -1;
		}
		
		//lastly check mathematical expression
		int pp = temp.lastIndexOf("+");
		int nn = temp.lastIndexOf("-");
		int mm = temp.lastIndexOf("*");
		int dd = temp.lastIndexOf("/");
		int pw = temp.lastIndexOf("^");
		int maxp_m = -1; char ct_m = ' ';
		if(pp>maxp_m) { maxp_m=pp; ct_m='+'; }
		if(nn>maxp_m) { maxp_m=nn; ct_m='-'; }
		if(maxp_m<0) {
			if(mm>maxp_m) { maxp_m=mm; ct_m='*'; }
			if(dd>maxp_m) { maxp_m=dd; ct_m='/'; }
		}
		if(maxp_m<0) {
			if(pw>maxp_m) { maxp_m=pw; ct_m='^'; }
		}
		if(maxp_m>=0) {
			int err = buildConditionTree(names, types, temp.substring(0, maxp_m), cond_list, err_msg);
			if(err>=0) return err+err_offset;
			String leftR = "r"+(cond_list.size()-1);
			err = buildConditionTree(names, types, temp.substring(maxp_m+1), cond_list, err_msg);
			if(err>=0) return err+maxp_m+1+err_offset;
			String rightR = "r"+(cond_list.size()-1);
			int numL = Integer.parseInt(leftR.substring(1));
			char dtL = cond_list.get(numL).split(" ")[1].charAt(0);
			int numR = Integer.parseInt(rightR.substring(1));
			char dtR = cond_list.get(numR).split(" ")[1].charAt(0);
			char dt_m = dtL=='L'&&dtR=='L' ? 'L' : 'D';
			cond_list.add(ct_m+" "+dt_m+" "+leftR+" "+rightR);
			return -1;
		}

		//additionally check for mathematical functions
		boolean foundFormula = false;
		char cf = '?'; int flen = 0;
		//"abs","acos","asin","atan","atan2","cbrt","cos","exp","log","nan","pow","sin","sqrt","tan"
		//A     C      S      T      U       Q      c     e     L     N     p     s     r      t
		if(temp.startsWith("abs["))   { foundFormula = true; cf = 'A'; flen = 3; }
		if(temp.startsWith("acos["))  { foundFormula = true; cf = 'C'; flen = 4; }
		if(temp.startsWith("asin["))  { foundFormula = true; cf = 'S'; flen = 4; }
		if(temp.startsWith("atan["))  { foundFormula = true; cf = 'T'; flen = 4; }
		if(temp.startsWith("atan2[")) { foundFormula = true; cf = 'U'; flen = 5; }
		if(temp.startsWith("cbrt["))  { foundFormula = true; cf = 'Q'; flen = 4; }
		if(temp.startsWith("cos["))   { foundFormula = true; cf = 'c'; flen = 3; }
		if(temp.startsWith("exp["))   { foundFormula = true; cf = 'e'; flen = 3; }
		if(temp.startsWith("log["))   { foundFormula = true; cf = 'L'; flen = 3; }
		if(temp.startsWith("nan["))   { foundFormula = true; cf = 'N'; flen = 3; }
		if(temp.startsWith("pow["))   { foundFormula = true; cf = 'p'; flen = 3; }
		if(temp.startsWith("sin["))   { foundFormula = true; cf = 's'; flen = 3; }
		if(temp.startsWith("sqrt["))  { foundFormula = true; cf = 'r'; flen = 4; }
		if(temp.startsWith("tan["))   { foundFormula = true; cf = 't'; flen = 3; }
		if(foundFormula) {
			int err = buildConditionTree(names, types, temp.substring(flen), cond_list, err_msg);
			if(err>=0) return err+flen+2+err_offset;
			String leftR = "r"+(cond_list.size()-1);
			int numL = Integer.parseInt(leftR.substring(1));
			char dtL = cond_list.get(numL).split(" ")[1].charAt(0);
			char dt_m = 'D';
			if(cf=='A' || cf=='p')
				dt_m = dtL=='L' ? 'L' : 'D';
			cond_list.add(cf+" "+dt_m+" "+leftR+" c0");
			return -1;
		}
		//TODO implement mathematical functions
		
		
		//this block can only be reached, if the expression contains a variable, a constant or the remnant of a bracket constellation
		if(temp.equals("true")) { cond_list.add("k B true true"); return -1; }
		if(temp.equals("false")) { cond_list.add("k B false false"); return -1; }
		int v_id = -2;
		switch(temp.charAt(0)) {
			case '@': v_id = DataHelper.strings_index(names, temp.substring(1));
				if(v_id<0) {
					err_msg = "variable @\""+temp.substring(1)+"\"not found";
					return err_offset;
				}
				char ct_a = types[v_id]==DataType.BOOL ? 'B' :
					types[v_id]==DataType.FLOAT || types[v_id]==DataType.DOUBLE ? 'D' : 'L';
				v_id += Constants.FIRST_IDX;
				cond_list.add("k "+ct_a+" $"+v_id+" $"+v_id);
				break;
			case '$': try {
						v_id = Integer.parseInt(temp.substring(1))-Constants.FIRST_IDX;
						char ct_d = types[v_id]==DataType.BOOL ? 'B' :
							types[v_id]==DataType.FLOAT || types[v_id]==DataType.DOUBLE ? 'D' : 'L';
						cond_list.add("k "+ct_d+" $"+v_id+" $"+v_id);
					} catch(NumberFormatException nfe) {
						//nfe.printStackTrace();
						err_msg = "cannot interprete variable id $\""+temp.substring(1)+"\"";
						return err_offset;
					} break;
			case 'r':
			case '[':
				String s_id = temp.charAt(0)=='r' ? temp.substring(1) : temp.substring(2, temp.length()-1);
				try {
					int sc_id_r = Integer.parseInt(s_id);
					String[] token_r = cond_list.get(sc_id_r).split(" ");
					char ct_r = '?';
					switch(token_r[0].charAt(0)) {
						case 'l': case 'g': case '<': case '=': case '>':
							ct_r = 'B'; break;
						case 'A': case 'k':
							ct_r = token_r[1].charAt(0); break;
						default:
							ct_r = 'D'; break;
					}
					cond_list.add("k "+ct_r+" r"+sc_id_r+" r"+sc_id_r);
				} catch(NumberFormatException nfe) {
					//nfe.printStackTrace();
					err_msg = "Parse error for getting sub-condition id";
					return err_offset;
				} break;
			case 'f': case 'F':
				cond_list.add("k B false false"); break;
			case 't': case 'T':
				cond_list.add("k B true true"); break;
			default:
				//interprete content as number:
				try {
					if(temp.contains(".")) {
						cond_list.add("k D c"+Double.parseDouble(temp)+" c0.0");
					} else {
						cond_list.add("k L c"+Long.parseLong(temp)+" c0");
					}
				} catch(NumberFormatException nfe) {
					nfe.printStackTrace();
					err_msg = "cannot interprete number \""+temp+"\"";
					return err_offset;
				}
		}
		
//		//analyse expression by type
//		if(maxp<0) { // no extra logic needed
//			int v_id = -2;
//			if(temp.charAt(0)=='@') {
//				v_id = DataHelper.strings_index(names, temp.substring(1));
//			}
//			if(temp.charAt(0)=='$') {
//				//if mathematical function...
//				int f_id = DataHelper.strings_index(mathfncts, temp.substring(1));
//				//TODO implement math...
//				//else it is numerical constant
//				if(f_id<0) try{
//					v_id = Integer.parseInt(temp.substring(1))-Constants.FIRST_IDX;
//				}catch(NumberFormatException nfe) {
//					err_msg = "could not find variable"; System.err.println(err_msg);
//					return 1+err_offset;
//				}
//			}
//			if(v_id==-1) {
//				err_msg = "variable not found"; System.err.println(err_msg);
//			}
//			if(v_id>=0) {
//				if(types[v_id]==DataType.BOOL) {
//					cond_list.add("k B $"+v_id+" t");
//				} else if(types[v_id]==DataType.BYTE || types[v_id]==DataType.SHORT || types[v_id]==DataType.INT || types[v_id]==DataType.LONG) {
//					cond_list.add("= L $"+v_id+" 0");
//					cond_list.add("n B r"+(cond_list.size()-1)+" false");
//				}
//			}
//			return -1;
//		} else if(ct=='n') {
//			if(maxp!=0) {
//				err_msg = "bad formated condition string, expect ~ as first char in sub condition"; System.err.println(err_msg);
//				return 0+err_offset;
//			}
//			//temp = temp.substring(1).trim();
//			String negSide = temp.substring(1).trim();
//			System.out.println("[COND-TREE]     analyse \""+negSide+"\" for negated expression..."); //TODO remove
//			int idx_logic_neg = Math.max(negSide.indexOf("<"), negSide.indexOf(">"));
//			idx_logic_neg = Math.max(idx_logic_neg, negSide.indexOf("="));
//			if(idx_logic_neg>=0) {
//				if(ct=='&' || ct=='|') {
//					err_msg = "malformed expression, only conditional expressions (<,<=,==,>=,>) or expressions within brackets can be negated"; System.err.println(err_msg);
//					return idx_logic_neg+temp.indexOf(negSide)+err_offset;
//				}
//				int negErr = buildConditionTree(names, types, negSide, cond_list, err_msg);
//				if(negErr!=-1) return 1+negErr+temp.indexOf(negSide)+err_offset;
//				cond_list.add("n B r"+(cond_list.size()-1)+" false");
//			} else {
//				if(negSide.equals("true")) cond_list.add("k B false false");
//				if(negSide.equals("false")) cond_list.add("k B true true");
//				int v_id = -2;
//				if(negSide.charAt(0)=='@') {
//					v_id = DataHelper.strings_index(names, negSide.substring(1));
//				}
//				if(negSide.charAt(0)=='$') try{
//					v_id = Integer.parseInt(negSide.substring(1))-Constants.FIRST_IDX;
//				}catch(NumberFormatException nfe) {
//					err_msg = "Could not identify variable"; System.err.println(err_msg);
//					return 1+err_offset;
//				}
//				if(v_id==-1) {
//					err_msg = "variable not found"; System.err.println(err_msg);
//					return 1+err_offset;
//				}
//				if(v_id>=0) {
//					if(types[v_id]==DataType.BOOL) {
//						cond_list.add("n B $"+v_id+" false");
//					} else if(types[v_id]==DataType.BYTE || types[v_id]==DataType.SHORT || types[v_id]==DataType.INT || types[v_id]==DataType.LONG) {
//						cond_list.add("= L $"+v_id+" 0");
//					}
//				}
//				if(hasBrackets && negSide.charAt(0)=='#') {
//					int replErr = buildConditionTree(names, types, replStr, cond_list, err_msg);
//					if(replErr!=-1) return replErr+temp.indexOf("###")+err_offset;
//					cond_list.add("n B r"+(cond_list.size()-1)+" false");
//				}
//			}
//		} else {
//			int leftEnd = maxp, rightBegin = maxp+1;
//			if(!(ct=='<' || ct=='>' || ct=='n')) rightBegin++;
//			String leftSide = temp.substring(0, leftEnd).trim();
//			if(hasBrackets && replPosB+3<=leftEnd) leftSide = leftSide.replace("###", replStr);
//			String rightSide = temp.substring(rightBegin).trim();
//			if(hasBrackets && replPosB>=rightBegin) rightSide = rightSide.replace("###", replStr);
//			System.out.println("[COND-TREE] splitted expression in left=\""+leftSide+"\" and right=\""+rightSide+"\""); //TODO remove
//			//analyse left side
//			int idx_logic_left = Math.max(leftSide.indexOf("<"), leftSide.indexOf(">"));
//			idx_logic_left = Math.max(idx_logic_left, Math.max(leftSide.indexOf("="), leftSide.indexOf("###")));
//			idx_logic_left = Math.max(idx_logic_left, Math.max(leftSide.indexOf("&"), leftSide.indexOf("|")));
//			idx_logic_left = Math.max(idx_logic_left, leftSide.indexOf("~"));
//			String leftCond = ""; DataType leftType = DataType.BOOL;
//			if(idx_logic_left>=0) {
//				if(!(ct=='&' || ct=='|')) {
//					err_msg = "malformed expression, expect <=> only in leafs of condition tree"; System.err.println(err_msg);
//					return leftEnd+err_offset;
//				}
//				int leftErr = buildConditionTree(names, types, leftSide, cond_list, err_msg);
//				if(leftErr!=-1) return leftErr+temp.indexOf(leftSide)+err_offset;
//				leftCond = "r"+(cond_list.size()-1);
//			} else {
//				if(leftSide.equals("true")) { leftCond = "true"; leftType = DataType.BOOL; }
//				if(leftSide.equals("false")) { leftCond = "false"; leftType = DataType.BOOL; }
//				int v_id = -2;
//				if(leftSide.charAt(0)=='@') {
//					v_id = DataHelper.strings_index(names, leftSide.substring(1));
//				}
//				if(leftSide.charAt(0)=='$') try{
//					v_id = Integer.parseInt(leftSide.substring(1))-Constants.FIRST_IDX;
//				}catch(NumberFormatException nfe) {
//					err_msg = "could not find variable"; System.err.println(err_msg);
//					return temp.indexOf(leftSide)+err_offset;
//				}
//				if(v_id==-1) {
//					err_msg = "variable not found"; System.err.println(err_msg);
//					return temp.indexOf(leftSide)+err_offset;
//				}
//				if(v_id>=0) {
//					leftType = types[v_id];
//					leftCond = "$"+v_id;
//				} else if(leftCond.length()<1) {
//					int idx_num_real = Math.max(leftSide.indexOf('.'), Math.max(leftSide.indexOf('e'), leftSide.indexOf('E')));
//					boolean read_successful = false;
//					if(idx_num_real<0) {
//						try {
//							leftCond = "c"+Long.parseLong(leftSide);
//							leftType = DataType.LONG;
//							read_successful = true;
//						}catch(NumberFormatException nfe) {
//							//do nothing
//						}
//					}
//					if(!read_successful) {
//						try {
//							leftCond = "c"+Double.parseDouble(leftSide);
//							leftType = DataType.DOUBLE;
//						}catch(NumberFormatException nfe) {
//							err_msg = "could not find variable"; System.err.println(err_msg);
//							return temp.indexOf(leftSide)+err_offset;
//						}
//					}
//				}
//			}
//			//analyse right side
//			int idx_logic_right = Math.max(rightSide.indexOf("<"), rightSide.indexOf(">"));
//			idx_logic_right = Math.max(idx_logic_right, Math.max(rightSide.indexOf("="), rightSide.indexOf("###")));
//			idx_logic_right = Math.max(idx_logic_right, Math.max(rightSide.indexOf("&"), rightSide.indexOf("|")));
//			idx_logic_right = Math.max(idx_logic_right, rightSide.indexOf("~"));
//			String rightCond = ""; DataType rightType = DataType.BOOL;
//			if(idx_logic_right>=0) {
//				if(!(ct=='&' || ct=='|')) {
//					err_msg = "malformed expression, expect <=> only in leafs of condition tree"; System.err.println(err_msg);
//					return leftEnd+err_offset;
//				}
//				int rightErr = buildConditionTree(names, types, rightSide, cond_list, err_msg);
//				if(rightErr!=-1) return rightErr+temp.indexOf(rightSide)+err_offset;
//				rightCond = "r"+(cond_list.size()-1);
//			} else {
//				if(rightSide.equals("true")) { rightCond = "true"; rightType = DataType.BOOL; }
//				if(rightSide.equals("false")) { rightCond = "false"; rightType = DataType.BOOL; }
//				int v_id = -2;
//				if(rightSide.charAt(0)=='@') {
//					v_id = DataHelper.strings_index(names, rightSide.substring(1));
//				}
//				if(rightSide.charAt(0)=='$') try{
//					v_id = Integer.parseInt(rightSide.substring(1))+Constants.FIRST_IDX;
//				}catch(NumberFormatException nfe) {
//					return temp.indexOf(rightSide)+err_offset;
//				}
//				if(v_id==-1) {
//					err_msg = "(numerical) variable not found"; System.err.println(err_msg);
//					return temp.indexOf(rightSide)+err_offset;
//				}
//				if(v_id>=0) {
//					rightType = types[v_id];
//					rightCond = "$"+v_id;
//				} else if(rightCond.length()<1) {
//					int idx_num_real = Math.max(rightSide.indexOf('.'), Math.max(rightSide.indexOf('e'), rightSide.indexOf('E')));
//					boolean read_successful = false;
//					if(idx_num_real<0) {
//						try {
//							rightCond = "c"+Long.parseLong(rightSide);
//							rightType = DataType.LONG;
//							read_successful = true;
//						}catch(NumberFormatException nfe) {
//							//do nothing
//						}
//					}
//					if(!read_successful) {
//						try {
//							rightCond = "c"+Double.parseDouble(rightSide);
//							rightType = DataType.DOUBLE;
//						}catch(NumberFormatException nfe) {
//							err_msg = "could not find (numerical) variable"; System.err.println(err_msg);
//							return temp.indexOf(rightSide)+err_offset;
//						}
//					}
//				}
//			}
//			String condVarType = "";
//			if(ct=='&' || ct=='|') {
//				if(leftType!=DataType.BOOL || rightType!=DataType.BOOL)
//					return leftEnd+err_offset;
//				condVarType = "B";
//			} else {
//				condVarType = "L";
//				if(leftType==DataType.FLOAT || leftType==DataType.DOUBLE) condVarType = "D";
//				if(rightType==DataType.FLOAT || rightType==DataType.DOUBLE) condVarType = "D";
//			}
//			cond_list.add(ct+" "+condVarType+" "+leftCond+" "+rightCond);
//		}
		return -1;
	}
	
	
	public static boolean evaluateConditionLines(DataFrame _df, int _i, String[][] conds_lowLua) {
		return evaluateConditionLinesInGeneral(_df, _i, _i, _i, conds_lowLua); }
	public static boolean evaluateConditionLines(DataFrame2D _df, int _i, int _j, String[][] conds_lowLua) {
		return evaluateConditionLinesInGeneral(_df, _i, _j, _j, conds_lowLua); }
	public static boolean evaluateConditionLines(DataFrame3D _df, int _i, int _j, int _k, String[][] conds_lowLua) {
		return evaluateConditionLinesInGeneral(_df, _i, _j, _k, conds_lowLua); }
	private static boolean evaluateConditionLinesInGeneral(Object dataframe, int _i, int _j, int _k, String[][] conds_lowLua) {
		double[] res = new double[conds_lowLua.length];
		boolean leftB=false, rightB=false;
		long leftL=0L, rightL=0L;
		double leftD=0d, rightD=0d;
		for(int s=0; s<res.length; s++) {
			char c0 = conds_lowLua[s][0].charAt(0);
			char c1 = conds_lowLua[s][1].charAt(0);
			char c2 = conds_lowLua[s][2].charAt(0);
			char c3 = conds_lowLua[s][3].charAt(0);
			switch(c1) {
				case 'B':
					if(c2=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
						leftB = getBoolean(dataframe, v_id, _i, _j, _k);
					} else
					if(c2=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][2].substring(1));
						leftB = res[sc_id] > 0.5d;
					} else { // conds_lowLua[s][2].charAt(0)=='c'
						leftB = Boolean.parseBoolean(conds_lowLua[s][2]);
					}
					if(c3=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
						rightB = getBoolean(dataframe, v_id, _i, _j, _k);
					} else
					if(c3=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][3].substring(1));
						rightB = res[sc_id] > 0.5d;
					} else { // conds_lowLua[s][3].charAt(0)=='c'
						rightB = Boolean.parseBoolean(conds_lowLua[s][3]);
					}
					break;
				case 'D':
					if(c2=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
						leftD = getDouble(dataframe, v_id, _i, _j, _k);
					} else
					if(c2=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][2].substring(1));
						leftD = res[sc_id];
					} else { // conds_lowLua[s][2].charAt(0)=='c'
						leftD = Double.parseDouble(conds_lowLua[s][2].substring(1));
					}
					if(c3=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
						rightD = getDouble(dataframe, v_id, _i, _j, _k);
					} else
					if(c3=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][3].substring(1));
						rightD = res[sc_id];
					} else { // conds_lowLua[s][3].charAt(0)=='c'
						rightD = Double.parseDouble(conds_lowLua[s][3].substring(1));
					}
					break;
				case 'L':
					if(c2=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][2].substring(1))+Constants.FIRST_IDX;
						leftL = getLong(dataframe, v_id, _i, _j, _k);
					} else
					if(c2=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][2].substring(1));
						leftL = (long) res[sc_id];
					} else { // conds_lowLua[s][2].charAt(0)=='c'
						leftL = Long.parseLong(conds_lowLua[s][2].substring(1));
					}
					if(c3=='$') {
						int v_id = Integer.parseInt(conds_lowLua[s][3].substring(1))+Constants.FIRST_IDX;
						rightL = getLong(dataframe, v_id, _i, _j, _k);
					} else
					if(c3=='r') {
						int sc_id = Integer.parseInt(conds_lowLua[s][3].substring(1));
						rightL = (long) res[sc_id];
					} else { // conds_lowLua[s][3].charAt(0)=='c'
						rightL = Long.parseLong(conds_lowLua[s][3].substring(1));
					}
					break;
				default:
					break;
			}
			switch(c0) {
				case '+': res[s] = c1=='L' ? leftL+rightL : leftD+rightD; break;
				case '-': res[s] = c1=='L' ? leftL-rightL : leftD-rightD; break;
				case '*': res[s] = c1=='L' ? leftL*rightL : leftD*rightD; break;
				case '/': res[s] = c1=='L' ? leftL/rightL : leftD/rightD; break;
				//{"abs","acos","asin","atan","atan2","cbrt","cos","exp","log","pow","sin","sqrt","tan"};
				case 'A': res[s] = c1=='L' ? Math.abs(leftL) : Math.abs(leftD); break;
				case 'C': res[s] = c1=='L' ? Math.acos(leftL) : Math.acos(leftD); break;
				case 'S': res[s] = c1=='L' ? Math.asin(leftL) : Math.asin(leftD); break;
				case 'T': res[s] = c1=='L' ? Math.atan(leftL) : Math.atan(leftD); break;
				case 'U': res[s] = c1=='L' ? Math.atan2(leftL, rightL) : Math.atan2(leftD, rightD); break;
				case 'Q': res[s] = c1=='L' ? Math.cbrt(leftL) : Math.cbrt(leftD); break;
				case 'c': res[s] = c1=='L' ? Math.cos(leftL) : Math.cos(leftD); break;
				case 'e': res[s] = c1=='L' ? Math.exp(leftL) : Math.exp(leftD); break;
				case 'L': res[s] = c1=='L' ? Math.log(leftL) : Math.log(leftD); break;
				case 'p': res[s] = c1=='L' ? Math.pow(leftL, rightL) : Math.pow(leftD, rightD); break;
				case 's': res[s] = c1=='L' ? Math.sin(leftL) : Math.sin(leftD); break;
				case 'r': res[s] = c1=='L' ? Math.sqrt(leftL) : Math.sqrt(leftD); break;
				case 't': res[s] = c1=='L' ? Math.tan(leftL) : Math.tan(leftD); break;
				case '^': res[s] = c1=='L' ? Math.pow(leftL, rightL) : Math.pow(leftD, rightD); break;
				case '<': switch(c1) {
						case 'D': res[s] = leftD<rightD ? 1d : 0d; break;
						case 'L': res[s] = leftL<rightL ? 1d : 0d; break;
						default:  res[s] = 0d;
					} break;
				case '>': switch(c1) {
						case 'D': res[s] = leftD>rightD ? 1d : 0d; break;
						case 'L': res[s] = leftL>rightL ? 1d : 0d; break;
						default:  res[s] = 0d;
					} break;
				case 'l': switch(c1) {
						case 'D': res[s] = leftD<=rightD ? 1d : 0d; break;
						case 'L': res[s] = leftL<=rightL ? 1d : 0d; break;
						default:  res[s] = 0d;
					} break;
				case 'g': switch(c1) {
						case 'D': res[s] = leftD>=rightD ? 1d : 0d; break;
						case 'L': res[s] = leftL>=rightL ? 1d : 0d; break;
						default:  res[s] = 0d;
					} break;
				case '=': switch(c1) {
						case 'B': res[s] = leftB==rightB ? 1d : 0d; break;
						case 'D': res[s] = leftD==rightD ? 1d : 0d; break;
						case 'L': res[s] = leftL==rightL ? 1d : 0d; break;
						default:  res[s] = 0d;
					} break;
				case '&': res[s] = leftB && rightB ? 1d : 0d; break;
				case '|': res[s] = leftB || rightB ? 1d : 0d; break;
				case 'n': res[s] = leftB ? 0d : 1d; break;
				case 'k': switch(c1) {
					case 'B': res[s] = leftB ? 1d : 0d; break;
					case 'D': res[s] = leftD; break;
					case 'L': res[s] = leftL; break;
				} break;
				default: res[s] = 0d; break;
			}
		}
		return res[res.length-1] > 0.5d;
	}
	
	private static boolean getBoolean(Object dataframe, int v_id, int i, int j, int k) {
		if(dataframe instanceof DataFrame) {
			if( ((DataFrame)dataframe).getVariableType(v_id) != DataType.BOOL)
				return false;
			return ((boolean[]) ((DataFrame)dataframe).getArray(v_id))[i];
		}
		if(dataframe instanceof DataFrame2D) {
			if( ((DataFrame2D)dataframe).getVariableType(v_id) != DataType.BOOL)
				return false;
			return ((boolean[][]) ((DataFrame2D)dataframe).getArray(v_id))[j][i];
		}
		if(dataframe instanceof DataFrame3D) {
			if( ((DataFrame3D)dataframe).getVariableType(v_id) != DataType.BOOL)
				return false;
			return ((boolean[][][]) ((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
		}
		return false;
	}
	private static long getLong(Object dataframe, int v_id, int i, int j, int k) {
		DataType dt = null; int dims = 0;
		if(dataframe instanceof DataFrame)   { dt = ((DataFrame)dataframe).getVariableType(v_id);   dims = 1; }
		if(dataframe instanceof DataFrame2D) { dt = ((DataFrame2D)dataframe).getVariableType(v_id); dims = 2; }
		if(dataframe instanceof DataFrame3D) { dt = ((DataFrame3D)dataframe).getVariableType(v_id); dims = 3; }
		if(dims==0) return 0L;
		if(dt==null) return 0L;
		switch(dt) {
			case BYTE: return
					dims==1 ? ((byte[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((byte[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((byte[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case SHORT: return
					dims==1 ? ((short[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((short[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((short[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case INT: return
					dims==1 ? ((int[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((int[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((int[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case LONG: return
					dims==1 ? ((long[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((long[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((long[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			default: return 0L;
		}
	}
	private static double getDouble(Object dataframe, int v_id, int i, int j, int k) {
		DataType dt = null; int dims = 0;
		if(dataframe instanceof DataFrame)   { dt = ((DataFrame)dataframe).getVariableType(v_id);   dims = 1; }
		if(dataframe instanceof DataFrame2D) { dt = ((DataFrame2D)dataframe).getVariableType(v_id); dims = 2; }
		if(dataframe instanceof DataFrame3D) { dt = ((DataFrame3D)dataframe).getVariableType(v_id); dims = 3; }
		if(dims==0) return 0d;
		switch(dt) {
			case BYTE: return
					dims==1 ? ((byte[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((byte[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((byte[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case SHORT: return
					dims==1 ? ((short[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((short[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((short[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case INT: return
					dims==1 ? ((int[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((int[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((int[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case LONG: return
					dims==1 ? ((long[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((long[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((long[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case FLOAT: return
					dims==1 ? ((float[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((float[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((float[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			case DOUBLE: return
					dims==1 ? ((double[])((DataFrame)dataframe).getArray(v_id))[i] :
					dims==2 ? ((double[][])((DataFrame2D)dataframe).getArray(v_id))[j][i] :
						((double[][][])((DataFrame3D)dataframe).getArray(v_id))[k][j][i];
			default: return 0d;
		}
	}
}
