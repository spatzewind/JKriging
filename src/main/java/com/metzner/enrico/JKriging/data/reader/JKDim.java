package com.metzner.enrico.JKriging.data.reader;

class JKDim {
	public String name;
	public int len;
	
	public JKDim() {
		name = ""; len = 0;
	}
	public JKDim(String n, int l) {
		name = n; len = l;
	}
	
	@Override
	public boolean equals(Object obj) {
//		System.out.println("[DEBUG] JKDim.equals(...) called");
		if(obj==null) return false;
		if(!(obj instanceof JKDim)) return false;
		if(this.name.equals(((JKDim)obj).name) && this.len==((JKDim)obj).len) return true;
		return super.equals(obj);
	}
	@Override
	public String toString() {
		return "JKDim@"+Integer.toHexString(this.hashCode())+"(name="+this.name+",len="+this.len+")";
	}
}
