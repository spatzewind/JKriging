package com.metzner.enrico.JKriging.kriging;

import com.metzner.enrico.JKriging.helper.DataHelper;
//import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.MathHelper;

public class KdTree {

	private KdTreeNode root;
	private int dimensions;
	private double progress;
	private double[][] identity_matrix;
	
	public KdTree() {
		root = null;
		dimensions = 0;
		progress = 0d;
		identity_matrix = new double[0][0];
	}
	
	public void build(boolean _polar, boolean _degree, double[]... coords) {
		for(int c=1; c<coords.length; c++) {
			if(coords[0].length != coords[c].length) {
				System.err.println("All dimension have to be the same length!. No tree is built!");
				DataHelper.printStackTrace(System.err);
				return;
			}
		}
		if(_polar) {
			if(coords.length<2) {
				System.err.println("For the use of spherical coordinats at least 2 dimension has to be given, but got "+coords.length+
						". So no tree is built!");
				DataHelper.printStackTrace(System.err);
				return;
			}
			System.out.println("INFO: if spherical coordinates are used, the first 2 dimensions are expected to be longitude and latitude respectivly!");
		}
		dimensions = coords.length;
		int npoints = coords[0].length;
		//(_polar?1:0) make a conversion from 2d spherical coordinates lat and lon to 3d kartesian coordinates possible
		// + 1 is required for indexing, so the index gets ordered parallel to the positions later in build-subroutine
		double[][] internalCoords = new double[dimensions + (_polar ? 1 : 0) + 1][coords[0].length];
		for(int c=0; c<dimensions; c++) {
			int cc = c + (_polar && c>=2 ? 1 : 0);
			for(int p=0; p<npoints; p++) {
				internalCoords[cc][p] = coords[c][p];
			}
		}
		if(_polar) {
			double d2r = Math.PI / 180d;
			for(int p=0; p<npoints; p++) {
				double lon = coords[0][p] * (_degree ? d2r : 1d);
				double lat = coords[1][p] * (_degree ? d2r : 1d);
				internalCoords[0][p] = Math.cos(lat) * Math.cos(lon);
				internalCoords[1][p] = Math.cos(lat) * Math.sin(lon);
				internalCoords[2][p] = Math.sin(lat);
			}
		}
		int indexdim = internalCoords.length-1;
		for(int p=0; p<npoints; p++) internalCoords[indexdim][p] = p + 0.001d;
		
		root = build(internalCoords, 0, npoints-1, 0);

		dimensions = internalCoords.length-1; //index column should not be part of later space transformation
		identity_matrix = new double[dimensions][dimensions];
		for(int j=0; j<dimensions; j++)
			for(int i=0; i<dimensions; i++)
				identity_matrix[j][i] = (i==j ? 1d : 0d);
	}
	private KdTreeNode build(double[][] coords, int index_start, int index_end, int axis) {
		if(index_start!=index_end) {
			int dimensions = coords.length - 1; // coords contain index column, not part of "position"
			double[] temp = new double[coords[axis].length];
			for(int v=0; v<temp.length; v++) temp[v] = coords[axis][v];
			DataHelper.sortem(temp, index_start, index_end, false, coords);
			int i = (index_start+index_end)/2; //<- this is may not be the original index, so the original index is carried in the last column of coords.
			double[] pos = new double[dimensions];
			for(int p=0; p<dimensions; p++) pos[p] = coords[p][i];
			int index = (int) coords[dimensions][i];
			KdTreeNode node = new KdTreeNode(index, pos);
			node.setChildren(index_start<i ? build(coords, index_start, i-1, (axis+1)%dimensions): null,
					         i<index_end ? build(coords, i+1, index_end, (axis+1)%dimensions) : null);
			return node;
		} else {
			double[] pos = new double[dimensions];
			for(int p=0; p<dimensions; p++) pos[p] = coords[p][index_start];
			int index = (int) coords[dimensions][index_start];
			KdTreeNode leaf = new KdTreeNode(index, pos);
			return leaf;
		}
	}
	
	public int getIndexOfClosestPointTo(boolean use_polar_coordinates, boolean is_in_degree, double... coords) {
		return getIndexOfNClosestPointsTo(1, Double.POSITIVE_INFINITY, use_polar_coordinates, is_in_degree, coords)[1];
	}
	public int[] getIndexOfNClosestPointsTo(int number_of_points, boolean use_polar_coordinates, boolean is_in_degree, double... coords) {
		return getIndexOfNClosestPointsTo(number_of_points, Double.POSITIVE_INFINITY, identity_matrix, use_polar_coordinates, is_in_degree, coords);
	}
	public int[] getIndexOfNClosestPointsTo(int number_of_points, double maximum_distance, boolean use_polar_coordinates, boolean is_in_degree, double... coords) {
		return getIndexOfNClosestPointsTo(number_of_points, maximum_distance, identity_matrix, use_polar_coordinates, is_in_degree, coords);
	}
	public int[] getIndexOfNClosestPointsTo(int number_of_points, double[][] rotation_matrix, boolean use_polar_coordinates, boolean is_in_degree, double... coords) {
		return getIndexOfNClosestPointsTo(number_of_points, Double.POSITIVE_INFINITY, rotation_matrix, use_polar_coordinates, is_in_degree, coords);
	}
	public int[] getIndexOfNClosestPointsTo(int number_of_points, double maximum_distance, double[][] rotation_matrix, boolean use_polar_coordinates, boolean is_in_degree, double... coords) {
		//System.out.println("[KDTREE] Search for Indices of N points:\n"+
		//                   "            N = "+number_of_points); //TODO DEBUG remove
		int[] indices = new int[number_of_points+1];
		for(int i=0; i<number_of_points; i++) indices[i+1] = -1;
		indices[0] = 0;
		if(root==null) {
			System.err.println("no k-d-tree was build or an other error occured! Return -1."); DataHelper.printStackTrace(System.err);
			return indices; }
		int rm_len = rotation_matrix.length;
		if(rm_len!=rotation_matrix[0].length) {
			System.err.println("Rotation matrix is not a square matrix. Return [-1,...]."); DataHelper.printStackTrace(System.err);
			return indices; }
		int dimcount = root.getPosition().length;
		if(rm_len!=dimcount) {
			System.err.println("Rotation matrix should be of shape "+dimcount+"x"+dimcount+", but found "+
					rm_len+"x"+rm_len+". Return [-1,...]."); DataHelper.printStackTrace(System.err);
			return indices; }
		//System.out.println("            rotmat = ");
		//FormatHelper.printMat(System.out, rotation_matrix);
		int pc = use_polar_coordinates ? 1 : 0;
		int givenCount = coords.length;
		if(dimcount!= givenCount+pc) {
			System.err.println("too "+(dimcount<givenCount+pc?"much":"few")+" points coordinates are given, "+
					           "expected "+(dimcount-pc)+" but got "+givenCount);
			DataHelper.printStackTrace(System.err);
			return indices;
		}
		double[] pivot = new double[dimcount];
		for(int p=0; p<givenCount; p++) {
			int c = p + (p>=2 ? pc : 0);
			pivot[c] = coords[p];
		}
		if(use_polar_coordinates) {
			double d2r = Math.PI / 180d;
			double lon = coords[0] * (is_in_degree ? d2r : 1d);
			double lat = coords[1] * (is_in_degree ? d2r : 1d);
			pivot[0] = Math.cos(lat) * Math.cos(lon);
			pivot[1] = Math.cos(lat) * Math.sin(lon);
			pivot[2] = Math.sin(lat);
		}
		
		double[][] transposed_inverse_rotation_matrix = MathHelper.transpose(MathHelper.inverse(rotation_matrix));
		//System.out.println("            trpinvrotmat = ");
		//FormatHelper.printMat(System.out, transposed_inverse_rotation_matrix);
		double max_dist_squared = maximum_distance * maximum_distance;
		for(int i=0; i<number_of_points; i++) {
			//System.out.println("            point "+(i+1)+":"); //TODO DEBUG remove
			//System.out.println("                -- min-dist = "+minimum_distance);
			KdTreeNode closest = closestNode(pivot, indices, max_dist_squared, rotation_matrix, transposed_inverse_rotation_matrix, root, 0);
			if(closest==null) break;
			indices[i+1] = closest.getIndex();
			indices[0]++;
		}
		return indices;
	}
	
	
	private KdTreeNode closestNode(double[] _pivot, int[] nodesToExclude, double max_sqr_dist, double[][] rotmat, double[][] trpinvrotmat, KdTreeNode current_branch, int axis) {
		if(current_branch == null) return null;
		
		KdTreeNode next_branch = null;
		KdTreeNode opp_branch = null;
		
		double[] brnc = current_branch.getPosition();
		if(_pivot[axis] < brnc[axis]) {
			next_branch = current_branch.getLeftChild();
			opp_branch  = current_branch.getRightChild();
		} else {
			next_branch = current_branch.getRightChild();
			opp_branch  = current_branch.getLeftChild();
		}
		
		KdTreeNode best = closer_point(
					_pivot,
					closestNode(_pivot, nodesToExclude, max_sqr_dist, rotmat, trpinvrotmat, next_branch, (axis+1)%_pivot.length),
					current_branch,
					nodesToExclude, max_sqr_dist, rotmat);
		double dist = max_sqr_dist;
		if(best!=null) {
			dist = MathHelper.sqdistKD(_pivot, best.getPosition(), rotmat);
		}
		//System.out.println("        inbetween: minimum/best distance: "+dist); //TODO DEBUG remove
		// calculate the distance between the pivot and the dividing plane from the k-d-tree
		//   as dot product of the normal vector of the transformed division plane
		//   and the vector to the division point on that plane
		double[] diff = new double[_pivot.length],
				 norm = new double[_pivot.length];
		for(int p=0; p<_pivot.length; p++) {
			diff[p] = brnc[p] - _pivot[p];
			norm[p] = (p==axis ? 1d : 0d);
		}
		diff = MathHelper.matmul(rotmat, diff);
		norm = MathHelper.matmul(trpinvrotmat, norm);
		double planedist=0d, normmag2=0d;
		for(int p=0; p<_pivot.length; p++) {
			planedist += diff[p] * norm[p];
			normmag2 += norm[p] * norm[p];
		}
		planedist *= planedist / normmag2;
		// if point is closer to the dividing plane as to the nearest neighbour so far,
		//   then the opposite branch could contain a point which is even closer
		if ( planedist < dist )
			best = closer_point(
					_pivot,
					closestNode(_pivot, nodesToExclude, max_sqr_dist, rotmat, trpinvrotmat, opp_branch, (axis+1)%_pivot.length),
					best,
					nodesToExclude, max_sqr_dist, rotmat);
		
		return best;
	}
	private KdTreeNode closer_point(double[] piv, KdTreeNode p1, KdTreeNode p2, int[] excluded_nodes, double maximum_squared_distance, double[][] rotation_matrix) {
		if(p2==null) return p1;
		if(p1==null) return p2;
		double[] p1p = p1.getPosition();
		double[] p2p = p2.getPosition();
		double d1=MathHelper.sqdistKD(piv, p1p, rotation_matrix),
			   d2=MathHelper.sqdistKD(piv, p2p, rotation_matrix);
		boolean p1toEx = false, p2toEx = false;
		for(int i=1; i<excluded_nodes.length; i++) {
			if(p1.getIndex()==excluded_nodes[i]) p1toEx = true;
			if(p2.getIndex()==excluded_nodes[i]) p2toEx = true;
		}
		if(p1toEx || d1>maximum_squared_distance)
			return (p2toEx || d2>maximum_squared_distance ? null : p2);
		if(p2toEx || d2>maximum_squared_distance)
			return p1;
		if(d1<d2) {
			return p1;
		} else {
			return p2;
		}
	}
	
	public KdTreeNode getRoot() { return root; }
	public double getProgress() { return progress; }
	
	
	
	public class KdTreeNode {
		private double[] pos;
		int index;
		KdTreeNode leftChild, rightChild;
		
		public KdTreeNode(int _idx, double... _pos) {
			pos = _pos;
			index = _idx;
			leftChild = null;
			rightChild = null;
		}
		
		public void setIndex(int _idx) { index = _idx; }
		public int getIndex() { return index; }
		public void setPosition(double... _pos) { pos = _pos; }
		public double[] getPosition() { return pos; }
		
		public void setChildren(KdTreeNode _left, KdTreeNode _right) { leftChild = _left; rightChild = _right; }
		public KdTreeNode getLeftChild() { return leftChild; }
		public KdTreeNode getRightChild() { return rightChild; }
		
		public String toJson() {
			String json = "{";
			json += "\"index\": "+index+", ";
			json += "\"point\": ["+pos[0];
			for(int c=1; c<pos.length; c++)
				json += ", "+pos[c];
			json += "]";
			if(leftChild!=null) json += ", \"left\": "+leftChild.toJson();
			if(rightChild!=null) json += ", \"right\": "+rightChild.toJson();
			json += "}";
			return json;
		}
	}
}
