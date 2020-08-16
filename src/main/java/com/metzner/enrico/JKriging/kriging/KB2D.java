package com.metzner.enrico.JKriging.kriging;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import com.metzner.enrico.JKriging.helper.FormatHelper;
import com.metzner.enrico.JKriging.helper.LinEquSolver;
import com.metzner.enrico.JKriging.probability.Covariance;

public class KB2D {
	
	private static final int    MAXNST    = 4;
	private static final double UNEST     = Double.NaN;
	private static final double EPSLON    = 1.0e-20d;
	private static final double VERSION   = 3.000d;
	
	
	private int[] it = new int[MAXNST];
	private double[] aa   = new double[MAXNST],
			         cc   = new double[MAXNST],
			         ang  = new double[MAXNST],
			         anis = new double[MAXNST];
	
	private int nx, ny, nd, ktype, nst, idbg;
	private double xmn,xsiz,
		           ymn,ysiz;
	private int ndmin,ndmax, nxdis,nydis;
	private double radius, skmea, c0;
	private double tmin,tmax;
	private String datafl, outfl, dbgfl;
	
	private double[] x,y,vr;
	int MAXSAM, MAXDIS, MAXKD, MAXKRG;
	
	private double[][] rotmat = new double[4][4];
	private double maxcov;
	
	public void readParam(String param_path) {
//        c-----------------------------------------------------------------------
//        c
//        c                  Initialization and Read Parameters
//        c                  **********************************
//        c
//        c The input parameters and data are read in from their files. Some quick
//        c error checking is performed and the statistics of all the variables
//        c being considered are written to standard output.
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    use       msflib
		int ixl,iyl,ivrl;
		//    use       geostat
		double[] var = new double[20];
//        c
//        c Unit numbers:
//        c
		//    lin  = 1
		//    lout = 2
		//    ldbg = 3
//        c
//        c Note VERSION number:
//        c
		System.out.println(" KB2D Version: "+VERSION);
//        c
//        c Get the name of the parameter file - try the default name if no input:
//        c
		if(param_path==null) {
			param_path = "res/kb2d.par";
		} else if(param_path.length()==0) {
			param_path = "res/kb2d.par";
		}
		File f = new File(param_path);
		if(!f.exists()) {
			System.err.println("ERROR - the parameter file does not exist,\n"
					          +"        check for the file and try again\n");
			if(param_path.equals("res/kb2d.par")) {
				System.out.println("        creating a blank parameter file");
				makepar();
			}
			return;
		}
		try(BufferedReader br = new BufferedReader(new FileReader(f))) {
			String line = br.readLine();
//        c
//        c Find Start of Parameters:
//        c
			while(!line.startsWith("STAR")) line = br.readLine();
//        c
//        c Read Input Parameters:
//        c
			line = br.readLine(); String[] parts = FormatHelper.splitBySpace(line);
			datafl = parts[0]; FormatHelper.chknam(datafl); System.out.println(" data file = "+datafl);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ixl = Integer.parseInt(parts[0].trim()); iyl = Integer.parseInt(parts[1].trim()); ivrl = Integer.parseInt(parts[2].trim());
			System.out.println(" columns for X,Y, VR = "+ixl+","+iyl+","+ivrl);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			tmin = Double.parseDouble(parts[0].trim()); tmax = Double.parseDouble(parts[1].trim());
			System.out.println(" trimming limits = "+tmin+" "+tmax);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			idbg = Integer.parseInt(parts[0].trim()); System.out.println(" debugging level = "+idbg);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			dbgfl = parts[0]; FormatHelper.chknam(dbgfl); System.out.println(" debugging file = "+dbgfl);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			outfl = parts[0]; FormatHelper.chknam(outfl); System.out.println(" output file = "+outfl);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nx = Integer.parseInt(parts[0].trim()); xmn = Double.parseDouble(parts[1].trim()); xsiz = Double.parseDouble(parts[2].trim());
			System.out.println(" nx, xmn, xsiz = "+nx+"  "+xmn+"  "+xsiz);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ny = Integer.parseInt(parts[0].trim()); ymn = Double.parseDouble(parts[1].trim()); ysiz = Double.parseDouble(parts[2].trim());
			System.out.println(" ny, ymn, ysiz = "+ny+"  "+ymn+"  "+ysiz);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nxdis = Integer.parseInt(parts[0].trim()); nydis = Integer.parseInt(parts[1].trim());
			System.out.println(" discretization = "+nxdis+"  "+nydis);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ndmin = Integer.parseInt(parts[0].trim()); ndmax = Integer.parseInt(parts[1].trim());
			if(ndmin<0) ndmin = 0;
			System.out.println(" min max data = "+ndmin+"  "+ndmax);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			radius = Double.parseDouble(parts[0].trim()); System.out.println(" isotropic radius = "+radius);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			ktype = Integer.parseInt(parts[0].trim()); skmea = Double.parseDouble(parts[1].trim());
			System.out.println(" ktype,skmea = "+ktype+" "+skmea);

			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nst = Integer.parseInt(parts[0].trim()); c0 = Double.parseDouble(parts[1].trim());
			System.out.println(" nst, nugget = "+nst+"  "+c0);

			if(nst<=0) {
				nst     = 1;
				it[0]   = 1;
				cc[0]   = 0.0d;
				ang[0]  = 0.0d;
				aa[0]   = 0.0d;
				anis[0] = 0.0d;
			} else {
				for(int i=0; i<nst; i++) {
					line = br.readLine(); parts = FormatHelper.splitBySpace(line);
					it[i]  = Integer.parseInt(parts[0].trim());
					cc[i]  = Double.parseDouble(parts[1].trim());
					ang[i] = Double.parseDouble(parts[2].trim());
					aa[i]  = Double.parseDouble(parts[3].trim());
					double a2 = Double.parseDouble(parts[4].trim());
					anis[i] = a2 / aa[i];
					System.out.println(" it,cc,ang,a_max,a_min = "+it[i]+", "+cc[i]+", "+ang[i]+", "+aa[i]+", "+a2);
					if(it[i]==4) {
						if(aa[i]<0.0d) throw new RuntimeException(" INVALID power variogram");
						if(aa[i]>2.0d) throw new RuntimeException(" INVALID power variogram");
					}
					if(it[i]==5) {
						System.out.println(" we do not support this variogram\n"
								          +" use kt3d instead!");
						return;
					}
				}
			}
		}catch(IOException io_e) {
			throw new RuntimeException("ERROR in parameter file!");
		}
//        c
//        c Find the needed parameters:
//        c
		MAXSAM = ndmax + 1;
		MAXDIS = nxdis * nydis;
		MAXKD = MAXSAM + 1;
		MAXKRG = MAXKD * MAXKD;
		double av,ss;
//        c
		if(nst>MAXNST) throw new RuntimeException("nst is too big .. modify .inc file");
//        c
//        c Check to make sure the data file exists, then either read in the
//        c data or write an error message and stop:
//        c
		f = new File(datafl);
		if(!f.exists()) {
			System.err.println("ERROR data file ',datafl,' does not exist!");
			return;
		}
//        c
//        c The data file exists so open the file and read in the header
//        c information. Initialize the storage that will be used to summarize
//        c the data found in the file:
//        c
		try(BufferedReader br = new BufferedReader(new FileReader(f))) {
			br.mark(1024*1024);
			String line = br.readLine();
			line = br.readLine(); String[] parts = FormatHelper.splitBySpace(line);
			int nvari = Integer.parseInt(parts[0].trim());
			for(int i=0; i<nvari; i++) {
				line = br.readLine();
			}
			int MAXDAT = 0;
			while(true) {
				line = br.readLine(); if(line==null) break;
				parts = FormatHelper.splitBySpace(line);
				for(int j=0; j<nvari; j++)
					var[j] = Double.parseDouble(parts[j].trim());
				if(var[ivrl-1]<tmin || var[ivrl-1]>=tmax) continue;
				MAXDAT++;
			}
//        c
//        c Allocate the needed memory:
//        c
			x = new double[MAXDAT];
//        c
			y = new double[MAXDAT];
//        c
			vr = new double[MAXDAT];
//        c
			br.reset();
			line = br.readLine();
			line = br.readLine(); parts = FormatHelper.splitBySpace(line);
			nvari = Integer.parseInt(parts[0].trim());
			av = 0d;
			ss = 0d;
			for(int i=0; i<nvari; i++) {
				line = br.readLine();
			}
//        c
//        c Read the data:
//        c
			nd = -1;
			while(true) {
				line = br.readLine(); if(line==null) break;
				parts = FormatHelper.splitBySpace(line);
				for(int j=0; j<nvari; j++)
					var[j] = Double.parseDouble(parts[j].trim());
				double vrt = var[ivrl-1];
				if(vrt<tmin || vrt>=tmax) continue;
				nd++;
				x[nd]  = var[ixl-1];
				y[nd]  = var[iyl-1];
				vr[nd] = vrt;
				av    += vrt;
				ss    += vrt*vrt;
			}
			nd++;
		} catch(IOException io_e) {
			throw new RuntimeException("ERROR in data file!");
		}
//        c
//        c Open the output files:
//        c
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfl)))) {
			bw.append("KB2D Output\n");
			int    nz    = 0;
			double zmn   = 0d;
			double zsiz  = 1d;
			int nreal = 1;
			bw.append(" 2 "+FormatHelper.nf(nx,4)+" "+FormatHelper.nf(ny,4)+" "+FormatHelper.nf(nz,4)
				     +" "+FormatHelper.nf(xmn,14,8)+" "+FormatHelper.nf(ymn,14,8)+" "+FormatHelper.nf(zmn,14,8)
				     +" "+FormatHelper.nf(xsiz,12,6)+" "+FormatHelper.nf(ysiz,12,6)+" "+FormatHelper.nf(zsiz,12,6)
				     +" "+FormatHelper.nf(nreal,4)+"\n");
			bw.append("Estimate\nEstimate Variance\n");
		 	bw.flush();
		}catch(IOException io_e) {
			io_e.printStackTrace();
			return;
		}

		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl)))) {
			bw.append("... init debug-file ...\n");
			bw.flush();
		} catch(IOException io_e) {
			io_e.printStackTrace();
			return;
		}
//        c
//        c Compute the averages and variances as an error check for the user:
//        c
		av /= Math.max(nd,1d);
		ss = (ss / Math.max(nd,1d)) - av * av;
//        c
//        c Write Some of the Statistics to the screen:
//        c
		System.out.println("   There are "+FormatHelper.nf(nd,8)+" data with:\n"
				          +"   mean value          = "+FormatHelper.nf(av,12,5)+"\n"
				          +"   standard deviation  = "+FormatHelper.nf(Math.sqrt(Math.max(ss,0d)),12,5));
		System.out.println("\n\nFINISHED READING\n");
		
		FormatHelper.printTable(10, x,y,vr);
	}

	public void kb2d() {
		System.out.println("[DEBUG]\n"
				+          "   '-> debug-file = "+dbgfl+"  (is null "+(dbgfl==null)+")");
//        c-----------------------------------------------------------------------
//        c
//        c           Ordinary/Simple Kriging of a 2-D Rectangular Grid
//        c           *************************************************
//        c
//        c This subroutine estimates point or block values of one variable by
//        c ordinary kriging.  All of the samples are rescanned for each block
//        c estimate; this makes the program simple but inefficient.  The data
//        c should NOT contain any missing values.  Unestimated points are
//        c returned as -1.0e21
//        c
//        c
//        c
//        c Original:  A.G. Journel                                           1978
//        c Revisions: B.E. Buxton                                       Apr. 1983
//        c-----------------------------------------------------------------------
		//    use geostat
		double[] xdb,ydb,xa,ya,vra,dist;
		double[] r,rr,s,a;
		int[] nums;
//        c
		boolean first = true;
		double PMX = 9999.0d;
		System.out.println("[DEBUG] MAXDIS="+MAXDIS+"  MAXSAM="+MAXSAM+"  MAXKD="+MAXKD+"  MAXKRG="+MAXKRG);
//        c
//        c Echo the input parameters if debugging flag is >2:
//        c
//        c Allocate the needed memory:
//        c
		xdb = new double[MAXDIS];
//        c
		ydb = new double[MAXDIS];
//        c
		xa = new double[MAXSAM];
//        c
		ya = new double[MAXSAM];
//        c
		vra = new double[MAXSAM];
//        c
		dist = new double[MAXSAM];
//        c
		nums = new int[MAXSAM];
//        c
		r = new double[MAXKD];
//        c
		rr = new double[MAXKD];
//        c
		s = new double[MAXKD];
//        c
		a = new double[MAXKRG];
//        c
		if(idbg>2) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				bw.append("KB2D Parameters\n\n");
				bw.append("Variogram Parameters for "+nst+" structures:\n"
						 +"  Nugget effect:            "+c0+"\n"
						 +"  Types of variograms:  ");
				for(int i=0; i<nst; i++) bw.append(" "+it[i]);
				bw.append("\n  Contribution cc       ");
				for(int i=0; i<nst; i++) bw.append(" "+cc[i]);
				bw.append("\n  Ranges:               ");
				for(int i=0; i<nst; i++) bw.append(" "+aa[i]);
				bw.append("\n  Angle for Continuity: ");
				for(int i=0; i<nst; i++) bw.append(" "+ang[i]);
				bw.append("\n  Anisotropy Factors:   ");
				for(int i=0; i<nst; i++) bw.append(" "+anis[i]);
				bw.append("\n\nGrid for Kriging:\n");
				bw.append("  Number of X and Y Blocks: "+nx+" "+ny+"\n");
				bw.append("  Origin of X and Y Blocks: "+xmn+" "+ymn+"\n");
				bw.append("  Size   of X and Y Blocks: "+xsiz+" "+ysiz+"\n\n");
				bw.append("Discretization of blocks:   "+nxdis+" "+nydis+"\n");
				bw.append("Search Radius:              "+radius+"\n");
				bw.append("Minimum number of samples:  "+ndmin+"\n");
				bw.append("Maximum number of samples:  "+ndmax+"\n");
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
//        c
//        c Echo the input data if debugging flag >1:
//        c
		if(idbg>=4) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				for(int id=0; id<nd; id++) {
					bw.append("Data: "+FormatHelper.nf(id+1,5)+" at "+FormatHelper.nf(x[id],12,3)+" "+FormatHelper.nf(y[id],12,3)+
							  " value: "+FormatHelper.nf(vr[id],12,5));
				}
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
//        c
//        c Set up the discretization points per block.  Figure out how many
//        c are needed, the spacing, and fill the xdb and ydb arrays with the
//        c offsets relative to the block center (this only gets done once):
//        c
		int ndb  = nxdis * nydis;
		if(ndb>MAXDIS) {
			throw new RuntimeException("ERROR KB2D: Too many discretization points \n"
					+                  "            Increase MAXDIS or lower n[xy]dis");
		}
		double xdis = xsiz  / Math.max(nxdis,1d);
		double ydis = ysiz  / Math.max(nydis,1d);
		double xloc = -0.5d*(xsiz+xdis);
		int ii   = -1;
		for(int ix=0; ix<nxdis; ix++) {
			xloc = xloc + xdis;
			double yloc = -0.5d*(ysiz+ydis);
			for(int iy=0; iy<nydis; iy++) {
				yloc = yloc + ydis;
				ii++;
				xdb[ii] = xloc;
				ydb[ii] = yloc;
			}
		}
		//System.out.println("[DEBUG] xdb ydb");
		//FormatHelper.printTable(30, xdb,ydb);
//        c
//        c Initialize accumulators:
//        c
		double cb;
		double cbb = 0d;
		double rad2 = radius*radius;
//        c
//        c Calculate Block Covariance. Check for point kriging.
//        c
		double cov = cova2(xdb[0],ydb[0],xdb[0],ydb[0],nst,c0,PMX,cc,aa,it,ang,anis,first);
//        c
//        c Keep this value to use for the unbiasedness constraint:
//        c
		double unbias = cov;
		first  = false;
		if (ndb<=1) {
			cbb = cov;
		} else {
			for(int i=0; i<ndb; i++) for(int j=0; j<ndb; j++) {
				cov = cova2(xdb[i],ydb[i],xdb[j],ydb[j],nst,c0,PMX,cc,aa,it,ang,anis,first);
				if(i==j) cov -= c0;
				cbb += cov;
			}
			cbb /= ndb*ndb;
		}
		if(idbg>1) {
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				bw.append("\nBlock Covariance: "+cbb+"\n");
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
//        c
//        c MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
//        c
		int    nk = 0;
		double ak = 0d,
			   vk = 0d;
		double est,estv;
		for(int iy=0; iy<ny; iy++) { // Loop ID: 4
			double yloc = ymn + iy*ysiz;
			for(int ix=0; ix<nx; ix++) { // Loop ID: 4
				xloc = xmn + ix*xsiz;
//        c
//        c Find the nearest samples within each octant: First initialize
//        c the counter arrays:
//        c
				int na = -1;
				for(int isam=0; isam<ndmax; isam++) {
					dist[isam] = 1.0e+20d;
					nums[isam] = 0;
				}
//        c
//        c Scan all the samples (this is inefficient and the user with lots of
//        c data should move to ktb3d):
//        c
				for(int id=0; id<nd; id++) { // Loop ID: 6
					double dx = x[id] - xloc;
					double dy = y[id] - yloc;
					double h2 = dx*dx + dy*dy;
					if(h2>rad2) continue;
//        c
//        c Do not consider this sample if there are enough close ones:
//        c
					//System.out.println("[TEST] na="+na+"  h2="+h2+"  dist[na]="+dist[na]);
					if(na==ndmax-1 && h2>dist[na]) continue;
//        c
//        c Consider this sample (it will be added in the correct location):
//        c
					if(na<ndmax-1) na++;
					nums[na]           = id;
					dist[na]           = h2;
					if(na==0) continue;
//        c
//        c Sort samples found thus far in increasing order of distance:
//        c
					//int n1 = na-1;
					//for(ii=0; ii<na; ii++) {
					//	int k=ii;
					//	if(h2<dist[ii]) {
					//		int jk = 0;
					//		for(int jj=k; jj<n1; jj++) {
					//			int j  = n1-jk;
					//			jk++;
					//			int j1 = j+1;
					//			dist[j1] = dist[j];
					//			nums[j1] = nums[j];
					//		}
					//		dist[k] = h2;
					//		nums[k] = id;
					//		break;
					//	}
					//}
					for(ii=na-1; ii>=0; ii--) {
						if(dist[ii]>dist[ii+1]) {
							int j = nums[ii];
							nums[ii] = nums[ii+1];
							nums[ii+1] = j;
							double tmpdst = dist[ii];
							dist[ii] = dist[ii+1];
							dist[ii+1] = tmpdst;
						}
					}
				}
				na++;
				System.out.println("    ... for cell "+(ix+1)+"|"+(iy+1)+" with pos "+xloc+" "+yloc+"  found "+na+" points");
//        c
//        c Is there enough samples?
//        c
				if(na<ndmin) {
					if(idbg>=2) {
						try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
							bw.append("Block "+(ix+1)+" "+(iy+1)+" not estimated\n");
							bw.flush();
						}catch(IOException io_e) {
							io_e.printStackTrace();
						}
					}
					System.out.println("Block "+(ix+1)+" "+(iy+1)+" not estimated due to too few points!");
					est  = UNEST;
					estv = UNEST;
					continue;
				}
//        c
//        c Put coordinates and values of neighborhood samples into xa,ya,vra:
//        c
				for(int ia=0; ia<na; ia++) {
					int jj      = nums[ia];
					xa[ia]  = x[jj];
					ya[ia]  = y[jj];
					vra[ia] = vr[jj];
				}
				//System.out.println("[DEBUG] xa, ya, vra, dist    for "+na+" points:");
				//FormatHelper.printTable(na, xa,ya,vra,dist);
//        c
//        c Handle the situation of only one sample:
//        c
				if(na==1) {
					double cb1 = cova2(xa[0],ya[0],xa[0],ya[0],nst,c0,PMX,cc,aa,it,ang,anis,first);
					double xx  = xa[0] - xloc;
					double yy  = ya[0] - yloc;
//        c
//        c Establish Right Hand Side Covariance:
//        c
					if(ndb<=1) {
						cb = cova2(xx,yy,xdb[0],ydb[0],nst,c0,PMX,cc,aa,it,ang,anis,first);
					} else {
						cb  = 0.0d;
						for(int i=0; i<ndb; i++) {
							cb += cova2(xx,yy,xdb[i],ydb[i],nst,c0,PMX,cc,aa,it,ang,anis,first);
							double dx = xx - xdb[i];
							double dy = yy - ydb[i];
							if(dx*dx+dy*dy<EPSLON) cb -= c0;
						}
						cb /= ndb;
					}
					if(ktype==0) {
						s[0] = cb/cbb;
						est  = s[0]*vra[0] + (1d-s[0])*skmea;
						estv = cbb - s[0] * cb;
					} else {
						est  = vra[0];
						estv = cbb - 2d*cb + cb1;
					}
				} else {
//        c
//        c Solve the Kriging System with more than one sample:
//        c
					int neq = na + ktype;
					//int nn  = (neq + 1)*neq/2;
//        c
//        c Set up kriging matrices:
//        c
					int in=-1;
					for(int j=0; j<na; j++) {
//        c
//        c Establish Left Hand Side Covariance Matrix:
//        c
						for(int i=0; i<=j; i++) {
							in++;
							a[in] = cova2(xa[i],ya[i],xa[j],ya[j],nst,c0,PMX,cc,aa,it,ang,anis,first);
						}
						double xx = xa[j] - xloc;
						double yy = ya[j] - yloc;
//        c
//        c Establish Right Hand Side Covariance:
//        c
						if(ndb<=1) {
							cb = cova2(xx,yy,xdb[0],ydb[0],nst,c0,PMX,cc,aa,it,ang,anis,first);
						} else {
							cb  = 0d;
							for(int j1=0; j1<ndb; j1++) {
								cb += cova2(xx,yy,xdb[j1],ydb[j1],nst,c0,PMX,cc,aa,it,ang,anis,first);
								double dx = xx - xdb[j1];
								double dy = yy - ydb[j1];
								if(dx*dx+dy*dy<EPSLON)
									cb -= c0;
							}
							cb /= ndb;
						}
						r[j]  = cb;
						rr[j] = r[j];
					}
//        c
//        c Set the unbiasedness constraint:
//        c
					if(ktype==1) {
						for(int i=0; i<na; i++) {
							in++;
							a[in] = unbias;
						}
						in++;
						a[in]   = 0d;
						r[neq-1]  = unbias;
						rr[neq-1] = r[neq-1];
					}

					//extra:
					//for(int ai=0; ai<neq*(neq+1)/2; ai++) a[ai] /= cbb;
					//for(int ri=0; ri<neq; ri++) r[ri] /= cbb;
//        c
//        c Write out the kriging Matrix if Seriously Debugging:
//        c
					if(idbg>=3) {
						try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
							bw.append("\nKriging Matrices for Node: "+FormatHelper.nf(ix+1,4)+" "+FormatHelper.nf(iy+1,4)+" RHS first\n");
							int is = 0;
							for(int i=0; i<neq; i++) {
								int ie = is + i;
								bw.append("  r("+FormatHelper.nf(i+1,2)+") = "+FormatHelper.nf(r[i],12,4)+"  a= ");
								for(int j=is; j<=ie; j++) bw.append(" "+FormatHelper.nf(a[j],12,4));
								bw.append("\n");
								is += i+1;
							}
						}catch(IOException io_e) {
							io_e.printStackTrace();
						}
					}
//        c
//        c Solve the Kriging System:
//        c
					s = LinEquSolver.ksol(1,neq,1,a,r);
//        c
//        c Write a warning if the matrix is singular:
//        c
					if(s==null) {
						try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
							bw.append("WARNING KB2D: singular matrix\n"
									+ "              for block "+(ix+1)+" "+(iy+1)+"\n");
						} catch(IOException io_e) {
							io_e.printStackTrace();
						}
						est  = UNEST;
						estv = UNEST;
						System.out.println("  WARNING! singular matrix or error in calculation!");
						continue;
					}
//        c
//        c Write the kriging weights and data if requested:
//        c
					if(idbg>=2) {
						try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
							bw.append("\nBLOCK: "+(ix+1)+" "+(iy+1)+"\n\n");
							if(ktype==1) bw.append("  Lagrange multiplier: "+(s[neq-1]*unbias)+"\n");
							bw.append("  BLOCK EST: x  y  vr  wt\n");
							for(int i=0; i<na; i++) {
								bw.append(FormatHelper.nf(xa[i],12,3)+" "+FormatHelper.nf(ya[i],12,3)+" "+
										  FormatHelper.nf(vra[i],12,3)+" "+FormatHelper.nf(s[i],12,3)+"\n");
							}
							bw.flush();
						} catch(IOException io_e) {
							io_e.printStackTrace();
						}
					}
//        c
//        c Compute the estimate and the kriging variance:
//        c
					est  = 0.0d;
					estv = cbb;
					double sumw = 0.0d;
					if(ktype==1) estv -= s[na+1]*unbias;
					for(int i=0; i<=na; i++) {
						sumw += s[i];
						est  += s[i]*vra[i];
						estv -= s[i]*rr[i];
					}
					if(ktype==0) est += (1.0-sumw)*skmea;
				}
				if(idbg>=2) {
					try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
						bw.append("  est  "+est+"\n"
								+ "  estv "+estv+"\n");
						bw.flush();
					}catch(IOException io_e) {
						io_e.printStackTrace();
					}
				}
//        c
//        c Write the result to the output file:
//        c
				try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfl),true))) {
					bw.append(FormatHelper.nf(est,14,8)+" "+FormatHelper.nf(estv,14,8)+"\n");
					bw.flush();
				}catch(IOException io_e) {
					io_e.printStackTrace();
				}
				if(est>UNEST) {
					nk++;
					ak += est;
					vk += est*est;
				}
//        c
//        c END OF MAIN LOOP OVER ALL THE BLOCKS:
//        c
			}
		}
		if(nk>=1) {
			ak /= nk;
			vk = vk/nk - ak*ak;
			try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dbgfl),true))) {
				String line = "  Estimated   "+FormatHelper.nf(nk,8)+" blocks\n"+
						      "  average     "+FormatHelper.nf(ak,9,4)+"\n"+
						      "  variance    "+FormatHelper.nf(vk,9,4);
				bw.append(line);
				System.out.println(line);
				bw.flush();
			}catch(IOException io_e) {
				io_e.printStackTrace();
			}
		}
		System.out.println("\n\nFINISHED!");
	}

	public double cova2(double x1, double y1, double x2, double y2, int nst, double c0,
			double PMX, double[] cc, double[] aa, int[] it, double[] ang, double[] anis, boolean first) {
//        c-----------------------------------------------------------------------
//        c
//        c              Covariance Between Two Points (2-D Version)
//        c              *******************************************
//        c
//        c This function returns the covariance associated with a variogram model
//        c that is specified by a nugget effect and possibly four different
//        c nested varigoram structures.  The anisotropy definition can be
//        c different for each of the nested structures (spherical, exponential,
//        c gaussian, or power).
//        c
//        c
//        c
//        c INPUT VARIABLES:
//        c
//        c   x1,y1            Coordinates of first point
//        c   x2,y2            Coordinates of second point
//        c   nst              Number of nested structures (max. 4).
//        c   c0               Nugget constant (isotropic).
//        c   PMX              Maximum variogram value needed for kriging when
//        c                      using power model.  A unique value of PMX is
//        c                      used for all nested structures which use the
//        c                      power model.  therefore, PMX should be chosen
//        c                      large enough to account for the largest single
//        c                      structure which uses the power model.
//        c   cc(nst)          Multiplicative factor of each nested structure.
//        c   aa(nst)          Parameter "a" of each nested structure.
//        c   it(nst)          Type of each nested structure:
//        c                      1. spherical model of range a;
//        c                      2. exponential model of parameter a;
//        c                           i.e. practical range is 3a
//        c                      3. gaussian model of parameter a;
//        c                           i.e. practical range is a*sqrt(3)
//        c                      4. power model of power a (a must be gt. 0  and
//        c                           lt. 2).  if linear model, a=1,c=slope.
//        c   ang(nst)         Azimuth angle for the principal direction of
//        c                      continuity (measured clockwise in degrees from Y)
//        c   anis(nst)        Anisotropy (radius in minor direction at 90 degrees
//        c                      from "ang" divided by the principal radius in 
//        c                      direction "ang")
//        c   first            A logical variable which is set to true if the
//        c                      direction specifications have changed - causes
//        c                      the rotation matrices to be recomputed.
//        c
//        c
//        c
//        c OUTPUT VARIABLES: returns "cova2" the covariance obtained from the
//        c                   variogram model.
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		double DTOR = Math.PI/180d,
			   EPSLON_S = 0.0000001d;
//        c
//        c The first time around, re-initialize the cosine matrix for the
//        c variogram structures:
//        c
		if(first) {
			maxcov = c0;
			for(int is=0; is<nst; is++) {
				double azmuth = (90d-ang[is])*DTOR;
				rotmat[0][is]  =  Math.cos(azmuth);
				rotmat[1][is]  =  Math.sin(azmuth);
				rotmat[2][is]  = -Math.sin(azmuth);
				rotmat[3][is]  =  Math.cos(azmuth);
				if(it[is]==4) {
					maxcov = maxcov + PMX;
				} else {
					maxcov = maxcov + cc[is];
				}
			}
		}
//        c
//        c Check for very small distance:
//        c
		double dx = x2-x1;
		double dy = y2-y1;
		if(dx*dx+dy*dy < EPSLON_S) {
			return maxcov;
		}
//        c
//        c Non-zero distance, loop over all the structures:
//        c
		double cova2 = 0d;
		for(int is=0; is<nst; is++) {
//        c
//        c Compute the appropriate structural distance:
//        c
			double dx1 =  dx*rotmat[0][is] + dy*rotmat[1][is];
			double dy1 = (dx*rotmat[2][is] + dy*rotmat[3][is])/anis[is];
			double h   = Math.sqrt(Math.max(dx1*dx1+dy1*dy1,0d));
			if(it[is]==Covariance.VARIOGRAM_SPHERICAL) {
//        c
//        c Spherical model:
//        c
				double hr = h/aa[is];
				if(hr<1d) cova2 += cc[is]*(1d-hr*(1.5d-.5d*hr*hr));
			} else if(it[is]==Covariance.VARIOGRAM_EXPONENTIAL) {
//        c
//        c Exponential model:
//        c
				cova2 += cc[is]*Math.exp(-3.0d*h/aa[is]);
			} else if(it[is]==Covariance.VARIOGRAM_GAUSSIAN) {
//        c
//        c Gaussian model:
//        c
				double hh=-3.0d*(h*h)/(aa[is]*aa[is]);
				cova2 += cc[is]*Math.exp(hh);
			} else {
//        c
//        c Power model:
//        c
				double cov1  = PMX - cc[is]*Math.pow(h,aa[is]);
				cova2 += cov1;
			}
		}
		return cova2;
	}

	public void makepar() {
//        c-----------------------------------------------------------------------
//        c
//        c                      Write a Parameter File
//        c                      **********************
//        c
//        c
//        c
//        c-----------------------------------------------------------------------
		//    lun = 99
		try(BufferedWriter bw = new BufferedWriter(new FileWriter(new File("res/kb2d.par")))) {
			bw.append("                  Parameters for KB2D\n"
					+ "                  *******************\n"
					+ "START OF PARAMETERS:\n"
					+ "res/cluster.dat              -file with data\n"
					+ "1   2   3                    -columns for X, Y, and variable\n"
					+ "-1.0e21   1.0e21             -trimming limits\n"
					+ "3                            -debugging level: 0,1,2,3\n"
					+ "res/kb2d.dbg                 -file for debugging output\n"
					+ "res/kb2d.out                 -file for kriged output\n"
					+ "2    0.0   25.00             -nx,xmn,xsiz\n"
					+ "2    0.0   25.00             -ny,ymn,ysiz\n"
					+ "1    1                       -x and y block discretization\n"
					+ "2    4                       -min and max data for kriging\n"
					+ "20.0                         -maximum search radius\n"
					+ "1    2.302                   -0=SK, 1=OK,  (mean if SK)\n"
					+ "1   2.0                      -nst, nugget effect\n"
					+ "1   8.0  0.0  10.0  10.0     -it, c, azm, a_max, a_min\n");
			bw.flush();
		} catch(IOException io_e) {
			io_e.printStackTrace();
		}
	}
}
