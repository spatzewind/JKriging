package com.metzner.enrico.JKriging.Tests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

import com.metzner.enrico.JKriging.data.DataFrame;
import com.metzner.enrico.JKriging.data.reader.DataReader;
import com.metzner.enrico.JKriging.error.UnknownFileFormatException;

@RunWith(Parameterized.class)
public class Test_SliceRead {
	
	@Parameter(0)
	public String dimvar;
	@Parameter(1)
	public int dimidx;
	
	
	@Parameters
	public static Collection<Object[]> data() {
        Object[][] data = new Object[][] {
        	{ "x", 0 },
        	{ "x", 2 },
        	{ "y", 0 },
        	{ "y", 1 }
        };
        List<Object[]> list = Arrays.asList(data);
        Collections.shuffle(list);
        return list;
    }
	
	private DataReader dr;
	
	@Test
	public void test_slice1d() {
		setup();
		
		try {
			System.out.println("JUnitTest: var="+dimvar+", idx="+dimidx);
			DataFrame df = dr.get1Dslice(dimvar, dimidx, "test2d");
			assertNotNull("DataFrame not created.", df);
			
			df.describe();
			df.describeDimension();
			System.out.println("\n");
			
		} catch (Exception iae) {
			iae.printStackTrace();
			fail(iae.getClass().getCanonicalName()+": "+iae.getMessage());
		}
		
		try {
			dr.close();
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getClass().getCanonicalName()+": "+e.getMessage());
		}
	}
	
	
	
	
	
	private void setup() {
		String file = "res/test2d.nc";
		try {
			dr = DataReader.openFile(file);
		} catch(IOException|UnknownFileFormatException e) {
//			e.printStackTrace();
			dr = null;
		}
		assertNotNull("File "+file+" not reachable or readable.", dr);
	}
}
