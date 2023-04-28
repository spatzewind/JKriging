package com.metzner.enrico.JKriging.error;

public class DimensionMismatchException extends RuntimeException {

	private static final long serialVersionUID = -4307944225412465625L;

	public DimensionMismatchException() {
	}
	
	public DimensionMismatchException(String message) {
		super(message);
	}
	
	public DimensionMismatchException(Throwable cause) {
		super(cause);
	}
	
	public DimensionMismatchException(String message, Throwable cause) {
		super(message, cause);
	}
}
