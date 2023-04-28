package com.metzner.enrico.JKriging.error;

public class UnknownFileFormatException extends Exception {
	
	public static final long serialVersionUID = 4910734366532162991L;
	
	public UnknownFileFormatException() {
		super();
	}
	public UnknownFileFormatException(String message) {
		super(message);
	}
	public UnknownFileFormatException(Throwable cause) {
		super(cause);
	}
	public UnknownFileFormatException(String message, Throwable cause) {
		super(message, cause);
	}
	
}
