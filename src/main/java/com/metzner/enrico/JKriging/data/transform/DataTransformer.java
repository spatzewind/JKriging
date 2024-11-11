package com.metzner.enrico.JKriging.data.transform;

import com.metzner.enrico.JKriging.data.DataFrame.DataType;

public interface DataTransformer {
	DataType sourceType();
	DataType destinationType();
	Object transformCell(Object in_cell);
}
