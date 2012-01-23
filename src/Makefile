
libsvmR.so: svm.cpp svm.h
	R CMD SHLIB -o libsvmR.so svm.cpp 
clean:
	rm -f *~ *.o *.so
tags:
	etags svm.cpp
