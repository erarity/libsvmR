
libsvmR.dll: svm.cpp svm.h
	R CMD SHLIB -o libsvmR.dll svm.cpp 
clean:
	rm -f *~ *.o *.so
tags:
	etags svm.cpp
