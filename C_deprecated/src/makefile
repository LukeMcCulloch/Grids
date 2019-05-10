#trimesh.C makefile
driver:driver.o tlmtrimesh.o TLMsearch.o make_nbrs.o List.o circletest.o
	g++ -g driver.o tlmtrimesh.o TLMsearch.o make_nbrs.o List.o circletest.o -o driver -lm

tlmtrimesh.o:tlmtrimesh.cpp trimesh.h
	g++ -g -c tlmtrimesh.cpp

TLMsearch.o:TLMsearch.cpp
	g++ -g -c TLMsearch.cpp

make_nbrs.o:make_nbrs.c
	g++ -g -c make_nbrs.c

List.o:List.cpp List.h
	g++ -g -c List.cpp

circletest.o:circletest.cpp circletest.h
	g++ -g -c circletest.cpp

.PHONY: clean


clean: 
	-rm -f  driver.o tlmtrimesh.o TLMsearch.o make_nbrs.o List.o circletest.o

realclean: 
	-rm -f  driver driver.o  tlmtrimesh.o TLMsearch.o make_nbrs.o List.o circletest.o
#