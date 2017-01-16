SRC = .
HOME = .
LIB_DIR = lib
PHOENIX_LIB_DIR = lib/phoenix-2.0.0/lib
INCLUDE_DIR = incl
PHOENIX = phoenix 
LIBS = -L$(HOME)/$(LIB_DIR) -L$(HOME)/$(PHOENIX_LIB_DIR) -l$(PHOENIX) -lpthread -lrt
LIB_PHOENIX = libphoenix.a
default: gpbf

mapred_gp.o: 
	mpic++ -D_LINUX_  -O3 -D__x86_64__  -I$(HOME)/$(INCLUDE_DIR) -c $(SRC)/mapred_gp.c  -o mapred_gp.o

gpbf.o:
	mpic++ -D_LINUX_   -O3 -D__x86_64__  -I$(HOME)/$(INCLUDE_DIR) -c $(SRC)/gpbf.cpp -o gpbf.o 

gpbf: gpbf.o
	mpic++ -D_LINUX_  -O3 -D__x86_64__ -o gpbf gpbf.o 
#	mpic++ -D_LINUX_  -O3 -D__x86_64__ -o gpbf gpbf.o mapred_gp.o $(LIBS)

clean:
	rm *.o gpbf

tag:
	ctags -R
sty:
	astyle --style=linux --indent=spaces=2 -p mvas_*.h
	astyle --style=linux --indent=spaces=2 -p *.c
