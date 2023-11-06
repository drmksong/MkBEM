CC = gcc
CXX = g++
CFLAGS = -c -g -frtti -DMKDEBUG  -DPROPOUT  # -I$(INC) -L$(MKLIB)  -D__STRICT_ANSI__
TFLAGS = -g -frtti -DMKDEBUG  -DPROPOUT -I$(INC) -L$(MKLIB)  -D__STRICT_ANSI__ # 
TARGET = thmbem.exe
INC = d:/MkLib_double
MKLIB = d:/MkLib_double
SRCS = mainunit.cpp inputunit.cpp bemunit.cpp
HEAD = mainunit.h inputunit.h bemunit.h
OBJS = mainunit.o inputunit.o bemunit.o 
MKOBJS = $(MKLIB)/MKmatrix.o $(MKLIB)/MKint.o 

.SUFFIXES : .cpp .o

$(TARGET) : $(OBJS)  
	$(CC) $(TFLAGS) $(OBJS) $(MKOBJS) -o $(TARGET) -I$(INC)  -lsupc++ -lm -lmkgcc

$(OBJS) : $(SRCS) $(HEAD)
	$(CC) $(CFLAGS) -I$(INC) $(SRCS)

.PHONY : all clean                            

clean :
	-rm $(OBJS) $(TARGET)
