OUT_DIR = ./Release

MAIN_OUT_DIR = ./Release

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main_DIR = $(PWD)
EXAM_DIR = ../

COMMON_DIR = ../../../common_src

export COMMON_DIR EXAM_DIR main_DIR
export OUT_DIR MAIN_OUT_DIR


 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#(1). rays
SRC_DIR = $(PWD)/src
MAIN_DIR = $(PWD)


export MAIN_DIR SRC_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 

#FC = gcc
FC = mpicc #gcc #arm-linux-gnu-gcc
export FC

bin = exe
 
export bin
 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sphclmotion: $(SRC_DIR)/*.c
	cd $(SRC_DIR) && make sphclmotion

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#.PHONY: clean
  
############################################################
clean:
	rm ./src/Release/*.o
	rm ../../sphclmotion
############################################################

 




