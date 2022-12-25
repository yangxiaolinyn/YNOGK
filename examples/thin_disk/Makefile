

OUT_DIR = ./Release

#Sample_DIR = $(PWD)/src/sampling
#image_DIR = $(PWD)/src/image

#MAIN_OUT_DIR = ../../../common_src/Release
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
disk: $(SRC_DIR)/*.c
	cd $(SRC_DIR) && make disk

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#.PHONY: clean
  
############################################################
clean:
	rm ./src/Release/*.o
	rm ../../disk
############################################################

 




