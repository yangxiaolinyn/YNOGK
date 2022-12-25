

OUT_DIR = Release
MOD_DIR = module
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
GRAYS_DIR = $(PWD)/src
GRAY_DIR = $(PWD)


export GRAY_DIR GRAYS_DIR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
 





#FC = gcc
FC = mpicc #gcc #arm-linux-gnu-gcc
export FC

bin = exe
 
export bin
 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grays: $(GRAYS_DIR)/*.c $(GRAYS_DIR)/*.h
	cd $(GRAYS_DIR) && make grays

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#.PHONY: clean
  
############################################################
clrays:
	rm ./src/$(OUT_DIR)/*.o
	rm ../../grays
############################################################

 




