# This is an commentary line in a makefile
# Start of the makefile


#
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
#
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#.SUFFIXES:
#.SUFFIXES: .f .F .o .a .f90 .F90

#SIESTA_ARCH=x86_64-ifort12-mpiifort

FPP=
FPP_OUTPUT=
#FC=mpiifort
#FC_SERIAL=ifort
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)
 
 

LIBS=$(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(NETCDF_LIBS)

  

  
OBJS = $(MAIN_OUT_DIR)/root34.o \
	$(MAIN_OUT_DIR)/1RandUtils.o \
	$(MAIN_OUT_DIR)/Carlsons.o \
	$(MAIN_OUT_DIR)/ellfunctions.o \
	$(MAIN_OUT_DIR)/ellCarlsons.o \
	$(MAIN_OUT_DIR)/particle.o \
	$(MAIN_OUT_DIR)/rmu_tps.o \
	$(MAIN_OUT_DIR)/BLcoordinates_new.o \
	$(MAIN_OUT_DIR)/ynogkini.o \
	$(MAIN_OUT_DIR)/ynogk_0aspin.o \
	$(MAIN_OUT_DIR)/ynogk_theta_part.o \
	$(MAIN_OUT_DIR)/ynogkBL.o \
	$(MAIN_OUT_DIR)/main.o
 


#*****************Making directions******************************* 
$(shell if [ ! -e $(MAIN_OUT_DIR) ];then mkdir -p $(MAIN_OUT_DIR); fi)
#$(shell if [ ! -e $(MOD_DIR) ];then mkdir -p $(MOD_DIR); fi) 
#*****************Making directions******************************* 

 

rays: $(OBJS)
	$(FC) $^ -o ../../../$@ -lm -g -O3 -Wall -std=c11


$(MAIN_OUT_DIR)/%.o: $(COMMON_DIR)/%.c
	$(FC) -c -O3 -lm -Wall -ld -g -std=c11 -I$(COMMON_DIR) -I$(MAIN_OUT_DIR) $< -o $@

$(MAIN_OUT_DIR)/%.o: ./%.c
	$(FC) -c -O3 -lm -Wall -ld -g -std=c11 -I$(COMMON_DIR) -I$(MAIN_OUT_DIR) $< -o $@


 
# End of the makefile
