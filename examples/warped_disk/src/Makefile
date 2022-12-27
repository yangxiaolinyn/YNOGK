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

  

  
OBJS =  $(OUT_DIR)/root34.o \
	$(OUT_DIR)/1RandUtils.o \
	$(OUT_DIR)/Carlsons.o \
	$(OUT_DIR)/ellfunctions.o \
	$(OUT_DIR)/ellCarlsons.o \
	$(OUT_DIR)/particle.o \
	$(OUT_DIR)/rmu_tps.o \
	$(OUT_DIR)/BLcoordinates_new.o \
	$(OUT_DIR)/ynogkini.o \
	$(OUT_DIR)/ynogk_0aspin.o \
	$(OUT_DIR)/ynogk_theta_part.o \
	$(OUT_DIR)/ynogkBL.o \
	$(OUT_DIR)/pemfindingUnit.o \
	$(OUT_DIR)/main.o



$(shell if [ ! -e $(OUT_DIR) ];then mkdir -p $(OUT_DIR); fi)
#$(shell if [ ! -e $(MOD_DIR) ];then mkdir -p $(MOD_DIR); fi) 

wpddisk: $(OBJS)
	$(FC) $^ -o ../../../$@ -lm -g -O3 -Wall -std=c11
 
 
$(OUT_DIR)/%.o: $(COMMON_DIR)/%.c Makefile 
	$(FC) -c -O3 -lm -Wall -ld -g -std=c11 -I$(COMMON_DIR) -I$(OUT_DIR) $< -o $@

$(OUT_DIR)/%.o: %.c
	$(FC) -c -O3 -lm -Wall -ld -g -std=c11 -I$(COMMON_DIR) -I$(OUT_DIR) $< -o $@

 
# End of the makefile
