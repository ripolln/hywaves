# ------------------------------------------------------------------------------
#                      Makefile for building SWASH program and documentation
# ------------------------------------------------------------------------------
#
# Before compilation, type "make config" first!
#
# To compile the serial executable type "make ser"
# To compile the OpenMP executable type "make omp"
# To compile the MPI executable type "make mpi"
#
# To remove compiled objects and modules: type "make clean"
#
# To remove compiled objects, modules and executable: type "make clobber"
#
# To compile the SWASH documentation type "make doc"
#
# To remove the PDF and HTML documents type "make cleandoc"
#
# Please do not change anything below, unless you're very sure what you're doing
# ------------------------------------------------------------------------------

include macros.inc

SWASH_EXE = swash.exe

SWASH_OBJS = \
ocpmod.$(EXTO) \
SwashModule1.$(EXTO) \
SwashModule2.$(EXTO) \
m_parall.$(EXTO) \
SwanGriddata.$(EXTO) \
SwanGridobjects.$(EXTO) \
SwanCompdata.$(EXTO) \
SwashFlowdata.$(EXTO) \
SwashSolvedata.$(EXTO) \
Swash.$(EXTO) \
SwashMain.$(EXTO) \
SwashInit.$(EXTO) \
SwashReadInput.$(EXTO) \
SwashReadTestpnts.$(EXTO) \
SwashInputGrid.$(EXTO) \
SwashInputField.$(EXTO) \
SwanReadGrid.$(EXTO) \
SwanReadADCGrid.$(EXTO) \
SwanReadTriangleGrid.$(EXTO) \
SwanReadEasymeshGrid.$(EXTO) \
SwanCheckGrid.$(EXTO) \
SwanCreateEdges.$(EXTO) \
SwanGridTopology.$(EXTO) \
SwanGridVert.$(EXTO) \
SwanGridCell.$(EXTO) \
SwanGridFace.$(EXTO) \
SwanPrintGridInfo.$(EXTO) \
SwashInitCompGrid.$(EXTO) \
SwashInitCompUgrid.$(EXTO) \
SwashInitCond.$(EXTO) \
SwashBackup.$(EXTO) \
SwashBounCond.$(EXTO) \
SwanBpntlist.$(EXTO) \
SwashBCspecfile.$(EXTO) \
SwashBCspectrum.$(EXTO) \
SwashBCshortwave.$(EXTO) \
SwashBCboundwave.$(EXTO) \
SwashIntWavgen.$(EXTO) \
SwashReqOutL.$(EXTO) \
SwashReqOutQ.$(EXTO) \
SwashCheckPrep.$(EXTO) \
SwashGeometrics.$(EXTO) \
SwashSpongeLayer.$(EXTO) \
SwanVertlist.$(EXTO) \
SwashPrintSettings.$(EXTO) \
SwashFlowDP.$(EXTO) \
SwashFloatObjects.$(EXTO) \
SwashUpdateData.$(EXTO) \
SwashReadBndval.$(EXTO) \
SwashUpdateFld.$(EXTO) \
SwashUpdFlowFlds.$(EXTO) \
SwashUpdPress.$(EXTO) \
SwashUpdateDepths.$(EXTO) \
SwashLayerIntfaces.$(EXTO) \
SwashInitSteady.$(EXTO) \
SwashInitBCtrans.$(EXTO) \
SwashServices.$(EXTO) \
SwashSolvers.$(EXTO) \
SwashBotFrict.$(EXTO) \
SwashLogLaw.$(EXTO) \
SwashPorousStruc.$(EXTO) \
SwashPorFricDep.$(EXTO) \
SwashPorFricLay.$(EXTO) \
SwashWindStress.$(EXTO) \
SwashVeget.$(EXTO) \
SwashDensity.$(EXTO) \
SwashHorzVisc.$(EXTO) \
SwashVertVisc.$(EXTO) \
SwashKepsMod1DH.$(EXTO) \
SwashKepsMod2DH.$(EXTO) \
SwashReynoldsStress.$(EXTO) \
SwashComputStruc.$(EXTO) \
SwashComputFlow.$(EXTO) \
SwashComputTrans.$(EXTO) \
SwashComputTurb.$(EXTO) \
SwashDryWet.$(EXTO) \
SwashPresFlow.$(EXTO) \
SwashBreakPoint.$(EXTO) \
SwashAverOutp.$(EXTO) \
SwashImpDep1DHflow.$(EXTO) \
SwashExpDep1DHflow.$(EXTO) \
SwashImpDepM1DHflow.$(EXTO) \
SwashExpDepM1DHflow.$(EXTO) \
SwashImpLay1DHflow.$(EXTO) \
SwashImpLayP1DHflow.$(EXTO) \
SwashExpLay1DHflow.$(EXTO) \
SwashExpLayP1DHflow.$(EXTO) \
SwashImpDep2DHflow.$(EXTO) \
SwashExpDep2DHflow.$(EXTO) \
SwashImpDepM2DHflow.$(EXTO) \
SwashExpDepM2DHflow.$(EXTO) \
SwashImpLay2DHflow.$(EXTO) \
SwashImpLayP2DHflow.$(EXTO) \
SwashExpLay2DHflow.$(EXTO) \
SwashExpLayP2DHflow.$(EXTO) \
SwashExpDep1DHtrans.$(EXTO) \
SwashExpLay1DHtrans.$(EXTO) \
SwashExpDep2DHtrans.$(EXTO) \
SwashExpLay2DHtrans.$(EXTO) \
SwashAntiCreep1DH.$(EXTO) \
SwashAntiCreep2DH.$(EXTO) \
SwashHDiffZplane1DH.$(EXTO) \
SwashHDiffZplane2DH.$(EXTO) \
SwanThreadBounds.$(EXTO) \
SwanFindPoint.$(EXTO) \
SwanPointinMesh.$(EXTO) \
SwashOutput.$(EXTO) \
SwashDecOutL.$(EXTO) \
SwashDecOutQ.$(EXTO) \
SwashCoorOutp.$(EXTO) \
SwashQuanOutp.$(EXTO) \
SwashHydroLoads.$(EXTO) \
SwashRunupHeight.$(EXTO) \
SwanInterpolatePoint.$(EXTO) \
SwanInterpolateOutput.$(EXTO) \
SwashCleanMem.$(EXTO) \
ocpids.$(EXTO) \
ocpcre.$(EXTO) \
ocpmix.$(EXTO) \
swanser.$(EXTO) \
swanout2.$(EXTO) \
swanparll.$(EXTO)

.SUFFIXES: .f .for .f90

.PHONY: help clean clobber

help:
	@echo "This Makefile supports the following:"
	@echo "make config    -- makes machine-dependent macros include file"
	@echo "make ser       -- makes the serial $(SWASH_EXE) executable"
	@echo "make omp       -- makes the OpenMP $(SWASH_EXE) executable"
	@echo "make mpi       -- makes the    MPI $(SWASH_EXE) executable"
	@echo "make doc       -- makes the SWASH documentation (PDF)"
	@echo "make clean     -- removes compiled objects and modules"
	@echo "make clobber   -- removes compiled objects, modules and $(SWASH_EXE)"
	@echo "make cleandoc  -- removes all SWASH documents"

config:
	@perl platform.pl

install:
	@perl platform.pl

ser:
	@perl switch.pl $(swch) *.ftn *.ftn90
	$(MAKE) FOR=$(F90_SER) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_SER)" \
                INCS="$(INCS_SER)" LIBS="$(LIBS_SER)" $(SWASH_EXE)

omp:
	@perl switch.pl $(swch) *.ftn *.ftn90
	$(MAKE) FOR=$(F90_OMP) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP)" \
	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_OMP)" \
                INCS="$(INCS_OMP)" LIBS="$(LIBS_OMP)" $(SWASH_EXE)

mpi:
	@perl switch.pl $(swch) -mpi *.ftn *.ftn90
	$(MAKE) FOR=$(F90_MPI) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI)" \
	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_MPI)" \
                INCS="$(INCS_MPI)" LIBS="$(LIBS_MPI)" $(SWASH_EXE)

doc:
	$(MAKE) -f Makefile.latex TARGET=swashuse doc
	$(MAKE) -f Makefile.latex TARGET=swashtech doc
	$(MAKE) -f Makefile.latex TARGET=swashimp doc

$(SWASH_EXE): $(SWASH_OBJS)
	$(FOR) $(SWASH_OBJS) $(FFLAGS) $(OUT)$(SWASH_EXE) $(INCS) $(LIBS)

.f.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.f90.o:
	$(FOR) $< -c $(FFLAGS90) $(INCS)

.for.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.for.obj:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.f90.obj:
	$(FOR) $< -c $(FFLAGS90) $(INCS)

clean:
	$(RM) *.$(EXTO) *.mod

clobber:
	$(RM) *.$(EXTO) *.mod *.f *.for *.f90 $(SWASH_EXE)

allclean:
	$(RM) *.$(EXTO) *.mod *.f *.for *.f90 $(SWASH_EXE)

cleandoc:
	$(MAKE) -f Makefile.latex TARGET=swashuse cleandoc
	$(MAKE) -f Makefile.latex TARGET=swashtech cleandoc
	$(MAKE) -f Makefile.latex TARGET=swashimp cleandoc
