# Set default shell
SHELL = /bin/sh

# Target document
TARGET = grid_generation

# External programs required (assumed to be in user path)
LATEX = latex
BIBTEX = bibtex
DVIPDF = dvipdf

# Make sure external programs are found
ifneq ($(MAKECMDGOALS),clean)
  ifeq ($(shell if [ -z "`which $(LATEX) 2>/dev/null`" ]; then echo 1; fi),1)
    $(error $(LATEX): command not found.  The $(LATEX) command is required)
  endif
  ifeq ($(shell if [ -z "`which $(BIBTEX) 2>/dev/null`" ]; then echo 1; fi),1)
    $(error $(BIBTEX): command not found.  The $(BIBTEX) command is required)
  endif
  ifeq ($(shell if [ -z "`which $(DVIPDF) 2>/dev/null`" ]; then echo 1; fi),1)
    $(error $(DVIPDF): command not found.  The $(DVIPDF) command is required)
  endif
endif

# Set base GRIDGEN directory
  GGEN = ../.

# IMDIR: location of images files
IMDIR ?= ${GGEN}/images
# IMGZ: name of gzipped images directory
IMGZ=images.tgz

print-%  : ; @echo $* = $($*)

# Make sure INMDIR exists (if not clean)
ifneq ($(MAKECMDGOALS),clean)
#  ifeq ($(shell if [ ! -e $(IMDIR) ]; then echo 1; fi),1)
#    IMDIR = ${GGEN}/images
#  endif
  ifeq ($(shell if [  -d $(IMDIR) ]; then echo 1; fi),1)
    $(shell tar -zxf $(IMGZ))
    ifeq ($(shell if [ ! -d $(IMDIR) ]; then echo 1; fi),1)
      $(error IMDIR=$(IMDIR) does not exist.  Attempt to create images \
        directory from images.tgz file failed.)
    endif
  endif
endif

# TeX files in this directory
TEXFILES = \
grid_generation.tex

# Bibliography inputs
BIBFILES = grid_generation.bib
BSTFILES = jas.bst

# Graphics inputs
GFXFILES = \
images/gshhs_area_distrib.eps \
images/flow_chart.eps \
images/grid_module.eps \
images/DBDB2_Bahamas_1deg_lattitude_transects.eps \
images/boundary_sketch1.eps \
images/boundary_sketch2.eps \
images/boundary_module.eps \
images/boundary_sketch3.eps \
images/global_boundaries.eps \
images/landmask_module.eps \
images/bahamas_landmask.eps \
images/wetcell_module.eps \
images/pacific_landmask.eps \
images/pacific_waterbodies.eps \
images/bahamas_landmask_withboundary.eps \
images/subgrid1.eps \
images/subgrid2.eps \
images/subgrid3.eps \
images/subgrid4.eps \
images/subgrid5.eps \
images/subgrid_module_1.eps \
images/subgrid_module_2.eps \
images/ec_lm_orig.eps \
images/ec_lm_final.eps \
images/Bahamas_basecont.eps \
images/Bahamas_2mincont.eps \
images/Bahamas_4mincont.eps \
images/Bahamas_8mincont.eps \
images/Bahamas_15mincont.eps \
images/Bahamas_30mincont.eps \
images/bahamas/gs_runs/ATLBAH_2min_noobstr.eps \
images/bahamas/newc/ATLBAH_2min_noobstr.eps \
images/bahamas/gs_runs2/ATLBAH_2min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_2min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_2min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_2min_obstr1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob1v2mnnoob.eps \
images/bahamas/gs_runs4/ATLBAH_4min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_8min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_15min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_30min_noobstr.eps \
images/bahamas/gs_runs4/ATLBAH_4min_obstr1.eps \
images/bahamas/gs_runs4/ATLBAH_8min_obstr1.eps \
images/bahamas/gs_runs4/ATLBAH_15min_obstr1.eps \
images/bahamas/gs_runs4/ATLBAH_30min_obstr1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob1v4mnob1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob1v8mnob1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob1v15mnob1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob1v30mnob1.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnnoobv4mnnoob.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnnoobv8mnnoob.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnnoobv15mnnoob.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnnoobv30mnnoob.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob0v4mnob0.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob0v8mnob0.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob0v15mnob0.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnob0v30mnob0.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnobrv4mnobr.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnobrv8mnobr.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnobrv15mnobr.eps \
images/bahamas/gs_runs4/ATLBAHcomp_2mnobrv30mnobr.eps \
images/bahamas/gs_runs3/ATLBAH_8min_obstr1.eps \
images/bahamas/gs_runs3/ATLBAH_8min_obstrr.eps \
images/bahamas/gs_runs3/ATLBAHcomp_8mnob1v8mnobr.eps \
images/Hawaii_basecont.eps \
images/Hawaii_2mincont.eps \
images/Hawaii_4mincont.eps \
images/Hawaii_8mincont.eps \
images/Hawaii_15mincont.eps \
images/Hawaii_30mincont.eps \
images/hawaii/gs_runs4/PACHAW_2min_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_2min_obstr1.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob1v2mnnoob.eps \
images/hawaii/gs_runs4/PACHAW_4min_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_8min_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_15min_noobstr.eps \
images/hawaii/gs_runs3/PACHAW_30min_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_4min_obstr1.eps \
images/hawaii/gs_runs4/PACHAW_8min_obstr1.eps \
images/hawaii/gs_runs4/PACHAW_15min_obstr1.eps \
images/hawaii/gs_runs4/PACHAW_30min_obstr1.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnnoobv4mnnoob.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnnoobv8mnnoob.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnnoobv15mnnoob.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnnoobv30mnnoob.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob0v4mnob0.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob0v8mnob0.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob0v15mnob0.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob0v30mnob0.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnobrv4mnobr.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnobrv8mnobr.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnobrv15mnobr.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnobrv30mnobr.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob1v4mnob1.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob1v8mnob1.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob1v15mnob1.eps \
images/hawaii/gs_runs4/PACHAWcomp_2mnob1v30mnob1.eps \
images/hawaii/gs_runs4/PACHAW_2minzoom_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_30minzoom_noobstr.eps \
images/hawaii/gs_runs4/PACHAW_30minzoom_obst0.eps \
images/hawaii/gs_runs4/PACHAW_30minzoom_obstr.eps \
images/FrenchPoly_basecont.eps \
images/FrenchPoly_2mincont.eps \
images/FrenchPoly_4mincont.eps \
images/FrenchPoly_8mincont.eps \
images/Hawaii_15mincont.eps \
images/Hawaii_30mincont.eps \
images/french_poly/gs_runs4/PACFPY_2min_noobstr.eps \
images/french_poly/gs_runs4/PACFPY_2min_obstr1.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob1v2mnnoob.eps \
images/french_poly/gs_runs4/PACFPY_4min_noobstr.eps \
images/french_poly/gs_runs4/PACFPY_8min_noobstr.eps \
images/french_poly/gs_runs4/PACFPY_15min_noobstr.eps \
images/french_poly/gs_runs4/PACFPY_30min_noobstr.eps \
images/french_poly/gs_runs4/PACFPY_4min_obstr1.eps \
images/french_poly/gs_runs4/PACFPY_8min_obstr1.eps \
images/french_poly/gs_runs4/PACFPY_15min_obstr1.eps \
images/french_poly/gs_runs4/PACFPY_30min_obstr1.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnnoobv4mnnoob.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnnoobv8mnnoob.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnnoobv15mnnoob.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnnoobv30mnnoob.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob0v4mnob0.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob0v8mnob0.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob0v15mnob0.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob0v30mnob0.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnobrv4mnobr.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnobrv8mnobr.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnobrv15mnobr.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnobrv30mnobr.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob1v4mnob1.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob1v8mnob1.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob1v15mnob1.eps \
images/french_poly/gs_runs4/PACFPYcomp_2mnob1v30mnob1.eps \
images/french_poly/frenchpoly_obsrect.eps \
images/french_poly/frenchpoly_obscurv.eps

# All dependencies
DEPFILES = $(TEXFILES) $(BIBFILES) $(BSTFILES) $(STYFILES) $(GFXFILES)

#------------------------------------------------------------------------------#
# Make targets
#------------------------------------------------------------------------------#

$(TARGET): $(TARGET).dvi $(TARGET).pdf

clean:
	\rm -f *.log */*.aux *.aux *.bbl *.blg *.toc *.out *~ */*~
	\rm -f *.dvi *.pdf
	\rm -rf images

$(TARGET).dvi: $(DEPFILES)
	$(LATEX) $(TARGET)
	$(BIBTEX) $(TARGET)
	-@sh bib_strip.sh $(TARGET)
	$(LATEX) $(TARGET)
	$(BIBTEX) $(TARGET)
	-@sh bib_strip.sh $(TARGET)
	$(LATEX) $(TARGET)
	$(BIBTEX) $(TARGET)
	-@sh bib_strip.sh $(TARGET)
	$(LATEX) $(TARGET)

$(TARGET).pdf: $(TARGET).dvi
	$(DVIPDF) $(TARGET)
