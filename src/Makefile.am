bin_PROGRAMS = vfgen
vfgen_SOURCES = vfgen.cpp vf.cpp vf_adolc.cpp vf_auto.cpp vf_boostodeint.cpp \
                vf_cvode.cpp vf_cvode7.cpp vf_dde23.cpp vf_dde_solver.cpp    \
                vf_ddebiftool.cpp vf_delay2ode.cpp vf_differentials.cpp      \
                vf_dstool.cpp vf_evf.cpp vf_gsl.cpp                          \
                vf_javamath.cpp vf_javascript.cpp vf_julia.cpp               \
                vf_latex.cpp vf_lsoda.cpp  vf_matcont.cpp vf_matlab.cpp      \
                vf_octave.cpp vf_pddecont.cpp vf_pydstool.cpp vf_pygsl.cpp   \
                vf_r.cpp vf_radau5.cpp vf_scilab.cpp vf_scipy.cpp            \
                vf_taylor.cpp vf_xml.cpp vf_xpp.cpp ginac_aux_functions.cpp  \
                strutils.cpp codegen_utils.cpp MyVec.cpp vf.h vf_help_text.h \
                vf_utils.h codegen_utils.h ginac_aux_functions.h             \
                ginac_declare_funcs.h strutils.h MyVec.h
vfgen_CPPFLAGS = -std=c++17 $(MXML_CFLAGS) $(GINAC_CFLAGS)
vfgen_LDADD = $(MXML_LIBS) $(GINAC_LIBS)
AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = CMakeLists.txt Makefile.vfgen LICENSE-GPL2
