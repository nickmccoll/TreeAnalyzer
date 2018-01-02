include $(ROOTSYS)/etc/Makefile.arch

#CXXFLAGS += -I../.. -I.
CXXFLAGS += -I.

ifeq ($(PLATFORM),macosx)
CXXFLAGS += -std=c++11
endif


TARGET=libTreeAnalysis.so
STARGET=libTreeAnalysis.a

ASTREEINT    = AnalysisSupport/TreeInterface/src
ASUTILITIES  = AnalysisSupport/Utilities/src
TREEANALYZER = TreeAnalyzer/src
TREEREADING  = TreeReaders/src
DATAFORMATS  = DataFormats/src
CORRECTIONS   = Processors/Corrections/src
GENTOOLS   = Processors/GenTools/src
VARIABLES   = Processors/Variables/src
EVENTSEL   = Processors/EventSelection/src
HCOMBI   = Utilities/HiggsCombineImport/src

SOURCE = $(wildcard $(ASTREEINT)/*.cc) $(wildcard $(ASUTILITIES)/*.cc) $(wildcard $(DATAFORMATS)/*.cc) $(wildcard $(CORRECTIONS)/*.cc) $(wildcard $(TREEREADING)/*.cc)  $(wildcard $(TREEANALYZER)/*.cc) $(wildcard $(GENTOOLS)/*.cc) $(wildcard $(VARIABLES)/*.cc) $(wildcard $(EVENTSEL)/*.cc) $(wildcard $(HCOMBI)/*.cc)
OBJ=$(join $(addsuffix ../obj/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cc=.o))) 
DEPENDS=$(join $(addsuffix ../.dep/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cc=.d)))




all:  $(TARGET) $(STARGET) 
	@true
clean:
	@-rm -f $(TARGET) $(STARGET) $(OBJ) $(DEPENDS) $(HCOMBI)/../obj/TreeAnalyzerDict.o $(HCOMBI)/../obj/TreeAnalyzerDict.cc
        
$(TARGET): $(OBJ) $(HCOMBI)/../obj/TreeAnalyzerDict.o
	@echo "============="
	@echo "Linking the target $@"
	@echo "============="
ifeq ($(PLATFORM),macosx)
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) -lGenVector -l RooFit -lRooFitCore -l RooStats -l Minuit -l Foam -lHistFactory
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS) $(ROOTLIBS)
endif
	@echo -- Link finished --



# ---- setup dictionairy for HiggsCombine import
$(HCOMBI)/../obj/TreeAnalyzerDict.cc: $(wildcard $(HCOMBI)/../interface/*.h) $(HCOMBI)/HiggsCombineImport_LinkDef.h
	@echo "Generating dictionary $@..."
	rootcling  -v -l -f $@ -c -I../.. $^ 
$(HCOMBI)/../obj/TreeAnalyzerDict.o: $(HCOMBI)/../obj/TreeAnalyzerDict.cc
	$(CXX)  $(CXXFLAGS) -c $<  -o $@	
#-----

$(STARGET): $(OBJ)
	@echo "============="
	@echo "Linking the target $@"
	@echo "============="
	ar rcs $@ $^
	@echo -- Link finished --
        
%.o : %.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@

$(ASTREEINT)/../obj/%.o : $(ASTREEINT)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(ASTREEINT)/../.dep/%.d: $(ASTREEINT)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(ASTREEINT)/../obj/$*.o^" > $@'
	
$(ASUTILITIES)/../obj/%.o : $(ASUTILITIES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(ASUTILITIES)/../.dep/%.d: $(ASUTILITIES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(ASUTILITIES)/../obj/$*.o^" > $@'
	
$(DATAFORMATS)/../obj/%.o : $(DATAFORMATS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(DATAFORMATS)/../.dep/%.d: $(DATAFORMATS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(DATAFORMATS)/../obj/$*.o^" > $@'
	
$(CORRECTIONS)/../obj/%.o : $(CORRECTIONS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(CORRECTIONS)/../.dep/%.d: $(CORRECTIONS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(CORRECTIONS)/../obj/$*.o^" > $@'

$(TREEREADING)/../obj/%.o : $(TREEREADING)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(TREEREADING)/../.dep/%.d: $(TREEREADING)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(TREEREADING)/../obj/$*.o^" > $@'

	
$(TREEANALYZER)/../obj/%.o : $(TREEANALYZER)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(TREEANALYZER)/../.dep/%.d: $(TREEANALYZER)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(TREEANALYZER)/../obj/$*.o^" > $@'
	
$(GENTOOLS)/../obj/%.o : $(GENTOOLS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(GENTOOLS)/../.dep/%.d: $(GENTOOLS)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(GENTOOLS)/../obj/$*.o^" > $@'
	
$(VARIABLES)/../obj/%.o : $(VARIABLES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(VARIABLES)/../.dep/%.d: $(VARIABLES)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(VARIABLES)/../obj/$*.o^" > $@'
	
	
$(EVENTSEL)/../obj/%.o : $(EVENTSEL)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(EVENTSEL)/../.dep/%.d: $(EVENTSEL)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(EVENTSEL)/../obj/$*.o^" > $@'
	
$(HCOMBI)/../obj/%.o : $(HCOMBI)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<  -o $@
$(HCOMBI)/../.dep/%.d: $(HCOMBI)/%.cc
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(HCOMBI)/../obj/$*.o^" > $@'

-include $(DEPENDS)
