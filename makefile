# Makefile

BINDIR := bin
LIBDIR := lib

CCFLAGS := -pedantic

CC := g++ -std=c++17

# src/ (declaracoes de funcoes, de classes + codigo)
# main/ (programas principais)
# bin/ (temporarios, .o, .exe)
# lib/ (bibliotecas) biblioteca FC

# making library
# - static: .a
# - shareable: .so

VPATH = main:src

ROOTLIB := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

SRC := $(wildcard src/*.C)
OBJ := $(patsubst %.C, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.h)

lib: $(LIBDIR)/libFC.a

$(LIBDIR)/libFC.a: $(OBJ) 
	@echo make lib...
	ar ruv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libFC.a 
	@echo compilink and linking... 
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l FC $(ROOTLIB) -lGeom -lEG -lpthread

$(BINDIR)/%.o: %.C | $(INC)
	@echo compiling... $<
	$(CC) -I src -I $(ROOTINC) -c -lpthread $< -o $@


######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(LIBDIR)/libFC.a	
