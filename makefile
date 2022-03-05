# Paths to the directories
SRCDIR = src
BINDIR = bin
OBJDIR = $(SRCDIR)/obj
MODDIR = $(SRCDIR)/mod

FILES = utils.f90 print_f.f90 main_sub.f90 magnetic.f90

# Create a list of source files
SRC=$(addprefix $(SRCDIR)/,$(FILES))

# Create a list of object files
OBJS=$(patsubst %.f90,%.o,$(SRC))
OBJ=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(OBJS))

# Compiler/Linker settings
FC = gfortran
ifeq ($(DEBUG),1)
FFLAGS = -J$(MODDIR) -Wall -Wextra -g3 -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow -fbounds-check -Warray-temporaries
else
FFLAGS = -J$(MODDIR) -Wall -fcheck=all
endif
LIBS = -lblas -llapack

# Program info
PROGRAM = magnetic.exe
	
# =====
# Rules
# =====

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@echo "\e[92;1mCOMPILING $<:\e[0m"
	@mkdir -p $(@D)
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(BINDIR)/$(PROGRAM): $(OBJ)
	@echo "\e[92;1mCOMPILING PROGRAM $(PROGRAM):\e[0m"
	@mkdir -p $(@D)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 
	@echo "\e[32;1mBUILD OK\e[0m"

debug:
	DEBUG=1 make clean print $(BINDIR)/$(PROGRAM)
	
print:
	@echo "\e[36;1mDIRECTORIES:\e[0m"
	@echo "SRCDIR = $(SRCDIR)"
	@echo "BINDIR = $(BINDIR)"
	@echo "OBJDIR = $(OBJDIR)"
	@echo "MODDIR = $(MODDIR)"
	@echo "\e[34;1mFILES:\e[0m"
	@echo "SRC = $(SRC)"
	@echo "OBJ = $(OBJ)"
	@echo "\e[31;1mCOMPILER FLAGS:\e[0m"
	@echo "FFLAGS = $(FFLAGS)"
	@echo "\e[35;1mPROGRAM:\e[0m"
	@echo "PROGRAM = $(PROGRAM)"

clean:
	@echo "\e[34;1mREMOVING FOLDERS:\e[0m"
	rm -fr $(BINDIR) $(OBJDIR) $(MODDIR)

