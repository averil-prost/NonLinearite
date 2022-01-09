############### Makefile (à laisser à la racine du projet)
# compile séparément tous les .cpp situés dans des dossiers de SRCDIR
# puis fait l'édition des liens avec le main
# les .o sont placés dans OBJDIR, l'exécutable dans BINDIR

CC = g++
OPT = -g

# à tester pour mac :
# CC = cl
# OPT = 

main = control
output = control

# ---------- répertoires
SRCDIR = src/
OBJDIR = obj/
BINDIR = bin/
DATADIR = data/

# ---------- librairies
LIB = 

# ---------- dépendances
### hypothèse : seul le main est à la racine, tout le reste dans des sous-répertoires
### de $(SRCDIR)
SOURCES = $(wildcard $(SRCDIR)*/*.cpp $(SRCDIR)*/*/*.cpp $(SRCDIR)*/*/*/*.cpp $(SRCDIR)*/*/*/*/*.cpp)
# $(info sources : $(SOURCES))
# DEP doit être de la forme $(OBJDIR)<même racine que le .cpp, mais dans OBJ>nomfile.obj
DEP = $(patsubst $(SRCDIR)%,$(OBJDIR)%.o,$(basename $(SOURCES)))
# $(info dep : $(DEP))
# chemins dans lesquels make doit chercher les .cpp
VPATH = $(dir $(SOURCES))
# $(info vpath : $(VPATH))

# ---------- exécutables
$(main) : $(DEP) $(SRCDIR)$(main).cpp
	$(CC) $(OPT) $^ -o $(BINDIR)$(output) $(LIB)

# ---------- compilation des dépendances

$(OBJDIR)%.o : $(SRCDIR)%.cpp
	mkdir -p $(dir $@)
	$(CC) $(OPT) $^ -c -o $@

# ---------- nettoyage 
.PHONY : clean
clean :
	rm $(BINDIR)$(output) $(OBJDIR)*.o $(OBJDIR)*/*.o $(OBJDIR)*/*/*.o $(OBJDIR)*/*/*/*.o $(OBJDIR)*/*/*/*/*.o 

.PHONY : rmvtk
rmvtk : 
	rm $(DATADIR)vtk/*.vtk