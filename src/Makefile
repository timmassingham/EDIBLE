CC = gcc
CFLAGS = -O3 -funroll-loops -DNDEBUG -fomit-frame-pointer -std=gnu99 -DHAVE_CONFIG_H
#CFLAGS = -static -O3 -funroll-loops -DNDEBUG -std=gnu99 -Wall
#CFLAGS =  -DNDEBUG -pg  -Wall -std=gnu99
#CFLAGS = -std=c99 -Wall -g -DWARNINGS -std=gnu99 -DHAVE_CONFIG_H
#CFLAGS = -O -funroll-loops  -std=gnu99 -DNDEBUG
LD = ld
LDFLAGS =  -lm -lc
INCFLAGS = 
#LDFLAGS = -lm
objects =  Meschach/meschach.a edible.o read.o utility.o tree.o random.o llh.o partial.o matrix.o options.o new_models.o gtr.o 


edible: $(objects)
	$(CC)  $(INCFLAGS) -o ../bin/$@ $(LDFLAGS) $(objects)

pcalc: pcalc/pcalc.o
	$(CC)  $(INCFLAGS) -o ../bin/$@ $< $(LDFLAGS)
	
.c.o:
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<
.f.o:
	g77 -O3 -o $@ -c $<
clean:
	rm -f *~ *.o  *.obj  *.exe *.lib  pcalc/pcalc.o && cd Meschach && make clean && rm -f meschach.a

Meschach/meschach.a:
	cd Meschach && ./configure && make basic

archive:
	make clean && cd ../.. && tar zcf edible.tgz --exclude=.git --exclude=bin --exclude=archive edible
