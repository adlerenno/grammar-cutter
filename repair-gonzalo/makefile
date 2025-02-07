
all: repair despair irepair idespair toint repair.tgz

repair.tgz: *.c *.h makefile COPYRIGHT README bal/*.c bal/*.h bal/makefile bal/README large/*.c large/*.h large/makefile large/README large/bal/*.c large/bal/*.h large/bal/makefile large/bal/README
	tar -czvf repair.tgz *.c *.h makefile COPYRIGHT README bal/*.c bal/*.h bal/makefile bal/README large/*.c large/*.h large/makefile large/README large/bal/*.c large/bal/*.h large/bal/makefile large/bal/README

repair: repair.o array.o hash.o heap.o records.o basics.o makefile
	gcc -O9 -o repair repair.o array.o hash.o heap.o records.o basics.o

repair.o: repair.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c repair.c

despair: despair.o basics.o makefile
	gcc -O9 -o despair despair.o basics.o

despair.o: despair.c basics.h makefile
	gcc -O9 -c despair.c

irepair: irepair.o array.o hash.o heap.o records.o basics.o makefile
	gcc -O9 -o irepair irepair.o array.o hash.o heap.o records.o basics.o

irepair.o: irepair.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c irepair.c

idespair: idespair.o basics.o makefile
	gcc -O9 -o idespair idespair.o basics.o

idespair.o: idespair.c basics.h makefile
	gcc -O9 -c idespair.c

array.o: array.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c array.c

hash.o: hash.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c hash.c

heap.o: heap.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c heap.c

records.o: records.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c records.c

basics.o: basics.c basics.h makefile
	gcc -O9 -c basics.c

toint: toint.c makefile
	gcc -O9 -o toint toint.c
