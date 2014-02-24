all: controller fitness optim

optim:
	gcc -Wall abc_optim.c ccABC.o -o ../bin/abc_optim
fitness:
	gcc -Wall tscontrol_fitness.c -lm -o ../bin/tscontrol_fitness
controller:
	gcc -Wall tscontroller.c fuzzycontroller.c poly.c -o ../bin/tscontroller
