

#
CC := g++
CFLAGS := -std=c++17 -Wall -g 


#
default: all


#
menu:
	@echo
	@echo "*** Project 4 Main Menu ***"
	@echo
	@echo "make menu            ==> This menu"
	@echo
	@echo "make all             ==> Run all targets"
	@echo "make test            ==> Run tests"
	@echo
	@echo "make test            ==> Build the maxdefense test"
	@echo "make maxdefense      ==> Build maxdefense"
	@echo


#
all: maxdefense test

test: maxdefense_test 
	./maxdefense_test

maxdefense_test: maxdefense.hh rubrictest.hh maxdefense_test.cc
	$(CC) $(CFLAGS) maxdefense_test.cc -o $@

maxdefense: maxdefense.hh timer.hh maxdefense_main.cc
	$(CC) $(CFLAGS) maxdefense_main.cc -o experiment

clean:
	-rm -f experiment maxdefense maxdefense_test


