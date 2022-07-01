all: main

main: hitchhiker.c helper.c
	@gcc -o main $^ -lisal