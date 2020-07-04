CC = g++
FLAGS = -Wall -O3
INCLUDES = -I./ -I../
LIB_NAME = librulkov.a

all: Base Rede Global Plasticidade SinalExterno 

Base:
	@$(CC) $(FLAGS) $(INCLUDES) -c src/Base.cpp -o src/Base.o
	@ar srv $(LIB_NAME) src/Base.o
	@rm src/Base.o

Rede: 
	@$(CC) $(FLAGS) $(INCLUDES) -c src/Rede.cpp -o src/Rede.o
	@ar srv $(LIB_NAME) src/Rede.o
	@rm src/Rede.o

Global: 
	@$(CC) $(FLAGS) $(INCLUDES) -c src/Global.cpp -o src/Global.o
	@ar srv $(LIB_NAME) src/Global.o
	@rm src/Global.o

Plasticidade: 
	@$(CC) $(FLAGS) $(INCLUDES) -c src/Plasticidade.cpp -o src/Plasticidade.o
	@ar srv $(LIB_NAME) src/Plasticidade.o
	@rm src/Plasticidade.o


SinalExterno:
	@$(CC) $(FLAGS) $(INCLUDES) -c src/SinalExterno.cpp -o src/SinalExterno.o
	@ar srv $(LIB_NAME) src/SinalExterno.o
	@rm src/SinalExterno.o


clean:
	@if [ -f $(LIB_NAME) ]; then rm $(LIB_NAME); fi

list:
	@cat makefile | grep : | grep -v cat
