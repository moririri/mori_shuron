all: mori_shuron solv

mori_shuron: ./src/mori_shuron.c
	gcc -o mori_shuron ./src/mori_shuron.c -std=c99 -lm -O3 -march=native

solv: ./src/solv.c
	gcc -o solv ./src/solv.c -std=c99 -lm

