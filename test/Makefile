
all: gcc_pop_para gcc_pop clang_pop_para clang_pop

gcc_pop_para:
	gcc -o gcc_pop_para sample_popcnt.c -O3 -mtune=native -march=native -fopenmp -lgomp -DPARA=1

gcc_pop:
	gcc -o gcc_pop sample_popcnt.c -O3 -mtune=native -march=native -fopenmp -lgomp

clang_pop_para:
	clang -o clang_pop_para sample_popcnt.c -O3 -mtune=native -march=native -fopenmp -DPARA=1

clang_pop:
	clang -o clang_pop sample_popcnt.c -O3 -mtune=native -march=native -fopenmp

clean:
	rm -rf gcc_pop_para gcc_pop clang_pop_para clang_pop
	

