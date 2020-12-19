all: mori_shuron mori_shuron_old mori_shuron_not_random mori_shuron_normal solv plot_check 

mori_shuron: ./src/mori_shuron.c
	gcc -o mori_shuron ./src/mori_shuron.c -std=c99 -lm -O3 -march=native -DENABLE_RANDOM=1

mori_shuron_old: ./src/mori_shuron_old.c
	gcc -o mori_shuron_old ./src/mori_shuron_old.c -std=c99 -lm -O3 -march=native

mori_shuron_not_random: ./src/mori_shuron.c
	gcc -o mori_shuron_not_random ./src/mori_shuron.c -std=c99 -lm -O3 -march=native -DENABLE_RANDOM=0

mori_shuron_normal: ./src/mori_shuron.c
	gcc -o mori_shuron_normal ./src/mori_shuron.c -std=c99 -lm -O3 -march=native

solv: ./src/solv.c
	gcc -o solv ./src/solv.c -std=c99 -lm

plot_check: ./src/plot_check.c
	gcc -o plot_check ./src/plot_check.c -std=c99 -lm

sample_str: ./test/sample_str.c
	gcc -o sample_str ./test/sample_str.c

clean:
	rm -rf mori_shuron mori_shuron_not_random mori_shuron_normal solv plot_check mori_shuron_old sample_str


