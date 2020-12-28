all: mori_shuron mori_shuron_clang mori_shuron_old mori_shuron_not_random mori_shuron_normal solv plot_check 

CC := clang
GCC := gcc
RM := rm

mori_shuron: ./src/mori_shuron.c
	$(GCC) -o mori_shuron ./src/mori_shuron.c -std=c99 -lm -O3 -mtune=native -march=native -mfpmath=both -DENABLE_RANDOM=1

mori_shuron_clang: ./src/mori_shuron_clang.c
	$(CC) -o mori_shuron_clang ./src/mori_shuron_clang.c -std=c99 -lm -O3 -mtune=native -march=native -DENABLE_RANDOM=1

mori_shuron_old: ./src/mori_shuron_old.c
	$(GCC) -o mori_shuron_old ./src/mori_shuron_old.c -std=c99 -lm -O3 -mtune=native -march=native -mfpmath=both 

mori_shuron_not_random: ./src/mori_shuron.c
	$(GCC) -o mori_shuron_not_random ./src/mori_shuron.c -std=c99 -lm -O3 -mtune=native -march=native -mfpmath=both -DENABLE_RANDOM=0

mori_shuron_normal: ./src/mori_shuron.c
	$(GCC) -o mori_shuron_normal ./src/mori_shuron.c -std=c99 -lm -O3 -mtune=native -march=native -mfpmath=both 

solv: ./src/solv.c
	$(GCC) -o solv ./src/solv.c -std=c99 -lm

plot_check: ./src/plot_check.c
	$(GCC) -o plot_check ./src/plot_check.c -std=c99 -lm

clean:
	$(RM) -rf mori_shuron mori_shuron_clang mori_shuron_not_random mori_shuron_normal solv plot_check mori_shuron_old


