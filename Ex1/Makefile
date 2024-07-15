main: build run clean
	@echo main done
build: main.c
	gcc -g -Wall -Wextra ./main.c -lm -o main

run:
	./main $(IN_DIR) $(OUT_DIR)

clean:
	rm main

test:
	gcc -Wall -Wextra ./test.c -lm -o test
	./test
	rm test