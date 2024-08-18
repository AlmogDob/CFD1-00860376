mesh_generate: build run clean
	@echo mesh_generate done
build: mesh_generate.c
	gcc -g -Wall -Wextra ./mesh_generate.c -lm -o mesh_generate

run:
	./mesh_generate

clean:
	rm mesh_generate

test:
	gcc -Wall -Wextra ./test.c -lm -o test
	./test
	rm test