main: build_main run_main clean_main  
	@echo main done
build_main: main.c
	gcc -g -Wall -Wextra ./main.c -lm -o main

run_main:
	./main $(IN_DIR) $(MESH_DIR)

clean_main:
	rm main

##########################################################
mesh_generate: build_mesh_generate run_mesh_generate clean_mesh_generate  
	@echo mesh_generate done
build_mesh_generate: mesh_generate.c
	gcc -g -Wall -Wextra ./mesh_generate.c -lm -o mesh_generate

run_mesh_generate:
	./mesh_generate

clean_mesh_generate:
	rm mesh_generate

test:
	gcc -Wall -Wextra ./test.c -lm -o test
	./test
	rm test