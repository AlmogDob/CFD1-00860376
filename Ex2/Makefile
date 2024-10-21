solver: build_solver run_solver clean_solver  
	@echo solver done
build_solver: solver.c
	gcc -g -Wall -Wextra ./solver.c -lm -o solver

run_solver:
	./solver $(IN_FILL) $(MESH_FILL) $(AUTO_DIR) $(NUM)

clean_solver:
	rm solver

##########################################################
mesh_generate: build_mesh_generate run_mesh_generate clean_mesh_generate  
	@echo mesh_generate done
build_mesh_generate: mesh_generate.c
	gcc -g -Wall -Wextra ./mesh_generate.c -lm -o mesh_generate

run_mesh_generate:
	./mesh_generate

clean_mesh_generate:
	rm mesh_generate

##########################################################
test:
	gcc -Wall -Wextra ./test.c -lm -o test
	./test
	rm test