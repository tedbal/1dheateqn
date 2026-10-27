main: src/main.cpp
	@echo "Building"
	@g++ src/main.cpp -I/lib/eigen-3.4.0/ -o bin/main 
	@echo "Done"




