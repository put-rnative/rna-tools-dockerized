all:
	@cd src; make; cd ..; if [ ! -d bin ]; then mkdir -v bin; fi; cp src/3dRNAscore/3dRNAscore-1.0 bin/3dRNAscore

clean:
	@rm -rf bin
	@cd src; make clean; cd ..


