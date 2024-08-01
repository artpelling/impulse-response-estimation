#!/usr/bin/env make

.PHONY: clean cleanlatex cleanplots document latex plots venv

clean : cleanlatex cleanplots cleanvenv

cleanlatex :
	-rm build.pdf
	-rm document.pdf
	-rm -rf ./latex/build

cleanplots :
	-rm ./latex/figures/*.pdf

cleanvenv :
	-rm -rf ./code/venv

pyfar :
	@$(MAKE) code/pyfar
	@echo 'Checking for upstream changes ...'
	@cd code/pyfar && git pull

code/pyfar :
	git clone -b toeplitz-deconvolution git@github.com:artpelling/pyfar.git code/pyfar

code/venv : code/requirements.txt
	python3 -m venv code/venv
	./code/venv/bin/python -m pip install --upgrade pip
	./code/venv/bin/python -m pip install -r ./code/requirements.txt
	@touch code/venv

venv : pyfar
	@$(MAKE) code/venv

document : cleanlatex plots
	cd latex && latexmk -pdf -f -view=none -interaction=nonstopmode -synctex=1 -shell-escape -bibtex -output-directory="build" "root.tex"
	cp latex/build/root.pdf document.pdf

latex : cleanlatex
	mkdir -p latex/build && touch latex/build/root.pdf && ln -s latex/build/root.pdf build.pdf
	cd latex && latexmk -pvc -pdf -view=none -interaction=nonstopmode -synctex=1 -shell-escape -bibtex -output-directory="build" "root.tex"

plots : venv
	./code/venv/bin/python ./code/plots.py
