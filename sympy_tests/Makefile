currSHELL=/usr/bin/
all: test1.sympyX.png test2.sympyX.png test3.sympyX.png test4.sympyX.png test5.sympyX.png test6.sympyX.png
test1: test1.sympyX.png
test2: test2.sympyX.png
test3: test3.sympyX.png
test4: test4.sympyX.png
test5: test5.sympyX.png
test6: test6.sympyX.png
test1.sympyX.png: test1.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test1.sympy.tex
test1.sympy.tex: test1.py
	$(currSHELL)pypy test1.py
test2.sympyX.png: test2.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test2.sympy.tex
test2.sympy.tex: test2.py
	$(currSHELL)pypy test2.py
test3.sympyX.png: test3.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test3 sympy
test3.sympy.tex: test3.py
	$(currSHELL)pypy test3.py
test4.sympyX.png: test4.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test4.sympy.tex
test4.sympy.tex: test4.py
	$(currSHELL)pypy test4.py
test5.sympyX.png: test5.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test5 sympy
test5.sympy.tex: test5.py
	$(currSHELL)pypy test5.py
test6.sympyX.png: test6.sympy.tex
	$(currSHELL)bash postprocessing.sympy.tex.sh test6 sympy
test6.sympy.tex: test6.py
	$(currSHELL)pypy test6.py
.PHONY: all test1 test2 test3 test4 test5 test6 clean
clean:
	$(currSHELL)rm -ifv *.aux *.log *.ps *.dvi *.toc *.los *.pdf *.tex *.pyc *.png *.bak *.so *sympy.c 1
	$(currSHELL)rm -frv build
	$(currSHELL)rm -ifv common/*.pyc
	$(currSHELL)rm -ifv redirect2sympy.txt redirect2maxima.max redirect2sympy.done.txt
	$(currSHELL)bash ./checkhangups.sh
