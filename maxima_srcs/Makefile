currSHELL=/usr/bin/
all: BevelBasisInducedSuscTensor.maximaX.png MoMElemSqs.maximaX.png DielecSuscNLTensInduced.maximaX.png
MoMElemSqs: MoMElemSqs.maximaX.png
MoMElemSqs.maximaX.png: MoMElemSqs.maxima.tex
DielecSuscNLTensInduced: DielecSuscNLTensInduced.maximaX.png
Bevel: BevelBasisInducedSuscTensor.maximaX.png
	$(currSHELL)bash postprocessing.maxima.tex.sh MoMElemSqs maxima
DielecSuscNLTensInduced.maximaX.png: DielecSuscNLTensInduced.maxima.tex
	$(currSHELL)bash postprocessing.maxima.tex.sh DielecSuscNLTensInduced maxima
BevelBasisInducedSuscTensor.maximaX.png: BevelBasisInducedSuscTensor*.maxima.tex
	$(currSHELL)bash postprocessing.maxima.tex.sh BevelBasisInducedSuscTensor* maxima
.PHONY: all MoMElemSqs DielecSuscNLTensInduced
clean:
	$(currSHELL)rm -ifv *.aux *.log *.ps *.dvi *.toc *.los *.pdf *.tex *.png *.bak 1
	$(currSHELL)bash ./checkhangups.sh
