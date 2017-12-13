from ROOT import gROOT, TCanvas, TH1, TLatex, TStyle
import FileIO

inputFile = TFile.Open('Efficiency.root')

effPlots = []
for key, obj in FileIO.FileIterator(inputFile):
    if 'TH1' in obj.ClassName():
        effPlots.append(obj)

canvas = TCanvas('canvas','',800,800)
effPlots[0].Draw()
