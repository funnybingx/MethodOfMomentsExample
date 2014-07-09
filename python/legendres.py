from ROOT import RooRealVar, TCanvas, RooLegendre
from ROOT import RooFit as rf

x = RooRealVar("x", "x", -1, 1)
frame = x.frame()

for i in range(10):
    pi = RooLegendre("P_"+str(i), "P_"+str(i)+"(x)", x, i)
    pi.plotOn(frame, rf.LineColor(i+1))

canvas = TCanvas()
frame.SetYTitle("legendre polynomial")
frame.Draw()
canvas.SaveAs("plots/legendres.pdf")
