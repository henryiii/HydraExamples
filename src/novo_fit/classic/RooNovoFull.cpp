#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>

#include <RooRealVar.h>
#include <RooNovosibirsk.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooMinuit.h>

// Compile with:
// `root-config --cxx --cflags --glibs` -lRooFit -lRooFitCore RooNovoFull.cpp 

int main() {
    // Main variable
    RooRealVar x{"x", "x", -10, 1};

    // Paramaters
    RooRealVar mean{"mean","mean", .4, -10, 10};
    RooRealVar sigma{"sigma","sigma", .6, 0, 1};
    RooRealVar tail{"tail","tail", 1.1, 0, 3};

    // PDF
    RooNovosibirsk novo{"novo", "novo", x, mean, sigma, tail};

    // Generating data
    RooDataSet *ds = novo.generate(RooArgSet(x), 1000);

    RooAbsReal *nll = novo.createNLL(*ds);

    RooMinuit m(*nll);
    m.setVerbose(kTRUE);
    m.migrad();


    // Plotting
    TApplication app{"myapp", 0, 0};
    TCanvas c{"c","canvas", 800, 600};
    RooPlot* xframe2 = x.frame(RooFit::Title("Novo Fit PDF"));
    ds->plotOn(xframe2);
    novo.plotOn(xframe2);
    xframe2->Draw();
    c.Draw();
    app.Run();

    return 0;
}

