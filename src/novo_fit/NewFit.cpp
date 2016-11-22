
#include <iostream>
#include <assert.h>
#include <time.h>
#include <string>
#include <vector>
#include <array>
#include <chrono>

//this lib
#include <hydra/Types.h>
#include <hydra/Vector4R.h>
#include <hydra/Containers.h>
#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Random.h>
#include <hydra/VegasState.h>
#include <hydra/Vegas.h>
#include <hydra/LogLikelihoodFCN.h>
#include <hydra/PointVector.h>
#include <hydra/Parameter.h>
#include <hydra/UserParameters.h>
#include <hydra/Pdf.h>
#include <hydra/AddPdf.h>

#include <thrust/transform.h>

#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

//root
#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TString.h>
#include <TStyle.h>

#include <Gauss.h>
#include <Novosibirsk.h>

using namespace std;
using namespace ROOT::Minuit2;
using namespace hydra;

struct {
    __host__ __device__
    GReal_t operator()(thrust::tuple<GReal_t> &val) const {
        return thrust::get<0>(val);
    }
} First;

GInt_t main(int argv, char** argc) {

	size_t  nentries           = 1e6;
	size_t  iterations         = 50000;
	GReal_t tolerance          = 1.0;
	GBool_t use_comb_minimizer = false;

	//Print::SetLevel(0);
	ROOT::Minuit2::MnPrint::SetLevel(3);
	//----------------------------------------------

	//Generator with current time count as seed.
	size_t seed = 0; //std::chrono::system_clock::now().time_since_epoch().count();
	Random<thrust::random::default_random_engine> Generator( seed  );

	//-------------------------------------------
	//range of the analysis
	std::array<GReal_t, 1>  min   ={ -10};
	std::array<GReal_t, 1>  max   ={ 1};

	//------------------------------------
	//parameters names
	std::string Mean("Mean");
	std::string Sigma("Sigma");
	std::string Tail("Tail");

	//fit paremeters
	//----------------------------------------------------------------------

	// 1) using named parameter idiom
	Parameter  mean_p  = Parameter::Create().Name("Mean").Value(.4) .Error(0.001).Limits(-10., 10.);
	Parameter  sigma_p = Parameter::Create().Name("Sigma").Value(.6).Error(0.001).Limits(0., 1.);
	Parameter  tail_p  = Parameter::Create().Name("Tail").Value(1.1).Error(0.001).Limits(0., 3.);

	//----------------------------------------------------------------------
    // registry the parameters
	UserParameters upar;

    upar.AddParameter(&mean_p);
    upar.AddParameter(&sigma_p);
    upar.AddParameter(&tail_p);

    //check all is fine
    upar.PrintParameters();

    //----------------------------------------------------------------------
	// create functor
    pdfs::Novosibirsk Novosibirsk{mean_p, sigma_p, tail_p, 0};

	//----------------------------------------------------------------------
	//get integration
    //Vegas state hold the resources for performing the integration
    VegasState<1> state = VegasState<1>( min, max); // nota bene: the same range of the analisys
	state.SetVerbose(-1);
	state.SetAlpha(1.75);
	state.SetIterations(5);
	state.SetUseRelativeError(1);
	state.SetMaxError(1e-3);

    //5,000 calls (fast convergence and precise result)
	Vegas<1> vegas(state, 10000);

	auto Novosibirsk_PDF   = make_pdf(Novosibirsk, &vegas);
	Novosibirsk_PDF.PrintRegisteredParameters();

	//----------------------------------------------------------------------
	//integrate with the current parameters just to test
	vegas.Integrate(Novosibirsk_PDF);
	cout << ">>> Novosibirsk intetgral prior fit "<< endl;
	cout << "Result: " << vegas.GetResult() << " +/- "
		 << vegas.GetAbsError() << " Chi2: "<< vegas.GetState().GetChiSquare() << endl;

	Novosibirsk_PDF.PrintRegisteredParameters();


	//--------------------------------------------------------------------
	//Generate data on the device with the original parameters
	//
    PointVector<device, GReal_t, 1> data_d(nentries);
    Generator.Sample<device>(Novosibirsk_PDF, min, max, data_d, nentries);


	//---------------------------
	//get data from device and fill histogram
	PointVector<host> data_h(data_d);

	TH1D hist_novo("novo", "", 1000, min[0], max[0]);

	for(auto point: data_h )
		hist_novo.Fill(point.GetCoordinate(0));

	//-------------------------------------------------
    //minimization

	//get the FCN
	auto modelFCN = make_loglikehood_fcn(Novosibirsk_PDF, data_d.begin(), data_d.end() );

    //print minuit parameters before the fit
	std::cout << upar << endl;

    //minimization strategy
	MnStrategy strategy(2);

	// create Migrad minimizer
	MnMigrad migrad(modelFCN, upar.GetState() ,  strategy);

	// create Minimize minimizer
	MnMinimize minimize(modelFCN,upar.GetState() ,  strategy);
	FunctionMinimum *minimum=0;

	// ... Minimize and profile the time
	auto start = std::chrono::high_resolution_clock::now();

	if(use_comb_minimizer){
		 minimum = new FunctionMinimum(minimize(iterations, tolerance));
	}
	else{
		 minimum = new FunctionMinimum(migrad(iterations, tolerance));
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;

	// output
	std::cout<<"minimum: "<<*minimum<<std::endl;

	//time
	std::cout << "-----------------------------------------"<<std::endl;
	std::cout << "| Time (ms) ="<< elapsed.count() <<std::endl;
	std::cout << "-----------------------------------------"<<std::endl;

	//------------------------------------------------------
	//Sampling the fitted model
        //Set the function with the fitted parameters
	Novosibirsk_PDF.SetParameters(minimum->UserParameters().Params());
	Novosibirsk_PDF.PrintRegisteredParameters();

	//sample fit function on the host nentries trials
	PointVector<host, GReal_t, 1> data2_h(0);
	Generator.SetSeed(std::chrono::system_clock::now().time_since_epoch().count()+1);//+1 because all can run very fast sometimes
	Generator.Sample(Novosibirsk_PDF, min, max, data2_h, nentries );

	TH1D hist_novo_fit("novo_fit", "", 100, min[0], max[0]);

	// histogram it for graphics representation
	for(auto point: data2_h )
			hist_novo_fit.Fill(point.GetCoordinate(0));

	TH1D hist_novo_plot("novo_plot", "", 100,  min[0], max[0]);
	for (size_t i=0 ; i<=100 ; i++) {
		GReal_t x = hist_novo_plot.GetBinCenter(i);
		hist_novo_plot.SetBinContent(i, Novosibirsk_PDF(&x) );
	}

	//scale
	hist_novo_fit.Scale(hist_novo.Integral()/hist_novo_fit.Integral() );
	hist_novo_plot.Scale(hist_novo.Integral()/hist_novo_plot.Integral() );



	/***********************************
	 * RootApp, Drawing...
	 ***********************************/

	TApplication *myapp=new TApplication("myapp",0,0);


	TCanvas canvas_gauss("canvas_gauss", "novo distribution", 500, 500);
	hist_novo.Draw("e0");
	hist_novo.SetMarkerSize(1);
	hist_novo.SetMarkerStyle(20);

	//sampled data after fit
	hist_novo_fit.Draw("barSAME");
	hist_novo_fit.SetLineColor(4);
	hist_novo_fit.SetFillColor(4);
	hist_novo_fit.SetFillStyle(3001);

	//original data
//	hist_novo.Draw("e0SAME");

	//plot
	hist_novo_plot.Draw("csame");
	hist_novo_plot.SetLineColor(2);
	hist_novo_plot.SetLineWidth(2);

	canvas_gauss.SaveAs("./fit.png");

	myapp->Run();

	delete minimum;

	return 0;


}



