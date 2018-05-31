////////////////Created by Philipp Kampmann and Philippe Cotte//////////////

#include <vector>

using namespace std;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

double Significance(double sigevents, double bgevents){
  double totevents = sigevents+bgevents;
  return TMath::Sqrt(2*(totevents*TMath::Log(totevents/bgevents)+bgevents-totevents));
}

pair<double,bool> pvalue(double N, double bckg, int nsigmas){
  bool Is_n_sigma = false;
  TF1 *poisson = new TF1("poisson","TMath::Poisson(x,[0])",0,5000);
  poisson->SetParameter(0,bckg);
  double pvalue =  1-poisson->Integral(0,N-1);
  double n_sigma =  poisson->Integral(0,bckg+nsigmas*TMath::Sqrt(bckg));
  if(1-pvalue > n_sigma){
    Is_n_sigma = true;
  }
  return make_pair(pvalue,Is_n_sigma);
}

double Pee2Flav(double energy, double dist){
  double sin22t12 = 0.85;
  double dm2 = 7.59e-5;
  double oscterm = TMath::Sin(1.27*dm2*dist/energy);
  return 1-sin22t12*oscterm*oscterm;
}

double ExpNuFlux(double dist, double MassTarget,  double Power){ //per min
  double constants = 6.257; 
  constants*=4; //Rough correction factor
  return constants*Power*MassTarget/dist/dist*60; //events/min
}

double ExpBGFlux(double MassTarget){
  return 1.e-5*MassTarget;
}

double speed(double power){
  return (power/3.)*1e3 /60.;
}

void PosReco_and_pvalue(){

  // Map: boat from x=50*1e3 to x=-50*1e3 in 1 min steps alors y=0, detectors at x=0, z=-1000, detectors spaced by 10 km each, from UK to Iceland (~800km so 80 modules). speed 50 km/h = 0.8333 km/min

  
  int Nmodules = 80;
  vector<double> Powers = {180,150,120,90,60,30};//in MW    ,120,90,60,30
  vector<int> GrColors = {632,600,417,401,800,433};
  double ModuleDy = 10*1e3; // in km. Distance between modules
  double zpos = 1e3; //1000 m overburden in the ocean (water)
  double dt = 10.; //min
  double MassTarget = 100; //in ktons
  double startx = -100*1e3;
  double endx = 100*1e3;
    
  vector <TGraphErrors*> approachgraph = {};
  TGraphErrors* approachgraphBG = new TGraphErrors();
  vector <TGraphErrors*> approachgrapSignificance = {};
  vector <TGraphErrors*> derivative_significance = {};
  
  int time_for_one_event = -1;
  TPolyMarker3D *modules_coordinates = new TPolyMarker3D(0);
  TPolyMarker3D *sub_coordinates = new TPolyMarker3D(0);
  bool first_iteration = true;
  TLegend *leg1 = new TLegend(0.7,0.1,0.9,0.3);
  
  for( int i = 0; i < Powers.size(); i++){
    double Power = Powers[i];
    double TotEventCount = 0.;
    double InstantEventCount = 0.;
    double BGEventCount = 1.;
    double InstantBGEventCount = 1.;
    double time = 0;
    double dx = speed(Power)*dt; //m
    approachgraph.push_back(new TGraphErrors());
    approachgrapSignificance.push_back(new TGraphErrors());
    derivative_significance.push_back(new TGraphErrors());
    for(double xpos = startx; xpos<endx; xpos+=dx){ // move in 1 min steps through the water
      if(first_iteration){
        sub_coordinates->SetNextPoint(xpos*1e-3,0,zpos*1e-3);
      }
      for(int i = 0; i<Nmodules;++i){ // Sum up many detectors...
        double ydet = 0;
        if(Nmodules % 2 == 0){
          ydet = (Nmodules*ModuleDy)/2.-i*ModuleDy - 5*1e3;
        }
        else{
          ydet = (Nmodules*ModuleDy)/2.-i*ModuleDy;
        }
        if(xpos == startx and first_iteration){
          modules_coordinates->SetNextPoint(0,ydet*1e-3,0);
        }
        double dist = TMath::Sqrt(xpos*xpos+ydet*ydet+zpos*zpos);
        InstantBGEventCount += ExpBGFlux(MassTarget);
        InstantEventCount += dt*ExpNuFlux(dist, MassTarget, Power)*Pee2Flav(3.6,dist);
        TotEventCount += dt*ExpNuFlux(dist, MassTarget, Power)*Pee2Flav(3.6,dist);
        BGEventCount += ExpBGFlux(MassTarget);
        if(TotEventCount+BGEventCount-1 > 1 and time_for_one_event == -1){
          time_for_one_event = time;
        }
      }
      if(Powers.size() > 1){
        approachgraph.back()->SetPoint(approachgraph.back()->GetN(),time,TotEventCount);
        approachgraphBG->SetPoint(approachgraphBG->GetN(),time,BGEventCount-1);
      }
      else{
        approachgraph.back()->SetPoint(approachgraph.back()->GetN(),xpos,TotEventCount);
        approachgraphBG->SetPoint(approachgraphBG->GetN(),xpos,BGEventCount-1);
      }
      approachgrapSignificance.back()->SetPoint(approachgrapSignificance.back()->GetN(),xpos, Significance(TotEventCount,BGEventCount));
      derivative_significance.back()->SetPoint(derivative_significance.back()->GetN(),xpos,Significance(InstantEventCount, InstantBGEventCount));
      InstantEventCount = 0.;
      InstantBGEventCount = 1.;
      time += dt;
    }
    cout << Power << " MW:" << endl;
    cout << "  First event seen after " << time_for_one_event << " minutes. Using that as t0" << endl;
 
    approachgraph.back()->SetMarkerColor(GrColors[i]);
    approachgraph.back()->SetLineColor(GrColors[i]);
    approachgraph.back()->SetMarkerStyle(2);
    approachgraph.back()->SetLineWidth(5);
    approachgrapSignificance.back()->SetMarkerColor(kGreen+i);
    approachgrapSignificance.back()->SetLineColor(kGreen+i);
    approachgrapSignificance.back()->SetMarkerStyle(2);
    approachgrapSignificance.back()->SetLineWidth(5);
    derivative_significance.back()->SetMarkerColor(kGreen+i);
    derivative_significance.back()->SetLineColor(kGreen+i);
    derivative_significance.back()->SetMarkerStyle(2);
    derivative_significance.back()->SetLineWidth(5);
    approachgraph.back()->SetTitle(to_string_with_precision(Power,3).data());
    if(first_iteration){
      approachgrapSignificance.back()->SetTitle("significance;x(km);# #sigma");
      derivative_significance.back()->SetTitle(string(";x(km);# #sigma / " + to_string_with_precision(dt,2) + " minutes").data());
      approachgraphBG->SetMarkerColor(kBlack);
      approachgraphBG->SetLineColor(kBlack);
      approachgraphBG->SetTitle("BG expected events");
      approachgraphBG->SetMarkerStyle(2);
      leg1->AddEntry(approachgraphBG);
    }
    leg1->AddEntry(approachgraph.back());
    first_iteration = false;
  }
  TMultiGraph *multi_approach = new TMultiGraph();
  for( int i = 0; i < approachgraph.size(); i++){
    multi_approach->Add(approachgraph[i]);
  }
  multi_approach->SetTitle("Signal events for different reactor power;x(km);events");
  
  TCanvas* events = new TCanvas("Moving_submarine_events","moving_submarine_events",1500,750);
  multi_approach->Draw("AC");
  approachgraphBG->Draw("same");
  leg1->Draw();
  events->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_events.png").data());
  
  TCanvas* significance_can = new TCanvas("Moving_submarine_pvalue","moving_submarine_pvalue",1500,750);
  approachgrapSignificance[0]->Draw("AC3");
  significance_can->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_pvalue.png").data());
  //significance_can->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_pvalue.root").data());
  
  TCanvas* significance_derivative_can = new TCanvas("significance_derivative_can","significance_derivative_can",1500,750);
  derivative_significance[0]->Draw("AC3");
  significance_derivative_can->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_timewindow.png").data());
  //significance_derivative_can->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_timewindow.root").data());
  
  
  TCanvas* c1 = new TCanvas("Sub_path","Sub_path",1500,750);
  TView3D *view = (TView3D*)TView::CreateView(1);
  TAxis3D *axis = new TAxis3D();
  axis->GetXaxis()->SetRangeUser(TMath::Min(-1.1*50,-1.1*Nmodules*ModuleDy*1e-3),TMath::Max(1.1*50,1.1*Nmodules*ModuleDy*1e-3));
  axis->GetYaxis()->SetRangeUser(TMath::Min(-1.1*50,-1.1*Nmodules*ModuleDy*1e-3),TMath::Max(1.1*50,1.1*Nmodules*ModuleDy*1e-3));
  axis->GetZaxis()->SetRangeUser(0,10*zpos*1e-3);
  axis->Draw();
  view->SetRange(TMath::Min(-1.1*50,-1.1*Nmodules*ModuleDy*1e-3),TMath::Min(-1.1*50,-1.1*Nmodules*ModuleDy*1e-3),0,TMath::Max(1.1*50,1.1*Nmodules*ModuleDy*1e-3),TMath::Max(1.1*50,1.1*Nmodules*ModuleDy*1e-3),10*zpos*1e-3);
  modules_coordinates->SetMarkerColor(kBlue);
  modules_coordinates->SetMarkerStyle(2);
  modules_coordinates->SetMarkerSize(2);
  sub_coordinates->SetMarkerColor(kRed);
  sub_coordinates->SetMarkerStyle(2);
  modules_coordinates->Draw();
  sub_coordinates->Draw();
  //c1->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_coordinates.png").data());
  //c1->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Powers[0],3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_coordinates.root").data());
  return;
}


