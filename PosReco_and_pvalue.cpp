using namespace std;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

pair<double,bool> pvalue(double N, double bckg, int nsigmas){
  bool Is_n_sigma = false;
  TF1 *poisson = new TF1("poisson","TMath::Poisson(x,[0])",0,5000);
  poisson->SetParameter(0,bckg);
  double pvalue =  poisson->Integral(N,5000);
  double n_sigma =  poisson->Integral(0,nsigmas*TMath::Sqrt(bckg));
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

void PosReco_and_pvalue(){

  // Map: boat from x=50*1e3 to x=-50*1e3 in 1 min steps alors y=0, detectors at x=0, z=-1000, detectors spaced by 10 km each, from UK to Iceland (~800km so 80 modules). speed 50 km/h = 0.8333 km/min
  
  TGraphErrors* approachgraph = new TGraphErrors();
  TGraphErrors* approachgraphBG = new TGraphErrors();
  TGraphErrors* approachgraphPVal = new TGraphErrors();
  
  int Nmodules = 80;
  double Power = 150*5.;//in MW    
  double ModuleDy = 10*1e3; // in km. Distance between modules
  double TotEventCount = 0.;
  double BGEventCount = 1.;
  double time = 0;
  double zpos = 1e3; //1000 m overburden in the ocean (water)
  double dt = 1.; //min
  double speed = (Power/3.)*1e3 /60.; //m/min
  double dx = speed*dt; //m
  double MassTarget = 100; //in ktons
  double startx = -100*1e3;
  double endx = 100*1e3;
  
  vector<int> time_sigmas;
  int time_for_one_event = -1;
  TPolyMarker3D *modules_coordinates = new TPolyMarker3D(0);
  TPolyMarker3D *sub_coordinates = new TPolyMarker3D(0);
  for(double xpos = startx; xpos<endx; xpos+=dx){ // move in 1 min steps through the water
    sub_coordinates->SetNextPoint(xpos*1e-3,0,zpos*1e-3);
    for(int i = 0; i<Nmodules;++i){ // Sum up many detectors...
      double ydet = 0;
      if(Nmodules % 2 == 0){
        ydet = (Nmodules*ModuleDy)/2.-i*ModuleDy - 5*1e3;
      }
      else{
        ydet = (Nmodules*ModuleDy)/2.-i*ModuleDy;
      }
      if(xpos == startx){
        modules_coordinates->SetNextPoint(0,ydet*1e-3,0);
      }
      double dist = TMath::Sqrt(xpos*xpos+ydet*ydet+zpos*zpos);
      TotEventCount+= 1*ExpNuFlux(dist, MassTarget, Power)*Pee2Flav(3.6,dist);
      BGEventCount+= ExpBGFlux(MassTarget);
      if(TotEventCount+BGEventCount-1 > 1 and time_for_one_event == -1){
        time_for_one_event = time;
      }
    }
    approachgraph->SetPoint(approachgraph->GetN(),time,TotEventCount);
    approachgraphBG->SetPoint(approachgraphBG->GetN(),time,BGEventCount);
    pair<double,bool> pvalue_and_sigma = pvalue(TotEventCount+BGEventCount,BGEventCount,time_sigmas.size()+1);
    approachgraphPVal->SetPoint(approachgraphPVal->GetN(),time,pvalue_and_sigma.first);
    if(pvalue_and_sigma.second){
      time_sigmas.push_back(time);
    }
    time += dt;
  }
  int sigma_lim = 5;
  vector<TLine*> sigma_lines;
  cout << "First event seen after " << time_for_one_event << " minutes. Using that as t0" << endl;
  if(time_for_one_event > 0){
    for(int i = 0; i<approachgraph->GetN(); i++){
      double t,n;
      approachgraphPVal->GetPoint(i,t,n);
      approachgraphPVal->SetPoint(i,t+time_for_one_event,n);
    }
    for(int i = 0; i<time_sigmas.size(); i++){
      time_sigmas[i] = time_sigmas[i] + time_for_one_event;
    }
  }
  if(time_sigmas.size() < 5){
    sigma_lim = time_sigmas.size();
    cout << "Could not find sub at more than " << time_sigmas.size() << " sigmas after " << time << " minutes." << endl;
  }
  
 
  vector<int> colors = {632,800,600,417};
  approachgraph->SetMarkerColor(kRed);
  approachgraph->SetLineColor(kRed);
  approachgraph->SetTitle("Signal events;t(min)");
  approachgraph->SetMarkerStyle(2);
  approachgraphBG->SetMarkerColor(kBlue);
  approachgraphBG->SetLineColor(kBlue);
  approachgraphBG->SetTitle("BG expected events;t(min)");
  approachgraphBG->SetMarkerStyle(2);
  approachgraphPVal->SetMarkerColor(kGreen+1);
  approachgraphPVal->SetLineColor(kGreen+1);
  approachgraphPVal->SetTitle("Pvalue;t(min)");
  approachgraphPVal->SetMarkerStyle(2);
  TLegend *leg1 = new TLegend(0.1,0.7,0.48,0.9);
  leg1->AddEntry(approachgraph);
  leg1->AddEntry(approachgraphBG);
  TLegend *leg2 = new TLegend(0.55,0.75,0.9,0.9);
  leg2->AddEntry(approachgraphPVal);
  TCanvas* c2 = new TCanvas("Moving_submarine","moving_submarine",1500,750);
  c2->Divide(2,1);
  c2->cd(1);
  approachgraph->Draw("AP"); 
  approachgraphBG->Draw("Psame");
  c2->Update();
  TLine *line_firstevent = new TLine(time_for_one_event,gPad->GetUymin(),time_for_one_event,gPad->GetUymax());
  line_firstevent->SetLineColor(kGreen+1);
  leg1->AddEntry(line_firstevent,"First event from submarine","l");
  line_firstevent->Draw();
  leg1->Draw();
  c2->cd(2);
  approachgraphPVal->Draw("AP");
  int n = approachgraphPVal->GetN();
  c2->Update();
  if(sigma_lim > 1){
    for(int i=1;i<sigma_lim;i++){
      cout << i+1 <<" sigmas reached after: " << time_sigmas[i] << " minutes" << endl;
      TLine *line = new TLine(time_sigmas[i],gPad->GetUymin(),time_sigmas[i],gPad->GetUymax());
      line->SetLineColor(colors[i-1]);
      leg2->AddEntry(line,TString::Format("%i sigmas", i+1),"l");
      line->Draw();
    }
  }
  leg2->Draw();
  c2->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Power,3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+".png").data());
  c2->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Power,3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+".root").data());
  
  
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
  c1->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Power,3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_coordinates.png").data());
  c1->SaveAs(string("modules"+to_string(Nmodules)+"_power"+to_string_with_precision(Power,3)+"_km"+to_string_with_precision(endx/1e3,3)+"_mass"+to_string_with_precision(MassTarget,3)+"_coordinates.root").data());
  return;
}


