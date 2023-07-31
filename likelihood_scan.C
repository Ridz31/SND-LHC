/*Evaluate the likelihood by scanning over the probability density functions of all the energies and Z position values*/
void likelihood_scan(int energy_true = 300, int Z_true = 291){

  std::string baseDir = "/eos/user/r/ridz01/";
  std::string pdfDir = baseDir+"Electron_Showers_pdfs/";

    bool debug = false;

    //int energy_test[10] = {50,100,200,300,400,500,750,1000,1250,1500};
    int energy_test[21] = {50,75,100,150,200,250,300,350,400,450,500,625,750,875,1000,1125,1250,1375,1500,1750,2000};
    int Z_test[8] = {285,287,289,291,293,295,297,299};

    int N_energy = 21;
    int N_Z= 8;

    std::vector<double> energy_bins;
    std::vector<double> Z_bins;

    for (int i_energy = 0; i_energy < N_energy; i_energy++) energy_bins.push_back(energy_test[i_energy]);
    energy_bins.push_back(energy_test[N_energy-1] + (energy_test[N_energy-1] - energy_test[N_energy-2]));

    for (int i_Z = 0; i_Z < N_Z; i_Z++) Z_bins.push_back(Z_test[i_Z]);
    Z_bins.push_back(Z_test[N_Z-1] + (Z_test[N_Z-1] - Z_test[N_Z-2]));

    // Lets get the probability histograms only once, at the start of the program, and store them in std::vectors
    std::vector<TH1D*> h_prob_ver_1;
    std::vector<TH1D*> h_prob_ver_2;
    std::vector<TH1D*> h_prob_ver_3;
    std::vector<TH1D*> h_prob_ver_4;
    std::vector<TH1D*> h_prob_ver_5;

    std::vector<TH1D*> h_prob_hor_1;
    std::vector<TH1D*> h_prob_hor_2;
    std::vector<TH1D*> h_prob_hor_3;
    std::vector<TH1D*> h_prob_hor_4;
    std::vector<TH1D*> h_prob_hor_5;

    std::vector<TH1D*> hist_SiPM_X_1;
    std::vector<TH1D*> hist_SiPM_X_2;
    std::vector<TH1D*> hist_SiPM_X_3;
    std::vector<TH1D*> hist_SiPM_X_4;
    std::vector<TH1D*> hist_SiPM_X_5;

    std::vector<TH1D*> hist_SiPM_Y_1;
    std::vector<TH1D*> hist_SiPM_Y_2;
    std::vector<TH1D*> hist_SiPM_Y_3;
    std::vector<TH1D*> hist_SiPM_Y_4;
    std::vector<TH1D*> hist_SiPM_Y_5;
    
    std::vector<TFile*> f_X;
    std::vector<TFile*> f_Y;

    std::vector<TFile*> f_X_SiPM;
    std::vector<TFile*> f_Y_SiPM;

    for ( int i=0; i<N_energy; i++){//energy
      for (int j=0; j<N_Z; j++){//Z
	f_X.push_back(new TFile((pdfDir+std::to_string(energy_test[i])+"GeV_"+std::to_string(Z_test[j])+"Z/Pdfs_X_hist.root").c_str(),"READ"));
	h_prob_ver_1.push_back((TH1D*)(f_X.back())->Get("Station 1;1"));
	h_prob_ver_2.push_back((TH1D*)(f_X.back())->Get("Station 2;1"));
	h_prob_ver_3.push_back((TH1D*)(f_X.back())->Get("Station 3;1"));
	h_prob_ver_4.push_back((TH1D*)(f_X.back())->Get("Station 4;1"));
	h_prob_ver_5.push_back((TH1D*)(f_X.back())->Get("Station 5;1"));
	
	f_Y.push_back(new TFile((pdfDir+std::to_string(energy_test[i])+"GeV_"+std::to_string(Z_test[j])+"Z/Pdfs_Y_hist.root").c_str(),"READ"));
	h_prob_hor_1.push_back((TH1D*)(f_Y.back())->Get("Station 1;1"));
	h_prob_hor_2.push_back((TH1D*)(f_Y.back())->Get("Station 2;1"));
	h_prob_hor_3.push_back((TH1D*)(f_Y.back())->Get("Station 3;1"));
	h_prob_hor_4.push_back((TH1D*)(f_Y.back())->Get("Station 4;1"));
	h_prob_hor_5.push_back((TH1D*)(f_Y.back())->Get("Station 5;1"));
	
    //Open SiPM files
	f_X_SiPM.push_back(new TFile((pdfDir+std::to_string(energy_test[i])+"GeV_"+std::to_string(Z_test[j])+"Z/SiPMs_X.root").c_str(),"READ"));
	hist_SiPM_X_1.push_back((TH1D*)(f_X_SiPM.back())->Get("Station 1 vertical;1"));
	hist_SiPM_X_2.push_back((TH1D*)(f_X_SiPM.back())->Get("Station 2 vertical;1"));
	hist_SiPM_X_3.push_back((TH1D*)(f_X_SiPM.back())->Get("Station 3 vertical;1"));
	hist_SiPM_X_4.push_back((TH1D*)(f_X_SiPM.back())->Get("Station 4 vertical;1"));
	hist_SiPM_X_5.push_back((TH1D*)(f_X_SiPM.back())->Get("Station 5 vertical;1"));
    
	f_Y_SiPM.push_back( new TFile((pdfDir+std::to_string(energy_test[i])+"GeV_"+std::to_string(Z_test[j])+"Z/SiPMs_Y.root").c_str(),"READ"));
	hist_SiPM_Y_1.push_back((TH1D*)(f_Y_SiPM.back())->Get("Station 1 horizontal;1"));
	hist_SiPM_Y_2.push_back((TH1D*)(f_Y_SiPM.back())->Get("Station 2 horizontal;1"));
	hist_SiPM_Y_3.push_back((TH1D*)(f_Y_SiPM.back())->Get("Station 3 horizontal;1"));
	hist_SiPM_Y_4.push_back((TH1D*)(f_Y_SiPM.back())->Get("Station 4 horizontal;1"));
	hist_SiPM_Y_5.push_back((TH1D*)(f_Y_SiPM.back())->Get("Station 5 horizontal;1"));
      }
    }
    // Similarly, let create the temporary histograms only once:
    // Create temporary histograms to loop through the hits and fill the corresponding h_hit histogram with the hit X or Y coordinate.
    TH1D * h_hit_ver_1 = new TH1D("XHits_Station_1", "h_X_Station_1;X(mm)", 100, -50, 0);
    TH1D * h_hit_ver_2 = new TH1D("XHits_Station_2", "h_X_Station_2;X(mm)", 100, -50, 0);
    TH1D * h_hit_ver_3 = new TH1D("XHits_Station_3", "h_X_Station_3;X(mm)", 100, -50, 0);
    TH1D * h_hit_ver_4 = new TH1D("XHits_Station_4", "h_X_Station_4;X(mm)", 100, -50, 0);
    TH1D * h_hit_ver_5 = new TH1D("XHits_Station_5", "h_X_Station_5;X(mm)", 100, -50, 0);

    TH1D * h_hit_hor_1 = new TH1D("YHits_Station_1", "h_Y_Station_1;Y(mm)", 120, 0, 60);
    TH1D * h_hit_hor_2 = new TH1D("YHits_Station_2", "h_Y_Station_2;Y(mm)", 120, 0, 60);
    TH1D * h_hit_hor_3 = new TH1D("YHits_Station_3", "h_Y_Station_3;Y(mm)", 120, 0, 60);
    TH1D * h_hit_hor_4 = new TH1D("YHits_Station_4", "h_Y_Station_4;Y(mm)", 120, 0, 60);
    TH1D * h_hit_hor_5 = new TH1D("YHits_Station_5", "h_Y_Station_5;Y(mm)", 120, 0, 60);

    // Open the simulation files
    // Directory for each event: /eos/user/r/ridz01/PG_electron_sim/{Energy}GeV/{Z}/sndLHC.PG_11-TGeant4_digCPP.root
    TFile * f_in = new TFile(("/eos/user/r/ridz01/PG_electron_sim/"+std::to_string(energy_true)+"GeV/"+std::to_string(Z_true)+"/sndLHC.PG_11-TGeant4_digCPP.root").c_str());

    // Now get the cbmsim TTree to access each event
    TTree * cbmsim = (TTree*)f_in->Get("cbmsim");

    // Scifi hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/sndScifiHit.h
    TClonesArray * DigiScifiHits = new TClonesArray("sndScifiHit");
    cbmsim->SetBranchAddress("Digi_ScifiHits", &DigiScifiHits);

    // Set up geometry
    TFile * fGeo = new TFile("/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V7_22November2022.root");
    TGeoManager * geo = (TGeoManager*) fGeo->Get("FAIRGeom");
    Scifi * ScifiDet = new Scifi("Scifi", kTRUE);

    // HACK HACK HACK! Set parameters by hand.
    ScifiDet->SetConfPar("Scifi/scifimat_length" , (Float_t) 39.0);
    ScifiDet->SetConfPar("Scifi/channel_width"   , (Float_t) (0.25/10));
    ScifiDet->SetConfPar("Scifi/epoxymat_z"      , (Float_t) 0.17);
    ScifiDet->SetConfPar("Scifi/nsipm_channels"  , (Int_t) 128);
    ScifiDet->SetConfPar("Scifi/nsipm_mat"       , (Int_t) 4);
    ScifiDet->SetConfPar("Scifi/nmats"           , (Int_t) 3);
    ScifiDet->SetConfPar("Scifi/sipm_edge"       , (Float_t) (0.17/10));
    ScifiDet->SetConfPar("Scifi/charr_width"     , (Float_t) (64 * 0.25/10));
    ScifiDet->SetConfPar("Scifi/charr_gap"       , (Float_t) (0.2/10));
    ScifiDet->SetConfPar("Scifi/sipm_diegap"     , (Float_t) (0.06/10));
    ScifiDet->SetConfPar("Scifi/firstChannelX"   , (Float_t) -19.528);
    // END of horrible horrible hack

    ScifiDet->SiPMmapping();

    // Simulation file of the true energy
    ofstream file;
    file.open((baseDir+"/Electron_Showers_Simulation_Plots/Likelihood_"+std::to_string(energy_true)+"GeV_"+std::to_string(Z_true)+"Z_100event_2Dscan.txt").c_str(), ios::app);

    // Initialising Histogram for the reconstruction of energy and Z
    //TH2D * hist_ZX = new TH2D(("Horizontal_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z").c_str(), "h_hor;z;x", 200, 280, 380, 200, -80, 80);
    TH1D * h_reco_energy = new TH1D("Reco_energy_test", "Energy_Reco; energy", N_energy, &energy_bins[0]);
    TH1D * h_reco_Z = new TH1D("Reco_Z_test", "Z_Reco; Z", N_Z, &Z_bins[0]);
    TH2D * h_reco = new TH2D("Reconstruction_Test", "Reco;energy;Z", N_energy, &energy_bins[0], N_Z, &Z_bins[0]);

    double epsilon = 1e-10; // to avoid zero values for log
    
    TVector3 a;
    TVector3 b;

    for (int i_event = 0; i_event < cbmsim->GetEntries(); i_event++){
      
        cout<<"Event Loop"<< i_event<< endl;
        cbmsim->GetEntry(i_event);
	    cout <<  DigiScifiHits->GetEntries() << endl;

	// Fill the temporary histograms only once per event
	// Reset the temporary histograms
	h_hit_ver_1->Reset();
	h_hit_ver_2->Reset();
	h_hit_ver_3->Reset();
	h_hit_ver_4->Reset();
	h_hit_ver_5->Reset();

	h_hit_hor_1->Reset();
	h_hit_hor_2->Reset();
	h_hit_hor_3->Reset();
	h_hit_hor_4->Reset();
	h_hit_hor_5->Reset();
                
	if (debug) cout<<"Temporary histograms reset" <<endl;

    std::vector<double> energy;
	std::vector<double> Z;
	std::vector<double> logL;

    // Loop through all the hits for each event and fill the histograms based on the Station 
    for (int scifihit = 0; scifihit<DigiScifiHits->GetEntries(); scifihit++){
	  if (debug) std::cout << "hit " << scifihit << " / " << DigiScifiHits->GetEntries() << std::endl;
            // Condition to check if it is a valid hit
            if (((sndScifiHit*)DigiScifiHits->At(scifihit))->isValid()){
                // Identify the station of the hit
                int hit_station = ((sndScifiHit*)DigiScifiHits->At(scifihit))->GetStation();
                if (debug) cout<< " Station "<< hit_station<<endl;
                // Find the SiPM id of the hit event 
                int id_hit = ((sndScifiHit*)DigiScifiHits->At(scifihit))->GetDetectorID();
                if (debug) cout<< "Hit id Hori  "<< id_hit<<endl;

                // To get the X and Y Coordinate
                ScifiDet->GetSiPMPosition(id_hit, a, b);
                // Fill the histograms based on the Station
                if (((sndScifiHit*)DigiScifiHits->At(scifihit))->isVertical()){
                    if (debug) cout<<"Vertical"<<endl;
                    if (hit_station==1){
                        h_hit_ver_1->Fill((a.X() + b.X())/2.);
                        if (debug) cout<<"Station 1 hist filled"<<endl;
                    }
                    if (hit_station==2){
                        h_hit_ver_2->Fill((a.X() + b.X())/2.);
                        if (debug) cout<<"Station 2 hist filled"<<endl;
                    }
                    if (hit_station==3){
                        h_hit_ver_3->Fill((a.X() + b.X())/2.);
                        if (debug) cout<<"Station 3 hist filled"<<endl;
                    }
                    if (hit_station==4){
                        h_hit_ver_4->Fill((a.X() + b.X())/2.);
                        if (debug) cout<<"Station 4 hist filled"<<endl;
                    }
                    if (hit_station==5){
                        h_hit_ver_5->Fill((a.X() + b.X())/2.);
                        if (debug) cout<<"Station 5 hist filled"<<endl;
                    }
                }
                else{
                    if (debug) cout<<"Horizontal"<<endl;
                    if (hit_station==1){
                        h_hit_hor_1->Fill((a.Y() + b.Y())/2.);
                        if (debug) cout<<"Station 1 hist filled"<<endl;
                    }
                    if (hit_station==2){
                        h_hit_hor_2->Fill((a.Y() + b.Y())/2.);
                        if (debug) cout<<"Station 2 hist filled"<<endl;
                    }
                    if (hit_station==3){
                        h_hit_hor_3->Fill((a.Y() + b.Y())/2.);
                        if (debug) cout<<"Station 3 hist filled"<<endl;
                    }
                    if (hit_station==4){
                        h_hit_hor_4->Fill((a.Y() + b.Y())/2.);
                        if (debug) cout<<"Station 4 hist filled"<<endl;
                    }
                    if (hit_station==5){
                        h_hit_hor_5->Fill((a.Y() + b.Y())/2.);
                        if (debug) cout<<"Station 5 hist filled"<<endl;
                    }
                }
            }// end of valid hit check
        }//end of hits loop

        // 2-D scan
        for ( int i=0; i<N_energy; i++){//energy
            for (int j=0; j<N_Z; j++){//Z

                //Initialize likelihood to 0 for each event, each energy and Z value
                double logLX = 0; double logLY =0;
		
                //Loop through all the bins for the X and Y
                for(int bin =1; bin<= h_prob_ver_1.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    if (debug) cout<< bin<< " "<< hist_SiPM_X_1.at(j+i*N_Z)->GetBinContent(bin)<< " " <<h_hit_ver_1->GetBinContent(bin)<< " " << h_prob_ver_1.at(j+i*N_Z)->GetBinContent(bin) <<endl;
                    double temp_log_X_1 = (hist_SiPM_X_1.at(j+i*N_Z)->GetBinContent(bin) - h_hit_ver_1->GetBinContent(bin))*log(1 - h_prob_ver_1.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_ver_1->GetBinContent(bin)*log(h_prob_ver_1.at(j+i*N_Z)->GetBinContent(bin) + epsilon);
                    if (debug) cout<< temp_log_X_1<< endl;
                    logLX +=  temp_log_X_1;
                }
                for(int bin =1; bin<= h_prob_ver_2.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    double temp_log_X_2 = (hist_SiPM_X_2.at(j+i*N_Z)->GetBinContent(bin) - h_hit_ver_2->GetBinContent(bin))*log(1 - h_prob_ver_2.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_ver_2->GetBinContent(bin)*log(h_prob_ver_2.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_X_2<< endl;
                    logLX +=  temp_log_X_2;
                }
                for(int bin =1; bin<= h_prob_ver_3.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    double temp_log_X_3 = (hist_SiPM_X_3.at(j+i*N_Z)->GetBinContent(bin) - h_hit_ver_3->GetBinContent(bin))*log(1 - h_prob_ver_3.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_ver_3->GetBinContent(bin)*log(h_prob_ver_3.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_X_3<< endl;
                    logLX +=  temp_log_X_3;
                }
                for(int bin =1; bin<= h_prob_ver_4.at(j+i*N_Z)->GetNbinsX(); bin++){
                    double temp_log_X_4 = (hist_SiPM_X_4.at(j+i*N_Z)->GetBinContent(bin) - h_hit_ver_4->GetBinContent(bin))*log(1 - h_prob_ver_4.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_ver_4->GetBinContent(bin)*log(h_prob_ver_4.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_X_4<< endl;
                    logLX +=  temp_log_X_4;
                }
                for(int bin =1; bin<= h_prob_ver_5.at(j+i*N_Z)->GetNbinsX(); bin++){
                    double temp_log_X_5 = (hist_SiPM_X_5.at(j+i*N_Z)->GetBinContent(bin) - h_hit_ver_5->GetBinContent(bin))*log(1 - h_prob_ver_5.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_ver_5->GetBinContent(bin)*log(h_prob_ver_5.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_X_5<< endl;
                    logLX +=  temp_log_X_5;
                }

                for(int bin =1; bin<= h_prob_hor_1.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    double temp_log_Y_1 = (hist_SiPM_Y_1.at(j+i*N_Z)->GetBinContent(bin) - h_hit_hor_1->GetBinContent(bin))*log(1 - h_prob_hor_1.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_hor_1->GetBinContent(bin)*log(h_prob_hor_1.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_Y_1<< endl;
                    logLY +=  temp_log_Y_1;
                }
                for(int bin =1; bin<= h_prob_hor_2.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    double temp_log_Y_2 = (hist_SiPM_Y_2.at(j+i*N_Z)->GetBinContent(bin) - h_hit_hor_2->GetBinContent(bin))*log(1 - h_prob_hor_2.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_hor_2->GetBinContent(bin)*log(h_prob_hor_2.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_Y_2<< endl;
                    logLY +=  temp_log_Y_2;
                }
                for(int bin =1; bin<= h_prob_hor_3.at(j+i*N_Z)->GetNbinsX(); bin++){ 
                    double temp_log_Y_3 = (hist_SiPM_Y_3.at(j+i*N_Z)->GetBinContent(bin) - h_hit_hor_3->GetBinContent(bin))*log(1 - h_prob_hor_3.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_hor_3->GetBinContent(bin)*log(h_prob_hor_3.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_Y_3<< endl;
                    logLY +=  temp_log_Y_3;
                }
                for(int bin =1; bin<= h_prob_hor_4.at(j+i*N_Z)->GetNbinsX(); bin++){
                    double temp_log_Y_4 = (hist_SiPM_Y_4.at(j+i*N_Z)->GetBinContent(bin) - h_hit_hor_4->GetBinContent(bin))*log(1 - h_prob_hor_4.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_hor_4->GetBinContent(bin)*log(h_prob_hor_4.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_Y_4<< endl;
                    logLY +=  temp_log_Y_4;
                }
                for(int bin =1; bin<= h_prob_hor_5.at(j+i*N_Z)->GetNbinsX(); bin++){
                    double temp_log_Y_5 = (hist_SiPM_Y_5.at(j+i*N_Z)->GetBinContent(bin) - h_hit_hor_5->GetBinContent(bin))*log(1 - h_prob_hor_5.at(j+i*N_Z)->GetBinContent(bin) + epsilon) + h_hit_hor_5->GetBinContent(bin)*log(h_prob_hor_5.at(j+i*N_Z)->GetBinContent(bin)+epsilon);
                    if (debug) cout<< temp_log_Y_5<< endl;
                    logLY +=  temp_log_Y_5;
                }

                file<< energy_test[i]<< "\t" << Z_test[j]<<"\t" <<logLX << "\t"<< logLY << "\t" <<logLX+logLY << endl;
                energy.push_back(energy_test[i]);
                Z.push_back(Z_test[j]);
                logL.push_back(logLX+logLY);
                if (debug) cout<< logLX << "\t"<< logLY << "\t" <<logLX+logLY << endl;
            }//end of Z value
        }// end of energy
        //find the maximum log value
        double  logL_max = *max_element(logL.begin(), logL.end());

        cout<<"Evaluating delta_nll"<<endl;
        std::vector<double> delta_nll;
        int min_pos = 0;
        double min = 2*(logL_max-logL.at(0));
        for(int i =0; i<logL.size(); i++){
            delta_nll.push_back(2*(logL_max-logL.at(i)));
            if (debug) cout<<energy.at(i)<<"\t"<<Z.at(i)<<"\t"<<delta_nll.at(i)<<endl; 
            if (delta_nll.at(i)<min){
                min = delta_nll.at(i);
                min_pos = i;
            }
        }
        cout<< " Minimum delta_nll " << min << "\t" << min_pos << "\t" << energy.at(min_pos)<< Z.at(min_pos) <<endl;
        // Filling reco histogram with the energy and Z value with the lowest delta_nll value
        h_reco -> Fill (energy.at(min_pos),Z.at(min_pos));
        h_reco_energy -> Fill (energy.at(min_pos));
        h_reco_Z -> Fill (Z.at(min_pos));
	file << endl;
    }//end of events loop 

    TCanvas * c_X_Prob = new TCanvas();
    h_reco->Draw("COLZ");
    c_X_Prob->SaveAs(("/afs/cern.ch/work/r/ridz01/public/Reco_2D_"+std::to_string(energy_true)+"_"+std::to_string(Z_true)+".png").c_str());

    TCanvas * c_1 = new TCanvas();
    h_reco_energy->Draw("COLZ");
    c_1->SaveAs(("/afs/cern.ch/work/r/ridz01/public/Reco_Energy_"+std::to_string(energy_true)+"_"+std::to_string(Z_true)+".png").c_str());

    TCanvas * c_2 = new TCanvas();
    h_reco_Z->Draw("COLZ");
    c_2->SaveAs(("/afs/cern.ch/work/r/ridz01/public/Reco_Z_"+std::to_string(energy_true)+"_"+std::to_string(Z_true)+".png").c_str());

    TFile * out_plots = new TFile(("/afs/cern.ch/work/r/ridz01/public/Reco_histograms_"+std::to_string(energy_true)+"_"+std::to_string(Z_true)+".root").c_str(), "RECREATE");
    h_reco->Write();
    h_reco_energy->Write();
    h_reco_Z->Write();
}



