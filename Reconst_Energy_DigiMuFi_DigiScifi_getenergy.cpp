/* Script to plot energy reconstruction plots based on the TGraphs of the no of hits(/Getenergy())
(digitised) v/s energy for Scifi and Muon Filter */

/* Find_Ev_Stat() - this method generate the event station based on the hits vector 
and the chosen threshold of 5% hits
std::vector <int> hit_sta - vector containing all the stations in which each event occured*/
int Find_Ev_Stat(std::vector <int> hit_sta, double thres){
    // hits_per_sta[] is an array with the no of hits in each station
    // ex: hits_per_sta[0] -> no of hits in station 0
    int hits_per_sta[5]={0};
    int valid_hits =0;
    /*for loop to count the no of hits for each station based on the hit_sta vector passed
    to the method*/  
    for (int i_hit =0; i_hit<hit_sta.size();i_hit++){
        cout<<"Hits station vector"<<hit_sta[i_hit]<<endl;
        if (hit_sta[i_hit]<=5&&hit_sta[i_hit]>0){
            hits_per_sta[hit_sta[i_hit]-1]++; 
            cout <<"Hits per station at station "<<hit_sta[i_hit]-1<<"is "<<hits_per_sta[hit_sta[i_hit]-1]<<endl;
            valid_hits++;
        }
        else {cout<<"WEIRD NUMBER"<<hit_sta[i_hit]<<endl;} 
    }
    cout<<"No of hits in station complete"<<endl;
    // calculate the fraction of the hits and determine the event station based on the 5% threshold value
    double f_hit =0;
    for(int find_sta=0;find_sta<5;find_sta++){
        f_hit = f_hit + hits_per_sta[find_sta]/(double)valid_hits;
        cout<<f_hit<<endl;
        if (f_hit>thres){ 
            cout << find_sta; 
            return find_sta;
        }
    }
    cout<<"Assignment of event station complete"<<endl;
    return -1;
}

void Reconst_Energy_DigiMuFi_DigiScifi_getenergy(){

    // STEP 1: READING MULTIPLE FILES
    // Instead of opening the TFile directly, create a "TChain" of TTrees:
     TChain * cbmsim = new TChain("cbmsim");
    // Now we can add the files one by one to the chain. I generated files in directories labeled 0 to 199 (some might be empty due to jobs which haven't finished yet or have failed).
    for (int i_job = 0; i_job < 200; i_job++){
      cbmsim->Add(("/eos/user/c/cvilela/SND_MC_June21/neutrino/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/"+std::to_string(i_job)+"/sndLHC.Genie-TGeant4_digCPP.root").c_str());
    }

    // STEP 2: DECLARE OBJECTS TO READ THE FILES
    // MC tracks. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipdata/ShipMCTrack.h
    TClonesArray * MCTracks = new TClonesArray("ShipMCTrack");
    cbmsim->SetBranchAddress("MCTrack", &MCTracks);
     
    // Scifi hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/sndScifiHit.h
    TClonesArray * DigiScifiHits = new TClonesArray("sndScifiHit");
    cbmsim->SetBranchAddress("Digi_ScifiHits", &DigiScifiHits);

    /*// ScifiPoint hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/ScifiPoint.h
    TClonesArray * ScifiHits = new TClonesArray("ScifiPoint");
    cbmsim->SetBranchAddress("ScifiPoint", &ScifiHits);*/
     
    // Muon filter hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/MuFilterHit.h
    TClonesArray * DigiMuFilterHits = new TClonesArray("MuFilterHit");
    cbmsim->SetBranchAddress("Digi_MuFilterHits", &DigiMuFilterHits);

    //Declare object of TDatabasePDG
    TDatabasePDG *pdg = new TDatabasePDG();

    // STEP 3: Reading the TGraphs
    // Scifi
    // Original TGraph 
    std::vector<TGraphErrors*> sci_en_NC;
    TFile* fin2_Sc = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Scifi_NC_energy_DigiScifi_irregular_231122.root","READ");// All Bins
    for (int i=0; i<5;i++){
        sci_en_NC.push_back((TGraphErrors*)fin2_Sc->Get(("Station_"+std::to_string(i)+"_sig;1").c_str()));
    }

    // TGraph of Errors
    std::vector<TGraphErrors*> sci_en_NC_err;
    TFile* fin2_err_Sc = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Scifi_NC_energy_errors_DigiScifi_irregular_231122.root","READ"); // All Bins
    for (int i=0; i<5;i++){
        sci_en_NC_err.push_back((TGraphErrors*)fin2_err_Sc->Get(("Station_"+std::to_string(i)+";1").c_str()));
    }

    //MuFi
    // TGraphs
    std::vector<TGraphErrors*> mufi_en_NC;
    //TFile* fin2_Mu = new TFile("/eos/user/r/ridz01/TGraph/TGraph_MuonFilter_NC_energy_DigiMufilter_051122.root","READ"); // All Bins","READ");
    TFile* fin2_Mu = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Mufi_NC_energy_DigiMufi_getenergy_irregular_241122.root","READ"); // All Bins","READ");
    for (int i=0; i<5;i++){
        mufi_en_NC.push_back((TGraphErrors*)fin2_Mu->Get(("Station_"+std::to_string(i)+"_sig;1").c_str()));
    }

    // TGraph of Errors
    std::vector<TGraphErrors*> mufi_en_NC_err;
    //TFile* fin2_err_Mu = new TFile("/eos/user/r/ridz01/TGraph/TGraph_MuonFilter_NC_energy_errors_DigiMufilter_051122.root","READ");
    TFile* fin2_err_Mu = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Mufi_NC_energy_DigiMufi_getenergy_errors_irregular_241122.root","READ");
    for (int i=0; i<5;i++){
        mufi_en_NC_err.push_back((TGraphErrors*)fin2_err_Mu->Get(("Station_"+std::to_string(i)+";1").c_str()));
    }

    // STEP 2: DECLARE HISTOGRAMS FOR THE ENERGY RESOLUTIONS

    // THID (1D plot) of the energy resolution (Avg_NC_1D_rec), true energy(Avg_NC_1D_etrue) and reconstructed energy(Avg_NC_1D_ereco)
    std::vector<TH1D*> Avg_NC_1D_rec;
    std::vector<TH1D*> Avg_NC_1D_etrue;
    std::vector<TH1D*> Avg_NC_1D_ereco;
    for (int i=0;i<5;i++){
        Avg_NC_1D_rec.push_back(new TH1D(("Avg_NC_Station_"+std::to_string(i)).c_str(), ";(E_rec-E_true)/E_true [GeV];Number of events", 100, -3, 3));
        Avg_NC_1D_etrue.push_back(new TH1D(("Avg_NC_Station_"+std::to_string(i)).c_str(), ";E_true [GeV];Number of events", 100, 0, 2500));
        Avg_NC_1D_ereco.push_back(new TH1D(("Avg_NC_Station_"+std::to_string(i)).c_str(), ";E_rec [GeV];Number of events", 100, 0, 2500));
    }
    // TH1D plot of energy resolution with all the Stations combined
    TH1D * Avg_NC_1D_rec_com = new TH1D("Avg_NC", ";(E_rec-E_true)/E_true [GeV];Number of events", 100, -3, 3);

    // TH2D 2-D plot of the energy resolution for each station where E_true is the x-axis and (E_rec-E_true)/E_true [GeV] is the y-axis
    std::vector<TH2D*> Avg_NC_rec;
    for (int i=0;i<5;i++){
        Avg_NC_rec.push_back(new TH2D(("Avg_NC"+std::to_string(i)).c_str(), "; E_true [GeV];(E_rec-E_true)/E_true [GeV]", 100, 0, 2500,100, -3, 3));
    }

    //TH2D 2-D plot for all stations combined
    TH2D * Avg_NC_rec_com = new TH2D("Avg_NC", "; E_true [GeV];(E_rec-E_true)/E_true [GeV]", 100, 0, 2500,100, -3, 3);

    // Main loop over the all the events
    for (int i_event=0; i_event<cbmsim->GetEntries(); i_event++){
        cout<<"Entering For loop for event"<< i_event<< endl;
        cbmsim->GetEntry(i_event);
        double enu = ((ShipMCTrack*)MCTracks->At(0))->GetEnergy();
        int Start_X = ((ShipMCTrack*)MCTracks->At(0))->GetStartX();
        int Start_Y = ((ShipMCTrack*)MCTracks->At(0))->GetStartY();
        int n_hits_Sci = 0;
        double enu_hadron = 0;
        double get_energy = 0;
        cout<<" Asignment of energies done" << endl;
        //Cutoff to consider the particles in the Scifi Filters
        if (((-8.0>Start_X)&&(Start_X>-47.0) && (15.5<Start_Y)&&(Start_Y<54.5))){
            // NC events selection condition
            if (TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==12||TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==14){
                std::vector <int> hit_station;
                cout<<"No of enteries for Scifi for event "<<i_event<<"  "<<DigiScifiHits->GetEntries()<<endl;
                // Count the total no of "valid" hits and generate the hit_station vector containing all the stations in which the event has occured
                for (int i=0; i<DigiScifiHits->GetEntries();i++){
                    if(((sndScifiHit*)DigiScifiHits->At(i))->isValid()){
                        hit_station.push_back(((sndScifiHit*)DigiScifiHits->At(i))->GetStation());
                        n_hits_Sci++;
                    }
                }
                cout<<"No of valid hits in Scifi"<< n_hits_Sci<<endl;
                if(hit_station.size()==0){continue;}
                // Determine the station at which the event is occuring using the Find_Ev_Stat function
                int i_hist = Find_Ev_Stat(hit_station, 0.05);
                cout<<"Event Station"<<i_hist<<endl; 
                // Determine the hits based on the getenergy() in MuFilter
                for (int i_hit=0; i_hit<DigiMuFilterHits->GetEntries(); i_hit++){
                    if (((MuFilterHit*)DigiMuFilterHits->At(i_hit))->GetSystem()==2) { 
                        get_energy = get_energy+((MuFilterHit*)DigiMuFilterHits->At(i_hit))->GetEnergy();
                    }//System 2
                }
                cout<<"Getenergy value in Mufi"<< get_energy<<endl;              
                enu_hadron = enu - ((ShipMCTrack*)MCTracks->At(1))->GetEnergy();// True energy
                cout<<"True Energy "<<enu_hadron<<endl;
                // E_sf is read from the TGraph using the Eval function
                double E_sf = sci_en_NC.at(i_hist)->Eval(n_hits_Sci); // i_hist is the station number 
                double dE_sf = sci_en_NC_err.at(i_hist)->Eval(n_hits_Sci);
                double E_mf = mufi_en_NC.at(i_hist)->Eval(get_energy); // 0 because only one plot(previous analysis). i_hist since we have 5 stations
                double dE_mf = mufi_en_NC_err.at(i_hist)->Eval(get_energy);
                //double Erec = (E_sf/pow(dE_sf,2)+E_us/pow(dE_us,2) + E_ds/pow(dE_ds,2))/(1/pow(dE_sf,2)+1/pow(dE_us,2) + 1/pow(dE_ds,2));
                double Erec = (E_sf/pow(dE_sf,2)+E_mf/pow(dE_mf,2))/(1/pow(dE_sf,2)+1/pow(dE_mf,2));
                cout<<E_sf<<" "<<dE_sf<<" "<<E_mf<<" "<<dE_mf<<" "<<(Erec-enu_hadron)/enu_hadron<<endl;
                Avg_NC_rec.at(i_hist)->Fill(enu_hadron,(Erec-enu_hadron)/enu_hadron);
                Avg_NC_1D_rec.at(i_hist)->Fill((Erec-enu_hadron)/enu_hadron);
                Avg_NC_rec_com->Fill(enu_hadron,(Erec-enu_hadron)/enu_hadron);
                Avg_NC_1D_rec_com->Fill((Erec-enu_hadron)/enu_hadron);
                Avg_NC_1D_ereco.at(i_hist)->Fill(Erec);
                Avg_NC_1D_etrue.at(i_hist)->Fill(enu_hadron);
            }
        }
    }

    TFile* fout2_Sc = new TFile("/eos/user/r/ridz01/Erec_NC_SF_Station_Digi_251122.root","RECREATE");
    for (int i=0; i<5; i++){        
        Avg_NC_1D_rec.at(i) ->Write(("Station_"+std::to_string(i)+"_1D").c_str());
        Avg_NC_1D_ereco.at(i) -> Write(("Station_"+std::to_string(i)+"_e_reco").c_str());
        Avg_NC_1D_etrue.at(i) -> Write(("Station_"+std::to_string(i)+"_e_true").c_str());
        Avg_NC_rec.at(i)->Write(("Station_"+std::to_string(i)+"_2D").c_str());  
    }
    fout2_Sc-> Close();

    TFile* fout2 = new TFile("/eos/user/r/ridz01/Erec_NC_SF_Digi_251122.root","RECREATE");
    Avg_NC_rec_com -> Write();
    Avg_NC_1D_rec_com -> Write();
    fout2-> Close();

    TCanvas *c_Avg_NC_rec_0 = new TCanvas();
    Avg_NC_rec.at(0)->Draw("COLZ");
    c_Avg_NC_rec_0->SaveAs("Rec_NC_Digi_00_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_0_1D = new TCanvas();
    Avg_NC_1D_rec.at(0)->Draw("COLZ");
    c_Avg_NC_rec_0_1D->SaveAs("Rec_NC_Digi_00_1D_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_1 = new TCanvas();
    Avg_NC_rec.at(1)->Draw("COLZ");
    c_Avg_NC_rec_1->SaveAs("Rec_NC_Digi_01_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_1_1D = new TCanvas();
    Avg_NC_1D_rec.at(1)->Draw("COLZ");
    c_Avg_NC_rec_1_1D->SaveAs("Rec_NC_Digi_01_1D_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_2 = new TCanvas();
    Avg_NC_rec.at(2)->Draw("COLZ");
    c_Avg_NC_rec_2->SaveAs("Rec_NC_Digi_02_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_2_1D = new TCanvas();
    Avg_NC_1D_rec.at(2)->Draw("COLZ");
    c_Avg_NC_rec_2_1D->SaveAs("Rec_NC_Digi_02_1D_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_3 = new TCanvas();
    Avg_NC_rec.at(3)->Draw("COLZ");
    c_Avg_NC_rec_3->SaveAs("Rec_NC_Digi_03_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_3_1D = new TCanvas();
    Avg_NC_1D_rec.at(3)->Draw("COLZ");
    c_Avg_NC_rec_3_1D->SaveAs("Rec_NC_Digi_03_1D_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_4 = new TCanvas();
    Avg_NC_rec.at(4)->Draw("COLZ");
    c_Avg_NC_rec_4->SaveAs("Rec_NC_Digi_04_251122(getenergy).png");

    TCanvas *c_Avg_NC_rec_4_1D = new TCanvas();
    Avg_NC_1D_rec.at(4)->Draw("COLZ");
    c_Avg_NC_rec_4_1D->SaveAs("Rec_NC_Digi_04_1D_251122(getenergy).png");

    TCanvas *c_Avg_NC = new TCanvas();
    Avg_NC_rec_com->Draw("COLZ");
    c_Avg_NC->SaveAs("Rec_NC_Digi_251122(getenergy).png");

    TCanvas *c_Avg_NC_1D_com = new TCanvas();
    Avg_NC_1D_rec_com->Draw();
    c_Avg_NC_1D_com->SaveAs("Rec_NC_Digi_1D_251122(getenergy).png");
}