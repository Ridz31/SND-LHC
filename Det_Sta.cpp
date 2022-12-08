/* Script to determine the station where the event occured*/

/* This function determines the station at which each event occured*/
int Find_Ev_Stat(std::vector <int> hit_sta, double thres){
    // hits_per_sta[] is an array with the no of hits in each station
    // ex: hits_per_sta[0] -> no of hits in station 0
    int hits_per_sta[5]={0};
    int valid_hits =0;
    // for loop to count the no of hits for each station
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
    // calculate the fraction of the hits and determine the event station based
    // on the threshold value
    double f_hit =0;
    for(int find_sta=0;find_sta<5;find_sta++){
        f_hit = f_hit + hits_per_sta[find_sta]/(double)valid_hits;
        cout<<f_hit<<endl;
        if (f_hit>thres){ 
            cout << find_sta; 
            return find_sta;}
    }
    cout<<"Assignment of event station complete"<<endl;
    return -1;
}

void Det_Sta(){
    
    // STEP 1: READING MULTIPLE FILES
    // Instead of opening the TFile directly, create a "TChain" of TTrees:
     TChain * cbmsim = new TChain("cbmsim");
    // Now we can add the files one by one to the chain. I generated files in directories labeled 0 to 199 (some might be empty due to jobs which haven't finished yet or have failed).
    for (int i_job = 0; i_job < 20; i_job++){
      cbmsim->Add(("/eos/user/c/cvilela/SND_MC_June21/neutrino/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/"+std::to_string(i_job)+"/sndLHC.Genie-TGeant4_digCPP.root").c_str());
    }

    // STEP 2: DECLARE OBJECTS TO READ THE FILES
    // MC tracks. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipdata/ShipMCTrack.h
    TClonesArray * MCTracks = new TClonesArray("ShipMCTrack");
    cbmsim->SetBranchAddress("MCTrack", &MCTracks);
     
    // Scifi hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/sndScifiHit.h
    TClonesArray * DigiScifiHits = new TClonesArray("sndScifiHit");
    cbmsim->SetBranchAddress("Digi_ScifiHits", &DigiScifiHits);

    //STEP 3: DECLARE HISTOGRAM
    int n_hist = 10;
    // 2 D histograms of the hit station
    std::vector<TH2D*> hist_scifi_st_CC;
    for (int i=0; i<n_hist; i++){
        hist_scifi_st_CC.push_back(new TH2D(("Scifi_station_CC"+std::to_string(i)).c_str(), ";Z; Station No", 100,250,380,100,-2,6));
    }
    std::vector<TH2D*> hist_scifi_st_NC;
    for (int i=0; i<n_hist;i++){
        hist_scifi_st_NC.push_back(new TH2D(("Scifi_station_NC"+std::to_string(i)).c_str(), ";Z; Station No", 100,250,380,100,-2,6));
    }
    // 1D histograms of the no of events
    std::vector<TH1D*> hist_scifi_st_CC_1D;
    for (int i=0; i<n_hist;i++){
        hist_scifi_st_CC_1D.push_back(new TH1D(("Scifi_station_CC"+std::to_string(i)).c_str(),";Z-St.No*dz; No of events", 100,250,380));
    }
    std::vector<TH1D*> hist_scifi_st_NC_1D;
    for (int i=0; i<n_hist;i++){
        hist_scifi_st_NC_1D.push_back(new TH1D(("Scifi_station_NC"+std::to_string(i)).c_str(),";Z-St.No*dz; No of events", 100,250,380));
    }

    cout<<"Declaration of objects and histograms complete"<<endl;

    //STEP 4: FILLING THE HISTOGRAM BASED ON THE EVENT STATION
    for (double thres=0.01; thres<=0.10; thres=thres+0.01){
        for (int i_event=0; i_event<cbmsim->GetEntries();i_event++){
            cbmsim->GetEntry(i_event);
            int Start_X = ((ShipMCTrack*)MCTracks->At(0))->GetStartX();
            int Start_Y = ((ShipMCTrack*)MCTracks->At(0))->GetStartY();
            if(((-8.0>Start_X)&&(Start_X>-47.0))&&((15.5<Start_Y)&&(Start_Y<54.5))){
                std::vector <int> hit_station;
                for (int i_hit_st=0; i_hit_st<DigiScifiHits->GetEntries();i_hit_st++){
                    hit_station.push_back(((sndScifiHit*)DigiScifiHits->At(i_hit_st))->GetStation());
                }
                cout<<"Station vector for event "<<i_event<<" initialised"<<endl;
                if(hit_station.size()==0){continue;}
                int ev_sta = Find_Ev_Stat(hit_station, thres); //the station is determined based on the no of hits vector and the chosen threshold
                cout << "For threshold"<<thres<<endl;
                cout<<"The starting station for the event is "<<ev_sta<<endl;
                if (TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==11){
                    hist_scifi_st_CC.at((thres*100-1)) -> Fill(((ShipMCTrack*)MCTracks->At(0))->GetStartZ(),ev_sta);
                    hist_scifi_st_CC_1D.at((thres*100-1)) -> Fill((((ShipMCTrack*)MCTracks->At(0))->GetStartZ())-ev_sta*13);
                }
                cout<<"Filled CC events hist"<<endl;
                if (TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==12||TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==14){
                    hist_scifi_st_NC.at((thres*100-1)) -> Fill(((ShipMCTrack*)MCTracks->At(0))->GetStartZ(),ev_sta);
                    hist_scifi_st_NC_1D.at((thres*100-1)) -> Fill((((ShipMCTrack*)MCTracks->At(0))->GetStartZ())-ev_sta*13);
                }
                cout<<"Filled NC events hits"<<endl;
            }
        }
    }

    //Storing the histograms in root files    
    TFile* fout1 = new TFile("/eos/user/r/ridz01/Hist_scifi_station_CC_test_04.root","RECREATE");
    for (int i_hist=0; i_hist<n_hist; i_hist++){
        hist_scifi_st_CC.at(i_hist)->Write();
        hist_scifi_st_CC_1D.at(i_hist) -> Write();
    }
    fout1->Close();
    cout<<"Histograms for CC events complete"<<endl;

    TFile* fout2 = new TFile("/eos/user/r/ridz01/Hist_scifi_station_NC_test_04.root","RECREATE");
    for (int i_hist=0; i_hist<n_hist; i_hist++){
        hist_scifi_st_NC.at(i_hist)->Write();
        hist_scifi_st_NC_1D.at(i_hist)->Write();
    }
    fout2->Close();
    cout<<"Histograms for NC events complete"<<endl;


    /*//STEP 5: PLOTTING THE HISTOGRAMS
    TCanvas *c_Scifi_st_CC = new TCanvas();
    Scifi_st_CC->Draw("COLZ");
    c_Scifi_st_CC->SaveAs("Scifi_station_CC.png");

    TCanvas *c_Scifi_st_NC = new TCanvas();
    Scifi_st_NC->Draw("COLZ");
    c_Scifi_st_NC->SaveAs("Scifi_station_NC.png");*/
}