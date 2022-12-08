/* Generate histograms of energies for each energy bin determined by adaptive binning
(Scifi hits true)*/

int Find_Ev_Stat(std::vector <int> hit_sta, double thres){
    // hits_per_sta[] is an array with the no of hits in each station
    // ex: hits_per_sta[0] -> no of hits in station 0
    int hits_per_sta[5]={0};
    int valid_hits =0;
    // for loop to count the no of hits for each station
    for (int i_hit =0; i_hit<hit_sta.size();i_hit++){
        if (hit_sta[i_hit]<=5&&hit_sta[i_hit]>0){
            hits_per_sta[hit_sta[i_hit]-1]++; 
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

void scifi_hist_energy_ScifiPoint(){
    
    cout<<"entered Script"<<endl;
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
     
    /*// Scifi hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/sndScifiHit.h
    TClonesArray * DigiScifiHits = new TClonesArray("sndScifiHit");
    cbmsim->SetBranchAddress("Digi_ScifiHits", &DigiScifiHits);*/

    // ScifiPoint hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/ScifiPoint.h
    TClonesArray * ScifiHits = new TClonesArray("ScifiPoint");
    cbmsim->SetBranchAddress("ScifiPoint", &ScifiHits);
     
    // Muon filter hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/MuFilterHit.h
    TClonesArray * DigiMuFilterHits = new TClonesArray("MuFilterHit");
    cbmsim->SetBranchAddress("Digi_MuFilterHits", &DigiMuFilterHits);

    //Declare object of TDatabasePDG
    TDatabasePDG *pdg = new TDatabasePDG();

    // STEP 3: DECALRE MUTLIPLE N HISTOGRAMS

    // Reading the binnings from the txt file
    std::vector<int> bin_0;
    std::vector<int> bin_1;
    std::vector<int> bin_2;
    std::vector<int> bin_3;
    std::vector<int> bin_4;

    int ctr_0=0;
    int ctr_1=0;
    int ctr_2=0;
    int ctr_3=0;
    int ctr_4=0;

    ifstream input_file("/eos/user/r/ridz01/EnergyBins_0_200.txt");
    while (!input_file.eof()) {
        int tmp;
        input_file >> tmp;
        bin_0.push_back(tmp);
        ctr_0++;  
    }

    ifstream input_file_1("/eos/user/r/ridz01/EnergyBins_1_200.txt");
    while (!input_file_1.eof()) {
        int tmp;
        input_file_1 >> tmp;
        bin_1.push_back(tmp);
        ctr_1++;  
    }

    ifstream input_file_2("/eos/user/r/ridz01/EnergyBins_2_200.txt");
    while (!input_file_2.eof()) {
        int tmp;
        input_file_2 >> tmp;
        bin_2.push_back(tmp);
        ctr_2++;  
    }

    ifstream input_file_3("/eos/user/r/ridz01/EnergyBins_3_200.txt");
    while (!input_file_3.eof()) {
        int tmp;
        input_file_3 >> tmp;
        bin_3.push_back(tmp);
        ctr_3++;  
    }

    ifstream input_file_4("/eos/user/r/ridz01/EnergyBins_4_200.txt");
    while (!input_file_4.eof()) {
        int tmp;
        input_file_4 >> tmp;
        bin_4.push_back(tmp);
        ctr_4++;  
    }
    // Number of bins would be no of edges -1
    int numbin_0 =ctr_0-1, numbin_1= ctr_1-1, numbin_2= ctr_2-1, numbin_3 = ctr_3-1, numbin_4 = ctr_4-1;

    // Histogram for Scifihits NC events
    std::vector<TH1D*> hist_vec_scifi_NC_en_00;
    std::vector<TH1D*> hist_vec_scifi_NC_en_01;
    std::vector<TH1D*> hist_vec_scifi_NC_en_02;
    std::vector<TH1D*> hist_vec_scifi_NC_en_03;
    std::vector<TH1D*> hist_vec_scifi_NC_en_04;
    // Fixed binning
    for (int i = 0; i<(ctr_0)-1; i++){
        hist_vec_scifi_NC_en_00.push_back(new TH1D(("h_Scifi_NC_en_00_"+std::to_string(i)).c_str(),";Energy[GeV];Number of events",10,0,2500)); 
    }
    for (int i = 0; i<(ctr_1)-1; i++){
        hist_vec_scifi_NC_en_01.push_back(new TH1D(("h_Scifi_NC_en_01_"+std::to_string(i)).c_str(),";Energy[GeV];Number of events",10,0,2500));
    } 
    for (int i = 0; i<(ctr_2)-1; i++){
        hist_vec_scifi_NC_en_02.push_back(new TH1D(("h_Scifi_NC_en_02_"+std::to_string(i)).c_str(),";Energy[GeV];Number of events",10,0,2500));
    } 
    for (int i = 0; i<(ctr_3)-1; i++){
        hist_vec_scifi_NC_en_03.push_back(new TH1D(("h_Scifi_NC_en_03_"+std::to_string(i)).c_str(),";Energy[GeV];Number of events",10,0,2500));
    }
    for (int i = 0; i<(ctr_4)-1; i++){
        hist_vec_scifi_NC_en_04.push_back(new TH1D(("h_Scifi_NC_en_04_"+std::to_string(i)).c_str(),";Energy[GeV];Number of events",10,0,2500));
    }

    int i_hist=0; // histogram no.
    // Main loop for the Scifi energies
    for (int i_event=0; i_event<cbmsim->GetEntries(); i_event++){
        cout<<"Entering For loop for event"<< i_event<< endl;
        cbmsim->GetEntry(i_event);
        double enu = ((ShipMCTrack*)MCTracks->At(0))->GetEnergy();
        int Start_X = ((ShipMCTrack*)MCTracks->At(0))->GetStartX();
        int Start_Y = ((ShipMCTrack*)MCTracks->At(0))->GetStartY();
        double enu_hadron = 0;
        cout<<" Asignment of energies done" << endl;
        //Cutoff to consider the particles in the Scifi Filters
        if (((-8.0>Start_X)&&(Start_X>-47.0) && (15.5<Start_Y)&&(Start_Y<54.5))){
            // NC events selection
            cout<<" Selecting the event if it is NC"<< endl;
            if (TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==12||TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==14){
                std::vector <int> hit_station;
                int n_hits_Sci=0; int charge_hadron=0;
                cout<<"Entering for loop for Scifihits"<<endl;
                cout<<"No of enteries for Scifi for event "<<i_event<<"  "<<ScifiHits->GetEntries()<<endl;
                for (int i=0; i<ScifiHits->GetEntries();i++){
                    int pdg_hadron = ((ScifiPoint*)ScifiHits->At(i))->PdgCode();
                    //cout<<"Pdg Hadron:"<<pdg_hadron<<endl;
                    if (pdg->GetParticle(pdg_hadron)){
                        charge_hadron = (pdg->GetParticle(pdg_hadron))->Charge();
                    }
                    else { charge_hadron = 1; }
                    //cout<<"Charge of the Hadron:"<<charge_hadron<<endl;
                    if (charge_hadron!=0){
                        n_hits_Sci++;
                        hit_station.push_back(((ScifiPoint*)ScifiHits->At(i))->station());
                    }
                }
                cout<<"No of valid hits "<< n_hits_Sci<<endl;
                if(hit_station.size()==0){continue;}
                int ev_sta = Find_Ev_Stat(hit_station, 0.05);
                cout<<"Event Station"<<ev_sta<<endl;               
                enu_hadron = enu - ((ShipMCTrack*)MCTracks->At(1))->GetEnergy();
                //int i_hist = (int)n_hits_Sci/delta_SF[ev_sta]; //Uniform binning

                if(ev_sta==0){ 
                    for( int i=0; i< numbin_0; i++){
                            if (bin_0.at(i)<=n_hits_Sci && bin_0.at(i+1)>n_hits_Sci) { 
                                cout << bin_0.at(i) << " " <<bin_0.at(i+1)<<" "<<i;
                                i_hist=i;}
                    }
                    cout<<i_hist<<endl;
                    cout<<hist_vec_scifi_NC_en_00.size()<<endl;
                    hist_vec_scifi_NC_en_00.at(i_hist)->Fill(enu_hadron);
                    cout<<"Hist station NC 0"<<i_hist<<endl;}
                if(ev_sta==1){                        
                    for( int i=0; i< numbin_1; i++){ 
                            if (bin_1.at(i)<=n_hits_Sci&& bin_1.at(i+1)>n_hits_Sci) { 
                                cout << bin_1.at(i) << " " <<bin_1.at(i+1)<<" "<<i;
                                i_hist=i;} 
                    }
                    cout<<i_hist<<endl;
                    cout<<hist_vec_scifi_NC_en_01.size()<<endl;
                    hist_vec_scifi_NC_en_01.at(i_hist)->Fill(enu_hadron);
                    cout<<"Hist station NC 1"<<i_hist<<endl;}
                if(ev_sta==2){
                    for( int i=0; i< numbin_2; i++){ 
                            if (bin_2.at(i)<=n_hits_Sci && bin_2.at(i+1)>n_hits_Sci) { 
                                cout << bin_2.at(i) << " " <<bin_2.at(i+1)<<" "<<i <<endl;
                                i_hist=i;} 
                    }
                    cout<<i_hist<<endl;
                    cout<<hist_vec_scifi_NC_en_02.size()<<endl;
                    hist_vec_scifi_NC_en_02.at(i_hist)->Fill(enu_hadron);
                    cout<<"Hist station NC 2"<<i_hist<<endl;}
                if(ev_sta==3){
                    for( int i=0; i< numbin_3; i++){ 
                            if (bin_3.at(i)<=n_hits_Sci && bin_3.at(i+1)>n_hits_Sci) { 
                                cout << bin_3.at(i) << " " <<bin_3.at(i+1)<<i;
                                i_hist=i;} 
                    }
                    cout<<i_hist<<endl;
                    cout<<hist_vec_scifi_NC_en_03.size()<<endl;
                    hist_vec_scifi_NC_en_03.at(i_hist)->Fill(enu_hadron);
                    cout<<"Hist station NC 3"<<i_hist<<endl;}
                if(ev_sta==4){
                    for( int i=0; i< numbin_4; i++){ 
                        if (bin_4.at(i)<=n_hits_Sci && bin_4.at(i+1)>n_hits_Sci) { 
                            cout << bin_4.at(i) << " " <<bin_4.at(i+1)<<i;
                            i_hist=i;} 
                    }
                    cout<<i_hist<<endl;
                    cout<<hist_vec_scifi_NC_en_04.size()<<endl;
                    hist_vec_scifi_NC_en_04.at(i_hist)->Fill(enu_hadron);
                    cout<<"Hist station NC 4"<<i_hist<<endl;}
            cout<<"Exiting Histogram assignment"<<endl;
            }
        }
    cout<<"Exiting for event "<<i_event<<endl;
    }

    // Saving Histograms
    TFile* fout3 = new TFile("/eos/user/r/ridz01/Hist_SciFi_NC_energy_ScifiPoint_irregular_101022_200files.root","RECREATE");
    for (int i_hist=0; i_hist<numbin_0; i_hist++){
        hist_vec_scifi_NC_en_00.at(i_hist)->Write();
    }
    for (int i_hist=0; i_hist<numbin_1; i_hist++){
        hist_vec_scifi_NC_en_01.at(i_hist)->Write();
    }
    for (int i_hist=0; i_hist<numbin_2; i_hist++){
        hist_vec_scifi_NC_en_02.at(i_hist)->Write();
    }
    for (int i_hist=0; i_hist<numbin_3; i_hist++){
        hist_vec_scifi_NC_en_03.at(i_hist)->Write();
    }
    for (int i_hist=0; i_hist<numbin_4; i_hist++){
         hist_vec_scifi_NC_en_04.at(i_hist)->Write();
    }

    cout<<"Filled scifi_NC hists"<<endl;
    fout3->Close(); 

}