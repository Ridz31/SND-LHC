/* Script to determine the binning edges to perform adaptive binning for digitised Scifi */
#include <algorithm>
// Determining the event station
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

void scifi_hist_energy_DigiScifi_binning(){
    
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

    //Declare object of TDatabasePDG
    TDatabasePDG *pdg = new TDatabasePDG();

    // No of counts per bin
    int countsPerBin = 100;
    
    std::vector<double> no_of_hits_0; 
    std::vector<double> no_of_hits_1;
    std::vector<double> no_of_hits_2;
    std::vector<double> no_of_hits_3;
    std::vector<double> no_of_hits_4;
    int ctr0=0, ctr1=0,ctr2=0,ctr3=0,ctr4=0;
    // Main loop for no of hits vectors
    for (int i_event=0; i_event<cbmsim->GetEntries(); i_event++){
        cbmsim->GetEntry(i_event);
        int Start_X = ((ShipMCTrack*)MCTracks->At(0))->GetStartX();
        int Start_Y = ((ShipMCTrack*)MCTracks->At(0))->GetStartY();
        int n_hits_Sci = 0;
        //Cutoff to consider the particles in the Scifi Filters
        if (((-8.0>Start_X)&&(Start_X>-47.0) && (15.5<Start_Y)&&(Start_Y<54.5))){
            // NC events slection
            if (TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==12||TMath::Abs(((ShipMCTrack*)MCTracks->At(1))->GetPdgCode())==14){
                std::vector <int> hit_station;
                // Determining the valid hits and counting the total number of hits for each event
                for (int i=0; i<DigiScifiHits->GetEntries();i++){
                    if (((sndScifiHit*)DigiScifiHits->At(i))->isValid()){
                        hit_station.push_back(((sndScifiHit*)DigiScifiHits->At(i))->GetStation());
                        n_hits_Sci++;
                    } 
                }
                if(hit_station.size()==0){continue;} 
                int ev_sta = Find_Ev_Stat(hit_station, 0.05); // The event station is returned and stored in ev_sta 
                // checking the station number
                if(ev_sta==0){
                    no_of_hits_0.push_back(n_hits_Sci);
                    ctr0++;
                }
                if(ev_sta==1){
                    no_of_hits_1.push_back(n_hits_Sci);
                    ctr1++;
                }
                if(ev_sta==2){
                    no_of_hits_2.push_back(n_hits_Sci);
                    ctr2++;               
                }
                if(ev_sta==3){
                    no_of_hits_3.push_back(n_hits_Sci);
                    ctr3++;
                }
                if(ev_sta==4){
                   no_of_hits_4.push_back(n_hits_Sci);
                   ctr4++;
                }
            }
        }
    }
    
    // No of bins per station
    int numbin_0 = ctr0/countsPerBin; cout<<numbin_0<<endl;
    int numbin_1 = ctr1/countsPerBin; cout<<numbin_1<<endl;
    int numbin_2 = ctr2/countsPerBin; cout<<numbin_2<<endl;
    int numbin_3 = ctr3/countsPerBin; cout<<numbin_3<<endl;
    int numbin_4 = ctr4/countsPerBin; cout<<numbin_4<<endl;
    
    //Sort the no of hits values
    std::sort(no_of_hits_0.begin(), no_of_hits_0.end());
    std::sort(no_of_hits_1.begin(), no_of_hits_1.end());
    std::sort(no_of_hits_2.begin(), no_of_hits_2.end());
    std::sort(no_of_hits_3.begin(), no_of_hits_3.end());
    std::sort(no_of_hits_4.begin(), no_of_hits_4.end());

    std::vector<double> binEdge_0;
    std::vector<double> binEdge_1;
    std::vector<double> binEdge_2;
    std::vector<double> binEdge_3;
    std::vector<double> binEdge_4;

    ofstream file0;
    file0.open ("/eos/user/r/ridz01/EnergyBins_0_200_Digi.txt");
    binEdge_0.push_back(no_of_hits_0.front());
    file0<<no_of_hits_0.front()<<endl;
    for (int edge=1; edge<numbin_0; edge++) {
        binEdge_0.push_back(no_of_hits_0.at(countsPerBin * edge));
        file0<<binEdge_0.at(edge)<<endl;
    }
    binEdge_0.push_back(no_of_hits_0.back());
    file0<<no_of_hits_0.back();
    file0.close();

    ofstream file1;
    file1.open ("/eos/user/r/ridz01/EnergyBins_1_200_Digi.txt");
    binEdge_1.push_back(no_of_hits_1.front());
    file1<<no_of_hits_1.front()<<endl;
    for (int edge=1; edge<numbin_1; edge++) {
        binEdge_1.push_back(no_of_hits_1.at(countsPerBin * edge));
        file1<<binEdge_1.at(edge)<<endl;
    }
    binEdge_1.push_back(no_of_hits_1.back());
    file1<<no_of_hits_1.back();
    file1.close();

    ofstream file2;
    file2.open ("/eos/user/r/ridz01/EnergyBins_2_200_Digi.txt");
    binEdge_2.push_back(no_of_hits_2.front());
    file2<<no_of_hits_2.front()<<endl;
    for (int edge=1; edge<numbin_2; edge++) {
        binEdge_2.push_back(no_of_hits_2.at(countsPerBin * edge));
        file2<<binEdge_2.at(edge)<<endl;   
    }
    binEdge_2.push_back(no_of_hits_2.back());
    file2<<no_of_hits_2.back();
    file2.close();

    ofstream file3;
    file3.open ("/eos/user/r/ridz01/EnergyBins_3_200_Digi.txt");
    binEdge_3.push_back(no_of_hits_3.front());
    file3<<no_of_hits_3.front()<<endl;
    for (int edge=1; edge<numbin_3; edge++) {
        binEdge_3.push_back(no_of_hits_3.at(countsPerBin * edge));
        file3<<binEdge_3.at(edge)<<endl;   
    }
    binEdge_3.push_back(no_of_hits_3.back());
    file3<<no_of_hits_3.back();
    file3.close();

    ofstream file4;
    file4.open ("/eos/user/r/ridz01/EnergyBins_4_200_Digi.txt");
    binEdge_4.push_back(no_of_hits_4.front());
    file4<<no_of_hits_4.front()<<endl;
      for (int edge=1; edge<numbin_4; edge++) {
        binEdge_4.push_back(no_of_hits_4.at(countsPerBin * edge));
        file4<<binEdge_4.at(edge)<<endl;   
    }
    binEdge_4.push_back(no_of_hits_4.back());
    file4<<no_of_hits_4.back();
    file4.close();

}