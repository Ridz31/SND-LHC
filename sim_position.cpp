/* Generate the hit pattern of the electron simulation event*/
void sim_position(){

    //int energy[10] = {50,100,200,300,400,500,750,1000,1250,1500};
    int energy[1]={1000};
    //int energy[19] = {50,75,100,150,200,250,300,350,400,450,500,625,750,875,1000,1125,1250,1375,1500};
    //int energy[9]={625,750,875,1000,1125,1250,1375,1500};
    //int Z[7] = {285,287,289,291,293,295,297};
    int Z[8] = {285,287,289,291,293,295,297,299};
    //int Z[1]={285};
    for(int j=0;j< 1;j++){
        for ( int i=0; i< 8; i++){
            
            //Reading Digitised simulation files
            TFile * f_in = new TFile(("/eos/user/r/ridz01/PG_electron_sim/"+std::to_string(energy[j])+"GeV/"+std::to_string(Z[i])+"/sndLHC.PG_11-TGeant4_digCPP.root").c_str());
            //TFile * f_in = new TFile(("/afs/cern.ch/work/r/ridz01/public/PG_electron_sim/"+std::to_string(energy[j])+"GeV/"+std::to_string(Z[i])+"/sndLHC.PG_11-TGeant4_digCPP.root").c_str());
            
            // Now get the cbmsim TTree to access each event
            TTree * cbmsim = (TTree*)f_in->Get("cbmsim");

            // Scifi hits. Documentation here: https://github.com/SND-LHC/sndsw/blob/master/shipLHC/sndScifiHit.h
            TClonesArray * DigiScifiHits = new TClonesArray("sndScifiHit");
            cbmsim->SetBranchAddress("Digi_ScifiHits", &DigiScifiHits);

            // Set up geometry of the detector to read the hits

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

            TH1D * h_X_1 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "X_Hits_Station_1;X(mm)", 100, -50, 0);
            TH1D * h_Y_1 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "Y_Hits_Station_1;Y(mm)", 120, 0, 60);

            TH1D * h_X_2 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "X_Hits_Station_2;X(mm)", 100, -50, 0);
            TH1D * h_Y_2 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "Y_Hits_Station_1;Y(mm)", 120, 0, 60);

            TH1D * h_X_3 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "X_Hits_Station_3;X(mm)", 100, -50, 0);
            TH1D * h_Y_3 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "Y_Hits_Station_3;Y(mm)", 120, 0, 60);

            TH1D * h_X_4 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "X_Hits_Station_4;X(mm)", 100, -50, 0);
            TH1D * h_Y_4 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "Y_Hits_Station_4;Y(mm)", 120, 0, 60);

            TH1D * h_X_5 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "X_Hits_Station_5;X(mm)", 100, -50, 0);
            TH1D * h_Y_5 = new TH1D((std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])).c_str(), "Y_Hits_Station_4;Y(mm)", 120, 0, 60);

            // Loop through each event
            for (int i_event = 0; i_event < cbmsim->GetEntries(); i_event++){
                cbmsim->GetEntry(i_event);
                //Initialize two TVector3 objects
                TVector3 a;
                TVector3 b;
                TVector3 mean_position;
                //Loop through Scifi hits:
                for (int scifihit =0; scifihit<DigiScifiHits->GetEntries(); scifihit++){
                    // Loop through all the hits for each event and fill the histograms based on the Station 
    
                    // Condition to check if it is a valid hit
                    if (((sndScifiHit*)DigiScifiHits->At(scifihit))->isValid()){
                        // Identify the station of the hit
                        int hit_station = ((sndScifiHit*)DigiScifiHits->At(scifihit))->GetStation();
                        
                        // Find the SiPM id of the hit event 
                        int id_hit = ((sndScifiHit*)DigiScifiHits->At(scifihit))->GetDetectorID();

                        // To get the X and Y Coordinate
                        ScifiDet->GetSiPMPosition(id_hit, a, b);
                        // Fill the histograms based on the Station

                        // If the hit is on the ZX plane, vertical fibres read the X coordinates
                        if (((sndScifiHit*)DigiScifiHits->At(scifihit))->isVertical()){
                            if (hit_station==1){
                                h_X_1->Fill((a.X() + b.X())/2.);
                            }
                            if (hit_station==2){
                                h_X_2->Fill((a.X() + b.X())/2.);
                            }
                            if (hit_station==3){
                                h_X_3->Fill((a.X() + b.X())/2.);
                            }
                            if (hit_station==4){
                                h_X_4->Fill((a.X() + b.X())/2.);
                            }
                            if (hit_station==5){
                                h_X_5->Fill((a.X() + b.X())/2.);
                            }
                        }
                        //If the hit is on the ZY plane, horizontal fibres read the Y coordinates
                        else{
                            if (hit_station==1){
                                h_Y_1->Fill((a.Y() + b.Y())/2.);
                            }
                            if (hit_station==2){
                                h_Y_2->Fill((a.Y() + b.Y())/2.);
                            }
                            if (hit_station==3){
                                h_Y_3->Fill((a.Y() + b.Y())/2.);
                            }
                            if (hit_station==4){
                                h_Y_4->Fill((a.Y() + b.Y())/2.);
                            }
                            if (hit_station==5){
                                h_Y_5->Fill((a.Y() + b.Y())/2.);
                            }
                        }
                    }// end of valid hit check
                }//end of Scifi hits loop
            }// end of events loop


            TFile* foutX = new TFile(("/eos/user/r/ridz01/Electron_Showers_Simulation_Plots/Hist_Scifi_ElectronShowers_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_1D_XPlane.root").c_str(),"RECREATE");
            h_X_1->Write("Station 1");
            h_X_2->Write("Station 2");
            h_X_3->Write("Station 3");
            h_X_4->Write("Station 4");
            h_X_5->Write("Station 5");
            foutX->Close();

            TFile* foutY = new TFile(("/eos/user/r/ridz01/Electron_Showers_Simulation_Plots/Hist_Scifi_ElectronShowers_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_1D_YPlane.root").c_str(),"RECREATE");
            h_Y_1->Write("Station 1");
            h_Y_2->Write("Station 2");
            h_Y_3->Write("Station 3");
            h_Y_4->Write("Station 4");
            h_Y_5->Write("Station 5");
            foutY->Close();

            ofstream file;
            file.open("/eos/user/r/ridz01/Electron_Showers_Simulation_Plots/Number_of_Entries.txt",ios::app);
            file<<energy[j]<<"\t"<<Z[i]<<"\t"<<cbmsim->GetEntries()<<endl;
        }
    }

}