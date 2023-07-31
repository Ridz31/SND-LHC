/* Generate the probability distribution function (pdf) using the previously generated hit pattern via sim_position.cpp*/
void hist_probability(){

    //int energy[10] = {50,10,200,300,400,500,750,1000,1250,1500};
    //int energy[4]={75,250,350,450};
    //int energy[11] = {50,75,100,150,200,250,300,350,400,450,500};
    //int energy[10] = {625,750,875,1000,1125,1250,1375,1500,1750,2000};
    int energy[1] = {2000};
    int Z[8] = {285,287,289,291,293,295,297,299};
    //int Z[2]={285,297};
    for(int j=0;j< 1;j++){
        for ( int i=0; i< 8; i++){

            // Open simulation files to calculate no of entries
            TFile * f_in = new TFile(("/eos/user/r/ridz01/PG_electron_sim/"+std::to_string(energy[j])+"GeV/"+std::to_string(Z[i])+"/sndLHC.PG_11-TGeant4_digCPP.root").c_str());
            TTree * cbmsim = (TTree*)f_in->Get("cbmsim");
            double nevents = cbmsim->GetEntries(); // Find the number of events simulated and eventually plotted in the hit histogram
            cout<<nevents<<endl;

            TFile* f_X = new TFile(("/eos/user/r/ridz01/Electron_Showers_Simulation_Plots/Hist_Scifi_ElectronShowers_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_1D_XPlane.root").c_str(),"READ");
            TH1D * hist_X_1 = (TH1D*)f_X->Get("Station 1;1");
            TH1D * hist_X_2 = (TH1D*)f_X->Get("Station 2;1");
            TH1D * hist_X_3 = (TH1D*)f_X->Get("Station 3;1");
            TH1D * hist_X_4 = (TH1D*)f_X->Get("Station 4;1");
            TH1D * hist_X_5 = (TH1D*)f_X->Get("Station 5;1");

            TFile* f_Y = new TFile(("/eos/user/r/ridz01/Electron_Showers_Simulation_Plots/Hist_Scifi_ElectronShowers_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_1D_YPlane.root").c_str(),"READ");
            TH1D * hist_Y_1 = (TH1D*)f_Y->Get("Station 1;1");
            TH1D * hist_Y_2 = (TH1D*)f_Y->Get("Station 2;1");
            TH1D * hist_Y_3 = (TH1D*)f_Y->Get("Station 3;1");
            TH1D * hist_Y_4 = (TH1D*)f_Y->Get("Station 4;1");
            TH1D * hist_Y_5 = (TH1D*)f_Y->Get("Station 5;1");

            cout<<"Read files"<<endl;

            // Scaling the histograms based on the number of entires in the histogram. hist_X -> X hits, hist_Y -> Y hits
            hist_X_1->Scale(1.0/nevents);
            hist_X_2->Scale(1.0/nevents);
            hist_X_3->Scale(1.0/nevents);
            hist_X_4->Scale(1.0/nevents);
            hist_X_5->Scale(1.0/nevents);

            hist_Y_1->Scale(1.0/nevents);
            hist_Y_2->Scale(1.0/nevents);
            hist_Y_3->Scale(1.0/nevents);
            hist_Y_4->Scale(1.0/nevents);
            hist_Y_5->Scale(1.0/nevents);

            /*TCanvas * c_Y_Scale = new TCanvas();
            hist_Y_1->Draw();
            c_Y_Scale->SaveAs("Test_Y_Scale.png");

            TCanvas * c_X_Scale = new TCanvas();
            hist_X_1->Draw();
            c_X_Scale->SaveAs("Test_X_Scale.png");*/

            //Saving the scaled histograms in a root file. hist_X -> X hits, hist_Y -> Y hits
            TFile* fout1X = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/ScalePlots_X.root").c_str(),"RECREATE");
            hist_X_1->Write("Station 1");
            hist_X_2->Write("Station 2");
            hist_X_3->Write("Station 3");
            hist_X_4->Write("Station 4");
            hist_X_5->Write("Station 5");
            fout1X->Close();

            TFile* fout1Y = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/ScalePlots_Y.root").c_str(),"RECREATE");
            hist_Y_1->Write("Station 1");
            hist_Y_2->Write("Station 2");
            hist_Y_3->Write("Station 3");
            hist_Y_4->Write("Station 4");
            hist_Y_5->Write("Station 5");
            fout1Y->Close();

        
            // Produce histograms of the number of SiPM chanels
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

            TVector3 a;
            TVector3 b;

            //Horizontal SiPM fibres read the Y-coordinate and Vertical SiPM fibres read the X-coodinate of the hit
            TH1D * h_SiPM_hor_1 = new TH1D("SiPM_horizontal", "Station_1;Y(mm)", 120, 0, 60);
            TH1D * h_SiPM_ver_1 = new TH1D("SiPM_vertical", "Station_1;X(mm)", 100, -50, 0);

            TH1D * h_SiPM_hor_2 = new TH1D("SiPM_horizontal", "Station_2;Y(mm)", 120, 0, 60);
            TH1D * h_SiPM_ver_2 = new TH1D("SiPM_vertical", "Station_2;X(mm)", 100, -50, 0);

            TH1D * h_SiPM_hor_3 = new TH1D("SiPM_horizontal", "Station_3;Y(mm)", 120, 0, 60);
            TH1D * h_SiPM_ver_3 = new TH1D("SiPM_vertical", "Station_3;X(mm)", 100, -50, 0);

            TH1D * h_SiPM_hor_4 = new TH1D("SiPM_horizontal", "Station_4;Y(mm)", 120, 0, 60);
            TH1D * h_SiPM_ver_4 = new TH1D("SiPM_vertical", "Station_4;X(mm)", 100, -50, 0);

            TH1D * h_SiPM_hor_5 = new TH1D("SiPM_horizontal", "Station_5;Y(mm)", 120, 0, 60);
            TH1D * h_SiPM_ver_5 = new TH1D("SiPM_vertical", "Station_5;X(mm)", 100, -50, 0);

            for (int i_sta = 1; i_sta <= 5; i_sta++){
                for (int i_mat = 0; i_mat < 3; i_mat++){
                    for (int i_sipm = 0; i_sipm < 4; i_sipm++){
                        for (int i_chan = 0; i_chan < 128; i_chan++){

                        int id = i_sta*1000000 + i_mat*10000 + i_sipm*1000 + i_chan;

                        ScifiDet->GetSiPMPosition(id, a, b);
                        if (i_sta ==1){
                            h_SiPM_hor_1->Fill((a.Y() + b.Y())/2.);
                            id += 100000; // Vertical
                            ScifiDet->GetSiPMPosition(id, a, b);
                            h_SiPM_ver_1->Fill((a.X() + b.X())/2.);
                        }

                        if (i_sta ==2){
                            h_SiPM_hor_2->Fill((a.Y() + b.Y())/2.);
                            id += 100000; // Vertical
                            ScifiDet->GetSiPMPosition(id, a, b);
                            h_SiPM_ver_2->Fill((a.X() + b.X())/2.);
                        }

                        if (i_sta ==3){
                            h_SiPM_hor_3->Fill((a.Y() + b.Y())/2.);
                            id += 100000; // Vertical
                            ScifiDet->GetSiPMPosition(id, a, b);
                            h_SiPM_ver_3->Fill((a.X() + b.X())/2.);
                        }

                        if (i_sta ==4){
                            h_SiPM_hor_4->Fill((a.Y() + b.Y())/2.);
                            id += 100000; // Vertical
                            ScifiDet->GetSiPMPosition(id, a, b);
                            h_SiPM_ver_4->Fill((a.X() + b.X())/2.);
                        }

                        if (i_sta ==5){
                            h_SiPM_hor_5->Fill((a.Y() + b.Y())/2.);
                            id += 100000; // Vertical
                            ScifiDet->GetSiPMPosition(id, a, b);
                            h_SiPM_ver_5->Fill((a.X() + b.X())/2.);
                        }	  
                        }
                    }
                }      
            }

            /*TCanvas * c_hor_SiPM = new TCanvas();
            h_hor_1->Draw();
            c_hor_SiPM->SaveAs("Test_X_SiPM_1.png");

            TCanvas * c_ver_SiPM = new TCanvas();
            h_ver_1->Draw();
            c_ver_SiPM->SaveAs("Test_Y_SiPM_1.png");*/

            TFile* fout2Y = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/SiPMs_Y.root").c_str(),"RECREATE");
            h_SiPM_hor_1->Write("Station 1 horizontal");
            h_SiPM_hor_2->Write("Station 2 horizontal");
            h_SiPM_hor_3->Write("Station 3 horizontal");
            h_SiPM_hor_4->Write("Station 4 horizontal");
            h_SiPM_hor_5->Write("Station 5 horizontal");
            fout2Y->Close();

            TFile* fout2X = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/SiPMs_X.root").c_str(),"RECREATE");
            h_SiPM_ver_1->Write("Station 1 vertical");
            h_SiPM_ver_2->Write("Station 2 vertical");
            h_SiPM_ver_3->Write("Station 3 vertical");
            h_SiPM_ver_4->Write("Station 4 vertical");
            h_SiPM_ver_5->Write("Station 5 vertical");
            fout2X->Close();

            // Divide the two histograms, the X hits are divided by the vertical SiPM histograms and vice versa
            hist_X_1->Divide(h_SiPM_ver_1);
            hist_X_2->Divide(h_SiPM_ver_2);
            hist_X_3->Divide(h_SiPM_ver_3);
            hist_X_4->Divide(h_SiPM_ver_4);
            hist_X_5->Divide(h_SiPM_ver_5);

            hist_Y_1->Divide(h_SiPM_hor_1);
            hist_Y_2->Divide(h_SiPM_hor_2);
            hist_Y_3->Divide(h_SiPM_hor_3);
            hist_Y_4->Divide(h_SiPM_hor_4);
            hist_Y_5->Divide(h_SiPM_hor_5);

            TFile* fout3X = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/Pdfs_X_hist.root").c_str(),"RECREATE");
            hist_X_1->Draw("HIST");
            hist_X_1->Write("Station 1");
            hist_X_2->Draw("HIST");
            hist_X_2->Write("Station 2");
            hist_X_3->Draw("HIST");
            hist_X_3->Write("Station 3");
            hist_X_4->Draw("HIST");
            hist_X_4->Write("Station 4");
            hist_X_5->Draw("HIST");
            hist_X_5->Write("Station 5");
            fout3X->Close();

            TFile* fout3Y = new TFile(("/eos/user/r/ridz01/Electron_Showers_pdfs/"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z/Pdfs_Y_hist.root").c_str(),"RECREATE");
            hist_Y_1->Draw("HIST");
            hist_Y_1->Write("Station 1");
            hist_Y_2->Draw("HIST");
            hist_Y_2->Write("Station 2");
            hist_Y_3->Draw("HIST");
            hist_Y_3->Write("Station 3");
            hist_Y_4->Draw("HIST");
            hist_Y_4->Write("Station 4");
            hist_Y_5->Draw("HIST");
            hist_Y_5->Write("Station 5");
            fout3Y->Close();
            
            /*TCanvas * c_X_Prob_1 = new TCanvas();
            hist_X_1->Draw("HIST");            
            c_X_Prob_1->SaveAs(("X_pdf_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_Station1.png").c_str());

            TCanvas * c_X_Prob_2 = new TCanvas();
            hist_X_2->Draw("HIST");            
            c_X_Prob_2->SaveAs(("X_pdf_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_Station2.png").c_str());

            TCanvas * c_X_Prob_3 = new TCanvas();
            hist_X_3->Draw("HIST");            
            c_X_Prob_3->SaveAs(("X_pdf_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_Station3.png").c_str());

            TCanvas * c_X_Prob_4 = new TCanvas();
            hist_X_4->Draw("HIST");            
            c_X_Prob_4->SaveAs(("X_pdf_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_Station4.png").c_str());

            TCanvas * c_X_Prob_5 = new TCanvas();
            hist_X_5->Draw("HIST");            
            c_X_Prob_5->SaveAs(("X_pdf_"+std::to_string(energy[j])+"GeV_"+std::to_string(Z[i])+"Z_Station5.png").c_str());*/
        }
    }
}