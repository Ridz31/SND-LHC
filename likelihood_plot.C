/*Plot the summary plots of the likelihood method*/
void likelihood_plot(){

    std::string baseDir = "/afs/cern.ch/work/r/ridz01/public/";

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

    TH2D * h_reco_energy_2D = new TH2D("Reco_energy_True_energy", "Energy_Reco versus Energy_True; energy_true; energy_reco", N_energy, &energy_bins[0], N_energy,&energy_bins[0]);
    TH2D * h_reco_Z_2D = new TH2D("Reco_Z_True_Z", "Z_Reco versus Z_True; Z_true; Z_reco", N_Z, &Z_bins[0], N_Z, &Z_bins[0]);


    std::vector<THStack*> h_reco_energy_1D;
    for (int j = 0; j < N_energy; j++ ){
        h_reco_energy_1D.push_back(new THStack(("Reco_energy_"+std::to_string(energy_test[j])).c_str(), " "));
    }

    for (int j = 0; j < N_energy; j++){
        for (int i_Z_true = 0; i_Z_true < N_Z; i_Z_true++){
            cout<<energy_test[j]<<"\t"<< Z_test[i_Z_true]<<endl;
            TFile * f_1d = new TFile((baseDir+"Reco_histograms_"+std::to_string(energy_test[j])+"_"+std::to_string(Z_test[i_Z_true])+".root").c_str(),"READ");
            TH1D * h_1d = (TH1D*)f_1d->Get("Reco_energy_test;1");
            
            //TH1D * temp = new TH1D(("Reco_energy_"+std::to_string(Z_test[i_Z_true])).c_str(), "Energy_Reco ; energy_reco", N_energy, &energy_bins[0]);
            for (int i = 1; i <= h_1d->GetNbinsX(); i++){
                h_reco_energy_2D->SetBinContent(i, j+1, h_reco_energy_2D->GetBinContent(j+1,i)+h_1d->GetBinContent(i));
            }
            h_1d->SetFillStyle(1001);
            h_1d->SetFillColor(40+i_Z_true);
            h_1d->SetStats(0);
            h_1d->SetTitle((std::to_string(Z_test[i_Z_true])+" cm").c_str());
            h_reco_energy_1D.at(j)->Add(h_1d);
            auto legend = new TLegend(0.6,0.7,0.9,0.9);
            legend->AddEntry(h_1d, (std::to_string(Z_test[i_Z_true])+" cm").c_str());
            legend->Draw();
        }
    }
    /*for (int i_Z_true = 0; i_Z_true < N_Z; i_Z_true++){
        TFile * f_1d = new TFile((baseDir+"Reco_histograms_2000_"+std::to_string(Z_test[i_Z_true])+".root").c_str(),"READ");
        TH1D * h_1d = (TH1D*)f_1d->Get("Reco_energy_test;1");

        for (int i = 0; i <= h_1d->GetNbinsX(); i++){
                h_reco_energy_2D->SetBinContent(i, 21, h_reco_energy_2D->GetBinContent(21,i)+h_1d->GetBinContent(i));
        }
    }*/

    for (int i_Z_true = 0; i_Z_true < N_Z; i_Z_true++){
        for (int j = 0; j < N_energy; j++){
            TFile * f_2d = new TFile((baseDir+"Reco_histograms_"+std::to_string(energy_test[j])+"_"+std::to_string(Z_test[i_Z_true])+".root").c_str(),"READ");
            TH1D * h_2d = (TH1D*)f_2d->Get("Reco_Z_test;1");
            for (int i = 1; i <= h_2d->GetNbinsX(); i++){
                h_reco_Z_2D->SetBinContent(i, i_Z_true+1, h_reco_Z_2D->GetBinContent(i_Z_true+1,i)+h_2d->GetBinContent(i));
            }
            /*for (int i = 0; i <= h_2d->GetNbinsX(); i++){
                h_reco_Z_2D->SetBinContent(i, 8, h_reco_Z_2D->GetBinContent(8,i)+h_2d->GetBinContent(i));
            }*/
       }
    }
    /*for (int j = 0; j < N_energy; j++){
            TFile * f_2d = new TFile((baseDir+"Reco_histograms_"+std::to_string(energy_test[j])+"_299.root").c_str(),"READ");
            TH1D * h_2d = (TH1D*)f_2d->Get("Reco_Z_test;1");
            for (int i = 0; i <= h_2d->GetNbinsX(); i++){
                h_reco_Z_2D->SetBinContent(i, 8, h_reco_Z_2D->GetBinContent(8,i)+h_2d->GetBinContent(i));
            }
    }*/


    TCanvas * c_1 = new TCanvas();
     h_reco_energy_2D->SetStats(0);
    h_reco_energy_2D->Draw("COLZ");
    c_1->SaveAs("/afs/cern.ch/work/r/ridz01/public/Reco_energy_true_energy.png");

    TCanvas * c_2 = new TCanvas();
    h_reco_Z_2D->SetStats(0);
    h_reco_Z_2D->Draw("COLZ");
    c_2->SaveAs("/afs/cern.ch/work/r/ridz01/public/Reco_Z_true_Z.png");

    /*TCanvas * c_2_1 = new TCanvas();
    h_reco_energy_1D.at(0)->Draw();
    c_2_1->SaveAs("Energy_Reco_Z_50GeV.png");

    TCanvas * c_2_2 = new TCanvas();
    h_reco_energy_1D.at(0)->Draw();
    c_2_2->SaveAs("Energy_Reco_Z_75GeV.png");

    TCanvas * c_2_3 = new TCanvas();
    h_reco_energy_1D.at(1)->Draw();
    c_2_3->SaveAs("Energy_Reco_Z_100GeV.png");

    TCanvas * c_2_3 = new TCanvas();
    h_reco_energy_1D.at(2)->Draw();
    c_2_3->SaveAs("Energy_Reco_Z_200GeV.png");

    TCanvas * c_2_4 = new TCanvas();
    h_reco_energy_1D.at(3)->Draw();
    c_2_4->SaveAs("Energy_Reco_Z_300GeV.png");

    TCanvas * c_2_5 = new TCanvas();
    h_reco_energy_1D.at(4)->Draw();
    c_2_5->SaveAs("Energy_Reco_Z_400GeV.png");

    TCanvas * c_2_6 = new TCanvas();
    h_reco_energy_1D.at(5)->Draw();
    c_2_6->SaveAs("Energy_Reco_Z_500GeV.png");*/

    /*TCanvas * c_2_7 = new TCanvas();
    h_reco_energy_1D.at(6)->Draw();
    c_2_7->SaveAs("Energy_Reco_Z_750GeV.png");

    TCanvas * c_2_8 = new TCanvas();
    h_reco_energy_1D.at(7)->Draw();
    c_2_8->SaveAs("Energy_Reco_Z_1000GeV.png");

    TCanvas * c_2_9 = new TCanvas();
    h_reco_energy_1D.at(8)->Draw();
    c_2_9->SaveAs("Energy_Reco_Z_1250GeV.png");

    TCanvas * c_2_10 = new TCanvas();
    h_reco_energy_1D.at(9)->Draw();
    c_2_10->SaveAs("Energy_Reco_Z_1500GeV.png");*/


    TFile * out_plots_2D = new TFile("/afs/cern.ch/work/r/ridz01/public/Reco_energy_True_energy.root", "RECREATE");
    h_reco_energy_2D->Write();
    h_reco_energy_1D.at(0)->Write();
    h_reco_energy_1D.at(1)->Write();
    h_reco_energy_1D.at(2)->Write();
    h_reco_energy_1D.at(3)->Write();
    h_reco_energy_1D.at(4)->Write();
    h_reco_energy_1D.at(5)->Write();
    h_reco_energy_1D.at(6)->Write();
    h_reco_energy_1D.at(7)->Write();
    h_reco_energy_1D.at(8)->Write();
    h_reco_energy_1D.at(9)->Write();
    h_reco_energy_1D.at(10)->Write();
    h_reco_energy_1D.at(11)->Write();
    h_reco_energy_1D.at(12)->Write();
    h_reco_energy_1D.at(13)->Write();
    h_reco_energy_1D.at(14)->Write();
    h_reco_energy_1D.at(15)->Write();
    h_reco_energy_1D.at(16)->Write();
    h_reco_energy_1D.at(17)->Write();
    h_reco_energy_1D.at(18)->Write();
    h_reco_energy_1D.at(19)->Write();
    h_reco_energy_1D.at(20)->Write();
    

}