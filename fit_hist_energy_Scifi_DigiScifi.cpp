/*Script to plot the TGraphs of energy v/s the no of hits for the digitised Scifi for each station*/
#include <fstream>
#include <cmath>
void fit_hist_energy_Scifi_DigiScifi(){

    // READING THE BIN EDGES FROM THE .txt FILES
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

    ifstream input_file("/eos/user/r/ridz01/EnergyBins_0_200_Digi.txt");
    while (!input_file.eof()) {
        int tmp;
        input_file >> tmp;
        bin_0.push_back(tmp);
        ctr_0++;  
    }

    ifstream input_file_1("/eos/user/r/ridz01/EnergyBins_1_200_Digi.txt");
    while (!input_file_1.eof()) {
        int tmp;
        input_file_1 >> tmp;
        bin_1.push_back(tmp);
        ctr_1++;  
    }

    ifstream input_file_2("/eos/user/r/ridz01/EnergyBins_2_200_Digi.txt");
    while (!input_file_2.eof()) {
        int tmp;
        input_file_2 >> tmp;
        bin_2.push_back(tmp);
        ctr_2++;  
    }

    ifstream input_file_3("/eos/user/r/ridz01/EnergyBins_3_200_Digi.txt");
    while (!input_file_3.eof()) {
        int tmp;
        input_file_3 >> tmp;
        bin_3.push_back(tmp);
        ctr_3++;  
    }

    ifstream input_file_4("/eos/user/r/ridz01/EnergyBins_4_200_Digi.txt");
    while (!input_file_4.eof()) {
        int tmp;
        input_file_4 >> tmp;
        bin_4.push_back(tmp);
        ctr_4++;  
    }
    
    // Creating Histogram for Scifihits NC events
    std::vector<TH1D*> hist_vec_scifi_NC_00;
    std::vector<TH1D*> hist_vec_scifi_NC_01;
    std::vector<TH1D*> hist_vec_scifi_NC_02;
    std::vector<TH1D*> hist_vec_scifi_NC_03;
    std::vector<TH1D*> hist_vec_scifi_NC_04;
    TFile* fin2 = new TFile("/eos/user/r/ridz01/Hist_SciFi_NC_energy_DigiScifi_irregular_231122_200files.root","READ");
    for (int i = 0; i<(ctr_0)-1; i++){
        hist_vec_scifi_NC_00.push_back((TH1D*)fin2->Get(("h_Scifi_NC_en_00_"+std::to_string(i)+";1").c_str()));
    }
    for (int i = 0; i<(ctr_1)-1; i++){
        hist_vec_scifi_NC_01.push_back((TH1D*)fin2->Get(("h_Scifi_NC_en_01_"+std::to_string(i)+";1").c_str()));
    } 
    for (int i = 0; i<(ctr_2)-1; i++){
         hist_vec_scifi_NC_02.push_back((TH1D*)fin2->Get(("h_Scifi_NC_en_02_"+std::to_string(i)+";1").c_str()));
    } 
    for (int i = 0; i<(ctr_3)-1; i++){
        hist_vec_scifi_NC_03.push_back((TH1D*)fin2->Get(("h_Scifi_NC_en_03_"+std::to_string(i)+";1").c_str()));
    }
    for (int i = 0; i<(ctr_4)-1; i++){
       hist_vec_scifi_NC_04.push_back((TH1D*)fin2->Get(("h_Scifi_NC_en_04_"+std::to_string(i)+";1").c_str()));
    }

    int numbin_0 =ctr_0-1, numbin_1= ctr_1-1, numbin_2= ctr_2-1, numbin_3 = ctr_3-1, numbin_4 = ctr_4-1;
    float x_i_Sc_NC_00[numbin_0],err_x_i_Sc_NC_00[numbin_0],err_y_i_Sc_NC_00[numbin_0];
    float mean_hist_Sc_NC_00[numbin_0], sig_hist_Sc_NC_00[numbin_0];
    for (int i_hist=0; i_hist<numbin_0; i_hist++){
        x_i_Sc_NC_00[i_hist] = (bin_0.at(i_hist) + bin_0.at(i_hist+1))/2;
        err_x_i_Sc_NC_00[i_hist]=(bin_0.at(i_hist+1) - bin_0.at(i_hist))/2;
        mean_hist_Sc_NC_00[i_hist]= hist_vec_scifi_NC_00.at(i_hist)->GetMean();
        sig_hist_Sc_NC_00[i_hist]=hist_vec_scifi_NC_00.at(i_hist)->GetRMS();
        err_y_i_Sc_NC_00[i_hist]=sig_hist_Sc_NC_00[i_hist]/ sqrt(hist_vec_scifi_NC_00.at(i_hist)->GetEntries());
        cout<<i_hist<<" "<<x_i_Sc_NC_00[i_hist]<<" "<<err_x_i_Sc_NC_00[i_hist]<<" "<<mean_hist_Sc_NC_00[i_hist]<<" "<< sig_hist_Sc_NC_00[i_hist]<<endl;
    }

    float x_i_Sc_NC_01[numbin_1],err_x_i_Sc_NC_01[numbin_1],err_y_i_Sc_NC_01[numbin_1];
    float mean_hist_Sc_NC_01[numbin_1], sig_hist_Sc_NC_01[numbin_1];
    for (int i_hist=0; i_hist<numbin_1; i_hist++){
        x_i_Sc_NC_01[i_hist] = (bin_1.at(i_hist) + bin_1.at(i_hist+1))/2;
        err_x_i_Sc_NC_01[i_hist]=(bin_1.at(i_hist+1) - bin_1.at(i_hist))/2;
        mean_hist_Sc_NC_01[i_hist]= hist_vec_scifi_NC_01.at(i_hist)->GetMean();
        sig_hist_Sc_NC_01[i_hist]=hist_vec_scifi_NC_01.at(i_hist)->GetRMS();
        err_y_i_Sc_NC_01[i_hist]=sig_hist_Sc_NC_01[i_hist]/ sqrt(hist_vec_scifi_NC_01.at(i_hist)->GetEntries());
        cout<<i_hist<<" "<<x_i_Sc_NC_01[i_hist]<<" "<<err_x_i_Sc_NC_01[i_hist]<<" "<<mean_hist_Sc_NC_01[i_hist]<<" "<< sig_hist_Sc_NC_01[i_hist]<<endl;
    }

    float x_i_Sc_NC_02[numbin_2],err_x_i_Sc_NC_02[numbin_2],err_y_i_Sc_NC_02[numbin_2];
    float mean_hist_Sc_NC_02[numbin_2], sig_hist_Sc_NC_02[numbin_2];
    for (int i_hist=0; i_hist<numbin_2; i_hist++){
        x_i_Sc_NC_02[i_hist] = (bin_2.at(i_hist) + bin_2.at(i_hist+1))/2;
        err_x_i_Sc_NC_02[i_hist]=(bin_2.at(i_hist+1) - bin_2.at(i_hist))/2;
        mean_hist_Sc_NC_02[i_hist]= hist_vec_scifi_NC_02.at(i_hist)->GetMean();
        sig_hist_Sc_NC_02[i_hist]=hist_vec_scifi_NC_02.at(i_hist)->GetRMS();
        err_y_i_Sc_NC_02[i_hist]=sig_hist_Sc_NC_02[i_hist]/ sqrt(hist_vec_scifi_NC_02.at(i_hist)->GetEntries());
        cout<<i_hist<<" "<<x_i_Sc_NC_02[i_hist]<<" "<<err_x_i_Sc_NC_02[i_hist]<<" "<<mean_hist_Sc_NC_02[i_hist]<<" "<< sig_hist_Sc_NC_02[i_hist]<<endl;
    }

    float x_i_Sc_NC_03[numbin_3],err_x_i_Sc_NC_03[numbin_3],err_y_i_Sc_NC_03[numbin_3];
    float mean_hist_Sc_NC_03[numbin_3], sig_hist_Sc_NC_03[numbin_3];
    for (int i_hist=0; i_hist<numbin_3; i_hist++){
        x_i_Sc_NC_03[i_hist] = (bin_3.at(i_hist) + bin_3.at(i_hist+1))/2;
        err_x_i_Sc_NC_03[i_hist]=(bin_3.at(i_hist+1) - bin_3.at(i_hist))/2;
        mean_hist_Sc_NC_03[i_hist]= hist_vec_scifi_NC_03.at(i_hist)->GetMean();
        sig_hist_Sc_NC_03[i_hist]=hist_vec_scifi_NC_03.at(i_hist)->GetRMS();
        err_y_i_Sc_NC_03[i_hist]=sig_hist_Sc_NC_03[i_hist]/ sqrt(hist_vec_scifi_NC_03.at(i_hist)->GetEntries());
        cout<<i_hist<<" "<<x_i_Sc_NC_03[i_hist]<<" "<<err_x_i_Sc_NC_03[i_hist]<<" "<<mean_hist_Sc_NC_03[i_hist]<<" "<< sig_hist_Sc_NC_03[i_hist]<<endl;
    }

    float x_i_Sc_NC_04[numbin_4],err_x_i_Sc_NC_04[numbin_4],err_y_i_Sc_NC_04[numbin_4];
    float mean_hist_Sc_NC_04[numbin_4], sig_hist_Sc_NC_04[numbin_4];
    for (int i_hist=0; i_hist<numbin_4; i_hist++){
        x_i_Sc_NC_04[i_hist] = (bin_4.at(i_hist) + bin_4.at(i_hist+1))/2;
        err_x_i_Sc_NC_04[i_hist]=(bin_4.at(i_hist+1) - bin_4.at(i_hist))/2;
        mean_hist_Sc_NC_04[i_hist]= hist_vec_scifi_NC_04.at(i_hist)->GetMean();
        sig_hist_Sc_NC_04[i_hist]=hist_vec_scifi_NC_04.at(i_hist)->GetRMS();
        err_y_i_Sc_NC_04[i_hist]=sig_hist_Sc_NC_04[i_hist]/sqrt(hist_vec_scifi_NC_04.at(i_hist)->GetEntries());
        cout<<i_hist<<" "<<x_i_Sc_NC_04[i_hist]<<" "<<err_x_i_Sc_NC_04[i_hist]<<" "<<mean_hist_Sc_NC_04[i_hist]<<" "<< sig_hist_Sc_NC_04[i_hist]<<endl;
    }
     
    
    auto Sci_NC_gr_hist_00_err = new TGraphErrors(numbin_0,x_i_Sc_NC_00,mean_hist_Sc_NC_00,err_x_i_Sc_NC_00,err_y_i_Sc_NC_00);
    auto Sci_NC_gr_hist_01_err = new TGraphErrors(numbin_1,x_i_Sc_NC_01,mean_hist_Sc_NC_01,err_x_i_Sc_NC_01,err_y_i_Sc_NC_01);
    auto Sci_NC_gr_hist_02_err = new TGraphErrors(numbin_2,x_i_Sc_NC_02,mean_hist_Sc_NC_02,err_x_i_Sc_NC_02,err_y_i_Sc_NC_02);
    auto Sci_NC_gr_hist_03_err = new TGraphErrors(numbin_3,x_i_Sc_NC_03,mean_hist_Sc_NC_03,err_x_i_Sc_NC_03,err_y_i_Sc_NC_03);
    auto Sci_NC_gr_hist_04_err = new TGraphErrors(numbin_4,x_i_Sc_NC_04,mean_hist_Sc_NC_04,err_x_i_Sc_NC_04,err_y_i_Sc_NC_04);

    auto Sci_NC_gr_hist_00_sig = new TGraphErrors(numbin_0,x_i_Sc_NC_00,mean_hist_Sc_NC_00,err_x_i_Sc_NC_00,sig_hist_Sc_NC_00);
    auto Sci_NC_gr_hist_01_sig = new TGraphErrors(numbin_1,x_i_Sc_NC_01,mean_hist_Sc_NC_01,err_x_i_Sc_NC_01,sig_hist_Sc_NC_01);
    auto Sci_NC_gr_hist_02_sig = new TGraphErrors(numbin_2,x_i_Sc_NC_02,mean_hist_Sc_NC_02,err_x_i_Sc_NC_02,sig_hist_Sc_NC_02);
    auto Sci_NC_gr_hist_03_sig = new TGraphErrors(numbin_3,x_i_Sc_NC_03,mean_hist_Sc_NC_03,err_x_i_Sc_NC_03,sig_hist_Sc_NC_03);
    auto Sci_NC_gr_hist_04_sig = new TGraphErrors(numbin_4,x_i_Sc_NC_04,mean_hist_Sc_NC_04,err_x_i_Sc_NC_04,sig_hist_Sc_NC_04);

    auto Sci_NC_err_hist_00 = new TGraph(numbin_0,x_i_Sc_NC_00,sig_hist_Sc_NC_00);
    auto Sci_NC_err_hist_01 = new TGraph(numbin_1,x_i_Sc_NC_01,sig_hist_Sc_NC_01);
    auto Sci_NC_err_hist_02 = new TGraph(numbin_2,x_i_Sc_NC_02,sig_hist_Sc_NC_02);
    auto Sci_NC_err_hist_03 = new TGraph(numbin_3,x_i_Sc_NC_03,sig_hist_Sc_NC_03);
    auto Sci_NC_err_hist_04 = new TGraph(numbin_4,x_i_Sc_NC_04,sig_hist_Sc_NC_04);

    Sci_NC_gr_hist_00_sig->SetTitle("Station 0; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_01_sig->SetTitle("Station 1; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_02_sig->SetTitle("Station 2; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_03_sig->SetTitle("Station 3; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_04_sig->SetTitle("Station 4; No of hits; Energy[GeV]");

    Sci_NC_gr_hist_00_err->SetTitle("Station 0; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_01_err->SetTitle("Station 1; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_02_err->SetTitle("Station 2; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_03_err->SetTitle("Station 3; No of hits; Energy[GeV]");
    Sci_NC_gr_hist_04_err->SetTitle("Station 4; No of hits; Energy[GeV]");

    ofstream file;
    file.open ("/eos/user/r/ridz01/TGraph/TGraph_DigiScifi_0_irregularbin_231122.txt");

    file<<"First points for each station:"<<endl;

    file<< "Station 0 point X=0,Y="<< Sci_NC_gr_hist_00_sig->Eval(0)<<"Y error="<<sig_hist_Sc_NC_00[0]<<"  "<<err_y_i_Sc_NC_00[0]<<endl;
    file<< "Station 1 point X=0,Y="<< Sci_NC_gr_hist_01_sig->Eval(0)<<"Y error="<<sig_hist_Sc_NC_01[0]<<"  "<<err_y_i_Sc_NC_01[0]<<endl;
    file<< "Station 2 point X=0,Y="<< Sci_NC_gr_hist_02_sig->Eval(0)<<"Y error="<<sig_hist_Sc_NC_02[0]<<"  "<<err_y_i_Sc_NC_02[0]<<endl;
    file<< "Station 3 point X=0,Y="<< Sci_NC_gr_hist_03_sig->Eval(0)<<"Y error="<<sig_hist_Sc_NC_03[0]<<"  "<<err_y_i_Sc_NC_03[0]<<endl;
    file<< "Station 4 point X=0,Y="<< Sci_NC_gr_hist_04_sig->Eval(0)<<"Y error="<<sig_hist_Sc_NC_04[0]<<"  "<<err_y_i_Sc_NC_04[0]<<endl;

    file<<"Correlation factor of the no of hits v/s energy for each station:"<<endl;

    file<< "Station 0 Correlation factor = "<< Sci_NC_gr_hist_00_sig->GetCorrelationFactor()<<endl;
    file<< "Station 1 Correlation factor ="<< Sci_NC_gr_hist_01_sig->GetCorrelationFactor()<<endl;
    file<< "Station 2 Correlation factor ="<< Sci_NC_gr_hist_02_sig->GetCorrelationFactor()<<endl;
    file<< "Station 3 Correlation factor ="<< Sci_NC_gr_hist_03_sig->GetCorrelationFactor()<<endl;
    file<< "Station 4 Correlation factor ="<< Sci_NC_gr_hist_04_sig->GetCorrelationFactor()<<endl;

    file.close();

    TFile* fout2 = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Scifi_NC_energy_DigiScifi_irregular_231122.root","RECREATE");    
    Sci_NC_gr_hist_00_err->Write("Station_0");
    Sci_NC_gr_hist_00_sig->Write("Station_0_sig");
    Sci_NC_gr_hist_01_err->Write("Station_1");
    Sci_NC_gr_hist_01_sig->Write("Station_1_sig");
    Sci_NC_gr_hist_02_err->Write("Station_2");
    Sci_NC_gr_hist_02_sig->Write("Station_2_sig");
    Sci_NC_gr_hist_03_err->Write("Station_3");
    Sci_NC_gr_hist_03_sig->Write("Station_3_sig");
    Sci_NC_gr_hist_04_err->Write("Station_4");
    Sci_NC_gr_hist_04_sig->Write("Station_4_sig");
    fout2->Close();

    TFile* fout2_err = new TFile("/eos/user/r/ridz01/TGraph/TGraph_Scifi_NC_energy_errors_DigiScifi_irregular_231122.root","RECREATE");
    Sci_NC_err_hist_00->Write("Station_0");
    Sci_NC_err_hist_01->Write("Station_1");
    Sci_NC_err_hist_02->Write("Station_2");
    Sci_NC_err_hist_03->Write("Station_3");
    Sci_NC_err_hist_04->Write("Station_4");
    fout2_err->Close();

    auto c0 = new TCanvas();
    Sci_NC_gr_hist_00_err->Draw("AP");
    c0->SaveAs("TGraph_00_02_DigiScifi.png");

    auto c1 = new TCanvas();
    Sci_NC_gr_hist_01_err->Draw("AP");
    c1->SaveAs("TGraph_01_02_DigiScifi.png");

    auto c2 = new TCanvas();
    Sci_NC_gr_hist_02_err->Draw("AP");
    c2->SaveAs("TGraph_02_02_DigiScifi.png");

    auto c3 = new TCanvas();
    Sci_NC_gr_hist_03_err->Draw("AP");
    c3->SaveAs("TGraph_03_02_DigiScifi.png");

    auto c4 = new TCanvas();
    Sci_NC_gr_hist_04_err->Draw("AP");
    c4->SaveAs("TGraph_04_02_DigiScifi.png");

}
