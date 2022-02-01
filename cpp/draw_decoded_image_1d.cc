#include <vector>
#include <string>
#include <iostream>
#include <fstream>


#include <TDirectoryFile.h>
#include <TLegend.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>

// typedef std::string<std::string<double>> S;

std::vector<std::vector<double>> read_sky_parameters(std::string filename, int Nsx, int Nsy) {
    std::ifstream infile;

    infile.open(filename);
    if (!infile) {
        std::cerr << " cannot open file: " << filename << std::endl;
        throw std::runtime_error("invalid file");
    }
    std::vector<std::vector<double>> ret(Nsx, std::vector<double>(Nsy, 0.0));
    double v;
    int x = 0;
    int y = 0;
    while (infile >> v) {
        ret[Nsx-x-1][y] = v;
        x++;
        if (x == Nsx) {
            x = 0;
            y++;
        }
    }
    infile.close();

    return ret;
}

void draw_image(int argc, char *argv[]) {
    int Nsx = std::atoi(argv[1]);
    int Nsy = std::atoi(argv[2]);
    double sky_pitch = std::stod(argv[3]);
    std::string outfilename = argv[4];

    std::vector<std::vector<std::vector<double>>> S_dists;
    std::vector<std::string> file_ids;
    for (int i=5; i<argc; i+=2){
        std::string filename = argv[i];
        std::string file_id = argv[i+1];
        std::vector<std::vector<double>> S_dist = read_sky_parameters(filename, Nsx, Nsy);
        S_dists.push_back(S_dist);
        file_ids.push_back(file_id);
    }
    std::cout << Nsx << " " << Nsy << " " << sky_pitch << std::endl;
    

    double xmin = -Nsx/2 * sky_pitch;
    double xmax = (Nsx-Nsx/2)* sky_pitch;
    double ymin = -Nsy/2 * sky_pitch;
    double ymax = (Nsy-Nsy/2)* sky_pitch;

    TCanvas* canvas = new TCanvas("c", "c", 800, 600);
    TLegend* leg = new TLegend(0.72,0.72,0.9,0.9);
    // TH2D *hist = new TH2D("image", "", Nsx, xmin, xmax, Nsy, ymin, ymax);
    for (int seq=0; seq<file_ids.size(); seq++){
        std::string file_id = file_ids[seq];

        TH1D *hist = new TH1D(file_id.c_str(), ";#theta_{x} [arcsec];#tilde{S}_{#theta_{x}, 0}", Nsx, xmin, xmax);
        for (int i=0; i<Nsx; i++){
            int j = Nsy/2 - 1;
            double v = S_dists[seq][i][j];
            double x = (i - Nsx/2) * sky_pitch;
            double y = (j - Nsy/2) * sky_pitch;
            // std::cout << x << " " << y << " " << v << std::endl;
            hist->Fill(x, v);
        }

        double zmax = hist->GetBinContent(hist->GetMaximumBin()) * 1.5;
        double zmin = 0;

        ((TGaxis*)hist->GetYaxis())->SetMaxDigits(2);
        // hist->GetZaxis()->SetLabelSize(0.02);
        // hist->GetYaxis()->SetNoExponent();
        hist->GetXaxis()->SetNoExponent();
        hist->SetLineColor(seq+1);

        // const double zmin = zmax / 50.0;
        gStyle->SetOptStat(0);
        gStyle->SetPalette(56);
        // canvas->cd(1);
        hist->GetYaxis()->SetRangeUser(zmin, zmax);
        leg->AddEntry(hist, file_id.c_str(), "l");
        hist->Draw("hist;same");
    
    }
    leg->Draw("same");
    canvas->SaveAs(outfilename.c_str());
}

void draw_polarization_image(char *argv[]) {
    std::string filename_Htype = argv[1];
    std::string filename_Vtype = argv[2];
    std::string outfilename = argv[3];
    int Nsx = std::atoi(argv[4]);
    int Nsy = std::atoi(argv[5]);
    double sky_pitch = std::stod(argv[6]);
    std::cout << Nsx << " " << Nsy << " " << sky_pitch << std::endl;
    std::vector<std::vector<double>> S_dist_Htype = read_sky_parameters(filename_Htype, Nsx, Nsy);
    std::vector<std::vector<double>> S_dist_Vtype = read_sky_parameters(filename_Vtype, Nsx, Nsy);

    double xmin = -Nsx/2 * sky_pitch;
    double xmax = (Nsx-Nsx/2)* sky_pitch;
    double ymin = -Nsy/2 * sky_pitch;
    double ymax = (Nsy-Nsy/2)* sky_pitch;    

    TCanvas* canvas = new TCanvas("c", "c", 1200, 600);
    canvas->Divide(2, 1);
    TH2D *hist_Htype = new TH2D("Htype", "Htype", Nsx, xmin, xmax, Nsy, ymin, ymax);
    for (int i=0; i<Nsx; i++){
        for (int j=0; j<Nsy; j++){
            double v = S_dist_Htype[i][j];
            double x = (i - Nsx/2) * sky_pitch;
            double y = (j - Nsy/2) * sky_pitch;
            // std::cout << x << " " << y << " " << v << std::endl;
            hist_Htype->Fill(x, y, v);
        }
    }
    TH2D *hist_Vtype = new TH2D("Vtype", "Vtype", Nsx, xmin, xmax, Nsy, ymin, ymax);
    for (int i=0; i<Nsx; i++){
        for (int j=0; j<Nsy; j++){
            double v = S_dist_Vtype[i][j];
            double x = (i - Nsx/2) * sky_pitch;
            double y = (j - Nsy/2) * sky_pitch;
            hist_Vtype->Fill(x, y, v);
        }
    }

    const double zmax_Htype = hist_Htype->GetBinContent(hist_Htype->GetMaximumBin());
    const double zmax_Vtype = hist_Vtype->GetBinContent(hist_Vtype->GetMaximumBin());
    const double zmax = std::max(zmax_Htype, zmax_Vtype);
    const double zmin = 0;

    ((TGaxis*)hist_Htype->GetZaxis())->SetMaxDigits(2);
    hist_Htype->GetZaxis()->SetLabelSize(0.02);
    hist_Htype->GetYaxis()->SetNoExponent();
    hist_Htype->GetXaxis()->SetNoExponent();

    ((TGaxis*)hist_Vtype->GetZaxis())->SetMaxDigits(2);
    hist_Vtype->GetZaxis()->SetLabelSize(0.02);
    hist_Vtype->GetYaxis()->SetNoExponent();
    hist_Vtype->GetXaxis()->SetNoExponent();
    // const double zmin = zmax / 50.0;
    gStyle->SetOptStat(0);
    canvas->cd(1);
    hist_Htype->GetZaxis()->SetRangeUser(zmin, zmax);
    hist_Htype->Draw("colz");
    canvas->cd(2);
    hist_Vtype->GetZaxis()->SetRangeUser(zmin, zmax);
    hist_Vtype->Draw("colz");
    canvas->SaveAs(outfilename.c_str());
}

int main(int argc, char *argv[]) {
    draw_image(argc, argv);
}