#include <vector>
#include <string>
#include <iostream>
#include <fstream>


#include <TDirectoryFile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPad.h>

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

void draw_image(char *argv[]) {
    std::string filename = argv[1];
    std::string outfilename = argv[2];
    int Nsx = std::atoi(argv[3]);
    int Nsy = std::atoi(argv[4]);
    double sky_pitch = std::stod(argv[5]);
    std::cout << Nsx << " " << Nsy << " " << sky_pitch << std::endl;
    std::vector<std::vector<double>> S_dist = read_sky_parameters(filename, Nsx, Nsy);

    double xmin = -Nsx/2 * sky_pitch;
    double xmax = (Nsx-Nsx/2)* sky_pitch;
    double ymin = -Nsy/2 * sky_pitch;
    double ymax = (Nsy-Nsy/2)* sky_pitch;    

    TCanvas* canvas = new TCanvas("c", "c", 600, 600);
    TH2D *hist = new TH2D("image", ";#theta_{x} [arcsec];#theta_{y} [arcsec];#tilde{S}_{#theta_{x}, #theta_{y}}", Nsx, xmin, xmax, Nsy, ymin, ymax);
    for (int i=0; i<Nsx; i++){
        for (int j=0; j<Nsy; j++){
            double v = S_dist[i][j];
            double x = (i - Nsx/2) * sky_pitch;
            double y = (j - Nsy/2) * sky_pitch;
            // std::cout << x << " " << y << " " << v << std::endl;
            hist->Fill(x, y, v);
        }
    }

    const double zmax = hist->GetBinContent(hist->GetMaximumBin());
    const double zmin = 0;

    ((TGaxis*)hist->GetZaxis())->SetMaxDigits(2);
    hist->GetZaxis()->SetLabelSize(0.02);
    hist->GetYaxis()->SetNoExponent();
    hist->GetXaxis()->SetNoExponent();

    // const double zmin = zmax / 50.0;
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.15);
    gPad->SetBottomMargin(0.15);
    // canvas->cd(1);
    hist->GetZaxis()->SetRangeUser(zmin, zmax);
    hist->Draw("colz");
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
    if (argc == 7) {
        draw_polarization_image(argv);
    } else if (argc == 6) {
        draw_image(argv);
    }
}