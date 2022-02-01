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
#include <TLatex.h>
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


void draw_polarization_image(char *argv[]) {
    std::string filename_Htype = argv[1];
    std::string filename_Vtype = argv[2];
    std::string outfilename_I = argv[3];
    std::string outfilename_P = argv[4];
    int Nsx = std::atoi(argv[5]);
    int Nsy = std::atoi(argv[6]);
    double sky_pitch = std::stod(argv[7]);
    double modulation_factor = std::stod(argv[8]);
    std::cout << Nsx << " " << Nsy << " " << sky_pitch << std::endl;
    std::vector<std::vector<double>> S_dist_Htype = read_sky_parameters(filename_Htype, Nsx, Nsy);
    std::vector<std::vector<double>> S_dist_Vtype = read_sky_parameters(filename_Vtype, Nsx, Nsy);

    double xmin = -Nsx/2 * sky_pitch;
    double xmax = (Nsx-Nsx/2)* sky_pitch;
    double ymin = -Nsy/2 * sky_pitch;
    double ymax = (Nsy-Nsy/2)* sky_pitch;    

    TCanvas* cI = new TCanvas("cI", "cI", 600, 600);
    TCanvas* cP = new TCanvas("cP", "cP", 600, 600);
    TH2D *hist_I = new TH2D("I", ";#theta_{x} [arcsec];#theta_{y} [arcsec];#tilde{I}_{#theta_{x}, #theta_{y}}", Nsx, xmin, xmax, Nsy, ymin, ymax);
    TH2D *hist_P = new TH2D("P", ";#theta_{x} [arcsec];#theta_{y} [arcsec];#tilde{#Pi}^{(Q)}_{#theta_{x}, #theta_{y}}", Nsx, xmin, xmax, Nsy, ymin, ymax);

    for (int i=0; i<Nsx; i++){
        for (int j=0; j<Nsy; j++){
            double sh = S_dist_Htype[i][j];
            double sv = S_dist_Vtype[i][j];
            double x = (i - Nsx/2) * sky_pitch;
            double y = (j - Nsy/2) * sky_pitch;
            // std::cout << x << " " << y << " " << v << std::endl;
            hist_I->Fill(x, y, sh+sv);
            hist_P->Fill(x, y, (sh-sv)/(sh+sv)/modulation_factor);
        }
    }

    cI->cd();
    const double zmax = hist_I->GetBinContent(hist_I->GetMaximumBin());
    const double zmin = 0;

    ((TGaxis*)hist_I->GetZaxis())->SetMaxDigits(2);
    hist_I->GetZaxis()->SetLabelSize(0.02);
    hist_I->GetYaxis()->SetNoExponent();
    hist_I->GetXaxis()->SetNoExponent();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hist_I->GetZaxis()->SetRangeUser(zmin, zmax);
    hist_I->Draw("colz");
    cI->SaveAs(outfilename_I.c_str());

    cP->cd();
    ((TGaxis*)hist_P->GetZaxis())->SetMaxDigits(2);
    hist_P->GetZaxis()->SetLabelSize(0.02);
    hist_P->GetYaxis()->SetNoExponent();
    hist_P->GetXaxis()->SetNoExponent();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hist_P->GetZaxis()->SetRangeUser(-1.2, 1.2);
    hist_P->Draw("colz");
    cP->SaveAs(outfilename_P.c_str());
}

int main(int argc, char *argv[]) {
    draw_polarization_image(argv);
}