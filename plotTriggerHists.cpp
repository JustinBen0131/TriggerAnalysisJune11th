#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <vector>

// Function to recursively save histograms as PNGs
void SaveHistograms(TDirectory* dir, const std::string& outputDir) {
    // Create the output directory if it doesn't exist
    gSystem->mkdir(outputDir.c_str(), true);

    // Get a list of keys in the current directory
    TIter next(dir->GetListOfKeys());
    TKey* key;
    
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        
        // If the object is a directory, recurse into it
        if (obj->InheritsFrom(TDirectory::Class())) {
            TDirectory* subdir = (TDirectory*)obj;
            std::string subDirName = subdir->GetName();
            std::string subDirOutputPath = outputDir + "/" + subDirName;
            SaveHistograms(subdir, subDirOutputPath);
        }
        // If the object is a 1D histogram, save it as a PNG
        else if (obj->InheritsFrom(TH1::Class())) {
            TH1* hist = (TH1*)obj;
            std::string histName = hist->GetName();
            std::string histOutputPath = outputDir + "/" + histName + ".png";
            TCanvas canvas;
            hist->GetYaxis()->SetTitle("dN/(dE*dt)");
            // Set the x-axis label based on histogram name
            if (histName.find("8x8") != std::string::npos) {
                hist->GetXaxis()->SetTitle("Maximum 8by8 Tower Energy Sum");
            } else {
                hist->GetXaxis()->SetTitle("Tower Energy");
            }
            
            hist->Draw();
            canvas.SaveAs(histOutputPath.c_str());
        }
    }
}
void OverlayLiveHistograms(TDirectory* liveDir, const std::string& outputDir) {
    // Define the histogram names and legend labels
    std::vector<std::string> histNames = {
        "trigger_28_live_max8x8",
        "trigger_29_live_max8x8",
        "trigger_30_live_max8x8",
        "trigger_31_live_max8x8"
    };
    std::vector<std::string> legendLabels = {
        "Photon 1",
        "Photon 2",
        "Photon 3",
        "Photon 4"
    };

    // Colors for the histograms
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta};

    TCanvas canvas;
    TLegend legend(0.7, 0.7, 0.9, 0.9);
    
    bool firstHist = true;
    double maxY = 0;

    // First loop to determine the maximum y-axis range
    for (const auto& histName : histNames) {
        TH1* hist = (TH1*)liveDir->Get(histName.c_str());
        if (hist) {
            std::cout << "Found histogram: " << histName << " with max Y value: " << hist->GetMaximum() << std::endl;
            double currentMaxY = hist->GetMaximum();
            if (currentMaxY > maxY) {
                maxY = currentMaxY;
            }
        } else {
            std::cout << "Histogram not found: " << histName << std::endl;
        }
    }
    // Second loop to draw the histograms
    for (size_t i = 0; i < histNames.size(); ++i) {
        TH1* hist = (TH1*)liveDir->Get(histNames[i].c_str());
        if (hist) {
            hist->SetLineColor(colors[i]);
            hist->SetLineWidth(2);
            hist->GetXaxis()->SetTitle("Maximum 8by8 Tower Energy Sum");
            hist->SetMaximum(maxY * 1.1);  // Set a bit higher than the maximum value to ensure visibility
            hist->SetStats(0);  // Turn off statistics box
            
            if (firstHist) {
                std::cout << "Drawing first histogram: " << histNames[i] << std::endl;
                hist->Draw("HIST");
                firstHist = false;
            } else {
                std::cout << "Drawing histogram: " << histNames[i] << std::endl;
                hist->Draw("HIST SAME");
            }
            
            legend.AddEntry(hist, legendLabels[i].c_str(), "l");
        } else {
            std::cout << "Histogram not found in second loop: " << histNames[i] << std::endl;
        }
    }

    legend.Draw();
    std::string overlayOutputPath = outputDir + "/overlay_live_max8x8.png";
    canvas.SaveAs(overlayOutputPath.c_str());
}


// Main macro function
void plotTriggerHists(const char* rootFilePath, const char* outputBaseDir) {
    // Open the root file
    TFile* file = TFile::Open(rootFilePath, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << rootFilePath << std::endl;
        return;
    }
    
    // Process each top-level directory
    std::vector<std::string> topDirs = {"raw", "live", "scaled"};
    for (const std::string& dirName : topDirs) {
        TDirectory* dir = (TDirectory*)file->Get(dirName.c_str());
        if (dir) {
            std::string outputDir = std::string(outputBaseDir) + "/" + dirName;
            SaveHistograms(dir, outputDir);

            // Overlay histograms if in the 'live' directory
            if (dirName == "live") {
                OverlayLiveHistograms(dir, outputDir);
            }
        } else {
            std::cerr << "Warning: Directory " << dirName << " not found in file " << rootFilePath << std::endl;
        }
    }
    
    // Close the file
    file->Close();
}
