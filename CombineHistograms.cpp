#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TSystem.h>
#include <TDirectory.h>

void CombineHistograms(const char* directoryPath) {
    // Set the output file name
    TString outputPath = TString(directoryPath) + "/HistCombined44686.root";

    // Open the output file
    TFile *outputFile = new TFile(outputPath, "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cerr << "Failed to create output file!" << std::endl;
        return;
    }

    // List all files in the specified directory ending with ".root"
    void *dirp = gSystem->OpenDirectory(directoryPath);
    const char *filename;
    std::map<std::string, TH1*> histogramMap;

    while ((filename = gSystem->GetDirEntry(dirp))) {
        TString filepath = TString(directoryPath) + "/" + TString(filename);
        if (!filepath.EndsWith(".root")) continue;

        TFile *file = TFile::Open(filepath);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            continue;
        }

        std::cout << "Processing file: " << filepath << std::endl;
        TList *list = file->GetListOfKeys();
        TIter iter(list);
        TKey *key;

        while ((key = (TKey*)iter())) {
            TObject *obj = key->ReadObj();
            if (obj->InheritsFrom("TH1")) {
                TH1 *hist = (TH1*)obj;
                std::cout << "Processing histogram: " << hist->GetName() << std::endl;

                // Only process histograms with entries
                if (hist->GetEntries() > 0) {
                    if (histogramMap.find(hist->GetName()) == histogramMap.end()) {
                        histogramMap[hist->GetName()] = (TH1*)hist->Clone();
                        histogramMap[hist->GetName()]->SetDirectory(0); // Remove from current file directory
                    } else {
                        histogramMap[hist->GetName()]->Add(hist);
                    }
                }
            }
        }
        file->Close();
    }

    // Write combined histograms to the output file
    outputFile->cd();
    for (auto &entry : histogramMap) {
        if (entry.second->GetEntries() > 0) { // Ensure it has entries before saving
            std::cout << "Writing histogram: " << entry.first << " with " << entry.second->GetEntries() << " entries." << std::endl;
            entry.second->Write();
        }
    }

    outputFile->Close();
    std::cout << "Histograms combined and saved to " << outputPath << std::endl;
}
