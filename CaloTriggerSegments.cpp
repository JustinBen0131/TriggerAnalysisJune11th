#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <bitset>
#include "TriggerDefs.h"

R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libtriggervalid.so)

/*
make list file: ./createListFile /sphenix/u/patsfan753/scratch/analysis/calotriggeremulator/output/44686 44686
 */

void CaloTriggerSegments(const char* filename, const char* output_filename, Long64_t nEvents = 0) {
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    std::cout << "ROOT file opened successfully: " << filename << std::endl;

    TTree *tree = (TTree*)file->Get("ttree");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'ttree' in the file." << std::endl;
        file->Close();
        return;
    }
    std::cout << "TTree 'ttree' loaded successfully." << std::endl;

    ULong64_t gl1_scaled[64], gl1_live[64];
    ULong64_t gl1_scaledvec;
    std::vector<float> *emcal_energy = nullptr;
    std::vector<float> *cluster_ecore = nullptr;

    tree->SetBranchAddress("gl1_scaled", gl1_scaled);
    tree->SetBranchAddress("gl1_scaledvec", &gl1_scaledvec);
    tree->SetBranchAddress("gl1_live", gl1_live);
    tree->SetBranchAddress("emcal_energy", &emcal_energy);
    tree->SetBranchAddress("cluster_ecore", &cluster_ecore);
    std::cout << "Branch addresses set successfully." << std::endl;

    TH1F *histograms_ecore[64];
    TH1F *histograms_energy[64];
    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i] = new TH1F(Form("trigger_%d_ecore", i), Form("Trigger %d Maximum Cluster ECore", i), 100, 0, 20);
        histograms_energy[i] = new TH1F(Form("trigger_%d_energy", i), Form("Trigger %d Maximum EMCAL Energy", i), 100, 0, 20);
    }
    std::cout << "Histograms initialized successfully." << std::endl;

    Long64_t nentries = tree->GetEntries();
    nentries = std::min(nentries, (nEvents > 0 ? nEvents : nentries));
    std::cout << "Number of entries to process: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        std::bitset<64> bits(gl1_scaledvec);
        std::cout << "Entry " << i << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;

        for (int j = 0; j < 64; ++j) {
            if (bits.test(j)) {
                float scaleFactor = 1.0 / (1.0 + static_cast<float>(gl1_live[j]) / static_cast<float>(gl1_scaled[j]));
                std::cout << "Trigger index " << j << " fired. Scale factor calculated: " << scaleFactor << std::endl;
                if (!cluster_ecore->empty()) {
                    float maxEcore = *std::max_element(cluster_ecore->begin(), cluster_ecore->end());
                    histograms_ecore[j]->Fill(maxEcore, scaleFactor);
                    std::cout << "Max ECore: " << maxEcore << " added to histogram." << std::endl;
                }
                if (!emcal_energy->empty()) {
                    float maxEnergy = *std::max_element(emcal_energy->begin(), emcal_energy->end());
                    histograms_energy[j]->Fill(maxEnergy, scaleFactor);
                    std::cout << "Max Energy: " << maxEnergy << " added to histogram." << std::endl;
                }
            }cat c
        }
    }

    TFile *outfile = new TFile(output_filename, "RECREATE");
    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i]->Write();
        histograms_energy[i]->Write();
    }
    std::cout << "Histograms written to file " << output_filename << std::endl;
    outfile->Close();
    file->Close();
    std::cout << "Files closed." << std::endl;
}
