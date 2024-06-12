#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <algorithm>

/*
make list file: ./createListFile /sphenix/u/patsfan753/scratch/analysis/calotriggeremulator/output/44686 44686
 */
void CaloTriggerSegments(const char* filename, const char* output_filename, Long64_t nEvents = 0) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    std::cout << "ROOT file opened successfully: " << filename << std::endl;
    
    // Get the TTree from the file
    std::cout << "Getting TTree from the file." << std::endl;
    TTree *tree = (TTree*)file->Get("ttree"); // Replace 'ttree' with the actual TTree name
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'ttree' in the file." << std::endl;
        return;
    }
    std::cout << "TTree loaded successfully." << std::endl;
    
    // Declare variables to hold branch data
    ULong64_t gl1_clock;
    ULong64_t gl1_scaled[64];
    ULong64_t gl1_live[64];
    ULong64_t gl1_raw[64];
    //    long gl1_rawvec;
    //    long gl1_livevec;
    //    long gl1_scaledvec;
    std::vector<float> *emcal_energy = 0;
    std::vector<float> *emcal_phibin = 0;
    std::vector<float> *emcal_etabin = 0;
    
    
    // Set branch addresses
    tree->SetBranchAddress("gl1_clock", &gl1_clock);
    tree->SetBranchAddress("gl1_scaled", gl1_scaled);
    tree->SetBranchAddress("gl1_live", gl1_live);
    tree->SetBranchAddress("gl1_raw", gl1_raw);
    //    tree->SetBranchAddress("gl1_rawvec", &gl1_rawvec);
    //    tree->SetBranchAddress("gl1_livevec", &gl1_livevec);
    //    tree->SetBranchAddress("gl1_scaledvec", &gl1_scaledvec);
    tree->SetBranchAddress("emcal_energy", &emcal_energy);
    tree->SetBranchAddress("emcal_phibin", &emcal_phibin);
    tree->SetBranchAddress("emcal_etabin", &emcal_etabin);
    std::cout << "Branch addresses set successfully." << std::endl;
    
    
    // Create histograms for each trigger and summed across all triggers
    std::cout << "Creating histograms for each trigger and summed across all triggers." << std::endl;
    TH1F *histograms_raw[64];
    TH1F *histograms_max8x8_raw[64];
    TH2F *histograms_2d_raw[64];
    TH1F *histogram_all_raw_emcal = new TH1F("all_raw_emcal", "All Triggers Raw EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_raw_max8x8 = new TH1F("all_raw_max8x8", "All Triggers Raw Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_raw_2d = new TH2F("all_raw_2d", "All Triggers Raw Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    TH1F *histograms_live[64];
    TH1F *histograms_max8x8_live[64];
    TH2F *histograms_2d_live[64];
    TH1F *histogram_all_live_emcal = new TH1F("all_live_emcal", "All Triggers Live EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_live_max8x8 = new TH1F("all_live_max8x8", "All Triggers Live Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_live_2d = new TH2F("all_live_2d", "All Triggers Live Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    TH1F *histograms_scaled[64];
    TH1F *histograms_max8x8_scaled[64];
    TH2F *histograms_2d_scaled[64];
    TH1F *histogram_all_scaled_emcal = new TH1F("all_scaled_emcal", "All Triggers Scaled EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_scaled_max8x8 = new TH1F("all_scaled_max8x8", "All Triggers Scaled Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_scaled_2d = new TH2F("all_scaled_2d", "All Triggers Scaled Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    for (int i = 0; i < 64; ++i) {
        histograms_raw[i] = new TH1F(Form("trigger_%d_raw_emcal", i), Form("Trigger %d Raw EMCal Energy", i), 100, 0, 20);
        histograms_max8x8_raw[i] = new TH1F(Form("trigger_%d_raw_max8x8", i), Form("Trigger %d Raw Max 8x8 EMCal Energy", i), 100, 0, 20);
        histograms_2d_raw[i] = new TH2F(Form("trigger_%d_raw_2d", i), Form("Trigger %d Raw Eta-Phi EMCal Energy", i), 96, 0, 96, 256, 0, 256);
        
        histograms_live[i] = new TH1F(Form("trigger_%d_live_emcal", i), Form("Trigger %d Live EMCal Energy", i), 100, 0, 20);
        histograms_max8x8_live[i] = new TH1F(Form("trigger_%d_live_max8x8", i), Form("Trigger %d Live Max 8x8 EMCal Energy", i), 100, 0, 20);
        histograms_2d_live[i] = new TH2F(Form("trigger_%d_live_2d", i), Form("Trigger %d Live Eta-Phi EMCal Energy", i), 96, 0, 96, 256, 0, 256);
        
        histograms_scaled[i] = new TH1F(Form("trigger_%d_scaled_emcal", i), Form("Trigger %d Scaled EMCal Energy", i), 100, 0, 20);
        histograms_max8x8_scaled[i] = new TH1F(Form("trigger_%d_scaled_max8x8", i), Form("Trigger %d Scaled Max 8x8 EMCal Energy", i), 100, 0, 20);
        histograms_2d_scaled[i] = new TH2F(Form("trigger_%d_scaled_2d", i), Form("Trigger %d Scaled Eta-Phi EMCal Energy", i), 96, 0, 96, 256, 0, 256);
    }
    std::cout << "Histograms created successfully." << std::endl;
    // Loop over all entries in the TTree
    // Determine number of entries to process
    Long64_t nentries = tree->GetEntries();
    if (nEvents > 0 && nEvents < nentries) {
        nentries = nEvents;
    }
    std::cout << "Number of entries to process: " << nentries << std::endl;
    

    for (Long64_t i = 0; i < nentries; i++) {
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            std::cout << "Processing entry " << i << std::endl;
            // Calculate the EMCal energy sum
            float emcal_energy_sum = 0;
            for (size_t j = 0; j < emcal_energy->size(); j++) {
                emcal_energy_sum += emcal_energy->at(j);
            }
            std::cout << "Entry " << i << ": EMCal energy sum calculated." << std::endl;
            // Calculate the maximum 8x8 EMCal energy sum
            float max8x8_energy_sum = 0;
            const int eta_bins = 96; // Number of eta bins (needs to be adjusted according to actual binning)
            const int phi_bins = 256; // Number of phi bins (needs to be adjusted according to actual binning)
            
            // Create a 2D array to hold the energy values for each tower
            float tower_energy[phi_bins][eta_bins] = {0};
            
            // Fill the 2D array with energy values
            for (size_t j = 0; j < emcal_energy->size(); j++) {
                int eta = static_cast<int>(emcal_etabin->at(j));
                int phi = static_cast<int>(emcal_phibin->at(j));
                tower_energy[phi][eta] = emcal_energy->at(j);
            }
            std::cout << "Entry " << i << ": 2D tower energy array filled." << std::endl;
            // Iterate over each tower to calculate the 8x8 energy sum
            for (int phi = 0; phi < phi_bins - 7; phi++) {
                for (int eta = 0; eta < eta_bins - 7; eta++) {
                    float sum8x8 = 0;
                    for (int dphi = 0; dphi < 8; dphi++) {
                        for (int deta = 0; deta < 8; deta++) {
                            sum8x8 += tower_energy[phi + dphi][eta + deta];
                        }
                    }
                    max8x8_energy_sum = std::max(max8x8_energy_sum, sum8x8);
                }
            }
            std::cout << "Entry " << i << ": Max 8x8 energy sum calculated." << std::endl;
            // Fill histograms based on raw triggers
            for (int j = 0; j < 64; ++j) {
                if (gl1_raw[j] >= 1) {
                    histograms_raw[j]->Fill(emcal_energy_sum);
                    histograms_max8x8_raw[j]->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_raw[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                    histogram_all_raw_emcal->Fill(emcal_energy_sum);
                    histogram_all_raw_max8x8->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_raw_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                }
            }
            std::cout << "Entry " << i << ": Raw trigger histograms filled." << std::endl;

            // Fill histograms based on live triggers
            for (int j = 0; j < 64; ++j) {
                if (gl1_live[j] >= 1) {
                    histograms_live[j]->Fill(emcal_energy_sum);
                    histograms_max8x8_live[j]->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_live[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                    histogram_all_live_emcal->Fill(emcal_energy_sum);
                    histogram_all_live_max8x8->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_live_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                }
            }
            std::cout << "Entry " << i << ": Live trigger histograms filled." << std::endl;            // Fill histograms based on scaled triggers
            for (int j = 0; j < 64; ++j) {
                if (gl1_scaled[j] >= 1) {
                    histograms_scaled[j]->Fill(emcal_energy_sum);
                    histograms_max8x8_scaled[j]->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_scaled[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                    histogram_all_scaled_emcal->Fill(emcal_energy_sum);
                    histogram_all_scaled_max8x8->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_scaled_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k));
                    }
                }
            }
            std::cout << "Entry " << i << ": Scaled trigger histograms filled." << std::endl;
        }
        
        // Save histograms to an output ROOT file
        TFile *outfile = new TFile(output_filename, "RECREATE");
        
        outfile->mkdir("raw");
        outfile->cd("raw");
        std::cout << "Writing raw histograms..." << std::endl;
        for (int i = 0; i < 64; ++i) {
            if (histograms_raw[i]) histograms_raw[i]->Write();
            if (histograms_max8x8_raw[i]) histograms_max8x8_raw[i]->Write();
            if (histograms_2d_raw[i]) histograms_2d_raw[i]->Write();
        }
        if (histogram_all_raw_emcal) histogram_all_raw_emcal->Write();
        if (histogram_all_raw_max8x8) histogram_all_raw_max8x8->Write();
        if (histogram_all_raw_2d) histogram_all_raw_2d->Write();
        std::cout << "Raw histograms written successfully." << std::endl;

        outfile->mkdir("live");
        outfile->cd("live");
        std::cout << "Writing live histograms..." << std::endl;
        for (int i = 0; i < 64; ++i) {
            if (histograms_live[i]) histograms_live[i]->Write();
            if (histograms_max8x8_live[i]) histograms_max8x8_live[i]->Write();
            if (histograms_2d_live[i]) histograms_2d_live[i]->Write();
        }
        if (histogram_all_live_emcal) histogram_all_live_emcal->Write();
        if (histogram_all_live_max8x8) histogram_all_live_max8x8->Write();
        if (histogram_all_live_2d) histogram_all_live_2d->Write();
        std::cout << "Live histograms written successfully." << std::endl;

        outfile->mkdir("scaled");
        outfile->cd("scaled");
        std::cout << "Writing scaled histograms..." << std::endl;
        for (int i = 0; i < 64; ++i) {
            if (histograms_scaled[i]) histograms_scaled[i]->Write();
            if (histograms_max8x8_scaled[i]) histograms_max8x8_scaled[i]->Write();
            if (histograms_2d_scaled[i]) histograms_2d_scaled[i]->Write();
        }
        if (histogram_all_scaled_emcal) histogram_all_scaled_emcal->Write();
        if (histogram_all_scaled_max8x8) histogram_all_scaled_max8x8->Write();
        if (histogram_all_scaled_2d) histogram_all_scaled_2d->Write();
        std::cout << "Scaled histograms written successfully." << std::endl;

        outfile->Close();
        std::cout << "Output ROOT file closed." << std::endl;;
        
        // Clean up
        try {
            delete file;
            std::cout << "Deleted file." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting file: " << e.what() << std::endl;
        }

        for (int i = 0; i < 64; ++i) {
            try {
                delete histograms_raw[i];
                std::cout << "Deleted histograms_raw[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_raw[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_max8x8_raw[i];
                std::cout << "Deleted histograms_max8x8_raw[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_max8x8_raw[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_2d_raw[i];
                std::cout << "Deleted histograms_2d_raw[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_2d_raw[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_live[i];
                std::cout << "Deleted histograms_live[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_live[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_max8x8_live[i];
                std::cout << "Deleted histograms_max8x8_live[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_max8x8_live[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_2d_live[i];
                std::cout << "Deleted histograms_2d_live[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_2d_live[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_scaled[i];
                std::cout << "Deleted histograms_scaled[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_scaled[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_max8x8_scaled[i];
                std::cout << "Deleted histograms_max8x8_scaled[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_max8x8_scaled[" << i << "]: " << e.what() << std::endl;
            }

            try {
                delete histograms_2d_scaled[i];
                std::cout << "Deleted histograms_2d_scaled[" << i << "]." << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error deleting histograms_2d_scaled[" << i << "]: " << e.what() << std::endl;
            }
        }

        try {
            delete histogram_all_raw_emcal;
            std::cout << "Deleted histogram_all_raw_emcal." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_raw_emcal: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_raw_max8x8;
            std::cout << "Deleted histogram_all_raw_max8x8." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_raw_max8x8: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_raw_2d;
            std::cout << "Deleted histogram_all_raw_2d." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_raw_2d: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_live_emcal;
            std::cout << "Deleted histogram_all_live_emcal." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_live_emcal: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_live_max8x8;
            std::cout << "Deleted histogram_all_live_max8x8." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_live_max8x8: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_live_2d;
            std::cout << "Deleted histogram_all_live_2d." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_live_2d: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_scaled_emcal;
            std::cout << "Deleted histogram_all_scaled_emcal." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_scaled_emcal: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_scaled_max8x8;
            std::cout << "Deleted histogram_all_scaled_max8x8." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_scaled_max8x8: " << e.what() << std::endl;
        }

        try {
            delete histogram_all_scaled_2d;
            std::cout << "Deleted histogram_all_scaled_2d." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error deleting histogram_all_scaled_2d: " << e.what() << std::endl;
        }
    }
}
