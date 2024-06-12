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
    long gl1_rawvec;
    long gl1_livevec;
    long gl1_scaledvec;
    float mbd_vertex_z = 0;
    
    std::vector<float> *emcal_energy = 0;
    std::vector<float> *cluster_ecore = 0;
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
    tree->SetBranchAddress("cluster_ecore", &cluster_ecore);
    tree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);  // Set branch address for z-vertex
    
    
    std::cout << "Branch addresses set successfully." << std::endl;
    
    
    // Create histograms for each trigger and summed across all triggers
    std::cout << "Creating histograms for each trigger and summed across all triggers." << std::endl;
    TH1F *histograms_raw[64]= {nullptr};
    TH1F *histograms_raw_ecore[64] = {nullptr};
    TH1F *histograms_max8x8_raw[64] = {nullptr};
    TH2F *histograms_2d_raw[64] = {nullptr};
    TH1F *histogram_all_raw_emcal = new TH1F("all_raw_emcal", "All Triggers Raw EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_raw_max8x8 = new TH1F("all_raw_max8x8", "All Triggers Raw Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_raw_2d = new TH2F("all_raw_2d", "All Triggers Raw Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    TH1F *histograms_live[64] = {nullptr};
    TH1F *histograms_live_ecore[64] = {nullptr};
    TH1F *histograms_max8x8_live[64] = {nullptr};
    TH2F *histograms_2d_live[64] = {nullptr};
    TH1F *histogram_all_live_emcal = new TH1F("all_live_emcal", "All Triggers Live EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_live_max8x8 = new TH1F("all_live_max8x8", "All Triggers Live Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_live_2d = new TH2F("all_live_2d", "All Triggers Live Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    TH1F *histograms_scaled[64] = {nullptr};
    TH1F *histograms_scaled_ecore[64] = {nullptr};
    TH1F *histograms_max8x8_scaled[64] = {nullptr};
    TH2F *histograms_2d_scaled[64] = {nullptr};
    TH1F *histogram_all_scaled_emcal = new TH1F("all_scaled_emcal", "All Triggers Scaled EMCal Energy", 100, 0, 20);
    TH1F *histogram_all_scaled_max8x8 = new TH1F("all_scaled_max8x8", "All Triggers Scaled Max 8x8 EMCal Energy", 100, 0, 20);
    TH2F *histogram_all_scaled_2d = new TH2F("all_scaled_2d", "All Triggers Scaled Eta-Phi EMCal Energy", 96, 0, 96, 256, 0, 256); // Adjust bins if necessary
    
    for (int i = 0; i < 64; ++i) {
        histograms_raw[i] = new TH1F(Form("trigger_%d_raw_emcal", i), Form("Trigger %d Raw EMCal Energy", i), 100, 0, 20);
        histograms_raw_ecore[i] = new TH1F(Form("trigger_%d_raw_cluster_ecore", i), Form("Trigger %d Raw Cluster Ecore", i), 100, 0, 20);
        histograms_max8x8_raw[i] = new TH1F(Form("trigger_%d_raw_max8x8", i), Form("Trigger %d Raw Max 8x8 EMCal Energy", i), 100, 0, 20);
        histograms_2d_raw[i] = new TH2F(Form("trigger_%d_raw_2d", i), Form("Trigger %d Raw Eta-Phi EMCal Energy", i), 96, 0, 96, 256, 0, 256);
        
        histograms_live[i] = new TH1F(Form("trigger_%d_live_emcal", i), Form("Trigger %d Live EMCal Energy", i), 100, 0, 20);
        histograms_live_ecore[i] = new TH1F(Form("trigger_%d_live_cluster_ecore", i), Form("Trigger %d Live Cluster Ecore", i), 100, 0, 20);
        histograms_max8x8_live[i] = new TH1F(Form("trigger_%d_live_max8x8", i), Form("Trigger %d Live Max 8x8 EMCal Energy", i), 100, 0, 20);
        histograms_2d_live[i] = new TH2F(Form("trigger_%d_live_2d", i), Form("Trigger %d Live Eta-Phi EMCal Energy", i), 96, 0, 96, 256, 0, 256);
        
        histograms_scaled[i] = new TH1F(Form("trigger_%d_scaled_emcal", i), Form("Trigger %d Scaled EMCal Energy", i), 100, 0, 20);
        histograms_scaled_ecore[i] = new TH1F(Form("trigger_%d_scaled_cluster_ecore", i), Form("Trigger %d Scaled Cluster Ecore", i), 100, 0, 20);
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
            std::cout << "Entry " << i << ": Total EMCal Energy calculated." << std::endl;
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
                // Iterate over eta bins, ensuring that we don't exceed array boundaries
                for (int eta = 0; eta < eta_bins - 7; eta++) {
                    float sum8x8 = 0;  // Initialize the sum for the current 8x8 block to zero
                    
                    // Iterate over the 8x8 block in the phi direction
                    for (int dphi = 0; dphi < 8; dphi++) {
                        // Iterate over the 8x8 block in the eta direction
                        for (int deta = 0; deta < 8; deta++) {
                            // Sum the energy in the current 8x8 block
                            sum8x8 += tower_energy[phi + dphi][eta + deta];
                        }
                    }
                    
                    // Update the maximum 8x8 energy sum if the current sum is larger
                    max8x8_energy_sum = std::max(max8x8_energy_sum, sum8x8);
                }
            }
            std::cout << "Entry " << i << ": Max 8x8 energy sum calculated." << std::endl;
            // Fill histograms based on raw triggers
            for (int j = 0; j < 64; ++j) {
                float weight_raw = gl1_raw[j]; // Using raw trigger counts as weights
                if (weight_raw >= 1) {
                    //if (std::abs(mbd_vertex_z) >= 10) continue;
                    histograms_raw[j]->Fill(emcal_energy_sum, weight_raw);
                    histograms_max8x8_raw[j]->Fill(max8x8_energy_sum, weight_raw);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_raw[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_raw);
                    }
                    histogram_all_raw_emcal->Fill(emcal_energy_sum, weight_raw);
                    histogram_all_raw_max8x8->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_raw_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_raw);
                    }
                    for (float ecore : *cluster_ecore) {
                        histograms_raw_ecore[j]->Fill(ecore, weight_raw);
                    }
                }
            }
            std::cout << "Entry " << i << ": Raw trigger histograms filled." << std::endl;

            // Fill histograms based on live triggers
            for (int j = 0; j < 64; ++j) {
                float weight_live = gl1_raw[j];
                if (weight_live >= 1) {
                    //if (std::abs(mbd_vertex_z) >= 10) continue;
                    histograms_live[j]->Fill(emcal_energy_sum, weight_live);
                    histograms_max8x8_live[j]->Fill(max8x8_energy_sum, weight_live);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_live[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_live);
                    }
                    histogram_all_live_emcal->Fill(emcal_energy_sum, weight_live);
                    histogram_all_live_max8x8->Fill(max8x8_energy_sum);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_live_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_live);
                    }
                    for (float ecore : *cluster_ecore) {
                        histograms_live_ecore[j]->Fill(ecore, weight_live);
                    }
                }
            }
            std::cout << "Entry " << i << ": Live trigger histograms filled." << std::endl;          // Fill histograms based on scaled triggers
            for (int j = 0; j < 64; ++j) {
                float weight_scaled = gl1_scaled[j];
                if (weight_scaled >= 1) {
                    //if (std::abs(mbd_vertex_z) >= 10) continue;
                    histograms_scaled[j]->Fill(emcal_energy_sum, weight_scaled);
                    histograms_max8x8_scaled[j]->Fill(max8x8_energy_sum, weight_scaled);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histograms_2d_scaled[j]->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_scaled);
                    }
                    histogram_all_scaled_emcal->Fill(emcal_energy_sum, weight_scaled);
                    histogram_all_scaled_max8x8->Fill(max8x8_energy_sum, weight_scaled);
                    for (size_t k = 0; k < emcal_energy->size(); k++) {
                        histogram_all_scaled_2d->Fill(emcal_etabin->at(k), emcal_phibin->at(k), emcal_energy->at(k) * weight_scaled);
                    }
                    for (float ecore : *cluster_ecore) {
                        histograms_scaled_ecore[j]->Fill(ecore, weight_scaled);
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
            if (histograms_raw_ecore[i]) histograms_raw_ecore[i]->Write();
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
            if (histograms_live_ecore[i]) histograms_live_ecore[i]->Write();
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
            if (histograms_scaled_ecore[i]) histograms_scaled_ecore[i]->Write();
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
            // Delete histograms_raw and safely handle exceptions
            if (histograms_raw[i] != nullptr) {
                delete histograms_raw[i];
                histograms_raw[i] = nullptr;
                std::cout << "Deleted histograms_raw[" << i << "]." << std::endl;
            }

            if (histograms_raw_ecore[i] != nullptr) {
                delete histograms_raw_ecore[i];
                histograms_raw_ecore[i] = nullptr;
                std::cout << "Deleted histograms_raw_ecore[" << i << "]." << std::endl;
            }

            if (histograms_max8x8_raw[i] != nullptr) {
                delete histograms_max8x8_raw[i];
                histograms_max8x8_raw[i] = nullptr;
                std::cout << "Deleted histograms_max8x8_raw[" << i << "]." << std::endl;
            }

            if (histograms_2d_raw[i] != nullptr) {
                delete histograms_2d_raw[i];
                histograms_2d_raw[i] = nullptr;
                std::cout << "Deleted histograms_2d_raw[" << i << "]." << std::endl;
            }

            if (histograms_live[i] != nullptr) {
                delete histograms_live[i];
                histograms_live[i] = nullptr;
                std::cout << "Deleted histograms_live[" << i << "]." << std::endl;
            }

            if (histograms_live_ecore[i] != nullptr) {
                delete histograms_live_ecore[i];
                histograms_live_ecore[i] = nullptr;
                std::cout << "Deleted histograms_live_ecore[" << i << "]." << std::endl;
            }

            if (histograms_max8x8_live[i] != nullptr) {
                delete histograms_max8x8_live[i];
                histograms_max8x8_live[i] = nullptr;
                std::cout << "Deleted histograms_max8x8_live[" << i << "]." << std::endl;
            }

            if (histograms_2d_live[i] != nullptr) {
                delete histograms_2d_live[i];
                histograms_2d_live[i] = nullptr;
                std::cout << "Deleted histograms_2d_live[" << i << "]." << std::endl;
            }

            if (histograms_scaled[i] != nullptr) {
                delete histograms_scaled[i];
                histograms_scaled[i] = nullptr;
                std::cout << "Deleted histograms_scaled[" << i << "]." << std::endl;
            }

            if (histograms_scaled_ecore[i] != nullptr) {
                delete histograms_scaled_ecore[i];
                histograms_scaled_ecore[i] = nullptr;
                std::cout << "Deleted histograms_scaled_ecore[" << i << "]." << std::endl;
            }

            if (histograms_max8x8_scaled[i] != nullptr) {
                delete histograms_max8x8_scaled[i];
                histograms_max8x8_scaled[i] = nullptr;
                std::cout << "Deleted histograms_max8x8_scaled[" << i << "]." << std::endl;
            }

            if (histograms_2d_scaled[i] != nullptr) {
                delete histograms_2d_scaled[i];
                histograms_2d_scaled[i] = nullptr;
                std::cout << "Deleted histograms_2d_scaled[" << i << "]." << std::endl;
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

//        try {
//            delete histogram_all_raw_2d;
//            std::cout << "Deleted histogram_all_raw_2d." << std::endl;
//        } catch (const std::exception& e) {
//            std::cerr << "Error deleting histogram_all_raw_2d: " << e.what() << std::endl;
//        }

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

