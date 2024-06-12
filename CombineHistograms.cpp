#include <iostream>
#include <map>
#include <functional>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TSystem.h>
#include <TList.h>
#include <TString.h>

void CombineHistograms(const char* directoryPath) {
    // Set the output file name
    TString outputPath = TString(directoryPath) + "/combined/Combined_TriggerHists.root";

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
    std::map<std::string, TDirectory*> directoryMap;

    while ((filename = gSystem->GetDirEntry(dirp))) {
        TString filepath = TString(directoryPath) + "/" + TString(filename);
        if (!filepath.EndsWith(".root")) continue;

        TFile *file = TFile::Open(filepath);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            continue;
        }

        std::cout << "Processing file: " << filepath << std::endl;

        // Function to recursively process directories
        std::function<void(TDirectory*, const std::string&)> processDirectory = [&](TDirectory* sourceDir, const std::string& parentPath) {
            TList *list = sourceDir->GetListOfKeys();
            TIter iter(list);
            TKey *key;

            while ((key = (TKey*)iter())) {
                TObject *obj = key->ReadObj();
                if (obj->InheritsFrom("TH1")) {
                    TH1 *hist = (TH1*)obj;
                    std::string histPath = parentPath + "/" + hist->GetName();
                    std::cout << "Processing histogram: " << histPath << std::endl;

                    if (hist->GetEntries() > 0) {
                        if (histogramMap.find(histPath) == histogramMap.end()) {
                            histogramMap[histPath] = (TH1*)hist->Clone();
                            histogramMap[histPath]->SetDirectory(0); // Remove from current file directory
                        } else {
                            histogramMap[histPath]->Add(hist);
                        }
                    }
                } else if (obj->InheritsFrom("TDirectory")) {
                    TDirectory *subDir = (TDirectory*)obj;
                    std::string subDirPath = parentPath + "/" + subDir->GetName();
                    if (directoryMap.find(subDirPath) == directoryMap.end()) {
                        // Navigate to the parent directory in the output file
                        TDirectory* parentDir = outputFile;
                        if (!parentPath.empty()) {
                            parentDir = directoryMap[parentPath];
                            if (!parentDir) {
                                std::cerr << "Parent directory not found for: " << parentPath << std::endl;
                                continue;
                            }
                        }
                        parentDir->cd();
                        TDirectory *newDir = parentDir->mkdir(subDir->GetName());
                        directoryMap[subDirPath] = newDir;
                    }
                    processDirectory(subDir, subDirPath);
                }
            }
        };

        processDirectory(file, "");

        file->Close();
    }

    // Write combined histograms to the output file
    for (auto &entry : histogramMap) {
        size_t pos = entry.first.find_last_of('/');
        std::string dirPath = entry.first.substr(0, pos);
        std::string histName = entry.first.substr(pos + 1);

        if (directoryMap.find(dirPath) != directoryMap.end()) {
            directoryMap[dirPath]->cd();
        } else {
            outputFile->cd();
        }

        if (entry.second->GetEntries() > 0) { // Ensure it has entries before saving
            std::cout << "Writing histogram: " << histName << " with " << entry.second->GetEntries() << " entries." << std::endl;
            entry.second->Write(histName.c_str());
        }
    }

    outputFile->Close();
    std::cout << "Histograms combined and saved to " << outputPath << std::endl;
}
