#include "plot_data.h" // Include the header file with function declarations
#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include "formatLabelName.h"
#include "TGraphErrors.h"


extern TTreeReader dataReader;
extern TTreeReader mcReader;
extern std::map<std::string, HistConfig> histConfigs;

#include <optional>

template<typename T>
void FillHistogram(TTreeReader& reader, const std::string& branchName, TH1D* hist,
        BaseKinematicCuts& kinematicCuts, int fitIndex, bool isMC) {
    TTreeReaderValue<T> val(reader, branchName.c_str());
    while (reader.Next()) {
        if (kinematicCuts.applyCuts(fitIndex, isMC)) {
            hist->Fill(*val);
        }
    }
}

template<typename T>
void FillHistogramForBins(TTreeReader& reader, const std::string& branchName, TH1D* hist,
      BaseKinematicCuts& kinematicCuts, int fitIndex, bool isMC,
      double varMin, double varMax) {
    TTreeReaderValue<T> val(reader, branchName.c_str());
    TTreeReaderValue<double> currentVariable(reader, propertyNames[fitIndex].c_str());

    while (reader.Next()) {
        if (kinematicCuts.applyCuts(fitIndex, isMC) && *currentVariable >= varMin && *currentVariable < varMax) {
            hist->Fill(*val);
        }
    }
}



template<typename T1, typename T2>
void createAndFillHistogram(TTreeReader& reader, TH2D* hist, const std::string& branchX, 
                            const std::string& branchY, KinematicCuts& kinematicCuts) {
    TTreeReaderValue<T1> valX(reader, branchX.c_str());
    TTreeReaderValue<T2> valY(reader, branchY.c_str());

    reader.Restart();
    while (reader.Next()) {
        if (kinematicCuts.applyCuts(0, false)) {
            hist->Fill(*valX, *valY);
        }
    }
}

void createIntegratedKinematicPlots() {
    const std::string outputDir = "output/integrated_plots/";
    const std::vector<std::string> branchesToSkip = {"helicity", "beam_pol", 
        "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"};

    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.05); // Increase the text size globally
    bool restart = true;
    for (Int_t i = 0; i < branches->GetEntries(); ++i) {

        TBranch* branch = (TBranch*)branches->At(i);
        std::string branchName = branch->GetName();
        if (branchName == "e_p" && restart) {
          // stupid hack to get it to do the runnum plot instead of it being blank 
          // due to reader restarts
          i = 0; 
          restart = false;
        }
        branch = (TBranch*)branches->At(i);
        branchName = branch->GetName();

        if (std::find(branchesToSkip.begin(), branchesToSkip.end(), branchName) != 
          branchesToSkip.end()) {
            continue; // Skip this branch
        }

        HistConfig config = {100, 0, 1}; // Default configuration
        if (histConfigs.find(branchName) != histConfigs.end()) {
            config = histConfigs[branchName];
        }

        TH1D* dataHist = new TH1D((branchName + "_data").c_str(), "", config.nBins, config.xMin, config.xMax);
        TH1D* mcHist = new TH1D((branchName + "_mc").c_str(), "", config.nBins, config.xMin, config.xMax);

        // Set x-axis title
        dataHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
        mcHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());

        dataHist->GetXaxis()->SetTitleSize(0.05); // Increase the x-axis title font size
        dataHist->GetYaxis()->SetTitleSize(0.05); // Increase the y-axis title font size
        mcHist->GetXaxis()->SetTitleSize(0.05);
        mcHist->GetYaxis()->SetTitleSize(0.05);

        // Set y-axis title and center it
        dataHist->GetYaxis()->SetTitle("Normalized Counts");
        dataHist->GetYaxis()->CenterTitle();
        mcHist->GetYaxis()->SetTitle("Normalized Counts");
        mcHist->GetYaxis()->CenterTitle();

        // Set y-axis title offset to make room for centering
        dataHist->GetYaxis()->SetTitleOffset(1.6);
        mcHist->GetYaxis()->SetTitleOffset(1.6);

        if (branchName == "runnum") {
          // Declare TTreeReaderValue for integers for dataReader
          TTreeReaderValue<int> dataVal(dataReader, branchName.c_str());
          // Fill histogram for dataReader
          FillHistogram<int>(dataReader, branchName, dataHist, *kinematicCuts, 0, false);
          // Check if the "runnum" branch exists in mcReader
          if (mcReader.GetTree()->GetBranch(branchName.c_str())) {
              // "runnum" branch exists, declare TTreeReaderValue for mcReader
              TTreeReaderValue<int> mcVal(mcReader, branchName.c_str());
              // Fill histogram for mcReader
              FillHistogram<int>(mcReader, branchName, mcHist, *mckinematicCuts, 0, true);
          } else {
              // "runnum" branch does not exist, use default value
              int defaultRunNum = 11;
              mcHist->Fill(defaultRunNum);
          }
        } else {
          // Declare TTreeReaderValue for doubles
          TTreeReaderValue<double> dataVal(dataReader, branchName.c_str());
          TTreeReaderValue<double> mcVal(mcReader, branchName.c_str());
          // Fill histograms for double values
          FillHistogram<double>(dataReader, branchName, dataHist, *kinematicCuts, 0, false);
          FillHistogram<double>(mcReader, branchName, mcHist, *mckinematicCuts, 0, true);
        }
        
        // Normalize the histograms
        dataHist->Scale(1.0 / dataHist->Integral());
        mcHist->Scale(1.0 / mcHist->Integral());

        // Find the maximum value for y-axis
        double maxY = 1.4*std::max(dataHist->GetMaximum(), mcHist->GetMaximum());
        dataHist->SetMaximum(1.1 * maxY);
        mcHist->SetMaximum(1.1 * maxY);

        // Create a canvas for drawing the histograms
        TCanvas* c = new TCanvas((branchName + "_canvas").c_str(), branchName.c_str(), 800, 600);
        // Adjust the margins to avoid cutting off labels
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.15);

        // Create a legend and adjust its font size
        TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9);
        leg->SetTextSize(0.04); // Increase the legend text size

        // Add entries to the legend with scientific notation for the number of entries
        dataHist->SetEntries(dataHist->GetEntries());
        mcHist->SetEntries(mcHist->GetEntries());
        leg->AddEntry(dataHist, (std::string("Data (") + std::to_string((int)dataHist->GetEntries()) + " entries)").c_str(), "l");
        leg->AddEntry(mcHist, (std::string("MC (") + std::to_string((int)mcHist->GetEntries()) + " entries)").c_str(), "l");

        // Set line colors for histograms
        dataHist->SetLineColor(kBlack);
        mcHist->SetLineColor(kRed);
        
        // Draw histograms on the canvas
        dataHist->Draw("HIST");
        mcHist->Draw("HISTSAME");
        leg->Draw();

        // Save the canvas to a file
        c->SaveAs((outputDir + branchName + ".png").c_str());
        // Clean up the created objects to avoid memory leaks
        delete dataHist;
        delete mcHist;
        delete c;
        delete leg;

        // Restart the TTreeReaders for the next branch
        dataReader.Restart();
        mcReader.Restart();
    }
}

void createIntegratedKinematicPlotsForBinsAndFits() {
    const std::string outputDir = "output/binned_plots/";
    const std::vector<std::string> branchesToSkip = {
        "helicity", "beam_pol", "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"
    };

    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.05); // Increase the text size globally

    // Loop over all sets of supplied kinematic variables
    for (size_t fitIndex = 0; fitIndex < allBins.size(); ++fitIndex) {
        std::string currentVariable = binNames[fitIndex]; // Assuming binNames is a vector<string> with the same size as allBins
        std::string branchVariable = propertyNames[fitIndex]; // Corresponding branch name for the current variable

        // Loop over all possible bins within the current set
        for (size_t binIndex = 0; binIndex < allBins[fitIndex].size() - 1; ++binIndex) {
            double binLowerEdge = allBins[fitIndex][binIndex];
            double binUpperEdge = allBins[fitIndex][binIndex + 1];

            // Format the bin edges with three decimal places
            std::ostringstream lowerEdgeStream, upperEdgeStream;
            lowerEdgeStream << std::fixed << std::setprecision(3) << binLowerEdge;
            upperEdgeStream << std::fixed << std::setprecision(3) << binUpperEdge;

            std::string binIndexLabel = "bin_" + std::to_string(binIndex + 1);

            // Now we iterate over all branches, except those we wish to skip
            for (Int_t i = 25; i < branches->GetEntries(); ++i) {
                TBranch* branch = (TBranch*)branches->At(i);
                std::string branchName = branch->GetName();

                if (std::find(branchesToSkip.begin(), branchesToSkip.end(), branchName) != branchesToSkip.end()) {
                    continue; // Skip this branch
                }

                // Determine histogram configuration, default or specific
                HistConfig config = {100, 0, 1}; // Default configuration
                if (histConfigs.find(branchName) != histConfigs.end()) {
                    config = histConfigs[branchName];
                } else {
                    std::cerr << "Warning: No specific histogram configuration found for " << branchName << ". Using default configuration." << std::endl;
                }

                // Create histogram title with formatted bin edges
                std::string formattedVariableName = formatLabelName(branchVariable);
                std::string plotTitle = lowerEdgeStream.str() + " < " + formattedVariableName + " < " + upperEdgeStream.str();

                // Create histograms with titles reflecting the bin edges and plotted variable
                std::string histName = currentVariable + "_" + branchName + "_" + binIndexLabel;
                TH1D* dataHist = new TH1D((histName + "_data").c_str(), plotTitle.c_str(), config.nBins, config.xMin, config.xMax);
                TH1D* mcHist = new TH1D((histName + "_mc").c_str(), plotTitle.c_str(), config.nBins, config.xMin, config.xMax);

                // Set histogram titles and styles
                dataHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                mcHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());

                dataHist->GetYaxis()->SetTitle("Normalized Counts");
                mcHist->GetYaxis()->SetTitle("Normalized Counts");

                dataHist->GetYaxis()->CenterTitle();
                mcHist->GetYaxis()->CenterTitle();

                dataHist->GetYaxis()->SetTitleOffset(1.6);
                mcHist->GetYaxis()->SetTitleOffset(1.6);

                dataHist->SetLineColor(kBlack);
                mcHist->SetLineColor(kRed);

                if (branchName == "runnum") {
                  TTreeReaderValue<int> dataVal(dataReader, branchName.c_str());
                  TTreeReaderValue<double> binVariable(dataReader, branchVariable.c_str());

                  if (mcReader.GetTree()->GetBranch(branchName.c_str())) {
                      TTreeReaderValue<int> mcVal(mcReader, branchName.c_str());
                      TTreeReaderValue<double> mcBinVariable(mcReader, branchVariable.c_str());
                      FillHistogramForBins<int>(dataReader, branchName, dataHist, *kinematicCuts, fitIndex, false, binLowerEdge, binUpperEdge);
                      FillHistogramForBins<int>(mcReader, branchName, mcHist, *mckinematicCuts, fitIndex, true, binLowerEdge, binUpperEdge);
                  } else {
                      int defaultRunNum = 11;
                      FillHistogramForBins<int>(dataReader, branchName, dataHist, *kinematicCuts, fitIndex, false, binLowerEdge, binUpperEdge);
                      mcHist->Fill(defaultRunNum);
                  }
                } else {
                  TTreeReaderValue<double> dataVal(dataReader, branchName.c_str());
                  TTreeReaderValue<double> binVariable(dataReader, branchVariable.c_str());
                  TTreeReaderValue<double> mcVal(mcReader, branchName.c_str());
                  TTreeReaderValue<double> mcBinVariable(mcReader, branchVariable.c_str());
                  FillHistogramForBins<double>(dataReader, branchName, dataHist, *kinematicCuts, fitIndex, false, binLowerEdge, binUpperEdge);
                  FillHistogramForBins<double>(mcReader, branchName, mcHist, *mckinematicCuts, fitIndex, true, binLowerEdge, binUpperEdge);
                }

                double dataScale = 0;
                double mcScale = 0;
                // Normalize the histograms
                if (dataHist->Integral() != 0) {
                    dataScale = dataHist->Integral();
                    dataHist->Scale(1.0 / dataScale);
                }
                if (mcHist->Integral() != 0) {
                    mcScale = mcHist->Integral();
                    mcHist->Scale(1.0 / mcScale);
                }

                // Find the maximum y-value between both histograms to set the y-axis range
                double maxY = std::max(dataHist->GetMaximum(), mcHist->GetMaximum());
                dataHist->SetMaximum(1.5 * maxY);
                mcHist->SetMaximum(1.5 * maxY);

                // Create a canvas for drawing the histograms
                TCanvas* c = new TCanvas((histName + "_canvas").c_str(), branchName.c_str(), 800, 600);
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                // Draw histograms on the canvas
                dataHist->Draw("HIST");
                mcHist->Draw("HIST SAME");

                // Create a legend for the histograms
                TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9);
                leg->AddEntry(dataHist, ("Data (" + std::to_string(static_cast<int>(dataHist->GetEntries())) + " entries)").c_str(), "l");
                leg->AddEntry(mcHist, ("MC (" + std::to_string(static_cast<int>(mcHist->GetEntries())) + " entries)").c_str(), "l");
                leg->Draw();

                // Save the canvas to a file
                std::string outputFileName = outputDir + histName + ".png";
                c->SaveAs(outputFileName.c_str());

                // // Now, create and handle the ratio plot for phi
                // if (branchName == "phi") { // make special rec/gen phi distribution
                //     // Assuming original histograms have a compatible binning scheme.
                //     const int targetBins = 24;
                //     double x[targetBins], y[targetBins], ex[targetBins], ey[targetBins];
                //     int originalBins = dataHist->GetNbinsX();
                //     int combineFactor = originalBins / targetBins; // Determine how many bins to combine.

                //     // Create new histograms for re-binned data
                //     TH1D* dataRebinned = new TH1D("data_rebinned", ";#phi;Normalized Counts", targetBins, 0, 2*3.1459);
                //     TH1D* mcRebinned = new TH1D("mc_rebinned", ";#phi;Normalized Counts", targetBins, 0, 2*3.1459);

                //     // Manually combine bins
                //     for (int i = 1; i <= targetBins; i++) {
                //         double dataSum = 0, mcSum = 0;
                //         for (int j = 0; j < combineFactor; j++) {
                //             int binIndex = (i - 1) * combineFactor + j + 1; // Calculate original bin index
                //             dataSum += dataHist->GetBinContent(binIndex);
                //             mcSum += mcHist->GetBinContent(binIndex);
                //         }
                //         dataRebinned->SetBinContent(i, dataSum);
                //         mcRebinned->SetBinContent(i, mcSum);
                //     }

                //     // filled dataRebinned and mcRebinned as before
                //     for (int i = 1; i <= targetBins; i++) {
                //         // Calculate the center of each bin as the x-value
                //         x[i-1] = dataRebinned->GetBinCenter(i);
                //         // y-values are the ratio
                //         double dataValue = dataRebinned->GetBinContent(i);
                //         double mcValue = mcRebinned->GetBinContent(i);
                //         y[i-1] = dataValue / mcValue;
                //         // Assuming no error in the x-direction (bin center)
                //         ex[i-1] = 0;

                //         // // Corrected calculation of error in y
                //         // if (dataValue > 0 && mcValue > 0) {
                //         //     double dataError = sqrt(dataValue); // For normalized, this should be scaled appropriately if needed
                //         //     double mcError = sqrt(mcValue); // Same as above
                //         //     // Calculate the relative errors and then the error on the ratio
                //         //     double relativeErrorData = dataError / dataValue;
                //         //     double relativeErrorMC = mcError / mcValue;
                //         //     // ey[i-1] = y[i-1] * sqrt(relativeErrorData * relativeErrorData + relativeErrorMC * relativeErrorMC);
                //         //     ey[i-1] = sqrt(mcError*mcError+dataError*dataError);
                //         // } else {
                //         //     ey[i-1] = 0; // Handle division by zero or negative values if necessary
                //         // }
                //     }


                    // Create a TGraphErrors with the arrays
                    TGraphErrors* graph = new TGraphErrors(targetBins, x, y, ex, ey);
                    graph->SetTitle(";#phi;Rec/Gen"); // Set title and axis labels
                    graph->SetMarkerStyle(21); // Choose a marker style
                    graph->SetLineColor(kBlue);

                    // Draw the TGraphErrors
                    TCanvas* graphCanvas = new TCanvas((histName + "_ratio_graph").c_str(), "Ratio Plot", 800, 600);
                    graph->Draw("APE"); // "AP" for drawing markers and lines, "E" for error bars

                    // Adjust Y-axis range manually
                    graph->GetHistogram()->SetMinimum(0.0); // Set minimum y-value
                    graph->GetHistogram()->SetMaximum(2.0); // Set maximum y-value

                    // Save the graph canvas
                    std::string graphOutputFileName = outputDir + histName + "_phi_ratio.png";
                    graphCanvas->SaveAs(graphOutputFileName.c_str());

                    // Cleanup
                    delete graph;
                    delete graphCanvas;
                }

                // Clean up
                delete dataHist;
                delete mcHist;
                delete c;
                delete leg;

                // Restart the TTreeReaders for the next branch
                dataReader.Restart();
                mcReader.Restart();
            }
        }
        // Increment the currentFits to process the next set of kinematic variables
        currentFits++;
    }
}

void createCorrelationPlotsforrunnum() {
    const std::string outputDir = "output/correlation_plots/";
    const std::vector<std::string> branchesToSkip = {"helicity", "beam_pol", "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"};

    // Assuming histConfigs is a global variable or it is accessible within this function's scope
    extern std::map<std::string, HistConfig> histConfigs;

    // Retrieve the list of branches
    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    // Prepare a vector with the names of the branches to be used
    std::vector<std::string> branchNames;
    for (Int_t i = 0; i < branches->GetEntries(); ++i) {
        std::string name = branches->At(i)->GetName();
        if (std::find(branchesToSkip.begin(), branchesToSkip.end(), name) == branchesToSkip.end()) {
            branchNames.push_back(name);
        }
    }

    // Initialize min and max runnum variables
    int minRunNum = std::numeric_limits<int>::max();
    int maxRunNum = std::numeric_limits<int>::min();

    // Loop over the dataReader to find the min and max runnum values
    TTreeReaderValue<int> runnumValue(dataReader, "runnum");
    while (dataReader.Next()) {
        minRunNum = std::min(minRunNum, *runnumValue);
        maxRunNum = std::max(maxRunNum, *runnumValue);
    }
    dataReader.Restart(); // Restart the dataReader after finding min and max

    // Adjust min and max for the plot range
    minRunNum -= 1;
    maxRunNum += 1;

    // Generate all possible pairs of branches to plot
    for (size_t i = 0; i < branchNames.size(); ++i) {
        const std::string& branchX = "runnum";
        const std::string& branchY = branchNames[i];

        TTreeReaderValue<int> valX(dataReader, "runnum");
        TTreeReaderValue<Double_t> valY(dataReader, branchY.c_str());

        // Get the configurations for each branch
        HistConfig configX = 
          histConfigs.count(branchX) ? histConfigs[branchX] : HistConfig{100, 0, 1};
        HistConfig configY = 
          histConfigs.count(branchY) ? histConfigs[branchY] : HistConfig{100, 0, 1};

        // Define the histogram for this pair
        std::string histName = branchX + "_vs_" + branchY;
        TH2D* hist = new TH2D(histName.c_str(), "",
          configX.nBins, minRunNum, maxRunNum, // Adjusted X-axis bins and range
          configY.nBins, configY.xMin, configY.xMax); // Y-axis bins and range

        // Set axis titles and center them
        hist->GetXaxis()->SetTitle(formatLabelName(branchX).c_str());
        hist->GetYaxis()->SetTitle(formatLabelName(branchY).c_str());
        hist->GetXaxis()->CenterTitle();
        hist->GetYaxis()->CenterTitle();

        // Increase the title font size for both axes
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetTitleSize(0.05);

        // Set the margins to avoid cutting off labels
        hist->GetXaxis()->SetTitleOffset(1.3);
        hist->GetYaxis()->SetTitleOffset(1.6);

        // Loop over dataReader to fill the histogram
        while (dataReader.Next()) {
            bool passedKinematicCuts = kinematicCuts->applyCuts(0, false);
            if (passedKinematicCuts) {
                hist->Fill(*valX, *valY);
            }
        }

        // Create a canvas for drawing the histogram
        TCanvas* c = new TCanvas(histName.c_str(), histName.c_str(), 800, 600);
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.15);
        c->SetRightMargin(0.15);

        // Draw the histogram on the canvas
        hist->Draw("COLZ"); // Draw as a 2D color plot

        // Save the canvas to a file
        c->SaveAs((outputDir + histName + ".png").c_str());

        // Clean up the created objects to avoid memory leaks
        delete hist;
        delete c;

        // Restart the TTreeReader for the next iteration
        dataReader.Restart();
    }
}

void createCorrelationPlots() {
    const std::string outputDir = "output/correlation_plots/";
    const std::vector<std::string> branchesToSkip = {"helicity", "beam_pol", "target_pol", "runnum", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"};

    // Assuming histConfigs is a global variable or it is accessible within this function's scope
    extern std::map<std::string, HistConfig> histConfigs;

    // Retrieve the list of branches
    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    // Prepare a vector with the names of the branches to be used
    std::vector<std::string> branchNames;
    for (Int_t i = 0; i < branches->GetEntries(); ++i) {
        std::string name = branches->At(i)->GetName();
        if (std::find(branchesToSkip.begin(), branchesToSkip.end(), name) == branchesToSkip.end()) {
            branchNames.push_back(name);
        }
    }

    // Generate all possible pairs of branches to plot
    for (size_t i = 0; i < branchNames.size(); ++i) {
        for (size_t j = i + 1; j < branchNames.size(); ++j) {
            const std::string& branchX = branchNames[i];
            const std::string& branchY = branchNames[j];

            TTreeReaderValue<int> runnum(dataReader, "runnum");
            TTreeReaderValue<Double_t> valX(dataReader, branchX.c_str());
            TTreeReaderValue<Double_t> valY(dataReader, branchY.c_str());

            // Get the configurations for each branch
            HistConfig configX = 
              histConfigs.count(branchX) ? histConfigs[branchX] : HistConfig{100, 0, 1};
            HistConfig configY = 
              histConfigs.count(branchY) ? histConfigs[branchY] : HistConfig{100, 0, 1};

            // Define the histogram for this pair
            std::string histName = branchX + "_vs_" + branchY;
            TH2D* hist = new TH2D(histName.c_str(), "",
              configX.nBins, configX.xMin, configX.xMax, // X-axis bins and range
              configY.nBins, configY.xMin, configY.xMax); // Y-axis bins and range

            // Set axis titles and center them
            hist->GetXaxis()->SetTitle(formatLabelName(branchX).c_str());
            hist->GetYaxis()->SetTitle(formatLabelName(branchY).c_str());
            hist->GetXaxis()->CenterTitle();
            hist->GetYaxis()->CenterTitle();

            // Increase the title font size for both axes
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetTitleSize(0.05);

            // Set the margins to avoid cutting off labels
            hist->GetXaxis()->SetTitleOffset(1.3);
            hist->GetYaxis()->SetTitleOffset(1.6);

            // Loop over dataReader to fill the histogram
            while (dataReader.Next()) {
                bool passedKinematicCuts = kinematicCuts->applyCuts(0, false);
                if (passedKinematicCuts) {
                    hist->Fill(*valX, *valY);
                }
            }

            // Create a canvas for drawing the histogram
            TCanvas* c = new TCanvas(histName.c_str(), histName.c_str(), 800, 600);
            c->SetLeftMargin(0.15);
            c->SetBottomMargin(0.15);
            c->SetRightMargin(0.15);

            // Draw the histogram on the canvas
            hist->Draw("COLZ"); // Draw as a 2D color plot

            // Save the canvas to a file
            c->SaveAs((outputDir + histName + ".png").c_str());

            // Clean up the created objects to avoid memory leaks
            delete hist;
            delete c;

            // Restart the TTreeReader for the next iteration
            dataReader.Restart();
        }
    }
}