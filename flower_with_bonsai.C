// C++
#include <algorithm> // std::sort
#include <iostream>
#include <math.h> // pow
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
// ROOT
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

// WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
// BONSAI
#include "WCSimBonsai.hh"
// FLOWER
#include "WCSimFLOWER.h"

// low energy reconstruction
int flower_with_bonsai(const char *detector, const char *filename = "../wcsim.root", const char *outfiletag = "a",
                       const int verbose = 0, const bool overwrite_nearest = false,
                       const double override_dark_rate = -99, const int start_event = 0, const int n_events = -1) {
    // sets the default nearest neighbour distances etc. Check DetectorEnumFromString inside
    // WCSimFLOWER.cpp for allowed values
    //  output file for input file.root will be in the format file_flower_bonsai<outfiletag>.root
    // if true, will overwrite the cached nearest neighbours file
    // if positive, will override the default dark rate (taken from WCSimRootOptions)

    // set up histogram
    TH1F *recEnergy = new TH1F("hErec", "Reconstructed Energy", 50, 0, 100);

    WCSimBonsai *bonsai = new WCSimBonsai();

    // Open the input ROOT file
    TFile *file = new TFile(filename, "read");
    if (!file->IsOpen()) {
        std::cout << "ERROR: Could not open input file " << filename << std::endl;
        return -1;
    }

    // Open the output ROOT file & create a TTree
    // Note: This is a flat tree that will be filled every trigger,
    // such that if you have an input file with multiple triggers,
    // some information (e.g. true_* variables) will be repeated
    // Note: If you have multiple true particles, the output tree won't account for this.
    // Currently it just saves the first true primary particle in the WCSimRootTrack TClonesArray
    TString outfilename(filename);
    outfilename.ReplaceAll(".root", TString::Format("_bonsai_flower%s%d.root", outfiletag, start_event));
    std::cout << "Output file: " << outfilename << std::endl;
    TFile *outfile = new TFile(outfilename, "RECREATE");
    TTree *out_tree = new TTree(
        "lowEreco", "FLOWER & BONSAI reconstruction (setup for events with 1 primary track & 1 trigger)");
    // setup true variables
    TVector3 true_pos, true_dir;
    double true_time, true_charge;
    float true_energy;
    out_tree->Branch("true_position", &true_pos);
    out_tree->Branch("true_direction", &true_dir);
    out_tree->Branch("true_time", &true_time, "true_time/D");
    out_tree->Branch("true_energy", &true_energy, "true_energy/F");
    out_tree->Branch("true_charge", &true_charge, "true_charge/D");
    // setup reco variables
    TVector3 reco_pos, reco_dir;
    double reco_time;
    float reco_energy, reco_neff, reco_neff2;
    out_tree->Branch("reco_position", &reco_pos);
    out_tree->Branch("reco_direction", &reco_dir);
    out_tree->Branch("reco_time", &reco_time, "reco_time/D");
    out_tree->Branch("reco_energy", &reco_energy, "reco_energy/F");
    out_tree->Branch("reco_neff", &reco_neff, "reco_neff/F");
    out_tree->Branch("reco_neff2", &reco_neff2, "reco_neff2/F");
    // setup hit variables
    int nhits, ndigits;
    out_tree->Branch("nhits", &nhits, "nhits/I");
    out_tree->Branch("ndigits", &ndigits, "ndigits/I");

    // Read geometry from WCSim file
    TTree *geotree = (TTree *)file->Get("wcsimGeoT");
    WCSimRootGeom *geo = 0;
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    if (geotree->GetEntries() == 0)
        exit(9); // exit if no geometry is defined in the WCSim file
    geotree->GetEntry(0);

    // Initalise BONSAI
    bonsai->Init(geo);

    // Read options tree from WCSim file
    TTree *opttree;
    file->GetObject("wcsimRootOptionsT", opttree);
    WCSimRootOptions *opt = 0;
    opttree->SetBranchAddress("wcsimrootoptions", &opt);
    if (opttree->GetEntries() == 0)
        exit(9); // exit if no options are defined in the WCSim file
    opttree->GetEntry(0);

    // Initialise FLOWER
    bool get_npmts_from_wcsimrootgeom = true;
    WCSimFLOWER *flower = new WCSimFLOWER(detector, geo, get_npmts_from_wcsimrootgeom, overwrite_nearest, verbose);
    if (override_dark_rate >= 0)
        flower->SetDarkRate(override_dark_rate);
    else
        flower->SetDarkRate(opt->GetPMTDarkRate("tank"));
    if (!flower->CheckNearestNeighbours()) {
        std::cerr << "Nearest neighbour distribution failure. Set nearest neighbour and try again" << std::endl;
        exit(8);
    }

    // Read event tree from WCSim file
    TTree *tree = (TTree *)file->Get("wcsimT");              // Get a pointer to the tree from the file
    WCSimRootEvent *event = new WCSimRootEvent();            // Create WCSimRootEvent to put stuff from the tree in
    tree->SetBranchAddress("wcsimrootevent", &event);        // Set branch address for reading from tree
    tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE); // Force deletion to prevent memory leak
    WCSimRootTrigger *trigger; // will contain triggers of event later (0: initial particle; 1..n: decay products)

    // Get the end event to reconstruct. If n_events is -1, simulate all events.
    int nWCSimEvents = tree->GetEntries();
    std::cout << "n_events is " << n_events << std::endl;
    if (n_events != -1) {
        nWCSimEvents = n_events + start_event;
    }

    std::cout << "nWCSimEvents is " << nWCSimEvents << std::endl;

    int ncherenkovdigihits, i;
    double x, y, z, eRec;

    // used by bonsai->BonsaiFit()
    float bsVertex[4], bsResult[6], bsGood[3];
    int *bsNhit;
    int bsNsel[2];

    // Loop over events
    std::cout << "Starting timer for the event loop" << std::endl;
    TStopwatch timer;

    for (int ev = start_event; ev < nWCSimEvents; ev++) {
        if (ev % 100 == 0 || ev == 0) {
            std::cout << "Event " << ev + 1 << " of " << nWCSimEvents << std::endl;
        }

        // reset tree variables
        true_pos.SetXYZ(-9999, -9999, -9999);
        true_dir.SetXYZ(-99, -99, -99);
        true_time = -9999;
        true_charge = -9999;
        true_energy = 0;
        reco_pos.SetXYZ(-9999, -9999, -9999);
        reco_dir.SetXYZ(-99, -99, -99);
        reco_time = -9999;
        reco_energy = -99;
        reco_neff = -99;
        reco_neff2 = -99;
        nhits = -99;
        ndigits = -99;

        if (verbose) {
            // red
            std::cout << "\033[1;31m";
            std::cout << "event number: " << ev << std::endl;
            // white
            std::cout << "\033[0m";
        }

        // Read the event from the tree into the WCSimRootEvent instance
        tree->GetEntry(ev);
        trigger = event->GetTrigger(0);

        // get the hits & digits
        nhits = trigger->GetNcherenkovhits();

        // See chapter 5 of doc/DetectorDocumentation.pdf in the WCSim repository
        // for more information on the structure of the root file.

        // get true vertex information
        true_pos.SetXYZ(trigger->GetVtx(0), trigger->GetVtx(1), trigger->GetVtx(2));
        true_time = trigger->GetHeader()->GetDate();

        // get true track information
        const int ntracks = trigger->GetNtrack();
        bool found_true_track = false;
        for (int itrack = 0; itrack < ntracks; itrack++) {
            TObject *element = (trigger->GetTracks())->At(itrack);
            if (!element)
                continue;
            WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack *>(element);

            if (wcsimroottrack->GetParenttype() == 0) {
                if (verbose > 1) {
                    // In green
                    std::cout << "\033[1;32m";

                    std::cout << "Track: " << itrack << std::endl << "  ";
                    // white
                    std::cout << "\033[0m";
                    int trackflag = wcsimroottrack->GetFlag();
                    if (trackflag == -1)
                        std::cout << "Primary neutrino track" << std::endl;
                    else if (trackflag == -2)
                        std::cout << "Neutrino target nucleus track" << std::endl;
                    else
                        std::cout << "Final state particle track" << std::endl;
                    printf("  Track ipnu (PDG code): %d\n", wcsimroottrack->GetIpnu());
                    printf("  PDG code of parent particle (0 for primary): %d\n", wcsimroottrack->GetParenttype());

                    std::cout << "  Track initial dir [unit 3-vector]: (" << wcsimroottrack->GetDir(0) << ", "
                              << wcsimroottrack->GetDir(1) << ", " << wcsimroottrack->GetDir(2) << ")"
                              << std::endl;
                    printf("  Track initial relativistic energy [MeV]: %f\n", wcsimroottrack->GetE());
                    printf("  Track initial momentum magnitude [MeV/c]: %f\n", wcsimroottrack->GetP());
                    printf("  Track mass [MeV/c2]: %f\n", wcsimroottrack->GetM());
                    printf("  Track ID: %d\n", wcsimroottrack->GetId());
                    printf("  Number of ID/OD crossings: %zu\n", wcsimroottrack->GetBoundaryPoints().size());
                    if (wcsimroottrack->GetBoundaryPoints().size() > 0)
                        printf(
                            "  First crossing position [mm]: (%f %f %f), KE [MeV]: %f, time [ns]: %f, type: %d\n",
                            wcsimroottrack->GetBoundaryPoints().at(0).at(0),
                            wcsimroottrack->GetBoundaryPoints().at(0).at(1),
                            wcsimroottrack->GetBoundaryPoints().at(0).at(2),
                            wcsimroottrack->GetBoundaryKEs().at(0), wcsimroottrack->GetBoundaryTimes().at(0),
                            wcsimroottrack->GetBoundaryTypes().at(0));
                } // verbose
            }

            // find the primary particle
            if (wcsimroottrack->GetParenttype() == 0 &&
                (abs(wcsimroottrack->GetIpnu()) == 11 || wcsimroottrack->GetIpnu() == 22)) {
                true_dir.SetXYZ(wcsimroottrack->GetDir(0), wcsimroottrack->GetDir(1), wcsimroottrack->GetDir(2));
                // red
                std::cout << "\033[1;31m";
                std::cout << "Primary track found!" << std::endl;
                std::cout << "  Track initial dir [unit 3-vector]: (" << wcsimroottrack->GetDir(0) << ", "
                          << wcsimroottrack->GetDir(1) << ", " << wcsimroottrack->GetDir(2) << ")" << std::endl;
                printf("  Track initial relativistic energy [MeV]: %f\n", wcsimroottrack->GetE());
                // white
                std::cout << "\033[0m";
                true_energy += wcsimroottrack->GetE();
                found_true_track = true;
                /* break; */
            }
        } // itrack
        if (!found_true_track) {
            std::cerr
                << "Couldn't find the true track! True energy & true direction in the output tree are set to "
                   "default values"
                << std::endl;
        }

        // Loop over triggers in the event
        for (int index = 0; index < event->GetNumberOfEvents(); index++) {
            trigger = event->GetTrigger(index);
            ncherenkovdigihits = trigger->GetNcherenkovdigihits();
            if (verbose)
                std::cout << "ncherenkovdigihits: " << ncherenkovdigihits << std::endl;
            if (ncherenkovdigihits == 0) {
                std::cout << "t, PID, 0, 0, 0, 0" << std::endl;
                continue;
            }
            std::vector<float> bsT(ncherenkovdigihits, 0);
            std::vector<float> bsQ(ncherenkovdigihits, 0);
            std::vector<int> bsCAB(ncherenkovdigihits, 0);

            bsNhit = &ncherenkovdigihits;

            // fill the tree variable
            ndigits = ncherenkovdigihits;

            // Get time, charge and PMT number for each hit
            for (i = 0; i < ncherenkovdigihits; i++) {
                TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
                WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit *>(element);

                bsT[i] = cherenkovdigihit->GetT();
                bsQ[i] = cherenkovdigihit->GetQ();
                bsCAB[i] = cherenkovdigihit->GetTubeId();

                if (bsCAB[i] == 0)
                    std::cout
                        << "Digit has tube ID 0. This shouldn't happen. WCSim TubeIds run from 1 to N. Has WCSim "
                           "changed?"
                        << std::endl;
            }

            // fit vertex position and direction using BONSAI
            int *bsCAB_a = &bsCAB[0]; // turn std::vector into a pointer for BONSAI
            float *bsT_a = &bsT[0];
            float *bsQ_a = &bsQ[0];

            if (verbose)
                std::cout << "Fitting event vertex with hk-BONSAI ..." << std::endl;

            try {
                bonsai->BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);
            } catch (int e) {
                std::cerr << "BONSAI threw an exception! Will fill the tree, then continue (not running FLOWER)"
                          << std::endl;
                out_tree->Fill();
                continue;
            }
            if (verbose) {
                /* std::cout << "Vertex found at:"; */
                for (int iv = 0; iv < 4; iv++) {
                    /* std::cout << " " << bsVertex[iv]; */
                    /* std::cout << std::endl << "BONSAI done!" << std::endl; */
                }
            }

            // Energy estimation for this trigger
            eRec = flower->GetEnergy(bsT, bsCAB, &bsVertex[0]);

            recEnergy->Fill(eRec);

            // reconstructed lepton direction (bsResult[0] is theta, bsResult[1] is phi)
            x = sin(bsResult[0]) * cos(bsResult[1]);
            y = sin(bsResult[0]) * sin(bsResult[1]);
            z = cos(bsResult[0]);

            // Print out reconstruction results in the format `time, PID, energy, direction (xyz), vertex (xyz)`
            // PID (i.e. electron vs. positron) and absolute time (as opposed to time relative to the trigger
            // window) are not available.
            /* std::cout << "t, PID, " << eRec << ", " << x << ", " << y << ", " << z << ", " << bsVertex[0] << ",
             * " */
            /* << bsVertex[1] << ", " << bsVertex[2] << std::endl; */

            // fill the reconstructed tree variables
            reco_pos.SetXYZ(bsVertex[0], bsVertex[1], bsVertex[2]);
            reco_dir.SetXYZ(x, y, z);
            reco_time = bsVertex[3];
            reco_energy = eRec;
            reco_neff = flower->RetrieveNEff();
            reco_neff2 = flower->RetrieveNEff2();

            true_charge = std::accumulate(bsQ.begin(), bsQ.end(), 0.0);

            // Free the memory used by these vectors
            std::vector<int>().swap(bsCAB);
            std::vector<float>().swap(bsT);
            std::vector<float>().swap(bsQ);

            out_tree->Fill();
        } // End of loop over triggers in event

        // reinitialize event between loops
        event->ReInitialize();
    } // End of loop over events
    timer.Print();
    std::cout << "Average times:"
              << "\tCPU   " << timer.CpuTime() / tree->GetEntries() << std::endl
              << "\tReal  " << timer.RealTime() / tree->GetEntries() << std::endl;

    outfile->cd();
    out_tree->Write();

    // display histogram(s)
    float winScale = 0.75;
    int nWide = 1;
    int nHigh = 1;
    TCanvas *c1 = new TCanvas("c1", "First canvas", 500 * nWide * winScale, 500 * nHigh * winScale);
    c1->Draw();
    c1->Divide(nWide, nHigh);
    c1->cd(1);
    recEnergy->Draw();
    c1->SaveAs("Erec.pdf");

    return 0;
}
