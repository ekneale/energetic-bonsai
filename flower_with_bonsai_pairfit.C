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
int flower_with_bonsai_pairfit(const char *detector, const char *filename = "../wcsim.root", const char *outfiletag = "a",
                       const int verbose = 1, const bool overwrite_nearest = false,
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
    // setup reco variables for positron
    TVector3 reco_dir, reco_pos;
    double reco_pos_x,reco_pos_y,reco_pos_z;
    double reco_time;
    float reco_energy, reco_neff, reco_neff2;
    //out_tree->Branch("reco_position", &reco_pos);
    out_tree->Branch("reco_position_x", &reco_pos_x);
    out_tree->Branch("reco_position_y", &reco_pos_y);
    out_tree->Branch("reco_position_z", &reco_pos_z);
    out_tree->Branch("reco_direction", &reco_dir);
    out_tree->Branch("reco_time", &reco_time, "reco_time/D");
    out_tree->Branch("reco_energy", &reco_energy, "reco_energy/F");
    out_tree->Branch("reco_neff", &reco_neff, "reco_neff/F");
    out_tree->Branch("reco_neff2", &reco_neff2, "reco_neff2/F");

    // setup reco variables for neutron
    TVector3 reco_dir_neutron, reco_pos_neutron;
    double reco_pos_x_neutron, reco_pos_y_neutron, reco_pos_z_neutron;
    double reco_time_neutron;
    float reco_energy_neutron, reco_neff_neutron, reco_neff2_neutron;
    //out_tree->Branch("reco_position", &reco_pos);
    out_tree->Branch("reco_position_x_neutron", &reco_pos_x_neutron);
    out_tree->Branch("reco_position_y_neutron", &reco_pos_y_neutron);
    out_tree->Branch("reco_position_z_neutron", &reco_pos_z_neutron);
    out_tree->Branch("reco_direction_neutron", &reco_dir_neutron);
    out_tree->Branch("reco_time_neutron", &reco_time_neutron, "reco_time_neutron/D");
    out_tree->Branch("reco_energy_neutron", &reco_energy_neutron, "reco_energy_neutron/F");
    out_tree->Branch("reco_neff_neutron", &reco_neff_neutron, "reco_neff_neutron/F");
    out_tree->Branch("reco_neff2_neutron", &reco_neff2_neutron, "reco_neff2_neutron/F");

    //setup reco variables for combined BONSAI fit
    TVector3 reco_dir_combined, reco_pos_combined;
    double reco_pos_x_combined, reco_pos_y_combined, reco_pos_z_combined;
    double reco_time_positron_combined, reco_time_neutron_combined;
    float reco_goodness_positron, reco_goodness_neutron;
    float reco_nhits100_positron, reco_nhits100_neutron;
    float reco_energy_combined, reco_neff_combined, reco_neff2_combined;
    out_tree->Branch("reco_position_x_combined", &reco_pos_x_combined);
    out_tree->Branch("reco_position_y_combined", &reco_pos_y_combined);
    out_tree->Branch("reco_position_z_combined", &reco_pos_z_combined);
    out_tree->Branch("reco_time_positron_combined", &reco_time_positron_combined, "reco_time_positron_combined/D");
    out_tree->Branch("reco_time_neutron_combined", &reco_time_neutron_combined, "reco_time_neutron_combined/D");
    out_tree->Branch("reco_goodness_positron", &reco_goodness_positron, "reco_goodness_positron/F");
    out_tree->Branch("reco_goodness_neutron", &reco_goodness_neutron, "reco_goodness_neutron/F");
    out_tree->Branch("reco_nhits100_positron", &reco_nhits100_positron, "reco_nhits100_positron/F");
    out_tree->Branch("reco_nhits100_neutron", &reco_nhits100_neutron, "reco_nhits100_neutron/F");
    out_tree->Branch("reco_energy_combined", &reco_energy_combined, "reco_energy_combined/F");
    out_tree->Branch("reco_neff_combined", &reco_neff_combined, "reco_neff_combined/F");
    out_tree->Branch("reco_neff2_combined", &reco_neff2_combined, "reco_neff2_combined/F");

    // setup hit variables for positron
    int nhits, ndigits;
    out_tree->Branch("nhits", &nhits, "nhits/I");
    out_tree->Branch("ndigits", &ndigits, "ndigits/I");

    //setup hit variables for neutron
    int nhits_neutron, ndigits_neutron;
    out_tree->Branch("nhits_neutron", &nhits_neutron, "nhits_neutron/I");
    out_tree->Branch("ndigits_neutron", &ndigits_neutron, "ndigits_neutron/I");

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
    // for neutron variables
    double x_neutron, y_neutron, z_neutron, eRec_neutron;
    //combined variables
    double x_combined, y_combined, z_combined, eRec_combined;

    // used by bonsai->BonsaiFit()
    float bsVertex[4], bsResult[6], bsGood[3];
    int *bsNhit;
    int bsNsel[2];

    //define neutron bsNhit
    float bsVertex_neutron[4], bsResult_neutron[6], bsGood_neutron[3];
    int *bsNhit_neutron;
    int bsNsel_neutron[2];

    //define BONSAIPairFit
    float pair_vertex[5], pair_goodness[2], pair_dir[2];
    int pair_nhits100[2];

    // Loop over events
    std::cout << "Starting timer for the event loop" << std::endl;
    TStopwatch timer;

    // for (int ev = start_event; ev < nWCSimEvents; ev++) {
     for (int ev = start_event; ev < 100; ev++) {
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
	reco_pos_x= -9999;
	reco_pos_y= -9999;
	reco_pos_z = -9999;
        reco_dir.SetXYZ(-99, -99, -99);
        reco_time = -9999;
        reco_energy = -99;
        reco_neff = -99;
        reco_neff2 = -99;
        nhits = -99;
        ndigits = -99;
	//tree variables for neutron
	   
        ndigits_neutron = -99;
        reco_pos_neutron.SetXYZ(-9999, -9999, -9999);
        reco_pos_x_neutron = -9999;
        reco_pos_y_neutron = -9999;
        reco_pos_z_neutron = -9999;
        reco_dir_neutron.SetXYZ(-99, -99, -99);
        reco_time_neutron = -9999;
        reco_energy_neutron = -99;
        reco_neff_neutron = -99;
        reco_neff2_neutron = -99;
        nhits_neutron = -99;

	//tree variables for combined BONSAI fit
	
	reco_pos_combined.SetXYZ(-9999, -9999, -9999);
	reco_pos_x_combined = -9999;
	reco_pos_y_combined = -9999;
	reco_pos_z_combined = -9999;
	reco_dir_combined.SetXYZ(-99, -99, -99);
	reco_time_neutron_combined = -9999;
	reco_time_positron_combined = -9999;
	reco_energy_combined = -99;
	reco_neff_combined = -99;
	reco_neff2_combined = -99;
	

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
        nhits_neutron = trigger->GetNcherenkovhits();

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
                if (verbose) {
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
                (abs(wcsimroottrack->GetIpnu()) == 11 )) {
                true_dir.SetXYZ(wcsimroottrack->GetDir(0), wcsimroottrack->GetDir(1), wcsimroottrack->GetDir(2));
                // red
                std::cout << "\033[1;31m";
                std::cout << "Primary track found!" << std::endl;
                std::cout << "  Track initial dir [unit 3-vector]: (" << wcsimroottrack->GetDir(0) << ", "
                          << wcsimroottrack->GetDir(1) << ", " << wcsimroottrack->GetDir(2) << ")" << std::endl;
                printf("  Track initial relativistic energy [MeV]: %f\n", wcsimroottrack->GetE());
                // white
                std::cout << "\033[0m";
                true_energy = wcsimroottrack->GetE();
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
	    std::cout << "true_energy " << true_energy << std::endl;
	    
	    if (true_energy > 50) continue;
	    
        // Loop over triggers in the event
        for (int index = 0; index < 1; index++) {
	    if (verbose)
	        std::cout << "trigger " << index << std::endl;
            trigger = event->GetTrigger(index);
            ncherenkovdigihits = trigger->GetNcherenkovdigihits();
            if (verbose)
                std::cout << "ncherenkovdigihits: " << ncherenkovdigihits << std::endl;
            if (ncherenkovdigihits == 0) {
                std::cout << "t, PID, 0, 0, 0, 0" << std::endl;
                continue;
            }
	    
	    //TODO temporarily skip higher-energy events
	    if (ncherenkovdigihits >2000) continue;
            
            //vectors for positron
            std::vector<float> bsT;
            std::vector<float> bsQ;
            std::vector<int> bsCAB;

	    //neutron vectors
	    std::vector<float> neutronT;
	    std::vector<float> neutronQ;
	    std::vector<int> neutronCAB;

            //bsNhit = &ncherenkovdigihits;

            // fill the tree variable
            //ndigits = ncherenkovdigihits;

            // Get time, charge and PMT number for each hit
	    
	    
            for (i = 0; i < ncherenkovdigihits; i++) {
                TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
                WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit *>(element);

		//using float and int to get the time, charge, and tube ID of the hit

		float time = cherenkovdigihit->GetT();
		float charge = cherenkovdigihit->GetQ();
		int tubeId = cherenkovdigihit->GetTubeId();

		if (i < 20) {
		  std::cout << "Hit" << i << "Time:" << time << "ns" << std::endl;
		}
		

		//skips hits that are 5000ns or more (filtering for positron hits)
		//so if hit time is >/= to 5000ns, the condition will skip to filter out late hits
		if (time < 2000) {

		  //update the bsT, bsQ, and bsCAB vectors with the filtered hits
                 //using .push_back command to add elements to the vector

                 bsT.push_back(time);
                 bsQ.push_back(charge);
                 bsCAB.push_back(tubeId);
		}

		// filtering neutron hits in the same way as positrons
		if (time >= 2000) {
		  neutronT.push_back(time);
		  neutronQ.push_back(charge);
		  neutronCAB.push_back(tubeId);
		}

		if (tubeId == 0) {
		  std::cout
		    <<"Digit has tube ID 0. This shouldn't happen. WCSim tube runs from 1 to N. Has WCSim changed?"
		    <<std::endl;
		}
	    }

	    //count the number of filtered hits (so less than 5000 ns)

	    int ncherenkovdigihits_filtered = bsT.size();
	    bsNhit = &ncherenkovdigihits_filtered;
	    ndigits = ncherenkovdigihits_filtered;

	    //count neutron filtered hits
	    int ncherenkovdigihits_neutron = neutronT.size();
	    bsNhit_neutron = &ncherenkovdigihits_neutron;
	    ndigits_neutron = ncherenkovdigihits_neutron;
	    

            // fit vertex position and direction using BONSAI
	   int *bsCAB_a = &bsCAB[0]; // turn std::vector into a pointer for BONSAI
           float *bsT_a = &bsT[0];
           float *bsQ_a = &bsQ[0];

	   //sliding window for neutron
	   //finding 200 ns window with most hits
	   if (neutronT.empty()) {
	     std::cout << "No neutron hits."
		       <<std::endl;
	     continue;
	   }
	   
	   int max_nhits = 0;
	   size_t max_start = 0;
	   size_t max_end = 0;

	   for (size_t start = 0; start<neutronT.size()-2; start++) {
	     for (size_t end = start + 1 ; end<neutronT.size()-1; end++) {
	       
	       float dt = neutronT[end]-neutronT[start];
	       if (dt >= 200) {
		 
		 int tmp_nhits = end - start;

		 if (tmp_nhits > max_nhits) {

		     max_nhits = tmp_nhits;
		     max_start = start;
		     max_end = end;
		   }
		 break;
	       }
	     }
	   }
	   

	     //define vectors to extract updated hits 
	     std::vector<float>neutronT_window;
	     std::vector<float>neutronQ_window;
	     std::vector<int>neutronCAB_window;

	     for (size_t i = max_start; i < max_end; ++i) {
	       
		 //push back command to add elements to vectors
		 neutronT_window.push_back(neutronT[i]);
		 neutronQ_window.push_back(neutronQ[i]);
		 neutronCAB_window.push_back(neutronCAB[i]);
	     }

	     //update neutron vectors 
	     neutronT = neutronT_window;
	     neutronQ = neutronQ_window;
	     neutronCAB = neutronCAB_window;

	     //update hits and vectors to BONSAI function
	     ncherenkovdigihits_neutron = neutronT.size();
	     bsNhit_neutron = &ncherenkovdigihits_neutron;
	     ndigits_neutron = ncherenkovdigihits_neutron;
	       
	     float *neutronT_a = &neutronT_window[0];
	     float *neutronQ_a = &neutronQ_window[0];
	     int *neutronCAB_a = &neutronCAB_window[0];

	     //add message to say no neutron hits
	     if (ncherenkovdigihits_neutron == 0) {
	       std::cout << "No neutron hits in 200 ns window"
			 <<std::endl;
	     }
	     
		 
		 

	    //vertex for neutron vectors
	    //int *neutronCAB_a = &neutronCAB[0];
	    //float *neutronT_a = &neutronT[0];
	    //float *neutronQ_a = &neutronQ[0];


	    
        //bonsai fit for positron
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

            //bonsai fit for neutron
            if (verbose)
                std::cout << "Fitting event vertex with hk-BONSAI ..." << std::endl;

            try {
                bonsai->BonsaiFit(bsVertex_neutron, bsResult_neutron, bsGood_neutron, bsNsel_neutron, bsNhit_neutron, neutronCAB_a, neutronT_a, neutronQ_a);
            }
            catch (int e) {
                std::cerr << "BONSAI threw an exception! Will fill the tree, then continue (not running FLOWER)"
                    << std::endl;
                out_tree->Fill();
                continue;
            }

	    //BONSAI combined fit
	    if (verbose)
                std::cout << "Fitting event vertex with hk-BONSAI ..." << std::endl;

            try {
	      bonsai->BonsaiPairFit(pair_vertex, pair_goodness, pair_dir, pair_nhits100, bsNhit, bsCAB_a, bsT_a, bsQ_a, bsNhit_neutron, neutronCAB_a, neutronT_a, neutronQ_a);
            }
            catch (int e) {
                std::cerr << "BONSAI PairFit threw an exception! Will fill the tree, then continue (not running FLOWER)"
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

	    std::cout << "Reconstructed event " << ev << ", trigger " << index << std::endl;

            // Energy estimation for this trigger
            eRec = flower->GetEnergy(bsT, bsCAB, &bsVertex[0]);

            recEnergy->Fill(eRec);

            // reconstructed lepton direction (bsResult[0] is theta, bsResult[1] is phi)
            x = sin(bsResult[0]) * cos(bsResult[1]);
            y = sin(bsResult[0]) * sin(bsResult[1]);
            z = cos(bsResult[0]);

            // filling reco energy for neutron
	
            eRec_neutron = flower->GetEnergy(neutronT, neutronCAB, &bsVertex_neutron[0]);
	    reco_energy_neutron = eRec_neutron;
	    reco_neff_neutron = flower->RetrieveNEff();
	    reco_neff2_neutron = flower->RetrieveNEff2();

            // reconstructed lepton direction (bsResult[0] is theta, bsResult[1] is phi) (for neutron)
            x_neutron = sin(bsResult_neutron[0]) * cos(bsResult_neutron[1]);
            y_neutron = sin(bsResult_neutron[0]) * sin(bsResult_neutron[1]);
            z_neutron = cos(bsResult_neutron[0]);

	    //filling reco energy for combined fit
	    eRec_combined = flower->GetEnergy(bsT, bsCAB, &pair_vertex[0]);
	    reco_energy_combined = eRec_combined;
	    reco_neff_combined = flower->RetrieveNEff();
	    reco_neff2_combined = flower->RetrieveNEff2();

	    //reconstruction lepton direction for combined fit
	    x_combined = sin(pair_dir[0]) * cos(pair_dir[1]);
	    y_combined = sin(pair_dir[0]) * sin(pair_dir[1]);
	    z_combined = cos(pair_dir[0]);


            // Print out reconstruction results in the format `time, PID, energy, direction (xyz), vertex (xyz)`
            // PID (i.e. electron vs. positron) and absolute time (as opposed to time relative to the trigger
            // window) are not available.
            /* std::cout << "t, PID, " << eRec << ", " << x << ", " << y << ", " << z << ", " << bsVertex[0] << ",
             * " */
            /* << bsVertex[1] << ", " << bsVertex[2] << std::endl; */

            // fill the reconstructed tree variables
            reco_pos.SetXYZ(bsVertex[0], bsVertex[1], bsVertex[2]);
	    reco_pos_x = bsVertex[0];
	    reco_pos_y = bsVertex[1];
	    reco_pos_z = bsVertex[2];
            reco_dir.SetXYZ(x, y, z);
            reco_time = bsVertex[3];
            reco_energy = eRec;
            reco_neff = flower->RetrieveNEff();
            reco_neff2 = flower->RetrieveNEff2();

            //fill reconstructed tree variables for neutron
            reco_pos_neutron.SetXYZ(bsVertex_neutron[0], bsVertex_neutron[1], bsVertex_neutron[2]);
            reco_pos_x_neutron = bsVertex_neutron[0];
            reco_pos_y_neutron = bsVertex_neutron[1];
            reco_pos_z_neutron = bsVertex_neutron[2];
            reco_dir_neutron.SetXYZ(x_neutron, y_neutron, z_neutron);
            reco_time_neutron = bsVertex_neutron[3];
            //reco_energy_neutron = eRec_neutron;
            //reco_neff_neutron = flower->RetrieveNEff();
            //reco_neff2_neutron = flower->RetrieveNEff2();

	    //fill reconstructed tree variables for BONSAI combined fit
	    reco_pos_combined.SetXYZ(pair_vertex[0], pair_vertex[1], pair_vertex[2]);
	    reco_pos_x_combined = pair_vertex[0];
	    reco_pos_y_combined = pair_vertex[1];
	    reco_pos_z_combined = pair_vertex[2];
	    reco_time_positron_combined = pair_vertex[3];
	    reco_time_neutron_combined = pair_vertex[4];
	    reco_goodness_positron = pair_goodness[0];
	    reco_goodness_neutron = pair_goodness[1];
	    reco_nhits100_positron = pair_nhits100[0];
	    reco_nhits100_neutron = pair_nhits100[1];
	    reco_dir_combined.SetXYZ(x_combined, y_combined, z_combined);

            true_charge = std::accumulate(bsQ.begin(), bsQ.end(), 0.0);

            // Free the memory used by these vectors
            std::vector<int>().swap(bsCAB);
            std::vector<float>().swap(bsT);
            std::vector<float>().swap(bsQ);

            // Free memory used by neutron vectors
            std::vector<int>().swap(neutronCAB);
            std::vector<float>().swap(neutronT);
            std::vector<float>().swap(neutronQ);

	    std::cout << "Filling event " << ev << ", trigger " << index << std::endl;
            out_tree->Fill();
	    std::cout << "Filled tree for event " << ev << ", trigger " << index << std::endl;
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
