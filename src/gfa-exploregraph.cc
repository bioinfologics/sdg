#include <iostream>
#include <fstream>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/WorkSpace.hpp>
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "bsg-workspace"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    if (argc <2){
        std::cout<<"Please specify one of: make, log, status, dump"<<std::endl;
        exit(1);
    }

    if (0==strcmp(argv[1],"repeat-explore")) {
        std::string filename;
        try {
            cxxopts::Options options("bsg-workspace kci-dump", "BSG workspace kci-dump");

            options.add_options()
                    ("help", "Print help")
                    ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename));

            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (result.count("workspace")==0) {
                throw cxxopts::OptionException(" please specify kmer spectra file");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        WorkSpace w;
        w.load_from_disk(filename);
        //auto min_size = 700;
        int LENGTH = 1000;
        SequenceGraph& sg(w.getGraph());
        for (auto n=1; n<sg.nodes.size(); ++n){
            // All components has to be long
            if (sg.get_bw_links(n).size()==2 and
              sg.get_fw_links(n).size()==2 and
              sg.nodes[llabs(sg.get_bw_links(n)[0].dest)].sequence.size()>=LENGTH and
              sg.nodes[llabs(sg.get_bw_links(n)[1].dest)].sequence.size()>=LENGTH and
              sg.nodes[llabs(sg.get_fw_links(n)[0].dest)].sequence.size()>=LENGTH and
              sg.nodes[llabs(sg.get_fw_links(n)[1].dest)].sequence.size()>=LENGTH){
               
                std::cout << "seq"<< llabs(sg.get_fw_links(n)[0].dest) << ",seq" << llabs(sg.get_fw_links(n)[1].dest) << ",seq" << n << ",seq" << llabs(sg.get_bw_links(n)[0].dest) << ",seq" << llabs(sg.get_bw_links(n)[1].dest) <<"," << std::endl;
                std::cout << sg.nodes[llabs(sg.get_fw_links(n)[0].dest)].sequence.size() << "," << sg.nodes[llabs(sg.get_fw_links(n)[1].dest)].sequence.size() << "," << n << "," << sg.nodes[llabs(sg.get_bw_links(n)[0].dest)].sequence.size() << "," << sg.nodes[llabs(sg.get_bw_links(n)[1].dest)].sequence.size() << std::endl;
             }


//            if (w.sg.get_bw_links(n).size()==2 and w.sg.get_fw_links(n).size()==2) {
//                auto fw_check=0;
//                for (auto fn: w.sg.get_fw_links(n)){
//                    if (w.sg.nodes[llabs(fn.dest)].sequence.size() >= LENGTH) ++fw_check;
//                }
//                auto bw_check=0;
//                for (auto bn:  w.sg.get_bw_links(n)){
//                    if (w.sg.nodes[llabs(bn.dest)].sequence.size() >= LENGTH) ++bw_check;
//                }
//
//                // All tips has to be long
//                if (bw_check == w.sg.get_bw_links(n).size() and fw_check == w.sg.get_fw_links(n).size()){
//                  std::cout << n << std::endl;
//                  for (auto bn:  w.sg.get_bw_links(n))
//                      std::cout << "seq"<< abs(bn.dest) << "(" << w.sg.nodes[llabs(bn.dest)].sequence.size() << "),";
//                  std::cout << " seq"<< n << ", ";
//                  for (auto fn:  w.sg.get_fw_links(n))
//                      std::cout << "seq"<< abs(fn.dest) << "(" << w.sg.nodes[llabs(fn.dest)].sequence.size() << "),";
//                  std::cout << std::endl;
//                }
//            }
        }

    } else {
        std::cout << "LALALALA..." << std::endl;
    }

}
