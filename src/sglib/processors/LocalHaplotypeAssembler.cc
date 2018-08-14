//
// Created by Bernardo Clavijo (EI) on 18/07/2018.
//

#include "LocalHaplotypeAssembler.hpp"
#include "GraphMaker.hpp"

void print_seq_kmer_location(const char * seq, std::unordered_map<uint64_t, graphPosition> & index){
    std::cout<<"Kmer to node analysis for "<<seq<<std::endl;
    std::vector<std::pair<uint64_t,bool>> readkmers;
    CStringKMerFactory cskf(31);
    cskf.create_kmers_direction(readkmers,seq);
    for (auto rki=readkmers.begin();rki<readkmers.end();++rki) {
        auto kp=index.find(rki->first);
        if (kp==index.end()) {std::cout <<" *";continue;} //kmer not found
        auto node=(rki->second ? -kp->second.node:kp->second.node);
        std::cout<<" "<<node;
    }
    std::cout<<std::endl;

}

void LocalHaplotypeAssembler::init_from_backbone( std::vector<sgNodeID_t> _backbone) {
    backbone=_backbone;
    //std::cout<<"Creating a LocalHaplotypeAssembler instance from backbone"<<std::endl;
    //std::cout<<"Backbone nodes:";
    //for (auto n:backbone) std::cout<<" "<<n;
    //std::cout<<std::endl;
    for (auto n:backbone) backbone_nodes.emplace_back(ws.sg.nodes[llabs(n)].sequence);

    //std::cout<<"Filling Candidate Tag Set..."<<std::endl;
    //Get tag reads in the nodes, and how many nodes with reads counts.
    std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
    for (auto &ln:backbone) {
        //partial counts in a single node
        std::map<bsg10xTag ,uint32_t> ntagcounts;
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(ln)]){
            auto tag=ws.linked_read_datastores[0].get_read_tag(rm.read_id);
            if (tag==0) continue;
            ++ntagcounts[tag];
        }
        //added to totals
        for (auto ntc:ntagcounts) {
            ++tagcounts[ntc.first].first;
            tagcounts[ntc.first].second+=ntc.second;
        }
    }
    uint64_t total_reads=0;
    for (auto tc:tagcounts) {
        auto & tag=tc.first;
        auto & tag_backbone_nodes=tc.second.first;
        auto & tag_backbone_reads=tc.second.second;
        auto tag_total_reads=ws.linked_read_datastores[0].get_tag_reads(tc.first).size();//TODO: this can be optimised by not creating the vector
        if (tc.second.first>1 or tag_total_reads/50<tag_backbone_reads and tag_total_reads>10) {
            tagSet.insert(tc.first);
            total_reads+=tag_total_reads;
        }
    }
    //std::cout<<"Local tags: "<<tagSet.size()<<" reads: "<<total_reads<<std::endl;
    //std::cout<<"Creating set of relevant paired reads..."<<std::endl;
    for (auto prl=0;prl<ws.paired_read_mappers.size();++prl) {
        paired_reads.emplace_back(std::make_pair(prl,std::vector<uint64_t>()));
        for (auto &ln:backbone) {
            auto nreads=ws.paired_read_mappers[prl].get_node_readpairs_ids(ln);
            paired_reads.back().second.insert(paired_reads.back().second.end(),nreads.begin(),nreads.end());
        }
        if (paired_reads.back().second.empty()) paired_reads.pop_back();
        //else std::cout<<paired_reads.back().second.size()<<" reads from "<<ws.paired_read_datastores[prl].filename<<std::endl;
    }
    for (auto lrl = 0; ws.long_read_mappers.size(); ++lrl) {
        long_reads.emplace_back(std::make_pair(lrl, std::vector<uint64_t>()));
        for (auto &ln:backbone) {
            auto nreads = ws.long_read_mappers[lrl].get_node_read_ids(ln);
            long_reads.back().second.insert(long_reads.back().second.end(), nreads.begin(), nreads.end());
        }
        if (long_reads.back().second.empty()) long_reads.pop_back();
    }
    //std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}

void LocalHaplotypeAssembler::init_from_file(std::string problem_file) {
    std::cout<<"Loading LocalHaplotypeAssembler problem from file: "<<problem_file<<std::endl;

    std::ifstream input_file(problem_file);
    uint64_t count;

    //load backbone;
    input_file.read((char *)&count,sizeof(count));
    backbone.resize(count);
    input_file.read((char *)backbone.data(),count*sizeof(backbone[0]));
    for (auto n:backbone) backbone_nodes.emplace_back(ws.sg.nodes[llabs(n)].sequence);

    //load 10x tags;
    input_file.read((char *)&count,sizeof(count));
    bsg10xTag t;
    while (count--) {
        input_file.read((char *) &t, sizeof(t));
        tagSet.insert(t);
    }

    //load paired read ids;
    input_file.read((char *)&count,sizeof(count));
    paired_reads.resize(count);
    for (auto &p:paired_reads){
        input_file.read((char *)&p.first,sizeof(p.first));
        input_file.read((char *)&count,sizeof(count));
        p.second.resize(count);
        input_file.read((char *)p.second.data(),count*sizeof(p.second[0]));
    }

    std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}

uint64_t LocalHaplotypeAssembler::expand_canonical_repeats_direct(int max_rep_size) {
    //A simpler canonical repeat expander: first look for the repeats, then find reads that contain them, then use only those
    std::cout<<"starting repeat expansion by direct read matching";
    //Find the repeat nodes
    std::vector<std::pair<sgNodeID_t, Node>> repeats;
    for (auto n=1;n<assembly.nodes.size();++n) {
        if (assembly.nodes[n].status==sgNodeDeleted or assembly.nodes[n].sequence.size()>max_rep_size) continue;
        auto bwl=assembly.get_bw_links(n);
        auto fwl=assembly.get_fw_links(n);
        if (bwl.size()!=2 or fwl.size()!=2) continue;
        //TODO: removing this may even solve some small loops.
        std::set<sgNodeID_t> nset;
        nset.insert(n);
        nset.insert(llabs(bwl[0].dest));
        nset.insert(llabs(bwl[1].dest));
        nset.insert(llabs(fwl[0].dest));
        nset.insert(llabs(fwl[1].dest));
        if (nset.size()<5) continue;
        repeats.emplace_back(n,assembly.nodes[n]);
    }
    std::cout<<"There are "<<repeats.size()<<" repeats to analyse"<<std::endl;
    //TODO: create strings for direct and rc versions of the repeat, and of all 4 possible resolutions, create a structure that holds those.
    std::vector<std::pair<sgNodeID_t, std::array<std::string,8>>> repeat_posibilities;
    for (auto &r:repeats) {
        //find the previous nucleotides
        auto bwl = assembly.get_bw_links(r.first);
        auto pn1 = assembly.nodes[llabs(bwl[0].dest)];
        if (bwl[0].dest > 0) pn1.make_rc();
        auto pb1 = pn1.sequence[pn1.sequence.size() + bwl[0].dist - 1];//XXX: assumes negative distance.
        auto pn2 = assembly.nodes[llabs(bwl[1].dest)];
        if (bwl[1].dest > 0) pn2.make_rc();
        auto pb2 = pn2.sequence[pn2.sequence.size() + bwl[1].dist - 1];//XXX: assumes negative distance.
        //find the next nucleotides
        auto fwl=assembly.get_fw_links(r.first);
        auto nn1=assembly.nodes[llabs(fwl[0].dest)];
        if (fwl[0].dest<0) nn1.make_rc();
        auto nb1=nn1.sequence[-fwl[0].dist];//XXX: assumes negative distance.
        auto nn2=assembly.nodes[llabs(fwl[1].dest)];
        if (fwl[1].dest<0) nn2.make_rc();
        auto nb2=nn2.sequence[-fwl[1].dist];//XXX: assumes negative distance.
        std::string saa,sab,sba,sbb;
        saa.push_back(pb1);
        saa+=r.second.sequence;
        saa.push_back(nb1);
        sab.push_back(pb1);
        sab+=r.second.sequence;
        sab.push_back(nb2);
        sba.push_back(pb2);
        sba+=r.second.sequence;
        sba.push_back(nb1);
        sbb.push_back(pb2);
        sbb+=r.second.sequence;
        sbb.push_back(nb2);
        auto saaN=Node(saa);
        saaN.make_rc();
        auto saar=saaN.sequence;
        auto sabN=Node(sab);
        sabN.make_rc();
        auto sabr=sabN.sequence;
        auto sbaN=Node(sba);
        sbaN.make_rc();
        auto sbar=sbaN.sequence;
        auto sbbN=Node(sbb);
        sbbN.make_rc();
        auto sbbr=sbbN.sequence;
        std::array<std::string,8> rp={saa,sab,sba,sbb,saar,sabr,sbar,sbbr};
        repeat_posibilities.emplace_back(std::make_pair(r.first,rp));
    }
    std::cout<<"Collecting votes..."<<std::endl;
    std::vector<uint64_t[4]> rvotes(repeat_posibilities.size());
    //For every single read: look
    for (auto t:tagSet) {
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            auto rseq=ws.linked_read_datastores[0].get_read_sequence(rid);
            for (auto i=0;i<repeat_posibilities.size();++i){
                for (auto p=0;p<8;++p) {
                    if (rseq.find(repeat_posibilities[i].second[p])<rseq.size()) ++rvotes[i][p%4];
                }
            }
        }
    }
    std::cout<<"Votes: "<<std::endl;
    for (auto i=0;i<repeat_posibilities.size();++i) {
        std::cout<<"Repead on node "<<repeat_posibilities[i].first<<" "<<rvotes[i][0]<<" "<<rvotes[i][1]<<" "<<rvotes[i][2]<<" "<<rvotes[i][3]<<std::endl;
    }

}

uint64_t LocalHaplotypeAssembler::expand_canonical_repeats() {
    std::vector<std::pair<sgNodeID_t ,std::pair<std::vector<std::vector<sgNodeID_t>>,std::vector<std::vector<sgNodeID_t>>>>> to_expand;
    const int min_cov_exp=5;
    const float min_diff_exp=3;
    for (auto n=0;n<assembly.nodes.size();++n){
        if (assembly.nodes[n].status==sgNodeDeleted) continue;
        auto bwl=assembly.get_bw_links(n);
        auto fwl=assembly.get_fw_links(n);
        if (bwl.size()!=2 or fwl.size()!=2) continue;
        std::set<sgNodeID_t> nset;
        nset.insert(n);
        nset.insert(llabs(bwl[0].dest));
        nset.insert(llabs(bwl[1].dest));
        nset.insert(llabs(fwl[0].dest));
        nset.insert(llabs(fwl[1].dest));
        if (nset.size()<5) continue;
        int v00=0,v11=0,v01=0,v10=0;

        for (auto p:linkedread_paths){
            if (p.size()<2) continue;
            //std::cout<<"analising path";
            int pb0=0,pb1=0,pn=0,pf0=0,pf1=0;
            for (auto &nip:p) {
                //std::cout<<" "<<nip;
                if (nip==-bwl[0].dest) pb0=1;
                else if (nip==bwl[0].dest) pb0=-1;
                else if (nip==-bwl[1].dest) pb1=1;
                else if (nip==bwl[1].dest) pb1=-1;
                else if (nip==n) pn=1;
                else if (nip==-n) pn=-1;
                else if (nip==fwl[0].dest) pf0=1;
                else if (nip==-fwl[0].dest) pf0=-1;
                else if (nip==fwl[1].dest) pf1=1;
                else if (nip==-fwl[1].dest) pf1=-1;
            }
            //std::cout<<std::endl;
            if (pb0==pf0 and pb0!=0 and pb1==0 and pf1==0) ++v00;
            if (pb0==pf1 and pb0!=0 and pb1==0 and pf0==0) ++v01;
            if (pb1==pf0 and pb1!=0 and pb0==0 and pf1==0) ++v10;
            if (pb1==pf1 and pb1!=0 and pb0==0 and pf0==0) ++v11;
            //if (pn!=0) std::cout<<"path with p, node votes "<<pb0<<" "<<pb1<<" "<<pn<<" "<<pf0<<" "<<pf1<<std::endl;
        }
        int lv00=0,lv11=0,lv01=0,lv10=0;
        for (auto p:pairedread_paths){
            if (p.size()<2) continue;
            //std::cout<<"analising path";
            int pb0=0,pb1=0,pn=0,pf0=0,pf1=0;
            for (auto &nip:p) {
                //std::cout<<" "<<nip;
                if (nip==-bwl[0].dest) pb0=1;
                else if (nip==bwl[0].dest) pb0=-1;
                else if (nip==-bwl[1].dest) pb1=1;
                else if (nip==bwl[1].dest) pb1=-1;
                else if (nip==n) pn=1;
                else if (nip==-n) pn=-1;
                else if (nip==fwl[0].dest) pf0=1;
                else if (nip==-fwl[0].dest) pf0=-1;
                else if (nip==fwl[1].dest) pf1=1;
                else if (nip==-fwl[1].dest) pf1=-1;
            }
            //std::cout<<std::endl;
            if (pb0==pf0 and pb0!=0 and pb1==0 and pf1==0) ++lv00;
            if (pb0==pf1 and pb0!=0 and pb1==0 and pf0==0) ++lv01;
            if (pb1==pf0 and pb1!=0 and pb0==0 and pf1==0) ++lv10;
            if (pb1==pf1 and pb1!=0 and pb0==0 and pf0==0) ++lv11;
            //if (pn!=0) std::cout<<"path with p, node votes "<<pb0<<" "<<pb1<<" "<<pn<<" "<<pf0<<" "<<pf1<<std::endl;
        }
//        std::cout<<"Repeat at node "<<n<<"( "<<assembly.nodes[n].sequence.size()<<"bp ), votes: "
//                  <<v00<<" "<<v11<<" "<<v10<<" "<<v01<<"  ("<<lv00<<" "<<lv11<<" "<<lv10<<" "<<lv01<<")"<<std::endl;
        if (v00>=min_cov_exp and v11>=min_cov_exp and std::min(v00,v11)>=min_diff_exp*std::max(v10,v01)) {
//            std::cout<<"solved AA!"<<std::endl;
            to_expand.push_back(std::make_pair(n,std::make_pair(std::vector<std::vector<sgNodeID_t>>({{bwl[0].dest},{bwl[1].dest}}),std::vector<std::vector<sgNodeID_t>>({{fwl[0].dest},{fwl[1].dest}}))));
        }
        if (v01>=min_cov_exp and v10>=min_cov_exp and std::min(v01,v10)>=min_diff_exp*std::max(v00,v11)) {
//            std::cout<<"solved AB!"<<std::endl;
            to_expand.push_back(std::make_pair(n,std::make_pair(std::vector<std::vector<sgNodeID_t>>({{bwl[0].dest},{bwl[1].dest}}),std::vector<std::vector<sgNodeID_t>>({{fwl[1].dest},{fwl[0].dest}}))));
        }
    }
    for (auto r:to_expand) assembly.expand_node(r.first,r.second.first,r.second.second);
    return to_expand.size();
}

uint64_t LocalHaplotypeAssembler::unroll_short_loops() {
    for (auto n = 0; n < assembly.nodes.size(); ++n) {
        if (assembly.nodes[n].status == sgNodeDeleted) continue;
        auto bwl = assembly.get_bw_links(n);
        auto fwl = assembly.get_fw_links(n);
        if (bwl.size() != 2 or fwl.size() != 2) continue;
        sgNodeID_t prev, next, loop;
        if (bwl[0].dest == -fwl[0].dest) {
            prev = -bwl[1].dest;
            loop = fwl[0].dest;
            next = fwl[1].dest;
        } else if (bwl[0].dest == -fwl[1].dest) {
            prev = -bwl[1].dest;
            loop = fwl[1].dest;
            next = fwl[0].dest;
        } else if (bwl[1].dest == -fwl[0].dest) {
            prev = -bwl[0].dest;
            loop = fwl[0].dest;
            next = fwl[1].dest;
        } else if (bwl[1].dest == -fwl[1].dest) {
            prev = -bwl[0].dest;
            loop = fwl[1].dest;
            next = fwl[0].dest;
        } else continue;
        if (assembly.get_fw_links(loop).size() != 1 or assembly.get_bw_links(loop).size() != 1) continue;
        //std::cout << "Loop detected on " << prev << " " << n << " " << loop << " " << n << " " << next << std::endl;
        int64_t rid = -1;
        for (auto &p:linkedread_paths) {
            rid += 2;
            if (p.size() < 3) continue;
            bool sprev=false,snext=false;
            for (auto &pp:p) {
                if (llabs(pp) == llabs(prev)) sprev = true;
                if (llabs(pp) == llabs(next)) snext = true;
                if(snext and sprev)
                {
                    for (auto &ppp:p) std::cout << ppp << " ";
                    std::cout << std::endl;
                    //print_seq_kmer_location(ws.linked_read_datastores[0].get_read_sequence(rid).c_str(),
                    //                        assembly.kmer_to_graphposition);
                    //print_seq_kmer_location(ws.linked_read_datastores[0].get_read_sequence(rid + 1).c_str(),
                    //                        assembly.kmer_to_graphposition);


                    continue;
                }
            }
        }
    }
}

void LocalHaplotypeAssembler::assemble(int k, int min_cov, bool tag_cov, bool simplify, std::string output_prefix){
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
    auto ltkmers128 = ws.linked_read_datastores[0].get_tags_kmers128(k, min_cov, tagSet, blrsg, tag_cov);
    GraphMaker gm(assembly);
    gm.new_graph_from_kmerset_trivial128(ltkmers128, k);
    gm.tip_clipping(200);
    if (simplify) {
        gm.remove_small_unconnected(500);
        //path_all_reads();
        assembly.create_index();
        path_linked_reads();
        if (!output_prefix.empty()) assembly.write_to_gfa(output_prefix + "pre_repex.gfa");
        while (expand_canonical_repeats() > 0) {
            assembly.join_all_unitigs();
            //path_all_reads();
            assembly.create_index();
            path_linked_reads();
        }
        //unroll_short_loops();
    }
}



void add_readkmer_nodes_lha(std::vector<sgNodeID_t> & kmernodes, std::vector<std::pair<uint64_t,bool>> & readkmers, std::unordered_map<uint64_t, graphPosition> & index, bool rev){
    //TODO allow for a minimum of kmers to count the hit?
    if (not rev) {
        for (auto rki=readkmers.begin();rki<readkmers.end();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? -kp->second.node:kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }
    else {
        for (auto rki=readkmers.rbegin();rki<readkmers.rend();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? kp->second.node:-kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }

}

void add_readkmer_nodes_lha128(std::vector<sgNodeID_t> & kmernodes, std::vector<std::pair<__uint128_t,bool>> & readkmers, std::unordered_map<__uint128_t, graphPosition> & index, bool rev){
    //TODO allow for a minimum of kmers to count the hit?
    if (not rev) {
        for (auto rki=readkmers.begin();rki<readkmers.end();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? -kp->second.node:kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }
    else {
        for (auto rki=readkmers.rbegin();rki<readkmers.rend();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? kp->second.node:-kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }

}


void LocalHaplotypeAssembler::path_linked_reads() {
    linkedread_paths.clear();
    linkedread_paths.reserve(1000000);//TODO: do this better!!!!
    //now populate the linked read paths first

    CStringKMerFactory cskf(31);
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
    std::vector<std::pair<uint64_t,bool>> read1kmers,read2kmers;
    std::vector<sgNodeID_t> kmernodes;

    //BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lib], 1000000, 1000);
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],100000,1000);
    for (auto t:tagSet) {
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;

            read1kmers.clear();
            read2kmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(read1kmers, blrsg.get_read_sequence(rid));
            cskf.create_kmers_direction(read2kmers, blrsg.get_read_sequence(rid + 1));
            //first put the kmers from read 1 in there;
            add_readkmer_nodes_lha(kmernodes, read1kmers, assembly.kmer_to_graphposition, false);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
            add_readkmer_nodes_lha(kmernodes, read2kmers, assembly.kmer_to_graphposition, true);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;

            linkedread_paths.emplace_back(kmernodes);
        }
    }
    sglib::OutputLog()<<linkedread_paths.size()<<" linked-reads paths created!"<<std::endl;

}


void LocalHaplotypeAssembler::path_linked_reads_informative_singles() {
    linkedread_paths.clear();
    linkedread_paths.reserve(1000000);//TODO: do this better!!!!
    //now populate the linked read paths first

    CStringKMerFactory128 cskf(63);
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
    std::vector<std::pair<__uint128_t,bool>> readkmers;
    std::vector<sgNodeID_t> kmernodes;
    //std::ofstream crf("chimeric_linkedreads.fasta");
    //BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lib], 1000000, 1000);
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],100000,1000);
    for (auto t:tagSet) {
//        auto chim=0;
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            readkmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(readkmers, blrsg.get_read_sequence(rid));

            add_readkmer_nodes_lha128(kmernodes, readkmers, assembly.k63mer_to_graphposition, false);
            if (kmernodes.size()>1) linkedread_paths.emplace_back(kmernodes);
//            if (kmernodes.size()==2 and kmernodes[0]==-kmernodes[1]) { //TODO quantify this by tag,
//                //std::cout<<"Read #"<<rid<<" has incoherent mapping kmernodes[0]==-32152 and kmernodes[1]==32152"<<std::endl;
//                //print_seq_kmer_location(blrsg.get_read_sequence(rid),assembly.kmer_to_graphposition);
//                crf<<">read_"<<rid<<"_tag_"<< t <<"_"<<llabs(kmernodes[0])<<std::endl<<blrsg.get_read_sequence(rid)<<std::endl;
//                ++chim;
//            }
        }
        //std::cout<<t<<","<<chim<<","<<ws.linked_read_datastores[0].get_tag_reads(t).size()<<","<<100.0*chim/ws.linked_read_datastores[0].get_tag_reads(t).size()<<std::endl;
    }
    //sglib::OutputLog()<<linkedread_paths.size()<<" informative single read paths created!"<<std::endl;
}

void LocalHaplotypeAssembler::path_paired_reads_informative_singles() {
    pairedread_paths.clear();
    pairedread_paths.reserve(1000000);//TODO: do this better!!!!
    //now populate the linked read paths first

    CStringKMerFactory cskf(31);
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
    std::vector<std::pair<uint64_t,bool>> readkmers;
    std::vector<sgNodeID_t> kmernodes;
    //std::ofstream crf("chimeric_pairedreads.fasta");
    //BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lib], 1000000, 1000);
    for (auto lpr:paired_reads) {
//        auto chim=0;
        BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lpr.first],100000,1000);
        for (auto rid : lpr.second) {
            readkmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(readkmers, bprsg.get_read_sequence(rid));

            add_readkmer_nodes_lha(kmernodes, readkmers, assembly.kmer_to_graphposition, false);
            if (kmernodes.size()>1) pairedread_paths.emplace_back(kmernodes);
//            if (kmernodes.size()==2 and kmernodes[0]==-kmernodes[1]) { //TODO quantify this by tag, then try to understand WTF is going on here. maybe tags from another haplotype?
//                //std::cout<<"Read #"<<rid<<" has incoherent mapping kmernodes[0]==-32152 and kmernodes[1]==32152"<<std::endl;
//                //print_seq_kmer_location(blrsg.get_read_sequence(rid),assembly.kmer_to_graphposition);
//                crf<<">read_"<<rid<<std::endl<<bprsg.get_read_sequence(rid)<<std::endl;
//                ++chim;
//            }
        }
//        std::cout<<lpr.first<<","<<chim<<","<<lpr.second.size()<<","<<100.0*chim/lpr.second.size()<<std::endl;
    }
    sglib::OutputLog()<<pairedread_paths.size()<<" informative single read paths created!"<<std::endl;
}

void LocalHaplotypeAssembler::path_all_reads() {
    linkedread_paths.clear();
    pairedread_paths.clear();
    linkedread_paths.reserve(1000000);//TODO: do this better!!!!
    pairedread_paths.reserve(1000000);//TODO: do this better!!!!

    //now populate the linked read paths first

    CStringKMerFactory cskf(31);
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
    std::vector<std::pair<uint64_t,bool>> read1kmers,read2kmers;
    std::vector<sgNodeID_t> kmernodes;

    //BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lib], 1000000, 1000);
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],100000,1000);
    for (auto t:tagSet) {
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;

            read1kmers.clear();
            read2kmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(read1kmers, blrsg.get_read_sequence(rid));
            cskf.create_kmers_direction(read2kmers, blrsg.get_read_sequence(rid + 1));
            //first put the kmers from read 1 in there;
            add_readkmer_nodes_lha(kmernodes, read1kmers, assembly.kmer_to_graphposition, false);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
            add_readkmer_nodes_lha(kmernodes, read2kmers, assembly.kmer_to_graphposition, true);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;

            linkedread_paths.emplace_back(kmernodes);
        }
    }
    sglib::OutputLog()<<linkedread_paths.size()<<" linked-reads paths created!"<<std::endl;

    //now do the same for each paired library
    for (auto lpr:paired_reads) {
        BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lpr.first],100000,1000);
        for (auto rid : lpr.second) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;

            read1kmers.clear();
            read2kmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(read1kmers, bprsg.get_read_sequence(rid));
            cskf.create_kmers_direction(read2kmers, bprsg.get_read_sequence(rid + 1));
            //first put the kmers from read 1 in there;
            add_readkmer_nodes_lha(kmernodes, read1kmers, assembly.kmer_to_graphposition, true);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
            add_readkmer_nodes_lha(kmernodes, read2kmers, assembly.kmer_to_graphposition, false);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;

            pairedread_paths.emplace_back(kmernodes);
        }
    }
    sglib::OutputLog()<<pairedread_paths.size()<<" paired-reads paths created!"<<std::endl;

}

void LocalHaplotypeAssembler::write_problem(std::string prefix) {
    std::ofstream output_file(prefix+".bsglhap");
    uint64_t count;

    //write down backbone;
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)backbone.data(),count*sizeof(backbone[0]));

    //write down 10x tags;
    count=tagSet.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &t:tagSet) output_file.write((char *)&t,sizeof(t));

    //write down paired read ids;
    count=paired_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paired_reads){
        output_file.write((char *)&p.first,sizeof(p.first));
        count=p.second.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write((char *)p.second.data(),count*sizeof(p.second[0]));
    }

}

void LocalHaplotypeAssembler::write_full(std::string prefix) {
    std::ofstream output_file(prefix+".bsglhapf");
    uint64_t count;
    //write down backbone;
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)backbone.data(),count*sizeof(backbone[0]));
    //write  backbone node's sequences
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto n:backbone_nodes){
        count=n.sequence.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write(n.sequence.c_str(),count);
    }
    //write down 10x tags;
    count=tagSet.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &t:tagSet) output_file.write((char *)&t,sizeof(t));

    //write the condensation of the 10x workspace
    ws.linked_read_datastores[0].write_selection(output_file,tagSet);
    //write the condensation of the paired reads workspace.
    count=paired_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paired_reads){
        ws.paired_read_datastores[p.first].write_selection(output_file,p.second);
    }

    count = long_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &l:long_reads){
        BufferedSequenceGetter sequenceGetter(ws.long_read_datastores[l.first]);
        sequenceGetter.write_selection(output_file,l.second);
    }
}


void LocalHaplotypeAssembler::init_from_full_file(std::string full_file) {
    std::cout<<"Loading LocalHaplotypeAssembler problem from file: "<<full_file<<std::endl;

    std::ifstream input_file(full_file);
    uint64_t count;
    //load backbone;
    input_file.read((char *)&count,sizeof(count));
    backbone.resize(count);
    input_file.read((char *)backbone.data(),count*sizeof(backbone[0]));
    //write  backbone node's sequences
    uint64_t bcount;
    input_file.read((char *)&bcount,sizeof(bcount));
    for (auto n=0;n<bcount;++n){
        std::string seq;
        input_file.read((char *)&count,sizeof(count));
        seq.resize(count);
        input_file.read((char *)seq.data(),count);
        backbone_nodes.emplace_back(seq);
    }
    //write down 10x tags;
    input_file.read((char *)&count,sizeof(count));
    for (auto n=0;n<count;++n) {
        bsg10xTag t;
        input_file.read((char *) &t, sizeof(t));
        tagSet.insert(t);
    }

    //write the condensation of the 10x workspace
    ws.linked_read_datastores.emplace_back();
    ws.linked_read_datastores[0].load_from_stream(full_file,input_file);
    //ws.linked_read_datastores[0].write_selection(output_file,tagSet);
    //write the condensation of the paired reads workspace.
    input_file.read((char *)&count,sizeof(count));
    for (auto n=0;n<count;++n){
        ws.paired_read_datastores.emplace_back();
        ws.paired_read_datastores.back().load_from_stream(full_file,input_file);
        std::vector<uint64_t> read_ids;
        read_ids.reserve(ws.paired_read_datastores.back().size());
        for(uint64_t i=0;i<ws.paired_read_datastores.back().size();++i) read_ids.emplace_back(i);
        paired_reads.emplace_back(std::make_pair(n,read_ids));
    }


    input_file.read((char *)&count,sizeof(count));
    for (auto n=0;n<count;++n){
        ws.long_read_datastores.emplace_back();
        ws.long_read_datastores.back().load_from_stream(full_file,input_file);
        std::vector<uint64_t> read_ids;
        read_ids.reserve(ws.long_read_datastores.back().size());
        for(uint64_t i=0;i<ws.long_read_datastores.back().size();++i) read_ids.emplace_back(i);
        long_reads.emplace_back(std::make_pair(n,read_ids));
    }
}

void LocalHaplotypeAssembler::write_anchors(std::string filename) {
    std::ofstream anchf(filename);
    for (auto n = 0; n < backbone.size(); ++n) {
        anchf << ">seq" << llabs(backbone[n]) << std::endl;
        anchf << backbone_nodes[n].sequence << std::endl;

    }
}

void LocalHaplotypeAssembler::write_gfa(std::string filename) {
    assembly.write_to_gfa(filename);
}

void LocalHaplotypeAssembler::construct_patches() {
    const size_t ENDS_SIZE = 200;
    patches.clear();
    std::vector<Node > all_nodes;
    for (auto &n:assembly.nodes) {
        if (n.sequence.size()<ENDS_SIZE*2) continue;
        all_nodes.emplace_back(n);
        all_nodes.emplace_back(n);
        all_nodes.back().make_rc();
    }
    //std::vector<std::pair<std::pair<sgNodeID_t , sgNodeID_t >, std::string>> patches;
    for (auto li = 0; li < backbone.size() - 1; ++li) {
        auto n1 = backbone_nodes[li];
        auto n2 = backbone_nodes[li+1];

        if (backbone[li] < 0) n1.make_rc();
        if (backbone[li + 1] < 0) n2.make_rc();
        if (n1.sequence.size() > ENDS_SIZE)
            n1.sequence = n1.sequence.substr(n1.sequence.size() - ENDS_SIZE - 1, ENDS_SIZE);
        if (n2.sequence.size() > ENDS_SIZE) n2.sequence.resize(ENDS_SIZE);
        std::vector<std::string> matches;
        for (auto &unitig:all_nodes) {
            auto n1pos = unitig.sequence.find(n1.sequence);
            auto n2pos = unitig.sequence.find(n2.sequence);
            if (n1pos < n2pos and n2pos < unitig.sequence.size()) {
                //std::cout << lines[i][li] << " and " << lines[i][li + 1] << " found on unitig " << n
                //          << std::endl;
                matches.emplace_back(unitig.sequence.substr(n1pos, n2pos + 2 * ENDS_SIZE - n1pos));
            }
        }
        //TODO: collapse unitigs that are equivalent.
        if (matches.size() == 1) {
            patches.emplace_back(std::make_pair(backbone[li], backbone[li + 1]), matches[0]);
        }
    }
}

void LocalHaplotypeAssembler::write_patches(std::string filename) {
    std::ofstream patchf(filename);
    for (auto &p:patches) {
        patchf << ">patch_" << -p.first.first << "_" << p.first.second << std::endl;
        patchf << p.second << std::endl;
    }
}

void LocalHaplotypeAssembler::write_patched_backbone(std::string filename) {
    std::ofstream patchf(filename);
    int i=0;
    for (auto &p:patched_backbone) {
        ++i;
        patchf << ">bb_" << backbone.front() << "_" << backbone.back() << "_" << i << std::endl;
        patchf << p << std::endl;
    }
}

void LocalHaplotypeAssembler::construct_patched_backbone(bool single_scaffold, bool extend_ends,
                                                         bool extend_internals) {
    const int ENDS_SIZE = 200;
    patched_backbone.clear();
    std::vector<std::string> start_seq,end_seq,anchor_seq;
    //for each anchor, in proper orientation, construct the start and end sequences.
    for (auto bn=0;bn<backbone_nodes.size();++bn){
        auto n = backbone_nodes[bn];
        if (backbone[bn] < 0) n.make_rc();
        anchor_seq.emplace_back(n.sequence);
        start_seq.emplace_back(n.sequence.substr(0,ENDS_SIZE));
        end_seq.emplace_back(n.sequence.substr(n.sequence.size() - ENDS_SIZE , ENDS_SIZE));
    }

    std::vector<std::string> all_nodes;
    //create a collection of all sequences and their reverses, for every node in the assembly >ENDS_SIZE
    for (auto &n:assembly.nodes) {
        if (n.sequence.size()<ENDS_SIZE*2) continue;
        all_nodes.emplace_back(n.sequence);
        auto n2=n;
        n2.make_rc();
        all_nodes.emplace_back(n2.sequence);
    }

    //now do the processing

    bool last_linked=0;
    std::string seq="";
    //std::cout<<std::endl;
    for (auto bn=0;bn<backbone_nodes.size();++bn) {

        {
            //find the start of this anchor
            std::string start_matching_seq = "";
            size_t start_p;
            for (auto &s:all_nodes) {
                auto p = s.find(start_seq[bn]);
                if (p < s.size()) {
                    if (!start_matching_seq.empty()) {
                        patched_backbone.clear();
                        return;
                    }
                    start_matching_seq = s;
                    start_p = p;
                    //std::cout<<"Start of anchor #"<<bn<<" ("<<backbone[bn]<<") found on pos "<<start_p<<" of unititg"<<std::endl;
                }
            }
            //if (start_matching_seq.empty()) std::cout<<"Start of anchor #"<<bn<<" ("<<backbone[bn]<<") NOT FOUND"<<std::endl;
            //pre-fill with previous sequence if asked for (i.e. pre-extension)
            if (!last_linked and !start_matching_seq.empty()) {
                if ((bn == 0 and extend_ends) or (bn != 0 and extend_internals)) {
                    //std::cout<<"pre-extending..."<<std::endl;
                    seq += start_matching_seq.substr(0, start_p);
                }
            }
        }

        //add node's sequence to seq
        {
            //std::cout<<"anchor #"<<bn<<" ("<<backbone[bn]<<") added to sequence"<<std::endl;
            seq+=anchor_seq[bn];
        }

        {
            //look for bn's end_sequence
            std::string end_matching_seq = "";
            size_t end_p;
            for (auto &s:all_nodes) {
                auto p = s.find(end_seq[bn]);
                if (p < s.size()) {
                    if (!end_matching_seq.empty()) {
                        patched_backbone.clear();
                        return;
                    }
                    end_matching_seq = s;
                    end_p = p+ENDS_SIZE;//THIS IS SO IF POINTS AFTER THE MATCH
                    //std::cout<<"End of anchor #"<<bn<<" ("<<backbone[bn]<<") found on pos "<<end_p<<" of unititg"<<std::endl;
                }
            }
            //if (end_matching_seq.empty()) std::cout<<"End of anchor #"<<bn<<" ("<<backbone[bn]<<") NOT FOUND"<<std::endl;
            last_linked=false;
            if (bn < backbone_nodes.size() - 1) {
                //check if next's node start_seq is on same contig, higher position
                auto np=end_matching_seq.find(start_seq[bn+1]);
                if (np<end_matching_seq.size()){
                    //std::cout<<"Start of next anchor #"<<bn+1<<" ("<<backbone[bn+1]<<") found on pos "<<np<<" of this anchor's end unititg, last_linked=true!"<<std::endl;
                    last_linked=true;
                    if (np>end_p) {
                        //true -> add "patch" sequence to seq; last_linked=true
                        seq += end_matching_seq.substr(end_p, np - end_p);
                        //std::cout<<"patch added!!!!"<<std::endl;
                    } else {
                        seq=seq.substr(0,seq.size()-(np-end_p));//TODO:test this one!!!!
                        //std::cout<<"overlap removed!!!!"<<std::endl;
                    }
                }
                //else std::cout<<"Start of next anchor #"<<bn+1<<" ("<<backbone[bn+1]<<") NOT FOUND on this anchor's end unititg"<<std::endl;
            }
            if (!last_linked) {
                if ((bn == backbone_nodes.size() - 1 and extend_ends) or
                    (bn != backbone_nodes.size() - 1 and extend_internals)) {
                    if (!end_matching_seq.empty()) {
                        //std::cout<<"Not last_linked, extending FW"<<std::endl;
                        seq+=end_matching_seq.substr(end_p);
                    }
                }
                if (!single_scaffold or bn == backbone_nodes.size() - 1) {
                    //std::cout<<"End of sequence, adding to collection"<<std::endl;
                    patched_backbone.emplace_back(seq);
                    seq = "";
                } else {
                    //std::cout<<"can't link fw but single scaffold mode, adding Ns"<<std::endl;
                    seq += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
                }
            }
        }
        //std::cout<<std::endl;

    }
}

std::vector<std::pair<std::string, std::string>> LocalHaplotypeAssembler::compute_metrics() {
    std::vector<std::pair<std::string, std::string>> metrics;
    construct_patches();
    metrics.emplace_back("Patches",std::to_string(patches.size()));
    metrics.emplace_back("Patches%",std::to_string(100.0*patches.size()/(backbone_nodes.size()-1)));
    int found_bb=0;
    for (auto bbn:backbone_nodes){
        auto bbnr=bbn;
        bbnr.make_rc();

        for (auto &n:assembly.nodes) {
            if (n.sequence.size() < bbn.sequence.size()) continue;
            if (n.sequence.find(bbn.sequence) < n.sequence.size() or n.sequence.find(bbnr.sequence) < n.sequence.size()) {
                ++found_bb;
                break;
            }
        }
    }
    metrics.emplace_back("InitialAnchors",std::to_string(backbone_nodes.size()));
    metrics.emplace_back("FullAnchors",std::to_string(found_bb));
    metrics.emplace_back("FullAnchors%",std::to_string(100.0*found_bb/backbone_nodes.size()));
    metrics.emplace_back("Unitigs",std::to_string(assembly.count_active_nodes()));
    //metrics.emplace_back("N50",std::to_string(assembly.computeNXX(50)));
    //metrics.emplace_back("Full_anchors",std::to_string(patches.size()));

    return metrics;
}

void LocalHaplotypeAssembler::problem_analysis(std::string prefix) {
    //std::ofstream logfile(prefix+"_analysis.log");

    write_anchors(prefix+"_anchors.fasta");

    //====== 1) Dump read files ===
    std::ofstream lr_r1file(prefix+"_lr_R1.fasta");
    std::ofstream lr_r2file(prefix+"_lr_R2.fasta");

    for (auto t:tagSet) {
//        auto chim=0;
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            if (rid%2==0) continue;
            lr_r1file<<">"<<t<<"_"<<rid<<""<<std::endl<<ws.linked_read_datastores[0].get_read_sequence(rid)<<std::endl;
            lr_r2file<<">"<<t<<"_"<<rid<<""<<std::endl<<ws.linked_read_datastores[0].get_read_sequence(rid+1)<<std::endl;
        }
    }
    auto lmplib=0;
    for (auto lpr:paired_reads) {
        ++lmplib;
        std::ofstream lmp_r1file(prefix+"_lmp"+std::to_string(lmplib)+"_R1.fasta");
        std::ofstream lmp_r2file(prefix+"_lmp"+std::to_string(lmplib)+"_R2.fasta");
        BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lpr.first],100000,1000);
        for (auto rid : lpr.second) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;
            lmp_r1file<<">"<<rid<<""<<std::endl<<bprsg.get_read_sequence(rid)<<std::endl;
            lmp_r2file<<">"<<rid<<""<<std::endl<<bprsg.get_read_sequence(rid+1)<<std::endl;
        }
    }

    //====== 2) Trivial 31-mer graph
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
    auto ltkmers128 = ws.linked_read_datastores[0].get_tags_kmers(31, 5, tagSet, blrsg);
    GraphMaker gm(assembly);
    gm.new_graph_from_kmerset_trivial(ltkmers128, 31);
    write_gfa(prefix+"_31-mer_DBG.gfa");

}