//
// Created by Bernardo Clavijo (EI) on 19/10/2020.
//

#include "CountFilter.hpp"

CountFilter::CountFilter(std::string _kcname, std::string _filter_count_name, int _filter_count_max, const std::vector<std::string> _value_count_names,const std::vector<int> _value_count_mins):
        kcname(_kcname),
        filter_count_name(_filter_count_name),
        filter_count_max(_filter_count_max),
        value_count_names(_value_count_names),
        value_count_mins(_value_count_mins){};

std::string CountFilter::get_pattern(NodeView nv) {
    if (patterns.count(llabs(nv.node_id()))==0){
        auto kcf=nv.kmer_coverage(kcname,filter_count_name);
        std::vector<int> fkpos;
        for (auto i=0;i<kcf.size();++i) {
            if (kcf[i]<=filter_count_max) fkpos.push_back(i);
        }
        if (fkpos.empty()) patterns[llabs(nv.node_id())]="";
        else {
            std::string p="";
            for (auto cni=0;cni<value_count_names.size();++cni){
                auto vcf=nv.kmer_coverage(kcname,value_count_names[cni]);
                for (auto i=0;i<vcf.size();++i) {
                    if (vcf[i]<value_count_mins[cni]) {
                        p.append("0");
                        break;
                    }
                }
                if (p.size()<=cni) p.append("1");
            }
            patterns[llabs(nv.node_id())]=p;
        }
    }
    return patterns[llabs(nv.node_id())];
}