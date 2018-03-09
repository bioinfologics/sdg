//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include "LinkedReadsDatastore.hpp"

void LinkedReadsDatastore::build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, LinkedReadsFormat format, int _rs, size_t chunksize) {

    //std::cout<<"Memory used by every read's entry:"<< sizeof(LinkedRead)<<std::endl;
    //read each read, put it on the index and on the appropriate tag
    readsize=_rs;
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Creating Datastore Index from "<<read1_filename<<" | "<<read2_filename<<std::endl;
    auto fd1=fopen(read1_filename.c_str(),"r");
    auto fd2=fopen(read2_filename.c_str(),"r");
    char readbuffer[1000];
    uint64_t r1offset,r2offset;
    uint64_t tagged_reads=0;
    //first, build an index of tags and offsets
    sglib::OutputLog()<<"Building tag sorted chunks of "<<chunksize<<" pairs"<<std::endl;
    std::vector<readData> readdatav;
    readdatav.reserve(chunksize);
    std::vector<std::ifstream> chunkfiles;
    std::vector<readData> next_in_chunk;
    readData currrent_read;
    //First, create the chunk files
    uint64_t pairs=0;
    while (!feof(fd1) and !feof(fd2)) {
        bsg10xTag newtag = 0;
        if (format==LinkedReadsFormat::UCDavis) {
            //LinkedRead r1,r2;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            //Tag to number from r1's name
            for (auto i = 1; i < 17; ++i) {
                newtag <<= 2;
                if (readbuffer[i] == 'C') newtag += 1;
                else if (readbuffer[i] == 'G') newtag += 2;
                else if (readbuffer[i] == 'T') newtag += 3;
                else if (readbuffer[i] != 'A') {
                    newtag = 0;
                    break;
                }
            }
            currrent_read.tag=newtag;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            currrent_read.seq1=std::string(readbuffer);
            if (currrent_read.seq1.back()=='\n') currrent_read.seq1.resize(currrent_read.seq1.size()-1);
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;

            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            currrent_read.seq2=std::string(readbuffer);
            if (currrent_read.seq2.back()=='\n') currrent_read.seq2.resize(currrent_read.seq2.size()-1);
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
        }
        else if (format==LinkedReadsFormat::seq){
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            //Tag to number from r1's name
            for (auto i = 0; i < 16; ++i) {
                newtag <<= 2;
                if (readbuffer[i] == 'C') newtag += 1;
                else if (readbuffer[i] == 'G') newtag += 2;
                else if (readbuffer[i] == 'T') newtag += 3;
                else if (readbuffer[i] != 'A') {
                    newtag = 0;
                    break;
                }
            }
            currrent_read.tag=newtag;
            currrent_read.seq1=std::string(readbuffer+16);
            if (currrent_read.seq1.back()=='\n') currrent_read.seq1.resize(currrent_read.seq1.size()-1);
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            currrent_read.seq2=std::string(readbuffer);
            if (currrent_read.seq2.back()=='\n') currrent_read.seq2.resize(currrent_read.seq2.size()-1);
        }
        if (0 != newtag) tagged_reads += 2;
        ++pairs;
        readdatav.push_back(currrent_read);
        if (readdatav.size()==chunksize){
            //sort
            std::sort(readdatav.begin(),readdatav.end());
            //dump
            std::ofstream ofile("sorted_chunk_"+std::to_string(chunkfiles.size())+".data");
            sglib::OutputLog()<<readdatav.size()<<" pairs dumping on chunk "<<chunkfiles.size()<<std::endl;
            //add file to vector of files
            char buffer[2*readsize+2];
            for (auto &r:readdatav){
                ofile.write((const char * ) &r.tag,sizeof(r.tag));
                bzero(buffer,2*readsize+2);
                memcpy(buffer,r.seq1.data(),(r.seq1.size()>readsize ? readsize : r.seq1.size()));
                memcpy(buffer+readsize+1,r.seq2.data(),(r.seq2.size()>readsize ? readsize : r.seq2.size()));
                ofile.write(buffer,2*readsize+2);
            }
            ofile.close();
            chunkfiles.emplace_back("sorted_chunk_"+std::to_string(chunkfiles.size())+".data");
            readdatav.clear();
            sglib::OutputLog()<<"dumped!"<<std::endl;
        }
    }
    if (readdatav.size()>0) {
        //sort
        std::sort(readdatav.begin(), readdatav.end());
        //dump
        std::ofstream ofile("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data");
        sglib::OutputLog() << readdatav.size() << " pairs dumping on chunk " << chunkfiles.size() << std::endl;
        //add file to vector of files
        char buffer[2 * readsize + 2];
        for (auto &r:readdatav) {
            ofile.write((const char *) &r.tag, sizeof(r.tag));
            bzero(buffer, 2 * readsize + 2);
            memcpy(buffer, r.seq1.data(), (r.seq1.size() > readsize ? readsize : r.seq1.size()));
            memcpy(buffer + readsize + 1, r.seq2.data(), (r.seq2.size() > readsize ? readsize : r.seq2.size()));
            ofile.write(buffer, 2 * readsize + 2);
        }
        ofile.close();
        chunkfiles.emplace_back("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data");
        readdatav.clear();
        sglib::OutputLog() << "dumped!" << std::endl;
    }
    sglib::OutputLog() << "performing merge from disk" << std::endl;
    //TODO: save space first for the tag index!!!
    std::ofstream output(output_filename.c_str());
    output.write((const char *) &readsize,sizeof(readsize));
    read_tag.resize(pairs+1);
    sglib::OutputLog() << "leaving space for " <<pairs<<" read_tag entries"<< std::endl;
    uint64_t rts=read_tag.size();
    output.write((const char *) &rts, sizeof(rts));
    output.write((const char *) read_tag.data(),sizeof(bsg10xTag)*read_tag.size());
    //multi_way merge of the chunks
    int openfiles=chunkfiles.size();
    bsg10xTag next_tags[chunkfiles.size()];
    char buffer[2 * readsize + 2];
    //read a bit from each file
    for (auto i=0;i<chunkfiles.size();++i) chunkfiles[i].read((char *)&next_tags[i],sizeof(bsg10xTag));
    read_tag.resize(1);
    while(openfiles){
        bsg10xTag mintag=UINT32_MAX;
        //find the minimum tag
        for (auto &t:next_tags) if (t<mintag) mintag=t;
        //copy from each file till the next tag is > than minimum
        for (auto i=0;i<chunkfiles.size();++i){
            while (!chunkfiles[i].eof() and next_tags[i]==mintag){
                //read buffer from file and write to final file
                chunkfiles[i].read(buffer,2*readsize+2);
                output.write(buffer,2*readsize+2);
                read_tag.push_back(mintag);
                //read next tag... eof? tag=UINT32_MAX
                chunkfiles[i].read((char *)&next_tags[i],sizeof(bsg10xTag));
                if (chunkfiles[i].eof()) {
                    --openfiles;
                    next_tags[i]=UINT32_MAX;
                    sglib::OutputLog() << "chunk "<<i<<" finished"<<std::endl;
                }
            }
        }
    }
    //go back to the beginning of the file and write the read_tag part again
    output.seekp(sizeof(readsize));
    sglib::OutputLog() << "writing down " <<pairs<<" read_tag entries"<< std::endl;
    rts=read_tag.size();
    output.write((const char *) &rts, sizeof(rts));
    output.write((const char *) read_tag.data(),sizeof(bsg10xTag)*read_tag.size());
    output.close();
    //delete all temporary chunk files
    for (auto &c:chunkfiles) c.close();
    for (auto i=0;i<chunkfiles.size();++i) ::unlink(("sorted_chunk_"+std::to_string(i)+".data").c_str());
    //DONE!
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Datastore with "<<(read_tag.size()-1)*2<<" reads, "<<tagged_reads<<" reads with tags"<<std::endl; //and "<<reads_in_tag.size()<<"tags"<<std::endl;
    filename=output_filename;
    fd=fopen(filename.c_str(),"r");
}

void LinkedReadsDatastore::read(std::ifstream &input_file) {
    //read filename
    uint64_t s;
    input_file.read((char *) &s, sizeof(s));
    filename.resize(s);
    input_file.read((char *) filename.data(), filename.size());
    load_index(filename);
}

void LinkedReadsDatastore::load_index(std::string _filename){
    uint64_t s;
    filename=_filename;
    fd=fopen(filename.c_str(),"r");
    fread( &readsize,sizeof(readsize),1,fd);
    fread(&s,sizeof(s),1,fd); read_tag.resize(s);
    fread(read_tag.data(),sizeof(read_tag[0]),read_tag.size(),fd);
    readpos_offset=ftell(fd);
}

void LinkedReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    uint64_t s=filename.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)filename.data(),filename.size());
}

std::string LinkedReadsDatastore::get_read_sequence(size_t readID) {
    char buffer[readsize+1];
    size_t read_offset_in_file=readpos_offset+(readsize+1)*(readID-1);
    fseek(fd,read_offset_in_file,SEEK_SET);
    fread(buffer,readsize+1,1,fd);
    return std::string(buffer);
}

std::vector<uint64_t> LinkedReadsDatastore::get_tag_reads(bsg10xTag tag) {
    std::vector<uint64_t> rids;
    rids.reserve(10000);
    for (auto n=std::lower_bound(read_tag.begin(),read_tag.end(),tag)-read_tag.begin();read_tag[n]==tag;++n) {
        rids.emplace_back( n * 2 + 1);
        rids.emplace_back( n * 2 + 2);
    }
    return rids;
}

bsg10xTag LinkedReadsDatastore::get_read_tag(size_t readID) {
    return read_tag[(readID-1)/2];
}


std::unordered_set<uint64_t> LinkedReadsDatastore::get_tags_kmers(int k, int min_tag_cov, std::unordered_set<bsg10xTag> tags, BufferedLRSequenceGetter & blrsg) {
    class StreamKmerFactory : public  KMerFactory {
    public:
        explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
        inline void produce_all_kmers(const char * seq, std::vector<uint64_t> &mers){
            // TODO: Adjust for when K is larger than what fits in uint64_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory skf(k);

    //reserve space by counting reads first, save only the integer, do not merge just count and insert in the set
    std::vector<uint64_t> all_kmers;
    std::vector<uint64_t> read_ids;
    for (auto t:tags) {
        auto reads=get_tag_reads(t);
        read_ids.insert(read_ids.end(),reads.begin(),reads.end());
    }
    all_kmers.reserve(read_ids.size()*(readsize-k+1));
    for (auto rid:read_ids){
            skf.produce_all_kmers(blrsg.get_read_sequence(rid),all_kmers);
     }
    std::sort(all_kmers.begin(),all_kmers.end());
    std::unordered_set<uint64_t> kset;
    auto ri=all_kmers.begin();
    auto nri=all_kmers.begin();
    while (ri<all_kmers.end()){
        while (nri<all_kmers.end() and *nri==*ri) ++nri;
        if (nri-ri>=min_tag_cov) kset.insert(*ri);
        ri=nri;
    }
    return std::move(kset);
}

const char* BufferedLRSequenceGetter::get_read_sequence(uint64_t readID) {
        size_t read_offset_in_file=datastore.readpos_offset+(datastore.readsize+1)*(readID-1);
        if (read_offset_in_file<buffer_offset or read_offset_in_file+chunk_size>buffer_offset+bufsize) {
            buffer_offset=read_offset_in_file;
            lseek(fd,read_offset_in_file,SEEK_SET);
            read(fd,buffer,bufsize);
        }
        return buffer+(read_offset_in_file-buffer_offset);
}
