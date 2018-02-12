//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include "LinkedReadsDatastore.hpp"

void LinkedReadsDatastore::build_from_fastq(std::string read1_filename, std::string read2_filename,
                                            LinkedReadsFormat format) {
    if (format==LinkedReadsFormat::raw) {
        std::cout<<"Raw 10x format not supported yet, datastore not populated!!!"<<std::endl;
        return;
    }
    //std::cout<<"Memory used by every read's entry:"<< sizeof(LinkedRead)<<std::endl;
    //read each read, put it on the index and on the appropriate tag
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Creating Datastore Index from "<<read1_filename<<" | "<<read2_filename<<std::endl;
    filename1=read1_filename;
    filename2=read2_filename;
    fd1=fopen(read1_filename.c_str(),"r");
    fd2=fopen(read2_filename.c_str(),"r");
    char readbuffer[1000];
    read_offset.resize(1);//leave read 0 empty.
    uint16_t r1offset,r2offset;
    uint64_t g1offset=0,g2offset=0;
    uint64_t readCountInFile=0;
    uint64_t tagged_reads=0;
    while (!feof(fd1) and !feof(fd2)){
        //LinkedRead r1,r2;
        bsg10xTag newtag=0;
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        //Tag to number from r1's name
        for(auto i=1;i<17;++i){
            newtag<<=2;
            if (readbuffer[i]=='C') newtag+=1;
            else if (readbuffer[i]=='G') newtag+=2;
            else if (readbuffer[i]=='T') newtag+=3;
            else if (readbuffer[i]!='A') {
                newtag=0;
                break;
            }
        }
        r1offset=ftell(fd1);
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        if (NULL == fgets(readbuffer, 999, fd1)) continue;

        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        r2offset=ftell(fd2);
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        if (NULL == fgets(readbuffer, 999, fd2)) continue;

        if (0==readCountInFile%group_size) {
            g1offset=r1offset;
            group_offset1.push_back(r1offset);
            g2offset=r2offset;
            group_offset2.push_back(r2offset);
        }
        read_offset.push_back(r1offset-g1offset);
        read_offset.push_back(r2offset-g2offset);
        ++readCountInFile;
        read_tag.push_back(newtag);
        if (0!=newtag) tagged_reads+=2;
    }
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Datastore with "<<read_offset.size()-1<<" reads, "<<tagged_reads<<" reads with tags"<<std::endl; //and "<<reads_in_tag.size()<<"tags"<<std::endl;

}

std::string LinkedReadsDatastore::get_read_sequence(size_t readID, FILE * file1, FILE * file2) {
    char buffer[1000];
    if (0==readID%2){
        auto pos_in_file=readID/2-1;
//        std::cout<<"Read from file2, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset2[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        fseek(file2,group_offset2[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        fread(buffer,1000,1,file2);
    } else {
        auto pos_in_file=readID/2;
//        std::cout<<"Read from file1, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset1[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        fseek(file1,group_offset1[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        fread(buffer,1000,1,file1);
    }
    for (auto &c:buffer) if (c=='\n') {c='\0'; break;}
    return std::string(buffer);
}

bsg10xTag LinkedReadsDatastore::get_read_tag(size_t readID) {
    return read_tag[(readID-1)/2];
}