#include <string>
#include <unordered_map>

//#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#include "cxxopts.hpp"
#include "sglib/filesystem/check_or_create_directory.h"
#include "sglib/readers/FileReader.h"
#include "sglib/factories/KMerIDXFactory.h"
#include "sglib/factories/ContigLink.h"
#include "sglib/SMR.h"

void log_kmer_density(const std::string &output_prefix, const std::vector<KmerIDX> &unique_kmers,
                      const std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> &sequences_fasta) {
    std::ofstream contig_kmer_density(output_prefix + "contig_kmer_density.log");
    for (const auto &seq: sequences_fasta){
        contig_kmer_density << seq.first << seq.second.second * 1.0f / unique_kmers.size();
    }
}

int main(int argc, char **argv) {

    const uint64_t GB(1024*1024*1024);
    unsigned int mem_limit(10);

    std::string asm_filename;
    std::string fastq_filename;
    uint16_t min_count(0);
    uint32_t max_count(100000);
    uint16_t max_coverage(1);
    uint32_t min_read_length(1000), min_contig_length(1000);
    uint32_t min_kmers_to_call_match(10);
    uint32_t min_seen_contig_to_write_output(1);

    std::string set_filelist;
    std::string output_prefix;

    uint8_t m(1),n(1),k(31);
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::string outdefault(
            std::to_string(tm.tm_year + 1900) + '-' + std::to_string(tm.tm_mon) + '-' + std::to_string(tm.tm_mday) +
            '_' + std::to_string(tm.tm_hour) + std::to_string(tm.tm_min));

    cxxopts::Options options("seq-sorter", "Sequence linking tool using long/linked reads.");
//@formatter:off
    options.add_options()("help", "Print help", cxxopts::value<std::string>()->implicit_value(""), "")
            ("a,assembly", "Sequence reference to link",cxxopts::value<std::string>(asm_filename), "FASTA - Sequence file")
            ("r,long_reads", "Reads to generate sequence-to-sequence links",cxxopts::value<std::string>(fastq_filename), "FASTQ - Reads")
            ("min_count","Minimum count to consider a link",cxxopts::value<uint16_t>(min_count)->default_value("10"),"uint")
            ("max_count","Maximum count to consider a link",cxxopts::value<uint32_t>(max_count)->default_value("100000"), "uint")
            ("max_coverage", "max coverage for a kmer to be considered",cxxopts::value<uint16_t>(max_coverage)->default_value("1"), "uint")
            ("min_read_length","minimum contig length",cxxopts::value<uint32_t>(min_read_length)->default_value("1000"), "uint")
            ("min_contig_length", "minimum read length",cxxopts::value<uint32_t>(min_contig_length)->default_value("1000"), "uint")
            ("min_kmers_to_call_match","minimum number of kmers to call a Read->Contig match",cxxopts::value<uint32_t>(min_kmers_to_call_match)->default_value("10"), "uint")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix)->default_value(outdefault),"prefix_dir");

    options.add_options("Output options")
            ("min_seen_contig_to_write_output","minimum number of seen contigs to report read on output",cxxopts::value<uint32_t>(min_seen_contig_to_write_output)->default_value("1"), "uint");

    options.add_options("Skip-mer shape (m every n, total k)")
            ("m,used_bases", "m", cxxopts::value<uint8_t>(m)->default_value("1"), "uint")
            ("n,skipped_bases", "n", cxxopts::value<uint8_t>(n)->default_value("1"), "uint")
            ("k,total_bases", "k",cxxopts::value<uint8_t>(k)->default_value("31"), "uint");

    options.add_options("Performance")
            ("mem_limit", "Memory limit in GB",cxxopts::value<unsigned int>(mem_limit)->default_value("10"), "uint");
//@formatter:on

    try
    {
        auto result (options.parse(argc, argv));

        if (0 != result.count("help")) {
            std::cout << options.help({"", "Performance"}) << std::endl;
            exit(0);
        }
        if (asm_filename.empty()) {
            throw cxxopts::OptionException("The assembly file parameter wasn't specified");
        } else if (!std::ifstream(asm_filename)) {
            cxxopts::OptionException("FASTA file doesn't exist");
        }
        if (fastq_filename.empty()) {
            cxxopts::missing_argument_exception("The read file parameter wasn't specified");
        } else if (!std::ifstream(fastq_filename)) {
            cxxopts::OptionException("FASTQ file doesn't exist");
        }
    } catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }

    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    std::string contig_file(asm_filename);
    std::string pacbio_reads_file(fastq_filename);

    uint64_t maxmem(mem_limit * GB);
    /*
     * Reads the reference/contigs file and generates the set of kmers and which entry number on the fasta they belong to
     * along with their orientation in the form of a int32_t with + for fwd and - for rev, the output is filtered
     * by min < coverage <= max
     */
    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            FastaReader<FastaRecord>,
            FastaRecord, FastxReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1}, {k}, {maxmem, 0, max_coverage,
                                                                              output_prefix});

    std::vector<KmerIDX> unique_kmers;

    // Get the unique_kmers from the file
    unique_kmers = kmerIDX_SMR.read_from_file(contig_file);

    std::vector<uint64_t > uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of " << int(k) << "-kmers that appear more than 0 < kmer_coverage <= " << max_coverage
              << " in assembly " << uniqKmer_statistics[1] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;

    /*
     * Store the lengths of the sequences longer than min_contig_length
     * and keep a set of the shorter ones to filter the unique_kmers map with
     */
    FastaReader<FastaRecord> reader({0}, asm_filename);
    std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> sequences_fasta;
    std::unordered_set<int> invalid_sequences;
    FastaRecord rec;
    std::ofstream out_links(output_prefix + "graph.gfa");
    out_links << "H\tVN:Z:1.0" << std::endl;
    while (reader.next_record(rec)) {
        if (rec.seq.size() > min_contig_length) {
            out_links << "S\t" << rec.id << "\t" << rec.seq << std::endl;
            sequences_fasta[rec.id] = std::make_pair(rec.seq.size(),0);
        } else {
            invalid_sequences.insert(rec.id);
        }
    }

    /*
     * Delete all the kmers in the map that come from "invalid" unitigs
     */
    std::remove_if(unique_kmers.begin(), unique_kmers.end(),
                   [&invalid_sequences] (const KmerIDX &kidx)
                   {
                       // Remove if I found the contigID on the invalid_sequences
                       return (invalid_sequences.find(kidx.contigID) != invalid_sequences.cend());
                   }
    );
    std::cout << "Number of " << int(k) << "-mers in contigs larger than " << min_contig_length << " " << unique_kmers.size() << std::endl;
    // Cleanup invalid_sequences
    invalid_sequences = std::unordered_set<int>();

    for (const auto &kdx: unique_kmers) {
        sequences_fasta[std::abs(kdx.contigID)].second++;
    }

    log_kmer_density(output_prefix, unique_kmers, sequences_fasta);

    /*
     * Instantiate a ContigLinkFactory that takes the unique_kmers as a parameter, a fastq file with long/linked reads
     * and generates a list of links between contigs.
     */
    SMR<ContigLink,
            ContigLinkFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, FilterSetParams> contigLink_SMR({min_read_length},   // ReaderParams
                                                                            {output_prefix, k,   // FilterParams
                                                                             unique_kmers, sequences_fasta,
                                                                             min_kmers_to_call_match,
                                                                             min_seen_contig_to_write_output},
                                                                            {maxmem, min_count, max_count,
                                                                             output_prefix});

    std::vector<ContigLink> links(contigLink_SMR.read_from_file(pacbio_reads_file));
    std::vector<uint64_t > link_statistics(contigLink_SMR.summaryStatistics());

    std::cout << "Number of links generated " << link_statistics[0] << std::endl;
    std::cout << "Number of links that appear more than " << min_count << " times but less than " << max_count << " " << link_statistics[1] << std::endl;
    std::cout << "Number of reads " << link_statistics[2] << std::endl;
    std::cout << "Number of reads larger than " << min_read_length << " " << link_statistics[3] << std::endl;


    // Output the count sorted links to a file in gfa friendly format.
    std::sort(links.begin(), links.end(), ContigLink::byCount());
    for(auto &l:links) {
        if (l.count >= min_count) {
            out_links << "L\t" << abs(l.contig1) << "\t" << ((l.contig1 > 0) ? '+' : '-') << "\t" << abs(l.contig2)
                      << "\t" << ((l.contig2 > 0) ? '+' : '-') << "\t0M\t" << "RC:i:" << l.count << std::endl;
        }
    }
    return 0;
}