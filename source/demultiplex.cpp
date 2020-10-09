#include <string>
#include <experimental/filesystem>
#include <iostream>
#include <seqan3/std/concepts>

#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/configuration/all.hpp> 
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>

#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/alphabet/nucleotide/concept.hpp>

std::vector<seqan3::dna5_vector> readinSpikeFasta(std::filesystem::path filename)
{
  seqan3::sequence_file_input inFasta{filename};
  std::vector<seqan3::dna5_vector> sequences;

  
     for (auto &[seq, id, qual] : inFasta)
    {
        seqan3::debug_stream << "ID:  " << id << '\n';
        seqan3::debug_stream << "SEQ: " << seq << '\n';
        seqan3::debug_stream << "EMPTY QUAL." << qual << '\n'; // qual is empty for FastA files
	sequences.push_back(seq);
    }
     //makes a copy so may want to do by reference for large file
     return sequences;
}
/*
void processFastqs(seqan3::sequence_file_input fin1, seqan3::sequence_file_input fin2)
{
  for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) // && is important!
    {                                                           // because seqan3::views::zip returns temporaries
    if (seqan3::get<seqan3::field::id>(rec1) != seqan3::get<seqan3::field::id>(rec2))
        throw std::runtime_error("Oh oh your pairs don't match.");
        seqan3::debug_stream << "ID: " << seqan3::get<seqan3::field::id>(rec1) << '\n';
    }

}
*/

template<typename T, std::integral ScoreType>
bool maxScore(T & readResults, ScoreType scoreVal){
  auto tempScore = 0;
  for(auto & res : readResults){
   //Todo
    tempScore = res.score();
    //for edit distance, a positive value is invalid
   if((tempScore >= scoreVal) && (tempScore < 1)){ 
      //          seqan3::debug_stream << res << std::endl;
	  return true;
	    }
	    }
 return false;
}

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

int main(int argc, char ** argv)
{
  seqan3::argument_parser myparser{"Spike_Processor", argc, argv}; // initialize
  
    std::filesystem::path spikeFastaName{};
    std::filesystem::path fastqOneName;
    std::filesystem::path fastqTwoName;
    
    //myparser.add_option(spikeFastaName, 'f', "spike","Fasta with Spike-in sequences",
    //                    seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fa","fasta"}});
    //only validates for .gz, could change
    myparser.add_option(fastqOneName, 'p', "fastq-1","Fastq One",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});
    myparser.add_option(fastqTwoName, 'q', "fastq-2","Fastq Two",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});
    
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "Spike-in processor - [PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    //read-in spike fasta

    //auto spikeSequence = "TTGACTCACCCATCAAGATACACAACCGCT"_dna5;
    //auto spikeRC = "AGCGGTTGTGTATCTTGATGGGTGAGTCAA"_dna5;


    //shorter spike sequence
    auto spikeSequence = "CTAATTTAAACTATTCTCTGTTC TTTCATGGGGAAGCAGA TTTGGGTACCACCCAAGTATTGACTCACCCATCAAGATACACAAC"_dna5;
    auto spikeRC = "GTTGTGTATCTTGATGGGTGAGTCAATACTTGGGTGGTACCCAAATCTGCTTCCCCATGAAAGAACAGAGAATAGTTTAAATTAG"_dna5;
    //paired fastqs
    std::cout << fastqOneName << std::endl;
    seqan3::sequence_file_input fin1{fastqOneName};
    seqan3::sequence_file_input fin2{fastqTwoName};
    //processFastqs(fastqOneName, fastqTwoName);

    std::filesystem::path outName1 = fastqOneName;
    outName1.replace_extension().replace_extension(".spike-removed.fq.gz");
    seqan3::sequence_file_output outFastq1{outName1};

    std::filesystem::path outName2 = fastqTwoName;
    outName2.replace_extension().replace_extension(".spike-removed.fq.gz");
    seqan3::sequence_file_output outFastq2{outName2};

    //out spike files
    std::filesystem::path outSpikeName1 = fastqOneName;
    outSpikeName1.replace_extension().replace_extension(".spike.fq.gz");
    seqan3::sequence_file_output outSpikeFastq1{outSpikeName1};

    std::filesystem::path outSpikeName2 = fastqTwoName;
    outSpikeName2.replace_extension().replace_extension(".spike.fq.gz");
    seqan3::sequence_file_output outSpikeFastq2{outSpikeName2};

   
    
    /*auto const config = seqan3::align_cfg::method_global |
      seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{0}}} |
      seqan3::align_cfg::result{seqan3::with_alignment} |
      seqan3::align_cfg::aligned_ends{seqan3::free_ends_first};*/

    auto const config = seqan3::align_cfg::edit |
      seqan3::align_cfg::aligned_ends{seqan3::free_ends_first};
      //seqan3::align_cfg::max_error{1u};
	  
    auto count = 0;
    auto spikeCount = 0;
    auto read1SpikeCount = 0;
    auto read2SpikeCount = 0;
    auto read1RCSpikeCount = 0;
    auto read2RCSpikeCount = 0;

    //need to record if a spike was found
    auto spikeBool = false;
    //now negative because edit distance
    auto minScore = -1;
    //iterate throught the fastq files
    for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) // && is important!
	{                                                           // because seqan3::views::zip returns temporaries
	  ++count;
	  if(count % 100000 == 0){
	    seqan3::debug_stream << count / 1000000.0 << " Million Reads processed" << std::endl;
	  }
	  

	  
          //alignments of of spike and RC to the reads and check if at least one meets the minimum score (hard coded to 29). Increment count
	  auto r1forward = seqan3::align_pairwise(std::tie(seqan3::get<seqan3::field::seq>(rec1), spikeSequence), config);
	  if(maxScore(r1forward, minScore)){
	    //seqan3::debug_stream << "Read 1 Spike" << std::endl;
	    read1SpikeCount++;
	    spikeCount++;
	    outSpikeFastq1.push_back(rec1);
	    outSpikeFastq2.push_back(rec2);
	    continue;
	  }

	  auto r2forward = seqan3::align_pairwise(std::tie(seqan3::get<seqan3::field::seq>(rec2), spikeSequence), config);
	  if(maxScore(r2forward, minScore)){
	    //	    seqan3::debug_stream << "Read 2 Spike" << std::endl;
	    read2SpikeCount++;
	    spikeCount++;

	    outSpikeFastq1.push_back(rec1);
	    outSpikeFastq2.push_back(rec2);
	    continue;
	  }

          auto r1reverse = seqan3::align_pairwise(std::tie(seqan3::get<seqan3::field::seq>(rec1),spikeRC), config);
	  if(maxScore(r1reverse, minScore)){
	    //	    seqan3::debug_stream << "Read 1 Reverse Complement Spike" << std::endl;
	    read1RCSpikeCount++;
	    spikeCount++;

	    outSpikeFastq1.push_back(rec1);
	    outSpikeFastq2.push_back(rec2);
	    continue;
	  }

	  auto r2reverse = seqan3::align_pairwise(std::tie(seqan3::get<seqan3::field::seq>(rec2), spikeRC), config);
	  if(maxScore(r2reverse, minScore)){
	    //	    seqan3::debug_stream << "Read 2 Reverse Complement Spike" << std::endl;
	    read2RCSpikeCount++;
	    spikeCount++;

	    outSpikeFastq1.push_back(rec1);
	    outSpikeFastq2.push_back(rec2);
	    continue;
	  }
	  
   	  //write to fastqs
	  outFastq1.push_back(rec1);
	  outFastq2.push_back(rec2);

	  
    if (seqan3::get<seqan3::field::id>(rec1).substr(0,37) != seqan3::get<seqan3::field::id>(rec2).substr(0,37)){
	    throw std::runtime_error("Oh oh your pairs don't match.");
	    }
	    }

    std::cout << "Number of Spikes: " << spikeCount << std::endl
              << "Number of Input Reads: " << count << std::endl
              << "R1 Spikes " << read1SpikeCount << std::endl
              << "R2 Spikes " << read2SpikeCount << std::endl
              << "R1 RC Spikes " << read1RCSpikeCount << std::endl
              << "R2 RC Spikes " << read2RCSpikeCount << std::endl;
        
    return 0;
}
