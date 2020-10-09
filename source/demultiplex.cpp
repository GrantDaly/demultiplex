#include <string>
<<<<<<< HEAD
#include <filesystem>
#include <iostream>
#include <fstream>
#include <seqan3/std/concepts>
#include <map>
#include <unordered_map>
#include <utility>
#include <any>
=======
#include <experimental/filesystem>
#include <iostream>
#include <seqan3/std/concepts>
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079

#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/configuration/all.hpp> 
<<<<<<< HEAD

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>

=======
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>

#include <seqan3/alignment/pairwise/all.hpp>
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
#include <seqan3/search/search.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/alphabet/nucleotide/concept.hpp>
<<<<<<< HEAD
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>


using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

struct Sample {
  std::string sampleName;
  std::string iSevenIndex;
  std::string ciFiveIndex;
};

std::istream& operator>>(std::istream& is, Sample& s)
   {
     
     is >> s.sampleName >> s.iSevenIndex >> s.ciFiveIndex;
     return is;
    }

std::ostream& operator<<(std::ostream& os, Sample& s)
{
  os << s.sampleName << "\t" << s.iSevenIndex << "\t" << s.ciFiveIndex;
  return os;
    }

std::vector<Sample> readinIndexes(std::filesystem::path indexFileName)
{

  std::ifstream indexFile{indexFileName};
  //indexFile.open(indexFileName, std::ios::in);

  //seqan3::sequence_file_input inFasta{filename};
  //  std::vector<seqan3::dna5_vector> sequences;
  std::vector<Sample> samples;
  
  if(indexFile) {
    while(indexFile.good()) {
      Sample temp;
      indexFile >> temp;
      //std::cout << temp << std::endl;
      samples.push_back(temp);

    }
  }
  else {
    std::cerr << "can't open barcode file" << std::endl;
  }

  //makes a copy so may want to do by reference for large file
  return samples;
}

// auto createOutFiles(std::vector<Sample> indexVec, std::string samplePrefix,
// 		    std::filesystem::path outDirectory, std::string samplePostfix) {

//   /*putting the relevant code in main for now because I'm having trouble putting both in unordered map */
//   std::cout << "sample prefix\t" << samplePrefix << "\tOutDirectory\t" << outDirectory.string() << std::endl;

//   //using out_type = seqan3::sequence_file_output<>;
//   std::unordered_map<std::string, seqan3::sequence_file_output<>> fileHandleMap;
//   for(auto s : indexVec) {
//     auto outFilePrefix = samplePrefix + s.sampleName;
//     std::cout << outFilePrefix << std::endl;

//     auto outName = outFilePrefix + samplePostfix;
//     //    auto outNameTwo = outFilePrefix + "_R2.fq.gz";
    
//     seqan3::sequence_file_output tmpOut{outName};
//     //seqan3::sequence_file_output tmpOutTwo{outNameTwo};
        
//     //auto filePair = std::make_pair(std::move(tmpOutOne), std::move(tmpOutTwo));
//     //auto filePair = 
//     std::string tempSample = s.sampleName;
//     //    fileHandleMap.insert(std::make_pair(std::move(s.sampleName), std::move(tmpOutTwo)));
//       auto filePair = std::make_pair(std::move(tempSample), std::move(tmpOut));
//       fileHandleMap.insert(filePair);
    
//   }

//   return fileHandleMap;

// }


seqan3::dna5_vector stringToDNA(std::string inString) {
  //took literal overload from dna5 source code
     seqan3::dna5_vector outVector;
     auto n = inString.length();
     outVector.resize(n);

     
     for (size_t i = 0; i < n; ++i)
     {
     //    outVector[i].assign_char(s[i]);
       outVector[i].assign_char(inString.at(i));
      
     }
       
       //seqan3::debug_stream << outVector << std::endl;
   return outVector;

}

seqan3::dna5_vector rcDNAvector(seqan3::dna5_vector dnaVector) {
  //auto vecSize = inVector.size();
  //seqan3::dna5_vector outVector{};
  //outVector.resize(vecSize);
  
  seqan3::dna5_vector reversedInput = dnaVector | 
    seqan3::views::complement | std::views::reverse |			 
    seqan3::views::to<seqan3::dna5_vector>;
  return reversedInput;

}

template<typename T, std::integral ScoreType>
bool anyGoodAlignments(T & readResults, ScoreType scoreVal){
=======

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
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
  auto tempScore = 0;
  for(auto & res : readResults){
   //Todo
    tempScore = res.score();
    //for edit distance, a positive value is invalid
<<<<<<< HEAD
   if(tempScore <= scoreVal){ 
=======
   if((tempScore >= scoreVal) && (tempScore < 1)){ 
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
      //          seqan3::debug_stream << res << std::endl;
	  return true;
	    }
	    }
 return false;
}

<<<<<<< HEAD

int main(int argc, char ** argv)
{
  seqan3::argument_parser myparser{"Demultiplexer", argc, argv}; // initialize
  
    std::filesystem::path indexFileName;
    std::filesystem::path fastqOneName;
    std::filesystem::path fastqTwoName;
    //    std::filesystem::directory_entry outDirectory;
    std::filesystem::path outDirectory;
    myparser.add_option(indexFileName, 'f', "barcodes","Barcode File",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"tsv"}});
=======
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
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
    //only validates for .gz, could change
    myparser.add_option(fastqOneName, 'p', "fastq-1","Fastq One",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});
    myparser.add_option(fastqTwoName, 'q', "fastq-2","Fastq Two",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});
<<<<<<< HEAD

    myparser.add_option(outDirectory, 'd', "out-dir","Output Directory",
    seqan3::option_spec::DEFAULT);

   try {
myparser.parse();
    }
   catch (seqan3::argument_parser_error const & ext) { // the user did something wrong 
        std::cerr << "Index Demultiplexer - [PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

   if(! std::filesystem::exists(outDirectory)) {
     std::filesystem::create_directory(outDirectory);
   }

    //read in sample and index file
    auto indexVec = readinIndexes(indexFileName);

    //    auto outFileMap = createOutFiles(std::vector<Sample> indexVec, samplePrefix);
    auto sampleNameOnly = fastqOneName.filename().string();
    auto samplePrefix = sampleNameOnly.substr(0, sampleNameOnly.find_last_of("_") + 1);
    //auto fileHandleMap = createOutFiles(indexVec, samplePrefix, outDirectory, "R1.fq.gz");
    

    for(auto s : indexVec) {
    auto outFilePrefix = samplePrefix + s.sampleName;

    seqan3::debug_stream << "Pulling out sample " << s.sampleName << std::endl;
    
    seqan3::dna5_vector sampleISeven = stringToDNA(s.iSevenIndex);
    seqan3::dna5_vector sampleCIFive = stringToDNA(s.ciFiveIndex);
    seqan3::dna5_vector sampleIFive = rcDNAvector(sampleCIFive);


    
    auto outNameOne = outFilePrefix + "_R1.fq.gz";
    auto outNameTwo = outFilePrefix + "_R2.fq.gz";

    //out fastqs
    seqan3::sequence_file_output outFastqOne{outDirectory / outNameOne};
    seqan3::sequence_file_output outFastqTwo{outDirectory / outNameTwo};
    									     
    //paired fastqs
        
    seqan3::sequence_file_input fin1{fastqOneName};
    seqan3::sequence_file_input fin2{fastqTwoName};

    long long int count = 0;
=======
    
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
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
    //iterate throught the fastq files
    for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) // && is important!
	{                                                           // because seqan3::views::zip returns temporaries
	  ++count;
<<<<<<< HEAD
	  if(count % 1000000 == 0){
	    seqan3::debug_stream << count / 1000000.0 << " Million  Reads processed" << std::endl;
	  }

	  auto headerOne = seqan3::get<seqan3::field::id>(rec1);
	  
	  //seqan3::debug_stream << headerOne << std::endl;

	  //get index sequences from read
	  std::string right = headerOne.substr(headerOne.find_last_of(" ") + 1, std::string::npos );
	  std::string readName = headerOne.substr(0, headerOne.find_last_of(" "));
	  //auto startOfIndex = headerOne.find_last_of(":") + 1;
	  auto indicesStart = right.find_last_of(":") + 1;
	  auto indexDelimiter = right.find("+");
	  	  
	  auto iSevenIndex = right.substr(indicesStart, 8);
	  
	  auto ciFiveIndex = right.substr(indexDelimiter +1, std::string::npos);
	  
	  //std::cout << "i7\t" << iSevenIndex << "\tci5\t" <<  ciFiveIndex << std::endl;
	  
	  seqan3::dna5_vector iSevenDNA = stringToDNA(iSevenIndex);
	  seqan3::dna5_vector ciFiveDNA = stringToDNA(ciFiveIndex);
	  seqan3::dna5_vector iFiveDNA = rcDNAvector(ciFiveDNA);
	  //seqan3::debug_stream << "i7\t" << iSevenDNA << "\tci5\t" <<  ciFiveIndex <<
	  // "\ti5\t" << iFiveDNA << std::endl;

    
	  
	  // Compute the alignment of a single pair.
	  auto const edit_config = seqan3::align_cfg::edit | seqan3::align_cfg::max_error{1u};
	  auto res1 = seqan3::align_pairwise(std::tie(sampleISeven, iSevenDNA), edit_config);

	  //seqan3::debug_stream << "The score: " << res1.score() << "\n";
	  if(anyGoodAlignments(res1, 1)) {

	      //if first index worked, check second index
	      auto res2  = seqan3::align_pairwise(std::tie(sampleCIFive, iFiveDNA), edit_config);
	      if(anyGoodAlignments(res2, 1)) {

		//seqan3::debug_stream << "********************" << std::endl;
		// seqan3::debug_stream << "Sample i7 " << sampleISeven  << "Sample i5 " << sampleCIFive << std::endl;
		// seqan3::debug_stream << " Read i7  " << iSevenDNA << std::endl << " Read i5  " << iFiveDNA << std::endl;
		//   //seqan3::debug_stream << "The score: " << *res2.score() << "\n";
		outFastqOne.push_back(rec1);
		outFastqTwo.push_back(rec2);
		}

	    }

	  

	  
	}
    }
=======
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
        
>>>>>>> 6c199f1c09578192530536361b175a63ce2e2079
    return 0;
}
