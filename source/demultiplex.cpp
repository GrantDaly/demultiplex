#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <seqan3/std/concepts>
#include <map>
#include <unordered_map>
#include <utility>
#include <any>

#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/configuration/all.hpp> 

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>

#include <seqan3/search/search.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/alphabet/nucleotide/concept.hpp>
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
  auto tempScore = 0;
  for(auto & res : readResults){
   //Todo
    tempScore = res.score();
    //for edit distance, a positive value is invalid
   if(tempScore <= scoreVal){ 
      //          seqan3::debug_stream << res << std::endl;
	  return true;
	    }
	    }
 return false;
}


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
    //only validates for .gz, could change
    myparser.add_option(fastqOneName, 'p', "fastq-1","Fastq One",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});
    myparser.add_option(fastqTwoName, 'q', "fastq-2","Fastq Two",
                        seqan3::option_spec::DEFAULT, seqan3::input_file_validator{{"fq", "fastq", "gz"}});

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
    //iterate throught the fastq files
    for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) // && is important!
	{                                                           // because seqan3::views::zip returns temporaries
	  ++count;
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
    return 0;
}
