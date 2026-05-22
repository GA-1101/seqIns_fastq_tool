#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <cctype>

// ============================================
// 1. Genetic code table
// ============================================
const std::unordered_map<std::string, char> CODON_TABLE = {
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
    {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};

// ============================================
// 2. Utility function: Convert string to uppercase
// ============================================
std::string toUpper(const std::string& s) {
    std::string result = s;
    for (char& c : result) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
    return result;
}

// ============================================
// 3. Core function: Translate DNA to protein
// ============================================
std::string translateDNA(const std::string& dna) {
    std::string protein;
    protein.reserve(dna.length() / 3 + 1);
    
    for (size_t i = 0; i + 2 < dna.length(); i += 3) {
        std::string codon = dna.substr(i, 3);
        auto it = CODON_TABLE.find(codon);
        if (it != CODON_TABLE.end()) {
            protein += it->second;
        } else {
            protein += 'X';  // Unknown codon represented by X
        }
    }
    return protein;
}

// ============================================
// 4. Core function: Extract sequence between two markers
// ============================================
std::string findInsertString(const std::string& seq, 
                              const std::string& seq1, 
                              const std::string& seq2) {
    // Find the position of seq1
    size_t pos1 = seq.find(seq1);
    if (pos1 == std::string::npos) {
        return "";  // seq1 not found
    }
    
    // Find the position of seq2 after seq1
    size_t start = pos1 + seq1.length();
    size_t pos2 = seq.find(seq2, start);
    if (pos2 == std::string::npos) {
        return "";  // seq2 not found
    }
    
    return seq.substr(start, pos2 - start);
}

// ============================================
// 5. File reading: Parsing FASTQ files
// ============================================
std::vector<std::string> readFASTQ(const std::string& filename) {
    std::vector<std::string> reads;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        std::exit(1);
    }
    
    std::string line;
    int lineNumber = 0;
    
    while (std::getline(file, line)) {
        // FASTQ format:Groups of 4 lines:
        // Line 1: @ (identifier)
        // Line 2: Sequence
        // Line 3: + (separator)
        // Line 4: Quality scores
        if (lineNumber % 4 == 1) {  // Line 2 is the sequence
            reads.push_back(line);
        }
        lineNumber++;
    }
    
    file.close();
    std::cout << "Successfully read " << reads.size() << " sequences" << std::endl;
    return reads;
}

// ============================================
// 6. Core processing: Analyze insert sequences in each read
// ============================================
struct InsertResult {
    std::string insertSeq;   // Inserted DNA sequence (or status marker)
    std::string insertPro;   // Translated protein sequence (or status marker)
    int insertLen;           // Length of the inserted sequence
};

std::vector<InsertResult> findInsertSeq(
    const std::vector<std::string>& reads,
    const std::string& seq1,
    const std::string& seq2) {
    
    std::vector<InsertResult> results;
    results.reserve(reads.size());
    
    for (const std::string& read : reads) {
        InsertResult result;
        
        // Check if both marker sequences are present
        bool hasSeq1 = (read.find(seq1) != std::string::npos);
        bool hasSeq2 = (read.find(seq2) != std::string::npos);
        
        if (hasSeq1 && hasSeq2) {
            std::string insert = findInsertString(read, seq1, seq2);
            
            if (insert.empty()) {
                // Found both markers, but no sequence between them
                result.insertSeq = "Not_Insert";
                result.insertPro = "Not_Insert";
                result.insertLen = 0;
            } else {
                result.insertSeq = insert;
                result.insertPro = translateDNA(insert);
                result.insertLen = insert.length();
            }
        } else {
            // Does not contain both markers
            result.insertSeq = "Not_Match";
            result.insertPro = "Not_Match";
            result.insertLen = -1;  // Let -1 represent a no-match.
        }
        
        results.push_back(result);
    }
    
    return results;
}

// ============================================
// 7. Statistical function: Count elements
// ============================================
std::unordered_map<std::string, int> countElements(
    const std::vector<std::string>& items) {
    std::unordered_map<std::string, int> counts;
    for (const std::string& item : items) {
        counts[item]++;
    }
    return counts;
}

// ============================================
// 8. Output function: Write to CSV file (sorted by Quantity in descending order)
// ============================================
void writeCSV(const std::unordered_map<std::string, int>& data,
              const std::string& filename,
              const std::string& keyHeader,
              const std::string& valHeader) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to create file " << filename << std::endl;
        return;
    }
    
    // Convert the unordered_map to a vector for sorting.
    std::vector<std::pair<std::string, int>> sortedData(data.begin(), data.end());
    
    // Sort by Quantity in descending order
    std::sort(sortedData.begin(), sortedData.end(),
              [](const std::pair<std::string, int>& a,
                 const std::pair<std::string, int>& b) {
                  return a.second > b.second;  // Descending order
              });
    
    file << keyHeader << "," << valHeader << "\n";
    for (const auto& pair : sortedData) {
        file << pair.first << "," << pair.second << "\n";
    }
    file.close();
    std::cout << "Generated: " << filename << " (Sorted by Quantity in descending order)" << std::endl;
}

// ============================================
// 9. Output function: Write to list file (TXT)
// ============================================
void writeList(const std::vector<std::string>& items,
               const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to create file " << filename << std::endl;
        return;
    }
    
    for (const std::string& item : items) {
        file << item << "\n";
    }
    file.close();
    std::cout << "Generated: " << filename << std::endl;
}

// ============================================
// 10. Output function: Write to JSON file
// ============================================
void writeJSON(const std::unordered_map<std::string, int>& data,
               const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to create file " << filename << std::endl;
        return;
    }
    
    file << "{\n";
    bool first = true;
    for (const auto& pair : data) {
        if (!first) file << ",\n";
        file << "  \"" << pair.first << "\": " << pair.second;
        first = false;
    }
    file << "\n}\n";
    file.close();
    std::cout << "Generated: " << filename << std::endl;
}

// ============================================
// 11. Main function
// ============================================
int main(int argc, char* argv[]) {
    // Record start time
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Check command line arguments
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] 
                  << " <fastq file> <sequence1> <sequence2> <output prefix>" << std::endl;
        std::cerr << "Example: " << argv[0] 
                  << " sample.fastq ATGCCG TAGGCA output" << std::endl;
        return 1;
    }
    
    // Parse arguments
    std::string seqsFile = argv[1];
    std::string sequence1 = toUpper(argv[2]);
    std::string sequence2 = toUpper(argv[3]);
    std::string outPrefix = argv[4];
    
    std::cout << "=============================================" << std::endl;
    std::cout << "FASTQ Insert Sequence Analysis Tool v1.0" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Input file: " << seqsFile << std::endl;
    std::cout << "Marker sequence 1: " << sequence1 << std::endl;
    std::cout << "Marker sequence 2: " << sequence2 << std::endl;
    std::cout << "Output prefix: " << outPrefix << std::endl;
    std::cout << std::endl;
    
    // 1. Read FASTQ file
    std::cout << "[1/5] Reading FASTQ file..." << std::endl;
    std::vector<std::string> reads = readFASTQ(seqsFile);
    
    // 2. Find insert sequences
    std::cout << "[2/5] Analyzing insert sequences..." << std::endl;
    std::vector<InsertResult> results = findInsertSeq(reads, sequence1, sequence2);
    
    // 3. Extract various data
    std::cout << "[3/5] Statistics..." << std::endl;
    std::vector<std::string> seqList, proList, lenList;
    for (const auto& r : results) {
        seqList.push_back(r.insertSeq);
        proList.push_back(r.insertPro);
        lenList.push_back(std::to_string(r.insertLen));
    }
    
    // 4. Calculate statistical information
    auto seqStats = countElements(seqList);
    auto proStats = countElements(proList);
    
    // Length statistics (requires special handling)
    std::unordered_map<std::string, int> lenStats;
    for (const auto& r : results) {
        if (r.insertLen > 0) {
            lenStats[std::to_string(r.insertLen)]++;
        } else if (r.insertLen == 0) {
            lenStats["0"]++;
        }
    }
    
    // 5. Generate output files
    std::cout << "[4/5] Generating output files..." << std::endl;
    
    writeCSV(seqStats, outPrefix + "_Insert_stats.csv", "Seq_Insert", "Quantity");
    writeCSV(lenStats, outPrefix + "_Insert_length_stats.csv", "Seq_Insert_Length", "Quantity");
    writeCSV(proStats, outPrefix + "_Insert_protein_stats.csv", "Seq_Insert_Protein", "Quantity");
    
    writeList(seqList, outPrefix + "_Insert_list.txt");
    writeList(lenList, outPrefix + "_Insert_length.txt");
    writeList(proList, outPrefix + "_Insert_protein.txt");
    
    writeJSON(seqStats, outPrefix + "_Insert_stats.json");
    writeJSON(lenStats, outPrefix + "_Insert_length_stats.json");
    writeJSON(proStats, outPrefix + "_Insert_protein_stats.json");
    
    // Print execution time
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        endTime - startTime).count();
    
    std::cout << "[5/5] Complete! Total time: " << duration / 1000.0 << " seconds" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    // Print statistical summary
    std::cout << "\n📊 Quick Statistics:" << std::endl;
    int matchCount = 0, noInsertCount = 0, notMatchCount = 0;
    for (const auto& r : results) {
        if (r.insertSeq == "Not_Match") notMatchCount++;
        else if (r.insertSeq == "Not_Insert") noInsertCount++;
        else matchCount++;
    }
    
    std::cout << "  Total sequences: " << results.size() << std::endl;
    std::cout << "  Successful matches: " << matchCount << std::endl;
    std::cout << "  No insert between markers: " << noInsertCount << std::endl;
    std::cout << "  Unmatched: " << notMatchCount << std::endl;
    
    return 0;
}
