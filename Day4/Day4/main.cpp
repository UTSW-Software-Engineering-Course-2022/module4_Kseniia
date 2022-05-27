/**
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <vector>
#include <chrono>


using namespace std::chrono;
using namespace std;

// TODO: support multiple chromosome
class FReference
{
public:
    // Default constructor
    FReference() {};

    // Constructor with fasta filename
    FReference(const std::string& InFilename) {
        LoadFromFasta(InFilename);
    }

    void LoadFromString(const std::string& InSequence)
    {
        Sequence = InSequence;
        if (Sequence.back() != '$') {
            Sequence.push_back('$');
        }
    }
    
    void LoadFromFasta(const std::string& InFilename)
    {
        std::ifstream fs(InFilename);
        if (!fs) {
            std::cerr << "Can't open file: " << InFilename << std::endl;
            return;
        }

        std::string buf;

        Sequence.clear();

        while(getline(fs, buf)) {
            if (buf.length() == 0) {
                // skip empty line
                continue;
            }

            if (buf[0] == '>') {
                // header line
                // TODO: save chromosome name
                continue;
            }

            Sequence.append(buf);
        }

        // append '$' as End of Sequence mark
        //Sequence.append("$");
    }

public:
    std::string Sequence;

};

//Load queries, create fector with quary sequences
std::vector<string>LoadQueryFromFasta(const std::string& InFilename)
{
    std::ifstream fs(InFilename);
    
    std::vector<string> queryVect;
    std::string buf;
    
    while(getline(fs, buf)) {
        if (buf.length() == 0) {
            // skip empty line
            continue;
        }
        if (buf[0] == '>') {
            //header line
            continue;
        }

        queryVect.push_back(buf);
        
    }
    
    return queryVect;
}
/*
 * Return filename from full path of file.
 *
 * InFilename   [In] Full path of file
 *
 * Return
 *   base filename
 *
 */
static std::string GetFilename(const std::string& InFilename)
{
    const size_t pos = InFilename.find_last_of("/\\");
    if (pos == std::string::npos) {
        return InFilename;
    }

    return InFilename.substr(pos + 1);
}

void PrintUsage(const std::string& InProgramName)
{
    std::cerr << "Invalid Parameters" << std::endl;
    std::cerr << "  " << InProgramName << " Reference_Fasta_File SuffixArray_File" << std::endl;
}

 
// Structure to store information of a suffix
struct suffix
{
    int index; // original index
    int rank[2]; // ranks
};
 

// Compares two pairs, returns 1 if first pair is smaller
int cmp(struct suffix a, struct suffix b)
{
    return (a.rank[0] == b.rank[0])? (a.rank[1] < b.rank[1] ?1: 0):
               (a.rank[0] < b.rank[0] ?1: 0);
}
 
class vector;
// build and return the suffix array for the given string
int *buildSuffixArray(char *txt, int n)
{
    // store suffixes and their indexes
    
    std::vector<struct suffix> suffixes(n);
    
    //struct suffix suffixes[n];
 
    for (int i = 0; i < n; i++)
    {
        suffixes[i].index = i;
        suffixes[i].rank[0] = txt[i] - 'a';
        suffixes[i].rank[1] = ((i+1) < n)? (txt[i + 1] - 'a'): -1;
    }
 
    // Sort the suffixes
    sort(suffixes.begin(), suffixes.end(), cmp);
    
    int *ind;
    ind = new int [n];
    //int ind[n];
    for (int k = 4; k < 2*n; k = k*2)
    {
        // Assigning rank and index values to first suffix
        int rank = 0;
        int prev_rank = suffixes[0].rank[0];
        suffixes[0].rank[0] = rank;
        ind[suffixes[0].index] = 0;
 
        // Assigning rank to suffixes
        for (int i = 1; i < n; i++)
        {
            // If first rank and next ranks are same as that of previous
            // suffix in array, assign the same new rank to this suffix
            if (suffixes[i].rank[0] == prev_rank &&
                    suffixes[i].rank[1] == suffixes[i-1].rank[1])
            {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = rank;
            }
            else // Otherwise increment rank and assign
            {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = ++rank;
            }
            ind[suffixes[i].index] = i;
        }
 
        // Assign next rank to every suffix
        for (int i = 0; i < n; i++)
        {
            int nextindex = suffixes[i].index + k/2;
            suffixes[i].rank[1] = (nextindex < n)?
                                  suffixes[ind[nextindex]].rank[0]: -1;
        }
 
        // Sort the suffixes according to first k characters
        sort(suffixes.begin(), suffixes.end(), cmp);
    }
 
    // Store indexes of all sorted suffixes in the suffix array
    int *suffixArr = new int[n];
    for (int i = 0; i < n; i++)
        suffixArr[i] = suffixes[i].index;
        
    free(ind);
    // Return the suffix array
    return  suffixArr;
}
 
// print array
void printArr(std::vector <int> arr, int n)
{
    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";
    cout << endl;
}

//save array
void saveArr(int arr[], int n, string output_path)
{
    ofstream myfile;
    myfile.open (output_path);
    
    for (int i = 0; i < n; i++)
        myfile << i+1<<" "<<arr[i]<<"\n";
    myfile.close();
    
}





std::vector <int> loadSuffixArray(const string& suffixarrayFile)
{
    string line;
    
    std::vector <int> index;
    std::ifstream fs(suffixarrayFile);
    if (!fs) {
        std::cerr << "Can't open file: " << suffixarrayFile << std::endl;
    }
    if (fs.is_open())
    {
        while (getline(fs, line))
        {
           index.push_back(stoi(line));
        }
        fs.close();
    }
    
    return index;
    
}

//maps character to an integer
size_t get_index(char const chr)
{
    switch (chr)
    {
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        default: throw std::logic_error{"wrong character"};
    }
}

//retuns count table with counts for each nucleotide (and $)
std::vector<uint16_t> CountTable(std::string const & bwt)
{
    std::vector<uint16_t> count_table(5);
 
    for (auto chr : bwt)
    {
        for (size_t i = get_index(chr) + 1; i < 5; ++i)
            ++count_table[i];
    }
 
    return count_table;
}

//defines compact dat atructure Bivecotrs
//this struct was published by www.fu-berlin.de
struct Bitvector
{
    std::vector<uint64_t> data;
    std::vector<uint16_t> blocks;
    std::vector<uint64_t> superblocks;
 
    uint64_t block_size;
    uint64_t superblock_size;
 
    Bitvector(size_t const count)
    {
        data.resize((count + 63) / 64); // the +63 are a trick to round up the fraction.
    }
 
    bool read(size_t const i) const
    {
        return (data[i / 64] >> (63 - (i % 64))) & 1;
    }
 
    void write(size_t const i, bool const value)
    {
        uint64_t mask = static_cast<uint64_t>(1) << (63 - (i % 64));
 
        if (value == true)
            data[i / 64] |= mask;
        else
            data[i / 64] &= ~mask;
    }
 
    void construct(size_t const new_block_size = 64, size_t const new_superblock_size = 512)
    {
        block_size = new_block_size;
        superblock_size = new_superblock_size;
 
        size_t number_of_bits = data.size() * 64;
 
        blocks.resize((number_of_bits + block_size - 1) / block_size, 0);
        superblocks.resize((number_of_bits + superblock_size - 1) / superblock_size, 0);
 
        size_t block_pos{0};
        size_t super_block_pos{0};
 
        uint16_t block_count{0};
        uint64_t super_block_count{0};
 
        for (size_t i = 0; i < number_of_bits; ++i)
        {
            if (i % block_size == 0)
            {
                if (i % superblock_size == 0)
                {
                    super_block_count += block_count; // update superblock count
 
                    superblocks[super_block_pos] = super_block_count;
 
                    ++super_block_pos; // move to the next position
                    block_count = 0;   // reset block count
                }
 
                blocks[block_pos] = block_count;
 
                ++block_pos; // move to the next position
            }
 
            if (read(i) == true)
                ++block_count;
        }
    }
 
    uint64_t rank(size_t const i) const
    {
        uint64_t rank{0};
 
        rank += superblocks[i / superblock_size];
        rank += blocks[i / block_size];
 
        for (size_t j = (i / block_size) * block_size; j < i; ++j)
        {
            rank += read(j);
        }
 
        return rank;
    }
};

//constructs occurance table - fills five bivectors using bwt
/*
 for ex.
   A T C C G $ T A C C A A
  -----------------------
$  0 0 0 0 0 1 0 0 0 0 0 0
A  1 0 0 0 0 0 0 1 0 0 1 1
G  0 0 0 0 1 0 0 0 0 0 0 0
T  0 1 0 0 0 0 1 0 0 0 0 0
C  0 0 1 1 0 0 0 0 1 1 0 0
 */
struct occurrence_table
{
    
    std::vector<Bitvector> data;
 
    occurrence_table(std::string const & bwt)
    {
       
        data.resize(5, Bitvector(bwt.size()));
 
        
        for (size_t i = 0; i < bwt.size(); ++i)
            data[get_index(bwt[i])].write(i, 1);
 
        for (Bitvector & bitv : data)
            bitv.construct(3, 6);
    }
 
    size_t read(char const chr, size_t const i) const
    {
        return data[get_index(chr)].rank(i + 1);
    }
};

//returns occurences of the queried motif
//by performing backward search
std::vector<int>  align_query(std::string const & query, std::string const & bwt, std::vector<uint16_t> const & C,
             occurrence_table const & OCC)
{   std::vector<int> suffixArray_range;
    int64_t i = query.size() - 1;
    size_t start = 0;
    size_t end = bwt.size() - 1;
 
    while ((start <= end) && (i >= 0))
    {
        char c = query[i];
        start = C[get_index(c)] + (start ? OCC.read(c, start - 1) : 0);
        end = C[get_index(c)] + OCC.read(c, end) - 1;
        
        i = i - 1;
    }
    
    
    if (end < start)
        return suffixArray_range;
    else
        //retyrb index of suffix array entries corresponding to the motif position
        for (int i2=start; i2<=end; i2++){
            suffixArray_range.push_back(i2);
        }
    return suffixArray_range;
    
}
/**
 *
 *
 */
int main(int argc, char* argv[])
{
    
        // load BW transformed sequence
        FReference bwt_arr(argv[1]);
        string bwt = bwt_arr.Sequence;
        
    
        //create count table with counts for each nucleotide (and $)
        std::vector<uint16_t> C = CountTable(bwt);
        
         
        //load suffix array file
        std::vector <int> suffixArr = loadSuffixArray(argv[2]);
        
        //compute occurance table using bwt
        occurrence_table OCC(bwt);

        
        //load fasta file with queries
        std::vector<string>queryVector;
        queryVector = LoadQueryFromFasta(argv[3]);
        
        int nq=queryVector.size();
        int ii=0;
        int lineNum=0;
        
        std::ofstream outfile;
        outfile.open(argv[4]);
                     
        while (ii < nq)
        {
            char quer[queryVector[ii].length()];
            strcpy(quer, queryVector[ii].c_str());
            ii = ii+1;
            
            //get a vecotr with suffix array indeces for each query occurance
            std::vector<int> suffixArray_pos = align_query(quer, bwt, C, OCC);
            if (suffixArray_pos.size()>0) {
                int st = suffixArray_pos.at(0);
                int i2=st;
                for (i2=st; i2<(st+suffixArray_pos.size()); i2++){
                    //save query match to the output file
                    outfile <<lineNum<<" "<<"query"<<ii<<" sample "<<suffixArr.at(i2)+1<<"\n";
                    lineNum = lineNum+1;
                   
                }
                
            }  else {
              continue;
            }
            
        
        }
       
    outfile.close();
    
    
    return 0;
}
