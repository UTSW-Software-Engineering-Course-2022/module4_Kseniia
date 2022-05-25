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
        Sequence.append("$");
    }

public:
    std::string Sequence;

};

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
void printArr(int arr[], int n)
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



// A suffix array based search function to search a given pattern
// 'pat' in given text 'txt' using suffix array suffArr[]

void search(char *pat, char *txt, int *suffArr, int n, string output_file, int queryNum)
{
    std::ofstream outfile;
    outfile.open("output_file", ios::app);
    int lineNum = 0;
    int m = strlen(pat); // get length of pattern, needed for strncmp()
    
    // Do simple binary search for the pat in txt
    
    // Initilize left and right indexes
    int l = 0, r = n-1;
    while (l <= r)
    {
        // Compare pat with the middle suffix in suffix array
        int mid = l + (r - l)/2;
        int res = strncmp(pat, txt+suffArr[mid], m);

        // If match found write a line into the output file
        if (res == 0)
        {
                        
            outfile <<lineNum<<" "<<"query"<<queryNum<<" sample "<<suffArr[mid]+1<<"\n";
            lineNum = lineNum+1;
            
        }
        
        outfile.close();
        // Move to left half if pattern is alphabtically less than the mid suffix
        if (res < 0) r = mid - 1;

        // Otherwise move to right half
        else l = mid + 1;
    }
    

}

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

// Structure to store data of a rotation
struct rotation {
    int index;
    char* suffix;
};
// Compares the rotations and
// sorts the rotations alphabetically
int cmpfunc(const void* x, const void* y)
{
    struct rotation* rx = (struct rotation*)x;
    struct rotation* ry = (struct rotation*)y;
    return strcmp(rx->suffix, ry->suffix);
}

// Takes suffix array and its size
// as arguments and returns the
// Burrows - Wheeler Transform of given text
char* BWT(char* input_text,
                   std::vector<int> suffix_arr, int n)
{
    // Iterates over the suffix array to find
    // the last char of each cyclic rotation
    char* bwt_arr = (char*)malloc(n * sizeof(char));
    int i;
    for (i = 0; i < n; i++) {
        // Computes the last char which is given by
        // input_text[(suffix_arr[i] + n - 1) % n]
        int j = suffix_arr[i] - 1;
        if (j < 0)
            j = j + n;
  
        bwt_arr[i] = input_text[j];
    }
  
    bwt_arr[i] = '\0';
  
    // Returns the computed Burrows - Wheeler Transform
    return bwt_arr;
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

/**
 *
 *
 */
int main(int argc, char* argv[])
{
    

    // Load reference from fasta file
    {
        // sa InReferenceFastaFile OutSuffixArrayFile
        if (argc < 3) {
            PrintUsage(GetFilename(argv[0]));
            return 1;
        }
          
        
       
        // load fasta file
        auto start = high_resolution_clock::now();
        FReference ref(argv[1]);
    
        //print first 100bp
        std::cout << "Reference sequence length: " << ref.Sequence.length() << std::endl;
        std::cout << ref.Sequence.substr(0, 100) << std::endl;
        int n = ref.Sequence.length();
        char *txt;
        txt = new char [n+1];
        strcpy(txt, ref.Sequence.c_str());
        
        //load suffix array file
        std::vector <int> suffixArr = loadSuffixArray(argv[2]);
        
        // Apply Burrows-Wheeler Transform algorithm
        char* bwt_arr = BWT(txt, suffixArr, n);
      
        //save output of Burrows-Wheeler Transform
        ofstream myfile;
        myfile.open (argv[3]);
        myfile <<bwt_arr<<"\n";
        myfile.close();

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << duration.count() << endl;
        cout<<"Done!"<<endl;
    
    }
    
    return 0;
}
