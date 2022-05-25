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



int main(int argc, char* argv[])
{
    

    // Load reference from fasta file
    {
        FReference ref("chr22.fa"); // load fasta file
        std::cout << "Reference sequence length: " << ref.Sequence.length() << std::endl;
        // print first 100bp
        std::cout << ref.Sequence.substr(0, 100) << std::endl;
       //std::cout << ref.Sequence<< std::endl;
        
        //create a suffix array for the sequence
        int n = ref.Sequence.length();
        char *txt;
        txt = new char [n+1];
        //char txt[n + 1];
        strcpy(txt, ref.Sequence.c_str());
        
//        int *suffixArr;
//        suffixArr = new int [n];
//        suffixArr = buildSuffixArray(txt,  n);
        int *suffixArr = buildSuffixArray(txt,  n);
        //cout << "Following is suffix array for " << txt << endl;
        //printArr(suffixArr, n);
        //save the result
        saveArr(suffixArr, n, "suffixarray.txt");
        cout<<"Done!"<<endl;
        
        
        
    }
    

    return 0;
}
