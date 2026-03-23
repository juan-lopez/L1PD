#include <limits.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include "cxxopts.hpp"
#include <omp.h>
#include <fstream>

using namespace std;
// This is the consolidated script so far work is still being done on the individual parts
// Outline error prone areas for future reference

/*
Precision & Recall

Authors:
    Juan O. Lopez (juano.lopez@upr. edu)
    Emanuel Martinez (emanuel.martinez8@upr.edu)
    Javier Quinones (javier.quinones3@upr.edu)

License: Creative Commons Attribution-ShareAlike 4.0
http://creativecommons.org/licenses/by-sa/4.0/1egalcode
*/
string get_slurm_dir()
{
    const char *slurm_job_id = std::getenv("SLURM_JOB_ID");
    if (slurm_job_id)
    {
        string command = "scontrol show job " + string(slurm_job_id) + " | awk -F= '/WorkDir=/{print $2}'";
        FILE *pipe = popen(command.c_str(), "r");
        if (!pipe)
            return "";
        char buffer[256];
        string result;
        if (fgets(buffer, sizeof(buffer), pipe) != nullptr)
        {
            result = buffer;
            result.erase(result.find_last_not_of(" \n\r\t") + 1);
        }
        pclose(pipe);
        return result;
    }
    else
    {
        char result[PATH_MAX];
        ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
        if (count != -1)
        {
            string path(result, count);
            size_t last_slash = path.find_last_of('/');
            if (last_slash != string::npos)
            {
                return path.substr(0, last_slash);
            }
            return path;
        }
        else
        {
            return "";
        }
    }
}

void executePythonScript(const string script_path, const string args, bool time)
{
    string command = "python3 " + script_path + " " + args;
    if (time)
        command = "time " + command;

    int result = system(command.c_str());
    if (result != 0)
        cerr << "Error executing: " << command << endl;
}

int executeCommand(const string command, string &output)
{
    FILE *p;
    int ch;
    output.clear();

    p = popen(command.c_str(), "r");
    if (p == NULL)
    {
        puts("Unable to open process");
        return 1;
    }

    while ((ch = fgetc(p)) != EOF)
        output += ch;

    return pclose(p);
}

bool fileExists(const std::string &filename)
{
    std::ifstream infile(filename);
    return infile.good();
}

int count_occurrences_in_file(const std::string &filepath, const std::string &search_term)
{
    std::ifstream file(filepath);
    if (!file)
    {
        std::cerr << "Error opening file: " << filepath << std::endl;
        return 0;
    }

    std::string line;
    int count = 0;
    while (std::getline(file, line))
    {
        // Convert line and search_term to lowercase for case-insensitive search
        std::string lower_line = line;
        std::string lower_term = search_term;
        std::transform(lower_line.begin(), lower_line.end(), lower_line.begin(), ::tolower);
        std::transform(lower_term.begin(), lower_term.end(), lower_term.begin(), ::tolower);

        // Check if search_term is in line
        if (lower_line.find(lower_term) != std::string::npos)
        {
            count++;
        }
    }
    file.close();
    return count;
}

int main(int argc, char *argv[])
{
    // Define command line options with descriptions matching the original Bash script flags
    cxxopts::Options options("precision_recall", "Precision and Recall computation script");

    options.add_options()
        // -g: Genome file
        ("g,genome", "Genome file.", cxxopts::value<string>())
        // -f: Fasta file
        ("f,fasta", "Fasta file with all the probes.", cxxopts::value<string>())
        // -s: SAM directory
        ("s,sam", "Directory where SAM files will be generated.", cxxopts::value<string>())
        // -c: Meta directory
        ("c,meta", "Directory with LINE-1 meta data file in .csv format.", cxxopts::value<string>())
        // -l: Entries count
        ("l,entries", "Entries for the entire genome.", cxxopts::value<int>())
        // -p: Print directory, default "result"
        ("p,print", "Directory where results and raw files will be saved in.", cxxopts::value<string>()->default_value("result"))
        // -e: Edit distance range, default "5, 20, 5"
        ("e,edit", "Vary the edit distance for MRFAST.", cxxopts::value<string>()->default_value("5,20,5"))
        // -t: Threshold range (note typo 'treshold' kept for exact match), default "25, 500, 25"
        ("t,treshold", "Vary threshold for Precision and Recall.", cxxopts::value<string>()->default_value("25,500,25"))
        // -m: Minimum kmers, default 2
        ("m,MinKmers", "Minimum amount of kmers to be used at once.", cxxopts::value<int>()->default_value("2"))
        // -G: Generate graph, default false
        ("G,graph", "Generates a graph from the results of precision and recall.", cxxopts::value<bool>()->default_value("false"))
        // -x: x-axis headers
        ("x,x_axis", "Specifies the header name for the x axis of the graph.", cxxopts::value<vector<string>>())
        // -y: y-axis headers
        ("y,y_axis", "Specifies the header name for the y axis of the graph.", cxxopts::value<vector<string>>())
        // -z: z-axis headers
        ("z,z_axis", "Specifies the header name for the z axis of the graph.", cxxopts::value<vector<string>>())
        // -d: 3D graph flag, default false
        ("d,three_d", "Graph will be three dimensional.", cxxopts::value<bool>()->default_value("false"))
        // -s: Stackable graphs, default false.
        // Note: duplicate short flag 's' with sam directory, short flag removed.
        ("stackable", "Graphs will be displayed on the same plain.", cxxopts::value<bool>()->default_value("false"))
        // -o: Output graph filename, default "graph.png".
        ("o,output", "Name of the output graph", cxxopts::value<string>()->default_value("graph.png"))
        // -q: Specify the genomeType to be evaluated
        ("q,sequence", "Sequence type to be evaluated.", cxxopts::value<string>()->default_value("Line1"))
        // -h: Print help menu.
        ("h,help", "Print usage");

    // Parse command line arguments
    auto result = options.parse(argc, argv);

    // If help requested, print usage and exit
    if (result.count("help"))
    {
        cout << options.help() << endl;
        return 0;
    }

    // Extract values from parsed options
    // Required arguments (no default values)
    string genome = result["genome"].as<string>();
    string fasta = result["fasta"].as<string>();
    string sam = result["sam"].as<string>();
    string meta = result["meta"].as<string>();
    int l1count = result["entries"].as<int>();
    string SequenceType = result["sequence"].as<string>();

    // Optional arguments with defaults
    string printdir = result["print"].as<string>();
    bool graph = result["graph"].as<bool>();
    bool three_d = result["three_d"].as<bool>();
    bool stackable = result["stackable"].as<bool>();
    string output = result["output"].as<string>();

    // Optional vector arguments, default to single headers if not provided
    vector<string> x_axis = result.count("x_axis") ? result["x_axis"].as<vector<string>>() : vector<string>{"Edit dist"};
    vector<string> y_axis = result.count("y_axis") ? result["y_axis"].as<vector<string>>() : vector<string>{"Threshold"};
    vector<string> z_axis = result.count("z_axis") ? result["z_axis"].as<vector<string>>() : vector<string>{"F1 score"};

    // Parse edit distance range string into integers
    string edit_range = result["edit"].as<string>();
    int EditStart, EditEnd, EditStep;
    std::sscanf(edit_range.c_str(), "%d,%d,%d", &EditStart, &EditEnd, &EditStep);

    // Parse threshold range string into integers
    string treshold_range = result["treshold"].as<string>(); // 'treshold' typo kept intentionally for exact match with original Bash script
    int ThresholdStart, ThresholdEnd, ThresholdStep;
    std::sscanf(treshold_range.c_str(), "%d,%d,%d", &ThresholdStart, &ThresholdEnd, &ThresholdStep);

    // Get SLURM working directory if applicable
    string DIR = get_slurm_dir();
    string pr_dir = printdir;
    string temp_output = ""; // we use temp output to extract the results of our execute commands function or the exit status if the command doesn't have a return

    // ##########################################
    // # Simplify access to parameter variables #
    // ##########################################

    string GENOME = genome;
    string FASTAfile = fasta;
    string SAMfile = sam;
    string DATA_DIR = meta;
    int L1Count = l1count;

    // ############################################
    // # Variable declaration for triple for loop #
    // ############################################
    // Extracted kmer size from the fasta file
    executeCommand("wc -L " + FASTAfile + " | sed 's/ .*//'", temp_output);
    int k = stoi(temp_output);
    // Extract and store total number of probes
    executeCommand("wc -l " + FASTAfile + " | sed 's/ .*//'", temp_output);
    int MaxKmers = stoi(temp_output) / 2;
    // As best F1 scores are obtained using more than 50% of the probes, run tests using 50% to 100% of probes
    int MinKmers = MaxKmers / 2;
    // Avoid using only 1 probe from 1 ORF in case there are only 2 probes
    if (MinKmers < 2)
    {
        MinKmers = 2;
    }
    // Backup of original MaxKmers & MinKmers value to be reset with each edit distance
    int OriginalMaxKmers = MaxKmers;
    int OriginalMinKmers = MinKmers;
    // Keep track of F1 score to only store values in best_F1_history.csv when an improvement is registered in a given kmer size
    int maxF1 = 0;
    // Delete the file if it exists (clears it before any writing happens)
    // Clear output files before initial run
    string allF1filename = pr_dir + "/" + std::to_string(k) + "mers_F1Scores.csv";
    string bestF1filename = pr_dir + "/" + to_string(k) + "mers_best_F1_history.csv";
    std::remove(allF1filename.c_str());
    std::remove(bestF1filename.c_str());

    // We want to use execute command function for this
    executeCommand("mkdir -p " + pr_dir, temp_output);

    std::ofstream bestF1File(pr_dir + "/" + to_string(k) + "mers_best_F1_history.csv");
    std::ofstream allF1File(pr_dir + "/" + to_string(k) + "mers_F1Scores.csv");

    std::string header = "K-mer size\tEdit dist\tThreshold\tK-mers min\tPrecision\tRecall\tF1 Score\n";
    bestF1File << header; // Add header to bestF1 file
    allF1File << header;  // Add header to allF1 file

    bestF1File.close();
    allF1File.close();

    // Since genome only has to be indexed once we check before entering the for loop if it is indexed
    // Check if genome is indexed if not then index it
    if (!fileExists(GENOME + ".index"))
    {
        executeCommand("mrfast --index " + GENOME, temp_output);
    }

// #################################################################################################
// # Vary edit distance (e), threshold(t) and ammount of probes used (m) to see effect on F1 Score #
// #################################################################################################

// Discuss in the next meeting
// Currently I am only parallelizing the outermost loop because each time we execute L1PD the sam file
// MUST exist beforehand and this way I can ensure the sam file is properly generated for each loop
// This implementation is best if the main bottleneck is the mrfast execution and we have a large edit distance
// Additionally the other implementation (tasks) doesn't help with I/O and command line execution which
// is mostly what our inner loops consist of so it will not be explored for now

// Try running wihtout the parallellization to see if thats what causing us issues
#pragma omp parallel for
    for (int e = EditStart; e <= EditEnd; e += EditStep)
    {
        // Add the edit distance to the name of the same file to avoid data race conditions
        std::string samFile_loop = SAMfile + "_" + std::to_string(e); // sam file doesn't end with.sam currently (fix it)
        std::string loop_out = "";                                    // Since this loop is parallelized we should use a different output variable to avoid overwriting
        std::string k_loop = std::to_string(k);
        double maxF1 = 0.0;
        // Generate sam file for current edit distance ensure name is unique for each edit distance to avoid

        // Is the sam file generated for mrfast unique for each edi distance or is it a single file for all??
        std::string searchCmd = "mrfast --search " + GENOME + " --seq " + FASTAfile + " -e " + std::to_string(e) + " -o " + samFile_loop;
        executeCommand(searchCmd, loop_out);

        // Reset variables
        int MaxKmers = OriginalMaxKmers;
        int MinKmers = OriginalMinKmers;
        int BestM = 0;
        double BestMF1 = 0;

        for (int t = ThresholdStart; t <= ThresholdEnd; t += ThresholdStep)
        {
            // Generate the CSV file with the established values
            // A file is created with "ofstream outFile" using the values from the for loops
            ofstream outFile(pr_dir + "/results_e" + to_string(e) + "_t" + to_string(t) + ".csv");

            // Checking if the file was created successfully
            if (outFile.is_open())
            {
                // Writing the header of the file
                outFile << "K-mers min\tPatterns\tTrue +\tFalse +\tFalse -\tPrecision\tRecall\tF1 Score\n";
            }

            // #######################################################################################################
            // #                                       MinKmer Runtime Reducer                                       #
            // # As a specific m value provides best F1 Score in the mayority of cases with any threshold in a       #
            // # given edit distance, do a test run using whole range of m values (50% to 100% of total probes)      #
            // # for each edit distance with first treshold value. After seeing which m value gave highest F1 value  #
            // # for current e value, repeat first and next treshold values with one m value for each edit distance. #
            // #######################################################################################################
            for (int m = MinKmers; m <= MaxKmers; m += 1)
            {
                string strE = to_string(e);
                string strT = to_string(t);
                string strM = to_string(m);
// Use critical to avoid data race conditions
#pragma omp critical
                {
                    cout << "edit " << strE << ", threshold " << strT << ", MinKmers " << strM << endl;
                    cout << "MinKmers: " << MinKmers << endl;
                    cout << "MaxKmers: " << MaxKmers << endl;
                }
                string rawfilename = "L1s_raw_e" + strE + "_t" + strT + "_m" + strM + ".csv";
                string filteredfilename = "L1s_filtered_e" + strE + "_t" + strT + "_m" + strM + ".csv";

                string seqTypeOption = "";
                if (SequenceType == "Alu")
                    seqTypeOption = " -q 0";

                // Generate a CSV with possible L1s
                executePythonScript(DIR + "/L1PD_files/L1PD.py", samFile_loop + " " + FASTAfile + " -t " + strT + " -m " + strM + " --data_dir " + DATA_DIR + seqTypeOption + " --csvoutput > " + pr_dir + "/" + rawfilename, true);
                // Compare possible L1s from L1PD against L1s from L1Base2 to filter L1s with no matches
                executePythonScript(DIR + "/precision_recall/filter_possible_sequences.py", pr_dir + "/" + rawfilename + " " + FASTAfile + " -t " + strT + " --data_dir " + DATA_DIR + seqTypeOption + " > " + pr_dir + "/" + filteredfilename, true);

                // Use wc t get line count, with sed's help to remove trailing file name
                // Line counts correspond to Positives and True Positives.
                // cambiar por grep (maybe)

                executeCommand("wc -l " + pr_dir + "/" + rawfilename + " | sed 's/ .*//'", loop_out);
                int PatCount = stoi(loop_out);

                if (PatCount > 0)
                {
                    // True positives, false positives, and false negatives
                    executeCommand("grep -v '^DUPLICATE' " + pr_dir + "/" + filteredfilename + " | wc -l | sed 's/ .*//'", loop_out);

                    int TPCount = stoi(loop_out);
                    int FPCount = PatCount - TPCount;
                    int FNCount = L1Count - TPCount;

                    // We assume the file with full-length intact L1s has 'FLI-L1'
                    // in its name (upper or lower case).
                    // executeCommand("grep -Fci 'FLI-L1' " + pr_dir + "/" + filteredfilename, loop_out);
                    // int FLIL1sFound = stoi(loop_out);
                    int FLIL1sFound = 0;
                    if (SequenceType == "Line1")
                    {
                        executeCommand("grep -Fci 'FLI-L1' " + pr_dir + "/" + filteredfilename, loop_out);
                        FLIL1sFound = stoi(loop_out);
                    }
                    else
                    {
                        executeCommand("grep -Fci 'alu' " + pr_dir + "/" + filteredfilename, loop_out);
                        FLIL1sFound = stoi(loop_out);
                    }

                    // Calculate precision, recall, and F1 Score.
                    float precision = static_cast<float>(TPCount) / PatCount;
                    float recall = static_cast<float>(TPCount) / L1Count;
                    float F1Score = 2 * precision * recall / (precision + recall);

                    string strPatCount = to_string(PatCount);
                    string strTPCount = to_string(TPCount);
                    string strFPCount = to_string(FPCount);
                    string strFNCount = to_string(FNCount);
                    string strPrecision = to_string(precision);
                    string strRecall = to_string(recall);
                    string strF1Score = to_string(F1Score);

                    ofstream resultsOutfile(pr_dir + "/results_e" + strE + "_t" + strT + "_m" + strM + ".txt");

                    if (resultsOutfile.is_open())
                    {
                        resultsOutfile << "FLI-L1s:   " + to_string(FLIL1sFound) + "\n";
                        resultsOutfile << "L1Count:   " + to_string(L1Count) + "\n";
                        resultsOutfile << "PatCount:  " + strPatCount + "\n";
                        resultsOutfile << "TPCount:   " + strTPCount + "\n";
                        resultsOutfile << "Precision: " + strPrecision + "\n";
                        resultsOutfile << "Recall:    " + strRecall + "\n";
                        resultsOutfile << "F1 Score:  " + strF1Score + "\n";
                    }

                    // Save data for best m score for e and t value
                    ofstream bestmOutfile(pr_dir + "/results_e" + strE + "_t" + strT + ".csv");

                    if (bestmOutfile.is_open())
                    {
                        bestmOutfile << strM + "\t" + strPatCount + "\t" + strTPCount + "\t" + strFPCount + "\t" + strFNCount + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
                    }
                    // If current run is part of real runs or only 2 probes are availible (no need to run test run)
                    if (MinKmers == MaxKmers)
                    {

// When writing to files with the same name we want to avoid data race conditions so we
// wrap it in an omp critical
#pragma omp critical
                        {
                            // Store all F1 scores with their respectively used parameters
                            // ofstream bestmOutfile(pr_dir + "/" + k_loop + "mers_F1Scores.csv");

                            string filename = pr_dir + "/" + k_loop + "mers_F1Scores.csv";
                            // Check if the file already exists and isn't empty
                            std::ofstream F1Outfile(filename, std::ios::app); // always appends after deletion
                            if (!F1Outfile.is_open())
                            {
                                std::cerr << "Failed to open file for appending: " << filename << std::endl;
                            }

                            // Diccionario de diccionarios
                            // guardando toda la informacion (k, e ,m , t precision, recall, f1 score)
                            if (F1Outfile.is_open())
                            {
                                F1Outfile << k_loop + "\t" + strE + "\t" + strT + "\t" + strM + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
                            }

                            // This error is referring to the parallelization of the edit distance which saves every score
                            // in kmers_best_F1_history.csv

                            //  The problem is for each edit distance they are the largest F1 score so all scores get stored
                            // Possible solution: Store all the data in an array then sort it by edit distance
                            // afterwards loop through the F1 scores and only send to the best F1 file
                            // when there's an improvement in F1 score this would be done outside of the loop

                            // We don't want to print to the file on each iteration we would rather store the value and only
                            // print it to the file at the end
                            // Example of command we could use: fgrep "F1 Score" results_e10_t100_m*.txt | sort -n -k 3
                            // fgrep "Recall" results_e10_t100_m*.txt | sort -n -k 2

                            // Max F1 is just grabbed at the end
                            // if (F1Score > maxF1)
                            // {

                            //     // Stores only F1 scores that show improvement in a given kmer value
                            //     ofstream bestF1Outfile(pr_dir + "/" + k_loop + "mers_best_F1_history.csv", std::ios::app);

                            //     if (bestF1Outfile.is_open())
                            //     {
                            //         bestF1Outfile << k_loop + "\t" + strE + "\t" + strT + "\t" + strM + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
                            //     }

                            //     // Stores current best F1 score to be used in if statmement comparison
                            //     maxF1 = F1Score;
                            // }
                        }
                        // As best m value has been selected, skip process of selecting best m value below
                        break;
                    }
                    // Check if current m value gave higher F1 score than previous m values in current edit distance during test runs
                    if (F1Score > BestMF1)
                    {
                        // Stores best m value for a given edit distance
                        BestM = m;
                        // Stores F1 Score for BestM to be compared in if statement
                        BestMF1 = F1Score;
                    }
                }

                else
                {
                    ofstream emptyPatCountOutfile(pr_dir + "/results_e" + strE + "_t" + strT + "_m" + strM + ".csv");

                    if (emptyPatCountOutfile.is_open())
                        emptyPatCountOutfile << "PatCount = 0";
                }
            }
            if (MinKmers != MaxKmers)
            {
                t = t - ThresholdStep;
                MinKmers = BestM;
                MaxKmers = BestM;
            }
        }
    }

    // Now lets compute the best F1 Score from all of our F1 scores
    // Did it with pandas as it would be cleaner, simpler and the execution delay should be minimal

    executePythonScript(DIR + "/precision_recall/best_f1.py", " -i " + pr_dir + "/" + to_string(k) + "mers_F1Scores.csv" + " -o " + pr_dir + "/" + to_string(k) + "mers_best_F1_history.csv", true);

    cout << "-i " + pr_dir + "/" + to_string(k) + "mers_F1Scores.csv" + " -o " + pr_dir + "/" + to_string(k) + "mers_best_F1_history.csv" << endl;
    // ofstream bestF1Outfile(pr_dir + "/" + k_loop + "mers_best_F1_history.csv", std::ios::app);

    // if (bestF1Outfile.is_open())
    // {
    //     bestF1Outfile << k_loop + "\t" + strE + "\t" + strT + "\t" + strM + "\t" + strPrecision + "\t" + strRecall + "\t" + strF1Score + "\n";
    // }
    // you MUST add the quotes for the variables if not it will seperate the variables due to spaces
    // Check if we want to make a graph
    if (graph)
    {
        // Graph command
        string cmd = "-i " + pr_dir + "/" + to_string(k) + "mers_F1Scores.csv -x";
        for (const auto &val : x_axis)
        {
            cmd += " '" + val + "'"; // Add a space before each value
        }
        cmd += " -y";
        for (const auto &val : y_axis)
        {
            cmd += " '" + val + "'"; // Add a space before each value
        }

        // Check if we want to make the graph 3D, 3D and stackable or 2D and Stackable
        if (three_d)
        {
            cmd += " -z";
            for (const auto &val : z_axis)
            {
                cmd += " '" + val + "'"; // Add a space before each value
            }

            cmd += " -d";
            if (stackable)
                cmd += " -s";
        }
        else if (stackable)
        {
            cmd += " -s";
        }

        cmd += " -o '" + output + "'";

        executePythonScript("read_csv.py", cmd, true);
    }
    return 0;
}
