#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// Function to split a string by a delimiter and return a vector of tokens
vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

int csv_to_latex_table() {
    string filename = "./csv/Worlddata_Ye_err.csv"; // Specify your CSV file name here
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Could not open the file!" << endl;
        return 1;
    }

    string line;
    vector<vector<string>> tableData;

    // Reading the CSV file line by line
    while (getline(file, line)) {
        // Remove any carriage return or extra spaces that might cause issues
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(remove(line.begin(), line.end(), '\t'), line.end());
        
        vector<string> row = split(line, ',');
        tableData.push_back(row);
    }

    file.close();

    // Outputting the LaTeX table
    cout << "\\begin{tabular}{|";
    for (size_t i = 0; i < tableData[0].size(); i++) {
        cout << "c|";
    }
    cout << "}" << endl;
    cout << "\\hline" << endl;

    for (size_t i = 0; i < tableData.size(); i++) {
        for (size_t j = 0; j < tableData[i].size(); j++) {
            cout << tableData[i][j];
            if (j < tableData[i].size() - 1) {
                cout << " & ";
            }
        }
        cout << " \\\\" << endl;
        cout << "\\hline" << endl;
    }

    cout << "\\end{tabular}" << endl;

    return 0;
}
