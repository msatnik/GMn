#ifndef CSVHANDLER_H
#define CSVHANDLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

class CSVHandler {
public:
    // Write values from a vector to a CSV file with configurable precision and values per line
    static void WriteToCSV(const std::string& filename, const std::vector<double>& data, int valuesPerLine = 10, int precision = 2) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to create file " << filename << std::endl;
            return;
        }

        int counter = 0;
        for (const auto& value : data) {
            outFile << std::fixed << std::setprecision(precision) << value << ",";

            counter++;
            if (counter == valuesPerLine) {
                outFile << std::endl;
                counter = 0;
            }
        }

        // Move the file pointer back to remove the last comma and replace it with a newline
        outFile.seekp(-1, std::ios_base::end);
        outFile << std::endl;

	cout<<"Vector written in CSV form to: "<<filename<<endl;

        outFile.close();
    }

    // Read values from a CSV file into a vector
    static void ReadFromCSV(const std::string& filename, std::vector<double>& data) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            while (std::getline(iss, token, ',')) {
                double value;
                std::istringstream(token) >> value;
                data.push_back(value);
            }
        }

	cout<<"CSV file read in from: "<<filename<<endl;

        file.close();
    }
};

#endif // CSVHANDLER_H
