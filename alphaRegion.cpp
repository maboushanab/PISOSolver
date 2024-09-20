#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <string>
std::string getLineByNumber(std::ifstream& file, int targetLineNumber) {
    std::string line;
    int currentLine = 0;

    // Reset the file stream back to the beginning
    file.clear(); // Clear any error flags
    file.seekg(0, std::ios::beg); // Rewind to the start of the file

    // Move to the target line by reading each line sequentially
    while (currentLine < targetLineNumber && std::getline(file, line)) {
        ++currentLine;
    }

    // Check if the target line number was reached
    if (currentLine == targetLineNumber) {
        return line;
    } else {
        return ""; // Return an empty string if the target line doesn't exist
    }
}
//only does rectangular shapes
//Read point region with alpha based on coordinates
std::vector<double> findPointRegion(std::ifstream& file, std::vector<std::vector<double>> regionCoords) {
    std::string CellStart = "# cells";
    std::string pointStart = "# points";
    std::string line;
    int cellLineNumber = 0;  // Initialize line counters
    int pointLineNumber = 0;

    if (file.is_open()) {
        // First, find the line starting with "# points"
        while (std::getline(file, line)) {
            ++pointLineNumber;
            if (line.find(pointStart) == 0) {
                break;
            }
        }
        while (std::getline(file, line)) {
            ++cellLineNumber;
            if (line.find(CellStart) == 0) {
                cellLineNumber += pointLineNumber;
                break;
            }
        }
        std::vector<std::tuple<double, double, double>> points;
        
        //point line is formatted as: id x y
        //look for point with x = region(0)(0)
        for (int i = pointLineNumber+1; i < cellLineNumber; i++) {
            line = getLineByNumber(file, i);
            // Create a string stream from the line
            std::istringstream iss(line);
            double num1, num2, num3;
            
            // Read three numbers from the line
            if (!(iss >> num1 >> num2 >> num3)) {
                // std::cerr << "Error reading line: " << line << std::endl;
                continue; // Skip lines with invalid format
            }
            
            // Add the tuple to the vector
            points.emplace_back(num1, num2, num3);
        }
        std::ofstream regionFile;
        regionFile.open("regionPoints.txt");
        for (int i = 0; i < points.size(); i++) {
            if (std::get<1>(points[i]) >= regionCoords[0][0] && std::get<1>(points[i]) <= regionCoords[1][0] && std::get<2>(points[i]) >= regionCoords[0][1] && std::get<2>(points[i]) <= regionCoords[2][1]) {
                // export list of points in region to a file
                regionFile << std::get<0>(points[i]) << std::endl;
            }
        }
        regionFile.close();
        file.close();
    } else {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Return some dummy vector for now since function returns a vector
    return std::vector<double>(); 
}

int main() {
    std::ifstream file("damBreak.mesh");
    std::vector<double> p1 = {0.0, 0.0};
    std::vector<double> p2 = {0.146, 0.0};
    std::vector<double> p3 = {0.0, 0.292};
    std::vector<double> p4 = {0.146, 0.292};
    std::vector<std::vector<double>> regionCoords = {p1, p2, p3, p4};
    findPointRegion(file, regionCoords);
    return 0;
}