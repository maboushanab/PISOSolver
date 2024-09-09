#include "data.h"
#include "matplotlibcpp.h"
#include <ctime>

// Function to generate a random color
std::string colorSequence(int i) {
    switch (i % 7) {
        case 0: return "r-";
        case 1: return "g-";
        case 2: return "b-";
        case 3: return "c-";
        case 4: return "m-";
        case 5: return "y-";
        case 6: return "k-"; 
    }
    return "r-";
}

// template<typename T>
// std::vector<T> stackToVector(std::stack<T> stack) {
//     std::vector<T> vec;
//     // Pop elements from the stack and push them into the vector
//     while (!stack.empty()) {
//         vec.push_back(stack.top());
//         stack.pop();
//     }
//     // Reverse the vector to maintain the original order
//     std::reverse(vec.begin(), vec.end());
//     return vec;
// }

void residualPlot(Data2D& data, std::string finalDirectoryName) {
    // create three graphs in the same plot for the vectors data.continuityResiduals with label (p), data.momentumXresiduals (u), data.momentumYresiduals (v)
    std::vector<double> x;
    for (int i = 0; i < data.continuityResiduals.size(); i++) {
        x.push_back(i);
    }
    std::vector<double> y1 = data.continuityResiduals;
    std::vector<double> y2 = data.momentumXResiduals;
    std::vector<double> y3 = data.momentumYResiduals;
    std::map<std::string, std::string> keywords;
    //y-log plot
    matplotlibcpp::named_semilogy("p", x, y1, "r-");
    matplotlibcpp::named_semilogy("u", x, y2, "g-");
    matplotlibcpp::named_semilogy("v", x, y3, "b-");
    matplotlibcpp::xlabel("Iteration");
    matplotlibcpp::ylabel("Residual");
    matplotlibcpp::title("Residuals");
    matplotlibcpp::legend();
    // matplotlibcpp::show();
    matplotlibcpp::grid(true);
    // matplotlibcpp::show();
    matplotlibcpp::save("../out/" + finalDirectoryName + "/residualPlot.png");
}