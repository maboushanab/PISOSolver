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

template<typename T>
std::vector<T> stackToVector(std::stack<T> stack) {
    std::vector<T> vec;
    // Pop elements from the stack and push them into the vector
    while (!stack.empty()) {
        vec.push_back(stack.top());
        stack.pop();
    }
    // Reverse the vector to maintain the original order
    std::reverse(vec.begin(), vec.end());
    return vec;
}

void residualPlot(Data2D& data, std::string finalDirectoryName) {
    if (data.stackOfContinuityResiduals.empty()) {
        std::cout << "No continuity residuals to plot." << std::endl;
        return;
    }
    int i = 0;
    while(!data.stackOfContinuityResiduals.empty()) {
        std::stack<double> stack = data.stackOfContinuityResiduals.top();
        std::vector<double> vStack = stackToVector(stack);
        std::vector<int> iterations(vStack.size());
        std::iota(iterations.begin(), iterations.end(), 0);
        std::string color = colorSequence(i); 
        matplotlibcpp::plot(iterations, vStack, color);
        data.stackOfContinuityResiduals.pop(); 
        i++;
    }
    matplotlibcpp::title("Continuity Residuals");
    matplotlibcpp::xlabel("Iteration");
    matplotlibcpp::ylabel("Residual");
    matplotlibcpp::grid(true);
    matplotlibcpp::legend();
    // matplotlibcpp::show();
    matplotlibcpp::save("../out/" + finalDirectoryName + "/residualPlot.png");
}