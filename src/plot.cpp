#include "data.h"
#include "matplotlibcpp.h"
#include <ctime>

void initResidualPlot(Data2D& data, std::string finalDirectoryName) {
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
    matplotlibcpp::xlim(0.0, static_cast<double>(x.size() - 1)); 
    matplotlibcpp::xlabel("Iteration");
    matplotlibcpp::ylabel("Residual");
    matplotlibcpp::title("Initial Residuals");
    matplotlibcpp::legend();
    // matplotlibcpp::show();
    matplotlibcpp::grid(true);
    // matplotlibcpp::show();
    matplotlibcpp::save("../out/" + finalDirectoryName + "/initialResidualPlot.png");
    matplotlibcpp::clf();
}
void outerLoopResidualPlot(Data2D& data, std::string finalDirectoryName) {
    matplotlibcpp::figure();  // Create a new figure
    // create three graphs in the same plot for the vectors data.pIterationRes with label (p), data.uIterationRes (u), data.vIterationRes (v)
    std::vector<double> x2;
    for (int i = 0; i < data.pIterationRes.size(); i++) {
        x2.push_back(i);
    }
    std::vector<double> y4 = data.pIterationRes;
    std::vector<double> y5 = data.uIterationRes;
    std::vector<double> y6 = data.vIterationRes;
    std::map<std::string, std::string> keywords;
    //y-log plot
    matplotlibcpp::named_semilogy("p", x2, y4, "r-");
    matplotlibcpp::named_semilogy("u", x2, y5, "g-");
    matplotlibcpp::named_semilogy("v", x2, y6, "b-");
    matplotlibcpp::xlim(0.0, static_cast<double>(x2.size() - 1)); 
    matplotlibcpp::xlabel("Iteration");
    matplotlibcpp::ylabel("Residual");
    matplotlibcpp::title("Outer-Loop Iteration Residuals");
    matplotlibcpp::legend();
    // matplotlibcpp::show();
    matplotlibcpp::grid(true);
    // matplotlibcpp::show();
    matplotlibcpp::save("../out/" + finalDirectoryName + "/outerLoopResidualPlot.png");
    matplotlibcpp::clf();
}