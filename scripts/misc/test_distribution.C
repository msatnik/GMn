#include <TROOT.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>

void test_distribution() {
    // Number of random samples to generate
    const int nSamples = 1000000;

    // Create a random number generator
    TRandom3 randGen;

    // Create a histogram for the results
    TH1F *hist = new TH1F("hist", "Product of Random Numbers;Product;Frequency", 100, 0, 9);

    // Generate random numbers and fill the histogram
    for (int i = 0; i < nSamples; ++i) {
        double num1 = randGen.Uniform(1, 3); // Random number between 0 and 1
        double num2 = randGen.Uniform(1, 3); // Random number between 1 and 2
        double product = num1 * num2;        // Multiply the two numbers

        hist->Fill(product); // Fill the histogram with the product
    }

    // Create a canvas to draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Random Number Product Histogram", 800, 600);
    hist->Draw();

    // Save the histogram as an image
    //c1->SaveAs("random_product_histogram.png");
}
