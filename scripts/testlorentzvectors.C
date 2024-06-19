#include <TLorentzVector.h>
#include <iostream>

int testlorentzvectors() {
    // Define four-momenta of two particles
    TLorentzVector particle1(10.0, 0.0, 0.0, 15.0); // ( px, py, pz,energy)
    TLorentzVector particle2(10.0, 0.0, 0.0, -10.0);

    // Add the four-momenta to get the total four-momentum
    TLorentzVector total = particle1 + particle2;

    // Print out the total energy and momentum
    std::cout << "Total energy: " << total.E() << std::endl;
    std::cout << "Total momentum: " << total.P() << std::endl;

    // Calculate invariant mass
    double invariantMass = total.M();
    std::cout << "Invariant mass: " << invariantMass << std::endl;

    double dotprod = particle1.Dot(particle2);
    std::cout<<"dot product: "<<dotprod<<std::endl;

    return 0;
}
