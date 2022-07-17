
#include "utilities.h"
#include "pressure.h"
#include "velocity.h"
#include <iostream>

int main() {
    if (!read_parameters() or !dump_simulation_parameters()) {
        std::cerr << "Error: Cannot handle simulation params.\n";
        return -1;
    }

//    if (!check_stability()){
//        cout<<"Stability conditions is not satisfied"<<endl;
//        return -1;
//    }

    initial_distribution();

    int step = 1;
    output_all_fields(step);
    while (true) {
//        cout << step << endl;

        calculate_auxilliary_velocity();

        calculate_pressure(step);
        calculate_velocity();

//        if (check_velocity_convergence()){
//            save_results(step);
//            output_all_fields(step);
//            return 0;
//        }

//        save_previous_velocity();
        step++;
        if (step < 100) {
            output_all_fields(step);
        }

        if (step < 1000 and step % 50 == 0) {
            output_all_fields(step);
        }

        if (step % 1000 == 0) {
            output_all_fields(step);
        }

        if (step % 100000000 == 0) {
            return 0;
        }

    }
}
