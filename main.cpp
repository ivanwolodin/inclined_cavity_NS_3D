
#include "utilities.h"
#include "pressure.h"
#include "velocity.h"


int main() {
    if (!read_parameters()) {
        std::cerr << "Error: Cannot read simulation params.\n";
        return -1;
    }

    if (!dump_simulation_parameters()) {
        std::cerr << "Error: Cannot handle simulation params.\n";
        return -1;
    }
//    if (!check_stability()){
//        cout<<"Stability conditions is not satisfied"<<endl;
//        return -1;
//    }

    initial_distribution();

    int step = 0;
    output_all_fields(step);
//    calculate_pressure(step);
    while (true) {
        step++;
        calculate_auxiliary_velocity();
        calculate_pressure(step);
        calculate_velocity();
        save_previous_velocity();

        if (check_velocity_convergence() and step > 100000000) {
//            save_results(step);
            output_all_fields(step);
            return 0;
        }

//        if (step < 300 or (step < 1000 and step % 100 == 0) or (step < 10000 and step % 1000 == 0) or
//            (step < 100000 and step % 10000 == 0)) {
//            output_all_fields(step);
//        }
        if ((step < 10000 and step % 100 == 0) or (step < 1000000 and step % 1000 == 0) or step % 1000000 == 0) {
            output_averaged_veloicty(step);
            output_all_fields(step);
        }
    }
}
