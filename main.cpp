#include <iostream>
#include <iomanip> // For formatted output

// Structure to represent the capacitor
struct Capacitor {
    double* time;       // time array
    double* voltage;    // voltage array
    double* current;    // current array
    double C;           // capacitance value
};

typedef struct Capacitor Capacitor;

// Function to allocate dynamic memory for arrays
void initialize_capacitor(Capacitor& cap, int timesteps) {
    cap.time = new double[timesteps];
    cap.voltage = new double[timesteps];
    cap.current = new double[timesteps];
}

// Function to simulate constant current charging of the capacitor
void simulate_constant_current(Capacitor& cap, int timesteps, double dt, double I_const) {
    for (int t = 1; t < timesteps; ++t) {
        cap.time[t] = cap.time[t-1] + dt;
        cap.voltage[t] = cap.voltage[t-1] + I_const * dt / cap.C;
    }
}

// Function to simulate constant voltage charging of the capacitor
void simulate_constant_voltage(Capacitor& cap, int timesteps, double dt, double R, double V0) {
    for (int t = 1; t < timesteps; ++t) {
        cap.time[t] = cap.time[t-1] + dt;
        cap.current[t] = cap.current[t-1] - (cap.current[t-1] * dt) / (R * cap.C); // Add missing semicolon here
        cap.voltage[t] = V0 - cap.current[t] * R;
    }
}

// Function to print results every 200 timesteps
void print_results(const Capacitor& cap, int timesteps) {
    for (int t = 0; t < timesteps; t += 200) {
        std::cout << "Time: " << std::setprecision(6) << cap.time[t] << " s, "
                  << "Voltage: " << cap.voltage[t] << " V, "
                  << "Current: " << cap.current[t] << " A" << std::endl;
    }
}

// Main function
int main() {
    int timesteps = 50000;
    double dt = 1e-10;      // change in time
    double final_time = 5e-6;
    double C = 100e-12;     // Capacitance (F)
    double R = 1e3;         // Resistance (ohms)
    double I_const = 1e-2;  // Constant current (A)
    double V0 = 10.0;       // Constant voltage source (V)

    Capacitor cap;
    cap.C = C;

    // Allocate memory for arrays
    initialize_capacitor(cap, timesteps);

    // Initialize time, voltage, and current at t = 0
    cap.time[0] = 0;
    cap.voltage[0] = 0; // Add missing semicolon here
    cap.current[0] = V0 / R; // For constant voltage, I(0) = V0/R

    // Simulate constant current charging
    simulate_constant_current(cap, timesteps, dt, I_const);

    // Simulate constant voltage charging
    simulate_constant_voltage(cap, timesteps, dt, R, V0);

    // Output the results
    print_results(cap, timesteps);

    // Free dynamically allocated memory
    delete[] cap.time;
    delete[] cap.voltage;
    delete[] cap.current;

    return 0;
}