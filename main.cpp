#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

#define CATCH_CONFIG_RUNNER  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

using namespace std;
using namespace arma;

// Opening object that writes to file
ofstream ofile;

// Declaring functions to do the diagonalization
void jacobi(double** A, double** R, int n);
void rotate(double** A, double** R, int k, int l, int n);
double maxvalue(double** A, int* k, int* l, int n);


// First unit test
TEST_CASE( "Maximum value is extracted", "[]" ) {
    int k,l;
    int m = 5;
    // Declaring dynamical matrix, with largest element = 10
    double **testmatrix = new double*[m];
    for(int i = 0; i < m; i++){
        testmatrix[i] = new double[m];
    }
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            testmatrix[i][j] = 5;
        }
    }
    testmatrix[1][4] = 10;

    // Doing the actual test
    double maximum_value_of_off_diagonal_elements = maxvalue(testmatrix,&k,&l,m);
    REQUIRE( maximum_value_of_off_diagonal_elements == 10 );

    // Deleting testmatrix
    for(int i = 0; i < m; ++i) {
        delete [] testmatrix[i];
    }
    delete [] testmatrix;
}

// End of first unit test
/*

TEST_CASE("Checking the Jacobi-matrix","[]"){
    // Declaring constants
    int m = 5;
    int k,l;
    double *egenverdier = new double[m];

    // Initializing matrices
    double **testmatrix = new double*[m];
    for(int i = 0; i < m; i++){
        testmatrix[i] = new double[m];
    }
    double **R = new double*[m];
    for(int i = 0; i < m; i++){
        R[i] = new double[m];
    }

    // Filling up matrices
    for(int i=0; i<m;i++){
        for(int j=0;j<m;j++){
            testmatrix[i][j] = 3*i*j*j;
        }
    }
    for(int i=0;i<m;i++){
        testmatrix[i][i] = 3*i*i;
    }
    testmatrix[m][m] = 18;

    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            if(i==j){
                R[i][j] = 1.0;
            }
            else{
                R[i][j] = 0.0;
            }
        }
    }

    // Calling function
    jacobi(testmatrix,R,m);

    // Sorting the eigenvalues in ascending values
    vec w = zeros<vec>(m);
    for(int i=0;i<m;i++){
        w(i) = testmatrix[i][i];
    }
    vec sorted_w = sort(w);

    for(int i=0;i<m;i++){
        egenverdier[i] = sorted_w(i);
    }

    // The test itself
    REQUIRE(egenverdier[0]==Approx(-246.7864));
    REQUIRE(egenverdier[1]==Approx(-80.8741));
    REQUIRE(egenverdier[2]==Approx(-18.4024));
    REQUIRE(egenverdier[3]==Approx(-0.5985));
    REQUIRE(egenverdier[4]==Approx(454.6613));

    // Deleting testmatrices
    for(int i = 0; i < m; ++i) {
        delete [] testmatrix[i];
    }
    delete [] testmatrix;

    for(int i = 0; i < m; ++i) {
        delete [] R[i];
    }
    delete [] R;
    delete [] egenverdier;

}
// End second unit test
*/
// Main program begins
int main(int argc, char *argv[])
{
    int variabel = 1;
    char* variabel2[1];
    int result = Catch::Session().run(variabel, variabel2);

    /*
    int variable3 = 1;
    char* variabel4[1];
    int result2 = Catch::Session().run(variable3, variabel4);*/


    // Declaring constants and pointers
    int n, indeks;
    double rho_zero, rho_max, e, konstant, h, hh, konstant2,  minimum_value;
    char *outfilename, *outfilename2;

    // Reading in outputfile. Abort if there are too few command line arguments
    if (argc<=3){
        cout << "Bad usage: " << argv[0] << "read also outputfile and 'n' on the same line"<<endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
        outfilename2=argv[2];
        n = atoi(argv[3]);
    }

    // Opening a file for the program
    ofile.open(outfilename);

    // Declaring dynamical matrixes
    double **A = new double*[n];
    for(int i = 0; i < n; i++)
        A[i] = new double[n];

    double **R = new double*[n];
    for(int i = 0; i < n; i++)
        R[i] = new double[n];

    double **data_points = new double*[n];
    for(int i = 0; i < n; i++)
        data_points[i] = new double[4];

    // Declaring dynamical arrays
    double *omega = new double[4];
    double *rho = new double[n];
    double *eigenvalues = new double[n];
    double *r = new double[n];
    double *psi = new double[n];
    double *r_values = new double[n];

    // Question b)
    // Determining values of constants
    rho_zero = 0.0;
    rho_max = 4.5;
    h = (rho_max - rho_zero)/((double) n+1);
    hh = h*h;
    e = -1.0/hh;
    konstant = 2.0/hh;

    // Filling up rho
    rho[0] = rho_zero + h;
    for(int i=1;i<n;i++){
        rho[i] = rho[0] + i*h;
    }

    // Filling up matrix
    for (int i =0; i<n; i++) {
        for (int j=0; j<n; j++) {
            A[i][j] = 0;
        }
    }
    for(int i=1;i<n-1;i++){
        A[i][i-1] = e;
        A[i][i] = rho[i]*rho[i] + konstant;
        A[i][i+1] = e;
    }


    // First and last terms of matrix
    A[0][0] = rho[0]*rho[0] + konstant;
    A[0][1] = e;
    A[n-1][n-2] = e;
    A[n-1][n-1] = rho[n-1]*rho[n-1] + konstant;

    // Armadillo-matrix - copying the A-matrix for use in Armadillo
    mat B = zeros<mat>(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            B(i,j) = A[i][j];
        }
    }

    // Here we start the time-taking
    clock_t start, finish;
    start = clock();

    // Using the Jacobi-method built into Armadillo to solve the eigenvalue problem
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, B);

    // Here we end the time-taking
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC ); // Calculating time
    cout << "Time Armadillo Jacobi-method in seconds: " << timeused << endl; // Printing out time spent on Jacobi

    // Sorting eigenvalues
    eigval = sort(eigval);

    // Printing eigenvalues produced by Armadillo
    cout<< "Armadillo eigenvalue 0 = " << eigval(0) << endl;
    cout<< "Armadillo eigenvalue 1 = " << eigval(1) << endl;
    cout<< "Armadillo eigenvalue 2 = " << eigval(2) << endl;

    // Then using my own implimintation of the Jacobi-method
    // Initializing the eigenvector-matrix
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                R[i][j] = 1.0;
            }
            else{
                R[i][j] = 0.0;
            }
        }
    }

    // Here we start the time-taking
    start = clock();

    // Calling function
    jacobi(A,R,n);

    // Here we end the time-taking
    finish = clock();
    double timeused_own_method = (double) (finish - start)/(CLOCKS_PER_SEC ); // Calculating time
    cout << "Time own Jacobi-method in seconds: " << timeused_own_method << endl; // Printing out time spent on Jacobi

    // Sorting the eigenvalues in ascending values
    vec v = zeros<vec>(n);
    for(int i=0;i<n;i++){
        v(i) = A[i][i];
    }
    vec sorted_v = sort(v);

    for(int i=0;i<n;i++){
        eigenvalues[i] = sorted_v(i);
    }

    // Printing out results
    cout<<"Eigenvalue 0 = " << eigenvalues[0] << endl;
    cout<<"Eigenvalue 1 = " << eigenvalues[1] << endl;
    cout<<"Eigenvalue 2 = " << eigenvalues[2] << endl;
    // End question b)

    cout<<"End question b)"<< endl;

    // Question d)
    // Determining values of constants
    rho_zero = 0.0;
    rho_max = 7.9;
    h = (rho_max - rho_zero)/((double) n+1);
    hh = h*h;
    e = -1.0/hh;
    konstant = 2.0/hh;

    // Filling up arrays
    rho[0] = rho_zero + h;
    for(int i=1;i<n;i++){
        rho[i] = rho[0] + i*h;
    }

    omega[0] = 0.01;
    omega[1] = 0.5;
    omega[2] = 1.0;
    omega[3] = 5.0;

    for(int k=0;k<4;k++){

        // Filling up matrix
        for(int i=1;i<n-1;i++){
            A[i][i-1] = e;
            A[i][i] = omega[k]*omega[k]*rho[i]*rho[i] + 1.0/rho[i]+ konstant;
            A[i][i+1] = e;
        }

        // First and last terms of matrix
        A[0][0] = omega[k]*omega[k]*rho[0]*rho[0] + 1.0/rho[0] + konstant;
        A[0][1] = e;
        A[n-1][n-2] = e;
        A[n-1][n-1] = omega[k]*omega[k]*rho[n-1]*rho[n-1] + 1.0/rho[n-1] + konstant;

        // Initializing the eigenvector-matrix
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i==j){
                    R[i][j] = 1.0;
                }
                else{
                    R[i][j] = 0.0;
                }
            }
        }

        // Calling function
        jacobi(A,R,n);

        // Sorting the eigenvalues in ascending values
        vec v = zeros<vec>(n);
        for(int i=0;i<n;i++){
            v(i) = A[i][i];
        }
        vec sorted_v = sort(v);

        for(int i=0;i<n;i++){
            eigenvalues[i] = sorted_v(i);
        }

        // Printing out results
        cout<<"Eigenvalue 0 = " << eigenvalues[0] <<" for omega = "<< omega[k] <<endl;


    }
    // End question d)
    cout<<"End question d)"<< endl;

    // Question e)
    // The r-vector
    konstant2 = ((6.58211899*6.58211899*9.0)/(0.510998910*1.44))*(1e-13);
    for(int i=0;i<n;i++){
        r[i] = konstant2*rho[i];
    }

    // First part of question e). Repulsion between the electrons.
    // Finding the eigenvectors for the four different omega's. These four eigenvectors are stored in the matrix "datapoints"
    for(int k=0;k<4;k++){
        // For k=0, omega = 0.01 and requires different value of rho_max
        if (k==0){
            // New rho_max
            rho_zero = 0.0;
            rho_max = 50;
            h = (rho_max - rho_zero)/((double) n+1);
            hh = h*h;
            e = -1.0/hh;
            konstant = 2.0/hh;

            // Filling up arrays
            rho[0] = rho_zero + h;
            for(int i=1;i<n;i++){
                rho[i] = rho[0] + i*h;
            }
            // The r-vector
            for(int i=0;i<n;i++){
                r_values[i] = konstant2*rho[i];
            }
            // Filling up matrix
            for(int i=1;i<n-1;i++){
                A[i][i-1] = e;
                A[i][i] = omega[k]*omega[k]*rho[i]*rho[i] + 1.0/rho[i]+ konstant;
                A[i][i+1] = e;
            }

            // First and last terms of matrix
            A[0][0] = omega[k]*omega[k]*rho[0]*rho[0] + 1.0/rho[0] + konstant;
            A[0][1] = e;
            A[n-1][n-2] = e;
            A[n-1][n-1] = omega[k]*omega[k]*rho[n-1]*rho[n-1] + 1.0/rho[n-1] + konstant;

            // Initializing the eigenvector-matrix
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    if(i==j){
                        R[i][j] = 1.0;
                    }
                    else{
                        R[i][j] = 0.0;
                    }
                }
            }

            // Calling function
            jacobi(A,R,n);

            // Finding the minimum eigenvalue with its corresponding eigenvector, psi
            indeks = 0;
            minimum_value = 1000.0;

            for(int i=0;i<n;i++){
                if(A[i][i] < minimum_value){
                    minimum_value = A[i][i];
                    indeks = i;
                }
            }

            for(int i=0;i<n;i++){
                psi[i] = R[i][indeks];
            }

            // Preserving data-points (the eigenvector) for plotting
            for(int i=0; i<n;i++){
                data_points[i][k] = psi[i];
            }
            // Switching back to original value of rho_max
            rho_zero = 0.0;
            rho_max = 7.9;
            h = (rho_max - rho_zero)/((double) n+1);
            hh = h*h;
            e = -1.0/hh;
            konstant = 2.0/hh;

            // Filling up arrays
            rho[0] = rho_zero + h;
            for(int i=1;i<n;i++){
                rho[i] = rho[0] + i*h;
            }
        }
        else{

        // Filling up matrix
        for(int i=1;i<n-1;i++){
            A[i][i-1] = e;
            A[i][i] = omega[k]*omega[k]*rho[i]*rho[i] + 1.0/rho[i]+ konstant;
            A[i][i+1] = e;
        }

        // First and last terms of matrix
        A[0][0] = omega[k]*omega[k]*rho[0]*rho[0] + 1.0/rho[0] + konstant;
        A[0][1] = e;
        A[n-1][n-2] = e;
        A[n-1][n-1] = omega[k]*omega[k]*rho[n-1]*rho[n-1] + 1.0/rho[n-1] + konstant;

        // Initializing the eigenvector-matrix
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i==j){
                    R[i][j] = 1.0;
                }
                else{
                    R[i][j] = 0.0;
                }
            }
        }

        // Calling function
        jacobi(A,R,n);

        // Finding the minimum eigenvalue with its corresponding eigenvector, psi
        indeks = 0;
        minimum_value = 1000.0;

        for(int i=0;i<n;i++){
            if(A[i][i] < minimum_value){
                minimum_value = A[i][i];
                indeks = i;
            }
        }

        for(int i=0;i<n;i++){
            psi[i] = R[i][indeks];
        }

        // Preserving data-points (the eigenvector) for plotting
        for(int i=0; i<n;i++){
            data_points[i][k] = psi[i];
        }
        }

    }

    // Writing to file for plotting (open file in MatLab and Import Data and plot in usual way).
    for(int i=0;i<n;i++){
        ofile << setw(15) << setprecision(8) << r[i] << "\t"
              << setw(15) << setprecision(8) << r_values[i] << "\t"
              << setw(15) << setprecision(8) << data_points[i][0] << "\t"
              << setw(15) << setprecision(8) << data_points[i][1] << "\t"
              << setw(15) << setprecision(8) << data_points[i][2] << "\t"
              << setw(15) << setprecision(8) << data_points[i][3] << endl;
    }

    ofile.close();
    ofile.open(outfilename2);

    // Second part of question e); no repulsion between the electrons
    // Finding the eigenvectors for the four different omega's. These four eigenvectors are stored in the matrix "datapoints"
    for(int k=0;k<4;k++){
        // For k=0, omega = 0.01 and requires different value of rho_max
        if(k==0){
            // New rho_max
            rho_zero = 0.0;
            rho_max = 50;
            h = (rho_max - rho_zero)/((double) n+1);
            hh = h*h;
            e = -1.0/hh;
            konstant = 2.0/hh;

            // Filling up arrays
            rho[0] = rho_zero + h;
            for(int i=1;i<n;i++){
                rho[i] = rho[0] + i*h;
            }
            // The r-vector
            for(int i=0;i<n;i++){
                r_values[i] = konstant2*rho[i];
            }
            // Filling up matrix
            for(int i=1;i<n-1;i++){
                A[i][i-1] = e;
                A[i][i] = omega[k]*omega[k]*rho[i]*rho[i] + konstant;
                A[i][i+1] = e;
            }

            // First and last terms of matrix
            A[0][0] = omega[k]*omega[k]*rho[0]*rho[0] + konstant;
            A[0][1] = e;
            A[n-1][n-2] = e;
            A[n-1][n-1] = omega[k]*omega[k]*rho[n-1]*rho[n-1] + konstant;

            // Initializing the eigenvector-matrix
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    if(i==j){
                        R[i][j] = 1.0;
                    }
                    else{
                        R[i][j] = 0.0;
                    }
                }
            }

            // Calling function
            jacobi(A,R,n);

            // Finding the minimum eigenvalue with its corresponding eigenvector, psi
            indeks = 0;
            minimum_value = 1000.0;

            for(int i=0;i<n;i++){
                if(A[i][i] < minimum_value){
                    minimum_value = A[i][i];
                    indeks = i;
                }
            }

            for(int i=0;i<n;i++){
                psi[i] = R[i][indeks];
            }

            // Preserving data-points (the eigenvector) for plotting
            for(int i=0; i<n;i++){
                data_points[i][k] = psi[i];
            }
            // Switching back to original value of rho_max
            rho_zero = 0.0;
            rho_max = 7.9;
            h = (rho_max - rho_zero)/((double) n+1);
            hh = h*h;
            e = -1.0/hh;
            konstant = 2.0/hh;

            // Filling up arrays
            rho[0] = rho_zero + h;
            for(int i=1;i<n;i++){
                rho[i] = rho[0] + i*h;
            }
        }
        else{
        // Filling up matrix
        for(int i=1;i<n-1;i++){
            A[i][i-1] = e;
            A[i][i] = omega[k]*omega[k]*rho[i]*rho[i] + konstant;
            A[i][i+1] = e;
        }

        // First and last terms of matrix
        A[0][0] = omega[k]*omega[k]*rho[0]*rho[0] + konstant;
        A[0][1] = e;
        A[n-1][n-2] = e;
        A[n-1][n-1] = omega[k]*omega[k]*rho[n-1]*rho[n-1] + konstant;

        // Initializing the eigenvector-matrix
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i==j){
                    R[i][j] = 1.0;
                }
                else{
                    R[i][j] = 0.0;
                }
            }
        }

        // Calling function
        jacobi(A,R,n);

        // Finding the minimum eigenvalue with its corresponding eigenvector, psi
        indeks = 0;
        minimum_value = 1000.0;

        for(int i=0;i<n;i++){
            if(A[i][i] < minimum_value){
                minimum_value = A[i][i];
                indeks = i;
            }
        }

        for(int i=0;i<n;i++){
            psi[i] = R[i][indeks];
        }

        // Preserving data-points (the eigenvector) for plotting
        for(int i=0; i<n;i++){
            data_points[i][k] = psi[i];
        }
        }
    }

    // Writing to file for plotting (open file in MatLab and Import Data and plot in usual way).
    for(int i=0;i<n;i++){
        ofile << setw(15) << setprecision(8) << r[i] << "\t"
              << setw(15) << setprecision(8) << r_values[i] << "\t"
              << setw(15) << setprecision(8) << data_points[i][0] << "\t"
              << setw(15) << setprecision(8) << data_points[i][1] << "\t"
              << setw(15) << setprecision(8) << data_points[i][2] << "\t"
              << setw(15) << setprecision(8) << data_points[i][3] << endl;
    }

    // Closing output file
    ofile.close();

    // Freeing memory
    delete [] rho;
    delete [] eigenvalues;
    delete [] omega;
    delete [] r;
    delete [] r_values;

    for(int i = 0; i < n; ++i) {
        delete [] A[i];
    }
    for(int i = 0; i < n; ++i) {
        delete [] R[i];
    }
    for(int i = 0; i < 4; ++i) {
        delete [] data_points[i];
    }

    delete [] A;
    delete [] R;
    delete [] data_points;

    // Success
    return 0;
}

// The Jacobi-function
void jacobi(double** A, double** R, int n){

    // Declaring local constants
    int k,l;
    double N = (double) n;
    double epsilon = 1.0e-8;
    double maximum_number_of_iterations = 0.01*N*N*N;
    int number_of_iterations = 0;
    double maximum_value_of_off_diagonal_elements = maxvalue(A,&k,&l,n); // maxvalue is a function that is given below

    // Doing the actual rotations
    while(maximum_value_of_off_diagonal_elements > epsilon && (double) number_of_iterations < maximum_number_of_iterations){
        rotate(A,R,k,l,n);
        maximum_value_of_off_diagonal_elements = maxvalue(A,&k,&l,n);
        number_of_iterations++;
    }

    cout<<"Number of similarity transformations = "<<number_of_iterations<<endl;
}

// Function that finds the maximum value of the off-diagonal elements of the current matrix A
double maxvalue(double** A, int* k, int* l, int n){
    double max = 0.0;
    double value = 0.0;
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            value = fabs(A[i][j]);
            if(value>max){
                max = value;
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}

// Function that does similarity transformation on the matrix A
void rotate(double** A, double** R, int k, int l, int n){

    double s,c; // Sine and cosine

    // Finding sine and cosine
    if(A[k][l] != 0.0){
        double t, tau, t1, t2;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        t1 = fabs(-tau + sqrt(tau*tau + 1));
        t2 = fabs(-tau - sqrt(tau*tau + 1));
        if(t1<t2){
            t = -tau + sqrt(tau*tau + 1);
        }
        else{
            t = -tau - sqrt(tau*tau + 1);
        }
        c = 1.0/sqrt(1 + t*t);
        s = t*c;
    }
    else{
        s = 0.0;
        c = 1.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il, a_kl; // Matrix elements

    // Preserving old values
    a_kk = A[k][k];
    a_ll = A[l][l];
    a_kl = A[k][l];

    // Doing the actual rotation
    // Changing all the matrix elements that are indexed with l and k
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = (a_kk - a_ll)*c*s + a_kl*(c*c - s*s);
    A[l][k] = A[k][l];


    // Changing all the matrix elements that are NOT indexed with l and k
    for(int i=0; i<n; i++){
        if(i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        // The new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
return;
}
