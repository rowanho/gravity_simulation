// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-3) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  int noBuckets = 10;
  
  int* bucket = new int[NumberOfBodies];
  
  double bucketThreshold = maxV /  noBuckets;
 
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];
  
  // bucket contains the bucket for each body
  // Sort the bodies into buckets
  for (int i=0; i<NumberOfBodies; i++){
      force0[i] = 0.0;
      force1[i] = 0.0;
      force2[i] = 0.0;
      if (t > 0.0){
          double tempV = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );      
          if (tempV >= maxV) {
              bucket[i] = noBuckets - 1;
          } else {
              bucket[i] = std::floor(tempV / bucketThreshold);
          }
      } else {
          bucket[i] = 0;
      }
  }
  
  
  
  
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  
  // Loop through buckets
  for (int b=0; b < 10; b++){
    int noSteps = std::pow(2, b);
    double timeStep = timeStepSize / noSteps;
    for (int step = 0; step < noSteps; step ++ ){
      for (int i = 0; i < NumberOfBodies; i++){
        if (bucket[i] != b){
          continue;
        }
      // Repeat for the number of timesteps
        force0[i] = 0.0;
        force1[i] = 0.0;
        force2[i] = 0.0;
        for (int j=0; j<NumberOfBodies; j++) {
          if (i == j) {
              continue;
          }
          const double distance = sqrt(
              (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
              (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
              (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
          );
          minDx = std::min( minDx,distance );                    
          const double massSquaredDist = mass[j]*mass[i] / distance / distance / distance ;
          force0[i] += (x[j][0]-x[i][0]) * massSquaredDist ;
          force1[i] += (x[j][1]-x[i][1]) * massSquaredDist ;
          force2[i] += (x[j][2]-x[i][2]) * massSquaredDist ;
        }
        
        
        double tempV = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
        if (tempV > maxV) {
            maxV = tempV;
        }
      }
      for (int i = 0; i < NumberOfBodies; i++){
        if (bucket[i] != b){
          continue;
        }
        x[i][0] += timeStep * v[i][0];
        x[i][1] += timeStep * v[i][1];
        x[i][2] += timeStep * v[i][2];
        v[i][0] += timeStep * force0[i] / mass[i];
        v[i][1] += timeStep * force1[i] / mass[i];
        v[i][2] += timeStep * force2[i] / mass[i];
      }
      // Deal with collisions
      for (int i=0; i<NumberOfBodies; i++) {
        if (bucket[i] != b){
          continue;
        }
        for (int j=i+1; j<NumberOfBodies; j++) {
            
          const double distance = sqrt(
            (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
            (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
            (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
          );
        
          if ( distance < 0.01) {
            // Merge object i and j        
            // Get mass of new object
            const double massSum = mass[i] + mass[j];
            const double massRatioI = mass[i] / massSum;
            const double massRatioJ = mass[j] / massSum;
            for (int k = 0; k < 3; k++){
              // Get veloity of new object
              v[i][k] = (massRatioI * v[i][k]) + ( massRatioJ * v[j][k]);
              // Merge object positions
              x[i][k] = (massRatioI * x[i][k]) + (massRatioJ * x[j][k])  ;
            }
            mass[i] =  massSum;
            std::cout << "Collided object position: (" << x[i][0] << ", "<< x[i][1] << ", " << x[i][2] <<  ")";

            for (int k = 0; k < 3; k++){
              v[j][k] = v[NumberOfBodies-1][k];
              x[j][k] = x[NumberOfBodies-1][k];
            }
            mass[j] = mass[NumberOfBodies-1];
            NumberOfBodies --;
          }
        }
      }

    }      
  }
  
  if (NumberOfBodies == 1){
    // Terminate the program
    tFinal = 0;
    std::cout << "Final object position: (" << x[0][0] << ", "<< x[0][1] << ", " << x[0][2] <<  ")";
    return;
  }
  
  delete[] bucket;
  
  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
              << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
