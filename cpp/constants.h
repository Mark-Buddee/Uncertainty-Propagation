#ifndef constants_H
#define constants_H

void removeDupWord(std::string str, double *xarray)
{
  int i = 0;
  istringstream ss(str);

  string word; // for storing each word
  while (ss >> word)
  {
    // print the read word
    double number = std::stod(word);
    xarray[i] = number;
    //cout << xarray[i] << endl;
    i = i + 1;
  }
}


struct parameters{

  double u0[3];
  double LU;
  double TU;
  double tof;
  double Tmax;
  double N;
  double mstart;
  double mu_earth;
  double Re;
  double dt;
  double J2  = 1.08262668e-3; 
  double method;
  double Isp;
  double g0;
  double meanel;
  double Rs; 
  double eclipses;
  double nu;
  double Area; 
  double Cd; 
  double MU;

};

parameters  createParams(){


  parameters params;
  ifstream infile;

  infile.open("params.dat"); // open a file to perform read operation using file object

    infile >> params.LU;	
    infile >> params.TU;
    infile >> params.tof;
    infile >> params.Tmax;
    infile >> params.N;
    infile >> params.mstart;
    infile >> params.mu_earth ;
    infile >> params.Re;
    infile >> params.dt;
    infile >> params.method;
    infile >> params.Isp;
    infile >> params.g0; 
    infile >> params.meanel; 
    infile >> params.Rs; 
    infile >> params.eclipses; 
    infile >> params.nu; 
    infile >> params.Area; 
    infile >> params.Cd; 
    infile >> params.MU;

//     cout <<  params.LU<< endl;
//     cout <<  params.TU<< endl;
//     cout <<  params.tof<< endl;
//     cout <<  params.Tmax<< endl;
//     cout <<  params.N<< endl;
//     cout <<  params.mstart<< endl;
//     cout <<  params.mu_earth<< endl;
//     cout <<  params.Re<< endl;
// 
//      cout <<  params.dt<< endl;
//     cout <<  params.method<< endl;
//     cout <<  params.Isp<< endl;
//      cout <<  params.g0<< endl;
//     cout <<  params.meanel<< endl;
//     cout <<  params.Rs<< endl;
//     cout <<  params.eclipses<< endl;
//     cout <<  params.nu<< endl;
//     cout <<  params.Area<< endl;
//         cout <<  params.Cd<< endl;
//     cout <<  params.MU<< endl;
// 


 infile.close(); // close the file object.
  


  //  cout <<  params.TU<< endl;
  //  cout <<  params.LU<< endl;
   

  return params;
}

#endif