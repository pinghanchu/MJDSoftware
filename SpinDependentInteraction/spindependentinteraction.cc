// Monte Carlo calculation by Pinghan Chu
// 2015.2.14
// modify for sigma.r/sigma.v/sigma.(vxr) for polarized 3He
//
// 2015.2.16
// modify the cell with curved surface
// mass still with flat surface
// 2015.3.28
// add other exotic forces
#include<iostream>
#include<iomanip>
#include<fstream>
#include<time.h>
#include<math.h>
#include<cmath>
#include<vector>
#include<cstdlib>
#include<ctime>
using namespace std;

//Input Parameters:******************************************************************************

//Nucleon density parameter for all the materials used as test substances per m^3
const long double unit = 1.0e0;
const long double Mn = 1.67e-27;
const long double Me = 9.11e-31;
const long double C = 299792458.*pow(unit,2);
const long double hbar = 1.055e-34*pow(unit,2);
const long double pi= 3.14159265359;
const long double particle_density= 4.29e30*pow(unit,-3); //BGO m^-3
//const long double particle_density= 6e28*pow(unit,-3); //AlNiCo m^-3; it can be spin density; PRL97, 021603
//const long double particle_density= 1e26*pow(unit,-3);//DyIG m^-3
//const long double normalization = 2.6518477625e-43*pow(unit,4);
//normalization=hbar^2/8pi/m_n  (kg^2.m^4/s^4/kg = kg.m^4/s^4)
//normalization=(1.055e-34*1.055e-34)/(25.1327*1.67e-27)
//hbar = 1.055e-34 m^2 kg/s
//mn = 	1.67e−27kg
//me = 9.11e-31kg
//const long double detector_density= 1.00e28*pow(unit,-3); //GGG m^-3
//const long double kappa=1.38e-23*pow(unit,2); //Boltzmann constant JK^-1 = [m^2 kg s^-2 K^-1]

const long double normalization2  =  pow(hbar,1)*C/(4*pi);
const long double normalization3  =  pow(hbar,3)/(4*pi*Me*Me*C);
const long double normalization4  = -pow(hbar,2)/(8.*pi*Me*C);
const long double normalization6  = -pow(hbar,2)/(4.*pi*Me*C);
const long double normalization8  =  pow(hbar,1)/(4.*pi*C);
const long double normalization9  =  pow(hbar,2)/(8.*pi*Me);
const long double normalization11 = -pow(hbar,2)/(4.*pi*Me);
const long double normalization12 =  pow(hbar,1)/(8.*pi);
const long double normalization14 =  pow(hbar,1)/(4.*pi);
const long double normalization15 = -pow(hbar,3)/(8.*pi*pow(C,2)*Me*Me);
const long double normalization16 = -pow(hbar,2)/(8.*pi*Me*C*C);
//hbar^3 /8pi c^2 /m_n m_e 
//////////////////////////////////////////////////////
const int MAXNUM = pow(2,8);//2^20
const long double cycles = 120;
const long double gyroratio = 4.396*pow(10,10);
const long double magsen = 3.81*pow(10,-15); 
const long double freqshift = gyroratio*magsen/sqrt(cycles);//Hz
const long double amp = 0.005*unit; // source mass amplitude; modulation amplitude 
const long double gap = 0.005*unit; // gap between source mass and detector (3He)
const long double omega = 2*pi*1; // 1Hz
const long double V = amp*omega;//~0.003 m/s
const long double shape1[3]={0.003*unit,0.003*unit,0.003*unit};
// shape of detector in cartesian coordinate
const long double shape2[3]={0.02*unit,0.02*unit,0.02*unit};
// shape of source mass in cartesian coordinate
const long double bounds1[6]={shape1[0]*0.5*unit,-shape1[0]*0.5*unit,shape1[1]*0.5*unit,-shape1[1]*0.5*unit,shape1[2]*0.5*unit,-shape1[2]*0.5*unit};
// boundaries of detector in cartesian coordinates
long double maxlength = sqrt(pow(shape2[0],2)*3)+gap;
//long double maxlength = gap + sqrt(pow(shape1[0],2)*2+pow(shape1[2]*2,2))+sqrt(pow(shape2[0],2)*3);
//long double maxlength = sqrt(pow(shape1[0]*0.5+shape2[0]*0.5,2)*2+pow(shape1[2]+shape1[0]*0.5+gap+shape2[2],2));
const long double bounds2[6]={maxlength*unit,gap,0.5*pi,0.,2.*pi,0.};
// boundaries of source mass in spherical coordinates (rho, theta, phi)
// rho = sqrt(0.05^2+0.05^2+0.05^2); the maximum distance in the mass
// rho + shape1[3]*2 + gap = 

const long double D[3] = {0.0*unit,0.0*unit, -(gap)};
// the distance between the mass and the detector due to cell thickness or other factor
//Functions defining the geometry and form of the potential:

long double f2(long double[],long double);

long double f3(long double[],long double);

long double f4(long double[],long double);

long double f6(long double[],long double);

long double f8(long double[],long double);

long double f9(long double[],long double);

long double f11(long double[],long double);

long double f12(long double[],long double);

long double f14(long double[],long double);

long double f15(long double[],long double);

long double f16(long double[],long double);

long double *dice(long double x[],long double a[],long double b[],int n){
  // x is the Monte Carlo generated points
  // a is the boundary of the lower limit
  // b is the boundary of the upper limit
  // n is the parameter number

  double max=RAND_MAX;
  for (int j=0;j<n;j++){
    x[j]=a[j] + (b[j]-a[j])*rand()/max;
  }
    
  return x;
}

int inside(long double lower[], long double upper[], long double point[],int n){
  // check point[] inside the box of lower[] and upper[]

  int x = 0;
  if(point[0] < lower[0] || point[0]>upper[0]){
    //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
    return x;
  }else if(point[1]< lower[1] || point[1]>upper[1]){
    //cout << "y " << point[1] << " " << lower[1] << " " << upper[1] <<endl;
    return x;
  }else if(point[2]<lower[2] || point[2]>upper[2]){
    //cout << "z " << point[2] << " " << lower[2] << " " << upper[2] << endl;
    return x;
  }else{
    x = 1;
    //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
    //cout << "y " << point[1] << " " << lower[1] << " " << upper[1] << endl;
    //if(point[2]>0.05){
    //cout << "z " << point[2] << " " << lower[2] << " " << upper[2] << endl;
    //}
          return x;
  }
}
int insidecylinder(long double lower[], long double upper[], long double point[],int n){
    // check point[] inside the box of lower[] and upper[]
  //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
  //cout << "y " << point[1] << " " << lower[1] << " " << upper[1] << endl;
  //cout << "z " << point[2] << " " << lower[2] << " " << upper[2] << endl;
    
  int x = 0;
  if(sqrt(pow(point[0],2)+pow(point[1],2))>upper[0] || sqrt(pow(point[0],2)+pow(point[1],2)) < lower[0]){
    x=0;
  }else if(point[2]<lower[2] || point[2]>upper[2]){
    x=0;
  }else{
    x = 1;
    //cout << "mass x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
    //cout << "mass y " << point[1] << " " << lower[1] << " " << upper[1] << endl;
    //cout << "mass z " << point[2] << " " << lower[2] << " " << upper[2] << endl;
    //return x;
  }
  return x;
}




int insidehalfcellcylinder(long double lower[], long double upper[], long double point[],int n){
    // check point[] inside the box of lower[] and upper[]
    //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
    //cout << "y " << point[1] << " " << lower[1] << " " << upper[1] << " " << sqrt(pow(upper[0],2)-pow(point[0],2)) << endl;
    //cout << "z " << point[2] << " " << lower[2] << " " << upper[2] << " " << upper[2]+sqrt(pow(upper[0],2)-pow(point[0],2)-pow(point[1],2)) << endl;
    
    int x = 0;
        
    if( (point[0]>lower[0] && point[0]<upper[0])){
        //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
        if(point[1]<sqrt(pow(upper[0],2)-pow(point[0],2))
        && point[1]>-sqrt(pow(upper[0],2)-pow(point[0],2))){
            if(point[2]< upper[2]+sqrt(pow(upper[0],2)-pow(point[0],2)-pow(point[1],2)) && point[2]> lower[2])
            {
                x=1;
            }else{
                x=0;
            }
        }
    }
    return x;
}

int insidecellcylinder(long double lower[], long double upper[], long double point[],int n){
    // check point[] inside the box of lower[] and upper[]
    //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
    //cout << "y " << point[1] << " " << lower[1] << " " << upper[1] << " " << sqrt(pow(upper[0],2)-pow(point[0],2)) << endl;
    //cout << "z " << point[2] << " " << lower[2] << " " << upper[2] << " " << upper[2]+sqrt(pow(upper[0],2)-pow(point[0],2)-pow(point[1],2)) << endl;
    
    int x = 0;
        
    if( (point[0]>lower[0] && point[0]<upper[0])){
        //cout << "x " << point[0] << " " << lower[0] << " " << upper[0] << endl;
        if(point[1]<sqrt(pow(upper[0],2)-pow(point[0],2))
        && point[1]>-sqrt(pow(upper[0],2)-pow(point[0],2))){
            if(point[2]< upper[2]+sqrt(pow(upper[0],2)-pow(point[0],2)-pow(point[1],2)) && point[2]> lower[2]-sqrt(pow(lower[0],2)-pow(point[0],2)-pow(point[1],2)))
            {
                x=1;
            }else{
                x=0;
            }
        }
    }
    return x;
}

int insidesphere(long double lower[], long double upper[], long double point[],int n){
    // check point[] inside the box of lower[] and upper[]
    
    int x = 0;
    if(sqrt(pow(point[0],2)+pow(point[1],2)+pow(point[2],2))>upper[0] || sqrt(pow(point[0],2)+pow(point[1],2)+pow(point[2],2)) < lower[0]){
        return x;
    }else{
        x = 1;
        return x;
    }
}




long double f2(long double x[], long double forcelength){
  //mass moving forward and backward

  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = \sigma_1\dot\sigma_2 (1/\rho) exp(-\rho/\lambda)d\rho
  //*pow(\rho,2)*sin(\theta)d\theta;
  //= -\rho \lambda dr^\prime \sin(\theta)d\theta

  long double dftot = -forcelength*rho*sin(theta);

  //do a transform : r^\prime = \exp(-\rho/\lambda)
  // d\rho = -\lambda/r^\prime dr^\prime
  return dftot;
 }

long double f3(long double x[], long double forcelength){
  // mass moving forward and backward

  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (\sigma_1\dot\sigma_2)(1/\lambda\rho^2+1/\rho^3)-(\sigma_1\cdot\hat{r})(\sigma_2\cdot\hat{r})(1/\lambda^2\rho+3/\lambda\rho^2+3/\rho^3) exp(-\rho/\lambda)d\rho
  //*pow(\rho,2)*sin(\theta)d\theta;
  //=( (1/\lambda+1/\rho)-\cos^2\theta(\rho/\lambda^2+3/\lambda+3/\rho))r^\prime (-\lambda/r^\prime dr^\prime)\sin\theta d\theta
  //=-((1+\lambda/\rho)-\cos^2\theta(\rho/\lambda+3+3\lambda/\rho))dr^\prime sin\theta\d\theta

  long double dftot = -((1+forcelength/rho)-cos(theta)*cos(theta)*(rho/forcelength+3+3*forcelength/rho))*sin(theta);

  return dftot;
}

long double f4(long double x[], long double forcelength){
  /////////////////////////////////////////////
  //consider the mass is spinning at frequency omega
  //the speed at a point will be r\cdot \omega/C 
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (\hat{\sigma} \cdot (\vec{v})\time \hat{rho}) (1/rho\lambda+1/rho^2)exp(-rho/\lambda) rho^2 \sin(\theta)
  //= (\rho\sin(\theta)\omega)(-\sin(\theta))(1/rho\lambda+1/rho^2)rho^2 exp(-rho/\lambda) dr \sin\thetad\theta
  //= (-\rho\omega)(\rho/\lambda+1)r^\prime(-\lambda/r^\prime dr^\prime) \sin^3\theta d\ehta
  //= (\rho\omega)(\rho+\lambda)dr^\prime \sin^3 d\theta
  long double dftot = (forcelength+rho)*(rho*omega)*pow(sin(theta),3);

  ///////////////////////////////////////
  //velocity is perpendicular to sigma_1
  //long double dftot = V*(forcelength+rho)*pow(sin(theta),2);
  //do a transform : r^\prime = \exp(-r/\lambda)
  // dr = -\lambda/r^\prime dr^\prime
  //////////////////////////////////////
  //cout  << forcelength << " " << rho <<" " << dftot << endl;


  return dftot;
}

long double f6(long double x[], long double forcelength){
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (\sigma_1\cdot v)(\sigma\cdot\hat{r})(1/\rho\lmabda + 1/\rho^2) e^{-\rho/\lambda} \rho^2 d\rho \sin\theta d\theta
  //=(v)(\cos\theta)(\rho/\lambda+1)r^\prime(-\lambda/r^\prime dr^\prime)\sin\theta d\theta
  //= -v (\rho+\lambda) \cos\theta\sin\theta dr^\prime d\theta

  long double dftot = -V*(rho+forcelength)*(sin(theta)*cos(theta));
  return dftot;
}

long double f8(long double x[], long double forcelength){
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (\sigma_1\cdot v)(\sigma_2\cdot v)(1/rho)*exp(-rho/forcelength);
  //*pow(rho,2)*sin(theta);
  //=v^2 (\rho)r^\prime(-\lambda/r^\prime dr^\prime) \sin\theta d\theta
  //=v^2 (-\rho\lambda)\sin\theta dr^\prime d\theta

  long double dftot = -V*V*(forcelength*rho)*sin(theta);

  //do a transform : r^\prime = \exp(-r/\lambda)
  // dr = -\lambda/r^\prime dr^\prime

  //cout  << forcelength << " " << rho <<" " << dftot << endl;
  return dftot;
}

long double f9(long double x[], long double forcelength){
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (1./pow(\rho,2)+1./(\rho\lambda))*cos(theta)*exp(-rho/forcelength);
  //*pow(rho,2)*sin(theta);
  //=(1+\rho/\lambda)(-\lambda dr^\prime)(\cos\theta\sin\theta d\theta)
  //=-(\lambda+\rho)(\cos\theta\sin\theta)dr^\prime d\theta

  long double dftot = -(forcelength+rho)*sin(theta)*cos(theta);

  //do a transform : r^\prime = \exp(-r/\lambda)
  // dr = -\lambda/r^\prime dr^\prime
  //cout  << forcelength << " " << rho <<" " << dftot << endl;
  return dftot;
}

long double f11(long double x[], long double forcelength){
  // the spin direction of mass is perpendicalar to detector Rb gas
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = ((\sigma_1\cross\simga_2)\cdot\hat{r})(1./pow(\rho,2)+1./(\rho\lambda))*exp(-rho/forcelength);
  //*pow(rho,2)*sin(theta);
  //=(-\sin\theta)(1+\rho/\lambda)(-\lambda dr^\prime)(\sin\theta d\theta)
  //=(\lambda+\rho)(\sin\theta\sin\theta)dr^\prime d\theta

  long double dftot = (forcelength+rho)*sin(theta)*sin(theta);

  //do a transform : r^\prime = \exp(-r/\lambda)
  // dr = -\lambda/r^\prime dr^\prime
  //cout  << forcelength << " " << rho <<" " << dftot << endl;
  return dftot;
}



long double f12(long double x[], long double forcelength){
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;
  //long double dftot = (\sigma\cdot v)(1/\rho)e^{-\rho/\lambda}\rho^2 d\rho \sin\theta d\theta
  //=v*(\rho)(r^\prime)(-\lambda/r^\prime dr^\prime)\sin\theta d\theta;
  //=v(-\rho\lambda)\sin\thetadr^\prime d\theta
    
  long double dftot = -V*(forcelength*rho)*sin(theta);

  //do a transform : r^\prime = \exp(-r/\lambda)
  // dr = -\lambda/r^\prime dr^\prime
  //cout  << forcelength << " " << rho <<" " << dftot << endl;
  return dftot;
}


long double f14(long double x[], long double forcelength){
  //mass spin is perpendicular to detector spin and velocity is perpendicular to these two spins.
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];                                                     
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;                                       
  long double dftot = -V*forcelength*rho*sin(theta);
  //cout << x[3] << " " << rho << " " << theta << " " << rho*(3+3*rho/forcelength+pow((rho/forcelength),2)) <<" " << dftot << endl;
  return dftot;
}

long double f15(long double x[], long double forcelength){
  //rotating mass
  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];                                                     
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;                                       
  long double dftot = 2*omega*pow(sin(theta),3)*cos(theta)*(rho*rho/forcelength+3*rho+3*forcelength);
  //velocity is prependicular to the sigma_1 and simga_2; simga_1 // sigma_2
  //long double dftot = 2*V*pow(sin(theta),2)*cos(theta)*(rho/forcelength+3+3*forcelength/rho);

  //cout << x[3] << " " << rho << " " << theta << " " << rho*(3+3*rho/forcelength+pow((rho/forcelength),2)) <<" " << dftot << endl;
  return dftot;
}

long double f16(long double x[], long double forcelength){
  //case 1 : sigma_2 and velocity are parallel and perpendicular to sigma_1
  //case 2 : sigma_1 and velocity are parallel and perpendicular to sigma_2

  long double rho = -forcelength*log(x[3]);
  //long double rho = x[3];                                                     
  long double theta = x[4];
  //cout << x[3] << " " << rho << endl;                                       
  long double dftot = V*V*sin(theta)*sin(theta)*(rho+forcelength);

  return dftot;
}

long double integration(long double range){


  int n = 6;

  long double up1[3];
  long double low1[3];
  long double up2[3];
  long double low2[3];
  //long double asinom = amp*sin(omega*1*0.25); // the distance because of moving; t=1/4, omega=2pi*1, maximum distance
  long double asinom = 0.;
////////////////////////////////////////////////////
  // setup the boundary corrdinate of detector and mass
  // detector
  up1[0] = shape1[0]/2.;
  low1[0] = -shape1[0]/2.;
  up1[1] = shape1[1]/2.;
  low1[1] = -shape1[1]/2.;
  up1[2] = shape1[2]/2.;
  low1[2] = -shape1[2]/2.;

  // mass
  up2[0] = shape2[0]/2.;
  low2[0] = -shape2[0]/2.;
  up2[1] = shape2[1]/2.;
  low2[1] = -shape2[1]/2.;  
  up2[2] = shape2[2]+gap+asinom+ shape1[2]/2;
  low2[2] = gap+asinom+shape1[2]/2;
  for(int i=0;i<3;i++){
    cout << "detector up="<<up1[i] << " detector low=" << low1[i] <<  " mass up=" << up2[i] << " mass low=" << low2[i] << endl;
  }
/////////////////////////////////////////////////  

  // generate Monte Carlo points; the first three is th coordinate
  // in the detector; the last three is the coordinate in the mass
  long double x[n];

  long double a[n],b[n];

  // boundary of Monte Carlo generting points
  // boundary of the detector
  for(int i = 0 ; i < n/2 ;i++){
    b[i] = bounds1[2*i]; 	
    a[i] = bounds1[2*i+1];
  }
// b[0] = bounds1[0] detector upper bound
// a[0] = bounds1[1].. detector lower bound
    
// boundary of the mass
  for(int i = n/2 ; i < n ;i++){
    b[i] = bounds2[2*i-6];
    a[i] = bounds2[2*i+1-6];
  }  
// b[3] = bounds2[0] mass upper bound
// a[3] = bounds2[1].. mass lower bound

  //for(int i=0;i<n;i++){
  //  cout << " up="<<b[i] << " low=" << a[i] << endl;
  //}
  //cout << a[3] << " " << b[3] << endl;    
  // transform mass to detector distance (radius) bound into exp
  b[3] = exp(-b[3]/range);
  a[3] = exp(-a[3]/range);
  //  for(int i=0;i<n;i++){
  //  cout << " up="<<b[i] << " low=" << a[i] << endl;
  //}
  // calculate the total volumn based on the boundary
  long double volumn = 1.;

  for(int i=3;i<6;i++){
    volumn = volumn*(b[i]-a[i]);
  }
//  volumn = pow(up2[0],2)*pi*(up2[2]-low2[2]);
  cout << "volumn=" << volumn<< endl;

  int in1=0;
  int in2=0;

  //initialize random #
  srand(time(NULL));

  long double dftot = 0.; 
  long double W = 0.;
  int count = 0 ;
  int totalcount = 0;
  while (count < MAXNUM){
    dice(x,a,b,n);
    totalcount++;
    long double point1[n/2],point2[n/2];
    for(int i = 0;i<n/2;i++){
      point1[i] = x[i]; 
    }
    //long double R = x[3];
    long double R = -range*log(x[3]);
    long double theta = x[4];
    long double phi = x[5];
    
    //    cout << x[3] << " " << x[4] << " " << x[5]<< " " << cos(x[5]) << " " << sin(x[5])<<endl;
    point2[0] = point1[0] + R*sin(theta)*cos(phi);
    point2[1] = point1[1] + R*sin(theta)*sin(phi);
    point2[2] = point1[2] + R*cos(theta);
    //cout << point2[0] << " " << point2[1] << "  "<< point2[2] << endl;
    
    //point2[0] = point1[0] - D[0] + R*sin(theta)*cos(phi);  
    //point2[1] = point1[1] - D[1] + R*sin(theta)*sin(phi);  
    //point2[2] = point1[2] - D[2] + R*cos(theta);  
    /*
    long double Xp = point2[0]-point1[0];
    long double Yp = point2[1]-point1[1];
    long double Zp = point2[2]-point1[2];
    long double Rp = sqrt(pow(Xp,2)+pow(Yp,2)+pow(Zp,2));
    long double Thetap = acos(Zp/Rp);
    long double Phip = atan(Yp/Xp);
    x[3] = exp(-Rp/range);
    x[4] = Thetap;
    x[5] = Phip;
    */
    //cout << x[3] << " "  << x[4] << " " << x[5]  << " " << cos(x[5]) << " " << sin(x[5]) <<endl;
    in1 = insidecylinder(low1,up1,point1,3);
    in2 = inside(low2,up2,point2,3);
    //long double xx = f(x,range);
    if(in1>0 && in2>0){
      long double xx = f4(x,range);
      dftot = dftot + xx;
      count = count + 1;
      //cout << "count="<<count << ", dftot=" << dftot << endl;
    }
  }

  W = normalization4*particle_density*(dftot/count)*(volumn*count/totalcount);

  return W;

  //  return 1;
}
//Main****************************************************************************************************
int main(){

  long double range = 0.;
  long double x= 0.;
  long double y= 0.1;
  ofstream fout1;
  fout1.open("potential.txt");
    ofstream fout2;
    //fout2.open("./data/coupling.4.1.opt.txt");
  fout2.open("./data/coupling.4.3.txt");

  for(int i = 0; i< 6;i++){
    for(int j = 1 ; j <=9;j++){
      range = 0.00001*pow(10,i)*j*pow(unit,1); // force interaction length (m)
      x=integration(range);
      y=abs(freqshift*hbar/x);
      fout1 << range << " " << x << endl;
      fout2 << range << " " << log10(y) << endl;
    }
  }
  fout1.close();
  fout2.close();
  return 0;
}

//End Program*********************************************************************************************

