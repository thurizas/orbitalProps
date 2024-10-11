// orbitalProps.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "constants.h"
#include "common.h"

#include "XGetopt.h"
#include "logger.h"


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

void CART2CLASS(psystemT, vec3, vec3, porbitalProp1T);        // converts cartesian (r and r') to to classical elements (a, e, M, i, W and w)
void CLASS2CART(psystemT, vec3*, vec3*, orbitalProp1T);       // convers classical elements (a, e, M, i, W and w) to cartesian (r and r')
void readFloat(const char*, double*, double conv = 1.0);
int readVector(const char*, vec3*);
double calDateToJulian(std::string);
double correctQuadrent(double, double, double);
void showUsage(const char*);

uint8_t LOC2CLASS = 0x01;
uint8_t CLASS2LOC = 0x02;
std::vector<orbitalProp1T>   refData;

double FNASN(double x) { return (atan(x / sqrt(-x * x + 1))); }
double FNACN(double x) { return (- atan(x / sqrt(-x * x + 1)) + (PI / 2)); }

void cmdOut(char* msg) { std::cout << msg; }

int main(int argc, char** argv)
{
    int choice = -1;
    int dbgLvl = 3;

    uint8_t mode = 0x00;

    while (-1 != (choice = getopt(argc, argv, "cpdh")))
    {
        switch (choice)
        {
        case 'c':
            mode |= LOC2CLASS;
            break;
        case 'p':
            mode |= CLASS2LOC;
            break;
        case 'd':
            dbgLvl--;
            if (dbgLvl < 1) dbgLvl = 1;
            if (dbgLvl > 6) dbgLvl = 6;
            break;
        case 'h':
            showUsage(argv[0]);
            exit(1);
        case '?':
        default:
            showUsage(argv[0]);
            exit(1);
        }
    }

    CLogger* pLogger = CLogger::getInstance();
    pLogger->regOutDevice(0, cmdOut);
    pLogger->outMsg(cmdLine, CLogger::level::NOTICE, "logging enabled\n");

    // read in solar system data for use in testing
    for (auto p: planets)
    {
        orbitalProp1T   refPlanet;
        
        strncpy(refPlanet.n, p.n, sizeof(p.n) - 1);
        refPlanet.a = p.a;
        refPlanet.e = p.e;
        refPlanet.m = p.m;
        refPlanet.i = p.i * DEG2RAD;                // modify values as appropriate.
        refPlanet.W = p.W * DEG2RAD;
        refPlanet.o = p.o;
        refPlanet.w = (p.w_bar - p.W) * DEG2RAD;    
        refPlanet.M = (p.L - p.w_bar) * DEG2RAD;

        refData.push_back(refPlanet);
    }

    if (mode == 0x00)
    {
        std::cerr << "[-] you must chose a mode of operation, exiting" << std::endl;
    }
    else if (mode == 0x03)
    {
        std::cerr << "[-] you must choose exactly one mode of operation, exiting" << std::endl;
    }
    else
    {
        // set up information about the overall system
        systemT  theSystem;
        theSystem.m_prim = 1.0;

        orbitalProp1T theOrbit = {0};              // classical properties of the orbit
        vec3 loc, vel;                             // cartesian properties of the orbit (position vector and time derivitive of position vector)


        if (mode == LOC2CLASS)
        {
            int converted;
            
            //  get the location and velocity from the user
            if (3 != (converted = readVector("Enter the position vector: ", &loc)))
            {
                std::cerr << "[-] input is malformed, require three components, only found " << converted << std::endl;
                exit(1);
            }
            if (3 != (converted = readVector("Enter the velocity vector: ", &vel)))
            {
                std::cerr << "[-] input is malformed, require three components, only found " << converted << std::endl;
                exit(1);
            }

            CART2CLASS(&theSystem, loc, vel, &theOrbit);

            std::cout << theOrbit;

        }
           
        if (mode == CLASS2LOC)
        {
            // get the orbital parameters from the user
            readFloat("enter the semi-major axis, in AU               : ", &theOrbit.a, static_cast<float>(AU2M));
            readFloat("enter the eccentricity of the orbit            : ", &theOrbit.e);
            readFloat("enter the inclination of the orbit, in degrees : ", &theOrbit.i, static_cast<float>(DEG2RAD));
            readFloat("enter the mean anomoly                         : ", &theOrbit.M);
            readFloat("enter the long. of ascending node, in degrees  : ", &theOrbit.W, static_cast<float>(DEG2RAD));
            readFloat("enter the argument of the periapsis, in degrees: ", &theOrbit.w, static_cast<float>(DEG2RAD));
            readFloat("enter the mass of the planet, in earth masses  : ", &theOrbit.m, static_cast<float>(MEARTH));

            std::cout << theOrbit;

            CLASS2CART(&theSystem, &loc, &vel, theOrbit);
            
            std::cout << "[+] the location vector is: " << loc;
            std::cout << "[+] the velocity vector is: " << vel;
        }
    }
 
    pLogger->outMsg(cmdLine, CLogger::level::NOTICE, "logging disabled");
    pLogger->delInstance();
    return 0;
}

void showUsage(const char* name)
{

}

void readFloat(const char* prompt, double* v, double conv)
{
    float val;
    std::cout << prompt;
    std::cin >> val;

    *v = static_cast<double>(val * conv);
}

int readVector(const char* prompt, vec3* vec)
{
    std::string line;
    std::string nbr;
    int converted = 0;                             // number of components read from the command line

    std::cout << prompt;
    std::getline(std::cin, line);

    bool processing = false;
    for(char c : line)
    {
        if (c == ',' || c == ' ' || c == '\n')
        {
            if(processing)
            {
                float f = static_cast<float>(atof(nbr.c_str()));
                vec->setComponent(converted++, f);
                nbr.erase();
                processing = false;
            }

            continue;
        }
        else if (isdigit(c) || c == '.' || c == '-' || c == '+')
        {
            processing = true;
            nbr.append(1, c);
        }
        else
        {
            std::cerr << "[-] vector a glyph other then a number, comma, decimal point or space - malformed entry" << std::endl;
            processing = false;
            break;
        }
    }

    // if processing is set and we got here, we have a number waiting to be converted and we hit end-of-line
    if(processing)
    { 
        float f = static_cast<float>(atof(nbr.c_str()));
        vec->setComponent(converted++, f);
        nbr.erase();
    }

    return converted;
} 

// expected format: year-month-day:hour-minute-second
double calDateToJulian(std::string calDate)
{
    int year = 0,month = 0,day = 0;        // y in range [1901,2099], m in range [1, 12], d in range [1, 31]
    int hour = 0,minute = 0,second = 0;    // h in range [0,24), m in range [0, 60), s in range [0, 60)
    double J0 = 0.0;                       // Julian date at 0H UTC on given date
    std::string date = "", time = "";
    std::string::size_type loc = std::string::npos;

    if (std::string::npos != (loc = calDate.find(':')))                  // found the colon delimiter
    {
        date = calDate.substr(0, loc);
        time = calDate.substr(loc + 1, calDate.length() - loc);
    }
    else                                                                 // no delimiter, assume all string is the date part.
    {
        date = calDate;
    }

    std::vector < std::string> dateParts;
    std::vector < std::string> timeParts;
    if (date != "")
    {
        std::string::size_type start = 0;
        do{
            if (std::string::npos != (loc = date.find("-", start)))    // 1985-August-3  first string 0:4 ; second string 5:10; third string 12:13
            {
                dateParts.push_back(date.substr(start, loc-start));
                start = loc + 1;
            }
            else
            {
                dateParts.push_back(date.substr(start, date.length() - start));
            }

        } while (loc != std::string::npos);
        
        // validate what we've read.. 
        if (dateParts.size() != 3) { std::cerr << "calendar date string is malformed, failed to find expected tokens in date part" << std::endl; exit(1); }
        else
        {
            year = stoi(dateParts[0]);
            if ((year < 1901) || (2099 < year)) { std::cerr << "year must be in the range [1901, 2099]" << std::endl; exit(1); }

            int n = 0;
            for (auto m : months)
            {
                n++;
                if (strcmp(m, dateParts[1].c_str()) == 0)
                    break;
            }
            month = n;
            if ((month < 1) || (12 < month)) { std::cerr << "month must be in the range [1, 12]" << std::endl; exit(1); }

            day = stoi(dateParts[2]);
            if ((day < 1) || (31 < day)) { std::cerr << "day must be in the range [1,31]" << std::endl; exit(1); }
        }
    }

    if (time != "")
    {
        std::string::size_type start = 0;
        do {
            if (std::string::npos != (loc = time.find("-", start)))    
            {
                timeParts.push_back(time.substr(start, loc - start));
                start = loc + 1;
            }
            else
            {
                timeParts.push_back(time.substr(start, time.length() - start));
            }

        } while (loc != std::string::npos);

        // validate what we've read... 
        if (timeParts.size() != 3) { std::cerr << "calendar date string is malformed, failed to find expected tokens in time part" << std::endl; exit(1); }
        else
        {
            hour = stoi(timeParts[0]);
            if ((hour < 0) || (23 < hour)) { std::cerr << "hour must be in the range [0,24)" << std::endl; exit(1); }

            minute = stoi(timeParts[1]);
            if ((minute < 0) || (59 < minute)) { std::cerr << "minute must be in the range [0, 60)" << std::endl; exit(1); }

            second = stoi(timeParts[2]);
            if ((second < 0) || (59 < second)) { std::cerr << "second must be in the range [0, 60)" << std::endl; exit(1); }
        }
    }

    J0 = (367 * year) - (floor((7 * (year + floor((month + 9) / 12))) / 4)) + (floor(275 * month / 9)) + day + 1721013.5;  // calculate J0 on date in question @ 0.0h UTC
    double UT = (double)hour + (((double)minute + ((double)second / 60.0)) / 60.0);                                        // convert h:m:s to h.m form.

     return J0 + (UT / 24.0);
}


/***********************************************************************************************************************
 * function  : 
 *
 * abstract  : converts cartesian (r and r') to to classical elements (a, e, M, i, W and w)
 *             derived from program 4.6.1 in Boulet, Dan. "Methods of Orbit Determination for the Micro Computer". 
 *             Willman-Bel:Richmond, Virginia : 1991
 *             line numbers refer to the code listing in the above refernece, pages 165-173
 * 
 * N.B. Initial positions (au) and velocities (au/day) for Sun and planets on Julian day (1969-June-28) 2440400.5
 *      see https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
 * 
 * object  vector        x-component             y-component                 z-component
 * Sun     position   0.00450250878464055477   0.00076707642709100705   0.00026605791776697764
           velocity  –0.00000035174953607552   0.00000517762640983341   0.00000222910217891203
Mercury    position   0.36176271656028195477  –0.09078197215676599295  –0.08571497256275117236
           velocity   0.00336749397200575848   0.02489452055768343341   0.01294630040970409203
Venus      position   0.61275194083507215477  –0.34836536903362219295  –0.19527828667594382236
           velocity   0.01095206842352823448   0.01561768426786768341   0.00633110570297786403
Mars       position  –0.11018607714879824523  –1.32759945030298299295  –0.60588914048429142236
           velocity   0.01448165305704756448   0.00024246307683646861  –0.00028152072792433877
Jupiter    position  –5.37970676855393644523  –0.83048132656339789295  –0.22482887442656542236
           velocity   0.00109201259423733748  –0.00651811661280738459  –0.00282078276229867897
Saturn     position   7.89439068290953155477   4.59647805517127300705   1.55869584283189997764
           velocity  –0.00321755651650091552   0.00433581034174662541   0.00192864631686015503
Uranus     position –18.26540225387235944523  –1.16195541867586999295  –0.25010605772133802236
           velocity   0.00022119039101561468  –0.00376247500810884459  –0.00165101502742994997
Neptune    position –16.05503578023336944523 –23.94219155985470899295  –9.40015796880239402236
           velocity   0.00264276984798005548  –0.00149831255054097759  –0.00067904196080291327
Pluto      position –30.48331376718383944523  –0.87240555684104999295   8.91157617249954997764
           velocity   0.00032220737349778078  –0.00314357639364532859  –0.00107794975959731297
 * 
 * 
 * 
 * 
 *
 * parameters:
 *
 * returns   :
 *
 * written   : (GKHuber)
 ***********************************************************************************************************************/
void CART2CLASS(psystemT psys, vec3 loc, vec3 vel, porbitalProp1T pprop)
{
    std::string calDate = "";
    std::string name = "";
    double obliquity = 0.0;
    double mass = 0.0;

    std::cout << "\n\n*************************************************************************\n\n";
    std::cout << __FUNCTION__ << " converting position and velocity to classical elements" << std::endl;

    // enter the calander date and calculate the julian time
    std::cout << "enter the calander date, year-month-day:hour-minute-second format: ";
    std::cin >> calDate;

    // enter the name of the planet
    std::cout << "enter the name of the object: ";
    std::cin >> name;

    // enter the obliquity of the ecliptic
    std::cout << "enter the obliquity of the ecliptic: "; 
    std::cin >> obliquity;
   
    bool bFound = false;
    for(auto p : refData)
    {
        if (strcmp(name.c_str(), p.n) == 0)
        {
            bFound = true;
            mass = p.m;
        }
    }
    if(!bFound)
    {
        // TODO : object is not in our reference data.
    }
 
    double mu = 0.0, jd=0.0;
    mu = psys->m_prim + (mass*MEARTH)/MSOL;             
    jd = calDateToJulian(calDate);
    
    // TODO : print out entered data:
    std::cout << "\n\n************************************************************\n\n" << std::endl;
    std::cout << "[+] supplied parameters" << std::endl;
    std::cout << "[+] name                            : " << name << std::endl;
    std::cout << "[+] epoch used                      : " << epoch << std::endl;
    std::cout << "[+] gaussian gravitational constant : " << std::setprecision(10) << K << std::endl;
    std::cout << "[+] obliquity of the ecliptic       : " << obliquity << std::endl;
    std::cout << "[+] combined mass                   : " << std::setprecision(10) << mu << std::endl;
    std::cout << "[+] Julian data                     : " << std::setprecision(10) << jd << " (" << calDate << ")" <<std::endl;
    std::cout << "[+] the location vector is          : " << loc;
    std::cout << "[+] the velocity vector is          : " << vel;
    
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] combined mass is: %8.6f\n", mu);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] t0 is %8.2f\n", jd);
    
    // convert obliquity of ecliptic to radian (line 1500)
    double EC = obliquity * DEG2RAD;
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] EC (line 1500) is %8.6f\n", EC);

    // convert position vector to ecliptic  (lines 1520 - 1540)
    vec3 R;
    R.x(loc.x());
    R.y(loc.y() * cos(EC) + loc.z() * sin(EC));
    R.z(loc.z() * cos(EC) - loc.y() * sin(EC));
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] R (lines 1520-1540) is <%8.6f,%8.6f,%8.6f>\n", R.x(), R.y(), R.z());
    std::cout << "cos is: " << cos(EC) << " sin is: " << sin(EC) << std::endl;
    // convert velocity vector to ecliptic (line 1560 - 1580)
    vec3 V;
    V.x(vel.x());
    V.y(vel.y() * cos(EC) + vel.z() * sin(EC));
    V.z(vel.z() * cos(EC) - vel.y() * sin(EC));
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] V (lines 1560-1580) is <%8.6f,%8.6f,%8.6f>\n", V.x(), V.y(), V.z());

    // calculate vectors E (eccentricity vector, eqn 4.17 and 4.69 ) ,H (angular momentum, eqn 4.72) and 
    // N (ascending node vector, eqn 4.74)
    double RMag = R.len();
    double VLen2 = V.len2();
    double RV = R.dot(V);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] length of R, is %8.6f\n", RMag);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] length^2 of V is %8.6f\n", VLen2);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] R dot V is %8.6f\n", RV);

    // eccentricity vector
    vec3 E;
    double F1 = ((VLen2 / mu) - (1 / RMag));
    double F2 = (RV / mu);
    for (int ndx = 0; ndx < 3; ndx++)
    {
        E.set(ndx, F1 * R.at(ndx) - F2 * V.at(ndx));
    }
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] factor, F1, is: %8.6f\n", F1);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] factor, F2, is: %8.6f\n", F2);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] eccentricity vector (eqn 4.17, 4.69) is <%8.6f,%8.6f,%8.6f>\n", E.x(), E.y(), E.z());

    // angular momentum vector
    vec3 H = R.cross(V);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] angular momentum vector (eqn 4.72) is <%8.6f,%8.6f,%8.6f>\n", H.x(), H.y(), H.z());

    // ascending node vector
    vec3 N;
    N.x(-H.y());
    N.y(H.x());
    N.z(0);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] ascending node vector (eqn 4.74) is <%8.6f,%8.6f,%8.6f>\n", N.x(), N.y(), N.z());

    // calculate values for  A recipercol of semi-major axis (eqn. 4.78), E eccentricity (eqn. 4.79) and 
    // Q perifocal distance (eqn. 4.82)

    double AI = (2 / RMag) - (VLen2 / mu);
    double EMag = E.len();
    double SP = H.len2() / mu;
    double Q = SP/(1 + EMag);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] inverse of semi-major axis (eqn 4.78) is: %8.6f\n", AI);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] eccentricity (eqn 4.79) is: %8.6f\n", EMag);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] factor, SP, is: %8.6f\n", SP);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] perifocal distance (eqn 4.82) Q, is: %8.6f\n", Q);

    // calculate values for I angle of inclination (eqn. 4.84), OO Longitude of ascending node (eqn. 4.86)
    // and W argument of perifocus (eqn. 4.88) which gives longitude of perifocus (w_bar = W + w)

    // angle of inclination
    double HLen = H.len();
    double I = acos((H.z() / HLen) / DEG2RAD);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] length of H is, %8.6f\n", HLen);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] angle of inclination (eqn 4.84), is %8.6f\n", I);

    // longitude of ascending node
    double NMag = N.len();
    double OO = acos((N.x() / NMag) / DEG2RAD);
    if (N.y() < 0) { OO = 2 * PI - OO; }
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] length of N is %8.6f\n", NMag);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] longitude of ascending node (eqn 4.86), is %8.6f\n", OO);

    // longitude of perifocus
    double NE = N.dot(E);
    double W = acos((NE / (EMag * NMag)) / DEG2RAD);
    if (E.z() < 0) { W = 2 * PI - W; }
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] N dot E  is %8.6f\n", NE);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] longitude of perifocus (eqn 4.88), is %8.6f\n", W);

    // x_bar (eqn 4.91), and y_bar (eqn. 4.92/)
    double XB = (SP - RMag) / EMag;
    double YB = RV * sqrt((SP / mu) / EMag);
    double PE = 0.0;
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] x-bar (eqn 4.91), is %8.6f\n", XB);
    CLogger::getInstance()->outMsg(cmdLine, CLogger::level::NOTICE, "[?] longitude of perifocus (eqn 4.92), is %8.6f\n", YB);
    // determine the type of orbit, based on the eccentricity
    double AQ = 0.0;
    double MM = 0.0;
    if (fabs(1 - EMag) < 0.001)            // we have parabolic motion
    {
        // line 2430
    }
    else if (EMag > 1)                      // we have hyperbolic motion
    {
        // line 2300
    }
    else if (EMag < 1)                      // we have elliptic motion line 2120
    {
        std::cout << "[+] calculating for an elliptical orbit" << std::endl;
        AQ = 1 / AI;
        double B = (1 / AI) * sqrt(1 - EMag * EMag);
        double CX = XB * AI + EMag;
        double SX = YB / B;
        double E1 =  correctQuadrent(asin(SX), SX, CX);  // line 11010
        MM = E1 - EMag * SX;                             // eqn 4.98
        double MT = MM / DEG2RAD;                        // convert mean anomaly to degrees
        double n = K*AI*sqrt(mu*AI);                     // mean motion eqn 4.99
        double TP = jd - MM / n;                         // time of perifocus passage, eqn 4.35
        PE = (2*PI/K)*sqrt(1/mu*AI*AI*AI);               // eqn 3.52

    }

    // line 2530
    std::cout << "[+] Classical elements: " << std::endl;
    std::cout << "    semi-major axis (A) : " << AQ << std::endl;
    std::cout << "    eccentricity (E)    : " << E << std::endl;
    std::cout << "    mean anomaly (MM)   : " << MM << std::endl;
    std::cout << "    inclination (I)     : " << I << std::endl;
    std::cout << "    longitude of ascending node : " << OO << std::endl;
    std::cout << "    argument of periapsis       : " << W << std::endl;




}


// convers classical elements (a, e, M, i, W and w) to cartesian (r and r')
void CLASS2CART(psystemT psys, vec3* ploc, vec3* pvec, orbitalProp1T prop)
{


}

double correctQuadrent(double X, double SX, double CX)
{
    if (fabs(SX) <= 0.707107) {X = FNASN(abs(SX)); }
    if (fabs(CX) <= 0.707107) {X = FNACN(abs(CX)); }
    if ((CX >= 0) && (SX >= 0)) { X = X; }
    if ((CX < 0) && (SX >= 0)) { X = 180 * DEG2RAD - X; }
    if ((CX < 0) && (SX < 0)) { X = 180 * DEG2RAD + X; }
    if ((CX >= 0) && (SX < 0)) { X = 180 * DEG2RAD - X; }
    return X;
}
