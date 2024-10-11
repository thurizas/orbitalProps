#ifndef _common_h_
#define _common_h_

#include "constants.h"
#include <cmath>

#include <iostream>

typedef struct _orbpro
{
    char     n[10];         // name of the object
    double   a;             // semi-major axis, AU
    double   e;             // eccentricity of orbit
    double   L;             // mean longitude, degrees
    double   i;             // angle of inclination of orbit, degrees
    double   W;             // longitude of ascending node, degrees
    double   w_bar;         // longitude of perihelion, degrees
    double   m;             // mass of the planet, earth masses
    double   o;             // obliquity of ecliptic, degrees
} orbitalPropT, * porbitalPropT;

typedef struct _orbpro1
{
    char     n[10];         // name of the object
    double    a;             // semi-major axis, AU
    double    e;             // eccentricity of orbit
    double    M;             // mean longitude, radians
    double    i;             // angle of inclination of orbit, radians
    double    W;             // longitude of ascending node, radians
    double    w;             // longitude of perihelion, degrees
    double    m;             // mass of the planet, kilograms
    double    o;             // obliquity of ecliptic, radians
} orbitalProp1T, * porbitalProp1T;

typedef struct _system
{
    float    m_prim;        // mass of the primary in solar masses
} systemT, *psystemT;

struct vec3
{
public:
    vec3() : m_comp{ 0.0f, 0.0f, 0.0f } { }
    vec3(float _x, float _y, float _z) : m_comp{_x, _y, _z} { } 

    float x() { return m_comp[0]; }
    void x(float _x) { m_comp[0] = _x; }
    float y() { return m_comp[1]; }
    void y(float _y) { m_comp[1] = _y; }
    float z() { return m_comp[2]; }
    void z(float _z) { m_comp[2] = _z; }

    vec3 operator+(const vec3& o) {vec3 res; for (int ndx = 0; ndx < 3; ndx++) { res.m_comp[ndx] = this->m_comp[ndx] + o.m_comp[ndx]; } return res;}
    vec3 operator-(const vec3& o) {vec3 res; for (int ndx = 0; ndx < 3; ndx++) { res.m_comp[ndx] = this->m_comp[ndx] + o.m_comp[ndx]; } return res;}

    bool operator==(const vec3& o) { bool ret = true; for (int ndx = 0; ndx < 3; ndx++) { ret &= (fabs(this->m_comp[ndx]-o.m_comp[ndx]) < 1E-6); } return ret; }

    float len2() { float sum = 0.0f;  for (int ndx = 0; ndx < 3; ndx++) { sum += m_comp[ndx] * m_comp[ndx]; } return sum; }
    float len() { return sqrt(len2()); }
    float dot(vec3& o) { float sum = 0.0f; for (int ndx = 0; ndx < 3; ndx++) { sum += m_comp[ndx] * o.getComp(ndx); } return sum; }
    vec3  cross(vec3& o) 
    { vec3 prdt; 
    float H_x = y() * o.z() - z() * o.y();
    float H_y = x() * o.z() - z() * o.x();
    float H_z = x() * o.y() - y() * o.x();
    prdt.x(H_x); prdt.y(-H_y); prdt.z(H_z); 
    return prdt; }

    float at(int ndx) { if ((ndx < 0) || (ndx > 2)) { throw std::runtime_error("index out of bounds"); } return m_comp[ndx]; }
    void set(int ndx, float val) { if ((ndx < 0) || (ndx > 2)) { throw std::runtime_error("index out of bounds"); } m_comp[ndx] = val;}
    
    // TODO : need better error handling here...
    void setComponent(int ndx, float val) { if ((0 <= ndx) && (ndx < 3)) m_comp[ndx] = val; }
    float getComp(int ndx) { float val = std::nanf("1");  if ((0 <= ndx) && (ndx < 3)) val = m_comp[ndx]; return val; }

private:
    float  m_comp[3];
};

std::ostream& operator<<(std::ostream& os, const vec3& vec)
{
    vec3 v = const_cast<vec3&>(vec);
    os << "<" << v.x() << "," << v.y() << "," << v.z() << ">" << std::endl;
    return os;
}


std::ostream& operator<<(std::ostream& os, const orbitalProp1T& prop)
{
    orbitalProp1T p = const_cast<orbitalProp1T&>(prop);
    os << "Classical Orbital Properties" << std::endl;
    os << "   semi-major axis, meters         : " << prop.a << std::endl;
    os << "   eccentricity of orbit           : " << prop.e << std::endl;
    os << "   mean anomoly                    : " << prop.M << std::endl;
    os << "   orbit inclination, radians      : " << prop.i << std::endl;
    os << "   long. of ascenging node, radians: " << prop.W << std::endl;
    os << "   argument of periapsis, radians  : " << prop.w << std::endl;
    os << "   mass of the planet              : " << prop.m << std::endl;
    return os;
}

#endif;
