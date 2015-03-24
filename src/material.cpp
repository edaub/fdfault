#include <cassert>
#include <cmath>
#include "material.hpp"

material::material() {
    // constructor
    
    // set all properties to default values
    // units can be set to any desired units, however the default values have one thing
    // of note: default density has funny units of MPa s^2 / km / m such that measuring
    // velocity in m/s, stress in MPa, and distance in km will lead to correct unit scaling
    // similarly, with these choices, measuring modulii in GPa (= MPa km / m) will give
    // correct values. However
    
    rho = 2.67;
    lambda = 32.04;
    g = 32.04;
    mu = 0.5735;
    beta = 0.2867;
    eta = 0.2775;
}

void material::set_rho(const double newrho) {
    // sets density
    assert(newrho > 0.);
    
    rho = newrho;
}

void material::set_lambda(const double newlambda) {
    // sets first lame parameter
    
    assert(newlambda > 0.);
    
    lambda = newlambda;
}

void material::set_g(const double newg) {
    // sets shear modulus
    
    assert(newg > 0.);
    
    g = newg;
}

void material::set_mu(const double newmu) {
    // sets plasticity internal friction
    
    assert(newmu > 0.);
    
    mu = newmu;
}

void material::set_beta(const double newbeta) {
    // sets plasticity dilatancy
    
    assert(newbeta > 0.);
    
    beta = newbeta;
}

    
void material::set_eta(const double neweta) {
    // sets plasticity viscosity
    
    assert(neweta > 0.);
    
    eta = neweta;
}

double material::get_rho() const {
    // returns density
    
    return rho;
}

double material::get_lambda() const {
    // returns first lame parameter
    
    return lambda;
}

double material::get_g() const {
    // returns shear modulus
    
    return g;
}

double material::get_mu() const {
    // returns plasticity internal friction
    
    return mu;
}

double material::get_beta() const {
    // returns plastic dilatancy
    
    return beta;
}

double material::get_eta() const {
    // returns plastic viscosity

    return eta;
}

double material::get_cs() const {
    // returns shear wave speed

    return sqrt(g/rho);
}

double material::get_cp() const {
    // returns dilatational wave speed
    
    return sqrt((lambda+2.*g)/rho);
}

double material::get_zs() const {
    // returns shear wave speed
    
    return rho*get_cs();
}

double material::get_zp() const {
    // returns dilatational wave speed
    
    return rho*get_cp();
}
