/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
    
    // Set the number of particles. //
    this->num_particles = 1000;
    
    double sigma_x = std[0];
    double sigma_y = std[1];
    double sigma_theta = std[2];
    
    normal_distribution<double> dist_x(x, sigma_x);
    normal_distribution<double> dist_y(y, sigma_y);
    normal_distribution<double> dist_theta(theta, sigma_theta);
    
    default_random_engine gen;
    
    for (int i = 0; i < this->num_particles; i++)
    {
	// Initialize all particles to first position (based on estimates of x, y, theta and their uncertainties from GPS). //
	// Add random Gaussian noise to each particle. //
	Particle par;
	par.id = i;
	par.x = dist_x(gen);
	par.y = dist_y(gen);
	par.theta = dist_theta(gen);
	// Initialize all particles weights to 1. //
	par.weight = 1.0;
	this->weights.push_back(par.weight);
	
	this->particles.push_back(par);
    }
    
    this->is_initialized = true;

}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
    
    double sigma_x = std_pos[0];
    double sigma_y = std_pos[1];
    double sigma_theta = std_pos[2];
    
    default_random_engine gen;
    
    // Add measurements to each particle and add random Gaussian noise. //
    for (int i = 0; i < this->num_particles; i++)
    {
	if (abs(yaw_rate) < pow(10,-6) )
	{
	    // Add measurement to x position //
	    particles[i].x += velocity*delta_t*cos(particles[i].theta);
	    // Add measurement to y position //
	    particles[i].y += velocity*delta_t*sin(particles[i].theta);
	}
	else
	{
	    // Add measurement to x position //
	    particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate*delta_t)) - sin(particles[i].theta));
	    // Add measurement to y position //
	    particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate*delta_t)));
	    // Add measurement to orientation //
	    particles[i].theta += yaw_rate*delta_t;
	}
	
	// Add random Gaussian noise. //
	normal_distribution<double> dist_x(particles[i].x, sigma_x);
	particles[i].x = dist_x(gen);
	normal_distribution<double> dist_y(particles[i].y, sigma_y);
	particles[i].y = dist_y(gen);
	normal_distribution<double> dist_theta(particles[i].theta, sigma_theta);
	particles[i].theta = dist_theta(gen);
    }
    
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

}


void ParticleFilter::Transform(LandmarkObs& original, Particle& reference, LandmarkObs& transformed)
{
    
    transformed.id = original.id;

    // Transform x coordinate //
    transformed.x = reference.x + (cos(reference.theta) * original.x) - (sin(reference.theta) * original.y);

    // transform y coordinate //
    transformed.y = reference.y + (sin(reference.theta) * original.x) + (cos(reference.theta) * original.y);

}


double ParticleFilter::DEucl(LandmarkObs& p_1, LandmarkObs& p_2)
{ 
    return sqrt( pow(p_1.x-p_2.x,2) + pow(p_1.y-p_2.y,2)); 
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) 
{
    unsigned int num_obs = observations.size();
    unsigned int num_pred = observations.size();
    unsigned int  best_match;
    double d_min;
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the observed measurement to this particular landmark. //
    for (int i=0; i < num_obs; i++)
    {
	d_min = DEucl(observations[i], predicted[0]);
	best_match = predicted[0].id;
	
	for (int i2=1; i2 < num_pred; i2++)
	{
	    double d = DEucl(observations[i], predicted[i2]);
	    if (d < d_min)
	    {
		d_min = d;
		best_match = predicted[0].id;
	    }
	}
	observations[i].id = best_match;
    }
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to implement this method and use it as a helper during the updateWeights phase.

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
    
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation 
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
    
//     double x = 6;
//     double y = 3;
//     double mu_x = 5;
//     double mu_y = 3;
//     double sigma_x = 0.3;
//     double sigma_y = 0.3;
// 
//     // calculate normalisation term //
//     double gauss_norm= (1/(2 * M_PI * sigma_x * sigma_y));
//     // calculate exponent //
//     double exponent= ((x - mu_x)**2)/(2 * sigma_x**2) + ((y - mu_y)**2)/(2 * sigma_y**2);
//     // calculate weight using normalisation terms and exponent //
//     double weight = gauss_norm * exp(exponent);
// 
//     double P = (1/(2*M_PI*sigma_x*sigma_y)) * exp(-( (pow(x-mu_x,2) / (2*sigma_x*sigma_x) ) + (pow(y-mu_y,2) / (2*sigma_y*sigma_y) ) ));
	
    
    
}


void ParticleFilter::resample() 
{
    
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}


Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
}


string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
