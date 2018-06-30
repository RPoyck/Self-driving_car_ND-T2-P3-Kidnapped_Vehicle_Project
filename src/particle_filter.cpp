/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
    
    std::cout << "!! Vehicle x: " << x << "\n";
    std::cout << "!! Vehicle y: " << y << "\n";
    std::cout << "!! Vehicle theta: " << theta << "\n";
    
    // Set the number of particles. //
    this->num_particles = 100;
    
    double sigma_x = std[0];
    double sigma_y = std[1];
    double sigma_theta = std[2];
    
    normal_distribution<double> dist_x(x, sigma_x);
    normal_distribution<double> dist_y(y, sigma_y);
    normal_distribution<double> dist_theta(theta, sigma_theta);
    
    default_random_engine gen;
    std::vector<double> weights_t;
    std::vector<Particle> particles_t;
    for (unsigned int i = 0; i < this->num_particles; i++)
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
	weights_t.push_back(par.weight);
	
	particles_t.push_back(par);
    }
    
    this->weights = weights_t;
    this->particles = particles_t;
    
    this->is_initialized = true;
    
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
    
    double sigma_x = std_pos[0];
    double sigma_y = std_pos[1];
    double sigma_theta = std_pos[2];
    
    normal_distribution<double> dist_x(0, sigma_x);
    normal_distribution<double> dist_y(0, sigma_y);
    normal_distribution<double> dist_theta(0, sigma_theta);
    
    default_random_engine gen;
    
    // Add measurements to each particle and add random Gaussian noise. //
    for (unsigned int i = 0; i < this->num_particles; i++)
    {
	if (fabs(yaw_rate) < 0.000001 )
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
	particles[i].x += dist_x(gen);
	particles[i].y += dist_y(gen);
	particles[i].theta += dist_theta(gen);
	
    }
    
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

}


// void ParticleFilter::Transform(LandmarkObs& original, Particle& reference, LandmarkObs& transformed)
// {
//     
//     transformed.id = original.id;
// 
//     // Transform x coordinate //
//     transformed.x = reference.x + (cos(reference.theta) * original.x) - (sin(reference.theta) * original.y);
// 
//     // transform y coordinate //
//     transformed.y = reference.y + (sin(reference.theta) * original.x) + (cos(reference.theta) * original.y);
// 
// }
// 
// 
// void ParticleFilter::Transform(std::vector<LandmarkObs> original, Particle reference, std::vector<LandmarkObs>& transformed)
// {
//     unsigned int num_originals = original.size();
//     for (unsigned int i = 0; i < num_originals; i++)
//     {
// 	transformed[i].id = original[i].id;
// 
// 	// Transform x coordinate //
// 	transformed[i].x = reference.x + (cos(reference.theta) * original[i].x) - (sin(reference.theta) * original[i].y);
// 
// 	// transform y coordinate //
// 	transformed[i].y = reference.y + (sin(reference.theta) * original[i].x) + (cos(reference.theta) * original[i].y);
//     }
// }


std::vector<LandmarkObs> ParticleFilter::Transform(std::vector<LandmarkObs> original, Particle reference)
{
    unsigned int num = original.size();
    std::vector<LandmarkObs> transformed(num);
    unsigned int num_originals = original.size();
    for (unsigned int i = 0; i < num_originals; i++)
    {
	transformed[i].id = original[i].id;

	// Transform x coordinate //
	transformed[i].x = reference.x + (cos(reference.theta) * original[i].x) - (sin(reference.theta) * original[i].y);

	// transform y coordinate //
	transformed[i].y = reference.y + (sin(reference.theta) * original[i].x) + (cos(reference.theta) * original[i].y);
	
    }
    
    return transformed;
}


double ParticleFilter::DEucl(LandmarkObs& p_1, LandmarkObs& p_2)
{ 
    return sqrt( pow(p_1.x-p_2.x,2) + pow(p_1.y-p_2.y,2)); 
}


double ParticleFilter::DEucl(LandmarkObs& p_1, Map::single_landmark_s& p_2)
{ 
    return sqrt( pow(p_1.x-p_2.x_f,2) + pow(p_1.y-p_2.y_f,2)); 
}


double ParticleFilter::DEucl(Particle& p_1, LandmarkObs& p_2)
{ 
    return sqrt( pow(p_1.x-p_2.x,2) + pow(p_1.y-p_2.y,2)); 
}


double ParticleFilter::DEucl(Particle& p_1, Map::single_landmark_s& p_2)
{
    return sqrt( pow(p_1.x-p_2.x_f,2) + pow(p_1.y-p_2.y_f,2)); 
}


void ParticleFilter::dataAssociation(Map predicted, std::vector<LandmarkObs>& observations, Particle& part, double sensor_range) 
{
    unsigned int num_obs = observations.size();
    unsigned int num_pred = predicted.landmark_list.size();
    
    // Find the predicted measurement that is closest to each observed measurement and assign the observed measurement to this particular landmark. //
    for (unsigned int i=0; i < num_obs; i++)
    {
	int  best_match = -1;
	double d_min = pow(10.0, 6);
	
	for (unsigned int i2=0; i2 < num_pred; i2++)
	{
// 	    if (sensor_range < DEucl(part, predicted.landmark_list[i2])) {continue;}
	    if (sensor_range < dist(part.x, part.y, predicted.landmark_list[i2].x_f, predicted.landmark_list[i2].y_f)) {continue;}
	    
// 	    double d = DEucl(observations[i], predicted.landmark_list[i2]);
	    double d = dist(observations[i].x, observations[i].y, predicted.landmark_list[i2].x_f, predicted.landmark_list[i2].y_f);
	    
	    if (d < d_min)
	    {
		d_min = d;
		best_match = i2;
	    }
	}
	observations[i].id = best_match;
    }
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to implement this method and use it as a helper during the updateWeights phase.
}	
	

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
    
    // Update the weights of each particle using a multi-variate Gaussian distribution. You can read more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution //
    
    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    // calculate normalisation term for multi-variate Gaussian distribution //
    double gauss_norm= 1/(2 * M_PI * sigma_x * sigma_y);
    
    for (unsigned int i=0; i < this->num_particles; i++)
    {
	// Transform all observations for this specific particle to map coordinate system //
	std::vector<LandmarkObs> observations_t = Transform(observations, this->particles[i]);
	
	// Associate the predicted landmarks for this particle with the observed landmarks //
	dataAssociation(map_landmarks, observations_t, this->particles[i], sensor_range);
	
	// Calculate the weight of the particle using multi-variate Gaussian distribution //
	this->particles[i].weight = 1.0;
	unsigned int num_obs = observations_t.size();
	for (unsigned int i2=0; i2 < num_obs; i2++)
	{
	    if (observations_t[i2].id == -1) {continue;}
	    double x = observations_t[i2].x;
	    double y = observations_t[i2].y;
	    double mu_x = map_landmarks.landmark_list[observations_t[i2].id].x_f;
	    double mu_y = map_landmarks.landmark_list[observations_t[i2].id].y_f;

	    // calculate exponent //
	    double exponent= ( pow((x - mu_x),2) / (2 * pow(sigma_x, 2)) ) + ( pow((y - mu_y), 2) / (2 * pow(sigma_y, 2)) );
	    // calculate weight using normalisation terms and exponent //
	    double weight = gauss_norm * exp(-1.0 * exponent);
	    
	    this->particles[i].weight *= weight;
	}
	this->weights[i] = this->particles[i].weight;
	
    }
    
    /* NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
      according to the MAP'S coordinate system. You will need to transform between the two systems.
      Keep in mind that this transformation requires both rotation AND translation (but no scaling).
      The following is a good resource for the theory:
      https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
      and the following is a good resource for the actual equation to implement (look at equation 
      3.33
      http://planning.cs.uiuc.edu/node99.html */
    
}


// Resample particles with replacement with probability proportional to their weight. //
void ParticleFilter::resample() 
{

    double max_w = this->weights[*max_element(std::begin(this->weights), std::end(this->weights))];
    
    // If the maximally weighted particle has a weight of (near) zero, all particles have probably drifted from the vehicle location too much and are all reinitialised //
    if (max_w < pow(10.0, -250.0)) {this->is_initialized = false;}
    
    std::default_random_engine re;
    discrete_distribution<int> unif(0.0, (2*max_w));
    
    double beta = 0.0;
    unsigned int i2 = 0;
    std::vector<Particle> particles_t;
    
    for (unsigned int i=0; i < this->num_particles; i++)
    {
	// Add a random number to beta between 0 and 2 times the max weight //
	// Subtract the weight of each particle from this beta until beta is smaller than the weight of the current particle //
	// Take this particle into the new particle list // 
	
	beta += unif(re);
	
	while (beta > this->weights[i2])
	{
	    beta -= this->weights[i2];
	    i2 = (i2+1) % this->num_particles;
	}
	particles_t.push_back(this->particles[i2]);
    }
    this->particles = particles_t;
	
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}


// Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, const std::vector<double>& sense_x, const std::vector<double>& sense_y)
// {
//     
//     //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
//     // associations: The landmark id that goes along with each listed association
//     // sense_x: the associations x mapping already converted to world coordinates
//     // sense_y: the associations y mapping already converted to world coordinates
// 
//     particle.associations= associations;
//     particle.sense_x = sense_x;
//     particle.sense_y = sense_y;
//     
// }


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
