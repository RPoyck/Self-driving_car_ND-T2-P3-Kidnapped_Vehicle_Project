/*
 * particle_filter.h
 *
 * 2D particle filter class.
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

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

    int id;
    double x;
    double y;
    double theta;
    double weight;
    std::vector<int> associations;
    std::vector<double> sense_x;
    std::vector<double> sense_y;
};



class ParticleFilter {

    // Number of particles to draw
    unsigned int num_particles; 

    // Flag, if filter is initialized
    bool is_initialized;

    // Vector of weights of all particles
    std::vector<double> weights;
    
public:
	
    // Set of current particles
    std::vector<Particle> particles;

    // Constructor
    // @param num_particles Number of particles
    ParticleFilter() : num_particles(0), is_initialized(false) {}


    // Destructor
    ~ParticleFilter() {}


    /**
    * init Initializes particle filter by initializing particles to Gaussian
    *   distribution around first position and all the weights to 1.
    * @param x Initial x position [m] (simulated estimate from GPS)
    * @param y Initial y position [m]
    * @param theta Initial orientation [rad]
    * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m] standard deviation of yaw [rad]]
    */
    void init(double x, double y, double theta, double std[]);


    /**
    * prediction Predicts the state for the next time step
    *   using the process model.
    * @param delta_t Time between time step t and t+1 in measurements [s]
    * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
    *   standard deviation of yaw [rad]]
    * @param velocity Velocity of car from t to t+1 [m/s]
    * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
    */
    void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	
    /**
    * Transform finds the transform of the original point within another reference frame.
    * @param original Pose of the points in the original ref. frame.
    * @param reference Pose of the original ref. frame within the goal ref. frame.
    * @param transformed Empty pose to write the resulting pose to.
    */
//     void Transform(LandmarkObs& original, Particle& reference, LandmarkObs& transformed);
    
    
    /**
    * Transform finds the transform of the original points within another reference frame.
    * @param original vector of poses of the points in the original ref. frame.
    * @param reference vector of poses of the original ref. frame within the goal ref. frame.
    * @param transformed Empty vector of poses to write the resulting poses to.
    */
//     void Transform(std::vector<LandmarkObs> original, Particle reference, std::vector<LandmarkObs>& transformed);
    
    std::vector<LandmarkObs> Transform(std::vector<LandmarkObs> original, Particle reference);

    
    /**
    * Returns the Euclidean distance between two points.
    * @param p_1 Point 1.
    * @param p_2 Point 2.
    */
    double DEucl(LandmarkObs& p_1, LandmarkObs& p_2);
    double DEucl(LandmarkObs& p_1, Map::single_landmark_s& p_2);
    double DEucl(Particle& p_1, LandmarkObs& p_2);
    double DEucl(Particle& p_1, Map::single_landmark_s& p_2);
    
    /**
    * dataAssociation Finds which observations correspond to which landmarks (likely by using
    *   a nearest-neighbors data association).
    * @param predicted Vector of predicted landmark observations
    * @param observations Vector of landmark observations
    */
//     void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
//     void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, Particle& part, double sensor_range);
    void dataAssociation(Map predicted, std::vector<LandmarkObs>& observations, Particle& part, double sensor_range);


    /**
    * updateWeights Updates the weights for each particle based on the likelihood of the 
    *   observed measurements. 
    * @param sensor_range Range [m] of sensor
    * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
    * @param observations Vector of landmark observations
    * @param map Map class containing map landmarks
    */
    void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks);


    /**
    * resample Resamples from the updated set of particles to form
    *   the new set of particles.
    */
    void resample();


    /*
    * Set a particles list of associations, along with the associations calculated world x,y coordinates
    * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
    */
//     Particle SetAssociations(Particle& particle, const std::vector<int>& associations, const std::vector<double>& sense_x, const std::vector<double>& sense_y);


    std::string getAssociations(Particle best);
    std::string getSenseX(Particle best);
    std::string getSenseY(Particle best);


    /**
    * initialized Returns whether particle filter is initialized yet or not.
    */
    const bool initialized() const { return is_initialized; }
	
};



#endif /* PARTICLE_FILTER_H_ */
