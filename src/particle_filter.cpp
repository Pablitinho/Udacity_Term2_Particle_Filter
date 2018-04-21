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
#include <limits>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles=150;
	weights.resize(num_particles);

	for (int i=0;i<num_particles;i++)
	{
		Particle particle_aux;
		particle_aux.id=i;

		particle_aux.x = dist_x(gen);
		particle_aux.y = dist_y(gen);
		particle_aux.theta = dist_theta(gen);
		particle_aux.weight = 1.0f;

		particles.push_back(particle_aux);

	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);

	//Avoid division by zero
	//TODO: Fix this issue
	if(fabsf(yaw_rate)<0.01)
	{
	   yaw_rate=0.001;
	}
	// Estimate the displacement
	for (int i=0;i<num_particles;i++)
	{
        particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
        particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
        particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
        //Debug
        if (isnan(particles[i].x))
        {
        	std::cout<<std::endl<<"yaw_rate: "<< yaw_rate;
        }
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	// NOT NEEDED
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	// Homogenous transformation
	float trans_x=0.0f;
	float trans_y=0.0f;
	float min_distance = std::numeric_limits<float>::max();;
	float distance = 0.0f;
	int id_near_lan=0;
	float weight=1.0;
	float dif_x=0.0f;
	float dif_y=0.0f;

	// Estimate std
	float sigma_xx = std_landmark[0]*std_landmark[0];
	float sigma_yy = std_landmark[1]*std_landmark[1];
	float sigma_xy = std_landmark[0]*std_landmark[1];

	for (int id_part = 0; id_part < num_particles; id_part++)
	{
		weight=1.0f;
		for (uint id_obs = 0; id_obs < observations.size(); id_obs++)
		{
			//Homogenous transformation VCS->MAP
			trans_x = observations[id_obs].x * cos(particles[id_part].theta) - observations[id_obs].y * sin(particles[id_part].theta) + particles[id_part].x;
			trans_y = observations[id_obs].x * sin(particles[id_part].theta) + observations[id_obs].y * cos(particles[id_part].theta) + particles[id_part].y;

			// Get the nearest landmark
			min_distance = std::numeric_limits<float>::max();

			for (uint id_lan = 0; id_lan < map_landmarks.landmark_list.size(); id_lan++)
			{
				distance = fabs(trans_x - map_landmarks.landmark_list[id_lan].x_f) + fabs(trans_y - map_landmarks.landmark_list[id_lan].y_f);

				if (distance<min_distance)
				{
					min_distance=distance;
					id_near_lan = id_lan;
				}
			}
			dif_x = trans_x- map_landmarks.landmark_list[id_near_lan].x_f;
			dif_y = trans_y- map_landmarks.landmark_list[id_near_lan].y_f;

			weight *= exp(-0.5*((dif_x * dif_x)/sigma_xx + (dif_y * dif_y)/sigma_yy))/(2*M_PI*sigma_xy);

		}

		if (weight<0.0001)
		{
			weight=0.0001;
		}

		//Assign the new weights
		particles[id_part].weight = weight;
		weights[id_part] = weight;

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<size_t> dist_index(weights.begin(), weights.end());

	vector<Particle> resampled_particles(particles.size());

	for (uint id_par = 0; id_par < particles.size(); id_par++)
	{
	  resampled_particles[id_par] = particles[dist_index(gen)];
	}

	particles = resampled_particles;

	/*Wheel version. Working worse*/
/*
	  // Get max weight.

	  float maxWeight = numeric_limits<float>::min();
	  for(int id_par = 0; id_par < num_particles; id_par++)
	  {
	   // weights.push_back(particles[i].weight);
	    if ( particles[id_par].weight > maxWeight )
	    {
	         maxWeight = particles[id_par].weight;
	    }
	  }

	  // Creating distributions.
	  uniform_real_distribution<float> dist_float(0.0, maxWeight);
	  uniform_int_distribution<int> distInt(0, num_particles - 1);
	  default_random_engine gen;

	  // Generating index.
	  int index = distInt(gen);

	  double beta = 0.0;

	  // Wheel algorithm
	  vector<Particle> resampledParticles;
	  for(int i = 0; i < num_particles; i++)
	  {
	    beta += dist_float(gen) * 2.0;
	    while( beta > weights[index])
	    {
	      beta -= weights[index];
	      index = (index + 1) % num_particles;
	    }
	    resampledParticles.push_back(particles[index]);
	  }

	  particles = resampledParticles;
	  */
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
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
