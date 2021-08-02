/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::sin;
using std::cos;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  std::random_device rd{};
  std::mt19937_64 gen{ rd() };
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i)
  {
	  Particle p;
	  p.x = dist_x(gen);
	  p.y = dist_y(gen);
	  p.theta = dist_theta(gen);
	  p.weight = 1.0;
	  particles.push_back(p);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) 
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	std::random_device rd{};
	std::mt19937_64 gen{ rd() };
	std::normal_distribution<double> dist_x(0.0, std_pos[0]);
	std::normal_distribution<double> dist_y(0.0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0.0, std_pos[2]);

	for (int i = 0; i < num_particles; ++i)
	{
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

    // TODO: Guard when yaw_rate is 0.
		double newX = dist_x(gen) + x + (velocity / yaw_rate) * (sin(theta + yaw_rate*delta_t) - sin(theta));
		double newY = dist_y(gen) + y + (velocity / yaw_rate) * (cos(theta) + cos(theta + yaw_rate * delta_t));
		double newTheta = dist_theta(gen) + theta + yaw_rate * delta_t;
		particles[i].x = newX;
		particles[i].y = newY;
		particles[i].theta = newTheta;
	}
}

void ParticleFilter::dataAssociation(
	vector<LandmarkObs> predicted, 
    vector<LandmarkObs>& observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	for (int i = 0; i < observations.size(); ++i)
	{
		LandmarkObs obs = observations[i];
		LandmarkObs closest = *std::min_element
		(
			predicted.begin(), 
			predicted.end(),
			[&obs]
			(
				const LandmarkObs& obs1, 
				const LandmarkObs& obs2
			) 
			{
				return dist(obs.x, obs.y, obs1.x, obs1.y) < dist(obs.x, obs.y, obs2.x, obs2.y);
			}
		);
		observations[i].id = closest.id;
	}
}

void ParticleFilter::updateWeights(
	double sensor_range, 
	double std_landmark[], 
    const vector<LandmarkObs> &observations, 
	const Map &map_landmarks) 
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  vector<LandmarkObs> landmarkObs;
  toLandmarkObs(map_landmarks, landmarkObs);

  for (int i = 0; i < particles.size(); i++) 
  {
    Particle p = particles[i];

    // mapBasedObs are observations in map/world coordinates
    vector<LandmarkObs> mapBasedObs;
    car2MapCoordinateTransform(p, observations, mapBasedObs);

    // For each observation find the closest landmark and copy its landmark id
    dataAssociation(landmarkObs, mapBasedObs);

	  double weight = 1.0;
    // Compute the prob for each observation
    for (auto obs : mapBasedObs)
    {
      double muX = obs.x;
      double muY  = obs.y;
      double x = map_landmarks.landmark_list[obs.id-1].x_f;
      double y = map_landmarks.landmark_list[obs.id-1].y_f;
      double sigmaX = std_landmark[0];
      double signmaY = std_landmark[1];
      weight *= multiv_prob(sigmaX, signmaY, x, y, muX, muY);
    }
    
	  p.weight = weight;
  }
}

void ParticleFilter::car2MapCoordinateTransform(
    const Particle& p, 
    const vector<LandmarkObs> &observations, 
    vector<LandmarkObs>& mapBasedObs)
{
  mapBasedObs.clear();
  for (auto obs : observations) 
  {
    LandmarkObs transformed;
    transformed.id = obs.id;
    transformed.x =  p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);
    transformed.y =  p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
    mapBasedObs.push_back(transformed);
  }
}

void ParticleFilter::toLandmarkObs(const Map& map, std::vector<LandmarkObs>& landmarkObs)
{
  landmarkObs.clear();
  for (auto lm : map.landmark_list)
  {
    LandmarkObs obs;
    obs.id = lm.id_i;
    obs.x = lm.x_f;
    obs.y = lm.y_f;
    landmarkObs.push_back(obs);
  }
}

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<double> weights;
	for (Particle p : particles)
	{
		weights.push_back(p.weight);
	}
	std::discrete_distribution<double> d(weights.begin(), weights.end());

	std::vector<Particle> newParticles;
	for (int i = 0; i < num_particles; ++i)
	{
		int idx = d(gen);
		newParticles.push_back(particles[idx]);
	}

	particles = newParticles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
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
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) 
{
  vector<double> v;

  if (coord == "X") 
  {
    v = best.sense_x;
  } 
  else 
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

double ParticleFilter::multiv_prob(
	double sig_x, 
	double sig_y, 
	double x_obs, 
	double y_obs,
	double mu_x, 
	double mu_y) 
{
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}