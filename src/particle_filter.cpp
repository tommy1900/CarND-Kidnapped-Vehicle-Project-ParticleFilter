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
#include "particle_filter.h"
#include "helper_functions.h"
using namespace std;

using std::string;
using std::vector;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  // init sensor noise
  normal_distribution<double> dist_x_init(0,std[0]);
  normal_distribution<double> dist_y_init(0,std[1]);
  normal_distribution<double> dist_theta_init(0,std[2]);

  // init particles
  Particle p;
  for (int i = 0; i<num_particles; i++){
    p.id = i;
    p.x = x + dist_x_init(gen);
    p.y = y + dist_y_init(gen);
    p.theta = theta + dist_theta_init(gen);
    p.weight = 1.0;

    //save into the current particle vector
    particles.push_back(p);
  }
  is_initialized = 1;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   normal_distribution<double> dist_x(0,std_pos[0]);
   normal_distribution<double> dist_y(0,std_pos[1]);
   normal_distribution<double> dist_theta(0,std_pos[2]);

   for (int i = 0; i < num_particles; i++){
     // predict the state
     if (fabs(yaw_rate) < 0.00001){ //check if yew rate is zero
       particles[i].x += velocity * delta_t * cos(particles[i].theta);
       particles[i].y += velocity * delta_t * sin(particles[i].theta);
     }
     else {
       particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
       particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
       particles[i].theta += yaw_rate * delta_t;
     }
    // with the noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
   }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
   // observation loop
   for (int i = 0; i<observations.size(); i++){
     // init the min_dist with the max distance possible
     double min_dist = numeric_limits<double>::max();
     // prediction loop O(mn)
     int associated_id = -1; // init the id
     for(int j = 0; j<predicted.size();j++){
       double cur_dist = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
       if(cur_dist < min_dist){
         min_dist = cur_dist;
         associated_id = predicted[j].id;
       }
     }
     observations[i].id = associated_id;
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
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

   //run thro all particles
   for(int i = 0; i < num_particles; i++){
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    //check for nearby landmarks
    vector<LandmarkObs> nearLandmarks;

    for(int j = 0; j < map_landmarks.landmark_list.size();j++){
      double landmarkx = map_landmarks.landmark_list[j].x_f;
      double landmarky = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      if(fabs(landmarkx - x) <= sensor_range && fabs(landmarky - y) <= sensor_range){
        nearLandmarks.push_back(LandmarkObs{ id, landmarkx, landmarky });
      }
    }

    // transform the observation
    vector<LandmarkObs> transObs;
    for(int j = 0; j < observations.size(); j++){
      double map_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double map_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      transObs.push_back(LandmarkObs{ observations[j].id, map_x, map_y });
    }
    // obs association with the landmarks
    dataAssociation(nearLandmarks,transObs); // now we know if the obs match up with the inranged lanmarks
    // reinit weight
    particles[i].weight = 1.0;

    //calculate for new weight
    for(auto mapObs:transObs){
      double pred_x, pred_y;
      for(auto mark:nearLandmarks){
        if(mark.id == mapObs.id){
          pred_x = mark.x;
          pred_y = mark.y;
          // if the inranged landmarks indeed match with the observations
          // we claim the prediction x y as the landmarked truth
        }
      }
      double noise_x = std_landmark[0];
      double noise_y = std_landmark[1];
      double obs_w = ( 1/(2*M_PI*noise_x*noise_y)) * exp( -( pow(pred_x-mapObs.x,2)/(2*pow(noise_x, 2)) + (pow(pred_y-mapObs.y,2)/(2*pow(noise_y, 2))) ) );
      // with some uncertainty, we can deliver the results:
      // how close the ith particle towards the identified lanmark
      particles[i].weight *= obs_w;
    }
   }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   //create a temp holder
   vector<Particle> new_particles;
   //w_o
   vector<double> w;
   for(int i = 0; i<num_particles; i++){
     w.push_back(particles[i].weight);
   }

   // use resampling wheel
   uniform_int_distribution<int> uniintdist(0, num_particles-1); //random engine among particles
   int index = uniintdist(gen);

   double beta = 0.0;
   double mw = *max_element(w.begin(),w.end());

   uniform_real_distribution<double> unidoubledist(0.0, mw); //random engine among particles

   for (int i = 0; i< num_particles; i++){
     beta += unidoubledist(gen) * 2.0;
     while(beta > w[index]){
       beta -= w[index];
       index = (index+1)%num_particles;
     }
     new_particles.push_back(particles[index]);
   }
   particles = new_particles;

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

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
