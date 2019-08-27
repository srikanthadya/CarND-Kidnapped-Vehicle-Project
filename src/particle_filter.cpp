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
using std::normal_distribution;
using namespace std;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if (is_initialized){
       return;
  }
  num_particles = 100;  // TODO: Set the number of particles
  normal_distribution<double> dist_x(0, std[0]);
  
  // TODO: Create normal distributions for y and theta
  normal_distribution<double> dist_y(0, std[1]);
  
  normal_distribution<double> dist_theta(0, std[2]);
  
  std::default_random_engine gen;

  for (int i = 0; i < num_particles; ++i) {
            Particle p;  
            p.id = i;
            p.weight  = 1.0;
            p.x = dist_x(gen);
            p.y = dist_y(gen);
            p.theta = dist_theta(gen);
            particles.push_back(p);
            }
  is_initialized = true;
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
   normal_distribution<double> dist_x(0, std_pos[0]);

   // TODO: Create normal distributions for y and theta
   normal_distribution<double> dist_y(0, std_pos[1]);

   normal_distribution<double> dist_theta(0, std_pos[2]);
   std::default_random_engine gen;
   for (int i = 0; i < particles.size(); ++i) {
       if (yaw_rate < 0.00001){
            particles[i].x += velocity*delta_t*cos(particles[i].theta);
            particles[i].y += velocity*delta_t*sin(particles[i].theta);
            }
       else {

            particles[i].x += velocity/yaw_rate *(sin(particles[i].theta + yaw_rate*delta_t) -sin(particles[i].theta)) ;
            particles[i].y += velocity/yaw_rate *(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) ;
            particles[i].theta += yaw_rate*delta_t;
            }
       // Adding Noise
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
  double min_dist = std::numeric_limits<double>::max(); 
  double measured_dist;
  for (int i=0;i<observations.size();i++){
    double x__ = observations[i].x;
    double y__ = observations[i].y;
    int closest = -1;
    for (int j=0;j<predicted.size();j++){
      double  measured_dist = dist(x__,y__,predicted[j].x,predicted[j].y);
      
      if ( measured_dist < min_dist) {
        closest = predicted[j].id; 
        min_dist = measured_dist;
      }
    }
    observations[i].id = closest;    
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
   *   My understanding is each observation is a ground truth measurement and each particle is the predicted vehicle position
   */
  
  double norm_weight = 0.0;
  for (int j=0;j<particles.size();j++){
    double prt_x = particles[j].x;
    double prt_y = particles[j].y;
    double prt_theta = particles[j].theta;
    vector<LandmarkObs> predictions;
    vector<LandmarkObs> transformed_obs;
    
    // Collecting only those landmarks within sensor range
    for (int k=0;k<map_landmarks.landmark_list.size();k++){
      
      double x_ = map_landmarks.landmark_list[k].x_f;
      double y_ = map_landmarks.landmark_list[k].y_f;
      
      
      if (dist(x_,y_,prt_x,prt_y) < sensor_range){
         int id_ = map_landmarks.landmark_list[k].id_i;
         predictions.push_back(LandmarkObs{id_,x_,y_});
      }
    }
      
    for (int i=0;i<observations.size();i++){
      LandmarkObs tobs;
      tobs.id = observations[i].id;
      tobs.x = prt_x + (observations[i].x * cos(prt_theta)) - observations[i].y * sin(prt_theta);
      tobs.y = prt_y + (observations[i].x * sin(prt_theta)) + observations[i].y * cos(prt_theta);
      transformed_obs.push_back(tobs);
    }
    dataAssociation(predictions,transformed_obs);
    //particles[j].weight = 1.0;
    for (int l=0;l < transformed_obs.size();l++){
      double pred_x,pred_y,obs_w;
      double tobs_x = transformed_obs[l].x;
      double tobs_y = transformed_obs[l].y;
      for (int m=0;m < predictions.size();m++){        
        if (predictions[m].id == transformed_obs[l].id){
          pred_x = predictions[m].x;
          pred_y = predictions[m].y;
        }
      }
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      //double obs_w =  (1./(2.*M_PI*s_x*s_y)) * exp( -1.*(pow(pred_x-tobs_x,2)/(2*pow(s_x, 2)) + (pow((pred_y-tobs_y,2)/(2*pow(s_y, 2))));
      obs_w = ( 1./(2.*M_PI*s_x*s_y)) * exp( -( pow(tobs_x-pred_x,2)/(2*pow(s_x, 2)) + (pow(tobs_y-pred_y,2)/(2*pow(s_y, 2))) ) );                                       
      //if ( obs_w == 0) {
        //particles[j].weight *= 0.00001;
      //}
      //else{

      particles[j].weight *= obs_w;
      //std::cout << "obs_w " << obs_w << std::endl;
      //}
    }

  
  norm_weight += particles[j].weight;
  
  }
  for (int i=0;i<particles.size();i++){
    //particles[i].weight /= norm_weight;
    weights.push_back(particles[i].weight);
  }
}
  
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  vector<Particle> resampled_particles;
  double max_wt = 2* *std::max_element(weights.begin(),weights.end());
  std::uniform_int_distribution<int> distribution1(0,particles.size()-1);
  std::uniform_real_distribution<double> distribution2(0.0,1.0);
  std::default_random_engine gen;
  int index = distribution1(gen);
  double beta = 0.0;
  for (int i=0;i<particles.size();i++){
    beta += distribution2(gen) * max_wt;
    while ( beta > weights[index]){
      beta -= weights[index];
      index = (index+1) % particles.size();
    }
    resampled_particles.push_back(particles[index]);
  }
  particles = resampled_particles;

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
