// File:          robot01.cpp
// Date:
// Description:
 // Author:
// Modifications:

// You may need to add webots include files such as
// <webots/DistanceSensor.hpp>, <webots/Motor.hpp>, etc.
// and/or to add some other includes
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/Keyboard.hpp>
#include <webots/PositionSensor.hpp>
#include <iostream>
#include <cstdint>

#define MAX_VELOCITY 5
#define MIN_VELOCITY -5
#define MAX_GEAR 5
#define MIN_GEAR -2

// All the webots classes are defined in the "webots" namespace
using namespace webots;

class Gripper{
  public:
    Gripper(Motor **fingers){ // right finger to be given as the first one
      this->fingers = fingers;
    }
    void set_velocities(){
      for(int i = 0; i < 2; i++){
        fingers[i]->setVelocity(0.1);
      }
    }
    bool is_open(){}
    bool is_closed(){}
    void close(){
      fingers[0]->setPosition(fingers[0]->getMaxPosition());
      fingers[1]->setPosition(fingers[1]->getMinPosition());
     this->set_velocities();
    }

    void open(){
       fingers[0]->setPosition(fingers[0]->getMinPosition());
       fingers[1]->setPosition(fingers[1]->getMaxPosition());
       this->set_velocities();
    }
    void release(){}
    void tighten(){
      double right_finger_pos = fingers[0]->getPositionSensor()->getValue();
      double left_finger_pos = fingers[1]->getPositionSensor()->getValue();
      std::cout << right_finger_pos << std::endl << left_finger_pos<<std::endl;
      fingers[0]->setPosition(right_finger_pos+0.01);
      fingers[1]->setPosition(left_finger_pos-0.01);
      this->set_velocities();
    }

  private:
    Motor **fingers;
};

float gear_velocities[8] = {-1,-0.5,0,0.3,0.5,0.7,1,1.2};

 double max(double a,double b){
   if(a > b) return a;
   return b;
 };
 
 double min(double a, double b){
   if(a < b) return a;
   return b;
 }
 //depracted
 int equalise(double velocity){
   if(velocity > MAX_VELOCITY) return MAX_VELOCITY;
   else if(velocity < MIN_VELOCITY) return MIN_VELOCITY;
   else return (int)velocity;
 }
 
 void set_angular_velocities(int gear,unsigned int turn,double& velocity_right, double& velocity_left){
   int index = gear+2;
    if(turn == 1){
      velocity_left = gear_velocities[index];
      velocity_right = gear_velocities[index]/2;  
    }
    else if(turn == 2){
      velocity_right = gear_velocities[index];
      velocity_left = gear_velocities[index]/2;  
    }
    else{
      velocity_right = gear_velocities[index];
      velocity_left = gear_velocities[index];  
    }
 }
 
 int set_arm1_velocity(unsigned int turn){
   if(turn == 0) return 0;
   else if(turn == 1) return 2;
   return -2;
 }
 
 int set_joint1_velocity(unsigned int turn){
   if(turn == 0) return 0;
   else if(turn == 1) return 2;
   return -2;
 }
  

// This is the main program of your controller.
// It creates an instance of your Robot instance, launches its
// function(s) and destroys it at the end of the execution.
// Note that only one instance of Robot should be created in
// a controller program.
// The arguments of the main function can be specified by the
// "controllerArgs" field of the Robot node
int main(int argc, char **argv) {
  // CREATE INSTANCES OF WEBOTS NODES/OBJECTS
  Keyboard kb;
  Robot *robot = new Robot();
  int timeStep = (int)robot->getBasicTimeStep();
  Motor *wheels[2];
  Motor *fingers[2];
  Motor *arm1;
  Motor *joint1;
  Motor *joint2;
  Motor *joint3;
  
  
  char finger_names[2][16] = {"right_finger","left_finger"};
  char wheels_names[2][8] = {"right","left"}; // "left/0"
  double velocity_left = 0.0;
  double velocity_right = 0.0;
  double velocity_arm1 = 0.0;
  double velocity_joint1 =0.0;
  double velocity_joint2 =0.0;
  double velocity_joint3 =0.0;
  unsigned int turn = 0;
  unsigned int arm1_turn = 0;
  unsigned int joint1_turn = 0;
  unsigned int joint2_turn = 0;
  unsigned int joint3_turn = 0;
  double turn_time_start = 0;
  double allowed_turn_time = 0.010;
  double arm1_time_start = 0;
  double joint1_time_start = 0;
  double joint2_time_start = 0;
  double joint3_time_start = 0;
  int gear = 0;
    

  
  for(int i = 0 ; i < 2 ; i++){
    wheels[i] = robot->getMotor(wheels_names[i]);
    wheels[i]->setPosition(INFINITY);
    wheels[i]->setVelocity(0.0);
  }
  
  for(int i = 0; i < 2; i++){
    fingers[i] = robot->getMotor(finger_names[i]);
    fingers[i]->setPosition(INFINITY);
    fingers[i]->setVelocity(0.0);
    fingers[i]->getPositionSensor()->enable(timeStep);
    std::cout << fingers[i]->getPositionSensor();
    std::cout << robot->getPositionSensor("glf_sensor");
  }
  
  Gripper grip(fingers);
  std::cout << "XD";
  
  arm1 = robot->getMotor("arm1");
  arm1->setPosition(INFINITY);
  arm1->setVelocity(0.0);
  
  joint1 = robot->getMotor("joint1");
  joint1->setPosition(INFINITY);
  joint1->setVelocity(0.0);
  
  joint2 = robot->getMotor("joint2");
  joint2->setPosition(INFINITY);
  joint2->setVelocity(0.0);
  
  joint3 = robot->getMotor("joint3");
  joint3->setPosition(INFINITY);
  joint3->setVelocity(0.0);
  
  // get the time step of the current world.

  kb.enable(timeStep);
  // You should insert a getDevice-like function in order to get the
  // instance of a device of the robot. Something like:
  //  Motor *motor = robot->getMotor("motorname");
  //  DistanceSensor *ds = robot->getDistanceSensor("dsname");
  //  ds->enable(timeStep);

  // Main loop:
  // - perform simulation steps until Webots is stopping the controller
  while (robot->step(timeStep) != -1) {
    if(robot->getTime() >= turn_time_start + allowed_turn_time && turn != 0)
      turn = 0;
    
    if(robot->getTime() >= arm1_time_start + allowed_turn_time && arm1_turn != 0)
      arm1_turn = 0;  
    
    if(robot->getTime() >= joint1_time_start + allowed_turn_time && joint1_turn != 0)
      joint1_turn = 0;
      
     if(robot->getTime() >= joint2_time_start + allowed_turn_time && joint2_turn != 0)
      joint2_turn = 0;
      
     if(robot->getTime() >= joint3_time_start + allowed_turn_time && joint3_turn != 0)
      joint3_turn = 0;
     
          
    int key = kb.getKey();
    if(key == kb.UP){
      gear++;
      if(gear > MAX_GEAR)
        gear = MAX_GEAR;
    }
    else if(key == kb.RIGHT){
      turn = 1;
      turn_time_start = robot->getTime();
    }
    else if(key == kb.LEFT){
      turn = 2;
      turn_time_start = robot->getTime();
    }
    else if(key == kb.DOWN){
       gear--;
       if(gear < MIN_GEAR)
          gear = MIN_GEAR;
    }
    else if(key == 'P'){ // press s to stop
      gear = 0;
    }
    
    else if(key == 'E'){
      arm1_turn = 2;
      arm1_time_start = robot->getTime();
    }
    
    else if(key == 'Q'){
      arm1_turn = 1;
      arm1_time_start = robot->getTime();
    }
    
     else if(key == 'S'){
      joint1_turn = 2;
      joint1_time_start = robot->getTime();
    }
    
    else if(key == 'W'){
      joint1_turn = 1;
      joint1_time_start = robot->getTime();
      }
      
        else if(key == '1'){
      joint2_turn = 2;
      joint2_time_start = robot->getTime();
    }
    
    else if(key == '2'){
      joint2_turn = 1;
      joint2_time_start = robot->getTime();
      }
         else if(key == '3'){
      joint3_turn = 2;
      joint3_time_start = robot->getTime();
    }
    
    else if(key == '4'){
      joint3_turn = 1;
      joint3_time_start = robot->getTime();
      }
    else if(key == 'G'){
       grip.close();
     }
     else if(key == 'O'){
       grip.open();
     }
     else if(key == 'T'){
       grip.tighten();
     }


      
    // Enter here functions to send actuator commands, like:
    //  motor->setPosition(10.0);
    
    set_angular_velocities(gear,turn,velocity_right,velocity_left);
    velocity_arm1 = set_arm1_velocity(arm1_turn);
    arm1->setVelocity(velocity_arm1);
    velocity_joint1 = set_joint1_velocity(joint1_turn);
    velocity_joint2 = set_joint1_velocity(joint2_turn);
    velocity_joint3 = set_joint1_velocity(joint3_turn);
    joint1->setVelocity(velocity_joint1);
    joint2->setVelocity(velocity_joint2);
    joint3->setVelocity(velocity_joint3);
    wheels[0]->setVelocity(velocity_left);
    wheels[1]->setVelocity(velocity_right);
  };

  // Enter here exit cleanup code.

  delete robot;
  delete arm1;
  delete joint1;
  delete joint2;
  delete joint3;
  delete[] wheels;
  return 0;
}
